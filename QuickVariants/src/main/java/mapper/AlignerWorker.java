package mapper;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Queue;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;

public class AlignerWorker extends Thread {
  static Logger silentLogger = Logger.NoOpLogger;

  public AlignerWorker(SequenceDatabase reference, int workerId, List<AlignmentListener> resultsListeners, AlignmentCache resultsCache, Queue completionListener) {
    this.sequenceDatabase = reference;
    this.workerId = "" + workerId;
    while (this.workerId.length() < 5) {
      this.workerId = " " + this.workerId;
    }
    this.resultsListeners = resultsListeners;
    this.completionListener = completionListener;
  }

  // queues the given work to execute later
  public void requestProcess(List<List<SamAlignment_Builder>> queryGroups, long startMillis, Logger alignmentLogger, Logger referenceLogger) {
    this.resetStatistics();

    this.logger = alignmentLogger;
    this.referenceLogger = referenceLogger;
    this.detailedAlignmentLogger = alignmentLogger.incrementScope();
    if (this.logger.getEnabled()) {
      log("\nOutput from worker " + this.workerId + ":");
    }
    this.groupedQueries = queryGroups;

    try {
      this.workQueue.put(true);
    } catch (InterruptedException e) {
      throw new IllegalArgumentException("Worker "  + this.workerId + " has no capacity for more work");
    }
  }

  // runs the given work in the current thread
  public void process(List<List<SamAlignment_Builder>> queryGroups, Logger logger) {
    requestProcess(queryGroups, 0, logger, logger);
    process();
  }

  private void resetStatistics() {
    numCacheHits = 0;
    numCacheMisses = 0;

    numCasesImmediatelyAcceptingFirstAlignment = 0;
  }

  public void noMoreQueries() {
    try {
      this.workQueue.put(false);
    } catch (InterruptedException e) {
      throw new IllegalArgumentException("Worker " + this.workerId + " is still working");
    }
  }

  @Override
  public void run() {
    while (true) {
      boolean succeeded = false;
      Boolean moreWork = false;
      try {
        moreWork = this.workQueue.take();
      } catch (InterruptedException e) {
        System.out.println("Interrupted");
        break;
      }
      if (!moreWork)
        break;
      //System.out.println("Worker.run got " + queries.size() + " queries");
      try {
        this.process();
        succeeded = true;
      } finally {
        // record error if any
        if (!succeeded)
          this.failed = true;
        this.completionListener.add(this);
      }
    }
  }

  private void process() {
    List<List<SamAlignment>> groupedQueries = this.buildQueries(this.groupedQueries);
    // List that for each query says where it aligns
    List<QueryAlignments> alignments = new ArrayList<QueryAlignments>();
    for (List<SamAlignment> queryAlignments: groupedQueries) {
      QueryAlignments alignmentsHere;
      try {
        alignmentsHere = this.align(queryAlignments);
      } catch (Exception e) {
        throw new RuntimeException("Failed to process " + queryAlignments, e);
      }
      // update some timing information
      // collect results
      alignments.add(alignmentsHere);
      if (alignmentsHere.getNumComponents() > 0) {
        if (this.logger.getEnabled()) {
          for (List<QueryAlignment> component: alignmentsHere.getAlignments()) {
            this.printAlignment(component);
          }
        }
      } else {
        if (this.logger.getEnabled()) {
          for (Sequence querySequence: alignmentsHere.getSequences()) {
            log("Unaligned    : " + querySequence.format());
          }
        }
      }
      if (this.logger.getEnabled()) {
        log(" ");
      }
    }
    this.sendResults(alignments);
  }

  public boolean tryComplete() throws InterruptedException {
    this.referenceLogger.flush();
    this.logger.flush();
    return !this.failed;
  }

  private void log(String message) {
    this.logger.log(message);
  }

  private List<List<SamAlignment>> buildQueries(List<List<SamAlignment_Builder>> queries) {
    List<List<SamAlignment>> groups = new ArrayList<List<SamAlignment>>();
    List<SamAlignment> currentGroup = null;
    for (List<SamAlignment_Builder> builderGroup: queries) {
      List<SamAlignment> group = new ArrayList<SamAlignment>(builderGroup.size());
      for (SamAlignment_Builder builder: builderGroup) {
        group.add(builder.build());
      }
      groups.add(group);
    }
    return groups;
  }

  // aligns to the unmodified reference we've been given
  public QueryAlignments align(List<SamAlignment> queries) {
    List<Sequence> queryComponents = new ArrayList<Sequence>();
    List<QueryAlignment> results = new ArrayList<QueryAlignment>(queries.size());
    for (SamAlignment query: queries) {
      QueryAlignment converted = tryConvertSamAlignment(query);
      if (converted == null)
        throw new IllegalArgumentException("Not a sam query: " + query);
      // If we have multiple alignments for the same query, and some alignments split the query into more pieces than others, we only keep the alignments using the most sequences
      boolean keep = converted.getComponents().size() >= queryComponents.size();
      boolean clear = converted.getComponents().size() != queryComponents.size();
      if (clear) {
        results.clear();
        queryComponents.clear();
        for (SequenceAlignment sequenceAlignment: converted.getComponents()) {
          queryComponents.add(sequenceAlignment.getSequenceA());
        }
      }
      if (keep) {
        results.add(converted);
      }
    }
    return QueryAlignments.singleComponent(queryComponents, results);
  }

  private QueryAlignment tryConvertSamAlignment(SamAlignment query) {
    List<SequenceAlignment> sequenceAlignments = new ArrayList<SequenceAlignment>(query.getNumSequences());
    double combinedScore = 0;
    for (Sequence sequence: query.getSequences()) {
      if (sequence instanceof SamRecord) {
        SamRecord samRecord = (SamRecord)sequence;
        SequenceAlignment sequenceAlignment = samRecord.toSequenceAlignment(this.sequenceDatabase);
        if (sequenceAlignment == null)
          return null;
        if (combinedScore == 0)
          combinedScore = samRecord.combinedScore;
        sequenceAlignments.add(sequenceAlignment);
      } else {
        return null;
      }
    }
    double combinedPenalty = -combinedScore;
    return new QueryAlignment(sequenceAlignments, 0, 0, 0, combinedPenalty, 0);
  }

  void printAlignment(List<QueryAlignment> alignments) {
    for (QueryAlignment alignment: alignments) {
      for (SequenceAlignment component : alignment.getComponents()) {
        this.printAlignment(component);
      }
    }
  }

  void printAlignment(SequenceAlignment alignment) {
    Sequence query = alignment.getSequenceA();

    int alignmentLength = alignment.getALength();

    String alignedQuery = alignment.getAlignedTextA();

    String alignedAncestralRef = alignment.getAlignedTextBHistory();
    String alignedUnmutatedRef = alignment.getAlignedTextB();

    String queryText = query.getText();
    String expectedAlignedText = queryText;
    if (alignment.isReferenceReversed()) {
      String originalQueryText = query.reverseComplement().getText();
      log("        Query: " + originalQueryText);
      log("     RC Query: " + queryText);
    } else {
      log("        Query: " + queryText);
    }

    if (!queryText.equals(alignedQuery)) {
      // If printing the aligned query is different from printing the query, then also print
      // the alignment of the query
      log("Aligned query: " + alignedQuery);
    }
    if (!alignedQuery.equals(alignedAncestralRef)) {
      StringBuilder differenceBuilder = new StringBuilder();
      differenceBuilder.append("Difference   : ");
      int max = Math.min(alignedQuery.length(), alignedAncestralRef.length());
      for (int i = 0; i < max; i++) {
        char c1 = alignedQuery.charAt(i);
        char c2 = alignedAncestralRef.charAt(i);
        if (c1 == c2) {
          differenceBuilder.append(" ");
        } else {
          if (Basepairs.canMatch(Basepairs.encode(c1), Basepairs.encode(c2))) {
            differenceBuilder.append("~");
          } else {
            differenceBuilder.append("!");
          }
        }
      }
      log(differenceBuilder.toString());
    }
    if (!alignedAncestralRef.equals(alignedUnmutatedRef)) {
      // If the ancestor analysis had an effect here, explain that too
      log("Ancestral ref: " + alignedAncestralRef + "(" + alignment.getSequenceBHistory().getName() + ", offset " + alignment.getStartOffset() + ")");
      log("Original ref : " + alignedUnmutatedRef + "(" + alignment.getSequenceB().getName() + ", offset " + alignment.getStartOffset() + ")");
    } else {
      log("Aligned ref  : " + alignedUnmutatedRef + "(" + alignment.getSequenceB().getName() + ", offset " + alignment.getStartOffset() + ")");
    }
  }

  private void sendResults(List<QueryAlignments> results) {
    for (AlignmentListener listener : this.resultsListeners) {
      listener.addAlignments(results);
    }
  }

  SequenceDatabase sequenceDatabase;
  List<AlignmentListener> resultsListeners;
  Logger logger;
  Logger detailedAlignmentLogger;
  Logger referenceLogger;
  String workerId;
  boolean failed = false;
  List<SequenceMatch> emptyMatchList = new ArrayList<SequenceMatch>(0);
  int numCacheHits;
  int numCacheMisses;

  int numCasesImmediatelyAcceptingFirstAlignment;
  Queue<AlignerWorker> completionListener;
  List<List<SamAlignment_Builder>> groupedQueries;

  BlockingQueue<Boolean> workQueue = new ArrayBlockingQueue<Boolean>(1);
}
