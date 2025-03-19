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
  public void requestProcess(List<List<QueryBuilder>> queryGroups, long startMillis, Logger alignmentLogger, Logger referenceLogger) {
    this.resetStatistics();

    this.startMillis = startMillis;
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
  public void process(List<List<QueryBuilder>> queryGroups, Logger logger) {
    requestProcess(queryGroups, 0, logger, logger);
    process();
  }

  private void resetStatistics() {
    numCacheHits = 0;
    numCacheMisses = 0;

    slowestAlignment = null;
    slowestAlignmentMillis = -1;
    millisSpentOnUnalignedQueries = 0;
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
    List<List<Query>> groupedQueries = this.buildQueries(this.groupedQueries);
    // List that for each query says where it aligns
    List<List<QueryAlignment>> alignments = new ArrayList<List<QueryAlignment>>(groupedQueries.size());
    List<Query> unalignedQueries = new ArrayList<Query>();
    for (List<Query> queryAlignments: groupedQueries) {
      long start = System.currentTimeMillis();
      List<QueryAlignment> alignmentsHere;
      try {
        alignmentsHere = this.align(queryAlignments);
      } catch (Exception e) {
        throw new RuntimeException("Failed to process " + queryAlignments, e);
      }
      // update some timing information
      long end = System.currentTimeMillis();
      long elapsed = end - start;
      if (elapsed > this.slowestAlignmentMillis) {
        this.slowestAlignmentMillis = (int)elapsed;
        this.slowestAlignment = alignmentsHere;
      }
      // collect results
      if (alignmentsHere.size() > 0) {
        alignments.add(alignmentsHere);
        if (this.logger.getEnabled()) {
          this.printAlignment(queryAlignments.get(0), alignmentsHere);
        }
      } else {
        unalignedQueries.addAll(queryAlignments);
        if (this.logger.getEnabled()) {
          for (Query query: queryAlignments) {
            log("Unaligned    : " + query.format());
          }
        }
      }
      if (this.logger.getEnabled()) {
        log(" ");
      }
    }
    this.sendResults(alignments, unalignedQueries);
  }

  public boolean tryComplete() throws InterruptedException {
    this.referenceLogger.flush();
    this.logger.flush();
    return !this.failed;
  }

  private void log(String message) {
    this.logger.log(message);
  }

  private List<List<Query>> buildQueries(List<List<QueryBuilder>> queries) {
    List<List<Query>> groups = new ArrayList<List<Query>>();
    List<Query> currentGroup = null;
    for (List<QueryBuilder> builderGroup: queries) {
      List<Query> group = new ArrayList<Query>(builderGroup.size());
      for (QueryBuilder builder: builderGroup) {
        group.add(builder.build());
      }
      groups.add(group);
    }
    return groups;
  }

  // aligns to the unmodified reference we've been given
  public List<QueryAlignment> align(List<Query> queries) {
    List<QueryAlignment> results = new ArrayList<QueryAlignment>(queries.size());
    for (Query query: queries) {
      QueryAlignment converted = tryConvertSamAlignment(query);
      if (converted == null)
        throw new IllegalArgumentException("Not a sam query: " + query);
      results.add(converted);
    }
    return results;
  }

  private QueryAlignment tryConvertSamAlignment(Query query) {
    List<SequenceAlignment> sequenceAlignments = new ArrayList<SequenceAlignment>(query.getNumSequences());
    for (Sequence sequence: query.getSequences()) {
      if (sequence instanceof SamRecord) {
        SamRecord samRecord = (SamRecord)sequence;
        SequenceAlignment sequenceAlignment = samRecord.toSequenceAlignment(this.sequenceDatabase);
        if (sequenceAlignment == null)
          return null;
        sequenceAlignments.add(sequenceAlignment);
      } else {
        return null;
      }
    }
    return new QueryAlignment(sequenceAlignments, 0, 0, 0);
  }

  void printAlignment(Query query, List<QueryAlignment> alignments) {
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

  private void sendResults(List<List<QueryAlignment>> results, List<Query> unalignedQueries) {
    for (AlignmentListener listener : this.resultsListeners) {
      listener.addAlignments(results);
      listener.addUnaligned(unalignedQueries);
    }
  }

  SequenceDatabase sequenceDatabase;
  List<AlignmentListener> resultsListeners;
  Logger logger;
  Logger detailedAlignmentLogger;
  Logger referenceLogger;
  String workerId;
  long startMillis;
  boolean failed = false;
  List<SequenceMatch> emptyMatchList = new ArrayList<SequenceMatch>(0);
  int numCacheHits;
  int numCacheMisses;

  Query slowestQuery;
  List<QueryAlignment> slowestAlignment;
  int slowestAlignmentMillis = -1;
  long millisSpentOnUnalignedQueries;
  int numCasesImmediatelyAcceptingFirstAlignment;
  Queue<AlignerWorker> completionListener;
  List<List<QueryBuilder>> groupedQueries;

  BlockingQueue<Boolean> workQueue = new ArrayBlockingQueue<Boolean>(1);
}
