package mapper;

import com.sun.management.HotSpotDiagnosticMXBean;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Queue;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import javax.management.MBeanServer;

public class Main {

  static Logger alignmentLogger;
  static Logger referenceLogger;
  static Logger silentLogger = Logger.NoOpLogger;

  public static void main(String[] args) throws IllegalArgumentException, FileNotFoundException, IOException, InterruptedException {
    MapperMetadata.setMainArguments(args);
    // load properties
    long startMillis = System.currentTimeMillis();
    Properties properties = new Properties();
    properties.load(Main.class.getResourceAsStream("/mapper.properties"));
    String version = properties.getProperty("mapper.version", "unknown");

    System.out.println("QuickVariants version " + version);

    // parse arguments
    List<String> referencePaths = new ArrayList<String>();
    List<GroupedQuery_Provider> queries = new ArrayList<GroupedQuery_Provider>();
    String outVcfPath = null;
    String outSamPath = null;
    String outUnalignedPath = null;
    String outAncestorPath = null;
    boolean vcfIncludeNonMutations = true;
    String outRefsMapCountPath = null;
    String outMutationsPath = null;
    MutationDetectionParameters mutationFilterParameters = MutationDetectionParameters.defaultFilter();
    MutationDetectionParameters vcfFilterParameters = MutationDetectionParameters.emptyFilter();
    int alignmentVerbosity = 0;
    int referenceVerbosity = 0;
    boolean allowNoOutput = false;
    boolean autoVerbose = false;

    double mutationPenalty = 1;
    double indelStart_penalty = 4.0/3.0;
    double indelExtension_penalty = 0.5;
    double maxErrorRate = 0.1;
    int maxNumMatches = Integer.MAX_VALUE;
    double max_penaltySpan = -1;

    int numThreads = 1;
    double queryEndFraction = 0.1;

    for (int i = 0; i < args.length; i++) {
      String arg = args[i];
      if ("--reference".equals(arg)) {
        referencePaths.add(args[i + 1]);
        i++;
        continue;
      }
      if ("--in-ordered-sam".equals(arg)) {
        String queryPath = args[i + 1];
        i++;
        GroupedQuery_Provider samParser = DataLoader.ParseSamAlignments(queryPath, false, false);
        queries.add(samParser);
        continue;
      }
      if ("--in-sam".equals(arg) || "--in-unordered-sam".equals(arg)) {
        String queryPath = args[i + 1];
        i++;
        GroupedQuery_Provider samParser = DataLoader.ParseSamAlignments(queryPath, false, true);
        queries.add(samParser);
        continue; 
      }
      if ("--out-vcf".equals(arg)) {
        outVcfPath = args[i + 1];
        i += 2;
        while (i < args.length) {
          arg = args[i];
          if ("--snp-threshold".equals(arg)) {
            vcfFilterParameters.minSNPTotalDepth = (float)Double.parseDouble(args[i + 1]);
            vcfFilterParameters.minSNPDepthFraction = (float)Double.parseDouble(args[i + 2]);
            i += 3;
            continue;
          }
          if ("--indel-start-threshold".equals(arg)) {
            vcfFilterParameters.minIndelTotalStartDepth = (float)Double.parseDouble(args[i + 1]);
            vcfFilterParameters.minIndelStartDepthFraction = (float)Double.parseDouble(args[i + 2]);
            i += 3;
            continue;
          }
          if ("--indel-continue-threshold".equals(arg)) {
            vcfFilterParameters.minIndelContinuationTotalDepth = (float)Double.parseDouble(args[i + 1]);
            vcfFilterParameters.minIndelContinuationDepthFraction = (float)Double.parseDouble(args[i + 2]);
            i += 3;
            continue;
          }
          if ("--indel-threshold".equals(arg)) {
            vcfFilterParameters.minIndelTotalStartDepth = (float)Double.parseDouble(args[i + 1]);
            vcfFilterParameters.minIndelContinuationTotalDepth = (float)Double.parseDouble(args[i + 1]);

            vcfFilterParameters.minIndelStartDepthFraction = (float)Double.parseDouble(args[i + 2]);
            vcfFilterParameters.minIndelContinuationDepthFraction = (float)Double.parseDouble(args[i + 2]);
            i += 3;
            continue;
          }
          i--;
          // maybe this argument is a top-level argument
          break;
        }

        continue;
      }
      if ("--out-sam".equals(arg)) {
        outSamPath = args[i + 1];
        i++;
        continue;
      }
      if ("--out-refs-map-count".equals(arg)) {
        outRefsMapCountPath = args[i + 1];
        i++;
        continue;
      }
      if ("--out-mutations".equals(arg)) {
        outMutationsPath = args[i + 1];
        i += 2;
        while (i < args.length) {
          arg = args[i];
          if ("--snp-threshold".equals(arg)) {
            mutationFilterParameters.minSNPTotalDepth = (float)Double.parseDouble(args[i + 1]);
            mutationFilterParameters.minSNPDepthFraction = (float)Double.parseDouble(args[i + 2]);
            i += 3;
            continue;
          }
          if ("--indel-start-threshold".equals(arg)) {
            mutationFilterParameters.minIndelTotalStartDepth = (float)Double.parseDouble(args[i + 1]);
            mutationFilterParameters.minIndelStartDepthFraction = (float)Double.parseDouble(args[i + 2]);
            i += 3;
            continue;
          }
          if ("--indel-continue-threshold".equals(arg)) {
            mutationFilterParameters.minIndelContinuationTotalDepth = (float)Double.parseDouble(args[i + 1]);
            mutationFilterParameters.minIndelContinuationDepthFraction = (float)Double.parseDouble(args[i + 2]);
            i += 3;
            continue;
          }
          if ("--indel-threshold".equals(arg)) {
            mutationFilterParameters.minIndelTotalStartDepth = (float)Double.parseDouble(args[i + 1]);
            mutationFilterParameters.minIndelContinuationTotalDepth = (float)Double.parseDouble(args[i + 1]);

            mutationFilterParameters.minIndelStartDepthFraction = (float)Double.parseDouble(args[i + 2]);
            mutationFilterParameters.minIndelContinuationDepthFraction = (float)Double.parseDouble(args[i + 2]);
            i += 3;
            continue;
          }
          i--;
          // maybe this argument is a top-level argument
          break;
        }
        continue;
      }
      if ("--out-ancestor".equals(arg)) {
        outAncestorPath = args[i + 1];
        i++;
        continue;
      }
      if ("--no-output".equals(arg)) {
        allowNoOutput = true;
        continue;
      }
      if ("--verbose".equals(arg) || "-v".equals(arg)) {
        alignmentVerbosity = Math.max(alignmentVerbosity, 1);
        continue;
      }
      if ("--verbose-alignment".equals(arg)) {
        alignmentVerbosity = Math.max(alignmentVerbosity, Integer.MAX_VALUE);
        continue;
      }
      if ("--verbose-reference".equals(arg)) {
        referenceVerbosity = Math.max(referenceVerbosity, 1);
        continue;
      }
      if ("-vv".equals(arg)) {
        alignmentVerbosity = Math.max(alignmentVerbosity, Integer.MAX_VALUE);
        referenceVerbosity = Math.max(referenceVerbosity, 1);
        continue;
      }
      if ("--verbosity-auto".equals(arg)) {
        autoVerbose = true;
        continue;
      }
      if ("--num-threads".equals(arg)) {
        String value = args[i + 1];
        numThreads = Integer.parseInt(value);
        i++;
        continue;
      }
      if ("--distinguish-query-ends".equals(arg)) {
        String value = args[i + 1];
        queryEndFraction = Double.parseDouble(value);
        i++;
        continue;
      }
      if ("--vcf-exclude-non-mutations".equals(arg)) {
        vcfIncludeNonMutations = false;
        continue;
      }

      if ("--snp-threshold".equals(arg) || "--indel-start-threshold".equals(arg) || "--indel-continue-threshold".equals(arg) || "--indel-threshold".equals(arg)) {
        usageError("" + arg + " is not a top-level argument: try --out-mutations <mutations.txt> " + arg + " <min total depth> <min supporting depth fraction>");
      }

      usageError("Unrecognized argument: " + arg);
    }
    if (referencePaths.size() < 1) {
      usageError("--reference is required");
    }
    if (queries.size() < 1) {
      usageError("--in-sam is required");
    }
    if (outVcfPath == null && outSamPath == null && outRefsMapCountPath == null && outUnalignedPath == null && outMutationsPath == null && !allowNoOutput) {
      usageError("No output specified. Try --out-vcf <output path>, or if you really don't want to generate an output file, --no-output");
    }
    alignmentLogger = new Logger(new PrintWriter(), 1, alignmentVerbosity);
    referenceLogger = new Logger(new PrintWriter(), 1, referenceVerbosity);

    if (numThreads < 1) {
      usageError("--num-threads must be >= 1");
    }
    if (queryEndFraction < 0 || queryEndFraction >= 1) {
      usageError("--distinguish-query-ends must be >= 0 and < 1");
    }

    if (max_penaltySpan < 0) {
      max_penaltySpan = mutationPenalty / 2;
    }
    // extending an insertion by some number of base pairs reduces the number of base pairs matching in the query by the same number of base pairs
    double additional_insertionExtension_penalty = maxErrorRate;

    System.out.println("" + referencePaths.size() + " reference files:");
    for (String referencePath: referencePaths) {
      System.out.println("Reference path = " + referencePath);
    }
    System.out.println("" + queries.size() + " sets of queries: ");
    for (GroupedQuery_Provider groupBuilder : queries) {
      System.out.println(groupBuilder.toString());
    }
    boolean successful = run(referencePaths, queries, outVcfPath, vcfIncludeNonMutations, outSamPath, outRefsMapCountPath, outMutationsPath, mutationFilterParameters, vcfFilterParameters, outUnalignedPath, numThreads, queryEndFraction, autoVerbose, outAncestorPath, startMillis);
    if (!successful) {
      System.exit(1);
    }
  }

  public static void usageError(String message) {
    System.err.println(
"\n" +
"Usage:\n"+
"  java -jar quick-variants.jar [--out-vcf <out.vcf>] [--out-mutations <out.txt>] [--out-sam <out.sam>] [--out-refs-map-count <counts.txt>] --reference <ref.fasta> --in-sam <input.sam> [options]\n" +
"\n" +
"    Converts a sam file to other formats, most notably .vcf\n" +
"\n" +
"  INPUT:\n" +
"\n" +
"    --reference <file> the reference to use. Should be in .fasta/.fa/.fna or .fastq/.fq format or a .gz of one of those formats.\n" +
"\n" +
"    --in-ordered-sam <file.sam> the input alignments. If there are Illumina-style paired-end reads, alignments for each mate should be adjacent in the sam file\n" +
"    --in-unordered-sam <file.sam> same as --in-ordered-sam but doesn't require alignments for matest to be listed in adjacent lines in the file\n" +
"      This argument is slower than --in-unordered-sam and requires more memory\n" +
"    --in-sam <file.sam> alias for --in-unordered-sam <file.sam>\n" +
"\n" +
"  OUTPUT FORMATS:\n" +
"\n" +
"    Summary by reference position\n" +
"\n" +
"      --out-vcf <file> output file to generate containing a description of mutation counts by position\n" +
"      --vcf-exclude-non-mutations if set, the output vcf file will exclude positions where no mutations were detected\n" +
"      --distinguish-query-ends <fraction> (default 0.1) In the output vcf file, we separately display which queries aligned at each position with <fraction> of the end of the query and which didn't.\n" +
"\n" +
"      --snp-threshold <min total depth> <min supporting depth fraction> (default 0, 0)\n" +
"        The minimum total depth and minimum supporting depth fraction required at a position to report the support for the mutation\n" +
"\n" +
"      --indel-start-threshold <min total depth> <min supporting depth fraction> (default 0, 0)\n" +
"        The minimum total (middle) depth and minimum supporting depth fraction required at a position to report support for the start of an insertion or deletion\n" +
"\n" +
"      --indel-continue-threshold <min total depth> <min supporting depth fraction> (default 0, 0)\n" +
"        The minimum total (middle) depth and minimum supporting depth fraction required at a position to report support for a continuation of an insertion or deletion\n" +
"      --indel-threshold <min total depth> <min supporting depth fraction>\n" +
"        Alias for --indel-start-threshold <min total depth> <min supporting depth frequency> and --indel-continue-threshold <min total depth> <min supporting depth frequency>\n" +
"\n" +
"    Summary by mutation\n" +
"\n" +
"      --out-mutations <file> output file to generate listing the mutations of the queries compared to the reference genome\n" +
"\n" +
"      --distinguish-query-ends <fraction> (default 0.1) When detecting indels, only consider the middle <fraction> of each query\n" +
"\n" +
"      --snp-threshold <min total depth> <min supporting depth fraction> (default 5, 0.9)\n" +
"        The minimum total depth and minimum supporting depth fraction required at a position to report it as a point mutation\n" +
"\n" +
"      --indel-start-threshold <min total depth> <min supporting depth fraction> (default 1, 0.8)\n" +
"        The minimum total (middle) depth and minimum supporting depth fraction required at a position to report it as the start of an insertion or deletion\n" +
"\n" +
"      --indel-continue-threshold <min total depth> <min supporting depth fraction> (default 1, 0.7)\n" +
"        The minimum total (middle) depth and minimum supporting depth fraction required at a position to report it as a continuation of an insertion or deletion\n" +
"      --indel-threshold <min total depth> <min supporting depth fraction>\n" +
"        Alias for --indel-start-threshold <min total depth> <min supporting depth frequency> and --indel-continue-threshold <min total depth> <min supporting depth frequency>\n" +
"\n" +
"    Summary by genome\n" +
"\n" +
"      --out-refs-map-count <file> the output file to summarize the number of reads mapped to each combination of references\n" +
"\n" +
"    Raw output\n" +
"\n" +
"      --out-sam <file> the output file in SAM format\n" +
"\n" +
"    --no-output if no output is requested, skip writing output rather than throwing an error\n" +
"\n" +
"    Debugging\n" +
"\n" +
"      -v, --verbose output diagnostic information\n" +
"\n" +
"      -vv more diagnostic information\n" +
"\n" +
"      --verbosity-auto set verbosity flags dynamically and automatically based on the data and alignment\n" +
"\n" +
"  Multiple output formats may be specified during a single run; for example:\n" +
"\n" +
"    --out-sam out.sam --out-vcf out.vcf\n" +
"\n" +
"  OTHER:\n" +
"\n" +
"    --num-threads <count> number of threads to use at once for processing. Higher values will run more quickly on a system that has that many CPUs available.\n"
);

    System.err.println(message);
    System.exit(1);
  }

  // performs alignment and outputs results
  public static boolean run(List<String> referencePaths, List<GroupedQuery_Provider> queriesList, String outVcfPath, boolean vcfIncludeNonMutations, String outSamPath, String outRefsMapCountPath, String outMutationsPath, MutationDetectionParameters mutationFilterParameters, MutationDetectionParameters vcfFilterParameters, String outUnalignedPath, int numThreads, double queryEndFraction, boolean autoVerbose, String outAncestorPath, long startMillis) throws IllegalArgumentException, FileNotFoundException, IOException, InterruptedException {
    VcfWriter vcfWriter = null;
    if (outVcfPath != null)
      vcfWriter = new VcfWriter(outVcfPath, vcfIncludeNonMutations, vcfFilterParameters);

    MutationsWriter mutationsWriter = null;
    if (outMutationsPath != null)
      mutationsWriter = new MutationsWriter(outMutationsPath, mutationFilterParameters);

    System.out.println("Loading reference");
    boolean keepQualityData = (outUnalignedPath != null);
    SequenceProvider reference = DataLoader.LoadFrom(referencePaths, false);
    List<Sequence> sortedReference = sortAndComplementReference(reference);
    SequenceDatabase sequenceDatabase = new SequenceDatabase(sortedReference);
    long now = System.currentTimeMillis();
    long elapsed = (now - startMillis) / 1000;

    GroupedQuery_Provider queries = new Composite_GroupedQuery_Provider(queriesList);

    List<AlignmentListener> listeners = new ArrayList<AlignmentListener>();
    MatchDatabase matchDatabase = new MatchDatabase(queryEndFraction);
    ReferenceAlignmentCounter referenceAlignmentCounter = new ReferenceAlignmentCounter();
    AlignmentCounter matchCounter = new AlignmentCounter();
    if (outVcfPath != null || outMutationsPath != null) {
      listeners.add(matchDatabase);
    }
    SamWriter samWriter = null;
    if (outSamPath != null) {
      samWriter = new SamWriter(sequenceDatabase, outSamPath);
      listeners.add(samWriter);
    }
    UnalignedQuery_Writer unalignedWriter = null;
    if (outUnalignedPath != null) {
      unalignedWriter = new UnalignedQuery_Writer(outUnalignedPath, queries.get_allReadsContainQualityInformation());
      listeners.add(unalignedWriter);
    }
    if (outRefsMapCountPath != null) {
      listeners.add(referenceAlignmentCounter);
    }
    listeners.add(matchCounter);
    AlignmentStatistics statistics = compare(sequenceDatabase, queries, startMillis, numThreads, queryEndFraction, listeners, autoVerbose);

    long numMatchingQuerySequences = matchCounter.getNumMatchingSequences();
    long numQuerySequences = matchCounter.getNumSequences();
    long matchPercent;
    if (numQuerySequences > 0)
      matchPercent = numMatchingQuerySequences * 100 / numQuerySequences;
    else
      matchPercent = 0;
    long totalAlignedQueryLength = matchCounter.getTotalAlignedQueryLength();

    // output referenceAlignmentCounter RefsMapCount
    if (outRefsMapCountPath != null) {
      long computationEnd = System.currentTimeMillis();
      long computationTime = (computationEnd - startMillis) / 1000;
      System.out.println("Writing RefsMapCount results at " + computationTime + "s");
      referenceAlignmentCounter.sumAlignments(outRefsMapCountPath);
      long writingEnd = System.currentTimeMillis();
      long writingTime = (writingEnd - startMillis) / 1000;
      System.out.println("Saved " + outRefsMapCountPath + " at " + writingTime + "s");
    }
    String displayCoverage = null;
    if (outVcfPath != null) {
      Map<Sequence, Alignments> alignments = matchDatabase.groupByPosition();
      long computationEnd = System.currentTimeMillis();
      double computationTime = (double)(computationEnd - startMillis) / 1000.0;
      System.out.println("Writing vcf results at " + computationTime + "s");
      vcfWriter.write(alignments, numThreads);
      System.out.println("Saved " + outVcfPath);
      long numMatchedPositions = vcfWriter.getNumReferencePositionsMatched();
      long numPositions = sequenceDatabase.getTotalForwardSize();

      // Format the coverage as a human-readable string
      double coverage = ((double)numMatchedPositions) / ((double)numPositions);
      displayCoverage = "" + (int)(coverage * 100) + "%";
      // If the coverage is less than 1% but more than 0%, emphasize that to make it easy to notice
      if (displayCoverage.equals("0%") && coverage > 0) {
        displayCoverage = "<1%";
      }
      displayCoverage = " Coverage                      : " + displayCoverage + " of the reference (" + numMatchedPositions + "/" + numPositions + ") was matched";
    }
    if (outMutationsPath != null) {
      Map<Sequence, Alignments> alignments = matchDatabase.groupByPosition();
      long computationEnd = System.currentTimeMillis();
      double computationTime = (double)(computationEnd - startMillis) / 1000.0;
      System.out.println("Writing mutation results at " + computationTime + "s");
      mutationsWriter.write(alignments, numThreads);
      System.out.println("Saved " + outMutationsPath);
    }
    // show statistics
    System.out.println("");
    System.out.println("Statistics: ");
    if (matchCounter.getNumMatchingSequences() != matchCounter.getNumAlignedQueries()) {
      // paired-end reads
      Distribution distance = matchCounter.getDistanceBetweenQueryComponents();
      System.out.println(" Query pair separation distance: avg: " + (float)distance.getMean() + " stddev: " + (float)distance.getStdDev());
    }
    System.out.println(" Alignment rate                : " + matchPercent + "% of query sequences (" + numMatchingQuerySequences + "/" + numQuerySequences + ")");
    if (displayCoverage != null) {
      System.out.println(displayCoverage);
    }
    //System.out.println(" Average penalty               : " + averagePenaltyPerBase + " per base (" + (long)totalAlignedPenalty + "/" + (long)totalAlignedQueryLength + ") in aligned queries");
    if (statistics != null) {
      System.out.println("");
      System.out.println("Timing:");

      Query slowestQuery = statistics.slowestQuery;

      if (slowestQuery != null) {
        String queryDisplayText = slowestQuery.format();
        String numAlignmentsText;
        int numAlignments = statistics.slowestQueryNumAlignments;
        if (numAlignments == 1)
          numAlignmentsText = "1 time";
        else
          numAlignmentsText = "" + numAlignments + " times";
        System.out.println(" Slowest query: #" + slowestQuery.getId() + " (" + statistics.slowestQueryMillis + "ms) : " + queryDisplayText + " aligned " + numAlignmentsText);
      }

      //int millisOnUnalignedQueries = (int)(statistics.cpuMillisSpentOnUnalignedQueries / 1000 / numThreads);
      //System.out.println(" Unaligned queries took        : " + statistics.cpuMillisSpentOnUnalignedQueries + " cpu-ms (" + millisOnUnalignedQueries + "s)");

      System.out.println(" Time reading queries          : " + statistics.millisReadingQueries + "ms");
      System.out.println(" Time launching workers        : " + statistics.millisLaunchingWorkers + "ms");
      System.out.println(" Time waiting for workers      : " + statistics.millisWaitingForWorkers + "ms");
    }
    boolean fullySuccessful = true;
    String successStatusMessage = "Done";
    if (statistics == null) {
      successStatusMessage = "Failed";
      fullySuccessful = false;
    } else {
      int numErrors = queries.getNumErrors();
      if (numErrors > 0) {
        successStatusMessage = "Completed with " + numErrors + " errors";
        fullySuccessful = false;
      }
    }
    if (samWriter != null)
      samWriter.close();
    if (unalignedWriter != null)
      unalignedWriter.close();
    long end = System.currentTimeMillis();
    double totalTime = ((double)end - (double)startMillis) / 1000.0;
    System.out.println("");
    System.out.println(successStatusMessage + " in " + totalTime + "s.");
    return fullySuccessful;
  }

  public static void dumpHeap() throws IOException {
    String outputPath = "mapper.hprof";
    System.out.println("dumping heap to " + outputPath);
    MBeanServer server = ManagementFactory.getPlatformMBeanServer();
    HotSpotDiagnosticMXBean bean =
        ManagementFactory.newPlatformMXBeanProxy(server,
        "com.sun.management:type=HotSpotDiagnostic", HotSpotDiagnosticMXBean.class);
    bean.dumpHeap("mapper.hprof", false);
    System.out.println("dumped heap to " + outputPath);
  }

  public static AlignmentStatistics compare(SequenceDatabase reference, GroupedQuery_Provider queries, long startMillis, int numThreads, double queryEndFraction, List<AlignmentListener> alignmentListeners, boolean autoVerbose) throws InterruptedException, IOException {
    long readingMillis = 0;
    long launchingMillis = 0;
    long waitingMillis = 0;
    int numQueriesForNextMessage = 1;
    long previousElapsedSeconds = 0;

    // Create some workers and assign some queries to each
    Set<AlignerWorker> workers = new HashSet<AlignerWorker>(numThreads);
    AlignmentCache alignmentCache = new AlignmentCache();

    long numQueriesLoaded = 0;
    int maxNumBasesPerJob = 500000;
    int workerIndex = 0;
    boolean doneReadingQueries = false;
    // pendingQueries[jobIndex][groupNumber][alignmentIndex] = the corresponding query builder
    List<List<List<QueryBuilder>>> pendingQueries = new ArrayList<List<List<QueryBuilder>>>();
    long lastPrintTime = 0;
    long nextCountToPrint = 0;
    long slowestAlignmentMillis = -1;
    Query slowestQuery = null;
    List<QueryAlignment> slowestAlignment = null;
    long millisSpentOnUnalignedQueries = 0;
    int numCacheHits = 0;
    int numCasesImmediatelyAcceptingFirstAlignment = 0;
    BlockingQueue<AlignerWorker> completedWorkers = new ArrayBlockingQueue<AlignerWorker>(numThreads);
    boolean everSaturatedWorkers = false;
    List<AlignerWorker> pendingWorkers = new ArrayList<AlignerWorker>();
    while (workers.size() > 0 || !doneReadingQueries || pendingQueries.size() > 0) {
      boolean progressed = false;
      if (workers.size() >= numThreads)
        everSaturatedWorkers = true;
      // Do we have to read more queries?
      if (!doneReadingQueries && (pendingQueries.size() < 1 || (workers.size() >= numThreads && pendingQueries.size() < numThreads * 10))) {
        long readStart = System.currentTimeMillis();
        int targetNumBases = maxNumBasesPerJob;
        List<List<QueryBuilder>> batch = new ArrayList<List<QueryBuilder>>();
        int totalLengthOfPendingQueries = 0;
        while (true) {
          if (totalLengthOfPendingQueries >= targetNumBases) {
            break;
          }
          List<QueryBuilder> group = queries.getNextGroup();
          if (group == null) {
            doneReadingQueries = true;
            break;
          }
          for (QueryBuilder queryBuilder: group) {
            numQueriesLoaded++;
            queryBuilder.setId(numQueriesLoaded);
            totalLengthOfPendingQueries += queryBuilder.getLength();
          }
          batch.add(group);
        }
        pendingQueries.add(batch);
        progressed = true;
        long readEnd = System.currentTimeMillis();
        readingMillis += (readEnd - readStart);
      }

      // if we haven't created all workers, create one
      if (workers.size() < numThreads) {
        long launchStart = System.currentTimeMillis();

        List<List<QueryBuilder>> queriesToProcess;
        // we have enough idle threads to spawn another worker
        if (pendingQueries.size() > 0) {
          // we have queries that haven't been assigned
          queriesToProcess = pendingQueries.get(pendingQueries.size() - 1);
          pendingQueries.remove(pendingQueries.size() - 1);
        } else {
          // each query has been given to some worker
          queriesToProcess = null;
        }
        if (queriesToProcess != null) {
          long now = System.currentTimeMillis();
          long elapsed = (now - startMillis) / 1000;

          // give these queries to this worker
          Logger workerAlignmentLogger = alignmentLogger.withWriter(new BufferedWriter());
          Logger workerReferenceLogger = referenceLogger.withWriter(new BufferedWriter());
          if (autoVerbose && workerIndex == 0) {
            workerAlignmentLogger = new Logger(new BufferedWriter(), 1, Integer.MAX_VALUE);
          }
          AlignerWorker worker;
          boolean workerAlreadyRunning;
          if (pendingWorkers.size() > 0) {
            worker = pendingWorkers.remove(pendingWorkers.size() - 1);
            workerAlreadyRunning = true;
          } else {
            worker = new AlignerWorker(reference, workerIndex, alignmentListeners, alignmentCache, completedWorkers);
            workerAlreadyRunning = false;
          }
          workers.add(worker);
          worker.requestProcess(queriesToProcess, startMillis, workerAlignmentLogger, workerReferenceLogger);
          workerIndex++;
          progressed = true;
          if (!workerAlreadyRunning)
            worker.start();

          // determine if enough queries have completed and enough time has passed for it to be worth issuing a status update
          if (numQueriesLoaded >= nextCountToPrint) {
            nextCountToPrint = determineNextCountToReport(numQueriesLoaded);
            if (elapsed != lastPrintTime) {
              // note elapsed != 0 because lastPrintTime starts at 0
              long queriesPerSecond = numQueriesLoaded / elapsed;
              System.out.println("Starting to process query " + numQueriesLoaded + " at " + elapsed + "s (" + queriesPerSecond + " q/s), " + workers.size() + " active workers, " + pendingQueries.size() + " ready jobs");
              checkMemoryUsage();
              lastPrintTime = elapsed;
            }
          }
        }
        long launchEnd = System.currentTimeMillis();
        launchingMillis += (launchEnd - launchStart);
      }

      long waitStart = System.currentTimeMillis();
      while (!progressed || (everSaturatedWorkers && completedWorkers.peek() != null)) {
        // process any workers that completed
        AlignerWorker worker = completedWorkers.take();
        boolean succeeded = worker.tryComplete();
        if (succeeded == false) {
          System.out.println("Worker failed; aborting");
          return null;
        }

        // remove this worker
        workers.remove(worker);
        pendingWorkers.add(worker);
        progressed = true;
      }
      long waitEnd = System.currentTimeMillis();
      waitingMillis += (waitEnd - waitStart);
    }
    for (AlignerWorker worker: pendingWorkers) {
      worker.noMoreQueries();
    }
    long doneAligningQueriesAt = System.currentTimeMillis();
    AlignmentStatistics result = new AlignmentStatistics();
    result.millisReadingQueries = readingMillis;
    result.millisLaunchingWorkers = launchingMillis;
    result.millisWaitingForWorkers = waitingMillis;
    result.cpuMillisSpentOnUnalignedQueries = millisSpentOnUnalignedQueries;
    result.numCasesImmediatelyAcceptingFirstAlignment = numCasesImmediatelyAcceptingFirstAlignment;
    result.numQueriesLoaded = numQueriesLoaded;

    return result;
  }

  private static void checkMemoryUsage() {
    Runtime runtime = Runtime.getRuntime();
    long maxAllowedMemory = runtime.maxMemory();
    long usedMemory = runtime.totalMemory() - runtime.freeMemory();
    double usageFraction = (double)usedMemory / (double)maxAllowedMemory;
    if (usageFraction >= 0.9) {
      long maxAllowedMemoryMB = maxAllowedMemory / 1024 / 1024;
      long usedMemoryMB = usedMemory / 1024 / 1024;
      int usagePercent = (int)(usageFraction * 100);
      System.out.println("Low memory! " + usagePercent + "% used (" + usedMemoryMB + "M/" + maxAllowedMemoryMB + "M). Try larger Xmx");
    }
  }

  public static List<Sequence> sortAndComplementReference(SequenceProvider provider) {
    Map<Integer, List<Sequence>> sequencesByLength = new TreeMap<Integer, List<Sequence>>();
    while (true) {
      SequenceBuilder builder = provider.getNextSequence();
      if (builder == null)
        break;
      Sequence sequence = builder.build();
      int key = sequence.getLength() * -1;
      List<Sequence> sequencesHere = sequencesByLength.get(key);
      if (sequencesHere == null) {
        sequencesHere = new ArrayList<Sequence>();
        sequencesByLength.put(key, sequencesHere);
      }
      sequencesHere.add(sequence);
      sequencesHere.add(sequence.reverseComplement());
    }
    List<Sequence> sorted = new ArrayList<Sequence>();
    for (int length : sequencesByLength.keySet()) {
      sorted.addAll(sequencesByLength.get(length));
    }
    return sorted;
  }

  // Given that we've completed <count> units of work, tells when to report again
  public static long determineNextCountToReport(long count) {
    // Only output numbers ending with lots of zeros and starting with up to two nonzero digits
    long multiplier = 1;
    while (count > 99) {
      count /= 10;
      multiplier *= 10;
    }
    return (count + 1) * multiplier;
  }
}
