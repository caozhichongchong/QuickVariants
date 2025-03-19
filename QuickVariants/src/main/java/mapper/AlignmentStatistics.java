package mapper;

public class AlignmentStatistics {

  long millisReadingQueries;
  long millisLaunchingWorkers;
  long millisWaitingForWorkers;

  long cpuMillisSpentOnUnalignedQueries;

  int slowestQueryNumAlignments;
  long slowestQueryMillis;

  long numCasesImmediatelyAcceptingFirstAlignment;
  long numQueriesLoaded;
}
