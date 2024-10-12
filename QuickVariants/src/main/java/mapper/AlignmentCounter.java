package mapper;

import java.util.List;
import java.util.Map;

public class AlignmentCounter implements AlignmentListener {
  public void addAlignments(List<QueryAlignments> alignments) {
    double newTotalAlignedPenalty = 0;
    long newTotalAlignedQueryLength = 0;
    int numNewAlignedQueries = 0;
    int numNewUnalignedQueries = 0;
    Distribution newTotalDistanceBetweenComponents = new Distribution();
    for (QueryAlignments queryAlignments : alignments) {
      for (List<QueryAlignment> choices: queryAlignments.getAlignments()) {
        if (choices.size() > 0) {
          numNewAlignedQueries++;

          newTotalAlignedPenalty += choices.get(0).getPenalty();
          newTotalAlignedQueryLength += choices.get(0).getALength();

          double currentTotalDistanceBetweenComponents = 0;
          for (QueryAlignment choice: choices) {
            if (choice.getNumSequences() > 1)
              newTotalDistanceBetweenComponents.add(choice.getTotalDistanceBetweenComponents(), (double)1.0 / (double)choices.size());
          }
        } else {
          if (queryAlignments.getNumComponents() == 1) {
            numNewUnalignedQueries += 1;
          } else {
            // we have a partially aligned query, which we don't count
          }
        }

      }
    }
    synchronized (this) {
      this.numAlignedQueries += numNewAlignedQueries;
      this.numUnalignedQueries += numNewUnalignedQueries;
      this.totalAlignedPenalty += newTotalAlignedPenalty;
      this.totalAlignedQueryLength += newTotalAlignedQueryLength;
      this.distanceBetweenQueryComponents = this.distanceBetweenQueryComponents.plus(newTotalDistanceBetweenComponents);
    }
  }

  public long getNumQueries() {
    return numUnalignedQueries + numAlignedQueries;
  }

  public long getTotalAlignedQueryLength() {
    return this.totalAlignedQueryLength;
  }

  public double getTotalAlignedPenalty() {
    return this.totalAlignedPenalty;
  }

  public long getNumAlignedQueries() {
    return numAlignedQueries;
  }

  public Distribution getDistanceBetweenQueryComponents() {
    return distanceBetweenQueryComponents;
  }

  long numAlignedQueries = 0;
  long numUnalignedQueries = 0;
  double totalAlignedPenalty;
  long totalAlignedQueryLength;
  Distribution distanceBetweenQueryComponents = new Distribution();
}
