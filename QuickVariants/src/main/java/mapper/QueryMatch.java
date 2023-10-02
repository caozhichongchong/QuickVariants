package mapper;

import java.util.List;

public class QueryMatch {
  public QueryMatch(List<SequenceMatch> components, int numDistinctHashBlockMismatches) {
    this.components = components;
    this.numDistinctHashBlockMismatches = numDistinctHashBlockMismatches;
  }
  public List<SequenceMatch> getComponents() {
    return this.components;
  }
  public int getNumSequences() {
    return this.components.size();
  }
  public int getNumDistinctHashBlockMismatches() {
    return this.numDistinctHashBlockMismatches;
  }

  public int getQueryTotalLength() {
    int total = 0;
    for (SequenceMatch match: this.components) {
      total += match.getSequenceA().getLength();
    }
    return total;
  }

  public int getTotalDistanceAcross() {
    SequenceMatch last = this.components.get(this.components.size() - 1);
    SequenceMatch first = this.components.get(0);
    if (this.getReversed())
      return first.getEndIndexB() - last.getStartIndexB();
    else
      return last.getEndIndexB() - first.getStartIndexB();
  }

  // Returns the total inner distance between subsequent pairs of components
  public int getTotalDistanceBetweenComponents() {
    int totalDistance = 0;
    SequenceMatch previousComponent = this.components.get(0);
    for (int i = 1; i < this.components.size(); i++) {
      SequenceMatch currentComponent = this.components.get(i);
      totalDistance += getDistance(previousComponent, currentComponent);
      previousComponent = currentComponent;
    }
    return totalDistance;
  }

  // Returns the distance between the two blocks
  // Can return a negative number if they overlap
  private int getDistance(SequenceMatch a, SequenceMatch b) {
    if (a.getSequenceB() != b.getSequenceB())
      return Integer.MAX_VALUE;
    int difference;
    if (this.getReversed())
      difference = a.getStartIndexB() - b.getEndIndexB();
    else
      difference = b.getStartIndexB() - a.getEndIndexB();
    return difference;
  }

  private boolean getReversed() {
    return this.components.get(0).getReversed();
  }

  private List<SequenceMatch> components;
  private int numDistinctHashBlockMismatches;
  private boolean reversed;
}
