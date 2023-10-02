package mapper;

import java.util.ArrayList;
import java.util.List;

// A QueryAlignment tells where a Query aligns
// Its like a SequenceAlignment but the query might consist of multiple sequences if it is two Illumina-style paired-end reads
public class QueryAlignment {

  public QueryAlignment(SequenceAlignment sequenceAlignment) {
    this.alignments = new ArrayList<SequenceAlignment>(1);
    this.alignments.add(sequenceAlignment);
  }

  public QueryAlignment(List<SequenceAlignment> sequenceAlignments, double spacingPenalty, double mirrorPenalty, double totalPenalty) {
    this.alignments = sequenceAlignments;
  }

  // list of alignments for each sequence
  public List<SequenceAlignment> getComponents() {
    return alignments;
  }

  public Sequence getSequenceB() {
    return this.alignments.get(0).getSection(0).getSequenceB();
  }

  public void putSequenceB(Sequence sequence) {
    for (SequenceAlignment component: this.alignments) {
      component.putSequenceB(sequence);
    }
  }

  public int getNumSequences() {
    return this.alignments.size();
  }

  public int getALength() {
    int total = 0;
    for (SequenceAlignment component: this.alignments) {
      total += component.getALength();
    }
    return total;
  }

  // Returns the total inner distance between subsequent pairs of components
  public int getTotalDistanceBetweenComponents() {
    int totalDistance = 0;
    SequenceAlignment previousComponent = this.alignments.get(0);
    for (int i = 1; i < this.alignments.size(); i++) {
      SequenceAlignment currentComponent = this.alignments.get(i);
      totalDistance += getDistance(previousComponent.getLastSection(), currentComponent.getFirstSection());
      previousComponent = currentComponent;
    }
    return totalDistance;
  }

  public int getNumAlignmentsCoveringIndexB(int referenceIndex) {
    if (this.minOverlap < 0)
      this.computeOverlap();
    if (referenceIndex < this.minOverlap)
      return 1;
    if (referenceIndex >= this.maxOverlap)
      return 1;
    return this.alignments.size();
  }

  private void computeOverlap() {
    for (SequenceAlignment alignment: this.alignments) {
      int min = alignment.getStartIndexB();
      int max = alignment.getEndIndexB();
      if (this.minOverlap < 0 || min >= this.minOverlap) {
        this.minOverlap = min;
      }
      if (this.maxOverlap < 0 || max <= this.maxOverlap) {
        this.maxOverlap = max;
      }
    }
  }

  // returns the distance between the two blocks
  private int getDistance(AlignedBlock a, AlignedBlock b) {
    if (a.getSequenceB() != b.getSequenceB())
      return Integer.MAX_VALUE;
    if (b.getSequenceB().getComplementedFrom() != null) {
      throw new IllegalArgumentException("QueryAlignment.getDistnace with a reversed: " + a.getSequenceB().getName());
    }
    return b.getStartIndexB() - a.getEndIndexB();
  }

  private List<SequenceAlignment> alignments;

  private int minOverlap = -1;
  private int maxOverlap = -1;

  // If multiple SequenceAlignments overlap, we multiply any penalties on non-overlapping portions of the alignments, so that we count them the same number of times (2) as overlapping portions of the alignments
}
