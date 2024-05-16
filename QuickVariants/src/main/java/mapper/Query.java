package mapper;

import java.util.ArrayList;
import java.util.List;

// A Query is a question that can be asked about where some sequences align
// It can model a single sequence or can model Illumina-style paired-end reads
public class Query {
  public Query(Sequence sequence) {
    this.sequences = new ArrayList<Sequence>(1);
    this.sequences.add(sequence);
    this.spacingDeviationPerUnitPenalty = 1;
  }

  public Query(List<Sequence> sequences, double expectedInnerDistance, double spacingDeviationPerUnitPenalty) {
    this.sequences = sequences;
    this.maxOffset = maxOffset;
    this.expectedInnerDistance = expectedInnerDistance;
    this.spacingDeviationPerUnitPenalty = spacingDeviationPerUnitPenalty;
  }

  public List<Sequence> getSequences() {
    return this.sequences;
  }

  public int getNumSequences() {
    return this.sequences.size();
  }

  public double getExpectedInnerDistance() {
    return this.expectedInnerDistance;
  }

  public double getSpacingDeviationPerUnitPenalty() {
    return this.spacingDeviationPerUnitPenalty;
  }

  public long getId() {
    return this.sequences.get(0).getId();
  }

  public int getLength() {
    int total = 0;
    for (Sequence sequence: this.sequences) {
      total += sequence.getLength();
    }
    return total;
  }

  public String format() {
    int totalSize = 0;
    for (Sequence sequence : this.sequences) {
       totalSize += sequence.getLength();
    }
    if (totalSize > 1000) {
      return "[" + this.sequences.size() + " sequences totalling " + totalSize + " base pairs]";
    }

    StringBuilder builder = new StringBuilder();
    for (int i = 0; i < this.sequences.size(); i++) {
      Sequence sequence = this.sequences.get(i);
      builder.append(sequence.getText());
      if (i < this.sequences.size() - 1) {
        builder.append(" / ");
      }
    }
    return builder.toString();
  }

  public boolean sameSequenceNames(Query other) {
    if (other.sequences.size() != this.sequences.size()) {
      return false;
    }
    for (int i = 0; i < this.sequences.size(); i++) {
      if (!this.sequences.get(i).getName().equals(other.sequences.get(i).getName())) {
        return false;
      }
    }
    return true;
  }

  @Override
  public String toString() {
    StringBuilder builder = new StringBuilder();
    for (Sequence sequence : this.sequences) {
      builder.append(sequence.getText() + " ");
    }
    return builder.toString();
  }

  private List<Sequence> sequences;
  private int maxOffset; // max distance betweeen where the sequences can be aligned
  private double expectedInnerDistance;
  private double spacingDeviationPerUnitPenalty;
}
