package mapper;

import java.util.ArrayList;
import java.util.List;

// A SamAlignment represents a specific way in which a specific query can align to the reference genome
// It can model a single sequence or can model Illumina-style paired-end reads
// If there are multiple alternate ways that the query could align to the reference genome, each should be represented by its own SamAlignment
public class SamAlignment {
  public SamAlignment(Sequence sequence) {
    this.sequences = new ArrayList<Sequence>(1);
    this.sequences.add(sequence);
  }

  public SamAlignment(List<Sequence> sequences) {
    this.sequences = sequences;
  }

  public List<Sequence> getSequences() {
    return this.sequences;
  }

  public int getNumSequences() {
    return this.sequences.size();
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

  public boolean sameSequenceNames(SamAlignment other) {
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
}
