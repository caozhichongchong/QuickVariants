package mapper;

import java.util.ArrayList;
import java.util.List;

// A SamAlignment represents a specific way in which a specific query can align to the reference genome
// It can model a single sequence or can model Illumina-style paired-end reads
// If there are multiple alternate ways that the query could align to the reference genome, each should be represented by its own SamAlignment
public class SamAlignment {
  public SamAlignment(SamRecord samRecord) {
    this.samRecords = new ArrayList<SamRecord>(1);
    this.samRecords.add((SamRecord)samRecord);
  }

  public SamAlignment(List<Sequence> sequences) {
    this.samRecords = new ArrayList<SamRecord>(sequences.size());
    for (Sequence sequence: sequences) {
      this.samRecords.add((SamRecord)sequence);
    }
  }

  public List<SamRecord> getSamRecords() {
    return this.samRecords;
  }

  public int getNumRecords() {
    return this.samRecords.size();
  }

  public long getId() {
    return this.samRecords.get(0).getId();
  }

  public int getLength() {
    int total = 0;
    for (Sequence sequence: this.samRecords) {
      total += sequence.getLength();
    }
    return total;
  }

  public String format() {
    int totalSize = 0;
    for (Sequence sequence : this.samRecords) {
       totalSize += sequence.getLength();
    }
    if (totalSize > 1000) {
      return "[" + this.samRecords.size() + " sequences totalling " + totalSize + " base pairs]";
    }

    StringBuilder builder = new StringBuilder();
    for (int i = 0; i < this.samRecords.size(); i++) {
      Sequence sequence = this.samRecords.get(i);
      builder.append(sequence.getText());
      if (i < this.samRecords.size() - 1) {
        builder.append(" / ");
      }
    }
    return builder.toString();
  }

  public boolean sameSequenceNames(SamAlignment other) {
    if (other.samRecords.size() != this.samRecords.size()) {
      return false;
    }
    for (int i = 0; i < this.samRecords.size(); i++) {
      if (!this.samRecords.get(i).getName().equals(other.samRecords.get(i).getName())) {
        return false;
      }
    }
    return true;
  }

  @Override
  public String toString() {
    StringBuilder builder = new StringBuilder();
    for (Sequence sequence : this.samRecords) {
      builder.append(sequence.getText() + " ");
    }
    return builder.toString();
  }

  private List<SamRecord> samRecords;
}
