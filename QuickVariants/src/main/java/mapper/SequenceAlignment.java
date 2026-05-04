package mapper;

import java.util.ArrayList;
import java.util.List;

// A SequenceAlignment says that two sequences resemble each other
// A SequenceAlignment models insertions and deletions
public class SequenceAlignment {
  public SequenceAlignment(AlignedBlock block, boolean referenceReversed, double totalPenalty, double alignedPenalty) {
    this.sections = new ArrayList<AlignedBlock>(1);
    this.sections.add(block);
    this.referenceReversed = referenceReversed;
    this.referenceContiguous = true;
    this.penalty = totalPenalty;
    this.alignedPenalty = alignedPenalty;
  }

  public SequenceAlignment(List<AlignedBlock> sections, boolean referenceReversed, double totalPenalty, double alignedPenalty) {
    this.sections = sections;
    this.referenceReversed = referenceReversed;
    this.computeContiguous();
    this.penalty = totalPenalty;
    this.alignedPenalty = alignedPenalty;
  }

  private void computeContiguous() {
    this.referenceContiguous = true;
    for (int i = 0; i < sections.size() - 1; i++) {
      if (this.sections.get(i).getEndIndexB() != this.sections.get(i + 1).getStartIndexB()) {
        this.referenceContiguous = false;
        return;
      }
    }
  }

  private void fail(String message) {
    throw new IllegalArgumentException(message);
  }

  public List<AlignedBlock> getSections() {
    return this.sections;
  }

  public int getNumSections() {
    return this.sections.size();
  }

  public AlignedBlock getSection(int index) {
    return this.sections.get(index);
  }

  public AlignedBlock getFirstSection() {
    return this.sections.get(0);
  }

  public AlignedBlock getLastSection() {
    return this.sections.get(this.sections.size() - 1);
  }

  public int getStartOffset() {
    return this.sections.get(0).getOffset();
  }

  public int getStartIndexB() {
    return this.sections.get(0).getStartIndexB();
  }

  public int getEndIndexB() {
    return this.sections.get(this.sections.size() - 1).getEndIndexB();
  }

  public int getStartIndexA() {
    return this.sections.get(0).getStartIndexA();
  }

  public int getEndIndexA() {
    return this.sections.get(this.sections.size() - 1).getEndIndexA();
  }

  public int getLengthA() {
    return this.getEndIndexA() - this.getStartIndexA();
  }

  public int getLengthB() {
    return this.getEndIndexB() - this.getStartIndexB();
  }

  public double getPenalty() {
    return this.penalty;
  }

  public double getAlignedPenalty() {
    return this.alignedPenalty;
  }

  public Sequence getSequenceA() {
    return this.sections.get(0).getSequenceA();
  }

  public Sequence getSequenceB() {
    return this.sections.get(0).getSequenceB(); 
  }

  public Sequence getSequenceBHistory() {
    return this.sections.get(0).getSequenceBHistory();
  }

  public String getAlignedTextA() {
    String result = "";
    for (AlignedBlock block: sections) {
      if (block.aLength > 0) {
        result += block.getTextA();
      } else {
        for (int i = 0; i < block.bLength; i++) {
          result += "-";
        }
      }
    }
    return result;
  }

  public String getAlignedTextBHistory() {
    String result = "";
    for (AlignedBlock block: sections) {
      if (block.bLength > 0) {
        result += block.getTextBHistory();
      } else {
        for (int i = 0; i < block.aLength; i++) {
          result += "-";
        }
      }
    }
    return result; 
  }

  public String getAlignedTextB() {
    String result = "";
    for (AlignedBlock block: sections) {
      if (block.bLength > 0) {
        result += block.getTextB();
      } else {
        for (int i = 0; i < block.aLength; i++) {
          result += "-";
        }
      }
    }
    return result;
  }

  // the length of the section of SequenceB that is included in the alignment
  public int getALength() {
    int total = 0;
    for (AlignedBlock block: sections) {
      total += block.aLength;
    }
    return total;
  }

  public boolean isReferenceReversed() {
    return referenceReversed;
  }

  public String format() {
    return this.getAlignedTextA() + "\n" + this.getAlignedTextB();
  }

  public void putSequenceB(Sequence sequence) {
    if (sequence == null) {
      throw new IllegalArgumentException("putSequenceB(null) for SequenceAlignment:\n" + this.format());
    }
    for (AlignedBlock block: this.sections) {
      block.putSequenceB(sequence);
    }
  }

  // Whether the sections of the reference described by this alignment are contiguous
  // This can be true even if there are indels, but not if there is a sam cigar alignment 'N' character
  public boolean isReferenceContiguous() {
    return this.referenceContiguous;
  }

  public boolean coversIndexB(int index) {
    for (AlignedBlock block: this.sections) {
      if (block.getStartIndexB() <= index && block.getEndIndexB() > index) {
        return true;
      }
    }
    return false;
  }

  @Override
  public int hashCode() {
    return this.sections.get(0).getOffset();
  }

  @Override
  public boolean equals(Object otherObject) {
    SequenceAlignment other = (SequenceAlignment)otherObject;
    if (this.sections.size() != other.sections.size())
      return false;
    if (this.referenceReversed != other.referenceReversed)
      return false;
    for (int i = 0; i < this.sections.size(); i++) {
      AlignedBlock ourBlock = this.sections.get(i);
      AlignedBlock theirBlock = other.sections.get(i);
      if (!ourBlock.equals(theirBlock)) {
        return false;
      }
    }
    return true;
  }

  public double weight = 1;
  private List<AlignedBlock> sections;
  private double penalty;
  private double alignedPenalty;
  private boolean referenceReversed;
  private boolean referenceContiguous;
}
