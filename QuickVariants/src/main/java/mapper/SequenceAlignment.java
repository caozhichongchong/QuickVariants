package mapper;

import java.util.ArrayList;
import java.util.List;

// A SequenceAlignment says that two sequences resemble each other
// A SequenceAlignment models insertions and deletions
public class SequenceAlignment {
  public SequenceAlignment(AlignedBlock block, boolean referenceReversed) {
    this.sections = new ArrayList<AlignedBlock>(1);
    this.sections.add(block);
    this.referenceReversed = referenceReversed;
  }

  public SequenceAlignment(List<AlignedBlock> sections, boolean referenceReversed) {
    this.sections = sections;
    this.referenceReversed = referenceReversed;
  }

  public void check() {
    for (int i = 0; i < this.sections.size() - 1; i++) {
      AlignedBlock prev = this.sections.get(i);
      AlignedBlock next = this.sections.get(i + 1);
      if ((next.getLengthA() > 0) != (next.getLengthB() > 0)) {
        // `next` is an indel
        if (next.getStartIndexA() != prev.getEndIndexA() && next.getEndIndexA() != prev.getStartIndexA()) {
          fail("next start a = " + next.getStartIndexA() + ", prev end a = " + prev.getEndIndexA() + ", next end a = " + next.getEndIndexA() + " prev start a = " + prev.getStartIndexA());
        }
        if (next.getStartIndexB() != prev.getEndIndexB() && next.getEndIndexB() != prev.getStartIndexB()) {
          fail("next start b = " + next.getStartIndexB() + ", prev end b = " + prev.getEndIndexB() + " next end b = " + next.getEndIndexB() + " prev start b = " + prev.getStartIndexB());
        }
      }
    }
    for (int i = 0; i < this.sections.size(); i++) {
      AlignedBlock block = this.sections.get(i);
      if (block.getStartIndexA() < 0) {
        fail("block start index " + block.getStartIndexA() + " < 0");
      }
      if (block.getEndIndexA() > block.getSequenceA().getLength()) {
        fail("block end index " + block.getEndIndexA() + " > sequence length " + block.getSequenceA().getLength() + " with text " + block.getSequenceA().getText());
      }
    }
  }

  private void fail(String message) {
    String n = null;
    System.out.println(message);
    System.out.println(n.toString());
  }

  public List<AlignedBlock> getSections() {
    return this.sections;
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
    for (AlignedBlock block: this.sections) {
      block.putSequenceB(sequence);
    }
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
  private boolean referenceReversed;
}
