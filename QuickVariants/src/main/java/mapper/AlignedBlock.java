package mapper;

// An AlignedBlock says that two sections of two sequences align with each other
// An AlignedBlock is like a Match but it refers to a subsection of the Sequence
public class AlignedBlock {
  public AlignedBlock(Sequence sequenceA, Sequence sequenceB, int aStartIndex, int bStartIndex, int aLength, int bLength) {
    this.sequenceA = sequenceA;
    this.sequenceB = sequenceB;
    this.sequenceBHistory = this.sequenceB;
    this.aStartIndex = aStartIndex;
    this.bStartIndex = bStartIndex;
    this.aLength = aLength;
    this.bLength = bLength;
    if (aLength != bLength) {
      if (aLength != 0 && bLength != 0) {
        throw new IllegalArgumentException("Attempted to create AlignedBlock with unequal query and reference lengths when neither is 0: block query start index = " + aStartIndex + ", block query length = " + aLength + ", block ref start index " + bStartIndex + " block ref length " + bLength);
      }
    }
  }
  public AlignedBlock withSequenceA(Sequence newSequenceA) {
    return new AlignedBlock(newSequenceA, sequenceB, aStartIndex, bStartIndex, aLength, bLength);
  }
  public String getTextA() {
    return this.sequenceA.getRange(this.aStartIndex, this.aLength);
  }
  public String getTextB() {
    return this.sequenceB.getRange(this.bStartIndex, this.bLength);
  }
  public String getTextBHistory() {
    return this.sequenceBHistory.getRange(this.bStartIndex, this.bLength);
  }

  public int getStartIndexA() {
    return aStartIndex;
  }
  public int getStartIndexB() {
    return bStartIndex;
  }
  public int getEndIndexA() {
    return aStartIndex + aLength;
  }
  public int getEndIndexB() {
    return bStartIndex + bLength;
  }
  public int getOffset() {
    return bStartIndex - aStartIndex;
  }
  public Sequence getSubsequenceA() {
    return this.sequenceA.getSubsequence(this.aStartIndex, this.aLength);
  }
  public Sequence getSequenceA() {
    return sequenceA;
  }

  public Sequence getSequenceB() {
    return sequenceB;
  }

  public Sequence getSequenceBHistory() {
    return sequenceBHistory;
  }

  public int getLengthA() {
    return aLength;
  }

  public int getLengthB() {
    return bLength;
  }

  public byte getFirstEncodedCharA() {
    return this.sequenceA.encodedCharAt(this.getStartIndexA());
  }
  public byte getFirstEncodedCharB() {
    return this.sequenceB.encodedCharAt(this.getStartIndexB());
  }
  public byte getLastEncodedCharA() {
    return this.sequenceA.encodedCharAt(this.getEndIndexA() - 1);
  }
  public byte getLastEncodedCharB() {
    return this.sequenceB.encodedCharAt(this.getEndIndexB() - 1);
  }

  public boolean equals(AlignedBlock other) {
    if (bStartIndex != other.bStartIndex)
      return false;
    if (aStartIndex != other.aStartIndex)
      return false;
    if (bLength != other.bLength)
      return false;
    if (aLength != other.aLength)
      return false;
    if (sequenceB != other.sequenceB)
      return false;
    if (sequenceA != other.sequenceA)
      return false;
    return true;
  }

  public int getIndelLength() {
    if (aLength == bLength)
      return 0;
    return Math.max(aLength, bLength);
  }

  // Modifies this alignment to refer to a new sequence
  // This can be useful if this alignment was computed in one way and applied to another sequence, for example, computed via an ancestor and applied to a child
  public void putSequenceB(Sequence newSequenceB) {
    this.sequenceB = newSequenceB;
  }

  // tells whether the indel status of this block (insertion, deletion, or neither) is the same as the other block
  public boolean sameIndelType(AlignedBlock other) {
    if ((this.aLength > 0) != (other.aLength > 0))
      return false;
    if ((this.bLength > 0) != (other.bLength > 0))
      return false;
    return true;
  }

  public boolean hasAmbiguousBasepairs() {
    for (int i = this.aStartIndex; i < this.aStartIndex + this.aLength; i++) {
      if (Basepairs.isAmbiguous(sequenceA.encodedCharAt(i)))
        return true;
    }
    for (int i = this.bStartIndex; i < this.bStartIndex + this.bLength; i++) {
      if (Basepairs.isAmbiguous(sequenceB.encodedCharAt(i)))
        return true;
    }
    return false;
  }

  public Sequence sequenceA;
  public Sequence sequenceB;
  public Sequence sequenceBHistory;
  public int aStartIndex;
  public int bStartIndex;
  public int aLength;
  public int bLength;
}
