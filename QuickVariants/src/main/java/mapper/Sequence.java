package mapper;

import java.util.Arrays;
import java.util.List;

// A Sequence is a list of genomic base pairs: A,C,G,T
public class Sequence {
  public Sequence(String name, String packedContents, int length, String path) {
    this.name = name;
    this.packedContents = packedContents;
    this.length = length;
    this.path = path;
  }

  // returns the name of this sequence
  public String getName() {
    return this.name;
  }

  // returns the name of the sequence that this sequence came from, or its own name if it didn't come from another sequence
  public String getSourceName() {
    return this.getName();
  }

  public String getPath() {
    return this.path;
  }

  // an internal unique identifier that is shorter than getName
  public long getId() {
    return this.identifier;
  }

  public void setId(long identifier) {
    this.identifier = identifier;
  }

  public String getRange(int startIndex, int count) {
    StringBuilder builder = new StringBuilder();
    int endIndex = startIndex + count;
    for (int i = startIndex; i < endIndex; i++) {
      builder.append(this.charAt(i));
    }
    return builder.toString();
  }

  public char charAt(int index) {
    byte basePair = encodedCharAt(index);
    return Basepairs.decode(basePair);
  }

  public byte encodedCharAt(int index) {
    int bitsPerBasePair = 4;
    int bitsPerCharacter = 16;
    int characterIndex = index * bitsPerBasePair / bitsPerCharacter;
    int offsetInCharacter = (index * bitsPerBasePair) % bitsPerCharacter;

    byte encoded = selectBits(characterIndex, offsetInCharacter, 4);
    return encoded;
  }

  private byte selectBits(int charIndex, int index, int count) {
    char character = this.packedContents.charAt(charIndex);
    int maxIndex = index + count;
    int shiftAmount = index;
    byte shifted = (byte)(character >> shiftAmount);
    byte inclusionBitmask = (byte)((1 << count) - 1);
    byte result = (byte)(shifted & inclusionBitmask);
    return result;
  }

  public Sequence getSubsequence(int startIndex, int count) {
    if (startIndex == 0 && count == this.getLength()) {
      return this;
    }
    return new Subsequence(this, startIndex, count);
  }

  public String getText() {
    return this.getRange(0, this.getLength());
  }

  public int getLength() {
    return this.length;
  }

  public Sequence reverseComplement() {
    return new ReverseComplementSequence(this);
  }

  // returns the Sequence that this one was created as the reverseComplement of, if any
  public Sequence getComplementedFrom() {
    return null;
  }

  public int compareTo(Sequence other) {
    return Long.compare(this.identifier, other.identifier);
  }

  private String name;
  private String packedContents;
  private long identifier;
  private int length;
  private String path;
}
