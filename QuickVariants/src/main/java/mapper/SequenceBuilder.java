package mapper;

import java.util.ArrayList;

public class SequenceBuilder {
  public SequenceBuilder() {
    this.stringBuilder = new StringBuilder();
  }
  public SequenceBuilder setName(String name) {
    this.name = name;
    return this;
  }
  public String getName() {
    return this.name;
  }
  public void setPath(String path) {
    this.path = path;
  }

  private String compress(String text) {
    this.currentValue = 0;
    this.currentCount = 0;
    StringBuilder compressed = new StringBuilder();
    for (int i = 0; i < text.length(); i++) {
      char basePair = text.charAt(i);
      byte newValue = Basepairs.encode(basePair);

      this.currentValue += ((int)newValue << this.currentCount);
      this.currentCount += 4;
      if (this.currentCount >= 16) {
        compressed.append(this.emitChar());
      }
    }
    if (currentCount > 0) {
      compressed.append(this.emitChar());
    }
    return compressed.toString();
  }

  public SequenceBuilder add(char text) {
    stringBuilder.append(text);
    this.length++;
    return this;
  }

  public SequenceBuilder add(String text) {
    stringBuilder.append(text);
    this.length += text.length();
    return this;
  }

  public void setId(long identifier) {
    this.identifier = identifier;
  }

  public Sequence build() {
    String compressed = this.compress(this.stringBuilder.toString().toUpperCase());
    if (this.buildSam) {
      SamAlignment result = new SamAlignment(name, compressed, length, path);
      result.referenceName = this.referenceName;
      result.referencePosition = this.referencePosition;
      result.setCigarString(this.cigarString);
      result.referenceReversed = this.referenceReversed;
      result.weight = this.alignmentWeight;
      result.setId(this.identifier);
      return result;
    }
    if (this.buildRead) {
      ReadSequence result = new ReadSequence(name, compressed, length, path);
      result.nameSuffix = this.nameSuffix;
      result.qualityString = this.qualityString;
      result.commentString = this.commentString;
      result.setId(this.identifier);
      return result;
    } else {
      Sequence result = new Sequence(name, compressed, length, path);
      result.setId(this.identifier);
      return result;
    }
  }

  public SequenceBuilder asRead(String nameSuffix, String qualityString, String commentString) {
    this.buildRead = true;
    this.nameSuffix = nameSuffix;
    this.qualityString = qualityString;
    this.commentString = commentString;
    return this;
  }

  public SequenceBuilder asAlignment(String referenceName, int referencePosition, String cigarString, boolean referenceReversed, boolean expectsMate) {
    this.buildSam = true;
    this.referenceName = referenceName;
    this.referencePosition = referencePosition;
    this.cigarString = cigarString;
    this.referenceReversed = referenceReversed;
    this.expectsMate = expectsMate;
    return this;
  }

  public boolean getExpectsMate() {
    return this.expectsMate;
  }

  public SequenceBuilder withAlignmentMate(SequenceBuilder other) {
    this.mate = other;
    return this;
  }

  public boolean getHasAlignmentMate() {
    return this.mate != null;
  }

  public int getLength() {
    return this.length;
  }

  public SequenceBuilder setAlignmentWeight(double weight) {
    this.alignmentWeight = weight;
    return this;
  }

  private char emitChar() {
    char newValue = (char)(this.currentValue % 65536);
    this.currentCount -= 16;
    this.currentValue = this.currentValue >> 16;
    return newValue;
  }

  private String name;
  private String path;
  StringBuilder stringBuilder;
  int currentValue;
  int currentCount;
  int length;
  boolean buildRead = false;
  String nameSuffix;
  String qualityString;
  String commentString;
  long identifier;

  boolean buildSam;
  String referenceName;
  int referencePosition;
  String cigarString;
  boolean referenceReversed;
  SequenceBuilder mate;
  boolean expectsMate;
  double alignmentWeight = 1;
}
