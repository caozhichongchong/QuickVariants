package mapper;

// A PositionDescriptor represents a position on a Sequence but doesn't have a reference to the sequence, only the name
public class PositionDescriptor {
  public PositionDescriptor(String sequenceName, int startIndex, boolean reverseComplemented) {
    this.sequenceName = sequenceName;
    this.startIndex = startIndex;
    this.reverseComplemented = reverseComplemented;
  }

  public String getSequenceName() {
    return this.sequenceName;
  }

  public int getStartIndex() {
    return this.startIndex;
  }

  public boolean isReverseComplemented() {
    return this.reverseComplemented;
  }

  public boolean equals(PositionDescriptor other) {
    if (!this.sequenceName.equals(other.sequenceName))
      return false;
    if (this.startIndex != other.startIndex)
      return false;
    if (this.reverseComplemented != other.reverseComplemented)
      return false;
    return true;
  }

  @Override
  public String toString() {
    String result = this.sequenceName;
    if (this.reverseComplemented)
      result += "-rev";
    result += "[" + startIndex + "]";
    return result;
  }

  private String sequenceName;
  private int startIndex;
  private boolean reverseComplemented;
}
