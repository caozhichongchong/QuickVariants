package mapper;

class VcfFormatRequest {
  public VcfFormatRequest(Sequence sequence, int startIndex, int length, FilteredAlignments alignments, int jobIndex) {
    this.sequence = sequence;
    this.startIndex = startIndex;
    this.length = length;
    this.alignments = alignments;
    this.jobIndex = jobIndex;
  }

  public Sequence sequence;
  public int startIndex;
  public int length;
  public FilteredAlignments alignments;
  public int jobIndex;

  @Override
  public String toString() {
    return "VcfFormatRequest " + jobIndex + " on " + sequence.getName() + " at " + startIndex + " length " + length;
  }
}
