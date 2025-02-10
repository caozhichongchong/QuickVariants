package mapper;

public class VcfFormatterWorker extends Thread {
  public VcfFormatterWorker(boolean includeNonMutations, boolean showSupportRead) {
    this.includeNonMutations = includeNonMutations;
    this.showSupportRead = showSupportRead;
  }

  public void request(VcfFormatRequest formatRequest) {
    this.formatRequest = formatRequest;
  }

  public String getResults() {
    return results;
  }

  public int getNumReferencePositionsMatched() {
    return this.numReferencePositionsMatched;
  }

  @Override
  public void run() {
    this.results = this.format(this.formatRequest);
  }

  private String format(VcfFormatRequest formatRequest) {
    StringBuilder stringBuilder = new StringBuilder();
    boolean showSupportRead = this.showSupportRead;
    String sequenceName = formatRequest.sequence.getName();
    FilteredAlignments bases = formatRequest.alignments;
    int endIndex = formatRequest.startIndex + formatRequest.length;
    for (int i = formatRequest.startIndex; i < endIndex; i++) {
      AlignmentPosition frequencies = bases.getPosition(i);
      int displayIndex = i + 1;
      if (frequencies.getCount() > 0) {
        if (this.includeNonMutations || frequencies.hasAlternates())
          writePosition(sequenceName, displayIndex, frequencies, stringBuilder, showSupportRead);
        this.numReferencePositionsMatched++;
      }
      int insertionIndex = 0;
      while (true) {
        AlignmentPosition insertion = bases.getInsertion(i, insertionIndex);
        if (!insertion.hasAlternates())
          break;
        writePosition(sequenceName, -1 * displayIndex, insertion, stringBuilder, showSupportRead);
        insertionIndex++;
      }
    }
    return stringBuilder.toString();
  }

  private void writePosition(String sequenceName, int rowNumber, AlignmentPosition frequencies, StringBuilder stringBuilder, boolean showSupportRead) {
    stringBuilder.append(sequenceName);
    stringBuilder.append('\t');
    stringBuilder.append(Integer.toString(rowNumber));
    stringBuilder.append('\t');

    char reference = frequencies.getReference();
    stringBuilder.append(reference);
    stringBuilder.append('\t');
    char[] alternates = frequencies.getNonzeroAlternates(reference);
    appendJoined(stringBuilder, ',', alternates);
    stringBuilder.append('\t');
    stringBuilder.append(frequencies.formatCount());
    stringBuilder.append('\t');

    writeDepths(frequencies, reference, alternates, false, stringBuilder);
    stringBuilder.append('\t');
    writeDepths(frequencies, reference, alternates, true, stringBuilder);

    if (showSupportRead) {
      stringBuilder.append('\t');
      boolean isFirst = true;
      for (char alternate : alternates) {
        // we found a mutation; let's show one sample read that mapped to this position
        Sequence sampleQuery = frequencies.getSampleAlternateSequence(alternate);
        if (isFirst) {
          isFirst = false;
        } else {
          stringBuilder.append(",");
        }
        if (sampleQuery != null) {
          int index = frequencies.getSampleAlternateIndex(alternate);
          if (frequencies.isSampleAlternateDeletion(alternate)) {
            stringBuilder.append(sampleQuery.getRange(0, index));
            stringBuilder.append("[-]");
            stringBuilder.append(sampleQuery.getRange(index, sampleQuery.getLength() - index));
          } else {
            stringBuilder.append(sampleQuery.getRange(0, index));
            stringBuilder.append('[');
            stringBuilder.append(sampleQuery.charAt(index));
            stringBuilder.append(']');
            stringBuilder.append(sampleQuery.getRange(index + 1, sampleQuery.getLength() - 1 - index));
          }
        }
      }
    }
    stringBuilder.append('\n');
  }

  private void writeDepths(AlignmentPosition frequencies, char reference, char[] alternates, boolean isQueryEnd, StringBuilder stringBuilder) {
    stringBuilder.append(frequencies.getCounts(reference, isQueryEnd));
    for (char alternate : alternates) {
      stringBuilder.append(';');
      stringBuilder.append(frequencies.getCounts(alternate, isQueryEnd));
    }
  }

  private void appendJoined(StringBuilder builder, char separator, char[] components) {
    boolean first = true;
    for (int i = 0; i < components.length; i++) {
      char component = components[i];
      if (first) {
        first = false;
      } else {
        builder.append(separator);
      }
      builder.append(component);
    }
  }

  VcfFormatRequest formatRequest;
  String results;
  int numReferencePositionsMatched;
  boolean includeNonMutations;
  boolean showSupportRead;
}
