package mapper;

import java.io.BufferedReader;
import java.io.IOException;

// A FastaParser parses .fasta files
public class FastaParser implements SequenceProvider {
  public FastaParser(BufferedReader reader, String path) {
    this.reader = reader;
    this.path = path;
  }

  public SequenceBuilder getNextSequence() {
    try {
      return doGetNextSequence();
    } catch(IOException e) {
      throw new RuntimeException(e);
    }
  }

  public boolean get_allReadsContainQualityInformation() {
    return false;
  }

  private SequenceBuilder doGetNextSequence() throws IOException {
    // find sequence start marker
    if (!hasReadASequence) {
      int first = reader.read();
      if (first < 0)
        return null;
      char firstChar = (char)first;
      if (firstChar != '>') {
        throw new IllegalArgumentException("Not in .fasta format: expected '>', not '" + firstChar + "'");
      }
      hasReadASequence = true;
    }

    // read sequence name
    String nameLine = reader.readLine();
    if (nameLine == null)
      return null;
    int spaceIndex = nameLine.indexOf(' ');
    if (spaceIndex > 0)
      nameLine = nameLine.substring(0, spaceIndex);
    SequenceBuilder builder = new SequenceBuilder();
    builder.setName(nameLine);
    builder.setPath(path);

    // read sequence content until next sequence start marker
    while (true) {
      int firstInLine = reader.read();
      if (firstInLine < 0) {
        break;
      }
      char firstCharInLine = (char)firstInLine;
      if (firstCharInLine == '>') {
        break;
      }
      builder.add(firstCharInLine);
      String line = this.reader.readLine();
      builder.add(line);
    }
    return builder;
  }

  public int getNumErrors() {
    // Failures to read a .fasta are considered fatal so if we get here there aren't any errors
    return 0;
  }

  @Override
  public String toString() {
    return this.path;
  }

  BufferedReader reader;
  String path;
  boolean hasReadASequence = false;
}
