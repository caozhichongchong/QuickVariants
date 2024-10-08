package mapper;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

// A SamReader reads .sam files and parses individual lines but doesn't associate related lines
public class SamReader implements SequenceProvider {
  public SamReader(BufferedReader reader, String path) {
    this.reader = reader;
    this.path = path;
  }

  public SequenceBuilder getNextSequence() {
    while (true) {
      String line = null;
      try {
        line = this.getNextNonCommentLine();
      } catch (IOException e) {
        throw new RuntimeException("Error reading " + path + ",", e);
      }
      if (line == null) {
        this.done();
        return null;
      }
      SequenceBuilder sequenceBuilder;
      try {
        sequenceBuilder = this.parseLine(line);
      } catch (Exception e) {
        throw new RuntimeException("Failed to parse sam line '" + line + "'", e);
      }
      if (sequenceBuilder != null)
        return sequenceBuilder;
    }
  }

  public boolean isReordered() {
    return false;
  }

  private void done() {
    int numWarnings = numSkippedSupplementalAlignments + numAlignmentsMissingQueryText + numAlignmentsWithoutIndelInformation;
    if (numWarnings > 0) {
      System.out.println("" + numWarnings + " warnings reading " + path);
    }
  }

  private String getNextNonCommentLine() throws IOException {
    while (true) {
      String line = this.reader.readLine();
      if (line == null)
        return line;
      if (line.startsWith("@"))
        continue;
      return line;
    }
  }

  private SequenceBuilder parseLine(String line) {
    SequenceBuilder sequenceBuilder = new SequenceBuilder();
    sequenceBuilder.setPath(this.path);

    String[] fields = line.split("\t");
    String name = fields[0];
    sequenceBuilder.setName(name);

    String samFlagString = fields[1];
    int samFlags = Integer.parseInt(samFlagString);
    boolean referenceReversed = (samFlags & 16) != 0;
    boolean isSupplementaryAlignment = (samFlags & 2048) != 0;
    // skipping supplementary alignments for now
    if (isSupplementaryAlignment) {
      if (numSkippedSupplementalAlignments < 1)
        System.out.println("Warning: skipping supplemental alignments, including " + line);
      numSkippedSupplementalAlignments++;
      return null;
    }

    String referenceContigName = fields[2];
    if ("*".equals(referenceContigName)) {
      // This indicates an unalighed query
      // We're not interested in them at the moment
      return null;
    }
    int spaceIndex = referenceContigName.indexOf(' ');
    if (spaceIndex != -1)
      referenceContigName = referenceContigName.substring(0, spaceIndex);
    int startPosition = Integer.parseInt(fields[3]);
    String cigarString = fields[5];
    if ("*".equals(cigarString)) {
      // Alignment doesn't tell where the indels are
      if (numAlignmentsWithoutIndelInformation < 1)
        System.out.println("Warning: skipping alignment where CIGAR string (indel information) is '" + cigarString + "', including " + line);
      numAlignmentsWithoutIndelInformation++;
      return null;
    }

    String mateContigName = fields[6];
    boolean hasMate = !mateContigName.equals("*");
    boolean hasUnalignedMate = (samFlags & 8) != 0;
    boolean expectsMateAlignment = hasMate && !hasUnalignedMate;

    sequenceBuilder.asAlignment(referenceContigName, startPosition, cigarString, referenceReversed, expectsMateAlignment);

    String queryText = fields[9];
    if ("*".equals(queryText)) {
      if (numAlignmentsMissingQueryText < 1)
        System.out.println("Warning: skipping alignments having query text '" + queryText + "', including " + line);
      numAlignmentsMissingQueryText++;
      return null;
    }

    sequenceBuilder.add(queryText);

    return sequenceBuilder;
  }

  @Override
  public String toString() {
    return this.path;
  }

  public int getNumErrors() {
    return 0;
  }

  public boolean get_allReadsContainQualityInformation() {
    return false;
  }

  BufferedReader reader;
  String path;

  int numSkippedSupplementalAlignments;
  int numAlignmentsMissingQueryText;
  int numAlignmentsWithoutIndelInformation;
  int numErrors;
}
