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
    int numWarnings = numAlignmentsMissingQueryText + numAlignmentsWithoutIndelInformation + numUnparseableSupplementaryAlignments;
    if (numWarnings > 0) {
      System.out.println("" + numWarnings + " warnings reading individual lines in " + path + ":");
      if (numAlignmentsMissingQueryText > 0) {
        System.out.println(" " + numAlignmentsMissingQueryText + " alignments missing query text");
      }
      if (numAlignmentsWithoutIndelInformation > 0) {
        System.out.println(" " + numAlignmentsWithoutIndelInformation + " alignments without indel information");
      }
      if (numUnparseableSupplementaryAlignments > 0) {
        System.out.println(" " + numUnparseableSupplementaryAlignments + " unparseable supplementary alignments");
      }
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

    String referenceContigName = fields[2];
    if ("*".equals(referenceContigName)) {
      // This indicates an unaligned query
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
    if ("=".equals(mateContigName))
      mateContigName = referenceContigName;
    boolean hasMate = !mateContigName.equals("*");
    boolean hasUnalignedMate = (samFlags & 8) != 0;
    boolean expectsMateAlignment = hasMate && !hasUnalignedMate;
    boolean mateReversed = (samFlags & 32) != 0;

    List<PositionDescriptor> otherComponentPositions = new ArrayList<PositionDescriptor>();
    if (expectsMateAlignment) {
      int matePosition = Integer.parseInt(fields[7]);
      otherComponentPositions.add(new PositionDescriptor(mateContigName, matePosition, mateReversed));
    }
    
    int queryTextIndex = 9;
    String queryText = fields[queryTextIndex];
    if ("*".equals(queryText)) {
      if (numAlignmentsMissingQueryText < 1)
        System.out.println("Warning: skipping alignments having query text '" + queryText + "', including " + line);
      numAlignmentsMissingQueryText++;
      return null;
    }

    sequenceBuilder.add(queryText);

    for (int i = queryTextIndex + 1; i < fields.length; i++) {
      String field = fields[i];
      List<PositionDescriptor> fieldAsNextComponentStarts = tryParseSupplementaryPosition(field);
      if (fieldAsNextComponentStarts != null)
        otherComponentPositions.addAll(fieldAsNextComponentStarts);
    }

    sequenceBuilder.asAlignment(referenceContigName, startPosition, cigarString, referenceReversed, otherComponentPositions);

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

  private List<PositionDescriptor> tryParseSupplementaryPosition(String text) {
    String prefix = "SA:Z:";
    if (!text.startsWith(prefix)) {
      return null;
    }
    String fieldValue = text.substring(prefix.length());
    String[] entries = fieldValue.split(";");
    List<PositionDescriptor> results = new ArrayList<PositionDescriptor>();
    for (String entry: entries) {
      String[] fields = entry.split(",");
      if (fields.length != 6) {
        if (this.numUnparseableSupplementaryAlignments < 1) {
          System.out.println("Warning: ignoring supplementary alignment text with unexpected number of comma-separated fields (" + fields.length + "): '" + text + "'");
        }
        this.numUnparseableSupplementaryAlignments++;
        return null;
      }
      String contigName = fields[0];
      int position = Integer.parseInt(fields[1]);
      boolean reversed = fields[2].equals("-");
      results.add(new PositionDescriptor(contigName, position, reversed));
    }
    return results;
  }
  
  BufferedReader reader;
  String path;

  int numAlignmentsMissingQueryText;
  int numAlignmentsWithoutIndelInformation;
  int numUnparseableSupplementaryAlignments;
}
