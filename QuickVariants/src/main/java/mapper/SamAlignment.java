package mapper;

import java.util.List;
import java.util.ArrayList;

// A SamAlignment is an alignment from a .sam file
public class SamAlignment extends Sequence {
  public SamAlignment(String name, String packedContents, int length, String path) {
    super(name, packedContents, length, path);
  }

  public void setCigarString(String cigarString) {
    this.cigarComponents = this.parseCigarComponents(cigarString);
    this.cigarString = cigarString;
  }

  public SequenceAlignment toSequenceAlignment(SequenceDatabase sequenceDatabase) {
    List<AlignedBlock> sections = new ArrayList<AlignedBlock>(this.cigarComponents.size());
    int queryStartIndex = 0;
    int referenceStartIndex = this.referencePosition - 1;
    Sequence query = this;
    Sequence reference = sequenceDatabase.getSequence(this.referenceName);
    if (reference == null) {
      throw new IllegalArgumentException("Reference contig named '" + this.referenceName + "' not found");
    }
    for (String component: cigarComponents) {
      if (component.equals("=")) {
        component = "" + this.getLength() + "M";
      }
      if (component.equals("*"))
        return null;

      AlignedBlock alignedBlock;
      char alignmentTypeChar = component.charAt(component.length() - 1);
      String alignmentLengthString = component.substring(0, component.length() - 1);
      int alignmentLength;
      try {
        alignmentLength = Integer.parseInt(alignmentLengthString);
      } catch (NumberFormatException e) {
        throw new IllegalArgumentException("Could not parse cigar component '" + component + "' from cigar string '" + this.cigarString + "'", e);
      }
      if (alignmentTypeChar == 'M') {
        alignedBlock = new AlignedBlock(query, reference, queryStartIndex, referenceStartIndex, alignmentLength, alignmentLength);
      } else {
        if (alignmentTypeChar == 'S' || (alignmentTypeChar == 'I' && queryStartIndex <= 0)) {
          queryStartIndex += alignmentLength;
          continue;
        }
        if (alignmentTypeChar == 'I') {
          alignedBlock = new AlignedBlock(query, reference, queryStartIndex, referenceStartIndex, alignmentLength, 0);
        } else {
          if (alignmentTypeChar == 'D') {
            alignedBlock = new AlignedBlock(query, reference, queryStartIndex, referenceStartIndex, 0, alignmentLength);
          } else {
            if (alignmentTypeChar == 'H') {
              continue;
            }
            throw new IllegalArgumentException("Unrecognized cigar character '" + alignmentTypeChar + "' from cigar component '" + component + "' from full cigar string '" + cigarString + "'");
          }
        }
      }
      queryStartIndex += alignedBlock.getLengthA();
      referenceStartIndex += alignedBlock.getLengthB();
      sections.add(alignedBlock);
    }
    SequenceAlignment result = new SequenceAlignment(sections, this.referenceReversed);
    result.weight = this.weight;
    return result;
  }

  private List<String> parseCigarComponents(String cigarString) {
    List<String> parsed = new ArrayList<String>();

    int prev = -1;
    for (int i = 0; i < cigarString.length(); i++) {
      char c = cigarString.charAt(i);
      if (c < '0' || c > '9' || i == cigarString.length() - 1) {
        String component = cigarString.substring(prev + 1, i + 1);
        parsed.add(component);
        //System.out.println("parsed cigar component '" + component + "'");
        prev = i;
      }
    }
    return parsed;
  }

  public String referenceName;
  public int referencePosition;
  public boolean referenceReversed;
  public double weight;
  private List<String> cigarComponents;
  private String cigarString;
}
