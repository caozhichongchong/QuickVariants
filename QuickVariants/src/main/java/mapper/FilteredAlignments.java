package mapper;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

// A FilteredAlignments tells how some subset of Sequences align to a particular reference contig
public class FilteredAlignments {
  public FilteredAlignments(Alignments alignments, MutationDetectionParameters filter) {
    this.alignments = alignments;
    this.filter = filter;
  }

  public int size() {
    return this.alignments.size();
  }

  public AlignmentPosition getPosition(int referenceIndex) {
    AlignmentPosition result = this.alignments.getPosition(referenceIndex);
    //result = filterSNPs(result, referenceIndex);
    //result = filterDeletions(result, referenceIndex);
    return result;
  }

  public AlignmentPosition getInsertion(int referenceIndex, int insertionIndex) {
    AlignmentPosition result = this.alignments.getInsertion(referenceIndex, insertionIndex);
    //result = filterSNPs(result, referenceIndex);
    //result = filterInsertions(result, referenceIndex, insertionIndex);
    return result;
  }

  private AlignmentPosition filterSNPs(AlignmentPosition position, int referenceIndex) {
    float totalCount = position.getCount();
    char[] nonzeroAlternates = position.getNonzeroAlternates();
    for (char alternate: nonzeroAlternates) {
      if (alternate != '-') {
        float alternateCount = position.getAlternateCount(alternate);
        if (!filter.supportsSNP(alternateCount, totalCount)) {
          position.replaceAlternateWithReference(alternate);
        }
      }
    }
    return position;
  }

  private AlignmentPosition filterDeletions(AlignmentPosition position, int referenceIndex) {
    // TODO: implement this
    return position;
  }

  private AlignmentPosition filterInsertions(AlignmentPosition position, int referenceIndex, int insertionIndex) {
    // TODO: implement this
    return position;
  }

  Alignments alignments;
  MutationDetectionParameters filter;
}
