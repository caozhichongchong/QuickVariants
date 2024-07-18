package mapper;

import java.util.ArrayList;
import java.util.List;
import org.junit.Assert;
import org.junit.Test;

public class MatchDatabase_Test {
  public MatchDatabase_Test() {
  }

  @Test
  public void testQueryEndingWithMismatch() {
    String queryText = "AACCACGT";
    String refText =   "AACCACGA";

    Sequence a = new SequenceBuilder().setName("a").add(queryText).build();
    Sequence b = new SequenceBuilder().setName("b").add(refText).build();
    SequenceAlignment sequenceAlignment = new SequenceAlignment(new AlignedBlock(a, b, 0, 0, queryText.length(), refText.length()), false);
    QueryAlignment alignment = new QueryAlignment(sequenceAlignment);
    MatchDatabase database = new MatchDatabase(0);
    List<QueryAlignment> alignmentList = new ArrayList<QueryAlignment>();
    alignmentList.add(alignment);
    List<List<QueryAlignment>> alignmentListList = new ArrayList<List<QueryAlignment>>();
    alignmentListList.add(alignmentList);
    database.addAlignments(alignmentListList);
    Alignments alignments = database.groupByPosition().get(b);
    for (int i = 0; i < refText.length(); i++) {
      AlignmentPosition position = alignments.getPosition(i);
      float count = position.getCount();
      if (count != 1) {
        fail("number of bases aligned at position " + i + " in reference equals " + count + " rather than 1");
      }
    }
  }

  @Test
  public void testNoncontiguousAlignment() { // which can happen if there is an 'N' in the sam cigar field
    String queryPrefix = "ACG";
    String refPrefix   = "AAG";

    String refMiddle = "TTTT";

    String querySuffix = "CCT";
    String refSuffix   = "CCT";

    String queryText = queryPrefix +             querySuffix;
    String refText   = refPrefix   + refMiddle + refSuffix;

    Sequence a = new SequenceBuilder().setName("a").add(queryText).build();
    Sequence b = new SequenceBuilder().setName("b").add(refText).build();

    List<AlignedBlock> alignedBlocks = new ArrayList<AlignedBlock>();
    alignedBlocks.add(new AlignedBlock(a, b, 0, 0, queryPrefix.length(), refPrefix.length()));
    alignedBlocks.add(new AlignedBlock(a, b, queryPrefix.length(), refPrefix.length() + refMiddle.length(), querySuffix.length(), refSuffix.length()));
    SequenceAlignment sequenceAlignment = new SequenceAlignment(alignedBlocks, false);
    QueryAlignment alignment = new QueryAlignment(sequenceAlignment);
    MatchDatabase database = new MatchDatabase(0);
    List<QueryAlignment> alignmentList = new ArrayList<QueryAlignment>();
    alignmentList.add(alignment);
    List<List<QueryAlignment>> alignmentListList = new ArrayList<List<QueryAlignment>>();
    alignmentListList.add(alignmentList);
    database.addAlignments(alignmentListList);
    Alignments alignments = database.groupByPosition().get(b);
    for (int i = 0; i < refText.length(); i++) {
      AlignmentPosition position = alignments.getPosition(i);
      float count = position.getCount();
      int expectedCount;
      if (i < refPrefix.length() || i >= refPrefix.length() + refMiddle.length())
        expectedCount = 1;
      else
        expectedCount = 0;

      if (count != expectedCount) {
        fail("At position " + i + " expected support of " + expectedCount + " but got " + count);
      }
    }
  }

  @Test
  public void testNoncontiguousPairedEndAlignment() { // which can happen if there is an 'N' in the sam cigar field
    String forwardQueryText = "AGCG" +"CCTT";
    String reverseQueryText =   "CGAT" +"TT";
    String referenceText    = "AACGTTTTCCTT";

    Sequence forward = new SequenceBuilder().setName("forward").add(forwardQueryText).build();
    Sequence reverse = new SequenceBuilder().setName("reverse").add(reverseQueryText).build();

    Sequence reference = new SequenceBuilder().setName("ref").add(referenceText).build();

    List<AlignedBlock> forwardBlocks = new ArrayList<AlignedBlock>();
    forwardBlocks.add(new AlignedBlock(forward, reference, 0, 0, 4, 4));
    forwardBlocks.add(new AlignedBlock(forward, reference, 4, 8, 4, 4));
    SequenceAlignment forwardAlignment = new SequenceAlignment(forwardBlocks, false);

    List<AlignedBlock> reverseBlocks = new ArrayList<AlignedBlock>();
    reverseBlocks.add(new AlignedBlock(forward, reference, 2, 2, 4, 4));
    reverseBlocks.add(new AlignedBlock(forward, reference, 6, 10, 2, 2));
    SequenceAlignment reverseAlignment = new SequenceAlignment(reverseBlocks, false);

    List<SequenceAlignment> sequenceAlignments = new ArrayList<SequenceAlignment>();
    sequenceAlignments.add(forwardAlignment);
    sequenceAlignments.add(reverseAlignment);

    QueryAlignment alignment = new QueryAlignment(sequenceAlignments, 0, 0, 0);

    MatchDatabase database = new MatchDatabase(0);
    List<QueryAlignment> alignmentList = new ArrayList<QueryAlignment>();
    alignmentList.add(alignment);
    List<List<QueryAlignment>> alignmentListList = new ArrayList<List<QueryAlignment>>();
    alignmentListList.add(alignmentList);
    database.addAlignments(alignmentListList);
    Alignments alignments = database.groupByPosition().get(reference);
    for (int i = 0; i < referenceText.length(); i++) {
      AlignmentPosition position = alignments.getPosition(i);
      float count = position.getCount();
      int expectedCount;
      // When paired-end alignments both align to a certain position, the support for that position should be 1
      if (i < 6 || i >= 8)
        expectedCount = 1;
      else
        expectedCount = 0;

      if (count != expectedCount) {
        fail("At position " + i + " expected support of " + expectedCount + " but got " + count);
      }
    }
  }

  private void fail(String message) {
    Assert.fail(message);
  }
}
