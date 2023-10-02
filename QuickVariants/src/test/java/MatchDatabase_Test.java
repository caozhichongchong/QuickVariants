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

  private void fail(String message) {
    Assert.fail(message);
  }
}
