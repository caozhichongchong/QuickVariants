package mapper;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import org.junit.Assert;
import org.junit.Test;

public class SamAlignmentBuilder_Test {

  @Test
  public void testEmpty() {
    SamAlignment_Builder builder = new SamAlignment_Builder();
    if (builder.isComplete()) {
      fail("empty SamAlignment_Builder.isComplete() should be true");
    }
  }

  @Test
  public void testUnpairedRead() {
    SamAlignment_Builder builder = new SamAlignment_Builder();
    SequenceBuilder read1 = new SequenceBuilder().setName("read1");
    if (!builder.accepts(read1)) {
      fail("builder should accept new unpaired read");
    }
    builder.add(read1);
    if (!builder.isComplete()) {
      fail("added unpaired read and builder is not complete");
    }
    SequenceBuilder read2 = new SequenceBuilder().setName("read2");
    if (builder.accepts(read2)) {
      fail("added unpaired read and builder still accepts more reads");
    }
  }

  @Test
  public void testPairedRead() {
    SamAlignment_Builder builder = new SamAlignment_Builder();
    List<PositionDescriptor> position2s = new ArrayList<PositionDescriptor>();
    position2s.add(new PositionDescriptor("ref", 0, true));
    List<PositionDescriptor> position1s = new ArrayList<PositionDescriptor>();
    position1s.add(new PositionDescriptor("ref", 0, false));

    SequenceBuilder mate1 = new SequenceBuilder().setName("query").asAlignment("ref", 0, null, false, position2s);
    SequenceBuilder mate2 = new SequenceBuilder().setName("query").asAlignment("ref", 0, null, true, position1s);

    SequenceBuilder wrongQuery = new SequenceBuilder().setName("otherQuery").asAlignment("ref", 0, null, true, position1s);
    
    if (!builder.accepts(mate1)) {
      fail("builder should accept mate1");
    }
    builder.add(mate1);
    if (builder.isComplete()) {
      fail("builder should not be complete after adding one mate");
    }
    if (builder.accepts(wrongQuery)) {
      fail("builder should not accept mate with different name");
    }
    if (!builder.accepts(mate2)) {
      fail("builder should accept mate with same name");
    }
    builder.add(mate2);

    if (!builder.isComplete()) {
      fail("builder should be complete after adding both mates");
    }
  }

  @Test
  public void testUnpairedSupplementaryAlignments() {
    SamAlignment_Builder builder = new SamAlignment_Builder();
    List<PositionDescriptor> position2s = new ArrayList<PositionDescriptor>();
    position2s.add(new PositionDescriptor("ref", 100, false));
    List<PositionDescriptor> position1s = new ArrayList<PositionDescriptor>();
    position1s.add(new PositionDescriptor("ref", 0, false));

    SequenceBuilder component1 = new SequenceBuilder().setName("query").asAlignment("ref", 0, null, false, position2s);
    SequenceBuilder component2 = new SequenceBuilder().setName("query").asAlignment("ref", 100, null, false, position1s);

    if (!builder.accepts(component1)) {
      fail("builder should accept component1");
    }
    builder.add(component1);
    if (builder.isComplete()) {
      fail("builder should not be complete after adding first component");
    }
    if (!builder.accepts(component2)) {
      fail("builder should accept component2");
    }
    builder.add(component2);
    if (!builder.isComplete()) {
      fail("builder should be complete after adding last component");
    }
  }

  @Test
  public void testMultipleSupplementaryAlignments() {
    SamAlignment_Builder builder = new SamAlignment_Builder();
    List<PositionDescriptor> position1Others = new ArrayList<PositionDescriptor>();
    position1Others.add(new PositionDescriptor("ref", 50, false));
    position1Others.add(new PositionDescriptor("ref", 100, false));
    List<PositionDescriptor> position2Others = new ArrayList<PositionDescriptor>();
    position2Others.add(new PositionDescriptor("ref", 0, false));
    position2Others.add(new PositionDescriptor("ref", 100, false));
    List<PositionDescriptor> position3Others = new ArrayList<PositionDescriptor>();
    position3Others.add(new PositionDescriptor("ref", 0, false));
    position3Others.add(new PositionDescriptor("ref", 50, false));

    SequenceBuilder component1 = new SequenceBuilder().setName("query").asAlignment("ref", 0, null, false, position1Others);
    SequenceBuilder component2 = new SequenceBuilder().setName("query").asAlignment("ref", 50, null, false, position2Others);
    SequenceBuilder component3 = new SequenceBuilder().setName("query").asAlignment("ref", 100, null, false, position3Others);

    if (!builder.accepts(component1)) {
      fail("builder should accept component1");
    }
    builder.add(component1);
    if (builder.isComplete()) {
      fail("builder should not be complete after adding first component");
    }
    if (!builder.accepts(component2)) {
      fail("builder should accept component2");
    }
    builder.add(component2);
    if (builder.isComplete()) {
      fail("builder should not be complete after adding second component");
    }
    if (!builder.accepts(component3)) {
      fail("builder should accept component3");
    }
    builder.add(component3);
    if (!builder.isComplete()) {
      fail("builder should be complete after adding last component");
    }
  }

  @Test
  public void testPairedSupplementaryAlignments() {
    SamAlignment_Builder builder = new SamAlignment_Builder();

    List<PositionDescriptor> mate1Others = new ArrayList<PositionDescriptor>();
    mate1Others.add(new PositionDescriptor("ref", 100, true));

    List<PositionDescriptor> mate2Component1Others = new ArrayList<PositionDescriptor>();
    mate2Component1Others.add(new PositionDescriptor("ref", 0, false));
    mate2Component1Others.add(new PositionDescriptor("ref", 150, true));
    List<PositionDescriptor> mate2Component2Others = new ArrayList<PositionDescriptor>();
    mate2Component2Others.add(new PositionDescriptor("ref", 100, true));

    SequenceBuilder mate1 = new SequenceBuilder().setName("query").asAlignment("ref", 0, null, false, mate1Others);
    SequenceBuilder mate2Component1 = new SequenceBuilder().setName("query").asAlignment("ref", 100, null, true, mate2Component1Others);
    SequenceBuilder mate2Component2 = new SequenceBuilder().setName("query").asAlignment("ref", 150, null, true, mate2Component2Others);

    if (!builder.accepts(mate1)) {
      fail("builder should accept mate1");
    }
    builder.add(mate1);
    if (builder.isComplete()) {
      fail("builder should not be complete after adding first mate");
    }
    if (!builder.accepts(mate2Component1)) {
      fail("builder should accept mate2Component1");
    }
    builder.add(mate2Component1);
    if (builder.isComplete()) {
      fail("builder should not be complete after adding mate2Component1");
    }
    if (!builder.accepts(mate2Component2)) {
      fail("builder should accept mate2Component2");
    }
    builder.add(mate2Component2);
    if (!builder.isComplete()) {
      fail("builder should be complete after adding mate2Component2");
    }
  }

  @Test
  public void testPairedSupplementaryAlignmentsSamePosition() {
    SamAlignment_Builder builder = new SamAlignment_Builder();

    List<PositionDescriptor> mate1Component1Others = new ArrayList<PositionDescriptor>();
    mate1Component1Others.add(new PositionDescriptor("ref", 150, false));
    mate1Component1Others.add(new PositionDescriptor("ref", 100, true));
    List<PositionDescriptor> mate1Component2Others = new ArrayList<PositionDescriptor>();
    mate1Component2Others.add(new PositionDescriptor("ref", 0, false));
 
    List<PositionDescriptor> mate2Component1Others = new ArrayList<PositionDescriptor>();
    mate2Component1Others.add(new PositionDescriptor("ref", 0, false));
    mate2Component1Others.add(new PositionDescriptor("ref", 150, true));
    List<PositionDescriptor> mate2Component2Others = new ArrayList<PositionDescriptor>();
    mate2Component2Others.add(new PositionDescriptor("ref", 100, true));

    SequenceBuilder mate1Component1 = new SequenceBuilder().setName("query").asAlignment("ref", 0, null, false, mate1Component1Others);
    SequenceBuilder mate1Component2 = new SequenceBuilder().setName("query").asAlignment("ref", 150, null, false, mate1Component2Others);
    SequenceBuilder mate2Component1 = new SequenceBuilder().setName("query").asAlignment("ref", 100, null, true, mate2Component1Others);
    SequenceBuilder mate2Component2 = new SequenceBuilder().setName("query").asAlignment("ref", 150, null, true, mate2Component2Others);

    if (!builder.accepts(mate1Component1)) {
      fail("builder should accept mate1");
    }
    builder.add(mate1Component1);
    if (builder.isComplete()) {
      fail("builder should not be complete after adding mate1Component1");
    }
    builder.add(mate1Component2);
    if (builder.isComplete()) {
      fail("builder should not be complete after adding mate1Component2");
    }
    if (!builder.accepts(mate2Component1)) {
      fail("builder should accept mate2Component1");
    }
    builder.add(mate2Component1);
    if (builder.isComplete()) {
      fail("builder should not be complete after adding mate2Component1");
    }
    if (!builder.accepts(mate2Component2)) {
      fail("builder should accept mate2Component2");
    }
    builder.add(mate2Component2);
    if (!builder.isComplete()) {
      fail("builder should be complete after adding mate2Component2");
    }
  }

  private void fail(String message) {
    Assert.fail(message);
  }
}
