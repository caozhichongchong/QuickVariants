package mapper;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import org.junit.Assert;
import org.junit.Test;

public class SamReader_Test {
  @Test
  public void testInvalidFile() {
    String corruptedText = "SampleText";
    String failureMessage = null;
    try {
      parse(corruptedText);
    } catch (Exception e) {
      failureMessage = e.toString();
    }
    if (failureMessage == null) {
      fail("Did not detect corrupted input");
    }
    if (failureMessage.indexOf(corruptedText) < 0) {
      fail("Did not mention corrupted line in failure message: '" + failureMessage + "'");
    }
  }

  @Test
  public void testCigarN() {
    String inputText = "Read1\t0\tcontig1\t0\t255\t3M4N3M\t*\t*\t*\tACGCCT";
    List<Sequence> sequences = parse(inputText);
    if (sequences.size() != 1) {
      Assert.fail("Expected 1 sequence, not " + sequences.size());
    }
    Sequence sequence = sequences.get(0);
    SamAlignment samAlignment = (SamAlignment)sequence;

    Sequence referenceContig = new SequenceBuilder().setName("contig1").add("ACGTTTTCCCT").build();
    SequenceDatabase reference = new SequenceDatabase(referenceContig);

    SequenceAlignment sequenceAlignment = samAlignment.toSequenceAlignment(reference);
    int numSections = sequenceAlignment.getNumSections();
    if (numSections != 2) {
      Assert.fail("Expected 2 sections, not " + numSections);
    }
  }

  private List<Sequence> parse(String text) {
    try {
      SamReader reader = buildReader(text);
      List<Sequence> sequences = new ArrayList<Sequence>();
      while (true) {
        SequenceBuilder builder = reader.getNextSequence();
        if (builder == null)
          break;
        sequences.add(builder.build());
      }
      return sequences;
    } catch (IOException e) {
      throw new RuntimeException("Failed to read from string", e);
    }
  }

  private SamReader buildReader(String text) {
    StringReader stringReader = new StringReader(text);
    BufferedReader bufferedReader = new BufferedReader(stringReader);
    String path = "test";
    return new SamReader(bufferedReader, path);
  }


  private void fail(String message) {
    Assert.fail(message);
  }
}
