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

  private List<Sequence> parse(String text) throws IOException {
    SamReader reader = buildReader(text);
    List<Sequence> sequences = new ArrayList<Sequence>();
    while (true) {
      SequenceBuilder builder = reader.getNextSequence();
      if (builder == null)
        break;
      sequences.add(builder.build());
    }
    return sequences;
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
