package mapper;

import java.util.ArrayList;
import java.util.List;
import org.junit.Assert;
import org.junit.Test;


public class SequenceDatabase_Test {
  public SequenceDatabase_Test() {
  }

  @Test
  public void testEncodingLargeSequences() {
    int sequenceSize = (int)Math.pow(2, 30);
    List<Sequence> reference = makeLargeSequences(16, sequenceSize);
    SequenceDatabase database = new SequenceDatabase(reference);
    for (int i = 0; i < reference.size(); i++) {
      testEncoding(database, reference.get(i), 0);
      testEncoding(database, reference.get(i), 100);
      testEncoding(database, reference.get(i), reference.get(i).getLength() - 100);
      testEncoding(database, reference.get(i), reference.get(i).getLength() - 1);
    }
  }

  @Test
  public void testEncodingManyLargeSequences() {
    int numSequences = (int)Math.pow(2, 13);
    int sequenceSize = (int)Math.pow(2, 21);
    List<Sequence> reference = makeLargeSequences(numSequences, sequenceSize);
    SequenceDatabase database = new SequenceDatabase(reference);
    for (int i = 0; i < reference.size(); i++) {
      testEncoding(database, reference.get(i), 0);
      testEncoding(database, reference.get(i), 100);
      testEncoding(database, reference.get(i), reference.get(i).getLength() - 100);
      testEncoding(database, reference.get(i), reference.get(i).getLength() - 1);
    }
  }

  private List<Sequence> makeLargeSequences(int numSequences, int sequenceLength) {
    // Make SequenceDatabase
    List<Sequence> reference = new ArrayList<Sequence>();
    for (int i = 0; i < numSequences; i++) {
      reference.add(makeLargeSequence(i, sequenceLength - i));
    }
    return reference;
  }

  private Sequence makeLargeSequence(int identifier, int length) {
    return new RepeatingSequence("" + identifier, 'A', length);
  }

  private void testEncoding(SequenceDatabase sequenceDatabase, Sequence sequence, int position) {
    long encoded = sequenceDatabase.encodePosition(sequence, position);
    SequencePosition decoded = sequenceDatabase.decodePosition(encoded);
    if (decoded.getSequence() != sequence || decoded.getStartIndex() != position) {
      fail("pre-decoded value " + sequence.getName() + "[" + position + "] != post-decoded value " + decoded.getSequence().getName() + "[" + decoded.getStartIndex() + "]");
    }
  }

  private void fail(String message) {
    Assert.fail(message);
  }
}
