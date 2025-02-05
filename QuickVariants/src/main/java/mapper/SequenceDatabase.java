package mapper;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;

public class SequenceDatabase {
  public SequenceDatabase(Collection<Sequence> sequences) {
    for (Sequence sequence : sequences) {
      add(sequence);
    }
    this.computeMetrics();
  }

  public SequenceDatabase(Collection<Sequence> sequences, boolean addReverseComplement) {
    for (Sequence sequence : sequences) {
      add(sequence);
      if (addReverseComplement) {
        add(sequence.reverseComplement());
      }
    }
    this.computeMetrics();
  }

  public SequenceDatabase(Sequence sequence) {
    this.add(sequence);
    this.computeMetrics();
  }

  private void add(Sequence sequence) {
    sequence.setId(this.sequences.size());
    this.sequences.add(sequence);
    if (sequence.getComplementedFrom() == null) {
      this.totalForwardSize += sequence.getLength();
    }
    this.totalForwardAndReverseSize += sequence.getLength();
  }
  public Sequence get(short sequenceId) {
    return this.sequences.get(sequenceId);
  }
  public List<Sequence> getAll() {
    return sequences;
  }
  private void computeMetrics() {
    computeEncodingSize();
  }
  private void computeEncodingSize() {
    this.numBitsToIdentifySequence = log2RoundUp(this.sequences.size());
    this.numBitsToIdentifyPosition = log2RoundUp(this.sequences.get(0).getLength());
    this.maxEncodableOffset = (1L << numBitsToIdentifyPosition) - 1;
    // we're going to pack this into a byte array, so we need to use at least 8 bits per item so we can determine how many items we have
    this.numBitsPerPosition = Math.max(8, this.numBitsToIdentifySequence + this.numBitsToIdentifyPosition);
    this.maxEncodableValue = (1L << numBitsPerPosition) - 1;
  }

  public int log2RoundUp(long value) {
    int numBits = 1;
    long exponent = 2;
    while (true) {
      if (exponent >= value) {
        return numBits;
      }
      numBits++;
      exponent *= 2;
    }
  }

  public int getDecodedLength(byte[] encoded) {
    return encoded.length * 8 / this.numBitsPerPosition;
  }

  public int getEncodedLength(int itemLength) {
    return (itemLength * this.numBitsPerPosition + 7) / 8;
  }

  public long encodePosition(Sequence sequence, int startIndex) {
    if (startIndex < 0 || startIndex >= sequence.getLength()) {
      throw new IllegalArgumentException("Invalid position " + startIndex + " for sequence " + sequence.getName() + " of length " + sequence.getLength());
    }
    long result = (sequence.getId() << this.numBitsToIdentifyPosition) + startIndex;
    if (result < 0 || result > this.maxEncodableValue) {
      throw new IllegalArgumentException("encoded " + sequence.getName() + "[" + startIndex + "] and received out-of-range result " + result);
    }
    return result;
  }

  public int getNumSequences() {
    return this.sequences.size();
  }

  public Sequence getSequence(int identifier) {
    return this.sequences.get(identifier);
  }

  public Sequence getSequence(String sequenceName) {
    if (sequenceName.equals("*")) {
      return null;
    }
    Sequence sequence = this.sequencesByName.get(sequenceName);
    if (sequence == null) {
      for (Sequence candidate: this.sequences) {
        if (candidate.getName().equals(sequenceName)) {
          sequence = candidate;
          break;
        }
      }
      if (sequence == null) {
        throw new IllegalArgumentException("Failed to find reference sequence with name '" + sequenceName + "'");
      }
      this.sequencesByName.put(sequenceName, sequence);
    }
    return sequence;
  }

  // the total size of all sequences in this database in the forward direction
  public long getTotalForwardSize() {
    return totalForwardSize;
  }

  public long getTotalForwardAndReverseSize() {
    return totalForwardAndReverseSize;
  }

  public Sequence getReverseComplement(Sequence sequence) {
    Sequence reverseComplement = sequence.getComplementedFrom();
    if (reverseComplement != null)
      return reverseComplement;
    int sequenceId = (int)sequence.getId();
    int complementId = sequenceId + 1 - (sequenceId % 2) * 2;
    reverseComplement = getSequence(complementId);
    if (reverseComplement.getComplementedFrom() != sequence) {
      throw new IllegalArgumentException("Reverse complement of " + sequence + " is unknown to this SequenceDatabase");
    }
    return reverseComplement;
  }

  public Map<String, String> getCacheKeys() {
    Map<String, String> keys = new HashMap<String, String>();
    keys.put("numSequences", "" + this.sequences.size());
    keys.put("totalSize", "" + totalForwardSize);
    keys.put("firstSequenceName", this.sequences.get(0).getName());
    keys.put("paths", this.getPathsString());
    if (this.ancestral)
      keys.put("ancestral", "true");
    return keys;
  }

  public void setAncestral() {
    this.ancestral = true;
  }

  private String getPathsString() {
    TreeSet<String> paths = new TreeSet<String>();
    for (Sequence sequence: sequences) {
      String path = sequence.getPath();
      if (path != null) {
        paths.add(path);
      }
    }
    StringBuilder builder = new StringBuilder();
    builder.append("{");
    boolean first = true;
    for (String path: paths) {
      if (first) {
        first = false;
      } else {
        builder.append(",");
      }
      File file = new File(path);
      String name = file.getName();
      builder.append(file.getName() + ":" + file.lastModified());
    }
    builder.append("}");
    return builder.toString();
  }

  private List<Sequence> sequences = new ArrayList<Sequence>();
  private long totalForwardSize;
  private long totalForwardAndReverseSize;
  int numBitsToIdentifySequence;
  int numBitsToIdentifyPosition;
  int numBitsPerPosition;
  long maxEncodableOffset;
  long maxEncodableValue;
  Map<String, Sequence> sequencesByName = new ConcurrentHashMap<String, Sequence>();
  private boolean ancestral;
}
