package mapper;

import java.util.ArrayList;
import java.util.List;

// An AlignmentPosition tells how many reads aligned to a specific position on a reference
public class AlignmentPosition {
  private static char[] emptyCharArray = new char[0];

  public AlignmentPosition(char referenceBase) {
    this.referenceBase = referenceBase;
  }

  public void replaceAlternateWithReference(char value) {
    if (value != this.referenceBase) {
      this.replaceAlternateWithReference(value, false, false);
      this.replaceAlternateWithReference(value, false, true);
      this.replaceAlternateWithReference(value, true, false);
      this.replaceAlternateWithReference(value, true, true);
    }
  }

  private void replaceAlternateWithReference(char value, boolean forward, boolean nearQueryEnd) {
    // clear the alternate count
    int alternate = getScaled(value, forward, nearQueryEnd);
    this.putScaled(value, 0, forward, nearQueryEnd);
    // add the old count to the reference count
    int reference = this.getScaled(this.referenceBase, forward, nearQueryEnd);
    this.putScaled(this.referenceBase, alternate + reference, forward, nearQueryEnd);
  }

  private int getScaled(char value, boolean forward, boolean nearQueryEnd) {
    AlignmentPosition_DirectionCounts container = getContainer(forward, nearQueryEnd);
    if (this.referenceBase == value) {
      return container.getScaledReference();
    } else {
      return container.getScaledAlternate(value);
    }
  }

  public void putScaled(char value, int scaledWeight, boolean forward, boolean nearQueryEnd) {
    AlignmentPosition_DirectionCounts container = getContainer(forward, nearQueryEnd);

    if (this.referenceBase == value) {
      container.putScaledReference(scaledWeight);
    } else {
      container.putScaledAlternate(value, scaledWeight);
    }
  }

  private AlignmentPosition_DirectionCounts getContainer(boolean forward, boolean nearQueryEnd) {
    AlignmentPosition_DirectionCounts container;
    if (forward) {
      if (nearQueryEnd) {
        return this.forwardEndCounts;
      } else {
        return this.forwardMiddleCounts;
      }
    } else {
     if (nearQueryEnd) {
        return this.reverseEndCounts;
      } else {
        return this.reverseMiddleCounts;
      }
    }
  }



  public void putSampleAlternateSequence(Sequence querySequence, int index, boolean isDeletion) {
    char alternate;
    if (isDeletion) {
      alternate = '-';
    } else {
      byte alternateEncoded = querySequence.encodedCharAt(index);
      if (Basepairs.isAmbiguous(alternateEncoded))
        return; // it should be uncommon for a query to have ambiguity so we ignore it for now
      alternate = Basepairs.decode(alternateEncoded);
    }
    int alternateIndex = this.indexForKey(alternate);

    if (this.sampleAlternateSequences == null) {
      int capacity = this.getAllKeys().length;
      this.sampleAlternateSequences = new Sequence[capacity];
      this.sampleAlternateIndices = new int[capacity];
    }

    this.sampleAlternateSequences[alternateIndex] = querySequence;
    int putIndex = index;
    if (isDeletion)
      putIndex *= -1;

    this.sampleAlternateIndices[alternateIndex] = putIndex;
  }

  public char getMostPopular() {
    float maxCount = this.getReferenceCount();
    char mostPopular = this.referenceBase;
    if (this.hasAlternates()) {
      for (int i = 0; i < getAllKeys().length; i++) {
        float count = this.getAlternateCount(i);
        if (count > maxCount) {
          maxCount = count;
          mostPopular = this.keyForIndex(i);
        }
      }
    }
    return mostPopular;
  }

  public char getMostPopularAlternate() {
    float maxCount = 0;
    char mostPopular = ' ';
    if (this.hasAlternates()) {
      for (int i = 0; i < getAllKeys().length; i++) {
        float count = this.getAlternateCount(i);
        if (count > maxCount || i == 0) {
          maxCount = count;
          mostPopular = this.keyForIndex(i);
        }
      }
    }
    return mostPopular;
  }

  public String getCounts(char b, boolean isQueryEnd) {
    Float forwardMiddleCount, forwardEndCount, reverseMiddleCount, reverseEndCount;
    if (this.referenceBase == b) {
      forwardMiddleCount = this.forwardMiddleCounts.getReferenceCount();
      forwardEndCount = this.forwardEndCounts.getReferenceCount();
      reverseMiddleCount = this.reverseMiddleCounts.getReferenceCount();
      reverseEndCount = this.reverseEndCounts.getReferenceCount();
    } else {
      forwardMiddleCount = this.forwardMiddleCounts.getAlternateCount(b);
      forwardEndCount = this.forwardEndCounts.getAlternateCount(b);
      reverseMiddleCount = this.reverseMiddleCounts.getAlternateCount(b);
      reverseEndCount = this.reverseEndCounts.getAlternateCount(b);
    }
    StringBuilder builder = new StringBuilder();
    if (isQueryEnd) {
      builder.append(formatNumber(forwardEndCount));
      builder.append(',');
      builder.append(formatNumber(reverseEndCount));
    } else {
      builder.append(formatNumber(forwardMiddleCount));
      builder.append(',');
      builder.append(formatNumber(reverseMiddleCount));
    }
    return builder.toString();
  }

  private float getReferenceCount(boolean forward, boolean nearQueryEnd) {
    if (forward) {
      if (nearQueryEnd)
        return this.forwardEndCounts.getReferenceCount();
      else
        return this.reverseEndCounts.getReferenceCount();
    } else {
      if (nearQueryEnd)
        return this.forwardMiddleCounts.getReferenceCount();
      else
        return this.reverseMiddleCounts.getReferenceCount();
    }
  }

  public float getMiddleReferenceCount() {
    return this.forwardMiddleCounts.getReferenceCount() + this.reverseMiddleCounts.getReferenceCount();
  }

  public static String formatNumber(float number) {
    int rounded = (int)number;
    if (rounded == number) {
      return Integer.toString(rounded);
    }
    return Float.toString(number);
  }

  public float getCount() {
    float total = this.getReferenceCount();
    if (this.hasAlternates()) {
      for (int i = 0; i < this.getAllKeys().length; i++) {
        total += this.getAlternateCount(i);
      }
    }
    return total;
  }

  public String formatCount() {
    return formatNumber(getCount());
  }

  public float getReferenceCount() {
    return this.forwardMiddleCounts.getReferenceCount() + this.forwardEndCounts.getReferenceCount() + this.reverseMiddleCounts.getReferenceCount() + this.reverseEndCounts.getReferenceCount();
  }

  public float getMiddleCount() {
    float total = this.getMiddleReferenceCount();
    if (this.hasAlternates()) {
      for (int i = 0; i < this.getAllKeys().length; i++) {
        total += this.getMiddleAlternateCount(i);
      }
    }
    return total;
  }

  public String formatAlternateCount(char c) {
    return formatNumber(getAlternateCount(c));
  }

  public float getAlternateCount(char c) {
    return getAlternateCount(indexForKey(c));
  }

  private float getAlternateCount(int index) {
    return this.forwardMiddleCounts.getAlternateCount(index) + this.forwardEndCounts.getAlternateCount(index) + this.reverseMiddleCounts.getAlternateCount(index) + this.reverseEndCounts.getAlternateCount(index);
  }

  public float getMiddleAlternateCount(char c) {
    return getMiddleAlternateCount(indexForKey(c));
  }

  public float getMiddleAlternateCount(int index) {
    return this.forwardMiddleCounts.getAlternateCount(index) + this.reverseMiddleCounts.getAlternateCount(index);
  }

  public boolean hasAlternates() {
    return this.forwardMiddleCounts.hasAlternates() || this.forwardEndCounts.hasAlternates() || this.reverseMiddleCounts.hasAlternates() || this.reverseEndCounts.hasAlternates();
  }

  private boolean hasAlternate(int index) {
    return this.forwardMiddleCounts.hasAlternate(index) || this.forwardEndCounts.hasAlternate(index) || this.reverseMiddleCounts.hasAlternate(index) || this.reverseEndCounts.hasAlternate(index);
  }

  public char[] getNonzeroAlternates() {
    return getNonzeroAlternates(getReference());
  }

  public char[] getNonzeroAlternates(char reference) {
    if (!this.hasAlternates()) {
      return emptyCharArray;
    }

    StringBuilder builder = new StringBuilder();

    char[] keys = getAllKeys();
    for (int i = 0; i < keys.length; i++) {
      if (hasAlternate(i)) {
        builder.append(keys[i]);
      }
    }
    return builder.toString().toCharArray();
  }

  public Sequence getSampleAlternateSequence(char havingValue) {
    if (this.sampleAlternateSequences == null)
      return null;
    return this.sampleAlternateSequences[indexForKey(havingValue)];
  }
  public int getSampleAlternateIndex(char havingValue) {
    return Math.abs(this.sampleAlternateIndices[indexForKey(havingValue)]);
  }
  public boolean isSampleAlternateDeletion(char havingValue) {
    return this.sampleAlternateIndices[indexForKey(havingValue)] < 0;
  }
  public char getReference() {
    return this.referenceBase;
  }

  private static char[] getAllKeys() {
    return AlignmentPosition_DirectionCounts.getAllKeys();
  }

  private static char keyForIndex(int index) {
    return AlignmentPosition_DirectionCounts.keyForIndex(index);
  }

  private static int indexForKey(char key) {
    return AlignmentPosition_DirectionCounts.indexForKey(key);
  }

  // The item on the reference
  private char referenceBase;
  // the number of forward or reverse matches/mismatches on the middle/end of a read
  private AlignmentPosition_DirectionCounts forwardMiddleCounts = new AlignmentPosition_DirectionCounts();
  private AlignmentPosition_DirectionCounts forwardEndCounts = new AlignmentPosition_DirectionCounts();
  private AlignmentPosition_DirectionCounts reverseMiddleCounts = new AlignmentPosition_DirectionCounts();
  private AlignmentPosition_DirectionCounts reverseEndCounts = new AlignmentPosition_DirectionCounts();

  private Sequence[] sampleAlternateSequences;
  private int[] sampleAlternateIndices;
}
