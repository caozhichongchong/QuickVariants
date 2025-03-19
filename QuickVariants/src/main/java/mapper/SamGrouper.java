package mapper;

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Queue;

// A SamGrouper reorders Sam records so that related records are next to each other
public class SamGrouper implements SequenceProvider {
  public SamGrouper(SequenceProvider reader) {
    this.reader = reader;
  }

  public SequenceBuilder getNextSequence() {
    while (true) {
      // check for a pending group
      if (this.readGroup != null) {
        SequenceBuilder nextInGroup = this.readGroup.poll();
        if (nextInGroup != null)
          return nextInGroup;
      }
      if (this.reader == null) {
        // if there is no pending group and nothing left to read, we're done
        return null;
      }

      // read next sequence
      SequenceBuilder sequence = this.reader.getNextSequence();
      if (sequence == null) {
        // we're done reading sequences, so we emit any remaining partial groups
        this.readGroup = this.flattenGroups(this.partiallyReadGroups);
        this.partiallyReadGroups = null;
        this.reader = null;
        this.printMetrics();
        continue;
      }

      // if this read doesn't expect a mate, we can just emit it now
      if (!sequence.getExpectsMate())
        return sequence;

      // if we found a read that expects a mate, we might need to wait for its mate
      String sequenceName = sequence.getName();
      SequenceBuilder mate = this.partiallyReadGroups.remove(sequenceName);
      if (mate != null) {
        // we completed a group so we can emit it now
        this.readGroup.add(mate);
        this.readGroup.add(sequence);
      } else {
        this.partiallyReadGroups.put(sequenceName, sequence);
        int numPendingGroups = this.partiallyReadGroups.size();
        if (this.maxNumPendingGroups < numPendingGroups) {
          this.maxNumPendingGroups = numPendingGroups;
        }
        if (numPendingGroups >= this.numPendingGroups_warningThreshold) {
          System.out.println("Many (>=" + numPendingGroups + ") reads not adjacent to their mates - this can reduce performance");
          this.numPendingGroups_warningThreshold = numPendingGroups * 4;
        }
      }
    }
  }

  public boolean isReordered() {
    return true;
  }

  public int getNumErrors() {
    return 0;
  }

  public boolean get_allReadsContainQualityInformation() {
    return false;
  }

  private void printMetrics() {
    if (this.maxNumPendingGroups > 1)
      System.out.println("Maximum number of sequences waiting for a mate at once: " + this.maxNumPendingGroups);
  }

  private Queue<SequenceBuilder> flattenGroups(Map<String, SequenceBuilder> groups) {
    Queue<SequenceBuilder> result = new ArrayDeque<SequenceBuilder>();
    for (SequenceBuilder sequence: groups.values()) {
      result.add(sequence);
    }
    return result;
  }

  SequenceProvider reader;
  Map<String, SequenceBuilder> partiallyReadGroups = new HashMap<String, SequenceBuilder>();
  Queue<SequenceBuilder> readGroup = new ArrayDeque<SequenceBuilder>();
  int maxNumPendingGroups;
  int numPendingGroups_warningThreshold = 20000;
}
