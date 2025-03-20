package mapper;

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Queue;

// A SamGroupOrderer reorders Sam records so that related records are next to each other
public class SamGroupOrderer implements SequenceProvider {
  public SamGroupOrderer(SequenceProvider reader) {
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

      // if we found a read that expects a mate, we might need to wait for its mate
      String sequenceName = sequence.getName();
      SamAlignment_Builder group = this.partiallyReadGroups.get(sequenceName);
      boolean groupAlreadyExisted;
      if (group == null) {
        groupAlreadyExisted = false;
        group = new SamAlignment_Builder();
      } else {
        groupAlreadyExisted = true;
      }
      group.add(sequence);
      if (group.isComplete()) {
        // we completed a group so we can emit it now
        this.readGroup.addAll(group.getComponents());
        if (groupAlreadyExisted)
          this.partiallyReadGroups.remove(sequenceName);
      } else {
        if (!groupAlreadyExisted) {
          this.partiallyReadGroups.put(sequenceName, group);
          int numPendingGroups = this.partiallyReadGroups.size();
          if (this.maxNumPendingGroups < numPendingGroups) {
            this.maxNumPendingGroups = numPendingGroups;
          }
          if (numPendingGroups >= this.numPendingGroups_warningThreshold) {
            this.numPendingGroups_warningThreshold = numPendingGroups * 4;
          }
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

  @Override
  public String toString() {
    return "reorder(" + this.reader.toString() + ")";
  }

  private void printMetrics() {
    if (this.maxNumPendingGroups > 1)
      System.out.println("Maximum number of sequences waiting for a mate at once: " + this.maxNumPendingGroups);
  }

  private Queue<SequenceBuilder> flattenGroups(Map<String, SamAlignment_Builder> groups) {
    Queue<SequenceBuilder> result = new ArrayDeque<SequenceBuilder>();
    for (SamAlignment_Builder alignmentBuilder: groups.values()) {
      result.addAll(alignmentBuilder.getComponents());
    }
    return result;
  }

  SequenceProvider reader;
  Map<String, SamAlignment_Builder> partiallyReadGroups = new HashMap<String, SamAlignment_Builder>();
  Queue<SequenceBuilder> readGroup = new ArrayDeque<SequenceBuilder>();
  int maxNumPendingGroups;
  int numPendingGroups_warningThreshold = 20000;
}
