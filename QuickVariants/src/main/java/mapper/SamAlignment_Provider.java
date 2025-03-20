package mapper;

import java.util.ArrayList;
import java.util.List;

// A SamAlignment_Provider converts from SequenceBuilder to SamAlignment_Builder
// It reads sequences, matches adjacent paired-end reads, and treats each result as a query
public class SamAlignment_Provider {
  public SamAlignment_Provider(SequenceProvider sequenceProvider, String path, boolean readerReordered) {
    this.sequenceProvider = sequenceProvider;
    this.path = path;
    this.readerReordered = readerReordered;
  }

  public SamAlignment_Builder getNextSamAlignment_Builder() {
    if (this.alreadyDone)
      return null;
    while (true) {
      SequenceBuilder sequenceBuilder = this.sequenceProvider.getNextSequence();
      if (sequenceBuilder == null) {
        SamAlignment_Builder result = this.pendingAlignment;
        this.pendingAlignment = null;
        if (result == null)
          this.done();
        return result;
      }
      if (!this.pendingAlignment.accepts(sequenceBuilder)) {
        SamAlignment_Builder result = this.pendingAlignment;
        this.pendingAlignment = new SamAlignment_Builder();
        this.pendingAlignment.add(sequenceBuilder);
        if (!result.isComplete()) {
          if (this.numIncompleteGroups < 1) {
            StringBuilder errorBuilder = new StringBuilder();
            errorBuilder.append("\n");
            if (this.readerReordered) {
              errorBuilder.append("Error reading " + this.path + ": " + result.explainIncompleteness());
            } else {
              errorBuilder.append("Error reading " + this.path + ": " + result.explainIncompleteness());
              errorBuilder.append("Previous alignment name:\n");
              errorBuilder.append(" " + result.getLastComponent().getName() + "\n");
              errorBuilder.append("Next alignment name:\n");
              errorBuilder.append(" " + sequenceBuilder.getName() + "\n");
              errorBuilder.append("You can consider replacing '--in-ordered-sam' with '--in-unordered-sam'");
            }
            errorBuilder.append("\n");
            System.out.println(errorBuilder.toString());
          }
          this.numIncompleteGroups++;
        }
        return result;
      }
      this.pendingAlignment.add(sequenceBuilder);
    }
  }

  public boolean get_allReadsContainQualityInformation() {
    return this.sequenceProvider.get_allReadsContainQualityInformation();
  }

  public int getNumErrors() {
    return this.numIncompleteGroups + this.sequenceProvider.getNumErrors();
  }

  @Override
  public String toString() {
    return this.sequenceProvider.toString();
  }

  private void done() {
    if (this.numIncompleteGroups > 0) {
      System.out.println("" + this.numIncompleteGroups + " incomplete groups (missing mates or supplementary alignments)");
    }
    this.alreadyDone = true;
  }

  private SamAlignment_Builder startNew(SequenceBuilder firstComponent) {
    SamAlignment_Builder result = this.pendingAlignment;
    this.pendingAlignment = new SamAlignment_Builder();
    this.pendingAlignment.add(firstComponent);
    return result;
  }
  
  private SequenceProvider sequenceProvider;
  private SamAlignment_Builder pendingAlignment = new SamAlignment_Builder();
  private int numIncompleteGroups;
  private String path;
  private boolean readerReordered;
  private boolean alreadyDone;
}
