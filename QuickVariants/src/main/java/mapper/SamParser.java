package mapper;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

// A SamParser processes a .sam file: checks that query sequence mates exist, and reports errors
public class SamParser implements SequenceProvider {
  public SamParser(SequenceProvider reader, String path, boolean readerReordered) {
    this.reader = reader;
    this.path = path;
    this.readerReordered = readerReordered;
  }

  public SequenceBuilder getNextSequence() {
    if (this.group.size() < 1) {
      try {
        this.readGroup();
      } catch(IOException e) {
        throw new RuntimeException(e);
      }
    }

    if (this.group.size() < 1) {
      return null;
    }
    SequenceBuilder result = this.group.get(0);
    this.group.remove(0);
    return result;
  }

  public boolean get_allReadsContainQualityInformation() {
    return false;
  }

  private void readGroup() throws IOException {
    if (this.group.size() < 1) {
      this.group = this.nextGroup;
      this.nextGroup = new ArrayList<SequenceBuilder>();
    }
    while (true) {
      SequenceBuilder mate1 = this.getNextSequenceBuilder();
      if (mate1 == null)
        return;
      boolean sameGroup = this.group.size() < 1 || mate1.getName().equals(this.group.get(0).getName());
      if (mate1.getExpectsMate()) {
        SequenceBuilder mate2 = this.peekNextSequenceBuilder();
        if (mate2 == null) {
          if (sameGroup) {
            // If this alignment is supposed to be paired but doesn't actually have a mate, ignore it
          } else {
            nextGroup.add(mate1);
          }
          return;
        }
        if (mate1.getName().equals(mate2.getName())) {
          mate1.withAlignmentMate(mate2);
          mate2.withAlignmentMate(mate1);
          this.pendingSequence = null;
        } else {
          if (sameGroup) {
            // If this alignment is supposed to be paired but doesn't actually have a mate, ignore it
            group.add(mate1);
          } else {
            nextGroup.add(mate1);
          }
          if (this.numErrors < 1) {
            StringBuilder errorBuilder = new StringBuilder();
            errorBuilder.append("\n");
            if (this.readerReordered) {
              errorBuilder.append("Error reading " + this.path + ": sequences specify mates, but no mates were found in the file for sequence:\n");
              errorBuilder.append(" " + mate1.getName() + "\n");
            } else {
              errorBuilder.append("Error reading " + this.path + ": sequences specify mates, but adjacent non-exempt lines have different names:\n");
              errorBuilder.append(" " + mate1.getName() + "\n");
              errorBuilder.append(" " + mate2.getName() + "\n");
              errorBuilder.append("You can consider replacing '--in-ordered-sam' with '--in-unordered-sam'");
            }
            errorBuilder.append("\n");
            System.out.println(errorBuilder.toString());
          }
          this.numErrors++;
          return;
        }
        if (sameGroup) {
          group.add(mate1);
          group.add(mate2);
        } else {
          nextGroup.add(mate1);
          nextGroup.add(mate2);
          return;
        }
      } else {
        if (sameGroup) {
          group.add(mate1);
        } else {
          nextGroup.add(mate1);
          return;
        }
      }
    }
  }

  private void readNextSequence() throws IOException {
    this.pendingSequence = this.reader.getNextSequence();
  }

  private SequenceBuilder peekNextSequenceBuilder() throws IOException {
    if (this.pendingSequence == null) {
      this.readNextSequence();
    }
    return this.pendingSequence;
  }

  private SequenceBuilder getNextSequenceBuilder() throws IOException {
    SequenceBuilder result = this.peekNextSequenceBuilder();
    this.pendingSequence = null;
    return result;
  }

  public int getNumErrors() {
    return this.numErrors;
  }

  @Override
  public String toString() {
    if (this.readerReordered) {
      return "reorder(" + this.path + ")";
    }
    return this.path;
  }

  SequenceProvider reader;
  boolean readerReordered = false;
  String path;
  private List<SequenceBuilder> group = new ArrayList<SequenceBuilder>();
  private List<SequenceBuilder> nextGroup = new ArrayList<SequenceBuilder>();
  private SequenceBuilder pendingSequence;
  int numErrors = 0;
}
