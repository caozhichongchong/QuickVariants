package mapper;

import java.util.ArrayList;
import java.util.List;

// A SamAlignment_Provider converts from SequenceBuilder to SamAlignment_Builder
// It reads sequences, matches adjacent paired-end reads, and treats each result as a query
public class SamAlignment_Provider {
  public SamAlignment_Provider(SequenceProvider sequenceProvider) {
    this.sequenceProvider = sequenceProvider;
  }

  public SamAlignment_Builder getNextSamAlignment_Builder() {
    SamAlignment_Builder alignmentBuilder = new SamAlignment_Builder();
    while (!alignmentBuilder.isComplete()) {
      SequenceBuilder sequenceBuilder = this.sequenceProvider.getNextSequence();
      if (sequenceBuilder == null)
        return null;
      alignmentBuilder.add(sequenceBuilder);
    }
    return alignmentBuilder;
  }

  public boolean get_allReadsContainQualityInformation() {
    return this.sequenceProvider.get_allReadsContainQualityInformation();
  }

  public int getNumErrors() {
    return this.sequenceProvider.getNumErrors();
  }

  @Override
  public String toString() {
    return this.sequenceProvider.toString();
  }
  
  private SequenceProvider sequenceProvider;
}
