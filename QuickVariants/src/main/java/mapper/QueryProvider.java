package mapper;

import java.util.ArrayList;
import java.util.List;

// A QueryProvider converts from SequenceBuilder to SamAlignment_Builder
// It reads sequences, matches adjacent paired-end reads, and treats each result as a query
public class QueryProvider {
  public QueryProvider(SequenceProvider sequenceProvider) {
    this.sequenceProvider = sequenceProvider;
  }

  public SamAlignment_Builder getNextSamAlignment_Builder() {
    SequenceBuilder builder = this.sequenceProvider.getNextSequence();
    if (builder == null) {
      return null;
    }
    if (builder.getHasAlignmentMate()) {
      SequenceBuilder builder2 = this.sequenceProvider.getNextSequence();
      if (builder2.getHasAlignmentMate()) {
        List<SequenceBuilder> sequenceBuilders = new ArrayList<SequenceBuilder>();
        sequenceBuilders.add(builder);
        sequenceBuilders.add(builder2);
        return new SamAlignment_Builder(sequenceBuilders, 1, 1);
      }
    }

    return new SamAlignment_Builder(builder);
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
