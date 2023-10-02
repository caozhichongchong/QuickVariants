package mapper;

import java.util.ArrayList;
import java.util.List;

// A PairedEndQueryProvider generates Query objects from Illumina-style paired-end reads
public class PairedEndQueryProvider implements QueryProvider {
  public PairedEndQueryProvider(SequenceProvider leftsProvider, SequenceProvider rightsProvider, double expectedInnerDistance, double spacingDeviationPerUnitPenalty) {
    this.sequenceProviders = new ArrayList<SequenceProvider>(1);
    this.sequenceProviders.add(leftsProvider);
    this.sequenceProviders.add(rightsProvider);
    this.expectedInnerDistance = expectedInnerDistance;
    this.spacingDeviationPerUnitPenalty = spacingDeviationPerUnitPenalty;
  }

  public QueryBuilder getNextQueryBuilder() {
    int firstLength = 0;
    boolean first = true;
    List<SequenceBuilder> components = new ArrayList<SequenceBuilder>(this.sequenceProviders.size());
    for (SequenceProvider provider : this.sequenceProviders) {
      SequenceBuilder builder = provider.getNextSequence();
      if (builder == null)
        return null;
      if (first) {
        firstLength = builder.getLength();
        first = false;
      }
      components.add(builder);
    }
    // Choose an upper bound on the maximum possible offset between sequences
    return new QueryBuilder(components, this.expectedInnerDistance, this.spacingDeviationPerUnitPenalty);
  }

  @Override
  public String toString() {
    return "paired queries: " + this.sequenceProviders.get(0).toString() + ", " + this.sequenceProviders.get(1).toString();
  }

  public boolean get_allReadsContainQualityInformation() {
    for (SequenceProvider sequenceProvider : this.sequenceProviders) {
      if (!sequenceProvider.get_allReadsContainQualityInformation()) {
        return false;
      }
    }
    return true;
  }

  public int getNumErrors() {
    int total = 0;
    for (SequenceProvider provider: this.sequenceProviders) {
      total += provider.getNumErrors();
    }
    return total;
  }

  private List<SequenceProvider> sequenceProviders;
  private double expectedInnerDistance;
  private double spacingDeviationPerUnitPenalty;
}
