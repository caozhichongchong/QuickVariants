package mapper;

import java.util.ArrayList;
import java.util.List;

// a SamAlignment_Builder builds a SamAlignment
public class SamAlignment_Builder {

  public SamAlignment_Builder() {
    this.sequenceProviders = new ArrayList<SequenceBuilder>();
  }

  public SamAlignment_Builder(SequenceBuilder sequenceProvider) {
    this.sequenceProviders = new ArrayList<SequenceBuilder>();
    this.add(sequenceProvider);
  }

  public SamAlignment_Builder(List<SequenceBuilder> sequenceProviders) {
    this.sequenceProviders = sequenceProviders;
  }

  public void add(SequenceBuilder sequenceBuilder) {
    this.sequenceProviders.add(sequenceBuilder);
  }

  public boolean isComplete() {
    if (this.sequenceProviders.size() < 1)
      return false;
    for (SequenceBuilder component: this.sequenceProviders) {
      PositionDescriptor missing = this.getMissingPosition(component);
      if (missing != null) {
        return false;
      }
    }
    return true;
  }

  public String explainIncompleteness() {
    for (SequenceBuilder component: this.sequenceProviders) {
      PositionDescriptor missing = this.getMissingPosition(component);
      if (missing != null) {
        StringBuilder errorBuilder = new StringBuilder();
        errorBuilder.append("\n in group of size " + this.sequenceProviders.size() + ":\n");
        for (SequenceBuilder other: this.sequenceProviders) {
          errorBuilder.append("  " + other.getName() + " at " + other.getPosition() + "\n");
        }
        errorBuilder.append(" " + component.getName() + " at " + component.getPosition() + " expects an alignment at " + missing + " but none was found.");
        PositionDescriptor existingReverseComplement = getExistingReverseComplementPosition(missing);
        if (existingReverseComplement != null) {
          errorBuilder.append("\n (reverse complement: " + existingReverseComplement + " was found, however)");
        }
        return errorBuilder.toString();
      }
    }
    return "Internal error: SamAlignment_Builder unable to explain incompleteness of group.";
  }

  // looks for an existing position that matches the reverse complement of this one, and if it exists, returns it
  private PositionDescriptor getExistingReverseComplementPosition(PositionDescriptor position) {
    for (SequenceBuilder sequence: this.sequenceProviders) {
      PositionDescriptor other = sequence.getPosition();
      PositionDescriptor complemented = new PositionDescriptor(other.getSequenceName(), other.getStartIndex(), !other.isReverseComplemented());
      if (position.equals(complemented))
        return complemented;
    }
    return null;
  }

  private PositionDescriptor getMissingPosition(SequenceBuilder component) {
    List<PositionDescriptor> expectations = component.getOtherComponentPositions();
    if (expectations != null) {
      for (PositionDescriptor expected: expectations) {
        if (!containsComponent(expected, component))
          return expected;
      }
    }
    return null;
  }

  private boolean containsComponent(PositionDescriptor expectation, SequenceBuilder exclude) {
    for (SequenceBuilder sequence: this.sequenceProviders) {
      if (sequence != exclude && sequence.getPosition().equals(expectation)) {
        return true;
      }
    }
    return false;
  }

  public boolean accepts(SequenceBuilder component) {
    if (this.sequenceProviders.size() < 1)
      return true;
    if (this.isComplete())
      return false;
    if (this.sequenceProviders.get(0).getName().equals(component.getName()))
      return true;
    return false;
  }

  public SamAlignment build() {
    List<Sequence> sequences = new ArrayList<Sequence>(this.sequenceProviders.size());
    boolean odd = false;
    for (SequenceBuilder sequenceProvider : this.sequenceProviders) {
      sequences.add(sequenceProvider.build());
    }
    return new SamAlignment(sequences);
  }

  public void setId(long id) {
    for (SequenceBuilder builder : this.sequenceProviders) {
      builder.setId(id);
    }
  }

  public int getLength() {
    int total = 0;
    for (SequenceBuilder builder : this.sequenceProviders) {
      total += builder.getLength();
    }
    return total;
  }

  public boolean sameSequenceNames(SamAlignment_Builder other) {
    if (other.sequenceProviders.size() != this.sequenceProviders.size())
      return false;
    for (int i = 0; i < this.sequenceProviders.size(); i++) {
      if (!this.sequenceProviders.get(i).getName().equals(other.sequenceProviders.get(i).getName()))
        return false;
    }
    return true;
  }

  public boolean hasSequenceName(String name) {
    for (SequenceBuilder sequenceProvider: this.sequenceProviders) {
      if (name.equals(sequenceProvider.getName()))
        return true;
    }
    return false;
  }

  public List<SequenceBuilder> getComponents() {
    return this.sequenceProviders;
  }

  public SequenceBuilder getLastComponent() {
    if (this.sequenceProviders.size() > 0)
      return this.sequenceProviders.get(this.sequenceProviders.size() - 1);
    return null;
  }

  public void setWeight(double weight) {
    for (SequenceBuilder component: this.sequenceProviders) {
      component.setAlignmentWeight(weight);
    }
  }

  public boolean isNonempty() {
    return this.sequenceProviders.size() > 0;
  }

  private List<SequenceBuilder> sequenceProviders;
  private int id;
}
