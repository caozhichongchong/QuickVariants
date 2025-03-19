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
    if (this.sequenceProviders.size() >= 2)
      return true;
    if (this.sequenceProviders.get(0).getExpectsMate())
      return false;
    return true;
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
