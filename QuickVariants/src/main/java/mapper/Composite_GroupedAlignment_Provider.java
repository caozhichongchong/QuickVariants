package mapper;

import java.util.ArrayList;
import java.util.List;

// A Composite_GroupedAlignment_Provider iterates over a list of GroupedAlignment_Provider
public class Composite_GroupedAlignment_Provider implements GroupedAlignment_Provider {
  public Composite_GroupedAlignment_Provider(GroupedAlignment_Provider provider) {
    this.providers = new ArrayList<GroupedAlignment_Provider>(1);
    this.providers.add(provider);
  }

  public Composite_GroupedAlignment_Provider(List<GroupedAlignment_Provider> providers) {
    this.providers = providers;
  }

  public List<SamAlignment_Builder> getNextGroup() {
    while (this.nextIndex < this.providers.size()) {
      List<SamAlignment_Builder> next = this.providers.get(this.nextIndex).getNextGroup();
      if (next != null) {
        return next;
      }
      this.nextIndex++;
    }
    return null;
  }

  public boolean get_allReadsContainQualityInformation() {
    for (GroupedAlignment_Provider provider : this.providers) {
      if (!provider.get_allReadsContainQualityInformation()) {
        return false;
      }
    }
    return true;
  }

  public int getNumErrors() {
    int total = 0;
    for (GroupedAlignment_Provider provider: this.providers) {
      total += provider.getNumErrors();
    }
    return total;
  }

  int nextIndex;
  List<GroupedAlignment_Provider> providers;
}
