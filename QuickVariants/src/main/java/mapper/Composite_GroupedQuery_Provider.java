package mapper;

import java.util.ArrayList;
import java.util.List;

// A Composite_GroupedQuery_Provider iterates over a list of GroupedQuery_Provider
public class Composite_GroupedQuery_Provider implements GroupedQuery_Provider {
  public Composite_GroupedQuery_Provider(GroupedQuery_Provider provider) {
    this.providers = new ArrayList<GroupedQuery_Provider>(1);
    this.providers.add(provider);
  }

  public Composite_GroupedQuery_Provider(List<GroupedQuery_Provider> providers) {
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
    for (GroupedQuery_Provider provider : this.providers) {
      if (!provider.get_allReadsContainQualityInformation()) {
        return false;
      }
    }
    return true;
  }

  public int getNumErrors() {
    int total = 0;
    for (GroupedQuery_Provider provider: this.providers) {
      total += provider.getNumErrors();
    }
    return total;
  }

  int nextIndex;
  List<GroupedQuery_Provider> providers;
}
