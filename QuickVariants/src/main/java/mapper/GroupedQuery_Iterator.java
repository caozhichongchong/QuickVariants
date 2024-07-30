package mapper;

import java.util.ArrayList;
import java.util.List;

// A GroupedQuery_Iterator provides a list of groups of queries
public class GroupedQuery_Iterator {
  public GroupedQuery_Iterator(GroupedQuery_Provider provider) {
    this.providers = new ArrayList<GroupedQuery_Provider>(1);
    this.providers.add(provider);
  }

  public GroupedQuery_Iterator(List<GroupedQuery_Provider> providers) {
    this.providers = providers;
  }

  public List<QueryBuilder> getNextGroup() {
    while (this.nextIndex < this.providers.size()) {
      List<QueryBuilder> next = this.providers.get(this.nextIndex).getNextGroup();
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
