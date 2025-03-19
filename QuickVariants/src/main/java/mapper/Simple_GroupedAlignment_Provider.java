package mapper;

import java.util.ArrayList;
import java.util.List;

// A Simple_GroupedAlignment_Provider combines adjacent alignments of the same query into groups
public class Simple_GroupedAlignment_Provider implements GroupedAlignment_Provider {
  public Simple_GroupedAlignment_Provider(SamAlignment_Provider queryProvider) {
    this.queryProvider = queryProvider;
  }

  public List<SamAlignment_Builder> getNextGroup() {
    List<SamAlignment_Builder> group = null;
    while (true) {
      SamAlignment_Builder query = this.peekNextQuery();
      if (query == null)
        break;
      if (group == null)
        group = new ArrayList<SamAlignment_Builder>();
      if (group.size() < 1 || query.sameSequenceNames(group.get(0))) {
        group.add(query);
        this.consumeQuery();
      } else {
        break;
      }
    }
    return group;
  }

  @Override
  public String toString() {
    return queryProvider.toString();
  }

  public boolean get_allReadsContainQualityInformation() {
    return this.queryProvider.get_allReadsContainQualityInformation();
  }

  public int getNumErrors() {
    return this.queryProvider.getNumErrors();
  }

  private SamAlignment_Builder peekNextQuery() {
    if (this.pendingQuery == null)
      this.pendingQuery = this.queryProvider.getNextSamAlignment_Builder();
    return this.pendingQuery;
  }

  private void consumeQuery() {
    this.pendingQuery = null;
  }

  private SamAlignment_Provider queryProvider;
  private SamAlignment_Builder pendingQuery;
}
