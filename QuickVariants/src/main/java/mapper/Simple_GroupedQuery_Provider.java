package mapper;

import java.util.ArrayList;
import java.util.List;

// A Simple_GroupedQuery_Provider combines adjacent queries with the same name into groups
public class Simple_GroupedQuery_Provider implements GroupedQuery_Provider {
  public Simple_GroupedQuery_Provider(QueryProvider queryProvider) {
    this.queryProvider = queryProvider;
  }

  public List<QueryBuilder> getNextGroup() {
    List<QueryBuilder> group = null;
    while (true) {
      QueryBuilder query = this.peekNextQuery();
      if (query == null)
        break;
      if (group == null)
        group = new ArrayList<QueryBuilder>();
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

  private QueryBuilder peekNextQuery() {
    if (this.pendingQuery == null)
      this.pendingQuery = this.queryProvider.getNextQueryBuilder();
    return this.pendingQuery;
  }

  private void consumeQuery() {
    this.pendingQuery = null;
  }

  private QueryProvider queryProvider;
  private QueryBuilder pendingQuery;
}
