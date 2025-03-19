package mapper;

import java.util.ArrayList;
import java.util.List;

// A GroupedQuery_Provider combines adjacent queries with the same name into groups
public interface GroupedQuery_Provider {
  List<QueryBuilder> getNextGroup();
  boolean get_allReadsContainQualityInformation();
  int getNumErrors();
}
