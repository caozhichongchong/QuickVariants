package mapper;

import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

class AlignmentCache {
  public AlignmentCache() {
  }

  public List<SequenceAlignment> get(Sequence query) {
    return null; //return this.cache.get(query.getText());
  }

  public void addAlignment(Sequence query, List<SequenceAlignment> alignments) {
    //this.cache.put(query.getText(), alignments);
  }

  public int getUsage() {
    return this.cache.size();
  }

  private ConcurrentHashMap<String, List<SequenceAlignment>> cache = new ConcurrentHashMap<String, List<SequenceAlignment>>();
}
