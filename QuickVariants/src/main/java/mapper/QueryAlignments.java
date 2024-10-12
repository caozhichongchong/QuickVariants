package mapper;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

// A QueryAlignments lists the places where each piece of a Query may align
// It's mostly a List<QueryAlignment> with some extra metadata about whether we were able to align each Sequence
public class QueryAlignments {

  // makes a QueryAlignments that says there was only query and only one alignment
  public static QueryAlignments singleChoice(QueryAlignment alignment) {
    List<QueryAlignment> choices = new ArrayList<QueryAlignment>(1);
    choices.add(alignment);
    List<Sequence> sequences = new ArrayList<Sequence>(1);
    sequences.add(alignment.getSequenceA());
    return singleComponent(sequences, choices);
  }

  // makes a QueryAlignments that says no alignments were found
  public static QueryAlignments unaligned(List<Sequence> querySequences) {
    return QueryAlignments.singleComponent(querySequences, new ArrayList<QueryAlignment>(0));
  }

  // makes a QueryAlignments that says that there was one query (potentially having multiple sequences) and these are its possible alignments
  public static QueryAlignments singleComponent(List<Sequence> querySequences, List<QueryAlignment> choices) {
    List<List<QueryAlignment>> components = new ArrayList<List<QueryAlignment>>(1);
    components.add(choices);
    return new QueryAlignments(querySequences, components);
  }

  // makes a QueryAlignments that says the query was split into these pieces and here are the possible alignments for each piece
  public QueryAlignments(List<Sequence> querySequences, List<List<QueryAlignment>> alignments) {
    this.sequences = querySequences;
    this.alignments = alignments;
  }

  public List<QueryAlignment> getTopLevelAlignments() {
    if (this.getNumQueries() != 1) {
      return new ArrayList<QueryAlignment>(0);
    }
    return this.alignments.get(0);
  }

  public List<List<QueryAlignment>> getAlignments() {
    return this.alignments;
  }

  public List<QueryAlignment> getAlignments(int index) {
    return this.alignments.get(index);
  }

  public int getNumComponents() {
    return this.alignments.size();
  }

  public int getTotalOfAllComponents() {
    int total = 0;
    for (List<QueryAlignment> value: this.alignments) {
      total += value.size();
    }
    return total;
  }

  // Returns the number of subqueries that this alignment represents
  // If our query wasn't a paired-end read, this number should be 1
  // If our query was a paired-end read:
  //  If neither mate aligned, this number should be 1
  //  If both mates aligned together, this number should be 1
  //  If one query aligned and one didn't, this number should be 2
  //  If both aligned to different places, this number should be 2
  public int getNumQueries() {
    return this.alignments.size();
  }

  // Returns the number of queries for which we found an alignment.
  // If our query wasn't a paired-end read:
  //  this number should be 0 (unaligned) or 1 (aligned)
  // If our query was a paired-end read:
  //  If neither mate aligned anywhere, this should return 0
  //  If both mates aligned together, this should return 1
  //  If one mate aligned and one didn't, this should return 1
  //  If both mates aligned but not together, this should return 2
  public int getNumQueriesHavingAlignments() {
    int count = 0;
    for (List<QueryAlignment> alignments: this.alignments) {
      if (alignments.size() > 0) {
        count++;
      }
    }
    return count;
  }

  public List<QueryAlignment> getFirstAlignments() {
    return this.alignments.get(0);
  }

  public Sequence getSequence(int index) {
    return this.sequences.get(index);
  }

  public List<Sequence> getSequences() {
    return this.sequences;
  }

  public Sequence getFirstSequence() {
    return this.sequences.get(0);
  }

  public int getQueryLength(int index) {
    if (this.alignments.size() == 1) {
      // This represents one query, so the total query length for any component is the total query length
      int total = 0;
      for (Sequence sequence: this.sequences) {
        total += this.sequences.get(index).getLength();
      }
      return total;
    }
    // This represents the alignments for multiple subqueries, so the total query length for any component is the length of that component's sequence
    return this.sequences.get(index).getLength();
  }

  private List<Sequence> sequences;
  private List<List<QueryAlignment>> alignments;
}
