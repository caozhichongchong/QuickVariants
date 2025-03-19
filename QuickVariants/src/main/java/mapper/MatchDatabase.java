package mapper;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

// Stores a collection of alignedBlocks
public class MatchDatabase implements AlignmentListener {
  public MatchDatabase(double queryEndFraction) {
    this.queryEndFraction = queryEndFraction;
    this.alignmentsBySequence = new HashMap<Sequence, Alignments>();
  }

  public void addAlignments(List<List<QueryAlignment>> alignments) {
    Map<Sequence, List<WeightedAlignment>> alignmentsByReference = this.groupByReference(alignments);

    for (Map.Entry<Sequence, List<WeightedAlignment>> job : alignmentsByReference.entrySet()) {
      Sequence reference = job.getKey();
      Alignments alignmentsHere;
      synchronized (this.alignmentsBySequence) {
        alignmentsHere = this.getOrCreateAlignments(reference);
      }
      alignmentsHere.add(job.getValue());
    }
  }

  public void addUnaligned(List<SamAlignment> unalignedQueries) {
  }

  private Map<Sequence, List<WeightedAlignment>> groupByReference(List<List<QueryAlignment>> alignments) {
    Map<Sequence, List<WeightedAlignment>> alignmentsByReference = new HashMap<Sequence, List<WeightedAlignment>>();

    for (List<QueryAlignment> queryAlignments: alignments) {
      for (QueryAlignment queryAlignment: queryAlignments) {
        for (SequenceAlignment alignment: queryAlignment.getComponents()) {
          float weight = (float)alignment.weight;
          List<AlignedBlock> blocks = alignment.getSections();
          if (blocks.size() > 0) {
            Sequence reference = blocks.get(0).getSequenceB();
            List<WeightedAlignment> alignmentsOnThisRef = alignmentsByReference.get(reference);
            if (alignmentsOnThisRef == null) {
              alignmentsOnThisRef = new ArrayList<WeightedAlignment>();
              alignmentsByReference.put(reference, alignmentsOnThisRef);
            }
            alignmentsOnThisRef.add(new WeightedAlignment(alignment, queryAlignment, weight));
          }
        }
      }
    }
    return alignmentsByReference;
  }


  // Map from name of sequence to Alignments on that sequence
  public Map<Sequence, Alignments> groupByPosition() {
    return alignmentsBySequence;
  }

  private Alignments getOrCreateAlignments(Sequence reference) {
    Alignments alignments = alignmentsBySequence.get(reference);
    if (alignments == null) {
      alignments = new Alignments(reference, this.queryEndFraction);
      alignmentsBySequence.put(reference, alignments);
    }
    return alignments;
  }

  private Map<Sequence, Alignments> alignmentsBySequence;
  private double queryEndFraction;
}
