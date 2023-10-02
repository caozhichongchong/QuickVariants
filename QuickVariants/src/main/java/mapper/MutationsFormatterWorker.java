package mapper;

import java.util.ArrayList;
import java.util.List;

public class MutationsFormatterWorker extends Thread {
  public MutationsFormatterWorker(MutationDetectionParameters parameters) {
    this.parameters = parameters;
  }

  public void request(MutationsFormatRequest formatRequest) {
    this.formatRequest = formatRequest;
  }

  public String getResults() {
    return results;
  }

  @Override
  public void run() {
    this.results = this.format(this.formatRequest);
  }

  private String format(MutationsFormatRequest formatRequest) {
    StringBuilder stringBuilder = new StringBuilder();

    String sequenceName = formatRequest.sequence.getName();
    Alignments bases = formatRequest.alignments;
    int endIndex = formatRequest.startIndex + formatRequest.length;
    List<AlignmentPosition> deletions = new ArrayList<AlignmentPosition>();
    int displayIndex = 1;
    for (int i = formatRequest.startIndex; i < endIndex; i++) {
      displayIndex = i + 1;
      AlignmentPosition frequencies = bases.getPosition(i);
      if (frequencies.getCount() > 0) {
        if (frequencies.hasAlternates()) {
          writeSNPs(sequenceName, displayIndex, frequencies, stringBuilder);
        }
        List<AlignmentPosition> candidateInsertions = new ArrayList<AlignmentPosition>();
        int insertionIndex = 0;
        while (true) {
          AlignmentPosition insertion = bases.getInsertion(i, insertionIndex);
          if (!insertion.hasAlternates())
            break;
          candidateInsertions.add(insertion);
          insertionIndex++;
        }
        if (candidateInsertions.size() > 0)
          writeInsertions(sequenceName, displayIndex, frequencies, candidateInsertions, stringBuilder);
      }
      boolean supportsDeletion;
      if (deletions.size() > 0) {
        supportsDeletion = supportsIndelContinuation(frequencies);
      } else {
        supportsDeletion = supportsIndelStart(frequencies);
      }
      if (supportsDeletion) {
        deletions.add(frequencies);
      } else {
        if (deletions.size() > 0) {
          int deletionDisplayIndex = displayIndex - (deletions.size() - 1);
          writeDeletions(sequenceName, deletionDisplayIndex, deletions, stringBuilder);
          deletions.clear();
        }
      }


    }
    // TODO: correctly handle deletions extending past the end of a block
    if (deletions.size() > 0) {
      int deletionDisplayIndex = displayIndex - (deletions.size() - 1);
      writeDeletions(sequenceName, deletionDisplayIndex, deletions, stringBuilder);
      deletions.clear();
    }
 
    return stringBuilder.toString();
  }

  private void writeSNPs(String sequenceName, int rowNumber, AlignmentPosition frequencies, StringBuilder stringBuilder) {
    char reference = frequencies.getReference();
    if (reference == 'N') {
      return;
    }
    String refString = "" + frequencies.getReference();
    char[] nonzeroAlternates = frequencies.getNonzeroAlternates();
    for (int i = 0; i < nonzeroAlternates.length; i++) {
      char alternate = nonzeroAlternates[i];
      if (alternate != '-') {
        float mutationDepth = frequencies.getAlternateCount(alternate);
        float totalDepth = frequencies.getCount();
        if (supportsSNP(mutationDepth, totalDepth))
          writeMutation(sequenceName, rowNumber, refString, "" + alternate, mutationDepth, totalDepth, stringBuilder);
      }
    }
  }

  private void writeMutation(String sequenceName, int rowNumber, String reference, String mutation, float mutationDepth, float totalDepth, StringBuilder stringBuilder) {
    stringBuilder.append(sequenceName);
    stringBuilder.append('\t');
    stringBuilder.append(Integer.toString(rowNumber));
    stringBuilder.append('\t');

    stringBuilder.append(reference);
    stringBuilder.append('\t');
    stringBuilder.append(mutation);
    stringBuilder.append('\t');
    String mutationDepthFormatted = AlignmentPosition.formatNumber(mutationDepth);
    stringBuilder.append(mutationDepthFormatted);
    stringBuilder.append('\t');
    String totalDepthFormatted = AlignmentPosition.formatNumber(totalDepth);
    stringBuilder.append(totalDepthFormatted);

    stringBuilder.append('\n');
  }

  private boolean supportsSNP(float mutationDepth, float totalDepth) {
    if (totalDepth < this.parameters.minSNPTotalDepth)
      return false;
    if (mutationDepth <= 0)
      return false;
    float mutationFraction = mutationDepth / totalDepth;
    if (mutationFraction < this.parameters.minSNPDepthFraction)
      return false;
    return true;
  }

  private boolean supportsIndelStart(AlignmentPosition frequencies) {
    float totalDepth = frequencies.getMiddleCount();
    if (totalDepth < this.parameters.minIndelTotalStartDepth) {
      return false;
    }
    float indelDepth;
    if (frequencies.getReference() == '-') {
      indelDepth = totalDepth - frequencies.getMiddleReferenceCount();
    } else {
      indelDepth = frequencies.getMiddleAlternateCount('-');
    }
    if (indelDepth <= 0)
      return false;
    float indelFraction = indelDepth / totalDepth;
    if (indelFraction < this.parameters.minIndelStartDepthFraction) {
      return false;
    }
    return true;
  }

  private boolean supportsIndelContinuation(AlignmentPosition frequencies) {
    float totalDepth = frequencies.getMiddleCount();
    if (totalDepth < this.parameters.minIndelContinuationTotalDepth)
      return false;

    float indelDepth;
    if (frequencies.getReference() == '-') {
      indelDepth = totalDepth;
    } else {
      indelDepth = frequencies.getMiddleAlternateCount('-');
    }
    if (indelDepth <= 0)
      return false;

    float indelFraction = indelDepth / totalDepth;
    if (indelFraction < this.parameters.minIndelContinuationDepthFraction)
      return false;
    return true; 
  }

  private void writeInsertions(String sequenceName, int rowNumber, AlignmentPosition frequencies, List<AlignmentPosition> insertions, StringBuilder stringBuilder) {
    StringBuilder reference = new StringBuilder();
    StringBuilder mutated = new StringBuilder();
    float totalDepth = frequencies.getMiddleCount();
    float insertionDepth = totalDepth;
    boolean insertionStarted = false;
    for (int i = 0; i < insertions.size(); i++) {
      AlignmentPosition insertion = insertions.get(i);
      boolean insert = false;
      if (insertionStarted) {
        insert = supportsIndelContinuation(insertion);
        if (!insert)
          break;
      } else {
        insert = supportsIndelStart(insertion);
        insertionStarted = insert;
      }
      if (insert) {
        reference.append('-');
        mutated.append(insertion.getMostPopularAlternate());
        insertionDepth = insertion.getMiddleCount() - insertion.getMiddleReferenceCount();
      }
    }
    if (insertionStarted) {
      writeMutation(sequenceName, rowNumber, reference.toString(), mutated.toString(), insertionDepth, totalDepth, stringBuilder);
    }
  }

  private void writeDeletions(String sequenceName, int displayIndex, List<AlignmentPosition> deletions, StringBuilder stringBuilder) {
    StringBuilder reference = new StringBuilder();
    StringBuilder mutated = new StringBuilder();
    float deletionDepth = -1;
    float totalDepth = -1;
    for (int i = 0; i < deletions.size(); i++) {
      AlignmentPosition deletion = deletions.get(i);
      reference.append(deletion.getReference());
      mutated.append('-');
      float totalDepthHere = deletion.getMiddleCount();
      float deletionDepthHere = deletion.getMiddleAlternateCount('-');
      if (deletionDepth < 0 || deletionDepth > deletionDepthHere)
        deletionDepth = deletionDepthHere;
      if (totalDepth < 0 || totalDepth > totalDepthHere)
        totalDepth = totalDepthHere;
    }
    writeMutation(sequenceName, displayIndex, reference.toString(), mutated.toString(), deletionDepth, totalDepth, stringBuilder);
  }

  MutationsFormatRequest formatRequest;
  String results;
  boolean includeNonMutations;
  MutationDetectionParameters parameters;
}
