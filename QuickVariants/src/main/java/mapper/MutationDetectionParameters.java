package mapper;

public class MutationDetectionParameters {
  public static MutationDetectionParameters defaultFilter() {
    MutationDetectionParameters result = new MutationDetectionParameters();
    result.minSNPTotalDepth = 5;
    result.minSNPDepthFraction = (float)0.9;

    result.minIndelTotalStartDepth = 1;
    result.minIndelStartDepthFraction = (float)0.8;

    result.minIndelContinuationTotalDepth = 1;
    result.minIndelContinuationDepthFraction = (float)0.7;

    return result;
  }

  public static MutationDetectionParameters emptyFilter() {
    return new MutationDetectionParameters();
  }

  public float minSNPTotalDepth;
  public float minSNPDepthFraction;

  public float minIndelTotalStartDepth;
  public float minIndelStartDepthFraction;

  public float minIndelContinuationTotalDepth;
  public float minIndelContinuationDepthFraction;

  public boolean supportsSNP(float mutationDepth, float totalDepth) {
    if (totalDepth < this.minSNPTotalDepth)
      return false;
    if (mutationDepth <= 0)
      return false;
    float mutationFraction = mutationDepth / totalDepth;
    if (mutationFraction < this.minSNPDepthFraction)
      return false;
    return true;
  }

  public boolean supportsIndelStart(AlignmentPosition frequencies) {
    float totalDepth = frequencies.getMiddleCount();
    if (totalDepth < this.minIndelTotalStartDepth) {
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
    if (indelFraction < this.minIndelStartDepthFraction) {
      return false;
    }
    return true;
  }

  public boolean supportsIndelContinuation(AlignmentPosition frequencies) {
    float totalDepth = frequencies.getMiddleCount();
    if (totalDepth < this.minIndelContinuationTotalDepth)
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
    if (indelFraction < this.minIndelContinuationDepthFraction)
      return false;
    return true;
  }


}
