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
}
