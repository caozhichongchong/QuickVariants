package mapper;

public class MutationDetectionParameters {
  public float minSNPTotalDepth = 5;
  public float minSNPDepthFraction = (float)0.9;

  public float minIndelTotalStartDepth = 1;
  public float minIndelStartDepthFraction = (float)0.8;

  public float minIndelContinuationTotalDepth = 1;
  public float minIndelContinuationDepthFraction = (float)0.7;
}
