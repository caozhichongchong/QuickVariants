package mapper;

public interface QueryProvider {
  QueryBuilder getNextQueryBuilder();
  boolean get_allReadsContainQualityInformation();
  public int getNumErrors();
}
