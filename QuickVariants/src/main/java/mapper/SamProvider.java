package mapper;

import java.io.IOException;

// A SamProvider lists lines in a .sam file
public interface SamProvider {
  SequenceBuilder getNextSequence() throws IOException;
  public boolean isReordered();
}
