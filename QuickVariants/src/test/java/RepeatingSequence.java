package mapper;

public class RepeatingSequence extends Sequence {
  public RepeatingSequence(String name, char repeatingContent, int length) {
    super(name, "", length, "");
    this.repeatingContent = Basepairs.encode(repeatingContent);
  }

  @Override
  public byte encodedCharAt(int i) {
    return repeatingContent;
  }

  private byte repeatingContent;
}
