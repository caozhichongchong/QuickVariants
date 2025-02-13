package mapper;

public class StderrWriter implements TextWriter {
  public StderrWriter() {
  }

  public void write(String message) {
    System.out.println(message);
  }

  public void flush() {
  }
}
