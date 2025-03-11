package mapper;

import java.util.ArrayList;
import java.util.List;

public class BufferedWriter implements TextWriter {
  public BufferedWriter(TextWriter writer) {
    this.writer = writer;
  }

  public void write(String message) {
    this.components.add(message);
  }

  public void write(List<String> messages) {
    this.components.addAll(messages);
  }

  public void flush() {
    this.writer.write(this.components);
    this.components.clear();
  }

  List<String> components = new ArrayList<String>();
  TextWriter writer;
}
