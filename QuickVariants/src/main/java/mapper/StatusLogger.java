package mapper;

import java.util.List;

// A StatusLogger is used for logging status messages
// A StatusLogger might choose not to log all messages if there are too many status updates consecutively
public class StatusLogger {
  public StatusLogger(Logger logger) {
    this.logger = logger;
    this.lastLoggedAt = 0;
  }

  // If no message was received recently or if this message is important, then this method will log it
  public void log(String message, boolean important) {
    long now = System.currentTimeMillis();
    if (now - lastLoggedAt > 1000 || important) {
      this.logger.log(message);
      this.lastLoggedAt = now;
    }
  }

  Logger logger;
  long lastLoggedAt;
}
