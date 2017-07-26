package au.edu.uq.imb.memesuite.io.html;

/**
 * An exception to indicate a problem parsing HTML data.
 */
public class HDataParseException extends Exception {
  public HDataParseException(String message) {
    super(message);
  }

  public HDataParseException(String message, Throwable cause) {
    super(message, cause);
  }

  public HDataParseException(Throwable cause) {
    super(cause);
  }
}
