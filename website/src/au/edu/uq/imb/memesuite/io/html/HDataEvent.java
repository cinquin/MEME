package au.edu.uq.imb.memesuite.io.html;

/**
* Created with IntelliJ IDEA.
* User: james
* Date: 24/06/15
* Time: 10:43 AM
* To change this template use File | Settings | File Templates.
*/
public class HDataEvent {
  protected HDataEventType type;
  protected String name;
  protected Object value;

  protected HDataEvent(HDataEventType type, String name, Object value) {
    this.type = type;
    this.name = name;
    this.value = value;
  }

  public HDataEventType getType() {
    return type;
  }

  public String getName() {
    if (name == null) throw new NullPointerException("Current event does not have a name!");
    return name;
  }

  public Object getValue() {
    if (type != HDataEventType.JSON_VALUE && type != HDataEventType.HIDDEN_FIELD)
      throw new NullPointerException("Current event does not have a value!");
    return value;
  }

  public boolean isNull() {
    if (type != HDataEventType.JSON_VALUE && type != HDataEventType.HIDDEN_FIELD)
      throw new NullPointerException("Current event does not have a value!");
    return value == null;
  }

  public boolean isBoolean() {
    if (type != HDataEventType.JSON_VALUE && type != HDataEventType.HIDDEN_FIELD)
      throw new NullPointerException("Current event does not have a value!");
    return value != null && value instanceof Boolean;
  }

  public boolean isString() {
    if (type != HDataEventType.JSON_VALUE && type != HDataEventType.HIDDEN_FIELD)
      throw new NullPointerException("Current event does not have a value!");
    return value != null && value instanceof String;
  }

  public boolean isNumber() {
    if (type != HDataEventType.JSON_VALUE && type != HDataEventType.HIDDEN_FIELD)
      throw new NullPointerException("Current event does not have a value!");
    return value != null && value instanceof Number;
  }

  public static HDataEvent createFileStart() {
    return new HDataEvent(HDataEventType.FILE_START, null, null);
  }
  public static HDataEvent createHiddenField(String name, String value) {
    return new HDataEvent(HDataEventType.HIDDEN_FIELD, name, value);
  }
  public static HDataEvent createJsonStartData(String name) {
    return new HDataEvent(HDataEventType.JSON_START_DATA, name, null);
  }
  public static HDataEvent createJsonEndData() {
    return new HDataEvent(HDataEventType.JSON_END_DATA, null, null);
  }
  public static HDataEvent createJsonStartObject() {
    return new HDataEvent(HDataEventType.JSON_START_OBJECT, null, null);
  }
  public static HDataEvent createJsonEndObject() {
    return new HDataEvent(HDataEventType.JSON_END_OBJECT, null, null);
  }
  public static HDataEvent createJsonStartList() {
    return new HDataEvent(HDataEventType.JSON_START_LIST, null, null);
  }
  public static HDataEvent createJsonEndList() {
    return new HDataEvent(HDataEventType.JSON_END_LIST, null, null);
  }
  public static HDataEvent createJsonStartProperty(String name) {
    return new HDataEvent(HDataEventType.JSON_START_PROPERTY, name, null);
  }
  public static HDataEvent createJsonEndProperty() {
    return new HDataEvent(HDataEventType.JSON_END_PROPERTY, null, null);
  }
  public static HDataEvent createJsonValue(Object value) {
    return new HDataEvent(HDataEventType.JSON_VALUE, null, value);
  }
  public static HDataEvent createJsonError() {
    return new HDataEvent(HDataEventType.JSON_ERROR, null, null);
  }
}
