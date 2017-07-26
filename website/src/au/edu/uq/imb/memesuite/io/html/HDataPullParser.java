package au.edu.uq.imb.memesuite.io.html;

import java.io.IOException;
import java.io.Reader;
import java.nio.CharBuffer;
import java.util.*;

/**
 *
 */
public class HDataPullParser {
  protected Deque<HDataEvent> events;
  protected Reader reader;
  protected CharBuffer buffer;
  protected HDataParser parser;
  protected HDataEvent current;

  protected class HDataEventGenerator implements HDataHandler {
    @Override
    public void htmlHiddenData(String name, String value) {
      events.add(HDataEvent.createHiddenField(name, value));
    }
    @Override
    public void jsonStartData(String name) {
      events.add(HDataEvent.createJsonStartData(name));
    }
    @Override
    public void jsonEndData() {
      events.add(HDataEvent.createJsonEndData());
    }
    @Override
    public void jsonStartObject() {
      events.add(HDataEvent.createJsonStartObject());
    }
    @Override
    public void jsonEndObject() {
      events.add(HDataEvent.createJsonEndObject());
    }
    @Override
    public void jsonStartList() {
      events.add(HDataEvent.createJsonStartList());
    }
    @Override
    public void jsonEndList() {
      events.add(HDataEvent.createJsonEndList());
    }
    @Override
    public void jsonStartProperty(String name) {
      events.add(HDataEvent.createJsonStartProperty(name));
    }
    @Override
    public void jsonEndProperty() {
      events.add(HDataEvent.createJsonEndProperty());
    }
    @Override
    public void jsonValue(Object value) {
      events.add(HDataEvent.createJsonValue(value));
    }
    @Override
    public void jsonError() {
      events.add(HDataEvent.createJsonError());
    }
  }

  public HDataPullParser(Reader in) {
    events = new ArrayDeque<HDataEvent>();
    reader = in;
    parser = new HDataParser();
    parser.setHandler(new HDataEventGenerator());
    current = HDataEvent.createFileStart();
    buffer = CharBuffer.allocate(100);
  }

  public void close() throws IOException {
    if (reader != null) {
      reader.close();
      reader = null;
    }
  }

  public HDataEventType getEventType() {
    return (current != null ? current.getType() : HDataEventType.FILE_END);
  }

  public String getName() {
    if (current == null) throw new NullPointerException("At end of file!");
    return current.getName();
  }

  public Object getValue() {
    if (current == null) throw new NullPointerException("At end of file!");
    return current.getValue();
  }

  public boolean isProperty(String name) {
    return getEventType() == HDataEventType.JSON_START_PROPERTY && !(name != null && !name.equals(getName()));
  }

  public void expect(HDataEventType type, String name, Class valueClass) throws HDataParseException {
    if (getEventType() != type) throw new HDataParseException("Expected state to be " + type.name() + " but actual state was " + getEventType().name());
    if (name != null && !name.equals(getName())) throw new HDataParseException("Expected " + type.name() + " to have name " + name + " however actual name was " + getName());
    if (valueClass != null && !valueClass.isInstance(getValue())) throw new HDataParseException("Expected " + type.name() +
        " to have a value that is an instance of class " + valueClass.getCanonicalName());
  }

  public void expect(String name) throws HDataParseException {
    expect(HDataEventType.JSON_START_PROPERTY, name, null);
  }

  public void expect(HDataEventType type, Class valueClass) throws HDataParseException {
    expect(type, null, valueClass);
  }

  public void expect(HDataEventType type) throws HDataParseException {
    expect(type, null, null);
  }

  public HDataEventType next() throws IOException {
    if (!events.isEmpty()) {
      current = events.remove();
      return current.getType();
    }
    if (reader == null) return HDataEventType.FILE_END;
    boolean error = true;
    try {
      while (events.isEmpty() && reader.read(buffer) != -1) {
        buffer.flip();
        parser.process(buffer);
        buffer.compact();
      }
      error = false;
    } finally {
      if (error) {
        try {
          reader.close();
        } catch (IOException e) { /* ignore */ }
        reader = null;
      }
    }
    current = events.poll();
    return (current != null ? current.getType() : HDataEventType.FILE_END);
  }

  public HDataEventType next(HDataEventType type, String name, Class valueClass) throws HDataParseException, IOException {
    expect(type, name, valueClass);
    return next();
  }

  public HDataEventType next(String name) throws HDataParseException, IOException {
    expect(name);
    return next();
  }

  public HDataEventType next(HDataEventType type, Class valueClass) throws HDataParseException, IOException {
    expect(type, valueClass);
    return next();
  }

  public HDataEventType next(HDataEventType type) throws HDataParseException, IOException {
    expect(type);
    return next();
  }

  public void skipObject() throws IOException, HDataParseException {
    next(HDataEventType.JSON_START_OBJECT);
    while (getEventType() == HDataEventType.JSON_START_PROPERTY) skipProperty();
    next(HDataEventType.JSON_END_OBJECT);
  }

  public Map<String,Object> nextObject() throws IOException, HDataParseException {
    Map<String,Object> properties = new TreeMap<String, Object>();
    next(HDataEventType.JSON_START_OBJECT);
    while (getEventType() == HDataEventType.JSON_START_PROPERTY) {
      String name = getName();
      Object value = nextProperty();
      properties.put(name, value);
    }
    next(HDataEventType.JSON_END_OBJECT);
    return properties;
  }

  public void skipList() throws IOException, HDataParseException {
    next(HDataEventType.JSON_START_LIST);
    while (getEventType() != HDataEventType.JSON_END_LIST) {
      switch (getEventType()) {
        case JSON_VALUE:
          next();
          break;
        case JSON_START_OBJECT:
          skipObject();
          break;
        case JSON_START_LIST:
          skipList();
          break;
        default:
          throw new HDataParseException("Expected state to be either JSON_VALUE, JSON_START_OBJECT, or JSON_START_LIST");
      }
    }
    next(HDataEventType.JSON_END_LIST);
  }

  public List<Object> nextList() throws IOException, HDataParseException {
    List<Object> values = new ArrayList<Object>();
    next(HDataEventType.JSON_START_LIST);
    while (getEventType() != HDataEventType.JSON_END_LIST) {
      switch (getEventType()) {
        case JSON_VALUE:
          values.add(getValue());
          next();
          break;
        case JSON_START_OBJECT:
          values.add(nextObject());
          break;
        case JSON_START_LIST:
          values.add(nextList());
          break;
        default:
          throw new HDataParseException("Expected state to be either JSON_VALUE, JSON_START_OBJECT, or JSON_START_LIST");
      }
    }
    next(HDataEventType.JSON_END_LIST);
    return values;
  }

  public void skipProperty() throws IOException, HDataParseException {
    expect(HDataEventType.JSON_START_PROPERTY);
    switch (next()) {
      case JSON_VALUE:
        next();
        break;
      case JSON_START_OBJECT:
        skipObject();
        break;
      case JSON_START_LIST:
        skipList();
        break;
      default:
        throw new HDataParseException("Expected state to be either JSON_VALUE, JSON_START_OBJECT, or JSON_START_LIST");
    }
    next(HDataEventType.JSON_END_PROPERTY);
  }

  public void skipProperty(String name) throws IOException, HDataParseException {
    expect(name);
    skipProperty();
  }

  public Object nextProperty() throws IOException, HDataParseException {
    expect(HDataEventType.JSON_START_PROPERTY);
    Object value;
    switch (next()) {
      case JSON_VALUE:
        value =  getValue();
        next();
        break;
      case JSON_START_OBJECT:
        value = nextObject();
        break;
      case JSON_START_LIST:
        value = nextList();
        break;
      default:
        throw new HDataParseException("Expected state to be either JSON_VALUE, JSON_START_OBJECT, or JSON_START_LIST");
    }
    next(HDataEventType.JSON_END_PROPERTY);
    return value;
  }

  public Object nextProperty(String name) throws IOException, HDataParseException {
    expect(name);
    return nextProperty();
  }

  public <T> List<T> nextList(Class<T> listType) throws IOException, HDataParseException {
    List<T> list = new ArrayList<T>();
    next(HDataEventType.JSON_START_LIST);
    while (getEventType() != HDataEventType.JSON_END_LIST) {
      expect(HDataEventType.JSON_VALUE, listType);
      list.add(listType.cast(getValue()));
      next();
    }
    next(HDataEventType.JSON_END_LIST);
    return list;
  }

  public <T> List<List<T>> nextMatrix(Class<T> listType) throws IOException, HDataParseException {
    List<List<T>> matrix = new ArrayList<List<T>>();
    next(HDataEventType.JSON_START_LIST);
    while (getEventType() != HDataEventType.JSON_END_LIST) {
      matrix.add(nextList(listType));
    }
    next(HDataEventType.JSON_END_LIST);
    return matrix;
  }

  public <T> T nextProperty(String propertyName, Class<T> type) throws IOException, HDataParseException {
    next(propertyName);
    expect(HDataEventType.JSON_VALUE, type);
    Object value = getValue();
    next();
    next(HDataEventType.JSON_END_PROPERTY);
    return type.cast(value);
  }

  public String nextStringProperty(String propertyName) throws IOException, HDataParseException {
    return nextProperty(propertyName, String.class);
  }

  public long nextIntegerProperty(String propertyName) throws IOException, HDataParseException {
    return nextProperty(propertyName, Long.class);
  }

  public double nextNumberProperty(String propertyName) throws IOException, HDataParseException {
    return nextProperty(propertyName, Number.class).doubleValue();
  }

  public boolean nextBooleanProperty(String propertyName) throws IOException, HDataParseException {
    return nextProperty(propertyName, Boolean.class);
  }

  public <T> List<T> nextListProperty(String propertyName, Class<T> listType) throws IOException, HDataParseException {
    List<T> list;
    next(propertyName);
    list = nextList(listType);
    next(HDataEventType.JSON_END_PROPERTY);
    return list;
  }

  public <T> List<List<T>> nextMatrixProperty(String propertyName, Class<T> listType) throws IOException, HDataParseException {
    List<List<T>> matrix;
    next(propertyName);
    matrix = nextMatrix(listType);
    next(HDataEventType.JSON_END_PROPERTY);
    return matrix;
  }

  public static boolean hasEntry(Object object, Object ... keys) {
    for (Object key : keys) {
      if (key instanceof Integer) {
        int index = (Integer)key;
        if (index < 0) throw new IllegalArgumentException("Integer keys must be non-negative");
        if (object instanceof List) {
          List list = (List)object;
          if (list.size() < index) return false;
          object = list.get(index);
        } else {
          return false;
        }
      } else if (key instanceof String) {
        String property = (String)key;
        if (object instanceof Map) {
          Map map = (Map)object;
          if (!map.containsKey(property)) return false;
          object = map.get(property);
        } else {
          return false;
        }
      } else if (key instanceof Class) {
        Class cls = (Class)key;
        return cls.isInstance(object);
      } else {
        throw new IllegalArgumentException("Keys must be string, integer or class");
      }
    }
    return true;
  }

  public static Object getEntry(Object object, Object ... keys) throws HDataParseException {
    for (Object key : keys) {
      if (key instanceof Integer) {
        int index = (Integer)key;
        if (index < 0) throw new IllegalArgumentException("Integer keys must be non-negative");
        if (object instanceof List) {
          List list = (List)object;
          if (list.size() < index) throw new HDataParseException("List did not contain index " + index);
          object = list.get(index);
        } else {
          throw new HDataParseException("Not a list");
        }
      } else if (key instanceof String) {
        String property = (String)key;
        if (object instanceof Map) {
          Map map = (Map)object;
          if (!map.containsKey(property)) throw new HDataParseException("Map did not contain property " + property);
          object = map.get(property);
        } else {
          throw new HDataParseException("Not a map");
        }
      } else if (key instanceof Class) {
        Class cls = (Class)key;
        if (cls.isInstance(object)) {
          return object;
        } else {
          throw new HDataParseException("Value was not of class " + cls);
        }
      } else {
        throw new IllegalArgumentException("Keys must be string, integer or class");
      }
    }
    return object;
  }

  public static Object getOptEntry(Object object, Object defVal, Object ... keys) throws HDataParseException {
    return (hasEntry(object, keys) ? getEntry(object, keys) : defVal);
  }
}
