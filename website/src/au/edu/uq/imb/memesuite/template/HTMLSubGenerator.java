package au.edu.uq.imb.memesuite.template;

import java.util.Iterator;

/**
 * Simple helper designed to make generating HTMLSub objects from a collection easier.
 */
public abstract class HTMLSubGenerator<E> implements Iterable<HTMLSub> {
  private Iterable<E> items;

  public HTMLSubGenerator(Iterable<E> items) {
    this.items = items;
  }

  @Override
  public Iterator<HTMLSub> iterator() {
    return new HTMLSubGeneratorIterator(items);
  }

  protected abstract HTMLSub transform(E item);

  private class HTMLSubGeneratorIterator implements Iterator<HTMLSub> {
    protected Iterator<E> itemIter;

    public HTMLSubGeneratorIterator(Iterable<E> items) {
      this.itemIter = items.iterator();
    }

    @Override
    public boolean hasNext() {
      return itemIter.hasNext();
    }

    @Override
    public HTMLSub next() {
      E item = itemIter.next();
      return transform(item);
    }

    @Override
    public void remove() {
      throw new UnsupportedOperationException("Remove is not supported by the HTMLSubGenerator");
    }
  }
}
