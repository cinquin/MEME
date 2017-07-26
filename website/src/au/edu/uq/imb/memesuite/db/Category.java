package au.edu.uq.imb.memesuite.db;

public class Category {
  private long id;
  private String name;
  private long count;
  private boolean priors;

  public Category(long id, String name, long count, boolean priors) {
    this.id = id;
    this.name = name;
    this.count = count;
    this.priors = priors;
  }

  public long getID() {
    return this.id;
  }

  public String getName() {
    return this.name;
  }

  public long getCount() {
    return this.count;
  }

  public boolean hasPriors() {
    return this.priors;
  }
}
