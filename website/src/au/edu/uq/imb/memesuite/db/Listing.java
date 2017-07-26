package au.edu.uq.imb.memesuite.db;


import au.edu.uq.imb.memesuite.data.AlphStd;

import java.util.EnumSet;

public class Listing {
  private long id;
  private String name;
  private String description;
  private EnumSet<AlphStd> alphabets;
  private boolean priors;

  public Listing(long id, String name, String description, EnumSet<AlphStd> alphabets, boolean priors) {
    this.id = id;
    this.name = name;
    this.description = description;
    this.alphabets = alphabets;
    this.priors = priors;
  }

  public long getID() {
    return this.id;  
  }

  public String getName() {
    return this.name;
  }

  public String getDescription() {
    return this.description;
  }

  public EnumSet<AlphStd> getAlphabets() {
    return this.alphabets;
  }

  public boolean hasPriors() {
    return priors;
  }
}
