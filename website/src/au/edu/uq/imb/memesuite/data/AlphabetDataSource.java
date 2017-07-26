package au.edu.uq.imb.memesuite.data;


import au.edu.uq.imb.memesuite.util.FileCoord;

import java.io.File;

/**
 *
 */
public class AlphabetDataSource extends NamedFileDataSource {
  private Alph alph;
  public AlphabetDataSource(File file, FileCoord.Name name, Alph alph) {
    super(file, name);
    this.alph = alph;
  }
  public Alph getAlph() {
    return alph;
  }
}
