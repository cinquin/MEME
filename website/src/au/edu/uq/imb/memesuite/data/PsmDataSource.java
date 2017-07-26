package au.edu.uq.imb.memesuite.data;

import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.File;
import java.io.IOException;

public class PsmDataSource extends NamedFileDataSource {

  public PsmDataSource(File file, FileCoord.Name name) {
    super(file, name);
  }
  
}
