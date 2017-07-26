package au.edu.uq.imb.memesuite.db;

import au.edu.uq.imb.memesuite.data.AlphStd;

import java.util.*;

/**
 * Object to store information related to a sequence database version.
 */
public class SequenceVersion {
  private EnumMap<AlphStd,SequenceDB> sequenceFileMap;
  private long edition;
  private String version;

  public SequenceVersion(long edition, String version, Collection<SequenceDB> sequenceDBs) {
    this.sequenceFileMap = new EnumMap<AlphStd, SequenceDB>(AlphStd.class);
    this.edition = edition;
    this.version = version;
    for (SequenceDB file : sequenceDBs) {
      sequenceFileMap.put(file.guessAlphabet(), file);
    }
  }

  public String getListingName() {
    return sequenceFileMap.values().iterator().next().getListingName();
  }

  public String getListingDescription() {
    return sequenceFileMap.values().iterator().next().getListingDescription();
  }

  public Set<AlphStd> getAlphabets() {
    return sequenceFileMap.keySet();
  }

  public long getEdition() {
    return edition;
  }

  public String getVersion() {
    return version;
  }

  public SequenceDB getSequenceFile(AlphStd alphabet) {
    return sequenceFileMap.get(alphabet);
  }


}
