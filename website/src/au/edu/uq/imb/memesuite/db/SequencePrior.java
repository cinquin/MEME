package au.edu.uq.imb.memesuite.db;

import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.IOException;
import java.net.URL;

/**
 *
 */
public class SequencePrior implements JsonWr.JsonValue {
  private long id;
  private long sequenceId;
  private String filePrior;
  private String fileDist;
  private String biosample;
  private String assay;
  private String source;
  private URL url;
  private String description;

  public SequencePrior(long id, long sequenceId, String filePrior, String fileDist, String biosample, String assay, String source, URL url, String description) {
    this.id = id;
    this.sequenceId = sequenceId;
    this.filePrior = filePrior;
    this.fileDist = fileDist;
    this.biosample = biosample;
    this.assay = assay;
    this.source = source;
    this.url = url;
    this.description = description;
  }

  public long getId() {
    return id;
  }

  public long getSequenceId() {
    return sequenceId;
  }

  public String getPriorsName() {
    return filePrior;
  }

  public String getPriorsDistName() {
    return fileDist;
  }

  public String getBiosample() {
    return biosample;
  }

  public String getAssay() {
    return assay;
  }

  public String getSource() {
    return source;
  }

  public URL getUrl() {
    return url;
  }

  public String getDescription() {
    return description;
  }

  public String toString() {
    return biosample + "; " + assay + "; " + source;
  }

  @Override
  public void outputJson(JsonWr out) throws IOException {
    out.startObject();
    out.property("id", id);
    out.property("biosample", biosample);
    out.property("assay", assay);
    out.property("source", source);
    out.property("url", url.toString());
    out.property("description", description);
    out.endObject();
  }
}
