package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.PsmDataSource;
import au.edu.uq.imb.memesuite.data.SequenceDataSource;
import au.edu.uq.imb.memesuite.data.SequenceInfo;
import au.edu.uq.imb.memesuite.db.SequenceDB;
import au.edu.uq.imb.memesuite.servlet.util.*;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import javax.activation.DataSource;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.Part;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.*;

public class Momo extends SubmitJob<Momo.Data> {
  private HTMLTemplate tmplMain;
  private HTMLTemplate tmplVerify;
  private ComponentHeader header;
  private ComponentSequences flanks;
  private ComponentJobDetails jobDetails;
  private ComponentAdvancedOptions advancedOptions;
  private ComponentSubmitReset submitReset;
  private ComponentFooter footer;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.Momo");

  protected class Data extends SubmitJob.JobData {
    public String email;
    public String description;
    public List<PsmDataSource> psms;
    public SequenceInfo flank;
    public Integer eliminateWidth;
    public Integer occurs;
    public Integer width;
    public String filterField;
    public String filterType;
    public Double filterThresh;
//    public Integer maxw;
    public boolean singlePerMass;
    
    public Integer countThreshold;
    public Double scoreThreshold;
    public Integer nmotifs;
    public String algorithm;
    public String fgFiletype;
    public String psmColumnName;

    public Data() {
      psms = new ArrayList<PsmDataSource>();
    }

    @Override
    public String email() {
      return email;
    }

    @Override
    public String description() {
      return description;
    }

    @Override
    public String emailTemplate() {
      return tmplVerify.getSubtemplate("message").toString();
    }

    @Override
    public List<DataSource> files() {
      ArrayList<DataSource> list = new ArrayList<DataSource>();
      if (flank != null && flank instanceof SequenceDataSource) {
        list.add(((SequenceDataSource) flank));
      }
      list.addAll(psms);
      return list;
    }

    @Override
    public String cmd() {
      StringBuilder args = new StringBuilder();
      if (singlePerMass) addArgs(args, "-single_per_mass");
      if (width != null) addArgs(args, "-width", width);
      if (eliminateWidth != null) addArgs(args, "-eliminate_width", eliminateWidth);
      if (filterField != null && filterType != null && filterThresh != null) {
        addArgs(args, "-filter_field", escapeOpal(filterField));
        addArgs(args, "-filter_type", filterType);
        addArgs(args, "-filter_thresh", filterThresh);
      }
      
      if (fgFiletype != null) addArgs(args, "-fg_filetype", fgFiletype);
      if (fgFiletype.equals("psm")) {
        if (psmColumnName != null) addArgs(args, "-psm_column_name", escapeOpal(psmColumnName));
      }
      
      if (flank != null) {
        if (flank instanceof SequenceDataSource) {
          addArgs(args, "-flank", ((SequenceDataSource) flank).getName());
        } else if (flank instanceof SequenceDB) {
          addArgs(args, "-flank", "db/" + ((SequenceDB) flank).getSequenceName());
        } else {
          throw new Error("Unexpected sequence type!");
        }
      }
      for (PsmDataSource psm : psms) {
        addArgs(args, "-psm", psm.getName());
      }
      
      if (algorithm.equals("alg_simp")) {
        if (occurs != null) addArgs(args, "-occurs", occurs);
      } else if (algorithm.equals("alg_mtfx")) {
        if (countThreshold != null) addArgs(args, "-count_threshold", countThreshold);
        if (scoreThreshold != null) addArgs(args, "-score_threshold", scoreThreshold);
      }
      if (algorithm != null) addArgs(args, "-algorithm", algorithm);
      
      return args.toString();
    }

    @Override
    public void cleanUp() {
      if (flank != null && flank instanceof SequenceDataSource) {
        if (!((SequenceDataSource) flank).getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              ((SequenceDataSource) flank).getFile());
        }
      }
      for (PsmDataSource psm : psms) {
        if (!psm.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " + psm.getFile());
        }
      }
    }

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      out.propertyAny("psms", psms);
      if (filterField != null && filterType != null && filterThresh != null) {
        out.property("filter_field", filterField);
        out.property("filter_type", filterType);
        out.property("filter_thresh", filterThresh);
      }
      
      if (fgFiletype != null) out.property("fg_filetype", fgFiletype);
      if (fgFiletype == "psm") {
        if (psmColumnName != null) out.property("psm_column_name", psmColumnName);
      }
      if (eliminateWidth != null) out.property("eliminate_width", eliminateWidth);
      if (width != null) out.property("width", width);
      
      out.property("single_per_mass", singlePerMass);
      
      if (algorithm == "alg_simp") {
        if (occurs != null) out.property("occurs", occurs);
      } else if (algorithm == "alg_mtfx") {
        if (countThreshold != null) out.property("count_threshold", countThreshold);
        if (scoreThreshold != null) out.property("score_threshold", scoreThreshold);
      }
      if (algorithm != null) out.property("algorithm", algorithm);
      
      out.endObject();
    }
  }

  public Momo() {
    super("MOMO", "MOMO");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the template
    this.tmplMain = cache.loadAndCache("/WEB-INF/templates/momo.tmpl");
    this.tmplVerify = cache.loadAndCache("/WEB-INF/templates/momo_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    flanks = new ComponentSequences(context, tmplMain.getSubtemplate("sequences"));
    jobDetails = new ComponentJobDetails(cache);
    advancedOptions = new ComponentAdvancedOptions(cache);
    submitReset = new ComponentSubmitReset(cache, jobTable.getCount(), jobTable.getDuration());
    footer = new ComponentFooter(cache, msp);
  }

  @Override
  protected void displayForm(HttpServletRequest request, HttpServletResponse response, long quotaMinWait) throws IOException {
    HTMLSub main = this.tmplMain.toSub();
    main.set("help", new HTMLSub[]{header.getHelp(),flanks.getHelp(),
        jobDetails.getHelp(), advancedOptions.getHelp(),
        submitReset.getHelp(), footer.getHelp()});
    main.set("header", header.getComponent());
    main.set("sequences", flanks.getComponent());
    main.set("job_details", jobDetails.getComponent());
    main.set("advanced_options", advancedOptions.getComponent());
    main.set("submit_reset", submitReset.getComponent(quotaMinWait));
    main.set("footer", footer.getComponent());
    response.setContentType("text/html; charset=UTF-8");
    main.output(response.getWriter());
  }

  @Override
  protected Data checkParameters(FeedbackHandler feedback, HttpServletRequest request) throws IOException, ServletException {
    // setup default file names
    FileCoord namer = new FileCoord();
    namer.createName("description");
    namer.createName("uuid");
    List<FileCoord.Name> psmNames = new ArrayList<FileCoord.Name>();
    for (int i = 1; request.getPart("psm_" + i) != null; i++) {
      psmNames.add(namer.createName("peptide_spectrum_matches_" + i + ".psm"));
    }
    FileCoord.Name flankName = namer.createName("match_flanks.faa");

    Data data = new Data();
    boolean error = true;
    try {
      // get the email
      data.email = jobDetails.getEmail(request, feedback);
      // get the description
      data.description = jobDetails.getDescription(request);
      // get the peptide-spectrum matches and their filters
      for (int i = 1; i <= psmNames.size(); i++) {
        FileCoord.Name psmName = psmNames.get(i-1);
        Part filePart = request.getPart("psm_" + i);
        if (filePart == null) throw new ServletException("Part psm_" + i + " is missing!");
        String fileName = getPartFilename(filePart);
        psmName.setOriginalName(fileName);
        if (filePart.getSize() == 0) {
          feedback.whine("PSM &quot;" + escapeForXML(fileName) + "&quot; is empty!");
          continue;
        }
        data.psms.add(new PsmDataSource(getPartTempFile(filePart), psmName));
      }
      
      if (paramBool(request, "filter_enable")) {
        String psmFilterField = null; String psmFilterType = null; Double psmFilterThreshold = null;
        psmFilterField = paramRequire(request, "filter_field");
        psmFilterType = paramChoice(request, "filter_type", "lt", "le", "eq", "ge", "gt");
        psmFilterThreshold = paramNumber(feedback, "filter threshold", request, "filter_thresh", null, null, 0.0);
        data.filterField = psmFilterField;
        data.filterType = psmFilterType;
        data.filterThresh = psmFilterThreshold;
      }

      data.fgFiletype = paramChoice(request, "fg_filetype", "psm", "prealigned", "fasta");
      data.psmColumnName = paramRequire(request, "psm_column_name");
      
      // get the flanking sequences for the PSMs
      data.flank = flanks.getSequences(flankName, request, feedback);
      // check if the flanking sequences need running through purge
      if (paramBool(request, "eliminate_enable")) {
        data.eliminateWidth = paramInteger(feedback, "eliminate repeats width", request, "eliminate_width", 0, null, 7);
      }
      // get the minimum number of occurrences needed to turn something into a motif
      data.occurs = paramInteger(feedback, "occurrences", request, "occurs", 1, null, 5);
      // get the minimum motif width
      data.width = paramInteger(feedback, "minimum width",
          request, "width", 2, 300, 7);
      data.singlePerMass = paramBool(request, "single_per_mass");
      
      data.countThreshold = paramInteger(feedback, "count_threshold", request, "count_threshold", 1, null, 20);
      data.scoreThreshold = paramNumber(feedback, "score_threshold", request, "score_threshold", 0.0, 1.0, 0.000001);
			
      data.algorithm = paramChoice(request, "algorithm", "alg_simp", "alg_mtfx", "alg_modl");
      
      error = false;
    } finally {
      if (error) data.cleanUp();
    }
    return data;
  }

  @Override
  public String title() {
    return tmplVerify.getSubtemplate("title").toString();
  }

  @Override
  public String subtitle() {
    return tmplVerify.getSubtemplate("subtitle").toString();
  }

  @Override
  public String logoPath() {
    return tmplVerify.getSubtemplate("logo").toString();
  }

  @Override
  public String logoAltText() {
    return tmplVerify.getSubtemplate("alt").toString();
  }
}
