package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.*;
import au.edu.uq.imb.memesuite.db.SequenceDB;
import au.edu.uq.imb.memesuite.db.SequencePrior;
import au.edu.uq.imb.memesuite.servlet.util.*;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.activation.DataSource;
import javax.servlet.*;
import javax.servlet.http.*;

public class Fimo extends SubmitJob<Fimo.Data> {
  private HTMLTemplate tmplMain;
  private HTMLTemplate tmplVerify;
  private ComponentHeader header;
  private ComponentSequenceAlphabet alphabet;
  private ComponentMotifs motifs;
  private ComponentSequences sequences;
  private ComponentJobDetails jobDetails;
  private ComponentAdvancedOptions advBtn;
  private ComponentSubmitReset submitReset;
  private ComponentFooter footer;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.fimo");
  
  protected class Data extends SubmitJob.JobData {
    public String email;
    public String description;
    public AlphabetDataSource alphabet;
    public MotifDataSource motifs;
    public SequenceInfo sequences;
    public SequencePrior priors;
    public double pthresh;
    public boolean norc;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      if (alphabet != null) out.property("alphabet", alphabet);
      out.property("motifs", motifs);
      out.property("sequences", sequences);
      out.property("priors", priors);
      out.property("pthresh", pthresh);
      out.property("norc", norc);
      out.endObject();
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
    public String cmd() {
  //    fimo_webservice [options] <motifs> <db seqs>
  //
  //      Options:
  //        -upseqs <file>    uploaded sequences
  //        -pvthresh <pv>    output p-value threshold
  //        -norc             scan given strand only
  //        -help             brief help message
      StringBuilder args = new StringBuilder();
      if (sequences instanceof SequenceDataSource) {
        addArgs(args, "-upseqs", ((SequenceDataSource) sequences).getName());
      }
      if (priors != null) {
        addArgs(args, "-psp", priors.getPriorsName());
        addArgs(args, "-prior-dist", priors.getPriorsDistName());
      }
      addArgs(args, "-pvthresh", pthresh);
      if (norc) addArgs(args, "-norc");
      addArgs(args, motifs.getName());
      if (sequences instanceof SequenceDB) {
        addArgs(args, ((SequenceDB) sequences).getSequenceName());
      }
      return args.toString();
    }
  
    @Override
    public List<DataSource> files() {
      ArrayList<DataSource> list = new ArrayList<DataSource>();
      if (alphabet != null) list.add(alphabet);
      if (motifs != null) list.add(motifs);
      if (sequences != null && sequences instanceof SequenceDataSource) {
        list.add((SequenceDataSource) sequences);
      }
      return list;
    }
  
    @Override
    public void cleanUp() {
      if (alphabet != null) {
        if (!alphabet.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              alphabet.getFile());
        }
      }
      if (sequences != null && sequences instanceof SequenceDataSource) {
        if (!((SequenceDataSource) sequences).getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              ((SequenceDataSource) sequences).getFile());
        }
      }
      if (motifs != null) {
        if (!motifs.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              motifs.getFile());
        }
      }
    }
  }

  public Fimo() {
    super("FIMO", "FIMO");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the templates
    tmplMain = cache.loadAndCache("/WEB-INF/templates/fimo.tmpl");
    tmplVerify = cache.loadAndCache("/WEB-INF/templates/fimo_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    alphabet = new ComponentSequenceAlphabet(context, tmplMain.getSubtemplate("alphabet"));
    motifs = new ComponentMotifs(context, tmplMain.getSubtemplate("motifs"));
    sequences = new ComponentSequences(context, tmplMain.getSubtemplate("sequences"));
    jobDetails = new ComponentJobDetails(cache);
    advBtn = new ComponentAdvancedOptions(cache);
    submitReset = new ComponentSubmitReset(cache, jobTable.getCount(), jobTable.getDuration());
    footer = new ComponentFooter(cache, msp);
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

  @Override
  protected void displayForm(HttpServletRequest request, HttpServletResponse response, long quotaMinWait) throws IOException {
    HTMLSub main = tmplMain.toSub();
    main.set("help", new HTMLSub[]{header.getHelp(), alphabet.getHelp(), motifs.getHelp(),
        sequences.getHelp(), jobDetails.getHelp(), advBtn.getHelp(),
        submitReset.getHelp(), footer.getHelp()});
    main.set("header", header.getComponent());
    main.set("alphabet", alphabet.getComponent());
    main.set("motifs", motifs.getComponent(request.getParameter("motifs_embed")));
    main.set("sequences", sequences.getComponent());
    main.set("job_details", jobDetails.getComponent());
    main.set("advanced_options", advBtn.getComponent());
    main.set("submit_reset", submitReset.getComponent(quotaMinWait));
    main.set("footer", footer.getComponent());
    response.setContentType("text/html; charset=UTF-8");
    main.output(response.getWriter());
  }

  @Override
  protected Data checkParameters(FeedbackHandler feedback,
      HttpServletRequest request) throws IOException, ServletException {
    // setup default file names
    FileCoord namer = new FileCoord();
    FileCoord.Name sequencesName = namer.createName("sequences.fa");
    FileCoord.Name motifsName = namer.createName("motifs.meme");
    namer.createName("description");
    namer.createName("uuid");
    Alph alph = null;
    Data data = new Data();
    // get the job details
    data.email = jobDetails.getEmail(request, feedback);
    data.description = jobDetails.getDescription(request);
    // get the motifs
    data.motifs =  (MotifDataSource)motifs.getMotifs(motifsName, request, feedback);
    if (data.motifs != null) alph = data.motifs.getAlphabet();
    // get the sequences
    data.sequences = sequences.getSequences(alph, sequencesName, request, feedback);
    data.priors = sequences.getPrior(request);
    // other options
    data.pthresh = WebUtils.paramNumber(feedback, "output p-value threshold",
        request, "output_pv", 0.0, 1.0, 1e-4);
    data.norc = WebUtils.paramBool(request, "norc");
    return data;
  }
}

