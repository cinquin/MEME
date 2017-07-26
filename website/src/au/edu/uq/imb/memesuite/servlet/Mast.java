package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.*;
import au.edu.uq.imb.memesuite.db.SequenceDB;
import au.edu.uq.imb.memesuite.servlet.util.*;
import au.edu.uq.imb.memesuite.template.*;
import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.activation.DataSource;
import javax.servlet.*;
import javax.servlet.http.*;

import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.*;

public class Mast extends SubmitJob<Mast.Data> {
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

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.mast");
  
  protected class Data extends SubmitJob.JobData {
    public String email;
    public String description;
    public AlphabetDataSource alphabet;
    public MotifDataSource motifs;
    public SequenceInfo sequences;
    public String strands;
    public boolean translateDNA;
    public boolean useSeqComp;
    public boolean scaleM;
    public boolean remCorr;
    public double seqEvalue;
    public Double motifEvalue;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      if (alphabet != null) out.property("alphabet", alphabet);
      out.property("motifs", motifs);
      out.property("sequences", sequences);
      out.property("translate_dna", translateDNA);
      out.property("strands", strands);
      out.property("use_seq_comp", useSeqComp);
      out.property("scale_m", scaleM);
      out.property("rem_corr", remCorr);
      out.property("seq_evalue", seqEvalue);
      out.property("motif_evalue", motifEvalue);
      out.endObject();
    }


    @Override
    public String email() {
      return email;
    }
  
    @Override
    public String description() {
      return description;  // generated code
    }

    @Override
    public String emailTemplate() {
      return tmplVerify.getSubtemplate("message").toString();
    }
  
    @Override
    public String cmd() {
  //    mast_webservice [options] <motifs file> [<sequence database>]
  //
  //      Options:
  //        -dna              search nucleotide database with protein motifs
  //        -seqp             Scale motif display threshold by sequence length
  //        -comp             Use individual sequence composition in E-value and p-value calculation
  //        -sep              treat the rc strand as a seperate sequence; not compatible with -norc
  //        -norc             don't process the rc strand; not compatible with -sep
  //        -remcorr          remove redundant motifs
  //        -ev <thresh>      display sequences with evalue below this threshold
  //        -mev <thresh>     Ignore motifs with evalue above this threshold
  //        -upload_db <file> uploaded sequence database
  //        -df <name>        Name to show for sequence database; underscores are converted to spaces
  //        -dl <dl>          in results use <dl> as link to search sequence
  //			      names; ignored when -dblist specified
  //        -help             brief help message
  
  
      StringBuilder args = new StringBuilder();
      if (translateDNA) addArgs(args, "-dna");
      if (scaleM) addArgs(args, "-seqp");
      if (remCorr) addArgs(args, "-remcorr");
      if (useSeqComp) addArgs(args, "-comp");
      if (strands.equals("separate")) addArgs(args, "-sep");
      if (strands.equals("ignore")) addArgs(args, "-norc");
      addArgs(args, "-ev", seqEvalue);
      if (motifEvalue != null) addArgs(args, "-mev", motifEvalue);
      addArgs(args, "-nseqs", sequences.getSequenceCount());
      // Add URL for wrapping protein sequence IDs in output.
      Alph alph = motifs.getAlphabet();
      if (alph != null && alph.getLikeAlphabet().toString() == "PROTEIN") {
        addArgs(args, "-dl", escapeOpal("http://www.uniprot.org/uniprot/?query=SEQUENCEID&sort=score"));
        strands = null;		// Suppress strands in status message.
      }
      if (sequences instanceof SequenceDB) {
        addArgs(args, "-df", escapeOpal(((SequenceDB) sequences).getSequenceName()));
      } else {
        addArgs(args, "-upload_db", ((SequenceDataSource) sequences).getName());
      }
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

  public Mast() {
    super("MAST", "MAST");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the templates
    tmplMain = cache.loadAndCache("/WEB-INF/templates/mast.tmpl");
    tmplVerify = cache.loadAndCache("/WEB-INF/templates/mast_verify.tmpl");
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
    FileCoord.Name alphName = namer.createName("alphabet.alph");
    FileCoord.Name sequencesName = namer.createName("sequences.fa");
    FileCoord.Name motifsName = namer.createName("motifs.meme");
    namer.createName("description");
    namer.createName("uuid");
    Alph alph = null;
    Data data = new Data();
    // translate dna
    data.translateDNA = WebUtils.paramBool(request, "translate_dna");
    // get the job details
    data.email = jobDetails.getEmail(request, feedback);
    data.description = jobDetails.getDescription(request);
    // Get the alphabet
    data.alphabet = alphabet.getAlphabet(alphName, request, feedback);
    if (data.alphabet != null) alph = data.alphabet.getAlph();
    // get the motifs
    data.motifs = (MotifDataSource)motifs.getMotifs(motifsName, request, feedback);
    if (data.motifs != null) {
      if (alph == null) alph = data.motifs.getAlphabet();
      if (data.translateDNA) {
        alph = AlphStd.PROTEIN.getAlph();
      }
    }
    // get the sequences
    data.sequences = sequences.getSequences(alph, sequencesName, request, feedback);
    // Options
    // strand handling
    data.strands = WebUtils.paramChoice(request, "strands", "combine", "separate", "ignore");
    // sequence E-value
    data.seqEvalue = WebUtils.paramNumber(feedback, "sequence E-value", request, "seq_evalue", 0.0, null, 10.0);
    // use sequence composition
    data.useSeqComp = WebUtils.paramBool(request, "use_seq_comp");
    // scale motif display threshold
    data.scaleM = WebUtils.paramBool(request, "scale_m");
    // remove redundant motifs 
    data.remCorr = WebUtils.paramBool(request, "rem_corr");
    // filter motifs by E-value
    if (WebUtils.paramBool(request, "motif_evalue_enable")) {
      data.motifEvalue = WebUtils.paramNumber(feedback, "motif E-value", request, "motif_evalue", 0.0, null, 0.05);
    } else {
      data.motifEvalue = null;
    }
    return data;
  }
}
