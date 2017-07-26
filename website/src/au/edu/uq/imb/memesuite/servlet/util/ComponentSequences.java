package au.edu.uq.imb.memesuite.servlet.util;

import au.edu.uq.imb.memesuite.data.*;
import au.edu.uq.imb.memesuite.db.*;
import au.edu.uq.imb.memesuite.io.fasta.FastaException;
import au.edu.uq.imb.memesuite.io.fasta.FastaParser;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLSubGenerator;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;
import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.Part;
import java.io.*;
import java.sql.SQLException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.*;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.SEQUENCE_DB_KEY;


/**
 * A component that is used for inputting FASTA sequences.
 */
public class ComponentSequences extends PageComponent {
  private ServletContext context;
  private HTMLTemplate tmplSequences;
  private AlphType type;
  private boolean shortOnly; // only short databases
  private String prefix;
  private String fieldName; // for feedback
  private HTMLTemplate title;
  private HTMLTemplate subtitle;
  private boolean loadDBs;
  private boolean loadPriors;
  private boolean enableNoSeq;
  private String registerFn;
  private DefaultOption defaultOption;
  // sequence options
  private boolean weights;
  private boolean mask;
  private boolean ambigs;
  private boolean gaps;
  private Integer maxNameLen;
  private Integer maxDescLen;
  private Integer minSeqLen;
  private Integer maxSeqLen;
  private Integer minSeqCount;
  private Integer maxSeqCount;
  private Integer maxSeqTotal;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.component.sequence");

  /**
   * Enum of the initial selection state.
   */
  public static enum DefaultOption {
    NOSEQ,
    TEXT,
    FILE,
    EMBED,
    DATABASE
  }

  private static class CategoryTemplate extends HTMLSubGenerator<Category> {
    private HTMLTemplate template;
    private boolean selectFirst;
    boolean selectable;
    boolean identifiable;
    boolean nameable;
    boolean first;

    public CategoryTemplate(DBList db, HTMLTemplate template, boolean selectFirst,
        boolean shortOnly, Set<AlphStd> allowedAlphabets) throws SQLException {
      super(db.getCategories(shortOnly, allowedAlphabets));
      this.template = template;
      this.selectFirst = selectFirst;
      selectable = template.containsSubtemplate("selected");
      identifiable = template.containsSubtemplate("id");
      nameable = template.containsSubtemplate("name");
      first = true;
    }

    @Override
    protected HTMLSub transform(Category item) {
      HTMLSub category = template.toSub();
      if (selectable && selectFirst && first) category.set("selected", "selected");
      if (identifiable) category.set("id", item.getID());
      if (nameable) category.set("name", item.getName());
      category.set("no_priors", item.hasPriors() ? "" : "no_priors");
      first = false;
      return category;
    }
  }

  public ComponentSequences(ServletContext context, HTMLTemplate info) throws ServletException {
    this.context = context;
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    tmplSequences = cache.loadAndCache("/WEB-INF/templates/component_sequences.tmpl");
    prefix = getText(info, "prefix", "sequences");
    fieldName = getText(info, "description", "sequences");
    registerFn = getText(info, "register", "nop");
    title = getTemplate(info, "title", null);
    subtitle = getTemplate(info, "subtitle", null);
    type = getEnum(info, "alph_type", AlphType.class, AlphType.ANY_ALPHABET);
    shortOnly = info.containsSubtemplate("short_only");
    loadDBs = info.containsSubtemplate("enable_db");
    loadPriors = loadDBs && info.containsSubtemplate("enable_priors");
    enableNoSeq = info.containsSubtemplate("enable_noseq");
    defaultOption = getEnum(info, "default", DefaultOption.class, DefaultOption.FILE);
    if (defaultOption == DefaultOption.NOSEQ) enableNoSeq = true;
    // sequence options
    mask = !info.containsSubtemplate("disable_masking");
    ambigs = !info.containsSubtemplate("disable_ambiguous");
    weights = info.containsSubtemplate("enable_weights");
    gaps = info.containsSubtemplate("enable_gaps");
    maxNameLen = getInt(info, "max_name_len", 50);
    maxDescLen = getInt(info, "max_desc_len", 500);
    minSeqLen = getInt(info, "min_seq_len", 1);
    maxSeqLen = getInt(info, "max_seq_len", null);
    minSeqCount = getInt(info, "min_seq_count", 1);
    maxSeqCount = getInt(info, "max_seq_count", null);
    maxSeqTotal = getInt(info, "max_seq_total", null);
  }

  public HTMLSub getComponent() {
    return getComponent(null, null);
  }

  public HTMLSub getComponent(String embedSeqs) {
    return getComponent(embedSeqs, null);
  }

  public HTMLSub getComponent(String embedSeqs, String embedName) {
    boolean selectFirstDB = false;
    DefaultOption defaultOption = this.defaultOption;
    HTMLSub sequences = tmplSequences.getSubtemplate("component").toSub();
    sequences.set("prefix", prefix);
    if (title != null) sequences.set("title", title);
    if (subtitle != null) sequences.set("subtitle", subtitle);
    if (!enableNoSeq) sequences.empty("noseq_option");
    if (embedSeqs != null) {
      sequences.getSub("embed_section").set("prefix", prefix).
          set("data", WebUtils.escapeForXML(embedSeqs));
      if (embedName != null) {
        sequences.getSub("embed_section").set("name", WebUtils.escapeForXML(embedName));
      }
      defaultOption = DefaultOption.EMBED;
    } else {
      sequences.empty("embed_option");
      sequences.empty("embed_section");
    }
    switch (defaultOption) {
      case NOSEQ:
        sequences.getSub("noseq_option").set("selected", "selected");
        break;
      case TEXT:
        sequences.getSub("text_option").set("selected", "selected");
        break;
      case FILE:
        sequences.getSub("file_option").set("selected", "selected");
        break;
      case EMBED:
        if (embedSeqs != null) sequences.getSub("embed_option").set("selected", "selected");
        break;
      case DATABASE:
        selectFirstDB = true;
    }
    SequenceDBList db = (SequenceDBList)context.getAttribute(SEQUENCE_DB_KEY);
    if (loadDBs && db != null) {
      try {
        HTMLTemplate cat_opt = tmplSequences.getSubtemplate("component").
            getSubtemplate("cat_options").getSubtemplate("cat_option");
        sequences.getSub("cat_options").set("cat_option",
            new CategoryTemplate(db, cat_opt, selectFirstDB, shortOnly, type.getStandardAlphabets()));
        sequences.getSub("db_section").set("prefix", prefix);
        if (loadPriors) {
          sequences.getSub("db_section").getSub("priors_section_2").set("prefix", prefix);
        } else {
          sequences.empty("priors_section_1");
          sequences.getSub("db_section").empty("priors_section_2");
        }
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error loading sequence categories", e);
        sequences.empty("priors_section_1");
        sequences.empty("cat_options");
        sequences.empty("db_section");
      }
    } else {
      sequences.empty("priors_section_1");
      sequences.empty("cat_options");
      sequences.empty("db_section");
    }

    StringWriter buf = new StringWriter();
    JsonWr jsonWr = new JsonWr(buf, 18);
    try {
      jsonWr.start();
      jsonWr.property("field", fieldName);
      jsonWr.property("alph_type", type.name());
      jsonWr.property("short_only", shortOnly);
      jsonWr.property("weights", weights);
      jsonWr.property("mask", mask);
      jsonWr.property("ambigs", ambigs);
      jsonWr.property("gaps", gaps);
      jsonWr.property("max_name_len", maxNameLen);
      jsonWr.property("max_desc_len", maxDescLen);
      jsonWr.property("min_seq_len", minSeqLen);
      jsonWr.property("max_seq_len", maxSeqLen);
      jsonWr.property("min_seq_count", minSeqCount);
      jsonWr.property("max_seq_count", maxSeqCount);
      jsonWr.property("max_seq_total", maxSeqTotal);
      jsonWr.end();
    } catch (IOException e) {
      // no IO exceptions should occur as this uses a StringBuffer
      throw new Error(e);
    }
    sequences.set("options", buf.toString());
    sequences.set("register_component", registerFn);
    return sequences;
  }

  public HTMLSub getHelp() {
    return tmplSequences.getSubtemplate("help").toSub();
  }

  private boolean checkSpec(Collection<Alph> restrictedAlphabets, SequenceStats stats, FeedbackHandler feedback) {
    boolean ok = true;
    if (minSeqCount != null && minSeqCount > stats.getSequenceCount()) {
      feedback.whine("Too few sequences for " + fieldName);
      ok = false;
    } else if (maxSeqCount != null && maxSeqCount < stats.getSequenceCount()) {
      feedback.whine("Too many sequences for " + fieldName);
      ok = false;
    }
    if (minSeqLen != null && minSeqLen > stats.getMinLength()) {
      feedback.whine(
          "There are one or more sequences that are too short in the "
          + fieldName);
      ok = false;
    }
    if (maxSeqLen != null && maxSeqLen < stats.getMaxLength()) {
      feedback.whine(
          "There are one or more " + fieldName + " that are too long");
      ok = false;
    }
    if (maxSeqTotal != null && maxSeqTotal < stats.getTotalLength()) {
      feedback.whine(
          "The combined length of all the " + fieldName + " is to long.");
      ok = false;
    }
    if (!stats.checkAlphabets(type, restrictedAlphabets)) {
      feedback.whine(
          "The alphabet of the " + fieldName + " seems to be " +
              stats.guessAlphabet() + " but it is not one of the allowed alphabets");
      ok = false;
    }
    return ok;
  }

  /**
   * Sequences come from 4 sources.
   * They can be typed, uploaded, embedded or be a pre-existing database.
   * Typed, embedded and uploaded sequences need to be preprocessed to calculate
   * statistics and to ensure they are valid.
   *
   * So what should happen if...
   * 1) An expected field is missing
   *    throw an exception and stop processing the request
   * 2) Sequences have not been sent
   *    whine complaining about the missing sequences
   *    return null
   * 3) Sequences contain a syntax error
   *    whine complaining about the syntax error
   *    return null
   * 4) Sequences violate some constraint
   *    whine complaining about the failed constraint
   *    return null
   * 5) Sequences are valid, pass all constraints
   *    return a sequence source
   *
   *
   * @param restrictedAlphabets specifies what alphabets are allowed after
   *                            other settings have restricted them.
   * @param name manages the name of user supplied sequences and avoids clashes.
   * @param request all the information sent to the webserver.
   * @param feedback an interface for providing error messages to the user.
   * @return a sequence source and return null when a source is not available.
   * @throws ServletException if request details are incorrect (like missing form fields).
   * @throws IOException if storing a parsed version of the sequences to file fails.
   */
  public SequenceInfo getSequences(List<Alph> restrictedAlphabets, FileCoord.Name name,
      HttpServletRequest request, FeedbackHandler feedback) throws ServletException, IOException {
    String source = paramRequire(request, prefix + "_source");
    if (source.equals("noseq")) return null;
    if (source.equals("text") || source.equals("file") || source.equals("embed")) {
      // get the reader from the part and use the file name if possible
      Part part = request.getPart(prefix + "_" + source);
      if (part == null || part.getSize() == 0) {
        feedback.whine("No " + fieldName + " provided.");
        return null; // no sequences submitted
      }
      if (source.equals("file")) {
        name.setOriginalName(getPartFilename(part));
      } else if (source.equals("embed")) {
        name.setOriginalName(request.getParameter(prefix + "_name"));
      }
      SequenceStats statistics = null;
      InputStream in = null;
      File file = null;
      OutputStream out = null;
      boolean success = false;
      try {
        in  = new BufferedInputStream(part.getInputStream());
        file = File.createTempFile("uploaded_sequences_", ".fa");
        file.deleteOnExit();
        out = new BufferedOutputStream(new FileOutputStream(file));
        FeedbackFastaWriter handler = new FeedbackFastaWriter(feedback, fieldName, out);
        FastaParser.parse(in, handler);
        try {out.close();} finally {out = null;}
        try {in.close();} finally {in = null;}
        // check the statistics
        statistics = handler.getStatsRecorder();
        success = checkSpec(restrictedAlphabets, statistics, feedback);
      } catch (FastaException e) {
        // ignore
      } finally {
        closeQuietly(in);
        closeQuietly(out);
        if (file != null && !success) {
          if (!file.delete()) file.deleteOnExit();
        }
      }
      if (success) return new SequenceDataSource(file, name, statistics);
      return null;
    } else { // database
      long listingId, edition;
      // databases only support standard alphabets so we have to convert any restriction alphabets into standard alphabets
      Set<AlphStd> dbAlphs;
      if (restrictedAlphabets != null && restrictedAlphabets.size() > 0) {
        dbAlphs = EnumSet.noneOf(AlphStd.class);
        for (Alph alph : restrictedAlphabets) {
          AlphStd alphStd = AlphStd.fromAlph(alph);
          if (alphStd != null && type.matches(alphStd)) dbAlphs.add(alphStd);
        }
      } else {
        dbAlphs = type.getStandardAlphabets();
      }
      try {
        listingId = Long.parseLong(paramRequire(request, prefix + "_db_listing"), 10);
      } catch (NumberFormatException e) {
        throw new ServletException("Expected the listing ID to be a number", e);
      }
      try {
        edition = Long.parseLong(paramRequire(request, prefix + "_db_version"), 10);
      } catch (NumberFormatException e) {
        throw new ServletException("Expected the version ID to be a number", e);
      }
      // now get a sequence database that matches the listing, version and alphabet
      SequenceDBList db = (SequenceDBList)context.getAttribute(SEQUENCE_DB_KEY);
      if (db == null) {
        throw new ServletException("Unable to access the sequence database.");
      }
      SequenceDB dbFile;
      try {
        dbFile = db.getSequenceFile(listingId, edition, dbAlphs);
      } catch (SQLException e) {
        throw new ServletException(e);
      }
      if (dbFile == null) {
        StringBuilder alphabetText = new StringBuilder();
        AlphStd[] alphabetList = dbAlphs.toArray(new AlphStd[dbAlphs.size()]);
        alphabetText.append(alphabetList[0]);
        for (int i = 1; i < alphabetList.length - 1; i++) {
          alphabetText.append(", ");
          alphabetText.append(alphabetList[i]);
        }
        if (alphabetList.length > 1) {
          alphabetText.append(" or ");
          alphabetText.append(alphabetList[alphabetList.length - 1]);
        }
          feedback.whine("There is no " + alphabetText + " variant of that database for the " + fieldName);
      }
      return dbFile;
    }
  }
  /**
   * Sequences come from 4 sources.
   * They can be typed, uploaded, embedded or be a pre-existing database.
   * Typed, embedded and uploaded sequences need to be preprocessed to calculate
   * statistics and to ensure they are valid.
   *
   * @param name manages the name of user supplied sequences and avoids clashes.
   * @param request all the information sent to the webserver.
   * @param feedback an interface for providing error messages to the user.
   * @return a sequence source and return null when a source is not available.
   * @throws ServletException if request details are incorrect (like missing form fields).
   * @throws IOException if storing a parsed version of the sequences to file fails.
   */
  public SequenceInfo getSequences(FileCoord.Name name, HttpServletRequest request, FeedbackHandler feedback) throws ServletException, IOException {
    return getSequences((List<Alph>)null, name, request, feedback);
  }

  /**
   * Sequences come from 4 sources.
   * They can be typed, uploaded, embedded or be a pre-existing database.
   * Typed, embedded and uploaded sequences need to be preprocessed to calculate
   * statistics and to ensure they are valid.
   *
   * @param restrictedAlphabet specifies which alphabet is allowed.
   * @param name manages the name of user supplied sequences and avoids clashes.
   * @param request all the information sent to the webserver.
   * @param feedback an interface for providing error messages to the user.
   * @return a sequence source and return null when a source is not available.
   * @throws ServletException if request details are incorrect (like missing form fields).
   * @throws IOException if storing a parsed version of the sequences to file fails.
   */
  public SequenceInfo getSequences(Alph restrictedAlphabet, FileCoord.Name name,
      HttpServletRequest request, FeedbackHandler feedback) throws ServletException, IOException {
    return getSequences((restrictedAlphabet != null ? Collections.singletonList(restrictedAlphabet) : null), name, request, feedback);
  }

  /**
   * Sequences come from 4 sources.
   * They can be typed, uploaded, embedded or be a pre-existing database.
   * Typed, embedded and uploaded sequences need to be preprocessed to calculate
   * statistics and to ensure they are valid.
   *
   * @param restrictedAlphabet specifies which alphabet is allowed.
   * @param name manages the name of user supplied sequences and avoids clashes.
   * @param request all the information sent to the webserver.
   * @param feedback an interface for providing error messages to the user.
   * @return a sequence source and return null when a source is not available.
   * @throws ServletException if request details are incorrect (like missing form fields).
   * @throws IOException if storing a parsed version of the sequences to file fails.
   */
  public SequenceInfo getSequences(AlphStd restrictedAlphabet, FileCoord.Name name,
      HttpServletRequest request, FeedbackHandler feedback) throws ServletException, IOException {
    return getSequences((restrictedAlphabet != null ? Collections.singletonList(restrictedAlphabet.getAlph()) : null), name, request, feedback);
  }

  /**
   * Load the selected prior.
   * @param request all the information sent to the webserver.
   * @return a sequence prior if a database was selected and a prior was selected.
   * @throws ServletException if request details are incorrect or there are problems accessing the database.
   */
  public SequencePrior getPrior(HttpServletRequest request) throws ServletException {
    String source = paramRequire(request, prefix + "_source");
    if (!loadPriors || source.equals("text") || source.equals("file") || source.equals("embed")) {
      return null; // only databases have priors currently
    } else {
      String priorId = request.getParameter(prefix + "_db_priors");
      if (priorId == null || priorId.isEmpty()) return null;
      // now get a prior for the id
      SequenceDBList db = (SequenceDBList)context.getAttribute(SEQUENCE_DB_KEY);
      if (db == null) {
        throw new ServletException("Unable to access the sequence database.");
      }
      SequencePrior prior;
      try {
        prior = db.getPrior(Integer.parseInt(priorId));
      } catch (SQLException e) {
        throw new ServletException(e);
      }
      return prior;
    }
  }
}
