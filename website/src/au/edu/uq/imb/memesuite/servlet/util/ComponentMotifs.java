package au.edu.uq.imb.memesuite.servlet.util;

import au.edu.uq.imb.memesuite.data.*;
import au.edu.uq.imb.memesuite.db.Category;
import au.edu.uq.imb.memesuite.db.DBList;
import au.edu.uq.imb.memesuite.db.MotifDB;
import au.edu.uq.imb.memesuite.db.MotifDBList;
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
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.*;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.MOTIF_DB_KEY;

public class ComponentMotifs extends PageComponent {
  private ServletContext context;
  private HTMLTemplate tmplMotifs;
  private AlphType alphType;
  private String prefix;
  private String fieldName;
  private String registerFn;
  private HTMLTemplate title;
  private HTMLTemplate subtitle;
  private boolean enableXAlph;
  private boolean enableDB;
  private boolean enableFilter;
  private DefaultOption defaultOption;

  private static final Pattern DB_ID_PATTERN = Pattern.compile("^\\d+$");
  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.component.motif");

  public static enum DefaultOption {
    TEXT,
    FILE,
    EMBED,
    DATABASE
  }

  private static class CategoryTemplate extends HTMLSubGenerator<Category> {
    private HTMLTemplate template;
    private boolean selectFirst;
    private boolean first;
    boolean selectable;
    boolean identifiable;
    boolean nameable;

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
      first = false;
      return category;
    }
  }

  public static class Selection {
    private static Pattern INDEX_RE = Pattern.compile("^[1-9][0-9]{0,2}$");
    private static Pattern ID_RE = Pattern.compile("^[a-zA-Z0-9:_\\.][a-zA-Z0-9:_\\.\\-]*");
    private String entry;
    private boolean isPos;
    public Selection(String entry) {
      this.entry = entry;
      this.isPos = INDEX_RE.matcher(entry).matches();
      if (!this.isPos && !ID_RE.matcher(entry).matches()) throw new IllegalArgumentException("Parameter entry is not valid");
    }
    public String getEntry() {
      return entry;
    }
    public boolean isPosition() {
      return isPos;
    }
  }

  public ComponentMotifs(ServletContext context, HTMLTemplate info) throws ServletException {
    this.context = context;
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    tmplMotifs = cache.loadAndCache("/WEB-INF/templates/component_motifs.tmpl");
    prefix = getText(info, "prefix", "motifs");
    fieldName = getText(info, "description", "motifs");
    registerFn = getText(info, "register", "nop");
    title = getTemplate(info, "title", null);
    subtitle = getTemplate(info, "subtitle", null);
    alphType = getEnum(info, "alph_type", AlphType.class, AlphType.ANY_ALPHABET);
    enableXAlph = info.containsSubtemplate("enable_xalph");
    enableDB = info.containsSubtemplate("enable_db");
    enableFilter = info.containsSubtemplate("enable_filter");
    defaultOption = getEnum(info, "default", DefaultOption.class, DefaultOption.FILE);
  }

  @Override
  public HTMLSub getComponent() {
    return getComponent(null, null);
  }

  public HTMLSub getComponent(String embedMotifs) {
    return getComponent(embedMotifs, null);
  }

  public HTMLSub getComponent(String embedMotifs, String embedName) {
    boolean selectFirstDB = false;
    DefaultOption defaultOption = this.defaultOption;
    HTMLSub motifs = tmplMotifs.getSubtemplate("component").toSub();
    motifs.set("prefix", prefix);
    if (title != null) motifs.set("title", title);
    if (subtitle != null) motifs.set("subtitle", subtitle);
    if (embedMotifs != null) {
      motifs.getSub("embed_section").set("prefix", prefix).set("data",
          WebUtils.escapeForXML(embedMotifs));
      if (embedName != null) {
        motifs.getSub("embed_section").set("name", WebUtils.escapeForXML(embedName));
      }
      defaultOption = DefaultOption.EMBED;
    } else {
      motifs.empty("embed_option");
      motifs.empty("embed_section");
    }
    if (enableFilter) {
      motifs.getSub("filter_section").set("prefix", prefix);
    } else {
      motifs.empty("filter_section");
    }
    switch (defaultOption) {
      case TEXT:
        motifs.getSub("text_option").set("selected", "selected");
        break;
      case FILE:
        motifs.getSub("file_option").set("selected", "selected");
        break;
      case EMBED:
        if (embedMotifs != null) motifs.getSub("embed_option").set("selected", "selected");
        break;
      case DATABASE:
        selectFirstDB = true;
    }
    MotifDBList db = (MotifDBList)context.getAttribute(MOTIF_DB_KEY);
    if (enableDB && db != null) {
      HTMLTemplate db_opt = tmplMotifs.getSubtemplate("component").
          getSubtemplate("cat_options").getSubtemplate("cat_option");
      try {
        motifs.getSub("cat_options").set("cat_option",
            new CategoryTemplate(db, db_opt, selectFirstDB, false, alphType.getStandardAlphabets()));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Unable to query motif listings", e);
        motifs.empty("cat_options");
      }
      motifs.getSub("db_section").set("prefix", prefix);
    } else {
      motifs.empty("cat_options");
      motifs.empty("db_section");
    }
    if (enableXAlph) {
      motifs.getSub("xalph_section").set("prefix", prefix);
    } else {
      motifs.empty("xalph_section");
    }

    StringWriter buf = new StringWriter();
    JsonWr jsonWr = new JsonWr(buf, 18);
    try {
      jsonWr.start();
      jsonWr.property("field", fieldName);
      jsonWr.property("alph_type", alphType.name());
      jsonWr.end();
    } catch (IOException e) {
      // no IO exceptions should occur as this uses a StringBuffer
      throw new Error(e);
    }
    motifs.set("options", buf.toString());
    motifs.set("register_component", registerFn);
    return motifs;
  }

  public HTMLSub getHelp() {
    return this.tmplMotifs.getSubtemplate("help").toSub();
  }

  public MotifInfo getMotifs(Set<Alph> restrictedAlphabets, FileCoord.Name name, HttpServletRequest request,
      FeedbackHandler feedback) throws ServletException, IOException {
    // determine the source
    String source = paramRequire(request, prefix + "_source");
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
      // check for a filter
      List<Selection> selections = new ArrayList<Selection>();
      if (enableFilter && source.equals("file") && paramBool(request, prefix + "_filter_on")) {
        String[] ids = paramRequire(request, prefix + "_filter").trim().split("\\s+");
        for (String id : ids) {
          try {
            selections.add(new Selection(id));
          } catch (IllegalArgumentException e) {
            feedback.whine("A motif identifier \"" + WebUtils.escapeForXML(id) +
                "\" in the filter field for the " + fieldName +
                " does not fit the allowed pattern.");
          }
        }
      }
      MotifStats statistics = null;
      InputStream in = null;
      File file = null;
      OutputStream out = null;
      boolean success = false;
      try {
        in  = new BufferedInputStream(part.getInputStream());
        file = File.createTempFile("uploaded_motifs_", ".fa");
        file.deleteOnExit();
        out = new BufferedOutputStream(new FileOutputStream(file));
        // copy to a temporary file
        byte[] buffer = new byte[10240]; // 10KB
        int len;
        while ((len = in.read(buffer)) != -1) {
            out.write(buffer, 0, len);
        }
        try {out.close();} finally {out = null;}
        try {in.close();} finally {in = null;}
        statistics = MotifValidator.validate(file);
        if (statistics == null) {
          feedback.whine("The " + fieldName + " did not pass validation.");
        } else if (!alphType.matches(statistics.getAlphabet())) {
          feedback.whine(
              "The alphabet of the " + fieldName + " seems to be " +
                  statistics.getAlphabet() + " but it is not one of the allowed alphabets");
        } else {
          success = true;
        }
      } finally {
        closeQuietly(in);
        closeQuietly(out);
        if (file != null && !success) {
          if (!file.delete()) {
            logger.log(Level.WARNING, "Failed to delete temporary file \"" + file +
                "\". A second attempt will be made at exit.");
          }
        }
      }
      if (success) return new MotifDataSource(file, name, statistics, selections);
      return null;
    } else {
      String listing = paramRequire(request, prefix + "_db_listing");
      Matcher m = DB_ID_PATTERN.matcher(listing);
      if (!m.matches()) {
        throw new ServletException("Parameter " + prefix + "_db_listing had a " +
            "value that did not match any of the allowed values.");
      }
      long dbId;
      try {
        dbId = Long.parseLong(listing, 10);
      } catch (NumberFormatException e) {
        throw new ServletException(e);
      }
      MotifDBList db = (MotifDBList)context.getAttribute(MOTIF_DB_KEY);
      if (db == null) {
        throw new ServletException("Unable to access the motif database.");
      }
      MotifDB motifDB;
      try {
        motifDB =  db.getMotifListing(dbId);
      } catch (SQLException e) {
        throw new ServletException(e);
      }
      if (!alphType.matches(motifDB.getAlphabet())) {
        feedback.whine(
            "The alphabet of the " + fieldName + " is " +
                motifDB.getAlphabet() + " but it is not one of the allowed alphabets");
      }
      return motifDB;
    }
  }

  public MotifInfo getMotifs(AlphStd restrictedAlphabet, FileCoord.Name name, HttpServletRequest request,
      FeedbackHandler feedback) throws ServletException, IOException {
    return getMotifs(Collections.singleton(restrictedAlphabet.getAlph()), name, request, feedback);
  }

  public MotifInfo getMotifs(Alph restrictedAlphabet, FileCoord.Name name, HttpServletRequest request,
      FeedbackHandler feedback) throws ServletException, IOException {
    return getMotifs(Collections.singleton(restrictedAlphabet), name, request, feedback);
  }

  public MotifInfo getMotifs(FileCoord.Name name, HttpServletRequest request,
      FeedbackHandler feedback) throws ServletException, IOException {
    return getMotifs((Set<Alph>)null, name, request, feedback);
  }

  public boolean isAlphabetExpansionAllowed(HttpServletRequest request) {
    return paramBool(request, prefix + "_xalph");
  }

}
