package au.edu.uq.imb.memesuite.servlet.util;

import au.edu.uq.imb.memesuite.data.Alph;
import au.edu.uq.imb.memesuite.data.AlphType;
import au.edu.uq.imb.memesuite.data.AlphabetDataSource;
import au.edu.uq.imb.memesuite.io.alph.AlphParser;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;
import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.Part;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;
import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.*;
import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.closeQuietly;

/**
 * A component used for inputting the alphabet of Fasta sequences.
 */
public class ComponentSequenceAlphabet extends PageComponent {
  private HTMLTemplate template;
  private String prefix;
  private String fieldName;
  private String registerFn;
  private HTMLTemplate title;
  private HTMLTemplate subtitle;
  private AlphType type;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.component.alphabet");

  public ComponentSequenceAlphabet(ServletContext context, HTMLTemplate info) throws ServletException {
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    template = cache.loadAndCache("/WEB-INF/templates/component_sequences.tmpl");
    prefix = getText(info, "prefix", "sequences");
    fieldName = getText(info, "description", "sequences");
    registerFn = getText(info, "register", "nop");
    title = getTemplate(info, "title", null);
    subtitle = getTemplate(info, "subtitle", null);
    type = getEnum(info, "alph_type", AlphType.class, AlphType.ANY_ALPHABET);
  }


  @Override
  public HTMLSub getComponent() {
    HTMLSub out = template.getSubtemplate("alphabet_selector").toSub();
    out.set("prefix2", prefix);
    out.set("title", title);
    out.set("subtitle", subtitle);
    out.set("register_component", registerFn);
    StringWriter buf = new StringWriter();
    JsonWr jsonWr = new JsonWr(buf, 18);
    try {
      jsonWr.start();
      jsonWr.property("field", fieldName);
      jsonWr.property("alph_type", type.name());
      jsonWr.end();
    } catch (IOException e) {
      // no IO exceptions should occur as this uses a StringBuffer
      throw new Error(e);
    }
    out.set("options2", buf.toString());
    return out;
  }

  @Override
  public HTMLSub getHelp() {
    return template.getSubtemplate("alphabet_selector_help").toSub();
  }

  public AlphabetDataSource getAlphabet(FileCoord.Name name, HttpServletRequest request, FeedbackHandler feedback) throws ServletException, IOException {
    if (!paramBool(request, prefix + "_custom")) {
      // default DNA, RNA or PROTEIN
      return null;
    }
    // Some custom alphabet
    // get the reader from the part and use the file name if possible
    Part part = request.getPart(prefix + "_file");
    if (part == null || part.getSize() == 0) {
      feedback.whine("No " + fieldName + " provided.");
      return null; // no sequences submitted
    }
    name.setOriginalName(getPartFilename(part));
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
      try {
        AlphabetDataSource dataSource = new AlphabetDataSource(file, name, AlphParser.parseFile(file.toPath()));
        success = true;
        return dataSource;
      } catch (AlphParser.AlphParseException e) {
        feedback.whine("The " + fieldName + " did not pass validation.");
      }
      return null;
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
  }
}
