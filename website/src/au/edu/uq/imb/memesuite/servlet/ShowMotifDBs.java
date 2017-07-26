package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.AlphStd;
import au.edu.uq.imb.memesuite.db.*;
import au.edu.uq.imb.memesuite.servlet.util.WebUtils;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLSubGenerator;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;

import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.EnumSet;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.MOTIF_DB_KEY;


/**
 * Display the available motif databases
 */
public class ShowMotifDBs extends HttpServlet {
  private ServletContext context;
  private HTMLTemplate template;
  private HTMLTemplate categoryTemplate;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.motifdb");

  public ShowMotifDBs() { }

  @Override
  public void init() throws ServletException {
    context = this.getServletContext();
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    template = cache.loadAndCache("/WEB-INF/templates/show_motif_dbs.tmpl");
    categoryTemplate = template.getSubtemplate("category");
  }

  @Override
  public void doGet(HttpServletRequest request, HttpServletResponse response)
      throws IOException, ServletException {
    if (request.getParameter("category") != null) {
      outputXmlListingsOfCategory(response, getId(request, "category"),
          isShortOnly(request), allowedAlphabets(request));
    } else if (request.getParameter("listing") != null) {
      outputXmlListing(response, getId(request, "listing"));
    } else {
      display(response);
    }
  }

  @Override
  public void doPost(HttpServletRequest request, HttpServletResponse response)
      throws IOException, ServletException {
    display(response);
  }

  private void display(HttpServletResponse response)
      throws IOException, ServletException {
    MotifDBList motifDBList = (MotifDBList)context.getAttribute(MOTIF_DB_KEY);
    response.setContentType("text/html; charset=UTF-8");
    HTMLSub out = template.toSub();
    if (motifDBList != null) {
      try {
        out.set("category", new AllCategory(motifDBList));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error reading motif categories", e);
        out.empty("category");
      }
    } else {
      out.empty("category");
    }
    out.output(response.getWriter());
  }

  private long getId(HttpServletRequest request, String name) throws ServletException {
    String value = request.getParameter(name);
    if (value == null) {
      throw new ServletException("Parameter '" + name + "' was not set.");
    }
    long id;
    try {
      id = Long.parseLong(value, 10);
    } catch (NumberFormatException e) {
      throw new ServletException("Parameter '" + name + "' is not a integer value.", e);
    }
    return id;
  }

  private boolean isShortOnly(HttpServletRequest request) throws ServletException{
    String value = request.getParameter("short");
    if (value == null) return false;
    boolean shortOnly;
    try {
      shortOnly = Integer.parseInt(value, 2) != 0;
    } catch (NumberFormatException e) {
      throw new ServletException("Parameter 'short' is not a binary value.", e);
    }
    return shortOnly;
  }

  private EnumSet<AlphStd> allowedAlphabets(HttpServletRequest request) throws ServletException {
    String value = request.getParameter("alphabets");
    if (value == null) return EnumSet.allOf(AlphStd.class);
    EnumSet<AlphStd> allowedAlphabets;
    try {
      allowedAlphabets = SQL.intToEnums(AlphStd.class, Integer.parseInt(value, 10));
    } catch (NumberFormatException e) {
      throw new ServletException("Parameter 'alphabets' is not a bitset.", e);
    } catch (IllegalArgumentException e) {
      throw new ServletException("Parameter 'alphabets' is not a bitset", e);
    }
    return allowedAlphabets;
  }

  private void outputXmlListingsOfCategory(HttpServletResponse response,
      long categoryId, boolean shortOnly, EnumSet<AlphStd> allowedAlphabets)
      throws IOException, ServletException {
    MotifDBList motifDBList = (MotifDBList)context.getAttribute(MOTIF_DB_KEY);
    response.setContentType("application/xml; charset=UTF-8");
    // turn off caching
//    response.setHeader("Cache-Control", "no-cache, no-store, must-revalidate"); // HTTP 1.1.
//    response.setHeader("Pragma", "no-cache"); // HTTP 1.0.
//    response.setDateHeader("Expires", 0); // Proxies.
    PrintWriter out = null;
    try {
      List<Listing> listingView = motifDBList.getListings(categoryId, shortOnly, allowedAlphabets);
      out = response.getWriter();
      out.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
      out.println("<listings category=\"" + categoryId + "\" short=\"" +
          (shortOnly ? "1" : "0") +"\" alphabets=\"" +
          SQL.enumsToInt(allowedAlphabets) + "\">");
      for (Listing listing : listingView) {
        out.println("<l i=\"" + listing.getID() + "\" n=\"" + WebUtils.escapeForXML(listing.getName()) +
            "\" a=\"" + SQL.enumsToInt(listing.getAlphabets()) + "\"/>");
      }
      out.println("</listings>");
    } catch (SQLException e) {
      throw new ServletException(e);
    } finally {
      if (out != null) out.close();
    }
  }


  private void outputXmlListing(HttpServletResponse response, long listingId)
      throws IOException {
    MotifDBList motifDBList = (MotifDBList)context.getAttribute(MOTIF_DB_KEY);
    response.setContentType("application/xml; charset=UTF-8");
    PrintWriter out = null;
    try {
      MotifDB motifDB = motifDBList.getMotifListing(listingId);
      out = response.getWriter();
      out.print("<motif_db id=\"" + listingId + "\" alphabet=\"" +
          motifDB.getAlphabetEn().name() + "\">");
      // name contains HTML, hence can't be escaped or it becomes useless
      out.print("<name><![CDATA[");
      out.print(motifDB.getName());
      out.println("]]></name>");
      // description contains HTML, hence can't be escaped or it becomes useless
      out.print("<description><![CDATA[");
      out.print(motifDB.getDescription());
      out.println("]]></description>");
      for (MotifDBFile file : motifDB.getMotifFiles()) {
        out.println("<file src=\"" +WebUtils.escapeForXML(file.getFileName()) +
            "\" count=\"" + file.getMotifCount() + "\" cols=\"" +
            file.getTotalCols() + "\" min=\"" + file.getMinCols() +
            "\" max=\"" + file.getMaxCols() + "\" />");
      }
      out.println("</motif_db>");
    } catch (SQLException e) {
      throw new IOException(e);
    } finally {
      if (out != null) out.close();
    }

  }

  private class AllCategory extends HTMLSubGenerator<Category> {
    private AllCategory(MotifDBList motifDBList) throws SQLException{
      super(motifDBList.getCategories(false, EnumSet.allOf(AlphStd.class)));
    }

    @Override
    protected HTMLSub transform(Category item) {
      HTMLSub out = categoryTemplate.toSub();
      out.set("id", item.getID());
      out.set("name", item.getName());
      out.set("cnt", item.getCount());
      return out;
    }
  }

}
