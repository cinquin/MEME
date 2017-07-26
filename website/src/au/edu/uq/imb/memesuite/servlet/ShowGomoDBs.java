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
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.GOMO_DB_KEY;

/**
 * Show the databases available to the GOMO program.
 */
public class ShowGomoDBs extends HttpServlet {
  private ServletContext context;
  private HTMLTemplate template;
  private HTMLTemplate categoryLiTemplate;
  private HTMLTemplate categoryTemplate;
  private HTMLTemplate listingTemplate;
  private HTMLTemplate secondaryTemplate;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web");

  public ShowGomoDBs() { }

  @Override
  public void init() throws ServletException {
    context = this.getServletContext();
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    template = cache.loadAndCache("/WEB-INF/templates/show_gomo_dbs.tmpl");
    categoryLiTemplate = template.getSubtemplate("content").getSubtemplate("category_li");
    categoryTemplate = template.getSubtemplate("content").getSubtemplate("category");
    listingTemplate = categoryTemplate.getSubtemplate("listing");
    secondaryTemplate = listingTemplate.getSubtemplate("secondaries").getSubtemplate("secondary");
  }

  @Override
  public void doGet(HttpServletRequest request, HttpServletResponse response)
      throws IOException, ServletException {
    if (request.getParameter("category") != null) {
      outputXmlListingsOfCategory(response, getId(request, "category"));
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
    GomoDBList gomoDBList = (GomoDBList)context.getAttribute(GOMO_DB_KEY);
    response.setContentType("text/html; charset=UTF-8");
    HTMLSub out = template.toSub();
    if (gomoDBList != null) {
      try {
        out.getSub("content").set("category_li", new AllCategoryList(gomoDBList));
        out.getSub("content").set("category", new AllCategories(gomoDBList));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error loading categories", e);
        out.empty("content");
      }
    } else {
      out.empty("content");
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

  private void outputXmlListingsOfCategory(HttpServletResponse response, long categoryId)
      throws IOException, ServletException {
    GomoDBList gomoDBList = (GomoDBList)context.getAttribute(GOMO_DB_KEY);
    response.setContentType("application/xml; charset=UTF-8");
    // turn off caching
//    response.setHeader("Cache-Control", "no-cache, no-store, must-revalidate"); // HTTP 1.1.
//    response.setHeader("Pragma", "no-cache"); // HTTP 1.0.
//    response.setDateHeader("Expires", 0); // Proxies.
    PrintWriter out = null;
    try {
      Iterator<Listing> listingView = gomoDBList.getListings(categoryId).iterator();
      out = response.getWriter();
      out.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
      out.println("<listings category=\"" + categoryId + "\">");
      while (listingView.hasNext()) {
        Listing listing = listingView.next();
        out.println("<l i=\"" + listing.getID() + "\" n=\"" + WebUtils.escapeForXML(listing.getName()) + "\"/>");
      }
      out.println("</listings>");
    } catch (SQLException e) {
      throw new ServletException("Error getting GOMO listings", e);
    } finally {
      if (out != null) out.close();
    }
  }

  private class AllCategoryList extends HTMLSubGenerator<Category> {
    private AllCategoryList(DBList gomoDBList) throws SQLException {
      super(gomoDBList.getCategories(false, EnumSet.allOf(AlphStd.class)));
    }

    @Override
    protected HTMLSub transform(Category item) {
      HTMLSub out = categoryLiTemplate.toSub();
      out.set("id", item.getID());
      out.set("name", item.getName());
      return out;
    }
  }

  private class AllCategories extends HTMLSubGenerator<Category> {
    private GomoDBList gomoDBList;
    private AllCategories(GomoDBList gomoDBList) throws SQLException {
      super(gomoDBList.getCategories(false, EnumSet.allOf(AlphStd.class)));
      this.gomoDBList = gomoDBList;
    }

    @Override
    protected HTMLSub transform(Category item) {
      HTMLSub out = categoryTemplate.toSub();
      out.set("id", item.getID());
      out.set("name", item.getName());
      try {
        out.set("listing", new AllListingsOfCategory(gomoDBList, item.getID()));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error loading GOMO listings", e);
        out.empty("listing");
      }
      return out;
    }
  }

  private class AllListingsOfCategory extends HTMLSubGenerator<Listing> {
    private GomoDBList gomoDBList;
    private AllListingsOfCategory(GomoDBList gomoDBList, long id) throws SQLException{
      super(gomoDBList.getListings(id));
      this.gomoDBList = gomoDBList;
    }

    @Override
    protected HTMLSub transform(Listing listing) {
      HTMLSub out = listingTemplate.toSub();
      out.set("name", listing.getName());
      out.set("description", listing.getDescription());
      try {
        GomoDB gomoDB = gomoDBList.getGomoListing(listing.getID());
        List<GomoDBSecondary> secondaries = gomoDB.getSecondaries();
        if (!secondaries.isEmpty()) {
          List<HTMLSub> list = new ArrayList<HTMLSub>();
          for (GomoDBSecondary secondary : secondaries) {
            list.add(secondaryTemplate.toSub().set("description", secondary.getDescription()));
          }
          out.getSub("secondaries").set("secondary", list);
        } else {
          out.empty("secondaries");
        }
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error querying GOMO listing",  e);
        out.empty("secondaries");
      }
      return out;
    }
  }

}
