package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.AlphStd;
import au.edu.uq.imb.memesuite.db.*;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLSubGenerator;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;
import au.edu.uq.imb.memesuite.util.JsonWr;

import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.SEQUENCE_DB_KEY;

/**
 * Display the available sequence databases.
 */
public class ShowSequenceDBs extends HttpServlet {
  private ServletContext context;
  private HTMLTemplate template;
  private HTMLTemplate categoryTemplate;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.sequencedb");

  public ShowSequenceDBs() { }

  @Override
  public void init() throws ServletException {
    this.context = this.getServletContext();
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    template = cache.loadAndCache("/WEB-INF/templates/show_sequence_dbs.tmpl");
    categoryTemplate = template.getSubtemplate("category");
  }

  @Override
  public void doGet(HttpServletRequest request, HttpServletResponse response)
      throws IOException, ServletException {
    if (request.getParameter("category") != null) {
      outputJsonListingsOfCategory(response, getId(request, "category"),
          isShortOnly(request), allowedAlphabets(request));
    } else if (request.getParameter("listing") != null) {
      outputJsonVersionsOfListing(response, getId(request, "listing"),
          isShortOnly(request), allowedAlphabets(request));
    } else if (request.getParameter("sequence") != null) {
      outputJsonPriorsOfSequence(response, getId(request, "sequence"));
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
    SequenceDBList sequenceDBList = (SequenceDBList)context.getAttribute(SEQUENCE_DB_KEY);
    response.setContentType("text/html; charset=UTF-8");
    HTMLSub out = template.toSub();
    if (sequenceDBList != null) {
      try {
        out.set("category", new AllCategory(sequenceDBList));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error loading sequence categories", e);
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

  private <E extends Enum> void outputEnumSetAsArray(JsonWr out, String property, Set<E> set) throws IOException{
    out.property(property);
    out.startArray();
    for (E value : set) {
      out.value(value.name());
    }
    out.endArray();
  }

  private void outputJsonListingsOfCategory(HttpServletResponse response,
      long categoryId, boolean shortOnly, EnumSet<AlphStd> allowedAlphabets)
      throws IOException, ServletException {
    SequenceDBList sequenceDBList = (SequenceDBList)context.getAttribute(SEQUENCE_DB_KEY);
    response.setContentType("application/json; charset=UTF-8");
    try {
      Iterator<Listing> listingView = sequenceDBList.getListings(categoryId, shortOnly, allowedAlphabets).iterator();
      JsonWr out = new JsonWr(response.getWriter());
      out.start();
      out.property("listings");
      out.startArray();
      while (listingView.hasNext()) {
        Listing listing = listingView.next();
        out.startObject();
        out.property("id", listing.getID());
        out.property("name", listing.getName());
        outputEnumSetAsArray(out, "alphabets", listing.getAlphabets());
        out.property("hasPriors", listing.hasPriors());
        out.endObject();
      }
      out.endArray();
      out.end();
    } catch (SQLException e) {
      throw new ServletException(e);
    }
  }

  private void outputJsonVersionsOfListing(HttpServletResponse response,
      long listingId, boolean shortOnly, EnumSet<AlphStd> allowedAlphabets)
      throws IOException, ServletException {
    SequenceDBList sequenceDBList = (SequenceDBList)context.getAttribute(SEQUENCE_DB_KEY);
    response.setContentType("application/json; charset=UTF-8");
    try {
      List<SequenceVersion> versions = sequenceDBList.getVersions(listingId, shortOnly, allowedAlphabets);
      JsonWr out = new JsonWr(response.getWriter());
      String listingName = "";
      String listingDescription = "";
      if (!versions.isEmpty()) {
        SequenceVersion ver = versions.get(0);
        listingName = ver.getListingName();
        listingDescription = ver.getListingDescription();
      }
      // calculate the combined alphabets of the versions
      EnumSet<AlphStd> availableAlphabets = EnumSet.noneOf(AlphStd.class);
      for (SequenceVersion version : versions) {
        availableAlphabets.addAll(version.getAlphabets());
      }

      out.start();
      out.property("name", listingName);
      out.property("description", listingDescription);
      outputEnumSetAsArray(out, "alphabets", availableAlphabets);
      out.property("versions");
      out.startArray();
      for (SequenceVersion version : versions) {
        out.startObject();
        out.property("id", version.getEdition());
        out.property("name", version.getVersion());
        outputEnumSetAsArray(out, "alphabets", version.getAlphabets());
        out.property("sequences");
        out.startObject();
        for (AlphStd alphabet : version.getAlphabets()) {
          SequenceDB file = version.getSequenceFile(alphabet);
          out.property(alphabet.name());
          out.startObject();
          out.property("id", file.getId());
          out.property("count", file.getSequenceCount());
          out.property("min", file.getMinLength());
          out.property("max", file.getMaxLength());
          out.property("avg", file.getAverageLength());
          out.property("total", file.getTotalLength());
          out.property("description", file.getDescription());
          out.property("priorCount", file.getPriorCount());
          out.endObject();
        }
        out.endObject();
        out.endObject();
      }
      out.endArray();
      out.end();
    } catch (SQLException e) {
      throw new ServletException(e);
    }
  }

  private void outputJsonPriorsOfSequence(HttpServletResponse response, long sequenceId) throws IOException, ServletException {
    SequenceDBList sequenceDBList = (SequenceDBList)context.getAttribute(SEQUENCE_DB_KEY);
    response.setContentType("application/json; charset=UTF-8");
    try {
      List<SequencePrior> priors = sequenceDBList.getPriors(sequenceId);
      JsonWr out = new JsonWr(response.getWriter());
      out.start();
      out.property("sequence", sequenceId);
      out.property("priors");
      out.startArray();
      for (SequencePrior prior : priors) {
        out.startObject();
        out.property("id", prior.getId());
        out.property("biosample", prior.getBiosample());
        out.property("assay", prior.getAssay());
        out.property("source", prior.getSource());
        out.property("url", prior.getUrl().toString());
        out.property("description", prior.getDescription());
        out.endObject();
      }
      out.endArray();
      out.end();
    } catch (SQLException e) {
      throw new ServletException(e);
    }

  }

  private class AllCategory extends HTMLSubGenerator<Category> {
    private AllCategory(SequenceDBList sequenceDBList) throws SQLException{
      super(sequenceDBList.getCategories(false, EnumSet.allOf(AlphStd.class)));
    }

    @Override
    protected HTMLSub transform(Category item) {
      HTMLSub out = categoryTemplate.toSub();
      out.set("id", item.getID());
      out.set("name", item.getName());
      out.set("count", item.getCount());
      out.set("no_priors", item.hasPriors() ? "" : "no_priors");
      return out;
    }
  }

}
