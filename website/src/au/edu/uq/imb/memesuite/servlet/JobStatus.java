package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.MemeSuiteProperties;
import au.edu.uq.imb.memesuite.servlet.util.WebUtils;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;

import edu.sdsc.nbcr.opal.AppServicePortType;
import edu.sdsc.nbcr.opal.StatusOutputType;
import org.globus.gram.GramJob;
import org.globus.gram.internal.GRAMConstants;

import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.rmi.RemoteException;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CONFIG_KEY;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;

/**
 * Displays the status of a job
 */
public class JobStatus extends HttpServlet {
  private HTMLTemplate tmplStatus;
  private MemeSuiteProperties msp;
  private Map<String, HTMLTemplate> reports;
  private RateLimit statusQueryLimiter;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.status");

  private static class RateEntry {
    public final String id;
    public final long timestamp;
    public final Long interval;
    protected int count;
    public RateEntry(String id, Long interval, long timestamp) {
      this.id = id;
      this.interval = interval;
      this.timestamp = timestamp;
      count = 1;
    }
    public int getCount() {
      return count;
    }
    public void incrementCount() {
      count++;
    }
  }
  private static class RateLoad {
    public final boolean accepted;
    public final double load;
    public final Long interval;
    public RateLoad(boolean accepted, double load, Long interval) {
      this.accepted = accepted;
      this.load = load;
      this.interval = interval;
    }
  }
  private static class RateLimit {
    public final int pollDelay;
    public final int maxPerSec;
    protected RateEntry[] buffer;
    protected int current;

    /**
     *
     * @param pollDelay the minimum delay allowed between checking the same ID
     * @param maxPerSec the maximum checks allowed per second
     */
    public RateLimit(int pollDelay, int maxPerSec) {
      this.pollDelay = pollDelay;
      this.maxPerSec = maxPerSec;
      buffer = new RateEntry[pollDelay * maxPerSec];
      current = 0;
      // fill the buffer with expired entries
      Arrays.fill(this.buffer, new RateEntry("", 0l, 0));
    }

    /**
     *
     * @param id the ID to check if allowed
     * @param interval the timer interval used by this source for polling
     * @return a rate load object
     */
    public synchronized RateLoad check(String id, Long interval) {
      long now = System.currentTimeMillis();
      long recent = now - 1000;
      long expired = now - (pollDelay * 1000);
      int recentAcceptedRequests = 0;
      int allRequests = 0;
      boolean seen = false;
      long intervalSum = 0;
      long intervalCount = 0;
      outer:
      {
        for (int i = current - 1; i >= 0; allRequests += buffer[i].getCount(), i--) {
          if (buffer[i].timestamp < expired) break outer;
          if (buffer[i].timestamp >= recent) recentAcceptedRequests++;
          if (buffer[i].id.hashCode() == id.hashCode() && buffer[i].id.equals(id)) seen = true;
          if (buffer[i].interval != null) {
            intervalSum += buffer[i].interval;
            intervalCount++;
          }
        }
        for (int i = buffer.length - 1; i >= current; allRequests += buffer[i].getCount(), i--) {
          if (buffer[i].timestamp < expired) break outer;
          if (buffer[i].timestamp >= recent) recentAcceptedRequests++;
          if (buffer[i].id.hashCode() == id.hashCode() && buffer[i].id.equals(id)) seen = true;
          if (buffer[i].interval != null) {
            intervalSum += buffer[i].interval;
            intervalCount++;
          }
        }
      }
      double load = (allRequests + 1.0) / buffer.length;
      Long avgInterval = (intervalCount > 0 ? intervalSum / intervalCount : null);
      if (seen) return new RateLoad(false, load, avgInterval); // reject, seen recently
      if (buffer[current].timestamp >= expired || recentAcceptedRequests > maxPerSec) { // reject, flooded
        buffer[(current + (buffer.length - 1)) % buffer.length].incrementCount();
        return new RateLoad(false, load, avgInterval);
      }
      buffer[current] = new RateEntry(id, interval, now);
      current = (current + 1) % buffer.length;
      return new RateLoad(true, load, avgInterval);
    }
  }

  @Override
  public void init() throws ServletException {
    ServletContext context = this.getServletContext();
    msp = (MemeSuiteProperties)context.getAttribute(CONFIG_KEY);
    if (msp == null) throw new ServletException("Failed to get MEME Suite properties");
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    tmplStatus = cache.loadAndCache("/WEB-INF/templates/status.tmpl");
    reports = new TreeMap<String, HTMLTemplate>();
    reports.put("AME", cache.loadAndCache("/WEB-INF/templates/ame_verify.tmpl"));
    reports.put("CENTRIMO", cache.loadAndCache("/WEB-INF/templates/centrimo_verify.tmpl"));
    reports.put("DREME", cache.loadAndCache("/WEB-INF/templates/dreme_verify.tmpl"));
    reports.put("FIMO", cache.loadAndCache("/WEB-INF/templates/fimo_verify.tmpl"));
    reports.put("GLAM2", cache.loadAndCache("/WEB-INF/templates/glam2_verify.tmpl"));
    reports.put("GLAM2SCAN", cache.loadAndCache("/WEB-INF/templates/glam2scan_verify.tmpl"));
    reports.put("GOMO", cache.loadAndCache("/WEB-INF/templates/gomo_verify.tmpl"));
    reports.put("MAST", cache.loadAndCache("/WEB-INF/templates/mast_verify.tmpl"));
    reports.put("MCAST", cache.loadAndCache("/WEB-INF/templates/mcast_verify.tmpl"));
    reports.put("MEME", cache.loadAndCache("/WEB-INF/templates/meme_verify.tmpl"));
    reports.put("MEMECHIP", cache.loadAndCache("/WEB-INF/templates/meme-chip_verify.tmpl"));
    reports.put("MOMO", cache.loadAndCache("/WEB-INF/templates/momo_verify.tmpl"));
    reports.put("SPAMO", cache.loadAndCache("/WEB-INF/templates/spamo_verify.tmpl"));
    reports.put("TOMTOM", cache.loadAndCache("/WEB-INF/templates/tomtom_verify.tmpl"));
    statusQueryLimiter = new RateLimit(10, 2);
  }

  protected static boolean urlExists(String url) {
    logger.log(Level.INFO, "Checking URL \"" + url + "\"");
    url = url.replaceFirst("https", "http");
    try {
      URLConnection connection = new URL(url).openConnection();
      if (connection instanceof HttpURLConnection) {
        HttpURLConnection http = (HttpURLConnection)connection;
        http.setRequestMethod("HEAD");
        int responseCode = http.getResponseCode();
        logger.log(Level.INFO, "Got response code " + responseCode + " for URL \"" + url + "\"");
        return responseCode == 200;
      } else {
        logger.log(Level.WARNING, "Unable to cast URLConnection to HttpURLConnection!");
        return false;
      }
    } catch (IOException e) {
      logger.log(Level.WARNING, "Failed to check URL \"" + url + "\"", e);
      return false;
    }
  }

  private static String jobStatus(StatusOutputType status) {
    switch (status.getCode()) {
      case GRAMConstants.STATUS_UNSUBMITTED:
      case GRAMConstants.STATUS_PENDING:
      case GRAMConstants.STATUS_STAGE_IN:
        return "pending";
      case GRAMConstants.STATUS_ACTIVE:
      case GRAMConstants.STATUS_STAGE_OUT:
        return "active";
      case GRAMConstants.STATUS_SUSPENDED:
        return "suspended";
      case GRAMConstants.STATUS_FAILED:
        return "failed";
      case GRAMConstants.STATUS_DONE:
        return "done";
      default:
        return "unknown";
    }
  }

  private static String jobUrl(StatusOutputType status) {
    String baseURL = status.getBaseURL().toString();
    return baseURL + (baseURL.endsWith("/") ? "" : "/") + "index.html";
  }

  @Override
  protected void doGet(HttpServletRequest req, HttpServletResponse resp) throws ServletException, IOException {
    String serviceName = req.getParameter("service");
    String jobId = req.getParameter("id");
    if (req.getParameter("xml") == null) {
      generateHtml(serviceName, jobId, resp);
    } else {
      String intervalStr = req.getParameter("interval");
      Long interval = null;
      try { interval = Long.parseLong(intervalStr, 10); } catch (NumberFormatException e) { /* ignore */ }
      generateXml(serviceName, jobId, interval, resp);
    }
  }

  protected void generateXml(String serviceName, String jobId, Long interval, HttpServletResponse resp) throws ServletException, IOException {
    RateLoad limiter = statusQueryLimiter.check(jobId, interval);
    StatusOutputType status = null;
    if (limiter.accepted && reports.containsKey(serviceName) && jobId != null) {
      try {
        AppServicePortType app = WebUtils.getOpal(msp, serviceName);
        status = app.queryStatus(jobId);
      } catch (RemoteException e) { /* ignore */  }
    }
    resp.setContentType("application/xml; charset=UTF-8");
    PrintWriter out = resp.getWriter();
    if (status != null) {
      String url = jobUrl(status);
      out.println("<job>");
      out.printf("<load>%.2f</load>\n", limiter.load);
      if (limiter.interval != null) {
        out.printf("<interval>%d</interval>\n", limiter.interval);
      }
      out.print("<status>");
      out.print(jobStatus(status));
      out.println("</status>");
      if (urlExists(url)) {
        out.print("<url>");
        out.print(url);
        out.println("</url>");
      }
      out.println("</job>");
    } else if (limiter.accepted) {
      out.println("<job>");
      out.printf("<load>%.2f</load>\n", limiter.load);
      if (limiter.interval != null) {
        out.printf("<interval>%d</interval>\n", limiter.interval);
      }
      out.println("<status>unknown</status>");
      out.println("</job>");
    } else {
      out.println("<job>");
      out.printf("<load>%.2f</load>\n", limiter.load);
      if (limiter.interval != null) {
        out.printf("<interval>%d</interval>\n", limiter.interval);
      }
      out.println("</job>");
    }
  }

  protected void generateHtml(String serviceName, String jobId, HttpServletResponse resp) throws ServletException, IOException {
    HTMLSub page = this.tmplStatus.toSub();
    // look up parameters
    // try to get info for service
    HTMLTemplate serviceInfo = reports.get(serviceName);
    if (serviceInfo != null) {
      // Yay, I've heard about this service!
      // input report information
      page.set("report", serviceInfo.getSubtemplate("message"));
      // setup page header
      page.empty("suite_header");
      page.getSub("program_header").
          set("title", serviceInfo.getSubtemplate("title")).
          set("subtitle", serviceInfo.getSubtemplate("subtitle")).
          set("logo", "../" + serviceInfo.getSubtemplate("logo")).
          set("alt", serviceInfo.getSubtemplate("alt"));
      page.set("service", serviceName);

      // see if we can discover information about the job
      StatusOutputType status = null;
      if (jobId != null) {
        try {
          AppServicePortType app = WebUtils.getOpal(msp, serviceName);
          status = app.queryStatus(jobId);
        } catch (RemoteException e) { /* ignore */  }
      }
      if (status != null) {
        // Yay, I know about this job!
        page.set("id", jobId);
        // set the displayed status message
        page.set("status", jobStatus(status));
        page.set("title", serviceInfo.getSubtemplate("title"));
        // display the output index if it has been created
        String fullURL = jobUrl(status);
        if (urlExists(fullURL)) {
          // note we randomize the ID so safari doesn't cache the content
          page.getSub("preview").set("url", fullURL).set("id", "IF_" + System.currentTimeMillis());
        } else {
          page.set("details", "expanded");
          page.empty("preview");
        }
      } else {
        // don't know about this job!
        page.set("status", "expired");
        page.empty("preview");
      }
    } else {
      // don't know about this service!
      page.set("report", "{}");
      page.empty("program_header");
      page.set("status", "expired");
    }
    resp.setContentType("text/html; charset=UTF-8");
    page.output(resp.getWriter());
  }
}
