package au.edu.uq.imb.memesuite.servlet.util;

import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;

import javax.servlet.ServletException;

/**
 * Component containing the submit and reset buttons.
 */
public class ComponentSubmitReset extends PageComponent {
  private HTMLTemplate tmplSubmitReset;
  private String submitTitle;
  private String resetTitle;
  private int quotaCount;
  private long quotaDuration;

  public ComponentSubmitReset(HTMLTemplateCache cache, int quotaCount, long quotaDuration) throws ServletException {
    tmplSubmitReset = cache.loadAndCache("/WEB-INF/templates/component_submit_reset.tmpl");
    submitTitle = "Start Search";
    resetTitle = "Clear Input";
    this.quotaCount = quotaCount;
    this.quotaDuration = quotaDuration;
  }

  public HTMLSub getComponent() {
    return getComponent(0);
  }

  public HTMLSub getComponent(long quotaMinWait) {
    HTMLSub sub = tmplSubmitReset.getSubtemplate("component").toSub();
    sub.set("submit_title", submitTitle);
    sub.set("reset_title", resetTitle);
    sub.set("submit_wait", quotaMinWait);
    sub.set("quota_count", quotaCount);
    sub.set("quota_duration", quotaDuration);
    return sub;
  }

  public HTMLSub getHelp() {
    return tmplSubmitReset.getSubtemplate("help").toSub();
  }
}
