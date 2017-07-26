
/*
 * $
 *
 * Shorthand function for getElementById
 */
function $(el) {
  return document.getElementById(el);
}


/*
 * See http://stackoverflow.com/a/5450113/66387
 * Does string multiplication like the perl x operator.
 */
function string_mult(pattern, count) {
    if (count < 1) return '';
    var result = '';
    while (count > 1) {
        if (count & 1) result += pattern;
        count >>= 1, pattern += pattern;
    }
    return result + pattern;
}

/*
 * See http://stackoverflow.com/questions/814613/how-to-read-get-data-from-a-url-using-javascript
 * Slightly modified with information from
 * https://developer.mozilla.org/en/DOM/window.location
 */
function parse_params() {
  "use strict";
  var search, queryStart, queryEnd, query, params, nvPairs, i, nv, n, v;
  search = window.location.search;
  queryStart = search.indexOf("?") + 1;
  queryEnd   = search.indexOf("#") + 1 || search.length + 1;
  query      = search.slice(queryStart, queryEnd - 1);

  if (query === search || query === "") return {};

  params  = {};
  nvPairs = query.replace(/\+/g, " ").split("&");

  for (i = 0; i < nvPairs.length; i++) {
    nv = nvPairs[i].split("=");
    n  = decodeURIComponent(nv[0]);
    v  = decodeURIComponent(nv[1]);
    // allow a name to be used multiple times
    // storing each value in the array
    if (!(n in params)) {
      params[n] = [];
    }
    params[n].push(nv.length === 2 ? v : null);
  }
  return params;
}

/*
 * coords
 *
 * Calculates the x and y offset of an element.
 * From http://www.quirksmode.org/js/findpos.html
 * with alterations to take into account scrolling regions
 */
function coords(elem) {
  var myX = myY = 0;
  if (elem.getBoundingClientRect) {
    var rect;
    rect = elem.getBoundingClientRect();
    myX = rect.left + ((typeof window.pageXOffset !== "undefined") ?
        window.pageXOffset : document.body.scrollLeft);
    myY = rect.top + ((typeof window.pageYOffset !== "undefined") ?
        window.pageYOffset : document.body.scrollTop);
  } else {
    // this fall back doesn't properly handle absolutely positioned elements
    // inside a scrollable box
    var node;
    if (elem.offsetParent) {
      // subtract all scrolling
      node = elem;
      do {
        myX -= node.scrollLeft ? node.scrollLeft : 0;
        myY -= node.scrollTop ? node.scrollTop : 0;
      } while (node = node.parentNode);
      // this will include the page scrolling (which is unwanted) so add it back on
      myX += (typeof window.pageXOffset !== "undefined") ? window.pageXOffset : document.body.scrollLeft;
      myY += (typeof window.pageYOffset !== "undefined") ? window.pageYOffset : document.body.scrollTop;
      // sum up offsets
      node = elem;
      do {
        myX += node.offsetLeft;
        myY += node.offsetTop;
      } while (node = node.offsetParent);
    }
  }
  return [myX, myY];
}

/*
 * position_popup
 *
 * Positions a popup relative to an anchor element.
 *
 * The avaliable positions are:
 * 0 - Centered below the anchor.
 */
function position_popup(anchor, popup, position) {
  "use strict";
  var a_x, a_y, a_w, a_h, p_x, p_y, p_w, p_h;
  var a_xy, spacer, margin, scrollbar, page_w;
  // define constants
  spacer = 5;
  margin = 15;
  scrollbar = 15;
  // define the positions and widths
  a_xy = coords(anchor);
  a_x = a_xy[0];
  a_y = a_xy[1];
  a_w = anchor.offsetWidth;
  a_h = anchor.offsetHeight;
  p_w = popup.offsetWidth;
  p_h = popup.offsetHeight;
  page_w = null;
  if (window.innerWidth) {
    page_w = window.innerWidth;
  } else if (document.body) {
    page_w = document.body.clientWidth;
  }
  // check the position type is defined
  if (typeof position !== "number") {
    position = 0;
  }
  // calculate the popup position
  switch (position) {
    case 1:
      p_x = a_x + a_w + spacer;
      p_y = a_y + (a_h / 2) - (p_h / 2);
      break;
    case 0:
    default:
      p_x = a_x + (a_w / 2) - (p_w / 2);
      p_y = a_y + a_h + spacer;
      break;
  }
  // constrain the popup position
  if (p_x < margin) {
    p_x = margin;
  } else if (page_w != null && (p_x + p_w) > (page_w - margin - scrollbar)) {
    p_x = page_w - margin - scrollbar - p_w;
  }
  if (p_y < margin) {
    p_y = margin;
  }
  // position the popup
  popup.style.left = p_x + "px";
  popup.style.top = p_y + "px";
}

function lookup_help_popup(popup_id) {
  var _body, pop, info;
  pop = document.getElementById(popup_id);
  if (pop == null) {
    _body = document.getElementsByTagName("body")[0];
    pop = document.createElement("div");
    pop.className = "pop_content";
    pop.id = popup_id;
    pop.style.backgroundColor = "#FFC";
    pop.style.borderColor = "black";
    info = document.createElement("p");
    info.style.fontWeight = "bold";
    info.appendChild(document.createTextNode("Error: No popup for topic \"" + popup_id + "\"."));
    pop.appendChild(info);
    // this might cause problems with the menu, but as this only happens
    // when something is already wrong I don't think that's too much of a problem
    _body.insertBefore(pop, _body.firstChild);
  }
  if (document.getElementsByTagName('body')[0].hasAttribute("data-autobtns")) {
    if (!/\bauto_buttons\b/.test(pop.className)) {
      pop.className += " auto_buttons";
      var back_btn_sec = document.createElement("div");
      back_btn_sec.className = "nested_only pop_back_sec";
      var back_btn = document.createElement("span");
      back_btn.className = "pop_back";
      back_btn.appendChild(document.createTextNode("<< back"));
      back_btn.addEventListener("click", function(e) {
        help_return();
      }, false);
      back_btn_sec.appendChild(back_btn);
      pop.insertBefore(back_btn_sec, pop.firstChild);
      var close_btn_sec = document.createElement("div");
      close_btn_sec.className = "pop_close_sec";
      var close_btn = document.createElement("span");
      close_btn.className = "pop_close";
      close_btn.appendChild(document.createTextNode("close"));
      close_btn.addEventListener("click", function(e) {
        help_popup();
      }, false);
      close_btn_sec.appendChild(close_btn);
      pop.appendChild(close_btn_sec);
    }
  }
  return pop;
}

/*
 * help_popup
 *
 * Moves around help pop-ups so they appear
 * below an activator.
 */
function help_popup(activator, popup_id) {
  "use strict";
  var pop;
  // set default values
  if (typeof help_popup.popup === "undefined") {
    help_popup.popup = [];
  }
  if (typeof help_popup.activator === "undefined") {
    help_popup.activator = null;
  }
  var last_pop = (help_popup.popup.length > 0 ? help_popup.popup[help_popup.popup.length - 1] : null);
  if (typeof(activator) == "undefined") { // no activator so hide
    if (last_pop != null) {
      last_pop.style.display = 'none';
      help_popup.popup = [];
    }
    return;
  }
  pop = lookup_help_popup(popup_id);
  if (pop == last_pop) {
    if (activator == help_popup.activator) {
      //hide popup (as we've already shown it for the current help button)
      last_pop.style.display = 'none';
      help_popup.popup = [];
      return; // toggling complete!
    }
  } else if (last_pop != null) {
    //activating different popup so hide current one
    last_pop.style.display = 'none';
  }
  help_popup.popup = [pop];
  help_popup.activator = activator;
  toggle_class(pop, "nested", false);
  //must make the popup visible to measure it or it has zero width
  pop.style.display = 'block';
  position_popup(activator, pop);
}

/*
 * help_refine
 * 
 * Intended for links within a help popup. Stores a stack of state so
 * you can go back.
 */
function help_refine(popup_id) {
  if (help_popup.popup == null || help_popup.popup.length == 0 || help_popup.activator == null) {
    throw new Error("Can not refine a help popup when one is not shown!");
  }
  var pop = lookup_help_popup(popup_id);
  var last_pop = help_popup.popup[help_popup.popup.length - 1];
  if (pop == last_pop) return; // slightly odd, but no real cause for alarm
  help_popup.popup.push(pop);
  toggle_class(pop, "nested", true);
  last_pop.style.display = "none";
  //must make the popup visible to measure it or it has zero width
  pop.style.display = "block";
  position_popup(help_popup.activator, pop);
}

/*
 * help_return
 * 
 * Intended for links within a help popup. Stores a stack of state so
 * you can go back.
 */
function help_return() {
  if (help_popup.popup == null || help_popup.popup.length == 0 || help_popup.activator == null) {
    throw new Error("Can not return to a earlier help popup when one is not shown!");
  }
  var last_pop = help_popup.popup.pop();
  last_pop.style.display = "none";
  var pop = (help_popup.popup.length > 0 ? help_popup.popup[help_popup.popup.length - 1] : null);
  if (pop != null) {
    toggle_class(pop, "nested", help_popup.popup.length > 1);
    pop.style.display = "block";
    position_popup(help_popup.activator, pop);
  } else {
    help_popup.activator = null;
  }
}

/*
 * update_scroll_pad
 *
 * Creates padding at the bottom of the page to allow
 * scrolling of anything into view.
 */
function update_scroll_pad() {
  var page, pad;
  page = (document.compatMode === "CSS1Compat") ? document.documentElement : document.body;
  pad = $("scrollpad");
  if (pad === null) {
    pad = document.createElement("div");
    pad.id = 'scrollpad';
    document.getElementsByTagName('body')[0].appendChild(pad);
  }
  pad.style.height = Math.abs(page.clientHeight - 100) + "px";
}

function substitute_classes(node, remove, add) {
  "use strict";
  var list, all, i, cls, classes;
  list = node.className.split(/\s+/);
  all = {};
  for (i = 0; i < list.length; i++) {
    if (list[i].length > 0) all[list[i]] = true;
  }
  for (i = 0; i < remove.length; i++) {
    if (all.hasOwnProperty(remove[i])) {
      delete all[remove[i]];
    }
  }
  for (i = 0; i < add.length; i++) {
    all[add[i]] = true;
  }
  classes = "";
  for (cls in all) {
    classes += cls + " ";
  }
  node.className = classes;
}

/*
 * toggle_class
 *
 * Adds or removes a class from the node. If the parameter 'enabled' is not 
 * passed then the existence of the class will be toggled, otherwise it will be
 * included if enabled is true.
 */
function toggle_class(node, cls, enabled) {
  var classes = node.className;
  var list = classes.replace(/^\s+/, '').replace(/\s+$/, '').split(/\s+/);
  var found = false;
  for (var i = 0; i < list.length; i++) {
    if (list[i] == cls) {
      list.splice(i, 1);
      i--;
      found = true;
    }
  }
  if (typeof enabled == "undefined") {
    if (!found) list.push(cls);
  } else {
    if (enabled) list.push(cls);
  }
  node.className = list.join(" ");
}

/*
 * find_child
 *
 * Searches child nodes in depth first order and returns the first it finds
 * with the className specified.
 * TODO replace with querySelector
 */
function find_child(node, className) {
  var pattern;
  if (node == null || typeof node !== "object") {
    return null;
  }
  if (typeof className === "string") {
    pattern = new RegExp("\\b" + className + "\\b");
  } else {
    pattern = className;
  }
  if (node.nodeType == Node.ELEMENT_NODE && 
      pattern.test(node.className)) {
    return node;
  } else {
    var result = null;
    for (var i = 0; i < node.childNodes.length; i++) {
      result = find_child(node.childNodes[i], pattern);
      if (result != null) break;
    }
    return result;
  }
}

/*
 * find_parent
 *
 * Searches parent nodes outwards from the node and returns the first it finds
 * with the className specified.
 */
function find_parent(node, className) {
  var pattern;
  pattern = new RegExp("\\b" + className + "\\b");
  do {
    if (node.nodeType == Node.ELEMENT_NODE && 
        pattern.test(node.className)) {
      return node;
    }
  } while (node = node.parentNode);
  return null;
}

/*
 * find_parent_tag
 *
 * Searches parent nodes outwards from the node and returns the first it finds
 * with the tag name specified. HTML tags should be specified in upper case.
 */
function find_parent_tag(node, tag_name) {
  do {
    if (node.nodeType == Node.ELEMENT_NODE && node.tagName == tag_name) {
      return node;
    }
  } while (node = node.parentNode);
  return null;
}

/*
 * __toggle_help
 *
 * Uses the 'topic' property of the this object to
 * toggle display of a help topic.
 *
 * This function is not intended to be called directly.
 */
function __toggle_help(e) {
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }

  help_popup(this, this.getAttribute("data-topic"));
}

function setup_help_button(button) {
  "use strict";
  var topic;
  if (button.hasAttribute("data-topic")) {
    topic = button.getAttribute("data-topic");
    if (document.getElementById(topic) != null) {
      button.tabIndex = "0"; // make keyboard selectable
      button.addEventListener("click", function() {
        help_popup(button, topic);
      }, false);
      button.addEventListener("keydown", function(e) {
        // toggle only on Enter or Spacebar, let other keys do their thing
        if (e.keyCode !== 13 && e.keyCode !== 32) return;
        // stop a submit or something like that
        e.preventDefault();
        help_popup(button, topic);
      }, false);
    } else {
      button.style.visibility = "hidden";
    }
  }
  button.className += " active";
}

/*
 * help_button
 *
 * Makes a help button for the passed topic.
 */
function help_button(topic) {
  var btn = document.createElement("div");
  btn.className = "help";
  btn.setAttribute("data-topic", topic);
  setup_help_button(btn);
  return btn;
}

/*
 * prepare_download
 *
 * Sets the attributes of a link to setup a file download using the given content.
 * If no link is provided then create one and click it.
 */
function prepare_download(content, mimetype, filename, link) {
  "use strict";
  // if no link is provided then create one and click it
  var click_link = false;
  if (!link) {
    link = document.createElement("a");
    click_link = true;
  }
  try {
    // Use a BLOB to convert the text into a data URL.
    // We could do this manually with a base 64 conversion.
    // This will only be supported on modern browsers,
    // hence the try block.
    var blob = new Blob([content], {type: mimetype});
    var reader = new FileReader();
    reader.onloadend = function() {
      // If we're lucky the browser will also support the download
      // attribute which will let us suggest a file name to save the link.
      // Otherwise it is likely that the filename will be unintelligible. 
      link.setAttribute("download", filename);
      link.href = reader.result;
      if (click_link) {
        // must add the link to click it
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
      }
    }
    reader.readAsDataURL(blob);
  } catch (error) {
    if (console && console.log) console.log(error);
    // probably an old browser
    link.href = "";
    link.visible = false;
  }
}

/*
 * add_cell
 *
 * Add a cell to the table row.
 */
function add_cell(row, node, cls, click_action) {
  var cell = row.insertCell(row.cells.length);
  if (node) cell.appendChild(node);
  if (cls && cls !== "") cell.className = cls;
  if (click_action) cell.addEventListener("click", click_action, false);
}

/*
 * add_header_cell
 *
 * Add a header cell to the table row.
 */
function add_header_cell(row, node, help_topic, cls, colspan) {
  var th = document.createElement("th");
  if (node) th.appendChild(node);
  if (help_topic && help_topic !== "") th.appendChild(help_button(help_topic));
  if (cls && cls !== "") th.className = cls;
  if (typeof colspan == "number" && colspan > 1) th.colSpan = colspan;
  row.appendChild(th);
}

/*
 * add_text_cell
 *
 * Add a text cell to the table row.
 */
function add_text_cell(row, text, cls, action) {
  var node = null;
  if (typeof(text) != 'undefined') node = document.createTextNode(text);
  add_cell(row, node, cls, action);
}

/*
 * add_text_header_cell
 *
 * Add a text header cell to the table row.
 */
function add_text_header_cell(row, text, help_topic, cls, action, colspan) {
  var node = null;
  if (typeof(text) != 'undefined') {
    var nbsp = (help_topic ? "\u00A0" : "");
    var str = "" + text;
    var parts = str.split(/\n/);
    if (parts.length === 1) {
      if (action) {
        node = document.createElement("span");
        node.appendChild(document.createTextNode(str + nbsp));
      } else {
        node = document.createTextNode(str + nbsp);
      }
    } else {
      node = document.createElement("span");
      for (var i = 0; i < parts.length; i++) {
        if (i !== 0) {
          node.appendChild(document.createElement("br"));
        }
        node.appendChild(document.createTextNode(parts[i]));
      }
    }
    if (action) {
      node.addEventListener("click", action, false);
      node.style.cursor = "pointer";
    }
  }
  add_header_cell(row, node, help_topic, cls, colspan);
}

function setup_help() {
  "use strict";
  var help_buttons, i;
  help_buttons = document.querySelectorAll(".help:not(.active)");
  for (i = 0; i < help_buttons.length; i++) {
    setup_help_button(help_buttons[i]);
  }
}

function setup_scrollpad() {
  "use strict";
  if (document.getElementsByTagName('body')[0].hasAttribute("data-scrollpad") && document.getElementById("scrollpad") == null) {
    window.addEventListener("resize", update_scroll_pad, false);
    update_scroll_pad();
  }
}

// anon function to avoid polluting global scope
(function() {
  "use strict";
  window.addEventListener("load", function load(evt) {
    window.removeEventListener("load", load, false);
    setup_help();
    setup_scrollpad();
  }, false);
})();

/*
 *  make_link
 *
 *  Creates a text node and if a URL is specified it surrounds it with a link.
 *  If the URL doesn't begin with "http://" it automatically adds it, as
 *  relative links don't make much sense in this context.
 */
function make_link(text, url) {
  var textNode = null;
  var link = null;
  if (typeof text !== "undefined" && text !== null) textNode = document.createTextNode(text);
  if (typeof url === "string") {
    if (url.indexOf("//") == -1) {
      url = "http://" + url;
    }
    link = document.createElement('a');
    link.href = url;
    if (textNode) link.appendChild(textNode);
    return link;
  }
  return textNode;
}
