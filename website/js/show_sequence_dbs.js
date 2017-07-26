"use strict";

/*****************************************************************************
 * 
 ****************************************************************************/
function make_priors_info(data) {
  var box, desc, table, tbody, tr, th, td, link;
  box = document.createElement("div");
  box.className = "priors_info";
  table = document.createElement("table");
  tbody = document.createElement("tbody");
  tr = document.createElement("tr");
  th = document.createElement("th");
  th.appendChild(document.createTextNode("Sample"));
  tr.appendChild(th);
  td = document.createElement("td");
  td.appendChild(document.createTextNode(data.biosample));
  tr.appendChild(td);
  tbody.appendChild(tr);
  tr = document.createElement("tr");
  th = document.createElement("th");
  th.appendChild(document.createTextNode("Assay"));
  tr.appendChild(th);
  td = document.createElement("td");
  td.appendChild(document.createTextNode(data.assay));
  tr.appendChild(td);
  tbody.appendChild(tr);
  tr = document.createElement("tr");
  th = document.createElement("th");
  th.appendChild(document.createTextNode("Source"));
  tr.appendChild(th);
  td = document.createElement("td");
  td.appendChild(document.createTextNode(data.source));
  tr.appendChild(td);
  tbody.appendChild(tr);
  tr = document.createElement("tr");
  th = document.createElement("th");
  th.appendChild(document.createTextNode("URL"));
  tr.appendChild(th);
  td = document.createElement("td");
  link = document.createElement("a");
  link.href = data.url;
  link.appendChild(document.createTextNode(data.url));
  td.appendChild(link);
  tr.appendChild(td);
  tbody.appendChild(tr);
  table.appendChild(tbody);
  box.appendChild(table);
  desc = document.createElement("div");
  desc.innerHTML = data.description;
  box.appendChild(desc);
  return box;
}


/*****************************************************************************
 * 
 ****************************************************************************/
function make_priors_section(container, data) {
  var hdr, i;
  container.innerHTML = "";
  if (data.priors.length > 0) {
    hdr = document.createElement("h4");
    hdr.appendChild(document.createTextNode("Available Priors"));
    container.appendChild(hdr);
    for (i = 0; i < data.priors.length; i++) {
      container.appendChild(make_priors_info(data.priors[i]));
    }
  }
}

/*****************************************************************************
 * 
 ****************************************************************************/
function load_priors(container, sequence_id) {
  "use strict";
  // now send the request
  var url = "sequences?sequence=" + sequence_id;
  var request = new XMLHttpRequest();
  request.addEventListener("load", function(evt) {
    make_priors_section(container, JSON.parse(request.responseText));
  }, false);
  request.open("GET", url, true);
  request.send();
}

/*****************************************************************************
 * Called when an alphabet or version button is clicked.
 * This checks to see that the currently selected opposite button (ie alphabet
 * for the version or vice-versa) can remain selected and if not it selects
 * one that can. It also toggles the look of the buttons to highlight the
 * selected one and gray out any options that aren't possible in combination
 * with the selected opposite.
 ****************************************************************************/
function update_selected_db(container, button, data) {
  var i, version, alphabet, active;
  if (button == null) return;
  var filter = document.getElementById("filter_priors").checked;
  // make sure the button we clicked is the active one
  if (!/\bactive\b/.test(button.className)) {
    // find existing active button
    active = container.querySelector("div.button.active." + (/\bversion\b/.test(button.className) ? "version" : "alphabet"));
    // deactivate existing active button
    if (active != null) toggle_class(active, "active", false);
    // make our button active
    toggle_class(button, "active", true);
  }
  // get a list of all the buttons
  var version_buttons = container.querySelectorAll("div.button.version");
  var alphabet_buttons = container.querySelectorAll("div.button.alphabet");
  // test if the currently selected opposite button can stay selected
  if (/\bversion\b/.test(button.className)) {
    version = data.versions[button.getAttribute("data-version")];
    // see if we can use the currently selected alphabet
    var active = container.querySelector("div.button.alphabet.active");
    if (active != null) {
      alphabet = active.getAttribute("data-alphabet");
      if (version.sequences[alphabet] == null || (filter && version.sequences[alphabet].priorCount == 0)) {
        // can not use this alphabet, so deselect it
        toggle_class(active, "active", false);
        active = null;
      }
    }
    if (active == null) {
      // try the other buttons
      for (i = 0; i < alphabet_buttons.length; i++) {
        active = alphabet_buttons[i];
        alphabet = active.getAttribute("data-alphabet");
        if (version.sequences[alphabet] != null && (!filter || version.sequences[alphabet].priorCount > 0)) {
          toggle_class(active, "active", true);
          break;
        }
      }
    }
  } else {
    alphabet = button.getAttribute("data-alphabet");
    // see if we can use the currently selected version
    var active = container.querySelector("div.button.version.active");
    if (active != null) {
      version = data.versions[active.getAttribute("data-version")];
      if (version.sequences[alphabet] == null || (filter && version.sequences[alphabet].priorCount == 0)) {
        // can not use this alphabet, so deselect it
        toggle_class(active, "active", false);
        active = null;
      }
    }
    if (active == null) {
      // try the other buttons
      for (i = 0; i < version_buttons.length; i++) {
        active = version_buttons[i];
        version = data.versions[active.getAttribute("data-version")];
        if (version.sequences[alphabet] != null && (!filter || version.sequences[alphabet].priorCount > 0)) {
          toggle_class(active, "active", true);
          break;
        }
      }
    }
  }
  // Now go over all the buttons and grey them out if they can't be used with
  // the currently selected opposite. Note that this is a purely visual effect,
  // they can still be selected but will of course deselect the currently
  // selected opposite.
  var btn, ver, alph;
  for (i = 0; i < version_buttons.length; i++) {
    btn = version_buttons[i];
    ver = data.versions[btn.getAttribute("data-version")];
    toggle_class(btn, "deactivated", (ver.sequences[alphabet] == null || (filter && ver.sequences[alphabet].priorCount == 0)));
  }
  for (i = 0; i < alphabet_buttons.length; i++) {
    btn = alphabet_buttons[i];
    alph = btn.getAttribute("data-alphabet");
    toggle_class(btn, "deactivated", (version.sequences[alph] == null || (filter && version.sequences[alph].priorCount == 0)));
  }
  // Set the description to match the current combination of version and alphabet
  container.querySelector("div.out").innerHTML = version.sequences[alphabet].description;
  // load priors
  load_priors(container.querySelector("div.priors"), version.sequences[alphabet].id);
}

/*****************************************************************************
 * When the filter is turned on update the selected entries of any listings
 * that are still visible.
 ****************************************************************************/
function update_all_selected_dbs() {
  var i, container, data, btn;
  var filter = document.getElementById("filter_priors").checked;
  // get all the listings that haven't been hidden
  var db_btn_containers = document.querySelectorAll(".listing" + (filter ? ":not(.no_priors)" : "") + " .info");
  for (i = 0; i < db_btn_containers.length; i++) {
    container = db_btn_containers[i];
    btn = container.querySelector(".button.version.active");
    if (btn != null) {
      // check that the currently selected button can remain selected
      if (filter && /\bno_priors\b/.test(btn.className)) {
        // select an alternative button
        btn = container.querySelector(".button.version:not(.no_priors)");
      }
      data = JSON.parse(container.getAttribute("data-info"));
      // update the state of the buttons
      update_selected_db(container, btn, data);
    }
  }
}

/*****************************************************************************
 * Create the version info section.
 ****************************************************************************/
function make_version_info(container, data) {
  var i, j, version, alph, file, btn, out, hasPriors, priors;
  var make_db_handler = function(container, button, data) {
    return function (evt) {
      update_selected_db(container, button, data);
    };
  }
  // clear the container
  container.innerHTML = "";
  container.setAttribute("data-info", JSON.stringify(data));
  var overview = document.createElement("p");
  overview.innerHTML = data.description;
  container.appendChild(overview);
  var hdr = document.createElement("h5");
  hdr.appendChild(document.createTextNode("Versions: "));
  container.appendChild(hdr);
  for (i = 0; i < data.versions.length; i++) {
    version = data.versions[i];
    hasPriors = false;
    for (j = 0; j < data.alphabets.length; j++) {
      if (version.sequences[data.alphabets[j]] != null && version.sequences[data.alphabets[j]].priorCount > 0) {
        hasPriors = true;
        break;
      }
    }
    btn = document.createElement("div");
    btn.className = "button version" + (hasPriors ? "" : " no_priors");
    btn.setAttribute("data-version", i);
    btn.appendChild(document.createTextNode(version.name));
    btn.addEventListener("click", make_db_handler(container, btn, data), false);
    container.appendChild(btn);
  }
  hdr = document.createElement("h5");
  hdr.appendChild(document.createTextNode("\u00A0\u00A0Alphabets: "));
  container.appendChild(hdr);
  var has_active = false;
  for (i = 0; i < data.alphabets.length; i++) {
    hasPriors = false;
    for (j = 0; j < data.versions.length; j++) {
      version = data.versions[j];
      if (version.sequences[data.alphabets[i]] != null && version.sequences[data.alphabets[i]].priorCount > 0) {
        hasPriors = true;
        break;
      }
    }
    alph = AlphStd[data.alphabets[i]];
    btn = document.createElement("div");
    btn.className = "button alphabet" + (hasPriors ? "" : " no_priors");
    btn.setAttribute("data-alphabet", data.alphabets[i]);
    btn.appendChild(document.createTextNode(alph.get_alphabet_name()));
    btn.addEventListener("click", make_db_handler(container, btn, data), false);
    container.appendChild(btn);
  }
  out = document.createElement("div");
  out.className = "out";
  container.appendChild(out);
  // create the priors section
  priors = document.createElement("div");
  priors.className = "priors";
  container.appendChild(priors);
  update_selected_db(container, container.querySelector("div.button.version" + (hasPriors ? ":not(.no_priors)" : "")), data);
}

/*****************************************************************************
 * Called when the versions of a listing need to be loaded.
 ****************************************************************************/
function load_versions(container, listing_id) {
  "use strict";
  // now send the request
  var url = "sequences?listing=" + listing_id;
  var request = new XMLHttpRequest();
  request.addEventListener("load", function(evt) {
    make_version_info(container, JSON.parse(request.responseText));
  }, false);
  request.open("GET", url, true);
  request.send();
}

/*****************************************************************************
 * Called when a listing is clicked.
 * This will expand/collapse the listing.
 * If it hasn't been loaded already it will do so.
 ****************************************************************************/
function toggle_listing(listing, button) {
  toggle_class(button, 'expanded');
  if (!/\bloading\b/.test(listing.className)) {
    listing.className += " loading";
    var container = listing.querySelector(".info");
    load_versions(container, listing.getAttribute("data-id"));
  }
}

/*****************************************************************************
 * Creates a listing.
 ****************************************************************************/
function create_listing(id, name, hasPriors) {
  var listing = document.createElement("div");
  listing.className = "listing" + (hasPriors ? "" : " no_priors");
  listing.setAttribute("data-id", id);
  var button = document.createElement("div");
  button.className = "btn";
  var heading = document.createElement("h4");
  heading.appendChild(document.createTextNode(name));
  button.appendChild(heading);
  button.appendChild(document.createTextNode("\u2002"));
  var more = document.createElement("span");
  more.className = "collapsed";
  more.appendChild(document.createTextNode("..."));
  button.appendChild(more);
  var less = document.createElement("span");
  less.className = "expanded";
  less.appendChild(document.createTextNode("\u25BC"));
  button.appendChild(less);
  listing.appendChild(button);
  var info = document.createElement("div");
  info.className = "info subcontent";
  info.appendChild(document.createTextNode("loading..."));
  listing.appendChild(info);
  button.addEventListener("click", function (evt) {
    toggle_listing(listing, button);
  }, false);
  button.addEventListener("keypress", function(evt) {
    if (evt.which == 32 || evt.keyCode == 32) {
      toggle_listing(listing, button);
    }
  }, false);
  return listing;
}

/*****************************************************************************
 * Called when the listings of a category need to be loaded after the
 * category has been expanded.
 *
 ****************************************************************************/
function load_listings(container, category_id) {
  "use strict";
  // now send the request
  var url = "sequences?category=" + category_id;
  var request = new XMLHttpRequest();
  request.addEventListener("load", function(evt) {
    var listings, listing, i;
    var data = JSON.parse(request.responseText);
    // clear the container
    container.innerHTML = "";
    // add the other options
    for (i = 0; i < data.listings.length; i++) {
      listing = data.listings[i];
      container.appendChild(create_listing(listing.id, listing.name, listing.hasPriors));
    }
    // re-enable the list
  }, false);
  request.open("GET", url, true);
  request.send();

}

/*****************************************************************************
 * Called when a category is clicked.
 * This loads the category or if it is already loaded it just toggles the
 * expanded state.
 ****************************************************************************/
function toggle_category(category, button) {
  toggle_class(button, 'expanded');
  if (!/\bloading\b/.test(category.className)) {
    category.className += " loading";
    var container = category.querySelector(".info");
    load_listings(container, category.getAttribute("data-id"));
  }
}


