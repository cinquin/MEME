
/******************************************************************************
 * Create a FASTA indicator handler.
 ******************************************************************************/
var FIH = function (container, indicator, popup, options, change_handler) {
  "use strict";
  this.container = container;
  this.indicator = indicator;
  this.popup = popup;
  this.change_handler = change_handler;
  this.update_timer = null;
  this.display_error = false;
  this.display_warn = false;
  FastaHandler.call(this, options);
};
FIH.prototype = Object.create(FastaHandler.prototype);
FIH.prototype.constructor = FIH;


/******************************************************************************
 * Clears all the indicators.
 ******************************************************************************/
FIH.prototype.clear_indicators = function() {
  "use strict";
  if (this.update_timer) {
    clearInterval(this.update_timer);
    this.update_timer = null;
  }
  this.popup.innerHTML = "";
  this.indicator.style.width = "0%";
  substitute_classes(this.container, ["good", "error", "warning"], []);
  this.display_type = 0;
};

/******************************************************************************
 *
 ******************************************************************************/
FIH.prototype._update_display = function () {
  "use strict";
  var table, summary, i, j, msg;
  var row, cell, i, ul, li;
  
  this.indicator.style.width = (this.fraction * 100) + "%";
  if (!this.updated) return;

  summary = this.summary();
  if (summary.error) {
    if (this.display_type != 3) {
      substitute_classes(this.container, ["good", "warning"], ["error"]); 
      this.display_type = 3;
    }
  } else if (summary.warning) {
    if (this.display_type != 2) {
      substitute_classes(this.container, ["good", "error"], ["warning"]);
      this.display_type = 2;
    }
  } else {
    if (this.display_type != 1) {
      substitute_classes(this.container, ["warning", "error"], ["good"]);
      this.display_type = 1;
    }
  }
  table = document.createElement("table");
  for (i = 0; i < summary.messages.length; i++) {
    msg = summary.messages[i];
    row = table.insertRow(table.rows.length);
    row.className = (msg.is_error ? "error" : "warning");
    cell = row.insertCell(row.cells.length);
    cell.appendChild(document.createTextNode(msg.is_error ? "\u2718" : "\u26A0"));
    cell = row.insertCell(row.cells.length);
    cell.appendChild(document.createTextNode(msg.message));
    if (typeof msg.reasons !== "undefined") {
      ul = document.createElement("ul");
      for (j = 0; j < msg.reasons.length; j++) {
        li = document.createElement("li");
        li.appendChild(document.createTextNode(msg.reasons[j]));
        ul.appendChild(li);
      }
      cell.appendChild(ul);
    }
  }
  
  this.popup.innerHTML = "";
  this.popup.appendChild(table);

};

/******************************************************************************
 *
 ******************************************************************************/
FIH.prototype.reset = function () {
  "use strict";
  FastaHandler.prototype.reset.apply(this, arguments);
  this.clear_indicators();
};

/******************************************************************************
 *
 ******************************************************************************/
FIH.prototype.begin = function () {
  "use strict";
  var me;
  FastaHandler.prototype.begin.apply(this, arguments);
  if (this.update_timer) clearInterval(this.update_timer);
  me = this; // reference 'this' inside closure
  this.update_timer = setInterval(function () { me._update_display();}, 100);
};

/******************************************************************************
 *
 ******************************************************************************/
FIH.prototype.end = function () {
  "use strict";
  FastaHandler.prototype.end.apply(this, arguments);
  if (this.update_timer) {
    clearInterval(this.update_timer);
    this.update_timer = null;
  }
  this._update_display();
  if (this.display_type == 1) this.clear_indicators();
  this.change_handler();
};

/******************************************************************************
 * Construct the manager of the alphabet input.
 ******************************************************************************/
var AlphabetInput = function(container) {
  "use strict";
  var me;
  // make 'this' accessable in inner scopes
  me = this;
  // get the parts
  this.container = container;
  // get radio buttons
  this.radio_custom_on = this.container.querySelector("input[type='radio'][value='1']");
  this.radio_custom_off = this.container.querySelector("input[type='radio'][value='0']");
  this.last_radio = null;
  // get the custom alphabet related parts
  this.cust_surround = this.container.querySelector("span.custom_alphabet");
  this.cust_input = this.cust_surround.querySelector("input");
  this.cust_popup = this.cust_surround.querySelector("div.popup");
  this.cust_error = false;
  this.cust_alphabet = null;
  // get the area used to display the alphabet
  this.alph_info = this.container.querySelector("span.alphabet_info");
  // set state
  this.cust_input.disabled = !this.radio_custom_on.checked;
  toggle_class(this.container, "custom", this.radio_custom_on.checked);
  this._alphabet_update();
  // detect radio button changes
  this.radio_custom_on.addEventListener("click", function(e) { me._radio_update(this); }, false);
  this.radio_custom_off.addEventListener("click", function(e) { me._radio_update(this); }, false);
  // detect changes in the custom alphabet
  this.cust_input.addEventListener("change", function(e) { me._alphabet_update() }, false);
  // detect form resets and reset properly
  if (this.radio_custom_on.form != null) {
    this.radio_custom_on.form.addEventListener("reset", function() {
      window.setTimeout(function () {
        me.reset();
      }, 50);
    }, false);
  }
};

/******************************************************************************
 * Resets the alphabet input to the default state.
 ******************************************************************************/
AlphabetInput.prototype.reset = function() {
  this.cust_input.value = "";
  this.cust_popup.innerHTML = "";
  substitute_classes(this.cust_surround, ["good", "error", "warning"], []);
  this.cust_error = false;
  this.cust_alphabet = null;
  this.radio_custom_on.checked = false;
  this.radio_custom_off.checked = true;
  this.cust_input.disabled = true;
  toggle_class(this.container, "custom", false);
  this._fire_alphabet_event();
};

/******************************************************************************
 * Returns the alphabet loaded from the file.
 ******************************************************************************/
AlphabetInput.prototype.get_custom_alphabet = function() {
  return this.cust_alphabet;
};

/******************************************************************************
 * Returns true if the option of a custom alphabet has been selected.
 ******************************************************************************/
AlphabetInput.prototype.has_custom_alphabet = function() {
  return this.radio_custom_on.checked;
};

/******************************************************************************
 * Returns the custom alphabet in a single element list or null.
 ******************************************************************************/
AlphabetInput.prototype.get_alphabets = function() {
  if (this.radio_custom_on.checked && this.cust_alphabet != null) {
    return [this.cust_alphabet];
  }
  return null;
};

/******************************************************************************
 * Fires an alphabet_changed event with the current status.
 ******************************************************************************/
AlphabetInput.prototype._fire_alphabet_event = function() {
  "use strict";
  try {
    // IE sometimes has problems with these lines.
    // I think they are related to the page not being fully loaded.
    var data = {"has_custom": this.radio_custom_on.checked, "alphabet": this.cust_alphabet};
    this.container.dispatchEvent(
        new CustomEvent("alphabet_changed", {detail: data}));
  } catch (e) {
    if (e.message && e.name && window.console) {
      console.log("Suppressed exception " + e.name + ": " + e.message);
    }
  }
};

/******************************************************************************
 * Displays the name of an alphabet truncated to 5em width with the details
 * shown on mouseover within a popup.
 ******************************************************************************/
AlphabetInput.prototype._show_alph_name = function() {
  "use strict";
  this.alph_info.innerHTML = "";
  var alphabet = this.cust_alphabet;
  if (alphabet != null) {
    var name = document.createElement("div");
    name.className = "alph_name";
    name.appendChild(document.createTextNode(alphabet.toString()));
    var title = document.createElement("h4");
    title.appendChild(document.createTextNode(alphabet.toString()));
    var popup = document.createElement("div");
    popup.className = "popup";
    popup.appendChild(title);
    popup.appendChild(alphabet.as_table());
    var info = document.createElement("div");
    info.className = "alph_info";
    info.appendChild(name);
    info.appendChild(popup);
    this.alph_info.appendChild(info);
  }
}

/******************************************************************************
 * CALLBACK HANDLER
 * Called when a problem is found with an alphabet file being parsed.
 * Creates an entry in the popup used to display warnings and errors.
 * Updates the error status.
 ******************************************************************************/
AlphabetInput.prototype._parse_error = function(is_error, message, reasons) {
  "use strict";
  var table, row, cell, i, ul, li;
  table = this.cust_popup.querySelector("table");
  if (table == null) {
    table = document.createElement("table");
    this.cust_popup.innerHTML = "";
    this.cust_popup.appendChild(table);
  }
  row = table.insertRow(table.rows.length);
  row.className = (is_error ? "error" : "warning");
  cell = row.insertCell(row.cells.length);
  cell.appendChild(document.createTextNode(is_error ? "\u2718" : "\u26A0"));
  cell = row.insertCell(row.cells.length);
  cell.appendChild(document.createTextNode(message));
  if (typeof reasons !== "undefined") {
    ul = document.createElement("ul");
    for (i = 0; i < reasons.length; i++) {
      li = document.createElement("li");
      li.appendChild(document.createTextNode(reasons[i]));
      ul.appendChild(li);
    }
    cell.appendChild(ul);
  }
  this.cust_error |= is_error;
  if (this.cust_error) {
    substitute_classes(this.cust_surround, ["good", "warning"], ["error"]);
  } else {
    substitute_classes(this.cust_surround, ["good"], ["warning"]);
  }
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called when an alphabet file has been successfully parsed.
 * Stores the alphabet and fires an alphabet event to notify listeners.
 ******************************************************************************/
AlphabetInput.prototype._parse_data = function(alphabet_data) {
  this.cust_alphabet = new Alphabet(alphabet_data);
  substitute_classes(this.cust_surround, ["good"], []);
  this._show_alph_name();
  this._fire_alphabet_event();
};

/******************************************************************************
 * EVENT HANDLER
 * Fires when the custom alphabet file is changed.
 * Parses any supplied alphabet file.
 ******************************************************************************/
AlphabetInput.prototype._alphabet_update = function() {
  var file, me;
  // reset before parsing
  this.cust_popup.innerHTML = "";
  this.cust_error = false;
  this.cust_alphabet = null;
  this._show_alph_name();
  // check for a file to load
  if (!(file = this.cust_input.files[0])) {
    // no file to load! Clear the status. 
    substitute_classes(this.cust_surround, ["good", "warning", "error"], []);
    this._fire_alphabet_event();
  } else {
    // set the status to "good"
    substitute_classes(this.cust_surround, ["warning", "error"], ["good"]);
    // start parsing
    me = this;
    var parser = new AlphabetParser({
      "error": function(is_error, message, reasons) {
        me._parse_error(is_error, message, reasons);
      },
      "data": function(alphabet_data) {
        me._parse_data(alphabet_data);
      }
    });
    parser.process_blob(file);
  }
};

/******************************************************************************
 * EVENT HANDLER
 * Fires when an alphabet source option is clicked.
 * Detects when the click actually changed the source between a standard
 * alphabet and a custom alphabet.
 * Ensures the alphabet file input is only enabled for custom alphabets.
 ******************************************************************************/
AlphabetInput.prototype._radio_update = function(radio_button) {
  if (this.last_radio !== radio_button) {
    this._fire_alphabet_event();
  }
  this.last_radio = radio_button;
  this.cust_input.disabled = !this.radio_custom_on.checked;
  toggle_class(this.container, "custom", this.radio_custom_on.checked);
};

/******************************************************************************
 * Construct the manager of the sequence input.
 ******************************************************************************/
var SequenceInput = function(container, options) {
  "use strict";
  var me, box;
  // make 'this' accessable in inner scopes
  me = this;
  // default fasta options
  // store the parameters
  this.container = container;
  this.options = options;
  this.alph_type = AlphType[options.alph_type];
  if (this.alph_type == null) throw new Error("Unknown alphabet type");
  this.custom_alphabet = null;
  this.expected_atype = null;
  // for making XMLHttpRequests
  this._xml_alph = this.alph_type.bitset();
  this._xml_short = (options.short_only ? 1 : 0);
  // lookup relevent components
  this.prior_filter = this.container.querySelector("input.prior_filter");
  this.source = this.container.querySelector("select.sequence_source");
  // get the text related parts
  this.text_surround = this.container.querySelector("div.sequence_text");
  this.text_indicator = this.text_surround.querySelector("span.indicator");
  this.text_control = this.text_surround.querySelector(".editor");
  this.text_area = this.text_control.querySelector("textarea");
  this.text_backdrop = this.text_control.querySelector("span");
  this.text_popup = this.text_surround.querySelector("div.popup");
  this.text_dbh = new FIH(this.text_surround, this.text_indicator, this.text_popup, options, function() {
    me._fire_sequences_checked_event();
  });
  // get the file related parts
  this.file_surround = this.container.querySelector("span.sequence_file");
  this.file_indicator = this.file_surround.querySelector("span.indicator");
  this.file_input = this.file_surround.querySelector("input");
  this.file_popup = this.file_surround.querySelector("div.popup");
  this.file_dbh = new FIH(this.file_surround, this.file_indicator, this.file_popup, options, function() {
    me._fire_sequences_checked_event();
  });
  // get the database related parts
  this.db_listing = container.querySelector("select.sequence_db.listing");
  this.db_version = container.querySelector("select.sequence_db.version");
  this.db_priors = container.querySelector("select.sequence_db.priors");
  // get the embed related parts
  this.embed_surround = this.container.querySelector("span.sequence_embed");
  if (this.embed_surround != null) {
    this.embed_name = this.embed_surround.querySelector("input.name");
    this.embed_data = this.embed_surround.querySelector("input.data");
    this.embed_handler = new FastaHandler();
    // override end function to hook in our event
    this.embed_handler.end = function() {
      FastaHandler.prototype.end.apply(this, arguments);
      me._fire_sequences_checked_event();
    };
  }
  // get the area used to display the alphabet
  this.alph_info = this.container.querySelector("span.sequence_alphabet_info");
  // disable the prior filter if no priors were provided
  if (this.source.querySelector("option:not(.no_priors)") == null) {
    box = container.querySelector("div.prior_filter_section");
    if (box != null) box.style.display = "none";
  }

  // create a list of submittable fields so we can disable the ones we are not using
  this.submittables = [this.text_area, this.file_input];
  if (this.db_listing != null) this.submittables.push(this.db_listing);
  if (this.db_version != null) this.submittables.push(this.db_version);
  if (this.db_priors != null) this.submittables.push(this.db_priors);
  if (this.embed_surround != null) {
    this.submittables.push(this.embed_name);
    this.submittables.push(this.embed_data);
  }
  // other things
  this.parser = null;
  this.timer = null;
  this.db_alphabets = [];
  // parse the embeded motifs
  if (this.embed_surround != null) {
    (new FastaChecker(this.embed_handler, this.alph_type.get_standard_alphabets())).process_blob(
        new Blob([this.embed_data.value], {"type": "text\/plain"}));
  }
  // initialise
  this._priors_filter_update();
  this._source_update();
  this._text_update();
  this._file_update();
  // add listeners
  if (this.prior_filter != null) {
    this.prior_filter.addEventListener('click', function() {
      me._priors_filter_update();
      me._source_update();
    }, false);
  }
  this.source.addEventListener('change', function() {
    me._source_update();
  }, false);
  // detect typing or pasting into motif text area
  this.text_area.addEventListener('input', function() {
    me._text_update(); 
  }, false);
  // detect file changes
  this.file_input.addEventListener('change', function() {
    me._file_update();
  }, false);
  // detect db changes
  if (this.db_listing != null) {
    this.db_listing.addEventListener('change', function() {
      me._load_versions(+me.db_listing.value);
    }, false);
  }
  if (this.db_version != null) {
    this.db_version.addEventListener('change', function() {
      me._update_db_alphabets();
      me._load_priors();
    },false);
  }
  // detect form resets and reset properly
  if (this.source.form != null) {
    this.source.form.addEventListener("reset", function() {
      window.setTimeout(function () {
        me.reset();
      }, 50);
    }, false);
  }
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype.set_custom_alphabet = function(use_custom, alphabet) {
  this.custom_alphabet = (use_custom ? alphabet : null);
  this._source_update();
};

/******************************************************************************
 * This is used by the database lookup to determine which sequence is used
 * and display the priors for it. If it is not set then priors will not be
 * available.
 ******************************************************************************/
SequenceInput.prototype.set_expected_alphabet = function(alphabet) {
  var alph_names, i;
  if (alphabet != null && typeof alphabet == "object" && alphabet instanceof Alphabet) {
    alph_names = ["DNA", "RNA", "PROTEIN"];
    for (i = 0; i < alph_names.length; i++) {
      if (AlphStd[alph_names[i]].equals(alphabet)) {
        this.expected_atype = alph_names[i];
      }
    }
  } else {
    this.expected_atype = null;
  }
  this._load_priors();
}

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype.check = function(restrict_alphabets) {
  "use strict";
  var source, id, exists, handler, summary, alphabet;
  // find out what source we're using
  source = this.source.value;
  if (source == "noseq") return true;
  if (/^\d+$/.test(source)) {
    id = parseInt(source, 10);
    source = "db";
    if (restrict_alphabets != null) {
      if (!AlphabetUtil.any_equal(this.db_alphabets, restrict_alphabets)) {
        alert("Please select a " + AlphabetUtil.display(restrict_alphabets) + 
            " database for the " + this.options.field + ".");
        return false;
      }
    }
    return true;
  } else {
    if (source == "text") {
      exists = /\S/.test(this.text_area.value);
      handler = this.text_dbh;
    } else if (source == "file") {
      exists = this.file_input.value.length > 0;
      handler = this.file_dbh;
    } else if (source == "embed") {
      exists = true;
      handler = this.embed_handler;
    } else {
      throw new Error("Unknown source");
    }
    if (!exists) {
      alert("Please input " + this.options.field + ".");
      return false;
    }
    summary = handler.summary();
    if (summary.error) {
      alert("Please correct errors in the " + this.options.field + ".");
      return false;
    }
    if (summary.warning) {
      if (!confirm("There are warnings for the " + this.options.field + ". Continue anyway?")) {
        return false;
      }
    }
    if (restrict_alphabets != null) {
      if (!AlphabetUtil.any_equal(handler.guess_alphabets(), restrict_alphabets)) {
        if (!confirm("The " + this.options.field + 
              " appear to use a different alphabet to other inputs. " +
              "Continue anyway?")) {
          return false;
        }
      }
    }
    return true;
  }
};

/******************************************************************************
 *
 ******************************************************************************/
// Has this field changed from the default "reset" state.
SequenceInput.prototype.changed = function() {
  "use strict";
  var source;
  if (!this.source.options[this.source.selectedIndex].defaultSelected) return true;
  source = this.source.value;
  if (/^\d+$/.test(source)) source = "db";
  if (source == "text") {
    if (this.text_area.value.length != 0) return true;
  } else if (source == "file") {
    if (this.file_input.value.length != 0) return true;
  } else if (source == "db") {
    if (!this.db_listing.options[this.db_listing.selectedIndex].defaultSelected) return true;
    if (!this.db_version.options[this.db_version.selectedIndex].defaultSelected) return true;
  }
  return false;
};

/******************************************************************************
 *
 ******************************************************************************/
// reset the fields to a consistant state
SequenceInput.prototype.reset = function() {
  "use strict";
  var i, opt;
  if (this.parser) this.parser.cancel();
  this.parser = null;
  if (this.timer) window.clearTimeout(this.timer);
  this.timer = null;
  this.file_input.value = "";
  this.file_dbh.reset();
  this.text_area.value = "";
  this.text_dbh.reset();
  if (this.db_listing != null) this._clear_select(this.db_listing);
  if (this.db_version != null) this._clear_select(this.db_version);
  for (i = 0; i < this.source.options.length; i++) {
    opt = this.source.options[i];
    opt.selected = opt.defaultSelected;
  }
  this._source_update();
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype.get_alphabets = function() {
  "use strict";
  var source, alphabet;
  source = this.source.value;
  if (/^\d+$/.test(source)) source = "db";
  if (source == "text") {
    return this.text_dbh.guess_alphabets();
  } else if (source == "file") {
    return this.file_dbh.guess_alphabets();
  } else if (source == "db") {
    return this.db_alphabets;
  } else if (source == "embed") {
    return this.embed_handler.guess_alphabets();
  }
  return [];
};


/******************************************************************************
 * Displays the name of an alphabet truncated to 5em width with the details
 * shown on mouseover within a popup.
 ******************************************************************************/
SequenceInput.prototype._show_alph_name = function() {
  "use strict";
  var i, alphabet;
  this.alph_info.innerHTML = "";
  var alphabets = this.get_alphabets();
  for (i = 0; i < alphabets.length; i++) {
    alphabet = alphabets[i];
    var name = document.createElement("div");
    name.className = "alph_name";
    name.appendChild(document.createTextNode(alphabet.toString()));
    var title = document.createElement("h4");
    title.appendChild(document.createTextNode(alphabet.toString()));
    var popup = document.createElement("div");
    popup.className = "popup";
    popup.appendChild(title);
    popup.appendChild(alphabet.as_table());
    var info = document.createElement("div");
    info.className = "alph_info";
    info.appendChild(name);
    info.appendChild(popup);
    this.alph_info.appendChild(info);
  }
}

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._source_update = function() {
  "use strict";
  var source, opt, id, i;
  if (this.prior_filter != null && this.prior_filter.checked &&
      /\bno_priors\b/.test(this.source.options[this.source.selectedIndex].className)) {
    // we can't keep the currently selected item so scan the other database entries looking for one we can select
    opt = this.source.querySelector("optgroup.db_options option:not(.no_priors)");
    if (opt == null) {
      // apparently there are no priors, oh well fallback on some other input
      opt = this.source.querySelector("option:not(.no_priors)");
    }
    if (opt != null) {
      this.source.selectedIndex = opt.index;
    } else {
      // our only choice is to turn off the filter
      this.prior_filter.checked = false;
    }
  }
  source = this.source.value;
  if (/^\d+$/.test(source)) {
    id = parseInt(source, 10);
    source = "db";
  }
  substitute_classes(this.container, ["noseq", "text", "file", "embed", "db"], [source]);
  for (i = 0; i < this.submittables.length; i++) {
    this.submittables[i].disabled = true;
  }
  if (this.parser) this.parser.cancel();
  if (source == "text") {
    this.text_area.disabled = false;
    this._text_validate();
  } else if (source == "file") {
    this.file_input.disabled = false;
    this._file_validate();
  } else if (source == "db") {
    // the db hidden inputs will be enabled by the _load_listings command
    this._load_listings(id);
  } else if (source == "embed") {
    this.embed_name.disabled = false;
    this.embed_data.disabled = false;
  }
  this._fire_sequences_checked_event();
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._fire_sequences_checked_event = function() {
  "use strict";
  var me;
  this._show_alph_name();
  me = this;
  try {
    // IE sometimes has problems with this line.
    // I think they are related to the page not being fully loaded.
    this.container.dispatchEvent(new CustomEvent("sequences_checked", {detail: {controler: me}}));
  } catch (e) {
    if (e.message && e.name && window.console) {
      console.log("Suppressed exception " + e.name + ": " + e.message);
    }
  }
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._text_validate = function() {
  var file;
  if (this.parser) this.parser.cancel();
  if (this.text_area.value.length == 0) {
    this.text_dbh.clear_indicators();
    substitute_classes(this.container, ["rna", "dna", "protein"], []);
    return;
  }
  file = new Blob([this.text_area.value], {"type": "text\/plain"});
  this.text_dbh.configure(this.options);
  var alphabets = (this.custom_alphabet != null ? [this.custom_alphabet] : this.alph_type.get_standard_alphabets());
  this.parser = new FastaChecker(this.text_dbh, alphabets);
  this.parser.process_blob(file);
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._text_update = function() {
  "use strict";
  var nodes, me;
  nodes = document.createTextNode(this.text_area.value);
  if (this.text_backdrop.firstChild) {
    this.text_backdrop.replaceChild(nodes, this.text_backdrop.firstChild);
  } else {
    this.text_backdrop.appendChild(nodes);
  }
  if (this.source.value == "text") {
    if (this.timer) window.clearTimeout(this.timer);
    me = this;
    this.timer = window.setTimeout(function() {
      me.timer = null;
      me._text_validate();
    }, 300);
  }
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._file_validate = function() {
  var file;
  if (this.parser) this.parser.cancel();
  if (!(file = this.file_input.files[0])) {
    this.file_dbh.clear_indicators();
    substitute_classes(this.container, ["rna", "dna", "protein"], []);
    return;
  }
  this.file_dbh.configure(this.options);
  var alphabets = (this.custom_alphabet != null ? [this.custom_alphabet] : this.alph_type.get_standard_alphabets());
  this.parser = new FastaChecker(this.file_dbh, alphabets);
  this.parser.process_blob(file);
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._file_update = function() {
  "use strict";
  if (this.source.value == "file") {
    this._file_validate();
  }
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._clear_select = function(list) {
  "use strict";
  var optgroup;
  // disable until update completes
  list.disabled = true;
  optgroup = list.getElementsByTagName("optgroup")[0];
  // clear previous values
  while (optgroup.childNodes.length > 0) {
    optgroup.removeChild(optgroup.lastChild);
  }
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._load_listings = function(category_id) {
  "use strict";
  var request, url, i, me, list, optgroup, selected_listing_id;
  // so we can access in the inner function
  me = this;
  list = this.db_listing;
  selected_listing_id = (list.value != null ? parseInt(list.value, 10) : null);
  optgroup = list.getElementsByTagName("optgroup")[0];
  this._clear_select(list);
  // now send the request
  url = "../db/sequences?category=" + category_id + "&alphabets=" + this._xml_alph + "&short=" + this._xml_short;
  request = new XMLHttpRequest();
  request.addEventListener("load", function(evt) {
    var data, listing, i, opt;
    data = JSON.parse(request.responseText);
    // ensure the list is empty
    me._clear_select(list);
    // add the other options
    for (i = 0; i < data.listings.length; i++) {
      listing = data.listings[i];
      opt = new Option(listing.name, listing.id);
      opt.className = (listing.hasPriors ? "" : "no_priors");
      optgroup.appendChild(opt);
      // try to retain old selection
      if (listing.id === selected_listing_id) opt.selected = true;
    }
    // check if we can keep the selected index
    if (me.prior_filter != null && me.prior_filter.checked &&
        /\bno_priors\b/.test(list.options[list.selectedIndex].className)) {
      opt = list.querySelector("option:not(.no_priors)");
      if (opt != null) list.selectedIndex = opt.index;
    }
    // re-enable the list
    list.disabled = false;
    if (list.value) {
      me._load_versions(parseInt(list.value, 10));
    }
  }, false);
  request.open("GET", url, true);
  request.send();
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._load_versions = function(listing_id) {
  "use strict";
  var me, request, url, i, list, optgroup, selected_version_id;
  // so we can access in the inner function
  me = this;
  list = this.db_version;
  selected_version_id = (list.value != null ? parseInt(list.value, 10) : null);
  optgroup = list.getElementsByTagName("optgroup")[0];
  this._clear_select(list);
  // now send the request
  url = "../db/sequences?listing=" + listing_id + "&alphabets=" + this._xml_alph + "&short=" + this._xml_short;
  request = new XMLHttpRequest();
  request.addEventListener("load", function(evt) {
    var data, version, i, j, hasPriors, opt;
    data = JSON.parse(request.responseText);
    // ensure the list is empty
    me._clear_select(list);
    // add the other options
    for (i = 0; i < data.versions.length; i++) {
      version = data.versions[i];
      hasPriors = false;
      for (j = 0; j < version.alphabets.length; j++) {
        if (version.sequences[version.alphabets[j]].priorCount > 0) {
          hasPriors = true;
          break;
        }
      }
      opt = new Option(version.name, version.id);
      opt.setAttribute("data-alphabets", JSON.stringify(version.alphabets));
      opt.setAttribute("data-sequences", JSON.stringify(version.sequences));
      opt.className = (hasPriors ? "" : "no_priors");
      optgroup.appendChild(opt);
      // try to retain old selection
      if (version.id === selected_version_id) opt.selected = true;
    }
    // check if we can keep the selected index
    if (me.prior_filter != null && me.prior_filter.checked &&
        /\bno_priors\b/.test(list.options[list.selectedIndex].className)) {
      opt = list.querySelector("option:not(.no_priors)");
      if (opt != null) list.selectedIndex = opt.index;
    }
    // re-enable the list
    list.disabled = false;
    me._update_db_alphabets();
    me._load_priors();
  }, false);
  request.open("GET", url, true);
  request.send();
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._update_db_alphabets = function() {
  "use strict";
  var filter, opt, sequences;
  filter = this.prior_filter != null && this.prior_filter.checked;
  opt = this.db_version.options[this.db_version.selectedIndex];
  if (opt != null) {
    sequences = JSON.parse(opt.getAttribute("data-sequences"));
    this.db_alphabets = JSON.parse(opt.getAttribute("data-alphabets")).filter(
          function (value) { return !filter || sequences[value].priorCount > 0; }
        ).map(
          function (value) { return AlphStd[value]; }
        );
  } else {
    this.db_alphabets = [];
  }
  this._fire_sequences_checked_event();
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._priors_filter_update = function() {
  "use strict";
  toggle_class(this.container, "filter_priors", this.prior_filter != null && this.prior_filter.checked);
};

/******************************************************************************
 *
 ******************************************************************************/
SequenceInput.prototype._load_priors = function() {
  "use strict";
  var list, optgroup, selected_prior_id, version_opt, sequences, alph, alphs, key, sequence_id, me, url, request;
  // check that priors are enabled
  list = this.db_priors;
  if (list == null) return;
  selected_prior_id = (list.value != null ? parseInt(list.value, 10) : null);
  optgroup = list.getElementsByTagName("optgroup")[0];
  this._clear_select(list);
  // only bother loading if the field is actually displayed
  if (this.prior_filter == null || !this.prior_filter.checked) return;
  // get the currently selected version
  version_opt = this.db_version.options[this.db_version.selectedIndex];
  // assuming there was a selected version then get the sequences associated with it
  sequences = (version_opt != null ? JSON.parse(version_opt.getAttribute("data-sequences")) : {});
  alph = this.expected_atype;
  if (alph == null) {
    // see if there is a single alphabet and automatically select it if there is
    alphs = [];
    for (key in sequences) {
      if (sequences.hasOwnProperty(key)) {
        if (sequences[key].priorCount > 0) {
          alphs.push(key);
        }
      }
    }
    if (alphs.length == 1) alph = alphs[0];
  }
  // check if there is a sequence for the currently expected alphabet
  if (alph != null && sequences[alph] != null && sequences[alph].priorCount > 0) {
    sequence_id = sequences[alph].id;
    me = this;
    url = "../db/sequences?sequence=" + sequence_id;
    request = new XMLHttpRequest();
    request.addEventListener("load", function() {
      var data, i, prior, opt;
      data = JSON.parse(request.responseText);
      // clear the list
      me._clear_select(list);
      // add the other options
      for (i = 0; i < data.priors.length; i++) {
        prior = data.priors[i];
        opt = new Option(prior.biosample + "; " + prior.assay + "; " + prior.source, prior.id);
        optgroup.appendChild(opt);
        // try to retain old selection
        if (prior.id === selected_prior_id) opt.selected = true;
      }
      list.disabled = false;
    }, false);
    request.open("GET", url, true);
    request.send();
  } else {
    optgroup.appendChild(new Option("Awaiting motifs with matching alphabet...", ""));
  }
}
