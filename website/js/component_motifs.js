/******************************************************************************
 * Takes the surrounding div around the motif input parts and sets up the
 * events that are required to make it interactive.
 ******************************************************************************/
var MotifInput = function(container, options) {
  "use strict";
  var me, i;
  me = this;
  // parameters
  this.container = container;
  this.options = options;
  // other settings
  this.alph_type = AlphType[options.alph_type];
  if (this.alph_type == null) throw new Error("Unknown alphabet type");
  // the alphabet used to display the text motifs
  this.alphabet = this.alph_type.get_standard_alphabets()[0];
  // information on the text motifs
  this.text_motif_bounds = [];
  this.text_motifs = [];
  this.last_update = Date.now() - 500;
  this.update_timer = null;
  // keeps track of the currently drawn motif logo in the popup
  this.displayed_motif = null;
  // keeps track of where the mouse was last
  this.last_x = 0;
  this.last_y = 0;
  // file parsing related stuff
  this.file_parser = null;
  this.file_motifs = [];
  this.file_meta = null;
  this.file_error = false;
  // embed parsing related stuff
  this.embed_motifs = [];
  this.embed_meta = null;
  this.embed_error = false;
  // db related 
  this.db_alphabet = null;
  this.db_motif_count = 0;
  // get the source selector
  this.source = this.container.querySelector("select.motif_source");
  // get the xalph option if available
  // this allows the motif alphabet to extend the set alphabet
  this.xalph_input = this.container.querySelector("div.motif_xalph input");
  // get the text related parts
  // this selects the alphabet of the text motif editor
  this.text_alphabet = this.container.querySelector("select.motif_text_alphabet");
  // this surrounds the text editor components
  this.text_surround = this.container.querySelector("div.motif_text");
  // hidden input field where the motifs are written
  this.text_hin = this.text_surround.querySelector("input");
  // the div that constrains the editor
  this.text_div = this.text_surround.querySelector("div");
  // the text input
  this.text_area = this.text_surround.querySelector("textarea");
  // the span where the text is replicated with additional formatting
  this.text_backdrop = this.text_surround.querySelector("span");
  this.text_backdrop.innerHTML = "";
  // get the custom alphabet related parts
  this.cust_surround = this.container.querySelector("span.motif_text_alphabet_file");
  this.cust_input = this.cust_surround.querySelector("input");
  this.cust_popup = this.cust_surround.querySelector("div.popup");
  this.cust_error = false;
  this.cust_alphabet = null;
  // get the file related parts
  this.file_surround = this.container.querySelector("span.motif_file");
  this.file_indicator = this.file_surround.querySelector("span.indicator");
  this.file_input = this.file_surround.querySelector("input");
  this.file_popup = this.file_surround.querySelector("div.popup");
  // get the embed related parts
  this.embed_surround = this.container.querySelector("span.motif_embed");
  this.embed_hin = (this.embed_surround != null ? this.embed_surround.querySelector("input.data") : null);
  // get the area used to display the alphabet
  this.alph_info = this.container.querySelector("span.motif_alphabet_info");
  // get the database related parts
  this.db_listing = this.container.querySelector("select.listing");
  // disable alphabets that are not allowed
  for (i = 0; i < this.text_alphabet.options.length; i++) {
    var alphabet_name = this.text_alphabet.options[i].value;
    if (alphabet_name == "custom") {
      this.text_alphabet.options[i].disabled = !this.alph_type.allow_custom_alphabets();
    } else {
      var alphabet = AlphStd[alphabet_name];
      this.text_alphabet.options[i].disabled = !this.alph_type.matches(alphabet);
    }
  }
  // check that the seleted alphabet is not disabled
  if (this.text_alphabet.options[this.text_alphabet.selectedIndex].disabled) {
    for (i = 0; i < this.text_alphabet.options.length; i++) {
      if (!this.text_alphabet.options[i].disabled) {
        this.text_alphabet.options[i].selected = true;
        break;
      }
    }
    if (i == this.text_alphabet.options.length) {
      throw new Error("No alphabets were allowed!");
    }
  }
  // parse the embeded motifs
  if (this.embed_hin) {
    (new MotifParser({
      "error": function(is_error, message, reasons) {
        me._embed_error(is_error, message, reasons); },
      "begin": function(size) { me._embed_begin(size); },
      "end": function () { me._embed_end(); },
      "meta": function (info) { me._embed_motif_meta(info); },
      "motif": function (motif) { me._embed_motif(motif); }
    })).process_blob(new Blob([this.embed_hin.value], {"type": "text/plain"}));
  }

  // create the popup
  this.popup = document.createElement("canvas");
  this.popup.className = "pop_logo";
  document.getElementsByTagName("body")[0].appendChild(this.popup);
  // activate motif editor
  this._source_update();
  this._file_update();
  this._cust_alphabet_update();
  this._text_alphabet_update();
  // add event listeners
  // detect changes in the motif source selection
  this.source.addEventListener("change", function() { me._source_update(); }, false);
  // detect changing of the uploaded file
  this.file_input.addEventListener("change", function(e) { me._file_update(); }, false);
  // detect changes in the alphabet selection for the text motif input
  this.text_alphabet.addEventListener("change", function() { me._text_alphabet_update(); }, false);
  // detect changes in the custom alphabet
  this.cust_input.addEventListener("change", function(e) { me._cust_alphabet_update(); }, false);
  // detect typing or pasting into motif text area
  this.text_area.addEventListener('input', function() { me._text_update(); }, false);
  // detect which text motif the mouse is over
  this.text_div.addEventListener("mousemove", function(e) { me._text_mousemove(e); }, false);
  // show the logo popup when over the text box
  this.text_div.addEventListener("mouseover", function(e) { me._text_mouseover(e); }, false);
  // hide the logo popup when mouse is not over the text box
  this.text_div.addEventListener("mouseout", function (e) { me._text_mouseout(e); }, false);
  // detect changes in the selected database
  if (this.db_listing != null) {
    this.db_listing.addEventListener("change", function (e) {
      me._db_update_listing(parseInt(me.db_listing.value, 10));
    }, false);
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
 * Check for constaint violations.
 * true = Ok, false = invalid.
 ******************************************************************************/
MotifInput.prototype.check = function(restrict_alphabets, allow_like) {
  "use strict";
  var source, i, j, tmotif, exists, warnings, errors, alphabet, xalph;
  source = this.source.value;
  warnings = false;
  errors = false;
  alphabet = null;
  // The xalph option allows the alphabet to be expanded to match the expected
  // alphabet as long as it contains all the core symbols.
  xalph = (this.xalph_input != null ? this.xalph_input.checked : false);
  if (source == "text") {
    // check we have got something other than just whitespace
    exists = /\S/.test(this.text_area.value);
    // scan the motifs for problems
    for (i = 0; i < this.text_motifs.length; i++) {
      tmotif = this.text_motifs[i];
      for (j = 0; j < tmotif.tokens.length; j++) {
        if ((tmotif.tokens[j].type & SimpleMotifTokEn.ERROR) != 0) {
          errors = true;
        } else if ((tmotif.tokens[j].type & SimpleMotifTokEn.WARN) != 0) {
          warnings = true;
        }
      }
    }
    // determine the alphabet
    if (this.text_alphabet.value == "custom") {
      alphabet = this.cust_alphabet;
    } else {
      alphabet = AlphStd[this.text_alphabet.value];
    }
  } else if (source == "file") {
    exists = this.file_input.value.length > 0;
    errors = this.file_error;
    if (this.file_meta != null) alphabet = this.file_meta.alphabet;
  } else if (source == "embed") {
    exists = true;
    errors = this.embed_error;
    if (this.embed_meta != null) alphabet = this.embed_meta.alphabet;
  } else {// db
    exists = true;
    alphabet = this.db_alphabet;
  }
  if (!exists) {
    alert("Please input " + this.options.field + ".");
    return false;
  }
  if (errors) {
    alert("Please correct errors in the " + this.options.field + ".");
    return false;
  }
  if (alphabet == null) {
    alert("Please input " + this.options.field + " alphabet.");
    return false;
  }
  if (warnings) {
    if (!confirm("There are warnings for the " + this.options.field + ". Continue anyway?")) {
      return false;
    }
  }
  if (restrict_alphabets != null && restrict_alphabets.length > 0 && alphabet != null) {
    if (!AlphabetUtil.any_equal([alphabet], restrict_alphabets)) {
      var test_alphabets = restrict_alphabets;
      // if we allow like alphabets then select them and try the expansion test
      if (!xalph && allow_like && alphabet.get_like() != null) {
        test_alphabets = [];
        for (i = 0; i < restrict_alphabets.length; i++) {
          if (alphabet.get_like() === restrict_alphabets[i].get_like()) {
            test_alphabets.push(restrict_alphabets[i]);
          }
        }
        xalph = true;
      }
      if (xalph) {
        var expandable = [];
        var expandable_with_warnings = [];
        var nonexpandable = [];
        for (i = 0; i < test_alphabets.length; i++) {
          switch (alphabet.check_core_subset(test_alphabets[i])) {
            case 0:
              nonexpandable.push(test_alphabets[i]);
              break;
            case -1:
              expandable_with_warnings.push(test_alphabets[i]);
              break;
            case 1:
              expandable.push(test_alphabets[i]);
              break;
          }
        }
        if (expandable.length == 0) {
          if (expandable_with_warnings.length > 0) {
            if (!confirm("Warning: the alphabet expansion for the " +
                  this.options.field + " converts the " + alphabet + " to the " +
                  AlphabetUtil.display(expandable_with_warnings) + 
                  " and requires changing the complementation " +
                  "rules. Continue anyway?")) {
              return false;
            }
          } else {
	    alert("The " + this.options.field + " use(s) the " + 
		alphabet + " alphabet which is not an expandable subset " +
		"of the " + AlphabetUtil.display(restrict_alphabets) +
		" used by the other inputs.");
	    return false;
          }
        }
      } else {
	var test_alphabets = restrict_alphabets;
	if (test_alphabets.length==1 && 
	    (
	      (alphabet=="DNA" && test_alphabets[0]=="RNA") 
	     || (alphabet=="RNA" && test_alphabets[0]=="DNA") 
	    )
	  ) { 
	  // DNA and RNA are compatible alphabets.
	} else {
	  alert("The " + this.options.field + " are in the " + 
	      alphabet + " alphabet but to be compatible " +
	      "with other inputs they should be in the " + 
	      AlphabetUtil.display(restrict_alphabets) + "." +
	      (this.xalph_input != null ?
	       "  If you think the two alphabets are compatible " +
	       "you can check the 'Allow alphabet expansion' option." : ""));
           return false;
        }
      }
    }
  }
  return true;
};

/******************************************************************************
 * Check if changed from default state.
 * true = changed, false = default.
 ******************************************************************************/
MotifInput.prototype.changed = function() {
  "use strict";
  var source;
  if (!this.source.options[this.source.selectedIndex].defaultSelected) return true;
  source = this.source.value;
  if (source == "text") {
    if (this.text_area.value.length != 0) return true;
  } else if (source == "file") {
    if (this.file_input.value.length != 0) return true;
  }
  return false;
};

/******************************************************************************
 * Reset to the default state.
 * Note that calls to changed() after a call to reset() should return false.
 ******************************************************************************/
MotifInput.prototype.reset = function() {
  "use strict";
  var i, opt;
  // stop actions in progress
  if (this.file_parser) {
    this.file_parser.cancel();
    this.file_parser = null;
  }
  if (this.update_timer != null) {
    clearTimeout(this.update_timer);
    this.update_timer = null;
  }
  // reset the text input
  this.text_area.value = "";
  for (i = 0; i < this.text_alphabet.options.length; i++) {
    opt = this.text_alphabet.options[i];
    opt.selected = opt.defaultSelected;
  }
  if (this.text_alphabet.options[this.text_alphabet.selectedIndex].disabled) {
    for (i = 0; i < this.text_alphabet.options.length; i++) {
      if (!this.text_alphabet.options[i].disabled) {
        this.text_alphabet.options[i].selected = true;
        break;
      }
    }
  }
  this.last_update = Date.now() - 500;
  this._text_alphabet_update();
  // reset the file input
  this.file_input.value = "";
  this._file_update();
  // reset the custom alphabet file input
  this.cust_input.value = "";
  this._cust_alphabet_update();
  // reset the selected source option
  for (i = 0; i < this.source.options.length; i++) {
    opt = this.source.options[i];
    opt.selected = opt.defaultSelected;
  }
  this._source_update();
};

/******************************************************************************
 * Get the alphabet used for the motifs.
 ******************************************************************************/
MotifInput.prototype.get_alphabet = function () {
  "use strict";
  var source;
  // determine the source
  source = this.source.value;
  if (source == "text") {
    if (this.text_alphabet.value == "custom") {
      return this.cust_alphabet;
    } else {
      return AlphStd[this.text_alphabet.value];
    }
  } else if (source == "file") {
    if (this.file_meta) return this.file_meta["alphabet"];
  } else if (source == "embed") {
    if (this.embed_meta) return this.embed_meta["alphabet"];
  } else { // db
    return this.db_alphabet;
  }
  return null;
};

/******************************************************************************
 * Get the alphabet used for the motifs. Returns in array form.
 ******************************************************************************/
MotifInput.prototype.get_alphabets = function() {
  var alph = this.get_alphabet();
  if (alph != null) return [alph];
  return null;
};

/******************************************************************************
 * Get the number of motifs.
 ******************************************************************************/
MotifInput.prototype.get_motif_count = function () {
  "use strict";
  var source;
  source = this.source.value;
  if (source == "text") {
    return this.text_motifs.length;
  } else if (source == "file") {
    return this.file_motifs.length;
  } else if (source == "embed") {
    return this.embed_motifs.length;
  } else { // db
    return this.db_motif_count;
  }
  return 0;
};

/******************************************************************************
 * Get the number of motifs.
 ******************************************************************************/
MotifInput.prototype._fire_motifs_loaded = function() {
  var me = this;
  try {
    this.container.dispatchEvent(new CustomEvent("motifs_loaded", {detail: {controler: me}}));
  } catch (e) {
    if (e.message && e.name && window.console) {
      console.log("Suppressed exception " + e.name + ": " + e.message);
    }
  }
};

/******************************************************************************
 * EVENT HANDLER
 * Fired by mouse movement over the text motif input field.
 * Repositions the logo popup and calls _update_popup(y) with the cooridinates.
 ******************************************************************************/
MotifInput.prototype._text_mousemove = function(e) {
  "use strict";
  var now_time, me;
  // position the popup offset 5px right to the mouse pointer
  this.popup.style.left = (e.pageX + 5) + "px";
  this.popup.style.top = (e.pageY) + "px";
  // display the motif which the mouse is over
  this._update_popup(e.clientY);
};

/******************************************************************************
 * EVENT HANDLER
 * Fired by mouse over the text motif input field.
 * Displays the logo popup (the content is not updated).
 ******************************************************************************/
MotifInput.prototype._text_mouseover = function(e) {
  toggle_class(this.popup, "mouseover", true);
};

/******************************************************************************
 * EVENT HANDLER
 * Fired by mouse leaving the text motif input field.
 * Hides the logo popup.
 ******************************************************************************/
MotifInput.prototype._text_mouseout = function (e) {
  toggle_class(this.popup, "mouseover", false);
};

/******************************************************************************
 * EVENT HANDLER
 * Fired when the motif file is changed.
 * Does the following:
 * Any previous motif file parse operation is canceled.
 * if no file is provided then
 *  -> all status indicators are cleared.
 * else if a file is provided then
 *  -> a motif parser is created and the file is passed to it.
 ******************************************************************************/
MotifInput.prototype._file_update = function () {
  var file, me;
  // cancel previous parsing operation
  if (this.file_parser) {
    this.file_parser.cancel();
    this.file_parser = null;
  }
  // check for a file to load
  if (!(file = this.file_input.files[0])) {
    // no file to load! Clear the status. 
    substitute_classes(this.file_surround, ["good", "warning", "error"], []);
    this.file_indicator.style.width = "0%";
    this.file_popup.innerHTML = "";
    this._show_alph_name();
    return;
  }
  me = this;
  this.file_parser = new MotifParser({
    "error": function(is_error, message, reasons) {
      me._file_error(is_error, message, reasons); },
    "begin": function(size) { me._file_begin(size); },
    "end": function (error_list) { me._file_end(error_list); },
    "progress": function (fraction) { me._file_progress(fraction); },
    "meta": function (info) { me._file_motif_meta(info); },
    "motif": function (motif) { me._file_motif(motif); }
  });
  this.file_parser.process_blob(file);
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called when the parser of the embeded motifs discovers errors.
 * Updates the state used by check() to decide if the field is valid.
 ******************************************************************************/
MotifInput.prototype._embed_error = function (is_error, message, reasons) {
  this.embed_error |= is_error;
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called when parsing the embeded motif begins.
 * Currently does nothing.
 ******************************************************************************/
MotifInput.prototype._embed_begin = function (size) {
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called when parsing the embeded motif ends.
 * When the source is set to the embedded motifs it fires the motifs_loaded event.
 ******************************************************************************/
MotifInput.prototype._embed_end = function () {
  if (this.source.value == "embed") this._fire_motifs_loaded();
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called when the embedded motif alphabet and background is known.
 * Stores the embeded motif meta information.
 ******************************************************************************/
MotifInput.prototype._embed_motif_meta = function (info) {
  "use strict";
  this.embed_meta = info;
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called for each motif found in the embeded motifs data.
 * Stores the embeded motif (though currently only used to know the motif count).
 ******************************************************************************/
MotifInput.prototype._embed_motif = function (motif) {
  "use strict";
  this.embed_motifs.push(motif);
};

/******************************************************************************
 * HELPER FUNCTION
 * Displays a warning/error message to a table in the passed popup box.
 ******************************************************************************/
MotifInput.prototype._add_error_entry = function (popup, is_error, message, reasons) {
  "use strict";
  var table, row, cell, i, ul, li;
  table = popup.querySelector("table");
  if (table == null) {
    table = document.createElement("table");
    popup.innerHTML = "";
    popup.appendChild(table);
  }
  row = table.insertRow(table.rows.length);
  row.className = (is_error ? "error" : "warning");
  cell = row.insertCell(row.cells.length);
  cell.appendChild(document.createTextNode(is_error ? "\u2718" : "\u26A0"));
  cell = row.insertCell(row.cells.length);
  cell.appendChild(document.createTextNode(message));
  if (reasons != null) {
    ul = document.createElement("ul");
    for (i = 0; i < reasons.length; i++) {
      li = document.createElement("li");
      li.appendChild(document.createTextNode(reasons[i]));
      ul.appendChild(li);
    }
    cell.appendChild(ul);
  }
};

/******************************************************************************
 * HELPER FUNCTION
 * Displays an warning/error message in a popup box that displays when hovering
 * over the motif file input field.
 * Changes the coloured border around the motif file input field to yellow
 * if only warnings have been reported, or red if any errors have been reported.
 * Updates the file_error status variable.
 ******************************************************************************/
MotifInput.prototype._file_error = function (is_error, message, reasons) {
  this._add_error_entry(this.file_popup, is_error, message, reasons);
  this.file_error |= is_error;
  if (this.file_error) {
    substitute_classes(this.file_surround, ["good", "warning"], ["error"]);
  } else {
    substitute_classes(this.file_surround, ["good"], ["warning"]);
  }
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called when parsing the motif file begins.
 * Clears all the errors and warnings from previous file parses.
 * Clears the displayed alphabet name.
 * Checks the size of the file to ensure it isn't too large to be considered.
 ******************************************************************************/
MotifInput.prototype._file_begin = function (size) {
  "use strict";
  substitute_classes(this.file_surround, ["warning", "error"], ["good"]);
  this.file_indicator.style.width = "0%";
  this.file_popup.innerHTML = "";
  this.file_error = false;
  if (typeof this.options["max_size"] === "number" && size > this.options["max_size"]) {
    this._file_error(true, "The file is larger than the maximum accepted size.",
        ["" + size + " > " +  this.options["max_size"]]);
  }
  this._show_alph_name(null);
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called when parsing the motif file ends.
 * If nothing went wrong then it removes the highlighting around the motif file
 * input field.
 * If the source is set to file it will fire a motifs_loaded event.
 ******************************************************************************/
MotifInput.prototype._file_end = function () {
  "use strict";
  this.file_indicator.style.width = "100%";
  if (!(this.file_error)) substitute_classes(this.file_surround, ["good"], []);
  if (this.source.value == "file") this._fire_motifs_loaded();
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called occasionally while parsing the motif file.
 * The background of the file is changed to indicate how much of the file has
 * been processed.
 ******************************************************************************/
MotifInput.prototype._file_progress = function (fraction) {
  this.file_indicator.style.width = (fraction * 100) + "%";
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called when the motif file's alphabet and background have been loaded.
 * Stores the loaded information.
 * Checks that the alphabet is valid for the alphabet type.
 ******************************************************************************/
MotifInput.prototype._file_motif_meta = function (info) {
  this.file_meta = info;
  this._show_alph_name();
  if (info.alphabet != null && info.alphabet instanceof Alphabet) {
    if (!this.alph_type.matches(info.alphabet)) {
      this._file_error(true, "" + info.alphabet + " motifs are not accepted by this tool.");
    }
  } else {
    throw new Error("Unknown alphabet");
  }
  // TODO this is for debugging, it can be removed later
  if (window.console && console.log) {
    console.log("version: " + info.version + " alphabet: " + info.alphabet + " norc: " + info.norc + " background: " + JSON.stringify(info.background));
  }
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called when the motif parser has found a motif in the file.
 * Stores the motif.
 ******************************************************************************/
MotifInput.prototype._file_motif = function (motif) {
  this.file_motifs.push(motif);
  // TODO this is for debugging, it can be removed later
  if (window.console && console.log) {
    console.log("motif id: " + motif.id);
  }
};

/******************************************************************************
 * Displays the name of an alphabet truncated to 5em width with the full name
 * as the title text.
 ******************************************************************************/
MotifInput.prototype._show_alph_name = function() {
  "use strict";
  this.alph_info.innerHTML = "";
  var alphabet = this.get_alphabet();
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
 * EVENT HANDLER
 * Fired when the custom alphabet for the text fields is changed.
 * Updates the alphabet used to display the text motifs. 
 * Also updates CSS rules related to automatic uppercasing of letters in the
 * text motif input as case-insensitive alphabets are displayed as uppercase.
 ******************************************************************************/
MotifInput.prototype._text_alphabet_update = function() {
  "use strict";
  if (this.text_alphabet.value == "custom") {
    if (this.cust_alphabet != null) {
      this.alphabet = this.cust_alphabet;
    }
  } else {
    var alphabet = AlphStd[this.text_alphabet.value];
    if (alphabet != null) {
      this.alphabet = alphabet;
    } else {
      throw new Error("Unknown alphabet constant");
    }
  }
  
  toggle_class(this.container, "custom", this.text_alphabet.value == "custom");
  toggle_class(this.container, "caseinsensitive", this.alphabet.is_case_insensitive());
  this._show_alph_name();
  this._text_update();
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called when the alphabet parser finds something wrong with the
 * alphabet definition.
 * Adds the warning/error message to a list displayed when the file input 
 * is moused over.
 ******************************************************************************/
MotifInput.prototype._cust_error = function (is_error, message, reasons) {
  "use strict";
  this._add_error_entry(this.cust_popup, is_error, message, reasons);
  this.cust_error |= is_error;
  if (this.cust_error) {
    substitute_classes(this.cust_surround, ["good", "warning"], ["error"]);
  } else {
    substitute_classes(this.cust_surround, ["good"], ["warning"]);
  }
};

/******************************************************************************
 * CALLBACK HANDLER
 * Called when the alphabet parser successfully parses an alphabet.
 ******************************************************************************/
MotifInput.prototype._cust_data = function (alphabet_data) {
  this.cust_alphabet = new Alphabet(alphabet_data);
  substitute_classes(this.cust_surround, ["good"], []);
  this._text_alphabet_update();
};

/******************************************************************************
 * EVENT HANDLER
 * Fired when the file selected for the custom alphabet is changed.
 * If no file is present then
 *    -> resets the file status indicators
 * else if a file is present then
 *    -> sets the file status to good and attempts to parse.
 ******************************************************************************/
MotifInput.prototype._cust_alphabet_update = function() {
  var file, me;
  // check for a file to load
  if (!(file = this.cust_input.files[0])) {
    // no file to load! Clear the status. 
    substitute_classes(this.cust_surround, ["good", "warning", "error"], []);
    this.cust_popup.innerHTML = "";
    return;
  }
  // reset before parsing
  substitute_classes(this.cust_surround, ["warning", "error"], ["good"]);
  this.cust_popup.innerHTML = "";
  this.cust_error = false;
  this.cust_alphabet = null;
  // start parsing
  me = this;
  var parser = new AlphabetParser({
    "error": function(is_error, message, reasons) {
      me._cust_error(is_error, message, reasons);
    },
    "data": function(alphabet_data) {
      me._cust_data(alphabet_data);
    }
  });
  parser.process_blob(file);
};

/******************************************************************************
 * Converts the parsed motif tokens into a backdrop with highlighting on
 * portions with errors.
 ******************************************************************************/
MotifInput.prototype._make_backdrop = function(simple_motif) {
  "use strict";
  var all, box, i, east, south, indicator;
  all = document.createElement("span");
  for (i = 0; i < simple_motif.tokens.length; i++) {
    if ((simple_motif.tokens[i].type & SimpleMotifTokEn.SPACE) !== 0) {
      all.appendChild(document.createTextNode(simple_motif.tokens[i].text));
    } else {
      box = document.createElement("span");
      box.className = "token";
      box.appendChild(document.createTextNode(simple_motif.tokens[i].text));
      if ((simple_motif.tokens[i].type & SimpleMotifTokEn.ERROR) !== 0) {
        box.className += "  error";
      } else if ((simple_motif.tokens[i].type & SimpleMotifTokEn.WARN) !== 0) {
        box.className += " warn";
      }
      east = ((simple_motif.tokens[i].type & SimpleMotifTokEn.EAST) !== 0);
      south = ((simple_motif.tokens[i].type & SimpleMotifTokEn.SOUTH) !== 0);
      if (east || south) {
        indicator = document.createElement("div");
        indicator.className = "indicator";
        indicator.className += (east ? " east" : " south");
        indicator.appendChild(document.createTextNode(simple_motif.tokens[i].indicator));
        box.appendChild(indicator);
      }
      all.appendChild(box);
    }
  }
  return all;
};

/******************************************************************************
 * Processes the content of the text area as motifs.
 * Highlights any problems found in the motifs.
 * Stores the Minimal MEME text format motif representation in a hidden field
 * ready for sending to the server.
 * Triggers the updating of the displayed logo (at most once every
 * 500 milliseconds).
 * Fires the motifs_loaded event.
 ******************************************************************************/
MotifInput.prototype._text_update = function() {
  "use strict";
  var nodes, text, line_starts, part, i, last, groups, box, bracket;
  var now_time, me;
  // get the text
  text = this.text_area.value;
  // 1) find positions of line starts
  line_starts = [0];
  for (i = text.indexOf('\n'); i != -1; i = text.indexOf('\n', i + 1)) {
    line_starts.push(i + 1);
  }
  line_starts.push(text.length);
  // 2) find starts and ends of contiguous non-empty lines as each of these
  // will be processed as a motif.
  for (i = 0, last = -1, groups = []; i < (line_starts.length - 1); i++) {
    part = text.substring(line_starts[i], line_starts[i+1]);
    if (/^\s*$/.test(part)) {
      // empty line
      if (last != -1) {
        groups.push({"start": line_starts[last], "end": line_starts[i] - 1});
        last = -1;
      }
    } else if (/^\s*>/.test(part)) {
      // begining of a new motif
      if (last != -1) {
        groups.push({"start": line_starts[last], "end": line_starts[i] - 1});
      }
      last = i;
    } else {
      //non-empty line
      if (last == -1) {
        last = i;
      }
    }
  }
  if (last != -1) {
    groups.push({"start": line_starts[last], "end": text.length});
  }
  // 3) Parse each motif and create a backdrop
  this.text_motifs = [];
  nodes = document.createElement("span");
  for (i = 0; i < groups.length; i++) {
    // add the whitespace between the motifs to the nodes list
    if (i == 0) {
      if (groups[0].start > 0) {
        nodes.appendChild(document.createTextNode(text.substring(0, groups[0].start)));
      }
    } else if (groups[i-1].end < groups[i].start)  {
      nodes.appendChild(document.createTextNode(text.substring(groups[i-1].end, groups[i].start)));
    }
    // parse the motif according to the current alphabet
    var tmotif = new SimpleMotif(text.substring(groups[i].start, groups[i].end), this.alphabet, i+1);
    // create a box to put the motif in
    box = document.createElement("div");
    box.className = "motif_box";
    // output the parsed motif tokens with highlighting to identify errors
    box.appendChild(this._make_backdrop(tmotif));
    // use CSS tricks to make it look like there are brackets surrounding it
    bracket = document.createElement("div");
    bracket.className = "motif_bracket_left";
    box.appendChild(bracket);
    bracket = document.createElement("div");
    bracket.className = "motif_bracket_right";
    box.appendChild(bracket);
    // add the motif to the nodes list
    nodes.appendChild(box);
    // store the parsed motif
    this.text_motifs[i] = tmotif;
  }
  // add the whitespace at the end to the nodes list
  if (groups.length == 0) {
    nodes.appendChild(document.createTextNode(text));
  } else if (groups[groups.length - 1].end < text.length) {
    nodes.appendChild(document.createTextNode(text.substring(groups[groups.length - 1].end)));
  }
  // replace the existing backdrop with the new one we just created
  if (this.text_backdrop.firstChild) {
    this.text_backdrop.replaceChild(nodes, this.text_backdrop.firstChild);
  } else {
    this.text_backdrop.appendChild(nodes);
  }
  // convert the parsed motifs into MEME text format.
  this.text_hin.value = this._text_as_meme();
  // determine the vertical divisions between motifs
  this._index_active_regions();
  // update the displayed logo though ensure that it only happens
  // once every 500 milliseconds.
  now_time = Date.now();
  me = this;
  if (this.update_timer != null) {
    clearTimeout(this.update_timer);
    this.update_timer = null;
  }
  if (now_time - this.last_update < 500) {
    this.update_timer = setTimeout(function() { me._update_popup(); }, 500 - (now_time - this.last_update));
  } else {
    this._update_popup(); 
    this.last_update = now_time;
  }
  // fire the motifs_loaded event
  this._fire_motifs_loaded();
};

/******************************************************************************
 * Processes information received from the server for a listing.
 * Extracts the listing alphabet and motif count.
 * Fires a motifs_loaded event.
 ******************************************************************************/
MotifInput.prototype._db_update_listing2 = function (doc) {
  "use strict";
  var motif_db = doc.querySelector("motif_db");
  var files = motif_db.querySelectorAll("file");
  var i, count;
  for (count = 0, i = 0; i < files.length; i++) {
    count += parseInt(files[i].getAttribute("count"), 10);
  }
  this.db_alphabet = AlphStd[motif_db.getAttribute("alphabet")];
  this.db_motif_count = count;
  this._show_alph_name();
  this._fire_motifs_loaded();
};

/******************************************************************************
 * Sends a request to the server for information on a listing (alphabet + count).
 ******************************************************************************/
MotifInput.prototype._db_update_listing = function (listing_id) {
  "use strict";
  var request, url, i, me;
  // so we can access in the inner function
  me = this;
  // now send the request
  url = "../db/motifs?listing=" + listing_id;
  request = new XMLHttpRequest();
  request.addEventListener("load", function(evt) { 
    me._db_update_listing2(request.responseXML); }, false);
  request.open("GET", url, true);
  request.send();
};

/******************************************************************************
 * Populate the list of motif databases for the category.
 ******************************************************************************/
MotifInput.prototype._db_update_category = function (category_id) {
  "use strict";
  var request, url, i, me, list, optgroup;
  if (this.db_listing == null) return;
  // so we can access in the inner function
  me = this;
  list = this.db_listing;
  list.disabled = true;
  optgroup = list.getElementsByTagName("optgroup")[0];
  // clear previous values
  while (optgroup.childNodes.length > 0) {
    optgroup.removeChild(optgroup.lastChild);
  }
  // now send the request
  url = "../db/motifs?category=" + category_id;
  request = new XMLHttpRequest();
  request.addEventListener("load", function(evt) { 
    var listings, all_l, listing, i, id, name, alphabets;
    listings = request.responseXML.firstChild;
    // add the other options
    all_l = listings.getElementsByTagName("l");
    for (i = 0; i < all_l.length; i++) {
      listing = all_l[i];
      id = listing.getAttribute("i");
      name = listing.getAttribute("n");
      //alphabets = listing.getAttribute("a");
      optgroup.appendChild(new Option(name, id));
    }
    // re-enable the list
    list.disabled = false;
    if (list.value) {
      me._db_update_listing(parseInt(list.value, 10));
    }
  }, false);
  request.open("GET", url, true);
  request.send();
};

/******************************************************************************
 * EVENT HANDLER
 * Called when the source is changed.
 * Disables the unused field so they are not sent to the server when the form is
 * submitted.
 * Changes the class on the container so css can hide/show required parts.
 * If a database is selected calls _db_update_category(source) to get the
 * required information.
 ******************************************************************************/
MotifInput.prototype._source_update = function() {
  "use strict";
  var source, match, classes;
  classes = this.container.className;
  source = this.source.value;
  // disable things to avoid sending data that is not needed
  this.text_hin.disabled = true;
  this.file_input.disabled = true;
  if (this.embed_hin != null) this.embed_hin.disabled = true;
  if (source == "text") {
    this.text_hin.disabled = false;
  } else if (source == "file") {
    this.file_input.disabled = false;
  } else if (source == "embed") {
    if (this.embed_hin != null) this.embed_hin.disabled = false;
  }
  // hide/show things
  if (/^(text|file|embed)$/.test(source)) {
    this.container.className = classes.replace(/(text|file|embed|db)/, this.source.value);
  } else { // db
    this.container.className = classes.replace(/(text|file|embed|db)/, "db");
    // do a remote request for more details on this database
    this._db_update_category(source);
  }
  this._show_alph_name();
};

/******************************************************************************
 * Calculate the vertical positions that are exactly between motifs and
 * serve as the cutoffs for displaying each motif when moving the mouse over
 * the text input box.
 ******************************************************************************/
MotifInput.prototype._index_active_regions = function() {
  var error_boxes, motif_boxes, i, node;
  var y, prev_end_y, limits;
  motif_boxes = this.text_backdrop.querySelectorAll(".motif_box");

  this.text_motif_bounds = [];
  for (i = 0; i < motif_boxes.length; i++) {
    node = motif_boxes[i];
    // first find the y position of the motif relative to the surrounding div
    y = 0;
    while (node != null && node != this.text_div) {
      y += node.offsetTop;
      node = node.offsetParent;
    }
    if (i > 0) {
      // calculate the position between 2 motif boxes (used to determine which
      // should be displayed)
      this.text_motif_bounds[i-1] = ((y - prev_end_y) / 2) + prev_end_y;
    }
    // calculate the end of the motif box
    prev_end_y = y + motif_boxes[i].offsetHeight;
  }
}

/******************************************************************************
 * Based on the position of the y coordinate (or a cached value)
 * determine which motif is closest to being under the mouse.
 * Paint that motif's logo onto the popup.
 ******************************************************************************/
MotifInput.prototype._update_popup = function(y) {
  "use strict";
  var rel_box, rel_y, motif;
  var imin, imax, imid;
  // allow using previous arguments
  if (arguments.length >= 1) {
    this.last_y = y;
  } else {
    y = this.last_y;
  }
  // calculate relative to the editor box
  rel_box = this.text_div.getBoundingClientRect();
  rel_y = y - rel_box.top;
  // check to see if we have any choice
  if (this.text_motifs.length == 0) {
    motif = null;
  } else {
    // binary search for closest motif
    imin = 0;
    imax = this.text_motifs.length - 1;
    while (imin < imax) {
      imid = imin + Math.floor((imax - imin) / 2);
      // note imid < imax
      if (this.text_motif_bounds[imid] < rel_y) {
        imin = imid + 1;
      } else {
        imax = imid;
      }
    }
    motif = this.text_motifs[imin];
  }
  // draw closest motif
  if (motif != null) {
    if (motif != this.displayed_motif) {
      if (Date.now() - this.last_update > 400) {
        draw_logo_on_canvas(logo_1(this.alphabet, "", 
              MotifUtils.as_pspm(motif)), this.popup, false, 0.5);
      }
      this.displayed_motif = motif;
    }
  }
  toggle_class(this.popup, "logo_exists", motif != null);
};

/******************************************************************************
 * Converts the motifs parsed from the text input into MEME text motif format.
 * Returns the MEME text format motifs.
 ******************************************************************************/
MotifInput.prototype._text_as_meme = function() {
  var i, out, motif, pspm;
  out = "";
  for (i = 0; i < this.text_motifs.length; i++) {
    motif = this.text_motifs[i];
    pspm = MotifUtils.as_pspm(motif);
    out += pspm.as_meme({"alphabet": this.alphabet, "with_header": (i === 0)});
  }
  return out;
};

