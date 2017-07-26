var current_motif = 0;
var meme_alphabet = new Alphabet(data.alphabet, data.background.freqs);

var DelayLogoTask = function(logo, canvas) {
  this.logo = logo;
  this.canvas = canvas;
};

DelayLogoTask.prototype.run = function () {
  draw_logo_on_canvas(this.logo, this.canvas, false);
};

function motif_pspm(index) {
  var motif, pwm, psm, name, ltrim, rtrim, nsites, evalue;
  // get motif
  motif = data["motifs"][index];
  // get motif paramters
  pwm = motif["pwm"]; 
  psm = motif["psm"];
  name = "" + (index + 1); ltrim = 0; rtrim = 0; 
  nsites = motif["nsites"]; evalue = motif["evalue"];
  // make pspm
  return new Pspm(pwm, name, ltrim, rtrim, nsites, evalue, psm);
}

function motif_count_matrix(index) {
  return motif_pspm(index).as_count_matrix();
}

function motif_prob_matrix(index) {
  return motif_pspm(index).as_probability_matrix();
}

function motif_minimal_meme(index) {
  return motif_pspm(index).as_meme({
    "with_header": true, 
    "with_pspm": true,
    "with_pssm": true,
    "version": data["version"],
    "alphabet": meme_alphabet,
    "strands": (meme_alphabet.has_complement() && data.options.revcomp ? 2 : 1)
  });
}

function motif_fasta(index) {
  "use strict";
  var motif, sites, site, seq, sequences, sequence, i, num, counter, out;
  counter = {};
  sequences = data["sequence_db"]["sequences"];
  motif = data["motifs"][index];
  sites = motif["sites"];
  out = "";
  for (i = 0; i < sites.length; i++) {
    site = sites[i];
    seq = site["seq"];
    sequence = sequences[seq];
    counter[seq] = (num = counter[seq]) ? (++num) : (num = 1); // inc counter
    if (i !== 0) {out += "\n";}
    out += ">" + sequence["name"] + "_site_" + num + " offset= " + site["pos"] + 
      (site["rc"] ? " RC\n" : "\n");
    out += site["match"];
  }
  return out;
}

function motif_raw(index) {
  "use strict";
  var sites, i, out;
  sites = data["motifs"][index]["sites"];
  out = "";
  for (i = 0; i < sites.length; i++) {
    if (i !== 0) {out += "\n";}
    out += sites[i]["match"];
  }
  return out;
}

function clone_template(template) {
  "use strict";
  var node, help_btns, i, button;
  node = $(template).cloneNode(true);
  toggle_class(node, "template", false);
  node.id = "";
  help_btns = node.querySelectorAll(".help");
  for (i = 0; i < help_btns.length; i++) {
    button = help_btns[i];
    if (button.hasAttribute("data-topic")) {
      button.tabIndex = "0";
      button.addEventListener("click", __toggle_help, false);
      button.addEventListener("keydown", __toggle_help, false);
    }
  }
  return node;
}

function set_tvar(template, tvar, value) {
  var node;
  node = find_child(template, tvar);
  if (node === null) {
    throw new Error("Template does not contain variable " + tvar);
  }
  node.innerHTML = "";
  if (typeof value !== "object") {
    node.appendChild(document.createTextNode(value));
  } else {
    node.appendChild(value);
  }
}

function make_logo(alphabet, pspm, rc, offset, className) {
  if (rc) pspm = pspm.copy().reverse_complement(alphabet);
  var logo = new Logo(alphabet, "");
  logo.add_pspm(pspm, offset);
  var canvas = document.createElement('canvas');
  canvas.height = 50;
  canvas.width = 0;
  canvas.className = className;
  size_logo_on_canvas(logo, canvas, false);
  add_draw_task(canvas, new DelayLogoTask(logo, canvas));
  return canvas;
}

function make_small_logo(alphabet, pspm, options) {
  if (typeof options === "undefined") options = {};
  if (options.rc) pspm = pspm.copy().reverse_complement(alphabet);
  var logo = new Logo(alphabet, {x_axis: false, y_axis: false});
  logo.add_pspm(pspm, (typeof options.offset === "number" ? options.offset : 0));
  var canvas = document.createElement('canvas');
  if (typeof options.className === "string") canvas.className = options.className;
  if (typeof options.width === "number" && options.width > 0) {
    canvas.height = 0;
    canvas.width = options.width;
    draw_logo_on_canvas(logo, canvas, false);
  } else {
    draw_logo_on_canvas(logo, canvas, false, 1/3);
  }
  return canvas;
}

function make_large_logo(alphabet, pspm, rc, offset, className) {
  if (rc) pspm = pspm.copy().reverse_complement(alphabet);
  var logo = new Logo(alphabet, "");
  logo.add_pspm(pspm, offset);
  var canvas = document.createElement('canvas');
  canvas.height = 200;
  canvas.width = 0;
  canvas.className = className;
  size_logo_on_canvas(logo, canvas, false);
  add_draw_task(canvas, new DelayLogoTask(logo, canvas));
  return canvas;
}

function make_sym_btn(symbol, title, action) {
  var box;
  box = document.createElement("div");
  box.tabIndex = 0;
  box.className = "sym_btn";
  box.appendChild(document.createTextNode(symbol));
  box.title = title;
  box.addEventListener('click', action, false);
  box.addEventListener('keydown', action, false);
  return box;
}

function make_seq(alphabet, seq) {
  var i, j, letter, lbox, sbox;
  sbox = document.createElement("span");
  for (i = 0; i < seq.length; i = j) {
    letter = seq.charAt(i);
    for (j = i+1; j < seq.length; j++) {
      if (seq.charAt(j) !== letter) {
        break;
      }
    }
    lbox = document.createElement("span");
    lbox.style.color = alphabet.get_colour(alphabet.get_index(letter));
    lbox.appendChild(document.createTextNode(seq.substring(i, j)));
    sbox.appendChild(lbox);
  }
  return sbox;
}

//
// make_pv_text
//
// Returns the string p-value, with the p italicised.
///
function make_pv_text() {
  var pv_text = document.createElement("span");
  var pv_italic_text = document.createElement("span");
  pv_italic_text.appendChild(document.createTextNode("p"));
  pv_italic_text.style.fontStyle = "italic";
  pv_text.appendChild(pv_italic_text);
  pv_text.appendChild(document.createTextNode("-value"));
  return pv_text;
}

function append_site_entries(tbody, motif, site_index, count) {
  "use strict";
  var i, end;
  var sites, site, sequences, sequence;
  var rbody;
  if (typeof count !== "number") {
    count = 20;
  }
  sequences = data["sequence_db"]["sequences"];
  sites = motif["sites"];
  end = Math.min(site_index + count, sites.length);
  for (i = site_index; i < end; i++) {
    site = sites[i];
    sequence = sequences[site["seq"]];

    rbody = tbody.insertRow(tbody.rows.length);
    add_text_cell(rbody, "" + (site["seq"] + 1) + ".", "site_num");
    add_text_cell(rbody, sequence["name"], "site_name");
    add_text_cell(rbody, site["rc"] ? "-" : "+", "site_strand");
    add_text_cell(rbody, site["pos"] + 1, "site_start");
    add_text_cell(rbody, site["pvalue"].toExponential(2), "site_pvalue");
    add_text_cell(rbody, site["lflank"], "site lflank");
    add_cell(rbody, make_seq(meme_alphabet, site["match"]), "site match");
    add_text_cell(rbody, site["rflank"], "site rflank");
  }
  return i;
}

function make_site_entries() {
  "use strict";
  var region;
  region = this;
  if (region.data_site_index >= region.data_motif["sites"].length) {
    // all sites created
    region.removeEventListener('scroll', make_site_entries, false);
    return;
  }
  // if there's still 100 pixels to scroll than don't do anything yet
  if (region.scrollHeight - (region.scrollTop + region.offsetHeight) > 100) {
    return;
  }

  region.data_site_index = append_site_entries(
      find_child(region, "sites_tbl").tBodies[0], 
      region.data_motif, region.data_site_index, 20
    ); 
}

function make_sites(motif) {
  "use strict";
  function add_site_header(row, title, nopad, help_topic, tag_class) {
    var div, divcp, th;
    th = document.createElement("th");
    div = document.createElement("div");
    div.className = "sites_th_inner";
    if (typeof title !== "object") {
      title = document.createTextNode("" + title);
    }
    div.appendChild(title);
    if (help_topic) {
      div.appendChild(document.createTextNode("\xA0"));
      div.appendChild(help_button(help_topic));
    }
    divcp = div.cloneNode(true);
    divcp.className = "sites_th_hidden";
    th.appendChild(div);
    th.appendChild(divcp);
    if (nopad) {
      th.className = "nopad";
    }
    if (tag_class) {
      th.className += " " + tag_class;
    }
    row.appendChild(th);
  }
  var outer_tbl, inner_tbl, tbl, thead, tbody, rhead;

  outer_tbl = document.createElement("div");
  outer_tbl.className = "sites_outer";

  inner_tbl = document.createElement("div");
  inner_tbl.className = "sites_inner";
  outer_tbl.appendChild(inner_tbl);

  tbl = document.createElement("table");
  tbl.className = "sites_tbl";
  inner_tbl.appendChild(tbl);

  thead = document.createElement("thead");
  tbl.appendChild(thead);
  tbody = document.createElement("tbody");
  tbl.appendChild(tbody);

  rhead = thead.insertRow(thead.rows.length);
  add_site_header(rhead, "", true);
  add_site_header(rhead, "Name", false, "pop_seq_name");
  add_site_header(rhead, "Strand", false, "pop_site_strand", "site_strand_title");
  add_site_header(rhead, "Start", false, "pop_site_start");
  add_site_header(rhead, make_pv_text(), false, "pop_site_pvalue");
  add_site_header(rhead, "", false);
  add_site_header(rhead, "Sites", true, "pop_site_match");
  add_site_header(rhead, "", false);

  inner_tbl.data_motif = motif;
  inner_tbl.data_site_index = append_site_entries(tbody, motif, 0, 20);
  if (inner_tbl.data_site_index < motif["sites"].length) {
    inner_tbl.addEventListener('scroll', make_site_entries, false);
  }
  return outer_tbl;
}

function make_motif_table_entry(row, alphabet, ordinal, motif, colw) {
  "use strict";
  function ev_sig(evalue_str) {
    "use strict";
    var ev_re, match, sig, exp, num;
    ev_re = /^(.*)e(.*)$/;
    if (match = ev_re.exec(evalue_str)) {
      sig = parseFloat(match[1]);
      exp = parseInt(match[2]);
      if (exp >= 0) {
        return false;
      } else if (exp <= -3) {
        return true;
      } else {
        return sig * Math.pow(10, exp) <= 0.05;
      }
    }
    return true;
  }
  function make_preview(alphabet, motif) {
    "use strict";
    var pspm, preview, preview_rc;
    var box, btn_box, logo_box, btn_plus, btn_minus;
    if (motif["preview_logo"]) {
      preview = motif["preview_logo"];
      preview_rc = motif["preview_logo_rc"];
    } else {
      pspm = new Pspm(motif["pwm"]);
      preview = make_logo(alphabet, pspm);
      motif["preview_logo"] = preview;
      if (alphabet.has_complement()) {
        preview_rc = make_logo(alphabet, pspm, true, 0, "logo_rc");
        motif["preview_logo_rc"] = preview_rc;
      }
    }
    if (preview_rc) {
      btn_plus = document.createElement("div");
      btn_plus.appendChild(document.createTextNode("+"));
      btn_plus.className = "preview_btn plus";
      btn_plus.tabIndex = "0";
      btn_plus.addEventListener("click", action_btn_rc, false);
      btn_plus.addEventListener("keydown", action_btn_rc, false);
      btn_minus = document.createElement("div");
      btn_minus.appendChild(document.createTextNode("-"));
      btn_minus.className = "preview_btn minus";
      btn_minus.tabIndex = "0";
      btn_minus.addEventListener("click", action_btn_rc, false);
      btn_minus.addEventListener("keydown", action_btn_rc, false);
      btn_box = document.createElement("div");
      btn_box.className = "preview_btn_box";
      btn_box.appendChild(btn_plus);
      btn_box.appendChild(btn_minus);
    }
    logo_box = document.createElement("div");
    logo_box.className = "preview_logo_box";
    logo_box.appendChild(preview);
    if (preview_rc) logo_box.appendChild(preview_rc);
    box = document.createElement("div");
    box.className = "preview_box";
    if (preview_rc) box.appendChild(btn_box);
    box.appendChild(logo_box);
    if (preview_rc) {
      if (motif["rc"]) {
        btn_minus.className += " active";
        logo_box.className += " show_rc_logo";
      } else {
        btn_plus.className += " active";
      }
    }
    return box;
  }
  var pspm, preview, preview_rc, c;
  row.data_motif = motif;
  row.data_ordinal = ordinal;
  if (!ev_sig(motif["evalue"])) {
    row.style.opacity = 0.4;
  }
  add_text_cell(row, "" + ordinal + ".", "motif_ordinal");
  add_cell(row, make_preview(alphabet, motif), "motif_logo");
  add_text_cell(row, motif["evalue"], "motif_evalue");
  add_text_cell(row, motif["nsites"], "motif_nsites");
  add_text_cell(row, motif["len"], "motif_width");
  add_cell(row, make_sym_btn("\u21A7", "Show more information.", 
        action_show_more), "motif_more");
  add_cell(row, 
      make_sym_btn("\u21E2", 
        "Submit the motif to another MEME Suite program or download it.",
        action_show_outpop), 
      "motif_submit");
  if (colw) {
    for (c = 0; c < row.cells.length; c++) {
      row.cells[c].style.minWidth = colw[c] + "px";
    }
  }
}

function make_motifs_table(alphabet, start_ordinal, motifs, colw, stop_reason) {
  var i, j;
  var tbl, thead, tbody, tfoot, row, preview;
  var motif, pspm;

  tbl = document.createElement("table");
  
  thead = document.createElement("thead");
  tbl.appendChild(thead);
  tbody = document.createElement("tbody");
  tbl.appendChild(tbody);
  tfoot = document.createElement("tfoot");
  tbl.appendChild(tfoot);

  row = thead.insertRow(thead.rows.length);
  add_text_header_cell(row, "", "", "motif_ordinal");
  add_text_header_cell(row, "Logo", "", "motif_logo");
  add_text_header_cell(row, "E-value", "pop_ev", "motif_evalue");
  add_text_header_cell(row, "Sites", "pop_sites", "motif_nsites");
  add_text_header_cell(row, "Width", "pop_width", "motif_width");
  add_text_header_cell(row, "More", "pop_more", "motif_more");
  add_text_header_cell(row, "Submit/Download", "pop_submit_dl", "motif_submit");

  for (i = 0; i < motifs.length; i++) {
    row = tbody.insertRow(tbody.rows.length);
    make_motif_table_entry(row, alphabet, start_ordinal + i, motifs[i], colw);
  }

  row = tfoot.insertRow(tfoot.rows.length);
  add_text_header_cell(row, stop_reason, "", "stop_reason", "", 6);

  return tbl;
}

function make_expanded_motif(alphabet, ordinal, motif, less_x, submit_x) {
  "use strict";
  var box, pspm, logo_box, large_logo, large_logo_rc, tab_logo, tab_logo_rc;
  var btn, offset, norc;

  box = clone_template("tmpl_motif_expanded");
  box.data_motif = motif;
  box.data_ordinal = ordinal;

  pspm = new Pspm(motif["pwm"]);
  if (typeof motif["rc"] !== "boolean") {
    motif["rc"] = false;
  }
  if (motif["large_logo"]) {
    large_logo = motif["large_logo"];
    large_logo_rc = motif["large_logo_rc"];
  } else {
    large_logo = make_large_logo(alphabet, pspm, false, 0);
    motif["large_logo"] = large_logo;
    if (alphabet.has_complement()) {
      large_logo_rc = make_large_logo(alphabet, pspm, true, 0, "logo_rc");
      motif["large_logo_rc"] = large_logo_rc;
    }
  }
  norc = (large_logo_rc == null);
  toggle_class(box, "norc", norc);

  logo_box = find_child(box, "tvar_logo");
  logo_box.appendChild(large_logo);
  if (large_logo_rc) logo_box.appendChild(large_logo_rc);
  toggle_class(logo_box, "show_rc_logo", motif["rc"]);

  tab_logo = find_child(box, "tvar_tab");
  tab_logo_rc = find_child(box, "tvar_tab_rc");

  toggle_class(tab_logo, "activeTab", !motif["rc"]);
  toggle_class(tab_logo_rc, "activeTab", motif["rc"]);

  tab_logo.addEventListener('click', action_rc_tab, false);
  tab_logo.addEventListener('keydown', action_rc_tab, false);
  tab_logo_rc.addEventListener('click', action_rc_tab, false);
  tab_logo_rc.addEventListener('keydown', action_rc_tab, false);

  set_tvar(box, "tvar_ordinal", ordinal); 
  set_tvar(box, "tvar_evalue", motif["evalue"]);
  set_tvar(box, "tvar_width", motif["len"]);
  set_tvar(box, "tvar_site_count", motif["nsites"]);
  set_tvar(box, "tvar_llr", motif["llr"]);
  set_tvar(box, "tvar_ic", motif["ic"]);
  set_tvar(box, "tvar_re", motif["re"]);
  set_tvar(box, "tvar_bt", motif["bt"]);
  set_tvar(box, "tvar_sites", make_sites(motif));

  offset = 32; // 1* 5px padding + 2 * 10px padding + 2 * 2px border + 3px ??

  btn = find_child(box, "tvar_less");
  btn.style.left = (less_x - offset) + "px";
  btn.addEventListener('click', action_show_less, false);
  btn.addEventListener('keydown', action_show_less, false);
  btn = find_child(box, "tvar_submit");
  btn.style.left = (submit_x - offset) + "px";
  btn.addEventListener('click', action_show_outpop, false);
  btn.addEventListener('keydown', action_show_outpop, false);
  return box;
}


//
//
///
function make_motifs() {
  "use strict";
  function pixel_value(str_in) {
    "use strict";
    var px_re, match;
    px_re = /^(\d+)px$/;
    if (match = px_re.exec(str_in)) {
      return parseInt(match[1], 10);
    }
    return 0;
  }
  var container, tbl;
  var colw, r, row, c, cell, cell_style, pad_left, pad_right;

  // make the motifs table
  container = $("motifs");
  container.innerHTML = ""; // clear content

  tbl = make_motifs_table(meme_alphabet, 1, data["motifs"], colw, data["stop_reason"]);
  container.appendChild(tbl);

  // measure table column widths
  colw = [];
  row = tbl.tBodies[0].rows[0];
  for (c = 0; c < row.cells.length; c++) {
    var padLeft, padRight;
    cell = row.cells[c];
    cell_style = window.getComputedStyle(cell, null);
    pad_left = pixel_value(cell_style.getPropertyValue("padding-left"));
    pad_right = pixel_value(cell_style.getPropertyValue("padding-right"));
    colw[c] = cell.clientWidth - pad_left - pad_right;
    if (typeof colw[c] !== "number" || colw[c] < 0) {
      colw[c] = 1;
    }
  }

  // set minimum table column widths on each row so later when we remove rows it still aligns
  for (r = 0; r < tbl.tBodies[0].rows.length; r++) {
    row = tbl.tBodies[0].rows[r];
    for (c = 0; c < row.cells.length; c++) {
      row.cells[c].style.minWidth = colw[c] + "px";
    }
  }

  // store the table column widths so we can create rows latter with the same minimums
  container.data_colw = colw;

  // calculate the x offset for the buttons
  row = tbl.tBodies[0].rows[0];
  container.data_more_x = coords(find_child(find_child(row, "motif_more"), "sym_btn"))[0];
  container.data_submit_x = coords(find_child(find_child(row, "motif_submit"), "sym_btn"))[0];

  draw_on_screen();
}

function make_meme_block(container, max_seq_len, is_scan, site) {
  "use strict";
  var motif = data.motifs[site.motif];
  var block = make_block(container, max_seq_len, site.pos, motif.len,
      site.pvalue, site.rc, site.motif, is_scan);
  var handler = (is_scan ?
      make_scan_popup(site, motif, block) :
      make_block_popup(site, motif, block));
  block.addEventListener("mouseover", handler, false);
  block.addEventListener("mouseout", handler, false);
}

function append_blocks_entries(tbody, seq_index, count) {
  "use strict";
  var i, end, j;
  var max_pvalue, max_block_height, max_seq_len, sequences;
  var sequence, sites, scans, scan;
  var container, plus, minus, rule, row;
  // define some constants
  max_seq_len = data.sequence_db.max_length;
  // determine how many to load
  end = Math.min(seq_index + count, data.sequence_db.sequences.length);
  for (i = seq_index; i < end; i++) {
    // get the sequence
    sequence = data.sequence_db.sequences[i];
    // make the containers for the block diagram
    container = make_block_container(meme_alphabet.has_complement(),
        data.options.revcomp, max_seq_len, sequence.length);
    // create blocks for the motif sites
    sites = sequence["sites"];
    for (j = 0; j < sites.length; j++)
      make_meme_block(container, max_seq_len, false, sites[j]);
    // create blocks for the scanned sites
    scan = data.scan[i];
    for (j = 0; j < scan.sites.length; j++)
      make_meme_block(container, max_seq_len, true, scan.sites[j]);
    // create a row for the sequence
    row = tbody.insertRow(tbody.rows.length);
    toggle_class(row, "empty_seq", sites.length == 0 && scan.sites.length == 0);
    toggle_class(row, "only_scan", sites.length == 0 && scan.sites.length > 0);
    add_text_cell(row, (i + 1) + ".", "blockdiag_num");
    add_text_cell(row, sequence["name"], "blockdiag_name");
    add_text_cell(row, scan["pvalue"].toExponential(2), "blockdiag_pvalue");
    add_cell(row, container, "block_td"); 
  }
  return end;
}

function make_blocks_entries() {
  "use strict";
  var region;
  region = this;
  if (region.data_blocks_index >= data["sequence_db"]["sequences"].length) {
    // all sites created
    region.removeEventListener('scroll', make_blocks_entries, false);
    return;
  }
  // if there's still 100 pixels to scroll than don't do anything yet
  if (region.scrollHeight - (region.scrollTop + region.offsetHeight) > 100) {
    return;
  }

  region.data_blocks_index = append_blocks_entries(
      find_child(region, "blocks_tbl").tBodies[0], 
      region.data_blocks_index, 20
    ); 
}

function make_blocks() {
  "use strict";
  function add_seqs_filter(container, id, checked, label_text, help_topic) {
    "use strict";
    var label, radio;
    radio = document.createElement("input");
    radio.type = "radio";
    radio.name = "seqs_display";
    radio.id = id;
    radio.checked = checked;
    radio.addEventListener('click', action_seqs_filter, false);
    label = document.createElement("label");
    label.appendChild(document.createTextNode(label_text));
    label.htmlFor = id;
    container.appendChild(radio);
    container.appendChild(label);
    if (help_topic) {
      container.appendChild(document.createTextNode("\xA0"));
      container.appendChild(help_button(help_topic));
    }
  }
  function add_blocks_header(row, title, nopad, help_topic) {
    "use strict";
    var div, divcp, th;
    th = document.createElement("th");
    div = document.createElement("div");
    div.className = "blocks_th_inner";
    if (typeof title !== "object") {
      title = document.createTextNode("" + title);
    }
    div.appendChild(title);
    if (help_topic) {
      div.appendChild(document.createTextNode("\xA0"));
      div.appendChild(help_button(help_topic));
    }
    divcp = div.cloneNode(true);
    divcp.className = "blocks_th_hidden";
    th.appendChild(div);
    th.appendChild(divcp);
    if (nopad) {
      th.className = "nopad";
    }
    row.appendChild(th);
  }
  var container;
  var page, view_height, outer_tbl, inner_tbl, tbl, thead, tbody, rhead;
  var in_view, i, seq_count;
  
  page = (document.compatMode === "CSS1Compat") ? document.documentElement : document.body;
  view_height = Math.max(page.clientHeight - 300, 300);

  container = $("blocks");
  toggle_class(container, "hide_empty_seqs", true);
  toggle_class(container, "hide_only_scan", true);
  container.innerHTML = "";
  add_seqs_filter(container, "rdo_sites_only", true, "Only Motif Sites", "pop_motif_sites");
  add_seqs_filter(container, "rdo_sites_and_scan", false, "Motif Sites+Scanned Sites", "pop_scanned_sites");
  add_seqs_filter(container, "rdo_all_seqs", false, "All Sequences", "pop_all_sequences");

  outer_tbl = document.createElement("div");
  outer_tbl.className = "blocks_outer";

  inner_tbl = document.createElement("div");
  inner_tbl.id = "blocks_scroll";
  inner_tbl.className = "blocks_inner";
  inner_tbl.style.maxHeight = view_height + "px";
  outer_tbl.appendChild(inner_tbl);

  tbl = document.createElement("table");
  tbl.className = "blocks_tbl";
  inner_tbl.appendChild(tbl);

  thead = document.createElement("thead");
  tbl.appendChild(thead);
  tbody = document.createElement("tbody");
  tbl.appendChild(tbody);

  rhead = thead.insertRow(thead.rows.length);
  add_blocks_header(rhead, "", true);
  add_blocks_header(rhead, "Name", false, "pop_seq_name");
  add_blocks_header(rhead, make_pv_text(), false, "pop_seq_pvalue");
  add_blocks_header(rhead, "Motif Location", false, "pop_motif_location");

  container.appendChild(outer_tbl);

  
  seq_count = data["sequence_db"]["sequences"].length;
  in_view = Math.max(Math.ceil(view_height / 25), 1);
  i = append_blocks_entries(tbody, 0, in_view);

  while (i < seq_count && inner_tbl.scrollHeight - (inner_tbl.scrollTop + inner_tbl.offsetHeight) < 400) {
    i = append_blocks_entries(tbody, i, 20);
  }
  inner_tbl.data_blocks_index = i;
  if (i < seq_count) {
    inner_tbl.addEventListener('scroll', make_blocks_entries, false);
  }
}

function make_scan_popup(site, motif) {
  return function (e) {
    "use strict";
    var pop, xy, padding, edge_padding, pop_left, pop_top, page_width;
    var lflank, match, rflank, pspm;
    if (!e) var e = window.event;
    pop = make_scan_popup.pop;
    if (e.type === "mouseover") {
      if (pop) return;
      pop = clone_template("tmpl_scan_info");
      pspm = new Pspm(motif.pwm);
      if (site.rc) pspm.reverse_complement(meme_alphabet);
      set_tvar(pop, "tvar_logo", make_small_logo(meme_alphabet, pspm, {"className": "scan_logo"}));
      set_tvar(pop, "tvar_motif", motif.id);
      set_tvar(pop, "tvar_pvalue", site.pvalue.toExponential(2));
      set_tvar(pop, "tvar_start", site.pos + 1);
      set_tvar(pop, "tvar_end", site.pos + motif.len);

      document.body.appendChild(pop);
      position_popup(this, pop);
      make_scan_popup.pop = pop;
    } else if (e.type === "mouseout") {
      if (pop) {
        pop.parentNode.removeChild(pop);
        make_scan_popup.pop = null;
      }
    }
  };
}

function make_block_popup(site, motif, block) {
  return function (e) {
    "use strict";
    var pop;
    var lflank, match, rflank, pspm, ruler, match_seq, match_width;
    if (!e) var e = window.event;
    pop = make_block_popup.pop;
    if (e.type === "mouseover") {
      if (pop) return;
      pop = clone_template("tmpl_block_info");
      pspm = new Pspm(motif.pwm);
      if (site.rc) { // must be dna
        pspm.reverse_complement(meme_alphabet);
        lflank = meme_alphabet.invcomp_seq(site.rflank);
        match = meme_alphabet.invcomp_seq(site.match);
        rflank = meme_alphabet.invcomp_seq(site.lflank);
      } else {
        lflank = site.lflank;
        match = site.match;
        rflank = site.rflank;
      }
      ruler = document.getElementById("measure_match");
      match_seq = make_seq(meme_alphabet, match);
      ruler.innerHTML = "";
      ruler.appendChild(match_seq);
      match_width = ruler.clientWidth;
      ruler.removeChild(match_seq);
      set_tvar(pop, "tvar_lflank", lflank);
      set_tvar(pop, "tvar_match", match_seq);
      set_tvar(pop, "tvar_rflank", rflank);
      set_tvar(pop, "tvar_logo_pad", lflank);
      set_tvar(pop, "tvar_logo", make_small_logo(meme_alphabet, pspm, {"width": match_width}));
      set_tvar(pop, "tvar_motif", motif.id);
      set_tvar(pop, "tvar_pvalue", site.pvalue.toExponential(2));
      set_tvar(pop, "tvar_start", site.pos + 1);
      set_tvar(pop, "tvar_end", site.pos + motif.len);

      document.body.appendChild(pop);
      position_popup(block, pop);
      make_block_popup.pop = pop;
    } else if (e.type === "mouseout") {
      if (pop) {
        pop.parentNode.removeChild(pop);
        make_block_popup.pop = null;
      }
    }
  };
}

function update_outpop_format(index) {
  switch(parseInt($("text_format").value)) {
    case 0: // count matrix
      $("outpop_text").value = motif_count_matrix(index);
      $("text_name").value = "motif_" + (index + 1) + "_counts.txt";
      break;
    case 1: // prob matrix
      $("outpop_text").value = motif_prob_matrix(index);
      $("text_name").value = "motif_" + (index + 1) + "_freqs.txt";
      break;
    case 2: // minimal meme
      $("outpop_text").value = motif_minimal_meme(index);
      $("text_name").value = "motif_" + (index + 1) + ".txt";
      break;
    case 3: // fasta
      $("outpop_text").value = motif_fasta(index);
      $("text_name").value = "motif_" + (index + 1) + "_fasta.txt";
      break;
    case 4: // raw
      $("outpop_text").value = motif_raw(index);
      $("text_name").value = "motif_" + (index + 1) + "_raw.txt";
      break;
    default:
      throw new Error("Unknown motif format");
  }
}

function update_outpop_motif(index) {
  "use strict";
  var motifs, motif, pspm, logo, canvas, num;
  motifs = data["motifs"];
  if (index < 0 || index >= motifs.length) {return;}
  current_motif = index;
  motif = motifs[index];
  pspm = new Pspm(motif["pwm"]);
  logo = new Logo(meme_alphabet, "");
  logo.add_pspm(pspm, 0);
  canvas = $("outpop_logo");
  canvas.width = canvas.width; // clear canvas
  draw_logo_on_canvas(logo, canvas, false);
  if (meme_alphabet.has_complement()) {
    pspm.reverse_complement(meme_alphabet);
    logo = new Logo(meme_alphabet, "");
    canvas = $("outpop_logo_rc");
    canvas.width = canvas.width; // clear canvas
    draw_logo_on_canvas(logo, canvas, false);
  }
  num = $("outpop_num");
  num.innerHTML = "";
  num.appendChild(document.createTextNode("" + (index + 1)));
  update_outpop_format(index);
}

//
// action_show_more
//
// Show more information on the motif.
///
function action_show_more(e) {
  var node, tr, tbody, table, container, motif, ordinal;
  var expanded_motif;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  // find the row that contains the cell
  node = this;
  do {
    if (node.tagName === "TR") break;
  } while (node = node.parentNode);
  if (!node) throw new Error("Expected to find row!?");
  tr = node;
  // get info
  motif = tr.data_motif;
  ordinal = tr.data_ordinal;
  // find tbody
  do {
    if (node.tagName === "TBODY") break;
  } while (node = node.parentNode);
  if (!node) throw new Error("Expected to find tbody!?");
  tbody = node;
  // find table
  do {
    if (node.tagName === "TABLE") break;
  } while (node = node.parentNode);
  if (!node) throw new Error("Expected to find table!?");
  table = node;
  // find container
  container = node.parentNode;
  // make a expanded motif
  motif["expanded"] = true;
  expanded_motif = make_expanded_motif(meme_alphabet, ordinal, motif, 
      container.data_more_x, container.data_submit_x);
  // now determine how to place it
  if (tbody.rows.length === 1) {
    // only us in the table so the table can be replaced
    container.replaceChild(expanded_motif, table);
  } else if (tbody.rows[0] === tr) {
    // first row, so remove and insert an expanded motif before
    table.deleteRow(tr.rowIndex);
    container.insertBefore(expanded_motif, table);
  } else if (tbody.rows[tbody.rows.length -1] === tr) {
    // last row, so remove and insert an expanded motif after
    table.deleteRow(tr.rowIndex);
    container.insertBefore(expanded_motif, table.nextSibling);
  } else {
    var table2, tbody2;
    table2 = table.cloneNode(false);
    table2.appendChild(table.tHead.cloneNode(true));
    tbody2 = table.tBodies[0].cloneNode(false);
    table2.appendChild(tbody2);
    container.insertBefore(table2, table.nextSibling);
    for (i = tbody.rows.length - 1; i >= 0; i--) {
      row = tbody.rows[i];
      row.parentNode.removeChild(row);
      if (row === tr) {
        break;
      }
      tbody2.insertBefore(row, tbody2.rows[0]);
    }
    container.insertBefore(expanded_motif, table2);
  }
  find_child(expanded_motif, "tvar_less").focus();
}

//
// action_show_less
//
// Show less information on the motif.
///
function action_show_less(e) {
  var btn;
  var expanded_motif, container, motif, ordinal, colw, focus_target;
  var table, tbody, tbody2, row, table_before, table_after;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  btn = this;
  // find expanded motif
  expanded_motif = find_parent(btn, "expanded_motif");
  if (!expanded_motif) throw new Error("Expected expanded motif.");
  // find the container
  container = expanded_motif.parentNode;
  // get data
  motif = expanded_motif.data_motif;
  ordinal = expanded_motif.data_ordinal;
  colw = container.data_colw;
  // get the table before
  table_before = expanded_motif.previousSibling;
  if (table_before && table_before.tagName !== "TABLE") {
    table_before = null;
  }
  // get the table after
  table_after = expanded_motif.nextSibling;
  if (table_after && table_after.tagName !== "TABLE") {
    table_after = null;
  }
  // see if there is a table below or above that we can put this in.
  // if there is a table both below and above then add this motif and
  // all ones below to the above table
  motif["expanded"] = false;
  if (table_before && table_after) {
    tbody = table_before.tBodies[0];
    row = tbody.insertRow(tbody.rows.length);
    make_motif_table_entry(row, meme_alphabet, ordinal, motif, colw);
    focus_target = find_child(row.cells[5], "sym_btn");
    container.removeChild(expanded_motif);
    tbody2 = table_after.tBodies[0];
    while (tbody2.rows.length > 0) {
      row = tbody2.rows[0];
      row.parentNode.removeChild(row);
      tbody.appendChild(row);
    }
    container.removeChild(table_after);
  } else if (table_before) {
    tbody = table_before.tBodies[0];
    row = tbody.insertRow(tbody.rows.length);
    make_motif_table_entry(row, meme_alphabet, ordinal, motif, colw);
    focus_target = find_child(row.cells[5], "sym_btn");
    container.removeChild(expanded_motif);
  } else if (table_after) {
    tbody = table_after.tBodies[0];
    row = tbody.insertRow(0);
    make_motif_table_entry(row, meme_alphabet, ordinal, motif, colw);
    focus_target = find_child(row.cells[5], "sym_btn");
    container.removeChild(expanded_motif);
  } else {
    //no table above or below!
    // make a new table
    table = make_motifs_table(meme_alphabet, ordinal, [motif], colw, data["stop_reason"]);
    focus_target = find_child(table.tBodies[0].rows[0].cells[5], "sym_btn");
    container.replaceChild(table, expanded_motif);
  }
  focus_target.focus();
}

function action_show_outpop(e) {
  "use strict";
  function init() {
    "use strict";
    var close_btn, next_btn, prev_btn, cancel_btn, do_btn;
    var tab1, tab2, tab3;
    var pnl1, pnl2, pnl3;
    var format_list;
    var tbl_submit, inputs, i, default_prog;
    close_btn = $("outpop_close");
    close_btn.addEventListener("click", action_hide_outpop, false);
    close_btn.addEventListener("keydown", action_hide_outpop, false);
    next_btn = $("outpop_next");
    next_btn.addEventListener("click", action_outpop_next, false);
    next_btn.addEventListener("keydown", action_outpop_next, false);
    prev_btn = $("outpop_prev");
    prev_btn.addEventListener("click", action_outpop_prev, false);
    prev_btn.addEventListener("keydown", action_outpop_prev, false);
    cancel_btn = $("outpop_cancel");
    cancel_btn.addEventListener("click", action_hide_outpop, false);
    do_btn = $("outpop_do");
    do_btn.addEventListener("click", action_outpop_submit, false);
    tab1 = $("outpop_tab_1");
    tab1.tabIndex = 0;
    tab1.addEventListener("click", action_outpop_tab, false);
    tab1.addEventListener("keydown", action_outpop_tab, false);
    tab2 = $("outpop_tab_2");
    tab2.tabIndex = 0;
    tab2.addEventListener("click", action_outpop_tab, false);
    tab2.addEventListener("keydown", action_outpop_tab, false);
    tab3 = $("outpop_tab_3");
    tab3.tabIndex = 0;
    tab3.addEventListener("click", action_outpop_tab, false);
    tab3.addEventListener("keydown", action_outpop_tab, false);
    pnl1 = $("outpop_pnl_1");
    pnl2 = $("outpop_pnl_2");
    pnl3 = $("outpop_pnl_3");
    toggle_class(tab1, "activeTab", true);
    toggle_class(tab2, "activeTab", false);
    toggle_class(tab3, "activeTab", false);
    pnl1.style.display = "block";
    pnl2.style.display = "none";
    pnl3.style.display = "none";
    format_list = $("text_format");
    format_list.addEventListener("change", action_outpop_format, false);
    // setup program selection
    tbl_submit = $("programs");
    // when not dna, hide the inputs for programs that require dna motifs
    toggle_class(tbl_submit, "alphabet_dna", meme_alphabet.has_complement());//TODO FIXME alphabet_dna is a bad name for a field when allowing custom alphabets
    // add a click listener for the radio buttons
    inputs = tbl_submit.querySelectorAll("input[type='radio']");
    for (i = 0; i < inputs.length; i++) {
      inputs[i].addEventListener("click", action_outpop_program, false);
    }
    // ensure that a default program option is selected for DNA and Protein
    default_prog = document.getElementById(meme_alphabet.has_complement() ? "submit_tomtom" : "submit_fimo"); //TODO FIXME Tomtom might require a more strict definition of DNA
    default_prog.checked = true;
    action_outpop_program.call(default_prog);
    // disable reverse-complement when not DNA
    $("logo_rc_option").disabled = !meme_alphabet.has_complement(); 
    // set errorbars on when ssc is on
    $("logo_ssc").addEventListener("change", action_outpop_ssc, false);
  }
  var node;
  // store the focused element
  action_hide_outpop.last_active = document.activeElement;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  // hide the help popup
  help_popup();
  // on first load initilize the popup
  if (!action_show_outpop.ready) {
    init();
    action_show_outpop.ready = true;
  }
  // load the motif logo
  node = this;
  do {
    if (/\bexpanded_motif\b/.test(node.className) || node.tagName === "TR") break;
  } while (node = node.parentNode);
  if (node === null) throw new Error("Expected node!");
  update_outpop_motif(node.data_ordinal - 1);
  // display the download popup
  $("grey_out_page").style.display = "block";
  $("download").style.display = "block";
  $("outpop_close").focus();
}

function action_hide_outpop(e) {
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  $("download").style.display = "none";
  $("grey_out_page").style.display = "none";
  if (typeof action_hide_outpop.last_active !== "undefined") {
    action_hide_outpop.last_active.focus();
  }
}

function action_outpop_next(e) {
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  update_outpop_motif(current_motif + 1);
}

function action_outpop_prev(e) {
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  update_outpop_motif(current_motif - 1);
}

function action_outpop_program() {
  "use strict";
  var table, tr, rows, i;
  tr = find_parent_tag(this, "TR");
  table = find_parent_tag(tr, "TABLE");
  rows = table.querySelectorAll("tr");
  for (i = 0; i < rows.length; i++) {
    toggle_class(rows[i], "selected", rows[i] === tr);
  }
}

function action_outpop_ssc() {
  "use strict";
  $("logo_err").value = $("logo_ssc").value;
}

function action_outpop_submit(e) {
  "use strict";
  var form, input, program, motifs;
  // find out which program is selected
  var radios, i;
  radios = document.getElementsByName("program");
  program = "fimo"; // default to fimo, since it works with all alphabet types
  for (i = 0; i < radios.length; i++) {
    if (radios[i].checked) program = radios[i].value;
  }

  motifs = motif_minimal_meme(current_motif);
  form = document.createElement("form");
  form.setAttribute("method", "post");
  form.setAttribute("action", site_url + "/tools/" + program);
  
  input = document.createElement("input");
  input.setAttribute("type", "hidden");
  input.setAttribute("name", "motifs_embed");
  input.setAttribute("value", motifs);
  form.appendChild(input);

  document.body.appendChild(form);
  form.submit();
  document.body.removeChild(form);
}

function action_outpop_download_motif(e) {
  $("text_form").submit();
}

function action_outpop_download_logo(e) {
  "use strict";
  $("logo_motifs").value = motif_minimal_meme(current_motif);
  $("logo_form").submit();
}

function action_btn_rc(e) {
  "use strict";
  var node, tr, motif, box, logo_box, tab_st, tab_rc, rc;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  node = this;
  do {
    if (node.tagName === "TR") break;
  } while (node = node.parentNode);
  if (!node) throw new Error("Expected to find row!?");
  tr = node;
  // get info
  motif = tr.data_motif;
  box = find_parent(this, "preview_box");
  logo_box = find_child(box, "preview_logo_box");
  tab_st = find_child(box, "plus");
  tab_rc = find_child(box, "minus");
  rc = (this === tab_rc);
  motif["rc"] = rc;
  toggle_class(logo_box, "show_rc_logo", rc);
  toggle_class(tab_st, "active", !rc);
  toggle_class(tab_rc, "active", rc);
}

function action_rc_tab(e) {
  "use strict";
  var box, logo_box, tab_st, tab_rc, rc;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  box = find_parent(this, "expanded_motif");
  logo_box = find_child(box, "tvar_logo");
  tab_st = find_child(box, "tvar_tab");
  tab_rc = find_child(box, "tvar_tab_rc");
  rc = (this === tab_rc);
  box.data_motif["rc"] = rc;
  toggle_class(logo_box, "show_rc_logo", rc);
  toggle_class(tab_st, "activeTab", !rc);
  toggle_class(tab_rc, "activeTab", rc);
}

function action_outpop_tab(e) {
  "use strict";
  var tab1, tab2, tab3, pnl1, pnl2, pnl3, do_btn;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  tab1 = $("outpop_tab_1");
  tab2 = $("outpop_tab_2");
  tab3 = $("outpop_tab_3");
  pnl1 = $("outpop_pnl_1");
  pnl2 = $("outpop_pnl_2");
  pnl3 = $("outpop_pnl_3");
  do_btn = $("outpop_do");

  toggle_class(tab1, "activeTab", (this === tab1));
  toggle_class(tab2, "activeTab", (this === tab2));
  toggle_class(tab3, "activeTab", (this === tab3));
  pnl1.style.display = ((this === tab1) ? "block" : "none");
  pnl2.style.display = ((this === tab2) ? "block" : "none");
  pnl3.style.display = ((this === tab3) ? "block" : "none");
  do_btn.value = ((this === tab1) ? "Submit" : "Download");
  do_btn.removeEventListener("click", action_outpop_submit, false);
  do_btn.removeEventListener("click", action_outpop_download_logo, false);
  do_btn.removeEventListener("click", action_outpop_download_motif, false);
  if (this === tab1) {
    do_btn.addEventListener("click", action_outpop_submit, false);
  } else if (this === tab2) {
    do_btn.addEventListener("click", action_outpop_download_motif, false);
  } else {
    do_btn.addEventListener("click", action_outpop_download_logo, false);
  }
}

function action_seqs_filter() {
  "use strict";
  var block_container;
  block_container = $("blocks");
  if ($("rdo_all_seqs").checked) {
    toggle_class(block_container, "hide_empty_seqs", false);
    toggle_class(block_container, "hide_only_scan", false);
  } else if ($("rdo_sites_and_scan").checked) {
    toggle_class(block_container, "hide_empty_seqs", true);
    toggle_class(block_container, "hide_only_scan", false);
  } else if ($("rdo_sites_only").checked) {
    toggle_class(block_container, "hide_empty_seqs", true);
    toggle_class(block_container, "hide_only_scan", true);
  }
}

function action_outpop_format() {
  update_outpop_format(current_motif);
}

//
// page_loaded
//
// Called when the page has loaded for the first time.
///
function page_loaded() {
  post_load_setup();
}

//
// page_loaded
//
// Called when a cached page is reshown.
///
function page_shown(e) {
  if (e.persisted) post_load_setup();
}

//
// page_loaded
//
// Called when the page is resized
///
function page_resized() {
  var page, blocks_scroll;
  update_scroll_pad();
  page = (document.compatMode === "CSS1Compat") ? document.documentElement : document.body;
  blocks_scroll = $("blocks_scroll");
  if (blocks_scroll) {
    blocks_scroll.style.maxHeight = Math.max(page.clientHeight - 300, 300) + "px";
  }
}

//
// pre_load_setup
//
// Run before the page is displayed
///
function pre_load_setup() {
  var start, hue, sat, light, divisions;
  var i, j, motifs, motif, sites, site, sequences, sequence;
  var max_seq_len;
  motifs = data["motifs"];
  sequences = data["sequence_db"]["sequences"];
  max_seq_len = 1;
  for (i = 0; i < sequences.length; i++) {
    sequence = sequences[i];
    sequence["sites"] = [];
    if (sequence["length"] > max_seq_len) {
      max_seq_len = sequence["length"];
    }
  }
  data["sequence_db"]["max_length"] = max_seq_len;
  // use hsl colours
  start = 0; //red
  sat = 100;
  light = 50;
  for (i = 0; i < motifs.length; i++) {
    motif = motifs[i];
    // give the motif a colour
    divisions = 1 << Math.ceil(Math.log(i + 1) / Math.LN2);
    hue = start + (360 / divisions) * ((i - (divisions >> 1)) * 2 + 1);
    motif["colour"] = "hsl(" + hue + ", " + sat + "%, " + light + "%)";
    // associate sites with sequences as well 
    // to make generating the block diagram easier
    sites = motif["sites"];
    for (j = 0; j < sites.length; j++) {
      site = sites[j];
      sequence = sequences[site["seq"]];
      // record the motif index
      site["motif"] = i;
      // add the site to the sequence
      sequence["sites"].push(site);
    }
  }
}

//
// post_load_setup
//
// Run when the page has loaded, or been reloaded.
//
function post_load_setup() {
  update_scroll_pad();
  if (data["motifs"].length > 0) {
    make_motifs();
    make_blocks();
  } else {
    $("motifs").innerHTML = "<p>No significant motifs found!</p>"; // clear content
    $("motifs").innerHTML += "<p><b>" + data["stop_reason"] + "</b></p>";
    $("blocks").innerHTML = "<p>No significant motifs found!</p>";
  }
}

pre_load_setup();
