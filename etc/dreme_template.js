var current_motif = 0;
var dreme_alphabet = new Alphabet(data.alphabet, data.control_db.freqs);

/*
 * Create a pspm for the given motif data
 */
function motif_pspm(m) {
  return new Pspm(m.pwm, m.id, 0, 0, m.nsites, m.evalue);
}

/*
 * Create a count matrix from the given motif data
 */
function motif_count_matrix(motif) {
  return motif_pspm(motif).as_count_matrix();
}

/*
 * Create a probablity matrix from the given motif data
 */
function motif_prob_matrix(motif) {
  return motif_pspm(motif).as_probability_matrix();
}

/*
 * Create a minimal meme format motif from the given motif data
 */
function motif_minimal_meme(motif) {
  return motif_pspm(motif).as_meme({
    "with_header": true, 
    "with_pspm": true,
    "with_pssm": false,
    "version": data["version"],
    "alphabet": dreme_alphabet,
    "strands": (data.options.revcomp ? 2 : 1)
  });
}

/*
 * Fill in a template variable
 */
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

/*
 * Make a canvas with the motif logo drawn on it. 
 */
function make_logo(motif, height, rc) {
  var pspm = new Pspm(motif["pwm"]);
  if (rc) pspm = pspm.copy().reverse_complement(dreme_alphabet);
  var logo = new Logo(dreme_alphabet);
  logo.add_pspm(pspm, 0);
  var canvas = document.createElement('canvas');
  canvas.height = height;
  canvas.width = 0;
  draw_logo_on_canvas(logo, canvas, false);
  return canvas;
}

/*
 * Create a button designed to contain a single symbol
 */
function make_sym_btn(symbol, title, action) {
  var box, sbox;
  box = document.createElement("div");
  box.tabIndex = 0;
  box.className = "sym_btn";
  sbox = document.createElement("span");
  if (typeof symbol == "string") {
    sbox.appendChild(document.createTextNode(symbol));
  } else {
    sbox.appendChild(symbol);
  }
  box.appendChild(sbox);
  box.title = title;
  box.addEventListener('click', action, false);
  box.addEventListener('keydown', action, false);
  return box;
}

/*
 * Create a pair of text spans with different classes.
 * This is useful when using CSS to only display one of them.
 */
function text_pair(txt1, cls1, txt2, cls2) {
  var container, part1, part2;
  container = document.createElement("span");
  part1 = document.createElement("span");
  part1.appendChild(document.createTextNode(txt1));
  part1.className = cls1;
  container.appendChild(part1);
  part2 = document.createElement("span");
  part2.appendChild(document.createTextNode(txt2));
  part2.className = cls2;
  container.appendChild(part2);
  return container;
}

/*
 * Make a colourised sequence.
 */
function make_seq(seq) {
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
    lbox.style.color = dreme_alphabet.get_colour(dreme_alphabet.get_index(letter));
    lbox.appendChild(document.createTextNode(seq.substring(i, j)));
    sbox.appendChild(lbox);
  }
  return sbox;
}

/*
 * Create a description element taking into account the newlines in the source text.
 */
function make_description(text) {
  var i, j, lines, p;
  var container = document.createElement("div");
  var paragraphs = text.split(/\n\n+/);
  for (i = 0; i < paragraphs.length; i++) {
    lines = paragraphs[i].split(/\n/);
    p = document.createElement("p");
    p.appendChild(document.createTextNode(lines[0]));
    for (j = 1; j < lines.length; j++) {
      p.appendChild(document.createElement("br"));
      p.appendChild(document.createTextNode(lines[j]));
    }
    container.appendChild(p);
  }
  return container;
}

/*
 * Make the table header for the discovered motifs.
 */
function make_motif_header() {
  var row = document.createElement("tr");
  add_text_header_cell(row, "", "", "motif_ordinal");
  add_text_header_cell(row, "Motif", "pop_motifs_word", "motif_word");
  add_text_header_cell(row, "Logo", "pop_motifs_logo", "motif_logo");
  if (data.options.revcomp) {
    add_text_header_cell(row, "RC Logo", "pop_motifs_rc_logo", "motif_logo");
  }
  add_text_header_cell(row, "E-value", "pop_motifs_evalue", "motif_evalue");
  add_text_header_cell(row, "Unerased E-value", "pop_motifs_uevalue", "motif_evalue");
  add_text_header_cell(row, "More", "pop_more", "motif_more");
  add_text_header_cell(row, "Submit/Download", "pop_submit_dl", "motif_submit");
  row.className = "more";
  return row;
}

/*
 * Make a compact motif summary row for the discovered motifs.
 */
function make_motif_row(tbody, ordinal, motif) {
  var row = document.createElement("tr");
  add_text_cell(row, "" + ordinal + ".", "motif_ordinal");
  add_text_cell(row, motif["id"], "motif_word");
  add_cell(row, make_logo(motif, 50, false), "motif_logo");
  if (data.options.revcomp) {
    add_cell(row, make_logo(motif, 50, true), "motif_logo");
  }
  add_text_cell(row, motif["evalue"], "motif_evalue");
  add_text_cell(row, motif["unerased_evalue"], "motif_evalue");
  add_cell(row, make_sym_btn(text_pair("\u21A7", "less", "\u21A5", "more"), "Show more information.", function(e) { toggle_class(tbody, "collapsed"); }, "\u21A5", ""), "motif_more");
  add_cell(row, make_sym_btn("\u21E2", "Submit the motif to another MEME Suite program or download it.", function(e) { action_show_outpop(e, ordinal); }), "motif_submit");
  return row;
}

/*
 * Make a sortable table of enriched matching rows.
 */
function make_motif_words(motif) {
  var row, i, match;
  var table = document.createElement("table");
  var thead = document.createElement("thead");
  row = document.createElement("tr");
  add_text_header_cell(row, "Word", "pop_match_word", "match_word", function(e) {sort_table(this, compare_words);});
  add_text_header_cell(row, "Positives", "pop_match_pos", "match_count", function(e) {sort_table(this, compare_counts);});
  add_text_header_cell(row, "Negatives", "pop_match_neg", "match_count", function(e) {sort_table(this, compare_counts);});
  add_text_header_cell(row, "P-value", "pop_match_pval", "match_evalue", function(e) {sort_table(this, compare_evalues);});
  add_text_header_cell(row, "E-value", "pop_match_eval", "match_evalue", function(e) {sort_table(this, compare_evalues);});
  thead.appendChild(row);
  table.appendChild(thead);
  var tbody = document.createElement("tbody");
  for (i = 0; i < motif.matches.length; i++) {
    match = motif.matches[i];
    row = document.createElement("tr");
    add_cell(row, make_seq(match.seq), "match_word");
    add_text_cell(row, match.p + " / " + data.sequence_db.count, "match_count");
    add_text_cell(row, match.n + " / " + data.control_db.count, "match_count");
    add_text_cell(row, match.pvalue, "match_evalue");
    add_text_cell(row, match.evalue, "match_evalue");
    tbody.appendChild(row);
  }
  table.appendChild(tbody);
  return table;
}

/*
 * Make an expanded view of a discovered motif.
 */
function make_motif_exp(tbody, ordinal, motif) {
  "use strict";
  var box, pspm, logo_box;
  box = $("tmpl_motif_expanded").cloneNode(true);
  toggle_class(box, "template", false);
  box.id = "";
  find_child(box, "tvar_logo").appendChild(make_logo(motif, 150, false));
  if (data.options.revcomp) {
    find_child(box, "tvar_rclogo").appendChild(make_logo(motif, 150, true));
  }
  set_tvar(box, "tvar_p", motif["p"]);
  set_tvar(box, "tvar_p_total", data.sequence_db.count);
  set_tvar(box, "tvar_n", motif["n"]);
  set_tvar(box, "tvar_n_total", data.control_db.count);
  set_tvar(box, "tvar_pvalue", motif["pvalue"]);
  set_tvar(box, "tvar_evalue", motif["evalue"]);
  set_tvar(box, "tvar_uevalue", motif["unerased_evalue"]);
  set_tvar(box, "tvar_words", make_motif_words(motif));
  var cell = document.createElement("td");
  cell.colSpan = 8;
  cell.appendChild(box);
  var row = document.createElement("tr");
  row.className = "more";
  row.appendChild(cell);
  return row;
}

/*
 * Convert a string containing a scientific number into the log of that number
 * without having an intermediate representation of the number.
 * This is intended to avoid underflow problems with the tiny evalues that
 * MEME and DREME can create.
 */
function sci2log(scinum) {
  "use strict";
  var ev_re, match, sig, exp;
  ev_re = /^(.*)e(.*)$/;
  if (match = ev_re.exec(scinum)) {
    sig = parseFloat(match[1]);
    exp = parseInt(match[2]);
    return Math.log(sig) + (exp * Math.log(10));
  }
  return 0;
}

/*
 * Create a table of discovered motifs. A fresh table body is used for each
 * motif to make hiding/showing rows with css easier.
 */
function make_motifs() {
  "use strict";
  var i, row, tbody, motif, ordinal;
  // make the motifs table
  var container = $("motifs");
  container.innerHTML = ""; // clear content
  var table = document.createElement("table");
  // add a header that is always shown
  var thead = document.createElement("thead");
  thead.appendChild(make_motif_header());
  table.appendChild(thead);
  for (i = 0; i < data.motifs.length; i++) {
    ordinal = i + 1;
    motif = data.motifs[i];
    tbody = document.createElement("tbody");
    tbody.className = "collapsed";
    tbody.appendChild(make_motif_row(tbody, ordinal, motif));
    tbody.appendChild(make_motif_exp(tbody, ordinal, motif));
    // create a following header for every row except the last one
    if ((i + 1) < data.motifs.length) tbody.appendChild(make_motif_header());
    table.appendChild(tbody);
  }
  container.appendChild(table);
}

/*
 * Create a table showing all the alphabet symbols, their names and frequencies.
 */
function make_alpha_bg(alph, freqs) {
  function colour_symbol(index) {
    var span = document.createElement("span");
    span.appendChild(document.createTextNode(alph.get_symbol(index)));
    span.style.color = alph.get_colour(index);
    span.className = "alpha_symbol";
    return span;
  }
  var table, thead, tbody, row, th, span, i;
  // create table
  table = document.createElement("table");
  table.className = "inputs";
  // create header
  thead = document.createElement("thead");
  table.appendChild(thead);
  row = thead.insertRow(thead.rows.length);
  if (alph.has_complement()) {
    add_text_header_cell(row, "Name", "pop_alph_name");
    add_text_header_cell(row, "Bg.", "pop_alph_control");
    add_text_header_cell(row, "");
    add_text_header_cell(row, "");
    add_text_header_cell(row, "");
    add_text_header_cell(row, "Bg.", "pop_alph_control");
    add_text_header_cell(row, "Name", "pop_alph_name");
  } else {
    add_text_header_cell(row, "");
    add_text_header_cell(row, "Name", "pop_alph_name");
    add_text_header_cell(row, "Bg.", "pop_alph_control");
  }
  // add alphabet entries
  tbody = document.createElement("tbody");
  table.appendChild(tbody);
  if (alph.has_complement()) {
    for (i = 0; i < alph.get_size_core(); i++) {
      var c = alph.get_complement(i);
      if (i > c) continue;
      row = tbody.insertRow(tbody.rows.length);
      add_text_cell(row, alph.get_name(i));
      add_text_cell(row, "" + freqs[i].toFixed(3));
      add_cell(row, colour_symbol(i)); 
      add_text_cell(row, "~");
      add_cell(row, colour_symbol(c)); 
      add_text_cell(row, "" + freqs[c].toFixed(3));
      add_text_cell(row, alph.get_name(c));
    }
  } else {
    for (i = 0; i < alph.get_size_core(); i++) {
      row = tbody.insertRow(tbody.rows.length);
      add_cell(row, colour_symbol(i)); 
      add_text_cell(row, alph.get_name(i));
      add_text_cell(row, "" + freqs[i].toFixed(3));
    }
  }
  return table;
}

/*
 * Updates the format download text in the popup.
 * This is called when either the format or current motif changes.
 */
function update_outpop_format(index) {
  var motif = data.motifs[index];
  var fn = [motif_count_matrix, motif_prob_matrix, motif_minimal_meme];
  var suffix = ["_counts.txt", "_freqs.txt", ".meme"];
  var format = parseInt($("text_format").value);
  var text = fn[format](motif);
  prepare_download(text, "text/plain", motif.id + suffix[format], $("outpop_text_dl"));
  $("outpop_text").value = text;
}

/*
 * Updates the motif logos and format download text in the popup.
 * This is called whenever the current motif changes.
 */
function update_outpop_motif(index) {
  "use strict";
  var motifs, motif, pspm, logo, canvas, num;
  motifs = data["motifs"];
  if (index < 0 || index >= motifs.length) {return;}
  current_motif = index;
  motif = motifs[index];
  pspm = new Pspm(motif["pwm"]);
  logo = new Logo(dreme_alphabet, "");
  logo.add_pspm(pspm, 0);
  canvas = $("outpop_logo");
  canvas.width = canvas.width; // clear canvas
  draw_logo_on_canvas(logo, canvas, false);
  canvas = $("outpop_logo_rc");
  canvas.width = canvas.width; // clear rc canvas
  if (data.options.revcomp) {
    pspm.reverse_complement(dreme_alphabet);
    logo = new Logo(dreme_alphabet, "");
    logo.add_pspm(pspm, 0);
    draw_logo_on_canvas(logo, canvas, false);
  }
  num = $("outpop_num");
  num.innerHTML = "";
  num.appendChild(document.createTextNode("" + (index + 1)));
  update_outpop_format(index);
}


/*
 * Initialise and display the download popup.
 */
function action_show_outpop(e, ordinal) {
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
    toggle_class(tbl_submit, "alphabet_dna", dreme_alphabet.has_complement());//TODO FIXME alphabet_dna is a bad name for a field when allowing custom alphabets
    // add a click listener for the radio buttons
    inputs = tbl_submit.querySelectorAll("input[type='radio']");
    for (i = 0; i < inputs.length; i++) {
      inputs[i].addEventListener("click", action_outpop_program, false);
    }
    // ensure that a default program option is selected for DNA and Protein
    default_prog = document.getElementById(dreme_alphabet.has_complement() ? "submit_tomtom" : "submit_fimo");
    default_prog.checked = true;
    action_outpop_program.call(default_prog);
    // disable reverse-complement when not DNA
    $("logo_rc_option").disabled = !dreme_alphabet.has_complement(); 
    // set errorbars on when ssc is on
    $("logo_ssc").addEventListener("change", action_outpop_ssc, false);
  }
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
  update_outpop_motif(ordinal - 1);
  // display the download popup
  $("grey_out_page").style.display = "block";
  $("download").style.display = "block";
  $("outpop_close").focus();
}

/*
 * Hide the download popup.
 */
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

/*
 * Show the next motif in the download popup.
 */
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

/*
 * Show the previous motif in the download popup.
 */
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

/*
 * Highlight the selected row in the program list.
 */
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

/*
 * Enable error bars when small sample correction is enabled.
 */
function action_outpop_ssc() {
  "use strict";
  $("logo_err").value = $("logo_ssc").value;
}

/*
 * Submit the motif to the selected program.
 */
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

  motifs = motif_minimal_meme(data.motifs[current_motif]);
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

/*
 * Download the format text.
 * Wire the link containing the data URI text to a download button so it looks
 * the same as the server submit stuff.
 */
function action_outpop_download_motif(e) {
  $("outpop_text_dl").click();
}

/*
 * Download the motif logo.
 * The EPS format can be calculated locally in Javascript
 */
function action_outpop_download_logo(e) {
  "use strict";
  var pspm, logo, eps;
  var logo_rc, logo_ssc, logo_width, logo_height;
  var motif = data.motifs[current_motif];
  if ($("logo_format").value == "0") { // EPS
    logo_rc = ($("logo_rc").value == "1");
    logo_ssc = ($("logo_ssc").value == "1");
    logo_width = parseFloat($("logo_width").value);
    if (isNaN(logo_width) || !isFinite(logo_width) || logo_width <= 0) logo_width = null;
    logo_height = parseFloat($("logo_height").value);
    if (isNaN(logo_height) || !isFinite(logo_height) || logo_height <= 0) logo_height = null;
    // create a PSPM from the motif
    pspm = motif_pspm(motif);
    if (logo_rc) pspm.reverse_complement(dreme_alphabet);
    logo = new Logo(dreme_alphabet);
    logo.add_pspm(pspm, 0);
    eps = logo.as_eps({"ssc": logo_ssc, "logo_width": logo_width, "logo_height": logo_height});
    prepare_download(eps, "application/postscript", motif.id + ".eps");
  } else {
    $("logo_motifs").value = motif_minimal_meme(motif);
    $("logo_form").submit();
  }
}

/*
 * Change the selected tab in the download popup.
 */
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

/*
 * Update the text in the download format popup.
 */
function action_outpop_format() {
  update_outpop_format(current_motif);
}

/*
 * Find all text nodes in the given container.
 */
function text_nodes(container) {
  var textNodes = [];
  var stack = [container];
  // depth first search to maintain ordering when flattened 
  while (stack.length > 0) {
    var node = stack.pop();
    if (node.nodeType == Node.TEXT_NODE) {
      textNodes.push(node);
    } else {
      for (var i = node.childNodes.length-1; i >= 0; i--) {
        stack.push(node.childNodes[i]);
      }
    }
  }
  return textNodes;
}

/*
 * Get the text out of a specific text node.
 */
function node_text(node) {
  if (node === undefined) {
    return '';
  } else if (node.textContent) {
    return node.textContent;
  } else if (node.innerText) {
    return node.innerText;
  } else {
    return '';
  }
}

/*
 * Get the text contained within the element.
 */
function elem_text(elem, separator) {
  if (separator === undefined) separator = '';
  return text_nodes(elem).map(node_text).join(separator);
}

/*
 * Sort all rows in the first table body based on the column of the given element and the comparison function.
 * The sort is not very fast and is intended for small tables only.
 */
function sort_table(colEle, compare_function) {
  //find the parent of colEle that is either a td or th
  var i, j;
  var cell = colEle;
  while (true) {
    if (cell == null) return;
    if (cell.nodeType == Node.ELEMENT_NODE && 
        (cell.tagName.toLowerCase() == "td" || cell.tagName.toLowerCase() == "th")) {
      break;
    }
    cell = cell.parentNode;
  }
  //find the parent of cell that is a tr
  var row = cell;
  while (true) {
    if (row == null) return;
    if (row.nodeType == Node.ELEMENT_NODE && row.tagName.toLowerCase() == "tr") {
      break;
    }
    row = row.parentNode;
  }
  //find the parent of row that is a table
  var table = row;
  while (true) {
    if (table == null) return;
    if (table.nodeType == Node.ELEMENT_NODE && table.tagName.toLowerCase() == "table") {
      break;
    }
    table = table.parentNode;
  }
  var column_index = cell.cellIndex;
  // do a bubble sort, because the tables are so small it doesn't matter
  var change;
  var trs = table.tBodies[0].getElementsByTagName('tr');
  var already_sorted = true;
  var reverse = false;
  while (true) {
    do {
      change = false;
      for (i = 0; i < trs.length -1; i++) {
        var v1 = elem_text(trs[i].cells[column_index]);
        var v2 = elem_text(trs[i+1].cells[column_index]);
        if (reverse) {
          var tmp = v1;
          v1 = v2;
          v2 = tmp;
        }
        if (compare_function(v1, v2) > 0) {
          exchange(trs[i], trs[i+1], table);
          change = true;
          already_sorted = false;
          trs = table.tBodies[0].getElementsByTagName('tr');
        }
      }
    } while (change);
    if (reverse) break;// we've sorted twice so exit
    if (!already_sorted) break;// sort did something so exit
    // when it's sorted one way already then sort the opposite way
    reverse = true;
  }
  // put arrows on the headers
  var headers = table.tHead.getElementsByTagName('tr');
  for (i = 0; i < headers.length; i++) {
    for (j = 0; j < headers[i].cells.length; j++) {
      var cell = headers[i].cells[j];
      var arrows = cell.getElementsByClassName("sort_arrow");
      var arrow;
      if (arrows.length == 0) {
        arrow = document.createElement("span");
        arrow.className = "sort_arrow";
        cell.insertBefore(arrow, cell.firstChild);
      } else {
        arrow = arrows[0];
      }
      arrow.innerHTML = "";
      if (j == column_index) {
        arrow.appendChild(document.createTextNode(reverse ? "\u25B2" : "\u25BC"));
      }
    }
  }
}

/*
 * Swap two rows in a table.
 */
function exchange(oRowI, oRowJ, oTable) {
  var i = oRowI.rowIndex;
  var j = oRowJ.rowIndex;
   if (i == j+1) {
    oTable.tBodies[0].insertBefore(oRowI, oRowJ);
  } if (j == i+1) {
    oTable.tBodies[0].insertBefore(oRowJ, oRowI);
  } else {
    var tmpNode = oTable.tBodies[0].replaceChild(oRowI, oRowJ);
    if(typeof(oRowI) != "undefined") {
      oTable.tBodies[0].insertBefore(tmpNode, oRowI);
    } else {
      oTable.appendChild(tmpNode);
    }
  }
}

/*
 * Compare two E-values which may be very small.
 */
function compare_evalues(v1, v2) {
  var e1 = sci2log(v1);
  var e2 = sci2log(v2);
  if (e1 < e2) return -1;
  else if (e1 > e2) return 1;
  return 0;
}

/*
 * Compare two counts.
 */
function compare_counts(v1, v2) {
  var re = /(\d+)\s*\/\s*\d+/;
  var m1 = re.exec(v1);
  var m2 = re.exec(v2);
  if (m1 == null && m2 == null) return 0;
  if (m1 == null) return -1;
  if (m2 == null) return 1;
  return parseInt(m2[1]) - parseInt(m1[1]);
}

/*
 * Compare two sequence words.
 */
function compare_words(v1, v2) {
  return v1.localeCompare(v2);
}


