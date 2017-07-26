var mcast_alphabet = new Alphabet(data.alphabet, data.background.freqs);

function redraw_annotated_sequence(tbody, first_display) {
  var annoboxcont, num_size, annobox, linelen, seqidx, maxlen, match, range;
  annoboxcont = tbody.querySelector(".annotated_sequence_container");
  num_size = $("num_ruler").offsetWidth + 20;
  linelen = Math.floor((annoboxcont.offsetWidth - (2 * num_size)) / $("ruler").offsetWidth * 5);
  annobox = tbody.querySelector(".annotated_sequence");
  annobox.style.left = num_size + "px";
  seqidx = parseInt(annobox.getAttribute("data-seqidx"), 10);
  if (seqidx < data.matches.length) {
    match = data.matches[seqidx];
  } else if (typeof data_aux != "undefined") {
    match = data_aux.matches[seqidx - data.matches.length];
  }
  if (match == null) {
    window.console && console.log("Unable to access match idx: " + seqidx);
    return;
  }
  if (first_display) {
    var start = match.segs[0].pos + match.start;
    var last_seg = match.segs[match.segs.length - 1];
    var end = last_seg.pos + last_seg.data.length + match.start;
    range = {"start": start, "end": end};
    set_block_needle_positions(tbody, range.start, range.end);
  } else {
    range = get_block_needle_positions(tbody);
  }
  annobox.innerHTML = "";
  annobox.appendChild(make_annotated_sequence(match, 0, range.start - match.start, range.end - match.start, linelen));
  annoboxcont.style.height = (annobox.offsetHeight + 10) + "px";
}

function redraw_annotated_sequences() {
  var all_expanded = document.querySelectorAll("table.seqs_table tbody.more");
  var i;
  for (i = 0; i < all_expanded.length; i++) {
    redraw_annotated_sequence(all_expanded[i], false);
  }
}

function calculate_wrap() {
  var page_width;
  if (window.innerWidth) {
    page_width = window.innerWidth;
  } else if (document.body) {
    page_width = document.body.clientWidth;
  } else {
    page_width = null;
  }
  var est_wrap = 0;
  if (page_width) {
    var ruler_width = document.getElementById("ruler").offsetWidth;
    est_wrap = Math.floor(page_width / (ruler_width / 5));
  }
  return est_wrap;
}

function rewrap() {
  var new_wrap = calculate_wrap();
  if (new_wrap != rewrap.wrap) {
    rewrap.wrap = new_wrap;
    redraw_annotated_sequences();
  }
}

function delayed_rewrap() {
  if (delayed_rewrap.timer != null) {
    clearTimeout(delayed_rewrap.timer);
  }
  delayed_rewrap.timer = setTimeout(rewrap, 1000);
}

/*
 * Make a colourised match.
 */
function make_seq(seq, alphabet) {
  var i, j, letter, letter_i, lbox, sbox;
  if (alphabet == null) alphabet = mcast_alphabet;
  sbox = document.createElement("span");
  for (i = 0; i < seq.length; i = j) {
    letter = seq.charAt(i);
    letter_i = mcast_alphabet.get_index(letter);
    for (j = i+1; j < seq.length; j++) {
      if (mcast_alphabet.get_index(seq.charAt(j)) !== letter_i) {
        break;
      }
    }
    if (letter_i != -1)  {
      lbox = document.createElement("span");
      lbox.style.color = mcast_alphabet.get_colour(letter_i);
      lbox.appendChild(document.createTextNode(seq.substring(i, j)));
      sbox.appendChild(lbox);
    } else {
      sbox.appendChild(document.createTextNode(seq.substring(i, j)));
    }
  }
  return sbox;
}

function make_sequence_databases_table() {
  var i, seq_sum, res_sum, db_stats, db_summary, row, db, name;
  seq_sum = 0; res_sum = 0;
  db_stats = $("sequence_db_stats");
  db_stats.innerHTML = "";
  db_summary = $("sequence_db_summary");
  db_summary.innerHTML = "";
  for (i = 0; i < data.sequence_dbs.length; i++) {
    row = document.createElement("tr");
    db = data.sequence_dbs[i];
    name = (db.name != null ? db.name : db.file);
    seq_sum += db.sequence_count;
    res_sum += db.residue_count;
    add_text_cell(row, name);
    add_text_cell(row, db.psp_file || "-", "priors");
    add_text_cell(row, db.dist_file || "-", "priors");
    add_text_cell(row, db.sequence_count);
    add_text_cell(row, db.residue_count);
    db_stats.appendChild(row);
  }
  row = document.createElement("tr");
  add_text_header_cell(row, "Total");
  add_text_cell(row, "", "priors");
  add_text_cell(row, "", "priors");
  add_text_cell(row, seq_sum);
  add_text_cell(row, res_sum);
  db_summary.appendChild(row);
  toggle_class($("sequence_db"), "has_priors", data.sequence_dbs.some(function (db) { return db.psp_file != null || db.dist_file != null; }));
}

function make_motif_databases_table() {
  var i, db_stats, row, db, name;
  db_stats = $("motif_db_stats");
  db_stats.innerHTML = "";
  for (i = 0; i < data.motif_dbs.length; i++) {
    row = document.createElement("tr");
    db = data.motif_dbs[i];
    name = (db.name != null ? db.name : db.file);
    add_text_cell(row, name);
    db_stats.appendChild(row);
  }
}

function action_btn_rc(e) {
  "use strict";
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  var rc = /\bminus\b/.test(this.className);
  var box = find_parent(this, "preview_box");
  toggle_class(box, "rc", rc);
}

function make_preview_logo(motif, alphabet, rc) {
  var pspm = new Pspm(motif.pwm);
  if (rc) pspm.reverse_complement(alphabet);
  var logo = new Logo(alphabet);
  logo.add_pspm(pspm);
  var canvas = document.createElement('canvas');
  canvas.className = "motif_box " + (rc ? "minus" : "plus");
  canvas.height = 50;
  canvas.width = 0;
  draw_logo_on_canvas(logo, canvas, false);
  return canvas;
}

function make_preview(alphabet, motif) {
  "use strict";
  if (alphabet.has_complement()) {
    var btn_plus = document.createElement("div");
    btn_plus.appendChild(document.createTextNode("+"));
    btn_plus.className = "preview_btn plus";
    btn_plus.tabIndex = "0";
    btn_plus.addEventListener("click", action_btn_rc, false);
    btn_plus.addEventListener("keydown", action_btn_rc, false);
    var btn_minus = document.createElement("div");
    btn_minus.appendChild(document.createTextNode("-"));
    btn_minus.className = "preview_btn minus";
    btn_minus.tabIndex = "0";
    btn_minus.addEventListener("click", action_btn_rc, false);
    btn_minus.addEventListener("keydown", action_btn_rc, false);
    var btn_box = document.createElement("div");
    btn_box.className = "preview_btn_box";
    btn_box.appendChild(btn_plus);
    btn_box.appendChild(btn_minus);
  }
  var logo_box = document.createElement("div");
  logo_box.className = "preview_logo_box";
  logo_box.appendChild(make_preview_logo(motif, alphabet, false));
  if (alphabet.has_complement()) {
    logo_box.appendChild(make_preview_logo(motif, alphabet, true));
  }
  var box = document.createElement("div");
  box.className = "preview_box";
  if (alphabet.has_complement()) box.appendChild(btn_box);
  box.appendChild(logo_box);
  return box;
}

function make_motif_table_entry(row, index, motif) {
  "use strict";
  add_text_header_cell(row, "" + (index + 1) + ".", "", "motif_ordinal");
  add_cell(row, make_preview(mcast_alphabet, motif), "motif_logo");
  add_text_cell(row, motif.id, "motif_id");
  add_text_cell(row, (motif.alt != null ? motif.alt : "") , "motif_alt");
  add_text_cell(row, motif.len, "motif_width");
}

function make_motifs_table() {
  var i, j;
  var tbl, thead, tbody, row;

  tbl = document.createElement("table");
  tbl.className = "motifs" + (data.settings.remove_correlated ? " remove_correlated" : "");
  
  thead = document.createElement("thead");
  tbl.appendChild(thead);
  tbody = document.createElement("tbody");
  tbl.appendChild(tbody);
  row = thead.insertRow(thead.rows.length);
  add_text_header_cell(row, "", "", "motif_ordinal");
  add_text_header_cell(row, "Logo", "", "motif_logo");
  add_text_header_cell(row, "Name", "pop_motif_id", "motif_id");
  add_text_header_cell(row, "Alt. Name", "pop_motif_alt", "motif_alt");
  add_text_header_cell(row, "Width", "pop_motif_width", "motif_width");

  for (i = 0; i < data.motifs.length; i++) {
    row = tbody.insertRow(tbody.rows.length);
    make_motif_table_entry(row, i, data.motifs[i]);
  }

  return tbl;
}

function make_other_settings() {
  $("opt_alpha").textContent = data.settings.alpha;
  $("opt_hard_mask").textContent = (data.settings.hard_mask ? "Yes" : "No");
  $("opt_max_gap").textContent = data.settings.max_gap;
  $("opt_max_stored_scores").textContent = data.settings.max_stored_scores;
  $("opt_motif_p_thresh").textContent = data.settings.motif_p_thresh;
  var thresh_type, thresh_value;
  if (data.settings.p_thresh != null) {
    thresh_type = "p-value";
    thresh_value = data.settings.p_thresh;
  } else if (data.settings.e_thresh != null) {
    thresh_type = "E-value";
    thresh_value = data.settings.e_thresh;
  } else if (data.settings.q_thresh != null) {
    thresh_type = "q-value";
    thresh_value = data.settings.q_thresh;
  }
  $("opt_threshold_type").textContent = thresh_type;
  $("opt_threshold_value").textContent = thresh_value;
  $("opt_parse_genomic_coord").textContent = (data.settings.genomic_coord ? "Yes" : "No");
  $("opt_synth").textContent = (data.settings.synth ? "Yes with seed " + data.settings.seed : "No");

  // adv opts (no way of setting via mcast command line)
  $("opt_min_match_score").textContent = data.settings.min_match_score;
  $("opt_cost_factor").textContent = data.settings.cost_factor;
  $("opt_gap_open_cost").textContent = data.settings.gap_open_cost;
  $("opt_gap_extend_cost").textContent = data.settings.gap_extend_cost;
}

function clone_template(template) {
  "use strict";
  var node, help_btns, i, button, new_button;
  node = $(template).cloneNode(true);
  toggle_class(node, "template", false);
  node.id = "";
  // replace help buttons (because they get broken)
  help_btns = node.querySelectorAll(".help");
  for (i = 0; i < help_btns.length; i++) {
    button = help_btns[i];
    if (button.hasAttribute("data-topic")) {
      new_button = help_button(button.getAttribute("data-topic"));
      button.parentNode.replaceChild(new_button, button);
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

function make_hit_popup(container, block, match, seg, hit) {
  return function(evt) {
    "use strict";
    var pop, xy, padding, edge_padding, pop_left, pop_top, page_width;
    var motif, pspm, parts;
    if (!evt) var evt = window.event;
    pop = make_hit_popup.pop;
    if (evt.type === "mouseover") {
      if (pop) return;
      pop = clone_template("t_hit_info");
      motif = data.motifs[hit.idx];
      var hit_len = motif.len;
      pspm = new Pspm(motif.pwm);
      if (hit.rc) {
        pspm.reverse_complement(mcast_alphabet);
      }
      var logo = new Logo(mcast_alphabet, {x_axis_hidden: true, y_axis: false});
      logo.add_pspm(pspm);
      /*
      */
      var fsize = 10;
      var lflank = seg.data.substring(Math.max(0, hit.pos - seg.pos - fsize), hit.pos - seg.pos);
      if (fsize > lflank.length && hit.pos > lflank.length) {
        lflank = string_mult("-", Math.min(hit.pos, fsize) - lflank.length) + lflank;
      }
      var seq = seg.data.substr(hit.pos - seg.pos, hit_len);
      var rflank = seg.data.substr(hit.pos - seg.pos + hit_len, fsize);
      if (fsize > rflank.length && match.length > (hit.pos + hit_len + rflank.length)) {
        rflank += string_mult("-", Math.min(match.length - (hit.pos + hit_len), fsize) - rflank.length);
      }
      var canvas = document.createElement('canvas');
      canvas.className = "motif_box";
      canvas.height = 0;
      canvas.width = hit_len * ($("ruler").offsetWidth / 5);
      draw_logo_on_canvas(logo, canvas, false);
      set_tvar(pop, "tvar_logo_pad", lflank);
      set_tvar(pop, "tvar_logo", canvas);
      pop.querySelector("div.xlate").style.display = "none";
      set_tvar(pop, "tvar_lflank", lflank);
      set_tvar(pop, "tvar_hit", make_seq(seq));
      set_tvar(pop, "tvar_rflank", rflank);
      set_tvar(pop, "tvar_motif", motif.id);
      set_tvar(pop, "tvar_pvalue", hit.pvalue.toExponential(1));
      set_tvar(pop, "tvar_start", match.start + hit.pos + 1);
      set_tvar(pop, "tvar_end", match.start + hit.pos + hit_len);
      document.body.appendChild(pop);
      position_popup(block, pop);
      make_hit_popup.pop = pop;
    } else if (evt.type === "mouseout") {
      if (pop) {
        pop.parentNode.removeChild(pop);
        make_hit_popup.pop = null;
      }
    }
  };
}

function make_seq_diagram(max_cluster_len, display_len, display_offset, match, range_handler) {
  var seg_i, seg, hit_i, hit, motif, colour_i, hit_len, elems;

  var container = make_block_container(
      mcast_alphabet.has_complement(),
      data.settings.norc,
      max_cluster_len, display_len, match.start + display_offset,
      range_handler
      );
  if (match.segs != null) {
    for (seg_i = 0; seg_i < match.segs.length; seg_i++) {
      seg = match.segs[seg_i];
      for (hit_i = 0; hit_i < seg.hits.length; hit_i++) {
        hit = seg.hits[hit_i];
        motif = data.motifs[hit.idx];
        hit_len = motif.pwm.length;
        block = make_block(container, max_cluster_len,
            hit.pos - display_offset, hit_len, hit.pvalue, hit.rc, motif.colour_index, false);
        popup = make_hit_popup(container, block, match, seg, hit);
        block.addEventListener("mouseover", popup, false);
        block.addEventListener("mouseout", popup, false);
      }
    }
  }
  return container;
}

function make_seq_legend(container) {
  var i, id, motif;
  container.innerHTML = "";
  for (i = 0; i < data.motifs.length; i++) {
    motif = data.motifs[i];
    if (data.settings.remove_correlated && motif.bad) continue;
    id = motif.id;
    container.appendChild(make_block_legend_entry((/^\d+$/.test(id) ? "Motif " + id : id), motif.colour_index));
  }
}

function make_seq_header() {
  var row = document.createElement("tr");
  add_text_header_cell(row, "Name", (data.settings.genome_coord ? "pop_seq_name" : "pop_seq_id"));
  add_text_header_cell(row, "Start", (data.settings.genome_coord ? "pop_genomic_start" : "pop_seq_start"), "number");
  add_text_header_cell(row, "Stop", (data.settings.genome_coord ? "pop_genomic_stop" : "pop_seq_stop"), "number");
  if (data.settings.stats) {
    add_text_header_cell(row, "E-value", "pop_cluster_evalue", "number");
  } else {
    add_text_header_cell(row, "Score", "pop_cluster_score", "number");
  }
  add_text_header_cell(row, ""); // header for the arrow
  add_text_header_cell(row, "Block Diagram", "pop_diagram");
  add_text_header_cell(row, ""); // header for the spacer
  var seqs_head = document.createElement("thead");
  seqs_head.appendChild(row);
  return seqs_head;
}

function make_seq_footer(max_seq_len) {
  var seqs_foot, row;
  row = document.createElement("tr");
  add_text_cell(row, ""); // sequence
  add_text_cell(row, ""); // start
  add_text_cell(row, ""); // stop
  add_text_cell(row, ""); // evalue/score
  add_text_cell(row, ""); // arrow
  add_cell(row, make_block_ruler(max_seq_len));
  add_text_cell(row, ""); // spacer
  seqs_foot = document.createElement("tfoot");
  seqs_foot.appendChild(row);
  return seqs_foot;
}

function make_seq_row1(max_cluster_len, match) {
  var ev1, ev2;

  // calculate the cluster length of this match
  var start, end, cluster_len, last_seg;
  start = match.segs[0].pos;
  last_seg = match.segs[match.segs.length - 1];
  end = start + last_seg.data.length;
  cluster_len = end - start;
  
  // determine the offset required
  var padding, padding_left, padding_right;
  padding = max_cluster_len - cluster_len;
  padding_left = Math.min(Math.floor(padding/2), start);
  padding_right = Math.min(padding - padding_left, match.length - end);
  padding_left = Math.min(padding - padding_right, start);

  var display_len = padding_left + cluster_len + padding_right;
  var display_offset = start - padding_left;

  var row = document.createElement("tr");
  // name
  add_text_cell(row, match.name);
  // start
  add_text_cell(row, match.start + display_offset + 1, "number");
  // stop
  add_text_cell(row, match.start + display_offset + display_len, "number");
  if (data.settings.stats) {
    // evalue
    add_text_cell(row, match.evalue != null ? match.evalue.toExponential(1) : "", "number");
  } else {
    // score
    add_text_cell(row, match.score, "number");
  }
  // more widget
  var down_arrow = document.createElement("span");
  down_arrow.className = "down";
  down_arrow.appendChild(document.createTextNode("\u21A7"));
  var up_arrow = document.createElement("span");
  up_arrow.className = "up";
  up_arrow.appendChild(document.createTextNode("\u21A5"));
  var more = document.createElement("span");
  more.appendChild(down_arrow);
  more.appendChild(up_arrow);
  more.className = "more_arrow";
  more.title = "Toggle additional information";
  more.addEventListener('click', function (e) {
    var tbody = row.parentNode;
    toggle_class(tbody, "more");
    if (/\bmore\b/.test(tbody.className)) redraw_annotated_sequence(tbody, true);
  }, false);
  add_cell(row, more);
  // block diagram
  add_cell(row, make_seq_diagram(max_cluster_len, display_len, display_offset, match,
    function (start, end, done) {
      if (done) redraw_annotated_sequence(row.parentNode, false);
    }
    ), "block_td");
  // add an invisible spacer column 
  add_text_cell(row, match.start + match.length, "spacer_num");
  return row;
}

function make_seq_row2(seqidx, match) {
  var row = document.createElement("tr");
  row.className = "more";
  var infobox = clone_template("t_infobox");
  var start = match.segs[0].hits[0].pos + 1 + match.start;
  var last_seg = match.segs[match.segs.length - 1];
  var last_hit = last_seg.hits[last_seg.hits.length - 1];
  var stop = last_hit.pos + data.motifs[last_hit.idx].len + match.start;
  infobox.querySelector(".tvar_start").textContent = start;
  infobox.querySelector(".tvar_stop").textContent = stop;
  infobox.querySelector(".tvar_score").textContent = match.score;
  infobox.querySelector(".tvar_pvalue").textContent = (match.pvalue != null ? match.pvalue.toExponential(1) : '-');
  infobox.querySelector(".tvar_evalue").textContent = (match.evalue != null ? match.evalue.toExponential(1) : '-');
  infobox.querySelector(".tvar_qvalue").textContent = (match.qvalue != null ? match.qvalue.toExponential(1) : '-');
  
  var annoseq = infobox.querySelector(".annotated_sequence");
  annoseq.setAttribute("data-seqidx", seqidx);
  var container = document.createElement("td");
  container.colSpan = 7;
  container.appendChild(infobox);
  row.appendChild(container);
  return row;
}

function make_seq_row_group(maxlen, seqidx, match) {
  var row_group;
  row_group = document.createElement("tbody");
  row_group.appendChild(make_seq_row1(maxlen, match));
  row_group.appendChild(make_seq_row2(seqidx, match));
  return row_group;
}

function make_seq_table() {
  // setup the match table
  var table = $("seqs_table");
  table.innerHTML = "";
  table.appendChild(make_seq_header());
  table.appendChild(make_seq_footer(data.calc.max_cluster_len));
  make_seq_legend($("seqs_legend_top"));
  make_seq_legend($("seqs_legend_bottom"));
  // create rows for each match
  var i;
  for (i = 0; i < data.matches.length; i++) {
    table.appendChild(make_seq_row_group(data.calc.max_cluster_len, i, data.matches[i]));
  }
  window.addEventListener("resize", delayed_rewrap, false);
}

function make_seq_table_extendable() {
  if (typeof data_aux == "undefined") {
    throw new Error("This function must not be called unless all the data has been loaded!");
  }
  var table = $("seqs_table");
  var index = data.matches.length;
  var total = data.matches.length + data_aux.matches.length;
  var more_results_btn = $("more_results");
  var add_rows = function (count) {
    var end = Math.min(index + count, total);
    for (; index < end; index++) {
      table.appendChild(make_seq_row_group(data.calc.max_cluster_len, index,
            data_aux.matches[index - data.matches.length]));
    }
    more_results_btn.textContent = "Display More Results ( " + end + " / " + total + " )";
    if (end == total) more_results_btn.style.display = "none";
  };
  more_results_btn.addEventListener("click", function (e) { add_rows(100); }, false);
  more_results_btn.textContent = "Display More Results ( " + index + " / " + total + " )";
  if (index == total) more_results_btn.style.display = "none";
}

function make_annoseq_hit(match, strand, line_start, line_end, hit_callback) {
  var line = document.createElement("div");
  line.className = "sequence";
  var i, j, seg, hit, motif, hit_len, hit_end_pos, box, label;
  var pos = line_start;
  for (i = 0; i < match.segs.length; i++) {
    seg = match.segs[i];
    if (seg.pos >= line_end) break;
    if ((seg.pos + seg.data.length) <= line_start) continue;
    for (j = 0; j < seg.hits.length; j++) {
      hit = seg.hits[j];
      if (hit.pos >= line_end) break;
      if (hit.rc ? strand > 0 : strand < 0) continue;
      motif = data.motifs[hit.idx];
      hit_len = motif.len;
      if ((hit.pos + hit_len) <= line_start) continue;
      if (pos < hit.pos) {
        line.appendChild(document.createTextNode(string_mult("\u00A0", hit.pos - pos)));
        pos = hit.pos;
      }
      hit_end_pos = Math.min(hit.pos + hit_len, line_end);
      line.appendChild(hit_callback(hit, motif, pos, hit_end_pos));
      pos = hit_end_pos;
    }
  }
  return line;
}

function make_annoseq_label(match, strand, line_start, line_end, label_callback) {
  return make_annoseq_hit(match, strand, line_start, line_end, function(hit, motif, pos_start, pos_end) {
    var box = document.createElement("span");
    box.className = "center_box";
    box.appendChild(document.createTextNode(string_mult("\u00A0", pos_end - pos_start)));
    var label = document.createElement("div");
    label.className = "center_box_label";
    label.appendChild(document.createTextNode(label_callback(hit, motif)));
    box.appendChild(label);
    return box;
  });
}

function make_annoseq_names(match, strand, line_start, line_end) {
  return make_annoseq_label(match, strand, line_start, line_end, function(hit, motif) {
    var show_strand = mcast_alphabet.has_complement() && !data.settings.norc;
    return (show_strand ? "(" + (hit.rc ? "-" : "+") + ") " : "") + 
        (/^\d+$/.test(motif.id) ? "Motif " : "") + motif.id;
  });
}

function make_annoseq_pvalues(match, strand, line_start, line_end) {
  return make_annoseq_label(match, strand, line_start, line_end, function(hit, motif) {
    return hit.pvalue.toExponential(1);
  });
}

function make_annoseq_motif(match, strand, line_start, line_end) {
  var line =  make_annoseq_hit(match, strand, line_start, line_end, function(hit, motif, pos_start, pos_end) {
    var pspm = new Pspm(motif.pwm);
    if (hit.rc) {
      pspm.reverse_complement(mcast_alphabet);
    }
    var logo = new Logo(mcast_alphabet, {x_axis_hidden: true, y_axis: false,
      xlate_nsyms: 1, xlate_start: pos_start - hit.pos, xlate_end: pos_end - hit.pos});
    logo.add_pspm(pspm);
    var canvas = document.createElement('canvas');
    canvas.className = "motif_box";
    canvas.height = 0;
    canvas.width = (pos_end - pos_start) * ($("ruler").offsetWidth / 5);
    draw_logo_on_canvas(logo, canvas, false);
    return canvas;
  });
  line.className += " motif_line";
  return line;
}

function make_annoseq_sequence(match, strand, line_start, line_end) {
  var line = document.createElement("div");
  line.className = "sequence sequence_line";
  var i, j, seg, hit, motif, hit_len, joiner, grey_seq, seg_end_pos, hit_end_pos;
  var start_pos = document.createElement("div");
  start_pos.className = "line_start";
  start_pos.appendChild(document.createTextNode(match.start + line_start + 1));
  line.appendChild(start_pos);
  var pos = line_start;
  for (i = 0; i < match.segs.length; i++) {
    seg = match.segs[i];
    if (seg.pos >= line_end) break;
    if ((seg.pos + seg.data.length) <= line_start) continue;
    if (pos < seg.pos) {
      // output dashes where we don't have sequence
      joiner = document.createElement("span");
      joiner.className = "joiner";
      joiner.appendChild(document.createTextNode(string_mult("-", seg.pos - pos)));
      line.appendChild(joiner);
      pos = seg.pos;
    }
    for (j = 0; j < seg.hits.length; j++) {
      hit = seg.hits[j];
      if (hit.pos >= line_end) break;
      if (hit.rc ? strand > 0 : strand < 0) continue;
      motif = data.motifs[hit.idx];
      hit_len = motif.len;
      if ((hit.pos + hit_len) <= line_start) continue;
      if (pos < hit.pos) {
        grey_seq = document.createElement("span");
        grey_seq.className = "joiner";
        grey_seq.appendChild(document.createTextNode(seg.data.substring(pos - seg.pos,  hit.pos - seg.pos)));
        line.appendChild(grey_seq);
        pos = hit.pos;
      }
      hit_end_pos = Math.min(hit.pos + hit_len, line_end);
      line.appendChild(make_seq(seg.data.substring(pos - seg.pos, hit_end_pos - seg.pos), mcast_alphabet));
      pos = hit_end_pos;
    }
    seg_end_pos = Math.min(seg.pos + seg.data.length, line_end);
    if (pos < seg_end_pos) {
      grey_seq = document.createElement("span");
      grey_seq.className = "joiner";
      grey_seq.appendChild(document.createTextNode(seg.data.substring(pos - seg.pos, seg_end_pos - seg.pos)));
      line.appendChild(grey_seq);
      pos = seg_end_pos;
    }
  }
  if (pos < line_end) {
    // output dashes where we don't have sequence
    joiner = document.createElement("span");
    joiner.className = "joiner";
    joiner.appendChild(document.createTextNode(string_mult("-", line_end - pos)));
    line.appendChild(joiner);
    pos = line_end;
  }
  var end_pos = document.createElement("div");
  end_pos.className = "line_end";
  end_pos.appendChild(document.createTextNode(match.start + line_end));
  line.appendChild(end_pos);
  return line;
}

function make_annoseq_bounds(match, strand, line_start, line_end) {
  return make_annoseq_hit(match, strand, line_start, line_end, function(hit, motif, pos_start, pos_end) {
    var hit_str = "\\" + string_mult("_", (motif.len) - 2) + "/";
    return document.createTextNode(hit_str.substring(pos_start - hit.pos, pos_end - hit.pos));
  });
}



function make_annotated_sequence(match, strand, start, end, line_length) {
  var line_start, bunch, line_end;
  var container = document.createElement("div");
  for (line_start = start; line_start < end; line_start += line_length) {
    line_end = Math.min(end, line_start + line_length);
    if (line_start >= line_end) break;
    bunch = document.createElement("div");
    bunch.className = "bunch";
    // line 1: Motif Names
    bunch.appendChild(make_annoseq_names(match, strand, line_start, line_end));
    // line 2: Motif P-values
    bunch.appendChild(make_annoseq_pvalues(match, strand, line_start, line_end));
    // line 3: Motif Best Match
    bunch.appendChild(make_annoseq_motif(match, strand, line_start, line_end));
    // line 6: Sequence
    bunch.appendChild(make_annoseq_sequence(match, strand, line_start, line_end));
    // line 7: Match bounds
    bunch.appendChild(make_annoseq_bounds(match, strand, line_start, line_end));
    container.appendChild(bunch);
  }
  return container;
}

// make some modifications to the data
(function() {
  var i;
  // assign each motif a colour index
  for (i = 0; i < data.motifs.length; i++) data.motifs[i].colour_index = i;
  // check for load
  window.addEventListener("load", make_seq_table_extendable, false);
})();
