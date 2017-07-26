var mast_alphabet = new Alphabet(data.alphabet, data.background.freqs);
var mast_seq_alphabet = (data.sequence_alphabet == null ? mast_alphabet : new Alphabet(data.sequence_alphabet));
var xnum = (data.xlate == null ? 1 : data.xlate);

function redraw_annotated_sequence(tbody, first_display) {
  var strands, annoboxcont, annobox, linelen, seqidx, maxlen, sequence, range;
  strands = parseInt(tbody.querySelector(".strand_selector").value, 10);
  annoboxcont = tbody.querySelector(".annotated_sequence_container");
  linelen = Math.floor((annoboxcont.offsetWidth - 80) / $("ruler").offsetWidth * 5);
  annobox = tbody.querySelector(".annotated_sequence");
  seqidx = parseInt(annobox.getAttribute("data-seqidx"), 10);
  sequence = data.sequences[seqidx];
  if (first_display) {
    var start = 0;
    range = {"start": start, "end": Math.min(start + linelen, sequence.length)};
    set_block_needle_positions(tbody, range.start, range.end);
  } else {
    range = get_block_needle_positions(tbody);
  }
  annobox.innerHTML = "";
  annobox.appendChild(make_annotated_sequence(sequence, strands, range.start, range.end, linelen));
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
 * Make a colourised sequence.
 */
function make_seq(seq, alphabet) {
  var i, j, letter, letter_i, lbox, sbox;
  if (alphabet == null) alphabet = mast_alphabet;
  sbox = document.createElement("span");
  for (i = 0; i < seq.length; i = j) {
    letter = seq.charAt(i);
    letter_i = mast_alphabet.get_index(letter);
    for (j = i+1; j < seq.length; j++) {
      if (mast_alphabet.get_index(seq.charAt(j)) !== letter_i) {
        break;
      }
    }
    if (letter_i != -1)  {
      lbox = document.createElement("span");
      lbox.style.color = mast_alphabet.get_colour(letter_i);
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
    name = (db.name != null ? db.name : db.source);
    seq_sum += db.sequence_count;
    res_sum += db.residue_count;
    add_text_cell(row, name);
    add_text_cell(row, db.sequence_count);
    add_text_cell(row, db.residue_count);
    add_text_cell(row, db.last_mod_date);
    db_stats.appendChild(row);
  }
  row = document.createElement("tr");
  add_text_header_cell(row, "Total");
  add_text_cell(row, seq_sum);
  add_text_cell(row, res_sum);
  add_text_cell(row, "\u00a0");
  db_summary.appendChild(row);
}

function make_motif_databases_table() {
  var i, db_stats, row, db, name;
  db_stats = $("motif_db_stats");
  db_stats.innerHTML = "";
  for (i = 0; i < data.motif_dbs.length; i++) {
    row = document.createElement("tr");
    db = data.motif_dbs[i];
    name = (db.name != null ? db.name : db.source);
    add_text_cell(row, name);
    add_text_cell(row, db.last_mod_date);
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
  if (data.sequence_alphabet == null && rc) pspm.reverse_complement(alphabet);
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
  var i, correlation;
  add_text_header_cell(row, "" + (index + 1) + ".", "", "motif_ordinal");
  add_cell(row, make_preview(mast_alphabet, motif), "motif_logo");
  add_text_cell(row, motif.id, "motif_id");
  add_text_cell(row, (motif.alt != null ? motif.alt : "") , "motif_alt");
  add_text_cell(row, motif.len, "motif_width");
  for (i = 0; i < data.motifs.length; i++) {
    if (i == index) {
      add_text_cell(row, "--", "motif_same");
    } else {
      add_text_cell(row, data.correlation[index][i], "motif_similarity" +
          (data.motifs[i].bad ? " remove" : "") +
          (data.correlation[index][i] > data.settings.max_correlation ? " bad" : ""));
    }
  }
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
  add_text_header_cell(row, "Motif Similarity Matrix", "pop_motif_similarity", "motif_similarity", null, data.motifs.length);

  row = thead.insertRow(thead.rows.length);
  add_text_header_cell(row, "");
  add_text_header_cell(row, "");
  add_text_header_cell(row, "");
  add_text_header_cell(row, "");
  add_text_header_cell(row, "");
  for (i = 1; i <= data.motifs.length; i++) {
    add_text_header_cell(row, "" + i + ".", "", "motif_ordinal" + (data.motifs[i-1].bad ? " remove" : ""));
  }

  for (i = 0; i < data.motifs.length; i++) {
    row = tbody.insertRow(tbody.rows.length);
    if (data.motifs[i].bad) row.className = "bad";
    make_motif_table_entry(row, i, data.motifs[i]);
  }

  return tbl;
}

function make_nos_diagram() {
  var i, row, container;
  if (data.nos == null) {
    $("nos_sec").style.display = "none";
    return;
  }
  var len = 0;
  for (i = 0; i < data.nos.length; i++) {
    if (i > 0) {
      len += data.nos[i].gap;
    }
    len += data.motifs[data.nos[i].idx].len * xnum;
  }
  var edge = Math.max(1, Math.round(len / 10));
  var full_len = len + (edge * 2);
  var table = $("nos_table");
  row = document.createElement("tr");
  container = make_block_container(
      mast_seq_alphabet.has_complement(),
      data.settings.strand_handling != "norc",
      full_len, full_len);
  var pos = edge;
  var motif, gap;
  for (i = 0; i < data.nos.length; i++) {
    if (i > 0) {
      gap = data.nos[i].gap;
      if (gap > 0) {
        make_block_label(container, full_len, pos, gap, "" + gap);
      }
      pos += gap;
    }
    motif = data.motifs[data.nos[i].idx];
    make_block(container, full_len,
        pos, motif.len * xnum, 1e-10, 
        false, motif.colour_index, false);
    pos += motif.len * xnum;
  }
  add_cell(row, container, "block_td");
  var tbody = document.createElement("tbody");
  tbody.appendChild(row);
  table.appendChild(tbody);
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

function make_hit_popup(container, block, sequence, seg, hit) {
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
      var hit_len = motif.len * xnum;
      pspm = new Pspm(motif.pwm);
      if (hit.rc) {
        if (data.sequence_alphabet == null) {
          pspm.reverse_complement(mast_alphabet);
        } else {
          pspm.reverse();
        }
      }
      var logo = new Logo(mast_alphabet, {x_axis_hidden: true, y_axis: false, xlate_nsyms: xnum});
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
      if (fsize > rflank.length && sequence.length > (hit.pos + hit_len + rflank.length)) {
        rflank += string_mult("-", Math.min(sequence.length - (hit.pos + hit_len), fsize) - rflank.length);
      }
      var match = hit.match.replace(/ /g, "\u00A0");
      if (xnum > 1) {
        parts = match.split('');
        parts.push(""); // need to add to last item too
        match = parts.join(string_mult("\u00A0", xnum - 1));
      }
      var canvas = document.createElement('canvas');
      canvas.className = "motif_box";
      canvas.height = 0;
      canvas.width = hit_len * ($("ruler").offsetWidth / 5);
      draw_logo_on_canvas(logo, canvas, false);
      set_tvar(pop, "tvar_logo_pad", lflank);
      set_tvar(pop, "tvar_logo", canvas);
      set_tvar(pop, "tvar_match_pad", lflank);
      set_tvar(pop, "tvar_match", match);
      if (data.sequence_alphabet != null) {
        var xlate = hit.translation;
        if (xnum > 1) {
          parts = xlate.split('');
          parts.push(""); // need to add to last item too
          xlate = parts.join(string_mult(".", xnum - 1));
        }
        set_tvar(pop, "tvar_xlate_pad", lflank);
        set_tvar(pop, "tvar_xlate", make_seq(xlate));
      } else {
        pop.querySelector("div.xlate").style.display = "none";
      }
      set_tvar(pop, "tvar_lflank", lflank);
      set_tvar(pop, "tvar_hit", make_seq(seq));
      set_tvar(pop, "tvar_rflank", rflank);
      set_tvar(pop, "tvar_motif", motif.id);
      set_tvar(pop, "tvar_pvalue", hit.pvalue.toExponential(1));
      set_tvar(pop, "tvar_start", hit.pos + 1);
      set_tvar(pop, "tvar_end", hit.pos + hit_len);
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

function make_seq_diagram(max_seq_len, sequence, range_handler) {
  var seg_i, seg, hit_i, hit, motif, colour_i, hit_len, elems;
  var genomic_coord = null;
  var container = make_block_container(
      mast_seq_alphabet.has_complement(),
      data.settings.strand_handling != "norc",
      max_seq_len, sequence.length, genomic_coord, range_handler
      );
  if (sequence.segs != null) {
    for (seg_i = 0; seg_i < sequence.segs.length; seg_i++) {
      seg = sequence.segs[seg_i];
      for (hit_i = 0; hit_i < seg.hits.length; hit_i++) {
        hit = seg.hits[hit_i];
        motif = data.motifs[hit.idx];
        hit_len = motif.psm.length * xnum;
        block = make_block(container, max_seq_len,
            hit.pos, hit_len, hit.pvalue, hit.rc, motif.colour_index, false);
        popup = make_hit_popup(container, block, sequence, seg, hit);
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
  add_text_header_cell(row, "Sequence", "pop_sequence");
  add_text_header_cell(row, "E-value", "pop_evalue", null, null, (data.settings.strand_handling == "separate" ? 3 : 1));
  add_text_header_cell(row, "");
  add_text_header_cell(row, "Block Diagram", "pop_diagram");
  var seqs_head = document.createElement("thead");
  seqs_head.appendChild(row);
  return seqs_head;
}

function make_seq_footer(max_seq_len) {
  var row = document.createElement("tr");
  add_text_cell(row, "");
  add_text_cell(row, "");
  if (data.settings.strand_handling == "separate") {
    add_text_cell(row, "");
    add_text_cell(row, "");
  }
  add_text_cell(row, "");
  add_cell(row, make_block_ruler(max_seq_len));
  var seqs_foot = document.createElement("tfoot");
  seqs_foot.appendChild(row);
  return seqs_foot;
}

function make_seq_row1(max_seq_len, sequence) {
  var ev1, ev2;
  var row = document.createElement("tr");
  var link = data.sequence_dbs[sequence.db].link;
  // Create name cell, wrapped in link to external DB, wrapped in span.
  var name = null;
  if (link == null) {
    //add_text_cell(row, sequence.name);
    name = document.createTextNode(sequence.name);
  } else {
    //add_cell(row, make_link(sequence.name, link.replace(/SEQUENCEID/, sequence.name)));
    name = make_link(sequence.name, link.replace(/SEQUENCEID/, sequence.name));
  }
  var span = document.createElement("span");
  if (sequence.comment == "") {
    span.title = "No description available.";
  } else {    
    span.title = sequence.comment;
  }
  span.appendChild(name);
  add_cell(row, span);

  // evalue
  if (data.settings.strand_handling != "separate") {
    add_text_cell(row, sequence.score[0].evalue.toExponential(1));
  } else {
    ev1 = (sequence.score[0] != null ? sequence.score[0].evalue : null);
    ev2 = (sequence.score[1] != null ? sequence.score[1].evalue : null);
    add_text_cell(row, (ev1 != null ? ev1.toExponential(1) : "--"), (ev1 == null || (ev2 != null && ev1 > ev2) ? "dim" : ""));
    add_text_cell(row, "/");
    add_text_cell(row, (ev2 != null ? ev2.toExponential(1) : "--"), (ev2 == null || (ev1 != null && ev2 > ev1) ? "dim" : ""));
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
  add_cell(row, make_seq_diagram(max_seq_len, sequence,
    function (start, end, done) {
      if (done) redraw_annotated_sequence(row.parentNode, false);
    }
    ), "block_td");
  return row;
}

function make_seq_row2(maxlen, seqidx, sequence) {
  var row = document.createElement("tr");
  row.className = "more";
  var infobox = clone_template("t_infobox");
  if (data.settings.strand_handling == "separate") infobox.className += " separate";
  infobox.querySelector(".comment").textContent = sequence.comment;
  var pv1_box = infobox.querySelector(".pv_1");
  var pv2_box = infobox.querySelector(".pv_2");
  if (data.settings.strand_handling != "separate") {
    pv1_box.textContent = sequence.score[0].combined_pvalue.toExponential(1);
  } else {
    var pv1 = (sequence.score[0] != null ? sequence.score[0].combined_pvalue : null);
    var pv2 = (sequence.score[1] != null ? sequence.score[1].combined_pvalue : null);
    pv1_box.textContent = (pv1 != null ? pv1.toExponential(1) : "--");
    pv1_box.className = (pv1 == null || (pv2 != null && pv1 > pv2) ? "dim" : "");
    pv2_box.textContent = (pv2 != null ? pv2.toExponential(1) : "--");
    pv2_box.className = (pv2 == null || (pv1 != null && pv2 > pv1) ? "dim" : "");
  }
  var strands = infobox.querySelector(".strand_selector");
  strands.addEventListener("click", function (evt) {
    redraw_annotated_sequence(row.parentNode, false);
  }, false);
  
  var annoseq = infobox.querySelector(".annotated_sequence");
  annoseq.setAttribute("data-maxlen", maxlen);
  annoseq.setAttribute("data-seqidx", seqidx);
  var container = document.createElement("td");
  container.colSpan = (data.settings.strand_handling == "separate" ? 6 : 4);
  container.appendChild(infobox);
  row.appendChild(container);
  return row;
}

function make_seq_table() {
  var i, row_group;
  // calculate the maximum sequence length
  var max_seq_len = 0;
  for (i = 0; i < data.sequences.length; i++) {
    if (data.sequences[i].length > max_seq_len) {
      max_seq_len = data.sequences[i].length;
    }
  }
  // setup the sequence table
  var table = $("seqs_table");
  table.innerHTML = "";
  table.appendChild(make_seq_header());
  table.appendChild(make_seq_footer(max_seq_len));
  // create rows for each sequence
  for (i = 0; i < data.sequences.length; i++) {
    row_group = document.createElement("tbody");
    row_group.appendChild(make_seq_row1(max_seq_len, data.sequences[i]));
    row_group.appendChild(make_seq_row2(max_seq_len, i, data.sequences[i]));
    table.appendChild(row_group);
  }
  make_seq_legend($("seqs_legend_top"));
  make_seq_legend($("seqs_legend_bottom"));

  window.addEventListener("resize", delayed_rewrap, false);
}

function make_annoseq_hit(sequence, strand, line_start, line_end, hit_callback) {
  var line = document.createElement("div");
  line.className = "sequence";
  var i, j, seg, hit, motif, hit_len, hit_end_pos, box, label;
  var pos = line_start;
  for (i = 0; i < sequence.segs.length; i++) {
    seg = sequence.segs[i];
    if (seg.pos >= line_end) break;
    if ((seg.pos + seg.data.length) <= line_start) continue;
    for (j = 0; j < seg.hits.length; j++) {
      hit = seg.hits[j];
      if (hit.pos >= line_end) break;
      if (hit.rc ? strand > 0 : strand < 0) continue;
      motif = data.motifs[hit.idx];
      hit_len = motif.len * xnum;
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

function make_annoseq_label(sequence, strand, line_start, line_end, label_callback) {
  return make_annoseq_hit(sequence, strand, line_start, line_end, function(hit, motif, pos_start, pos_end) {
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

function make_annoseq_names(sequence, strand, line_start, line_end) {
  return make_annoseq_label(sequence, strand, line_start, line_end, function(hit, motif) {
    var show_strand = data.settings.strand_handling != "unstranded";
    return (show_strand ? "(" + (hit.rc ? "-" : "+") + ") " : "") + 
        (/^\d+$/.test(motif.id) ? "Motif " : "") + motif.id;
  });
}

function make_annoseq_pvalues(sequence, strand, line_start, line_end) {
  return make_annoseq_label(sequence, strand, line_start, line_end, function(hit, motif) {
    return hit.pvalue.toExponential(1);
  });
}

function make_annoseq_motif(sequence, strand, line_start, line_end) {
  var line =  make_annoseq_hit(sequence, strand, line_start, line_end, function(hit, motif, pos_start, pos_end) {
    var pspm = new Pspm(motif.pwm);
    if (hit.rc) {
      if (data.sequence_alphabet == null) {
        pspm.reverse_complement(mast_alphabet);
      } else {
        pspm.reverse();
      }
    }
    var logo = new Logo(mast_alphabet, {x_axis_hidden: true, y_axis: false,
      xlate_nsyms: xnum, xlate_start: pos_start - hit.pos, xlate_end: pos_end - hit.pos});
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

function make_annoseq_match(sequence, strand, line_start, line_end) {
  return make_annoseq_hit(sequence, strand, line_start, line_end, function(hit, motif, pos_start, pos_end) {
    var match = hit.match.replace(/ /g, "\u00A0");
    if (xnum > 1) {
      parts = match.split('');
      parts.push(""); // need to add to last item too
      match = parts.join(string_mult("\u00A0", xnum - 1));
    }
    return document.createTextNode(match.substring(pos_start - hit.pos, pos_end - hit.pos));
  });
}

function make_annoseq_translation(sequence, strand, line_start, line_end) {
  return make_annoseq_hit(sequence, strand, line_start, line_end, function(hit, motif, pos_start, pos_end) {
    var translation = hit.translation;
    if (xnum > 1) {
      var parts = translation.split('');
      parts.push("");
      translation = parts.join(string_mult(".", xnum - 1));
    }
    return make_seq(translation.substring(pos_start - hit.pos, pos_end - hit.pos), mast_alphabet);
  });
}

function make_annoseq_sequence(sequence, strand, line_start, line_end) {
  var line = document.createElement("div");
  line.className = "sequence sequence_line";
  var i, j, seg, hit, motif, hit_len, joiner, grey_seq, seg_end_pos, hit_end_pos;
  var start_pos = document.createElement("div");
  start_pos.className = "line_start";
  start_pos.appendChild(document.createTextNode(line_start + 1));
  line.appendChild(start_pos);
  var pos = line_start;
  for (i = 0; i < sequence.segs.length; i++) {
    seg = sequence.segs[i];
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
      hit_len = motif.len * xnum;
      if ((hit.pos + hit_len) <= line_start) continue;
      if (pos < hit.pos) {
        grey_seq = document.createElement("span");
        grey_seq.className = "joiner";
        grey_seq.appendChild(document.createTextNode(seg.data.substring(pos - seg.pos,  hit.pos - seg.pos)));
        line.appendChild(grey_seq);
        pos = hit.pos;
      }
      hit_end_pos = Math.min(hit.pos + hit_len, line_end);
      line.appendChild(make_seq(seg.data.substring(pos - seg.pos, hit_end_pos - seg.pos), mast_seq_alphabet));
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
  end_pos.appendChild(document.createTextNode(line_end));
  line.appendChild(end_pos);
  return line;
}

function make_annoseq_bounds(sequence, strand, line_start, line_end) {
  return make_annoseq_hit(sequence, strand, line_start, line_end, function(hit, motif, pos_start, pos_end) {
    var hit_str = "\\" + string_mult("_", (motif.len * xnum) - 2) + "/";
    return document.createTextNode(hit_str.substring(pos_start - hit.pos, pos_end - hit.pos));
  });
}



function make_annotated_sequence(sequence, strand, start, end, line_length) {
  var line_start, bunch, line_end;
  var container = document.createElement("div");
  for (line_start = start; line_start < end; line_start += line_length) {
    line_end = Math.min(end, line_start + line_length);
    if (line_start >= line_end) break;
    bunch = document.createElement("div");
    bunch.className = "bunch";
    // line 1: Motif Names
    bunch.appendChild(make_annoseq_names(sequence, strand, line_start, line_end));
    // line 2: Motif P-values
    bunch.appendChild(make_annoseq_pvalues(sequence, strand, line_start, line_end));
    // line 3: Motif Best Match
    bunch.appendChild(make_annoseq_motif(sequence, strand, line_start, line_end));
    // line 4: Match Strength
    bunch.appendChild(make_annoseq_match(sequence, strand, line_start, line_end));
    if (data.sequence_alphabet != null) {
      // line 5: Match Sequence Translation (Optional)
      bunch.appendChild(make_annoseq_translation(sequence, strand, line_start, line_end));
    }
    // line 6: Sequence
    bunch.appendChild(make_annoseq_sequence(sequence, strand, line_start, line_end));
    // line 7: Match bounds
    bunch.appendChild(make_annoseq_bounds(sequence, strand, line_start, line_end));
    container.appendChild(bunch);
  }
  return container;
}

function make_other_settings() {
  $("opt_strand_handling").innerHTML = (function() {
    var sh = data.settings.strand_handling;
    if (sh === "combine") {
      return "The result of scanning both strands is <b>combined</b>. When matches overlap the non-overlapping combination with the best p-value is shown.";
    } else if (sh === "separate") {
      return "The two <b>separate</b> strands are scanned as if they were different sequences and only recombined to display.";
    } else if (sh === "norc") { 
      return "Only the given strand is scanned.";
    } else if (sh === "unstranded") {
      return "The alphabet is unstranded.";
    }
  })();
  $("opt_max_correlation").innerHTML = (function() {
    var max_corr = data.settings.max_correlation;
    if (max_corr < 1.0) {
      return "Motifs with a correlation greater than <b>" + max_corr + "</b> are marked for potential removal dependant on the <code>--remcorr</code> option.";
    } else {
      return "Motifs never reach the maximum correlation because it is set to 1 or larger.";
    }
  })();
  $("opt_remove_correlated").innerHTML = (function() {
    if (data.settings.remove_correlated) {
      return "Correlated motifs exceeding the threshold are <b>removed</b>.";
    } else {
      return "Correlated motifs exceeding the threshold are <b>highlighted</b> and their removal is recommended.";
    }
  })();
  $("opt_seq_evalue").innerHTML = "Sequences with an <i>E</i>-value less than <b>" +
    data.settings.max_seq_evalue + "</b> are included in the output.";
  $("opt_adjust_hit").innerHTML = (function() {
    if (data.settings.adjust_hit) {
      return "The hit <i>p</i>-value is <b>adjusted</b> for the length of the sequence so the same hit in a longer sequence will be considered less important.";
    } else {
      return "The hit <i>p</i>-value is <b>not adjusted</b> for the length of the sequence.";
    }
  })();
  $("opt_max_hit_pvalue").innerHTML = "The <i>p</i>-value of a hit must be less than <b>" + data.settings.max_hit_pvalue + "</b> to be shown in the output.";
  $("opt_weak_hit_pvalue").innerHTML = (function() {
    if (data.settings.max_weak_pvalue > data.settings.max_hit_pvalue) {
      return "The <i>p</i>-value of a weak hit msut be less than <b>" + data.settings.max_weak_pvalue + "</b> to be shown in the output.";
    } else {
      return "Weak hits are <b>not displayed</b>.";
    }
  })();
}

function psm2pwm(psm, bg) {
  "use strict";
  var i, j, score, freq, freqs, min, max, sum;
  var pwm = [];
  for (i = 0; i < psm.length; i++) {
    freqs = [];
    min = Number.POSITIVE_INFINITY;
    for (j = 0; j < psm[i].length; j++) {
      score = psm[i][j];
      freq = Math.pow(2, score / 100) * bg[j];
      if (freq < min) min = freq;
      freqs[j] = freq;
    }
    // shift to make positive
    if (min < 0) {
      for (j = 0; j < freqs.length; j++) freqs[j] += -min;
    }
    // scale so they sum to 1
    sum = 0;
    for (j = 0; j < freqs.length; j++) sum += freqs[j];
    for (j = 0; j < freqs.length; j++) freqs[j] /= sum;
    pwm[i] = freqs;
  }
  return pwm;
}

// make some modifications to the data
(function() {
  var motif;
  var bad_motifs = 0;
  var colour_index = 0;
  for (i = 0; i < data.motifs.length; i++) {
    motif = data.motifs[i];
    motif.pwm = psm2pwm(motif.psm, data.motif_dbs[motif.db].bg);
    if (typeof motif.bad !== "boolean") motif.bad = false;
    if (motif.bad) bad_motifs++;
    if (!data.settings.remove_correlated || !motif.bad) motif.colour_index = colour_index++;
    else motif.colour_index = null;
  }
  data.bad_motif_count = bad_motifs; 
})();
