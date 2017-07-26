"use strict";

var MAX_INT = 9007199254740992;
var alphabet = new Alphabet(data.alphabet, data.background);
var orient_list = [[0],[1],[2],[3],[0, 2], [0, 3], [2, 1], [1, 3], [0, 2, 1, 3]];
var orient_names = ["up+", "dn+", "up-", "dn-", 
    "up+/up-", "up+/dn-", "up-/dn+", "dn+/dn-", "all"];
var orient_desc = ["upstream / same strand", "downstream / same strand",
    "upstream / opposite strand", "downstream / opposite strand", 
    "upstream / secondary palindromic", "upstream / primary palindromic", 
    "downstream / primary palindromic", "downstream / secondary palindromic", 
    "all / both palindromic"];
if (!alphabet.has_complement()) {
  orient_names = ["up", "dn"];
  orient_desc = ["upstream", "downstream"];
}

var DelayLogoTask = function(logo, canvas) {
  "use strict";
  canvas.width = canvas.width; // clear canvas
  this.logo = logo;
  this.canvas = canvas;
};

DelayLogoTask.prototype.run = function () {
  "use strict";
  this.canvas.width = this.canvas.width; // clear canvas
  draw_logo_on_canvas(this.logo, this.canvas);
};

function make_logo(canvas, motif) {
  var pspm;
  pspm = new Pspm(motif.pwm, motif.id, motif.ltrim, motif.rtrim, motif.nsites);
  var logo = new Logo(alphabet);
  logo.add_pspm(pspm, 0);
  canvas.width = 0; // allow any width
  size_logo_on_canvas(logo, canvas);
  add_draw_task(canvas, new DelayLogoTask(logo, canvas));
  return canvas;
}

function make_pwm_logo(canvas, motif_link, eps_link, pwm, dots, name, shift) {
  if (pwm == null) {
    canvas.style.visibility = "hidden";
    if (motif_link) motif_link.style.visibility = "hidden";
    if (eps_link) eps_link.style.visibility = "hidden";
    return 0;
  }
  var pspm;
  pspm = new Pspm(pwm, name);
  var logo = new Logo(alphabet, {x_axis: !dots, y_axis: !dots});
  logo.add_pspm(pspm, shift);
  if (motif_link) {
    prepare_download(pspm.as_meme({"alphabet": alphabet, "with_header": true}), "text/plain", (typeof name == "string" ? name.replace(/ /g, "_") + ".meme" : "motif.meme"), motif_link);
    motif_link.style.visibility = "visible";
  }
  canvas.width = 0; // clear canvas
  size_logo_on_canvas(logo, canvas);
  if (eps_link) {
    prepare_download(logo.as_eps(), "application/postscript", (typeof name == "string" ? name.replace(/ /g, "_") + ".eps" : "logo.eps"), eps_link);
    eps_link.style.visibility = "visible"
  }
  canvas.style.visibility = "visible";
  add_draw_task(canvas, new DelayLogoTask(logo, canvas));
  return canvas.width;
}

function make_minimal_spacing_diagram(spacing, start, length) {
  "use strict";
  var i, j, canvas, ctx, entry, colours;
  var w, h, g;
  colours = ["#DDD", "#CCC", "#BBB", "#AAA", "#999", "#888", "#777", "#666",
          "#555", "#444", "#333", "#222", "#111", "#000"];
  w = 1;
  h = 2;
  g = 0;
  canvas = document.createElement("canvas");
  canvas.height = 9 * (h + g);
  canvas.width = length * (w+g);
  ctx = canvas.getContext('2d');
  ctx.fillStyle = "#EEE";
  ctx.fillRect(0, 0, canvas.width, canvas.height); 
  ctx.fillStyle = "#000"
  for (i = 0; i < spacing.length; i++) {
    entry = spacing[i];
    if (entry.bin < start || entry.bin >= (start + length)) continue;
    ctx.fillStyle = (i != 0 ? 
        colours[
        Math.min(
          Math.max(
            0, 
            Math.round(-(Math.log(entry.pvalue)/Math.LN10))),
          colours.length - 1)
        ] :
        "red");
    ctx.fillRect((entry.bin - start) * (w + g), entry.orient * (h + g), w, h);
  }
  return canvas;
}

function make_compact_graph(secondary) {
  var i, j;
  var h = 30;
  var w = data.options.margin;
  var bsp = secondary.spacings[0]; // best spacing
  var counts = [];
  for (i = 0; i < w; i++) counts[i] = 0;
  var orients = orient_list[bsp.orient];
  for (i = 0; i < orients.length; i++) {
    for (j = 0; j < secondary.counts[orients[i]].length; j++) {
      counts[j] += secondary.counts[orients[i]][j];
    }
  }
  var count_max = Math.max.apply(null, counts);
  var spacing_bins = {};
  for (i = 0; i < secondary.spacings.length; i++) {
    if (secondary.spacings[i].orient != bsp.orient) continue;
    spacing_bins[secondary.spacings[i].bin] = true;
  }
  var canvas = document.createElement("canvas");
  canvas.width = w;
  canvas.height = h;
  var ctx = canvas.getContext('2d');
  ctx.fillStyle = "white";
  ctx.fillRect(0, 0, canvas.width, canvas.height); 
  ctx.fillStyle = "black";
  for (i = 0; i < w; i++) {
    if (spacing_bins[i]) continue;
    ctx.fillRect(i, h, 1, -(h * (counts[i] / count_max)));
  }
  ctx.fillStyle = "red";
  for (i = 1; i < secondary.spacings.length; i++) {
    var spacing = secondary.spacings[i];
    if (spacing.orient != bsp.orient) continue;
    ctx.fillRect(spacing.bin, h, 1, -(h * (counts[spacing.bin] / count_max)));
  }
  ctx.fillRect(bsp.bin, h, 1, -(h * (counts[bsp.bin] / count_max)));
  return canvas;
}

function make_alignment_eps(manager) {
  "use strict";
  var i, alignments, logo, pspm;
  alignments = manager.spacing.alignment_pwm;
  logo = new Logo(alphabet, {x_axis: false, y_axis: false});
  for (i = 0; i < alignments.length; i++) {
    if (alignments[i] != null) {
      pspm = new Pspm(alignments[i], orient_desc[i]);
      logo.add_pspm(pspm, 0);
    }
  }
  var seqs_desc = "alignseqs_" + manager.primary.motif.id + "_with_" +
    manager.secondary.motif.id + "_g" + manager.spacing.bin + "_o" +
    manager.spacing.orient;
  prepare_download(logo.as_eps(), "application/postscript", seqs_desc,
      manager.block.querySelector(".sa_align_eps"));
}

function make_secondary_pair_eps(manager) {
  "use strict";
  var logo, motif, pspm, eps_content;

  logo = new Logo(alphabet);
  motif = manager.secondary.motif;
  pspm = new Pspm(motif.pwm, motif.id, motif.ltrim, motif.rtrim, motif.nsites);
  logo.add_pspm(pspm, 0);
  pspm = new Pspm(manager.spacing.inferred_pwm);
  logo.add_pspm(pspm, 0);
  eps_content = logo.as_eps({"title": motif.id, "xaxislabel": "Inferred motif"});
  prepare_download(eps_content, "application/postscript", "inferred.eps",
      manager.block.querySelector(".sa_slogo_eps"));
}

function draw_overview_graph(manager) {
  "use strict";
  var canvas, hl, spacings, spacing, orients, i, j, counts, max_count, spamo_graph;
  canvas = manager.block.querySelector(".sa_overview_graph")
  hl = [{},{},{},{}];
  // set pink highlights for all bins involved in a significant spacing
  if (manager.highlight_all) {
    spacings = manager.secondary.spacings;
    for (i = 0; i < spacings.length; i++) {
      spacing = spacings[i];
      orients = orient_list[spacing.orient];
      for (j = 0; j < orients.length; j++) {
        hl[orients[j]][spacing.bin] = "pink";
      }
    }
  }
  // set red highlight for all bins involved in the selected spacing
  if (manager.highlight_selected) {
    orients = orient_list[manager.spacing.orient];
    for (j = 0; j < orients.length; j++) {
      hl[orients[j]][manager.spacing.bin] = "red";
    }
  }

  counts = manager.secondary.counts;
  max_count = Math.max.apply(null, counts.map( function (list) { return Math.max.apply(null, list);} ));
  spamo_graph = new SpamoQuadGraph(data.options.margin, data.options.bin_size, max_count, manager.secondary.motif.len, counts, hl);
  canvas.width = canvas.width;
  var ctx = canvas.getContext('2d');
  var eps_ctx = new EpsContext(ctx, canvas.width, canvas.height);
  eps_ctx.register_font("9px Helvetica", "Helvetica", 9);
  eps_ctx.register_font("12px Helvetica", "Helvetica", 12);
  eps_ctx.register_font("14px Helvetica", "Helvetica", 14);
  eps_ctx.register_font("18px Helvetica", "Helvetica", 18);
  spamo_graph.draw_graph(ctx, canvas.width, canvas.height);
  var filename = "overview_" + manager.primary.motif.id + "_to_" +
      manager.secondary.motif.id + "_g" + manager.spacing.bin + "_o" +
      manager.spacing.orient + ".eps";
  prepare_download(eps_ctx.eps(), "application/postscript", filename,
      manager.block.querySelector(".sa_overview_graph_dl"));
}

function best_window(orient, range, req_bin, max_bin, spacings) {
  var i;
  var left = Math.max(req_bin - range + 1, 0);
  var right = Math.min(req_bin + range - 1, max_bin);
  var size = right - left + 1;
  var starts = [];
  for (i = 0; i < size; i++) starts[i] = 0;
  var peak_i = 0;
  var peak_pvalue = 1;
  for (i = 0; i < spacings.length; i++) {
    var spacing = spacings[i];
    if (spacing.orient == orient && spacing.bin >= left && spacing.bin <= right) {
      starts[spacing.bin - left] = Math.log(spacing.pvalue);
      if (spacing.pvalue < peak_pvalue) {
        peak_pvalue = spacing.pvalue;
        peak_i = spacing.bin - left;
      }
    }
  }
  for (i = 1; i < size; i++) starts[i] += starts[i - 1];
  var best_log_pvalues = 0;
  var best_peak_dist = range;
  var best_i = [];
  for (i = 0; i < (size - range + 1); i++) {
    var log_pvalues = starts[i + range - 1] - starts[i];
    var peak_dist = Math.abs((range / 2) - peak_i + i);
    if (log_pvalues < best_log_pvalues || 
        (log_pvalues == best_log_pvalues && peak_dist < best_peak_dist)) {
      best_log_pvalues = log_pvalues;
      best_peak_dist = peak_dist;
      best_i = [i];
    } else if (log_pvalues == best_log_pvalues && peak_dist == best_peak_dist) {
      best_i.push(i);
    }
  }
  return (best_i.length > 0 ? left + best_i[0] : req_bin - Math.round(range / 2));
}

function draw_sg(manager) {
  "use strict";
  var canvas, ctx;
  var i, hl;
  var min_count, max_count, metrics;
  canvas = manager.block.querySelector(".sa_selected_graph");
  hl = {};
  if (manager.highlight_all) {
    for (i = 0; i < manager.secondary.spacings.length; i++) {
      if (manager.secondary.spacings[i].orient == manager.spacing.orient) {
        hl[manager.secondary.spacings[i].bin] = "pink";
      }
    }
  }
  if (manager.highlight_selected) {
    hl[manager.spacing.bin] = "red";
  }
  var range = 20;
  var max_bin = data.options.margin - manager.secondary.motif.len;
  var start = best_window(manager.spacing.orient, range, manager.spacing.bin, max_bin, manager.secondary.spacings);
  var bin_size = data.options.bin_size;
  var graph = new SpamoOrientGraph(manager.secondary.counts, orient_list[manager.spacing.orient], start, range, bin_size, hl);
  canvas.width = canvas.width; // clear the canvas
  ctx = canvas.getContext('2d');
  var eps_ctx = new EpsContext(ctx, canvas.width, canvas.height);
  eps_ctx.register_font("9px Helvetica", "Helvetica", 9);
  eps_ctx.register_font("12px Helvetica", "Helvetica", 12);
  eps_ctx.register_font("14px Helvetica", "Helvetica", 14);
  eps_ctx.register_font("18px Helvetica", "Helvetica", 18);
  graph.draw_graph(ctx, canvas.width, canvas.height);
  var filename = "selorient_" + manager.primary.motif.id + "_to_" + manager.secondary.motif.id + "_g" + manager.spacing.bin + "_o" + manager.spacing.orient + ".eps"
  prepare_download(eps_ctx.eps(), "application/postscript", filename,
      manager.block.querySelector(".sa_selected_graph_dl"));
}

function make_name_txt(motif, nbsp) {
  var name;
  name = motif.id + (typeof motif.alt === "string" ? (nbsp ? "\u00A0" : " ") + "(" + motif.alt + ")" : "");
  return name;
}

function make_name(motif, nbsp, should_link) {
  var name, link;
  name = motif.id + (typeof motif.alt === "string" ? (nbsp ? "\u00A0" : " ") + "(" + motif.alt + ")" : "");
  if (should_link && typeof motif.url === "string") {
    link = document.createElement("a");
    link.href = motif.url;
    link.appendChild(document.createTextNode(name));
    return link;
  } else {
    return document.createTextNode(name);
  }
}

function make_pri_table() {
  "use strict";
  function _make_list_handler(primaryi, primary, secondary) {
    return function(e) {
      var block, manager;
      // get the block
      block = document.getElementById("primary_" + primaryi);
      if (!block) return;
      // get the manager
      manager = block.data_manager;
      manager.primary = primary;
      manager.secondary = secondary;
      manager.spacing = secondary.spacings[0];
      update_secondaries_list(manager);
      update_selected_secondary(manager);
      block.scrollIntoView();
      e.preventDefault();
    };
  }
  var pri_table, pri_row, pri_tbody;
  var primary, secondaries, secondary;
  var row, row_mlist, mlink;
  var i, j;

  pri_table = $("pri_tbl");
  pri_row = pri_table.querySelector(".pri_row");
  pri_tbody = pri_row.parentNode;
  pri_tbody.removeChild(pri_row);
  for (i = 0; i < data.primaries.length; i++) {
    primary = data.primaries[i].motif;
    secondaries = data.primaries[i].secondaries;
    row = pri_row.cloneNode(true);
    row.querySelector(".pri_db").appendChild(document.createTextNode(data.primary_dbs[primary.db].name));
    row.querySelector(".pri_name").appendChild(make_name(primary, false, true));
    make_logo(row.querySelector(".pri_logo"), primary);
    row.querySelector(".pri_nlist").appendChild(document.createTextNode(secondaries.length));
    row_mlist = row.querySelector(".pri_list");
    for (j = 0; j < secondaries.length; j++) {
      if (j > 0) row_mlist.appendChild(document.createTextNode(",\u2003 "));// EM Space
      secondary = secondaries[j];
      mlink = document.createElement("a");
      mlink.href = "?pdb=" + primary.db + "&pid=" + primary.id + "&sdb=" + secondary.motif.db + "&sid=" + secondary.motif.id;
      mlink.className = (i % 2 == 0 ? "ml1" : "ml2");
      mlink.appendChild(make_name(secondary.motif, true, false));
      mlink.addEventListener("click", _make_list_handler(i, data.primaries[i], secondary), false);
      //mlink.href = "#match_" + i + "_" + secondary.idx;
      row_mlist.appendChild(mlink);
    }
    pri_tbody.appendChild(row);
  }
}

function make_seq_table() {
  "use strict";
  var seq_table, seq_row, seq_tbody;
  var seq_db, used;
  var row;
  var i;

  seq_table = $("seq_tbl");
  seq_row = seq_table.querySelector(".seq_row");
  seq_tbody = seq_row.parentNode;
  seq_tbody.removeChild(seq_row);
  for (i = 0; i < data.sequence_dbs.length; i++) {
    seq_db = data.sequence_dbs[i];
    row = seq_row.cloneNode(true);
    row.querySelector(".seq_name").appendChild(document.createTextNode(seq_db.name));
    row.querySelector(".seq_last_modified").appendChild(document.createTextNode(seq_db.last_modified));
    row.querySelector(".seq_loaded").appendChild(document.createTextNode(seq_db.loaded));
    row.querySelector(".seq_too_short").appendChild(document.createTextNode(seq_db.excluded_too_short));
    row.querySelector(".seq_ambiguous").appendChild(document.createTextNode(seq_db.excluded_ambigs));
    row.querySelector(".seq_no_primary").appendChild(document.createTextNode(seq_db.excluded_no_match));
    row.querySelector(".seq_too_similar").appendChild(document.createTextNode(seq_db.excluded_similar));
    used = seq_db.loaded - seq_db.excluded_too_short - seq_db.excluded_ambigs - seq_db.excluded_no_match - seq_db.excluded_similar;
    row.querySelector(".seq_used").appendChild(document.createTextNode(used));

    seq_tbody.appendChild(row);
  }
}

function make_sdb_table() {
  "use strict";
  var sdb_table, sdb_row, sdb_tbody;
  var sdb_db, smotif, pm_count, rm_count;
  var row;
  var i, j;

  // count how many motifs
  pm_count = [];
  rm_count = [];
  for (i = 0; i < data.sequence_dbs.length; i++) {
    pm_count[i] = 0;
    rm_count[i] = 0;
  }
  for (i = 0; i < data.secondary_motifs.length; i++) {
    smotif = data.secondary_motifs[i];
    if (smotif.nonredundant) {
      pm_count[smotif.db]++;
    } else {
      rm_count[smotif.db]++;
    }
  }

  sdb_table = $("sdb_tbl");
  sdb_row = sdb_table.querySelector(".sdb_row");
  sdb_tbody = sdb_row.parentNode;
  sdb_tbody.removeChild(sdb_row);
  for (i = 0; i < data.secondary_dbs.length; i++) {
    sdb_db = data.secondary_dbs[i];
    row = sdb_row.cloneNode(true);
    row.querySelector(".sdb_name").appendChild(document.createTextNode(sdb_db.name));
    row.querySelector(".sdb_last_modified").appendChild(document.createTextNode(sdb_db.last_modified));
    row.querySelector(".sdb_loaded").appendChild(document.createTextNode(sdb_db.loaded - sdb_db.excluded));
    row.querySelector(".sdb_sig_motifs").appendChild(document.createTextNode(pm_count[i]));
    row.querySelector(".sdb_redundant_motifs").appendChild(document.createTextNode(rm_count[i]));
    sdb_tbody.appendChild(row);
  }
}

function update_sequence_list(manager) {
  "use strict";
  var csc_span, textarea, format, seqs, mimetype, filename, desc, seqs_desc;
  csc_span = manager.block.querySelector(".sa_contr_seq_count");
  textarea = manager.block.querySelector(".sa_contr_seqs");
  format = manager.block.querySelector(".sa_contr_seqs_format").value;
  seqs = manager.spacing.seqs;
  seqs = seqs.map( function (seqidx) { return data.sequence_names[seqidx]; } );
  mimetype = "text/plain";
  desc =  "_g" + manager.spacing.bin + "_o" + manager.spacing.orient;
  seqs_desc = manager.primary.motif.id + "_with_" + manager.secondary.motif.id + desc;
  filename = "seqs_" + seqs_desc + ".txt";
  if (format == 1) {
    seqs = seqs.map( function (seqid) { return seqid.replace(/^([^:]+):(\d+)-(\d+)$/, "$1\t$2\t$3"); } );
    mimetype = "text/x-bed";
    filename = "seqs_" + seqs_desc + ".bed";
  }
  csc_span.textContent = seqs.length;
  textarea.value = seqs.join("\n");
  prepare_download(textarea.value, mimetype, filename,
      manager.block.querySelector(".sa_contr_seqs_dl"));
}

//FIXME_TLB
function isArray(x) {
    return x.constructor.toString().indexOf("Array") > -1;
}

function update_selected_spacing(manager) {
  var slogo_out, desc, seqs_desc, lwidth;
  slogo_out = manager.block.querySelector(".sa_slogo_out");
  slogo_out.width = 0;
  desc =  "_gap_" + manager.spacing.bin + "_orientation_" + manager.spacing.orient;
  seqs_desc = manager.primary.motif.id + "_with_" + manager.secondary.motif.id + desc;
  var shift = manager.secondary.motif.ltrim;
  make_pwm_logo(slogo_out, manager.block.querySelector(".sa_slogo_out_meme"),
      manager.block.querySelector(".sa_slogo_out_eps"),
      manager.spacing.inferred_pwm, false,
      manager.secondary.motif.id + "_near_" + manager.primary.motif.id + desc, shift);
  lwidth = 0;
  lwidth = make_pwm_logo(manager.block.querySelector(".sa_align_up_pos"),
      manager.block.querySelector(".sa_align_up_pos_meme"),
      manager.block.querySelector(".sa_align_up_pos_eps"),
      manager.spacing.alignment_pwm[SpamoQuadGraph.SAME_LEFT], true,
      "alignseqs_up+_" + seqs_desc, 0);
  lwidth = Math.max(make_pwm_logo(manager.block.querySelector(".sa_align_up_neg"),
      manager.block.querySelector(".sa_align_up_neg_meme"),
      manager.block.querySelector(".sa_align_up_neg_eps"),
      manager.spacing.alignment_pwm[SpamoQuadGraph.OPPOSITE_LEFT], true,
      "alignseqs_up-_" + seqs_desc, 0), lwidth);
  lwidth = Math.max(make_pwm_logo(manager.block.querySelector(".sa_align_dn_pos"),
      manager.block.querySelector(".sa_align_dn_pos_meme"),
      manager.block.querySelector(".sa_align_dn_pos_eps"),
      manager.spacing.alignment_pwm[SpamoQuadGraph.SAME_RIGHT], true,
      "alignseqs_dn+_" + seqs_desc, 0), lwidth);
  lwidth = Math.max(make_pwm_logo(manager.block.querySelector(".sa_align_dn_neg"),
      manager.block.querySelector(".sa_align_dn_neg_meme"),
      manager.block.querySelector(".sa_align_dn_neg_eps"),
      manager.spacing.alignment_pwm[SpamoQuadGraph.OPPOSITE_RIGHT], true,
      "alignseqs_dn-_" + seqs_desc, 0), lwidth);
  manager.block.querySelector(".sa_align_scroll").style.width = lwidth + "px";
  var scroller = manager.block.querySelector(".sa_align_scrollbox");
  var scroll = Math.max(scroller.scrollWidth - scroller.offsetWidth, 0) / 2;
  scroller.scrollLeft = scroll;
  manager.block.querySelector(".sa_scroll_up_pos").scrollLeft = scroll;
  manager.block.querySelector(".sa_scroll_dn_pos").scrollLeft = scroll;
  manager.block.querySelector(".sa_scroll_up_neg").scrollLeft = scroll;
  manager.block.querySelector(".sa_scroll_dn_neg").scrollLeft = scroll;
  make_alignment_eps(manager);
  make_secondary_pair_eps(manager);
  update_sequence_list(manager);
  draw_overview_graph(manager);
  draw_sg(manager);
}

function make_spacing_handler(manager, spacing, row, table) {
  if (manager.spacing == spacing) {
    manager.spacing_row = row;
    toggle_class(row, "active", true);
    substitute_classes(table, ["orient_0", "orient_1", "orient_2", "orient_3",
        "orient_4", "orient_5", "orient_6", "orient_7", "orient_8"],
        ["orient_" + spacing.orient]);
    update_selected_spacing(manager);
  }
  return function(e) {
    // check if already selected
    if (manager.spacing === spacing) return; // ignore multiple clicks
    // update selection
    if (typeof manager.spacing_row !== "undefined") toggle_class(manager.spacing_row, "active", false);
    toggle_class(row, "active", true);
    substitute_classes(table, ["orient_0", "orient_1", "orient_2", "orient_3",
        "orient_4", "orient_5", "orient_6", "orient_7", "orient_8"],
        ["orient_" + spacing.orient]);
    manager.spacing = spacing;
    manager.spacing_row = row;
    update_selected_spacing(manager);
  };
}

function update_selected_secondary(manager) {
  function _set_text(node, text, defval) {
    if (typeof defval === "undefined") defval = "";
    node.innerHTML = "";
    node.appendChild(document.createTextNode((typeof text !== "undefined" ? text : defval)));
  }
  var name, table, tbody, tr, i;
  var spacings = manager.secondary.spacings;
  var best_sp = (spacings.length > 0 ? spacings[0] : null);
  // fill in details of the secondary
  name = manager.block.querySelector(".sa_name");
  name.innerHTML = "";
  name.appendChild(make_name(manager.secondary.motif));
  //_set_text(manager.block.querySelector(".sa_cluster"), manager.secondary.cluster.motif.id);
  _set_text(manager.block.querySelector(".sa_cluster"), make_name_txt(manager.secondary.cluster.motif));
  _set_text(manager.block.querySelector(".sa_evalue"), (best_sp.pvalue * manager.nmotifs).toExponential(2));
  _set_text(manager.block.querySelector(".sa_gap"), best_sp.bin);
  _set_text(manager.block.querySelector(".sa_orient"), orient_desc[best_sp.orient]);
  // fill in the secondary logo
  make_logo(manager.block.querySelector(".sa_slogo_in"), manager.secondary.motif);
  // get the spacings table
  table = manager.block.querySelector(".sa_spacings");
  // clear the table
  for (i = table.tBodies.length - 1; i >= 0; i--) {
    table.removeChild(table.tBodies[i]);
  }
  // now add all spacings
  tbody = document.createElement("tbody");
  for (i = 0; i < spacings.length; i++) {
    var spacing = spacings[i];
    tr = document.createElement("tr");
    tr.className = "sa_sp_row orient_" + spacing.orient;
    add_text_cell(tr, spacing.bin, "sa_sp_gap"); // spacing gap
    add_text_cell(tr, orient_names[spacing.orient], "sa_sp_orient"); // spacing orientation
    add_text_cell(tr, spacing.pvalue.toExponential(2), "sa_sp_pvalue"); // spacing p-value
    // check for clicks
    tr.addEventListener("click", make_spacing_handler(manager, spacing, tr, table), false);
    tbody.appendChild(tr);
  }
  table.appendChild(tbody);
  
}

function make_secondary_handler(manager, secondary, row) {
  if (manager.secondary === secondary) {
    manager.secondary_row = row;
    toggle_class(row, "active", true);
  }
  return function(e) {
    // check if already selected
    if (manager.secondary === secondary) return; // ignore multiple clicks
    // update selection
    if (typeof manager.secondary_row !== "undefined") toggle_class(manager.secondary_row, "active", false);
    if (typeof manager.spacing_row !== "undefined") toggle_class(manager.spacing_row, "active", false);
    toggle_class(row, "active", true);
    manager.secondary = secondary;
    manager.secondary_row = row;
    manager.spacing = secondary.spacings[0];
    delete manager.spacing_row;
    update_selected_secondary(manager);
  };
}

function make_hover_handler(secondary) {
  if (typeof make_hover_handler.timer === "undefined") make_hover_handler.timer = null;
  return function(e) {
    move_logo(e);
    this.addEventListener('mousemove', move_logo, false);
    if (make_hover_handler.motif_idx === secondary.idx) { 
      $("logo_popup").style.display = "block";
    } else {
      if (make_hover_handler.timer) clearTimeout(make_hover_handler.timer);
      make_hover_handler.timer = setTimeout(function() {
        var motif, pspm, logo, canvas;
        // create the motif pspm
        motif = secondary.motif;
        pspm = new Pspm(motif.pwm, motif.id, motif.ltrim, motif.rtrim, motif.nsites);
        // draw given motif
        canvas = $("logo_popup_canvas");
        logo = logo_1(alphabet, "", pspm);
        draw_logo_on_canvas(logo, canvas, false, 0.5);
        // draw RC motif
        canvas = $("logo_popup_canvas_rc");
        if (alphabet.has_complement()) {
          pspm.reverse_complement(alphabet);
          logo = logo_1(alphabet, "", pspm);
          draw_logo_on_canvas(logo, canvas, false, 0.5);
        } else {
          canvas.style.display = "none"; 
        }
        // record the displayed motif
        make_hover_handler.motif_idx = secondary.idx;
        make_hover_handler.timer = null;
        // show the popup
        $("logo_popup").style.display = "block";
      }, 200);
    }
  };
}

function dehover_handler() {
  var popup;
  if (make_hover_handler.timer) {
    clearTimeout(make_hover_handler.timer);
    make_hover_handler.timer = null;
  }
  popup = $("logo_popup");
  popup.style.display = "none";
  this.removeEventListener('mousemove', move_logo, false);
}

/*
 * move_logo
 * 
 * keeps the motif logo at a set distance from the cursor.
 */
function move_logo(e) {
  var popup = $("logo_popup");
  popup.style.left = (e.pageX + 20) + "px";
  popup.style.top = (e.pageY + 20) + "px";
}

function make_lock_handler(secondary) {
  return function(e) {
    secondary.locked = this.checked;
    e.stopPropagation();
  };
}

function create_secondary_row(manager, secondary) {
  // get the best spacing (if one exists)
  var best_sp = (secondary.spacings.length > 0 ? secondary.spacings[0] : null);
  // create the string describing the evalue
  var evalue_str = (best_sp.pvalue * manager.nmotifs).toExponential(2);
  // create the string describing the best gap
  var gap = (best_sp.bin * data.options.bin_size);
  var gap_end = gap + data.options.bin_size - 1;
  var gap_str = "" + gap + (gap_end > gap ? " - " + gap_end : "");
  // create the string describing the orientation
  var orient_str = orient_desc[best_sp.orient];
  // create the row
  var tr = document.createElement("tr");
  tr.className = "sa_row";
  var chk_lock = document.createElement("input");
  chk_lock.type = "checkbox";
  chk_lock.checked = secondary.locked;
  chk_lock.addEventListener("click", make_lock_handler(secondary), false);
  add_cell(tr, chk_lock, "sa_lock");
  add_text_cell(tr, secondary.motif.id, "sa_id"); // secondary motif id
  add_text_cell(tr, secondary.motif.alt, "sa_id"); // secondary motif alt
  add_text_cell(tr, secondary.cluster.motif.id, "sa_cluster"); // cluster
  add_text_cell(tr, evalue_str, "sa_evalue");// secondary best spacing E-value
  add_text_cell(tr, gap_str, "sa_gap");
  add_text_cell(tr, orient_str);
  add_cell(tr, make_compact_graph(secondary), "sa_spacings"); // spacing diagram // make_minimal_spacing_diagram(secondary.spacings, 0, 150)
  tr.addEventListener("click", make_secondary_handler(manager, secondary, tr), false);
  tr.addEventListener("mouseover", make_hover_handler(secondary), false);
  tr.addEventListener("mouseout", dehover_handler, false);
  return tr;
}

function update_secondaries_list(manager) {
  var table, tbody, tr, i, last;
  var secondaries;
  // now get the elements we're manipulating
  table = manager.block.querySelector(".sa_table");
  // clear the table
  for (i = table.tBodies.length - 1; i >= 0; i--) {
    table.removeChild(table.tBodies[i]);
  }
  // now make a filtered copy of the secondaries
  secondaries = manager.primary.secondaries.filter(function (elem) {
    if (elem === manager.secondary) return true;
    if (elem.locked) return true;
    if (manager.filter_id != null && !manager.filter_id.test(elem.motif.id)) return false;
    if (manager.filter_name != null && !manager.filter_name.test(typeof elem.motif.alt === "string" ? elem.motif.alt : "")) return false;
    if (manager.filter_cluster != null && !manager.filter_cluster.test(elem.cluster.motif.id)) return false;
    if (manager.filter_pvalue != null && elem.spacings[0].pvalue > manager.filter_pvalue) return false;
    if (manager.filter_bins != null) {
      var pv = (manager.filter_pvalue == null ? 1 : manager.filter_pvalue);
      var i, spacing;
      for (i = 0; i < elem.spacings.length; i++) {
        spacing = elem.spacings[i];
        // exit loop if any spacings pass filter
        if (manager.filter_bins[spacing.bin] && spacing.pvalue < pv) break;
      }
      if (i == elem.spacings.length) return false;
    }
    return true;
  });
  // sort
  secondaries.sort(manager.sort);
  // limit to top N secondaries keeping any locked items that would otherwise be excluded
  if (manager.filter_top != null && secondaries.length > manager.filter_top) {
    last = secondaries.length;
    for (i = secondaries.length - 1; i >= manager.filter_top; i--) {
      if (secondaries[i].locked || secondaries[i] === manager.secondary) {
        if ((i+1) < last) secondaries.splice(i+1, last - (i+1));
        last = i;
      }
    }
    if ((i+1) < last) secondaries.splice(i+1, last - (i+1));
  }
  tbody = document.createElement("tbody");
  for (i = 0; i < secondaries.length; i++) {
    tbody.appendChild(create_secondary_row(manager, secondaries[i]));
  }
  table.appendChild(tbody);
}

function make_contr_seqs_format_handler(manager) {
  return function(e) {
    update_sequence_list(manager);
  };
}

function make_alignment_scroll_handler(scroller, up_pos, dn_pos, up_neg, dn_neg) {
  return function(e) {
    var left = scroller.scrollLeft;
    up_pos.scrollLeft = left;
    dn_pos.scrollLeft = left;
    up_neg.scrollLeft = left;
    dn_neg.scrollLeft = left;
  };
}

function parse_filter_bins(field) {
  if (field.disabled) return null;
  var value = field.value;
  if (!/^[0-9,\- ]*$/.test(value)) return null;
  try {
    // set all the bins to off
    var bins = [];
    var i;
    for (i = 0; i < data.options.bin_pvalue_calc_range; i += data.options.bin_size) bins.push(false);
    // split the input into parts
    var re = /(\d+|-|,)/g;
    var m;
    var items = [];
    while (m = re.exec(value)) items.push(m[1]);
    // process the parts
    var last = -1; // track the last used item
    for (i = 0; i < items.length; i++) {
      var num1, num2;
      if (items[i] == "-") { // found a dash indicating a range
        // check that there is a number before that hasn't been used already
        if ((i - 1) <= last) throw new Error("Range start already used");
        if (!/^\d+$/.test(items[i - 1])) throw new Error("Range start not a number");
        num1 = parseInt(items[i - 1]);
        // check that there is a number after
        if ((i + 1) >= items.length) throw new Error("Range end not available");
        if (!/^\d+$/.test(items[i + 1])) throw new Error("Range end not a number");
        num2 = parseInt(items[i + 1]);
        // swap them if the first is larger than the second
        if (num1 > num2) {
          var temp = num1;
          num1 = num2;
          num2 = temp;
        }
        // update the last used item and skip over the second number we used
        last = ++i;
      } else if (items[i] == ",") { // found a comma which we can safely skip
        // update the last used item
        last = i;
        continue;
      } else { // by a process of elimination this must be a number
        // check to see if there is a next item and if it is a dash indicating a range
        if ((i + 1) < items.length && items[i + 1] == "-") continue; // range upcomming so don't process this one yet
        // parse the number as a range of size one
        num1 = parseInt(items[i]);
        num2 = num1;
        // update the last used item
        last = i;
      }
      // mark any bins touched by the gap range
      for (i = num1; i <= num2; i++) {
        bins[Math.floor(i / data.options.bin_size)] = true;
      }
    }
    // return the filtered bins
    return bins;
  } catch (err) {
    return null;
  }
}

function make_sort_handler(manager) {
  return function(e) {
    var sort_input, filter_top, filter_id, filter_name, filter_cluster,
        filter_evalue, value;
    var sci_num_re = /^\s*([+]?\d+(?:\.\d+)?)(?:[eE]([+-]\d+))?\s*$/;
    // set number of results to return (excluding locked)
    filter_top = manager.block.querySelector(".sa_filter_ipt_top");
    if (!filter_top.disabled && /^\s*\d+\s*$/.test(value = filter_top.value)) {
      manager.filter_top = parseInt(value);
    } else {
      manager.filter_top = null;
    }
    // set filter on ID
    filter_id = manager.block.querySelector(".sa_filter_ipt_id");
    if (!filter_id.disabled) {
      try { manager.filter_id = new RegExp(filter_id.value, "i"); }
      catch (err) { manager.filter_id = null; }
    } else {
      manager.filter_id = null;
    }
    // set filter on name
    filter_name = manager.block.querySelector(".sa_filter_ipt_name");
    if (!filter_name.disabled) {
      try { manager.filter_name = new RegExp(filter_name.value, "i"); }
      catch (err) { manager.filter_name = null; }
    } else {
      manager.filter_name = null;
    }
    // set filter on cluster
    filter_cluster = manager.block.querySelector(".sa_filter_ipt_cluster");
    if (!filter_cluster.disabled) {
      try { manager.filter_cluster = new RegExp(filter_cluster.value, "i"); }
      catch (err) { manager.filter_cluster = null; }
    } else {
      manager.filter_cluster = null;
    }
    // set filter on evalue
    filter_evalue = manager.block.querySelector(".sa_filter_ipt_ev");
    if (!filter_evalue.disabled && sci_num_re.test(value = filter_evalue.value)) {
      manager.filter_pvalue = parseFloat(value) / manager.nmotifs;
    } else {
      manager.filter_pvalue = null;
    }
    // set filter on bins
    manager.filter_bins = parse_filter_bins(manager.block.querySelector(".sa_filter_ipt_gap"));
    // set sort function
    var sort_fns = [cmp_secondaries_by_id, cmp_secondaries_by_alt,
             cmp_secondaries_by_cluster, make_cmp_secondaries_by_evalue, 
             cmp_secondaries_by_gap, cmp_secondaries_by_orientation,
             cmp_secondaries_by_spacing];
    sort_input = manager.block.querySelector(".sa_sort");
    manager.sort = sort_fns[parseInt(sort_input.value, 10)](manager.filter_bins);
    // update the displayed table of results
    update_secondaries_list(manager);
  };
}

function setup_opt_input(manager, id, chk_idclass, lbl_idclass, ipt_idclass) {
  var chk, lbl, ipt;
  chk = manager.block.querySelector("." + chk_idclass);
  lbl = manager.block.querySelector("." + lbl_idclass);
  ipt = manager.block.querySelector("." + ipt_idclass);
  ipt.id = id;
  lbl.htmlFor = ipt.id;
  // add listeners
  chk.addEventListener("click", function() {
    ipt.disabled = !chk.checked;
  }, false);
  lbl.addEventListener("click", function() {
    if (!chk.checked) chk.click();
  }, false);
  // handle form resets
  if (chk.form != null) {
    chk.form.addEventListener("reset", function() {
      window.setTimeout(function() {
        ipt.disabled = !chk.checked;
      }, 50);
    }, false);
  }
  // set to current state
  ipt.disabled = !chk.checked;
}

function setup_filter_sort(manager, primaryi) {
  var sort_input, sort_label, update_button;
  // setup sort selector
  sort_input = manager.block.querySelector(".sa_sort");
  sort_label = manager.block.querySelector(".sa_sort_lbl");
  sort_input.id = "sort_" + primaryi;
  sort_label.htmlFor = sort_input.id;
  // filters
  setup_opt_input(manager, "filter_top_" + primaryi, "sa_filter_chk_top", "sa_filter_lbl_top", "sa_filter_ipt_top");
  setup_opt_input(manager, "filter_id_" + primaryi, "sa_filter_chk_id", "sa_filter_lbl_id", "sa_filter_ipt_id");
  setup_opt_input(manager, "filter_name_" + primaryi, "sa_filter_chk_name", "sa_filter_lbl_name", "sa_filter_ipt_name");
  setup_opt_input(manager, "filter_cluster_" + primaryi, "sa_filter_chk_cluster", "sa_filter_lbl_cluster", "sa_filter_ipt_cluster");
  setup_opt_input(manager, "filter_ev_" + primaryi, "sa_filter_chk_ev", "sa_filter_lbl_ev", "sa_filter_ipt_ev");
  // both the pvalue and gap fields work on the same checkbox
  setup_opt_input(manager, "filter_gap_" + primaryi, "sa_filter_chk_pv_gap", "sa_filter_lbl_gap", "sa_filter_ipt_gap");
  // update button
  update_button = manager.block.querySelector(".sa_update_filter_sort");
  update_button.addEventListener("click", make_sort_handler(manager), false);
}

function setup_highlight(manager, primaryi) {
  var all_chk, sel_chk, all_lbl, sel_lbl;
  all_chk = manager.block.querySelector(".sa_hl_all_chk");
  all_chk.id = "hl_all_" + primaryi;
  all_lbl = manager.block.querySelector(".sa_hl_all_lbl");
  all_lbl.htmlFor = all_chk.id;
  all_chk.addEventListener("click", function(e) {
    manager.highlight_all = this.checked;
    draw_overview_graph(manager);
    draw_sg(manager);
  }, false);
  sel_chk = manager.block.querySelector(".sa_hl_sel_chk");
  sel_chk.id = "hl_sel_" + primaryi;
  sel_lbl = manager.block.querySelector(".sa_hl_sel_lbl");
  sel_lbl.htmlFor = sel_chk.id;
  sel_chk.addEventListener("click", function(e) {
    manager.highlight_selected = this.checked;
    draw_overview_graph(manager);
    draw_sg(manager);
  }, false);
}

function make_spacing_analysis() {
  "use strict";
  var template, container, block, checkbox, label;
  var primary, orient, i, nmotifs;
  var manager;
  // hide things when no reverse complement avaliable
  if (alphabet.has_complement()) {
    document.getElementsByTagName('body')[0].className += " revcomp";
  }
  // count total motifs so we can do E-value conversions
  nmotifs = 0;
  for (i = 0; i < data.secondary_dbs.length; i++) {
    nmotifs += (data.secondary_dbs[i].loaded - data.secondary_dbs[i].excluded);
  }
  // get the template and hide it
  template = document.getElementById("tmpl_primary");
  // get the container
  container = template.parentNode;
  // empty the container of everything but the primary template
  // note that the primary template is hidden by css based on its ID
  container.removeChild(template);
  container.innerHTML = "";
  container.appendChild(template);
  // loop over primary motifs
  for (i = 0; i < data.primaries.length; i++) {
    primary = data.primaries[i];
    // make a copy of the template
    block = template.cloneNode(true);
    // set the id
    block.id = "primary_" + i;
    // create an object to track the selected objects
    manager = {};
    manager.nmotifs = nmotifs;
    manager.block = block;
    manager.primary = primary;
    manager.secondary = primary.secondaries[0];
    manager.spacing = primary.secondaries[0].spacings[0];
    manager.sort = cmp_secondaries_by_evalue;
    manager.filter_top = null;
    manager.filter_id = null;
    manager.filter_name = null;
    manager.filter_cluster = null;
    manager.filter_evalue = null;
    manager.filter_bins = null;
    manager.highlight_selected = true;
    manager.highlight_all = true;
    // store the manager in the block for easy access
    block.data_manager = manager;
    // set the name
    block.querySelector(".sa_primary_name").appendChild(make_name(primary.motif));
    // make the logo
    make_logo(block.querySelector(".sa_plogo_in"), primary.motif);
    // setup options
    setup_filter_sort(manager, i);
    setup_highlight(manager, i);
    // add event listener to format selector
    block.querySelector(".sa_contr_seqs_format").addEventListener("click", make_contr_seqs_format_handler(manager), false);
    block.querySelector(".sa_align_scrollbox").addEventListener("scroll",
        make_alignment_scroll_handler(
          block.querySelector(".sa_align_scrollbox"),
          block.querySelector(".sa_scroll_up_pos"),
          block.querySelector(".sa_scroll_dn_pos"),
          block.querySelector(".sa_scroll_up_neg"),
          block.querySelector(".sa_scroll_dn_neg")
          ), false);
    container.appendChild(block);
    // populate the list of secondaries
    update_secondaries_list(manager);
    update_selected_secondary(manager);
  }
}

function make_program_summary() {
  "use strict";
  $("version").appendChild(document.createTextNode(data.version));
  $("release").appendChild(document.createTextNode(data.release));
  $("cmd").value = data.cmd.join(" ");
  $("runtime").appendChild(document.createTextNode(data.run_time.real));
}

/*
 * cmp_spacing_by_pvalue
 */
function cmp_spacing_by_pvalue(a, b) {
  if (a == null || b == null) {
    if (a == null && b == null) {
      return 0;
    } else if (a == null) {
      return 1;
    } else if (b == null) {
      return -1;
    }
  }
  if (a.pvalue < b.pvalue) {
    return -1;
  } else if (a.pvalue > b.pvalue) {
    return 1;
  }
  if (a.bin < b.bin) {
    return -1;
  } else if (a.bin > b.bin) {
    return 1;
  }
  if (a.orient < b.orient) {
    return -1
  } else if (a.orient > b.orient) {
    return 1
  }
  return 0;
}

/*
 * cmp_secondaries_by_evalue
 */
function cmp_secondaries_by_evalue(a, b) {
  var cmp;
  if (a === b) return 0;
  if ((cmp = cmp_spacing_by_pvalue(a.spacings[0], b.spacings[0])) != 0) return cmp;
  if (a.idx < b.idx) {
    return -1;
  } else if (a.idx > b.idx) {
    return 1;
  }
  return 0;
}

/*
 * cmp_secondaries_by_id
 */
function cmp_secondaries_by_id(allowed_bins) {
  var ev_cmp = make_cmp_secondaries_by_evalue(allowed_bins);
  return function(a, b) {
    var a_alt, b_alt;
    if (a.motif.id < b.motif.id) {
      return -1;
    } else if (a.motif.id > b.motif.id) {
      return 1;
    }
    a_alt = (typeof a.motif.alt === "string" ? a.motif.alt : "");
    b_alt = (typeof b.motif.alt === "string" ? b.motif.alt : "");
    if (a_alt < b_alt) {
      return -1;
    } else if (a_alt > b_alt) {
      return 1;
    }
    return ev_cmp(a, b);
  };
}

/*
 * cmp_secondaries_by_alt
 */
function cmp_secondaries_by_alt(allowed_bins) {
  var ev_cmp = make_cmp_secondaries_by_evalue(allowed_bins);
  return function(a, b) {
    var a_alt, b_alt;
    a_alt = (typeof a.motif.alt === "string" ? a.motif.alt : "");
    b_alt = (typeof b.motif.alt === "string" ? b.motif.alt : "");
    if (a_alt < b_alt) {
      return -1;
    } else if (a_alt > b_alt) {
      return 1;
    }
    if (a.motif.id < b.motif.id) {
      return -1;
    } else if (a.motif.id > b.motif.id) {
      return 1;
    }
    return ev_cmp(a, b);
  };
}

/*
 * cmp_secondaries_by_cluster
 */
function cmp_secondaries_by_cluster(allowed_bins) {
  var ev_cmp = make_cmp_secondaries_by_evalue(allowed_bins);
  return function(a, b) {
    var cluster_cmp;
    cluster_cmp = ev_cmp(a.cluster, b.cluster);
    if (cluster_cmp != 0) return cluster_cmp;
    return ev_cmp(a, b);
  };
}

/*
 * cmp_secondaries_by_gap
 */
function cmp_secondaries_by_gap(allowed_bins) {
  var ev_cmp = make_cmp_secondaries_by_evalue(allowed_bins);
  return function(a, b) {
    var a_space, b_space;
    if (a === b) return 0;
    a_space = a.spacings[0];
    b_space = b.spacings[0];
    if (a_space.bin < b_space.bin) {
      return -1;
    } else if (a_space.bin > b_space.bin) {
      return 1;
    }
    return ev_cmp(a, b);
  };
}

/*
 * cmp_secondaries_by_orientation
 */
function cmp_secondaries_by_orientation(allowed_bins) {
  var ev_cmp = make_cmp_secondaries_by_evalue(allowed_bins);
  return function(a, b) {
    var a_space, b_space;
    if (a === b) return 0;
    a_space = a.spacings[0];
    b_space = b.spacings[0];
    if (a_space.orient < b_space.orient) {
      return -1
    } else if (a_space.orient > b_space.orient) {
      return 1
    }
    return ev_cmp(a, b);
  };
}

/*
 * cmp_secondaries_by_spacing
 */
function cmp_secondaries_by_spacing(allowed_bins) {
  return function(a, b) {
    var a_space, b_space;
    if (a === b) return 0;
    a_space = a.spacings[0];
    b_space = b.spacings[0];
    if (a_space.bin < b_space.bin) {
      return -1;
    } else if (a_space.bin > b_space.bin) {
      return 1;
    }
    if (a_space.orient < b_space.orient) {
      return -1
    } else if (a_space.orient > b_space.orient) {
      return 1
    }
    if (a_space.pvalue < b_space.pvalue) {
      return -1;
    } else if (a_space.pvalue > b_space.pvalue) {
      return 1;
    }
    if (a.idx < b.idx) {
      return -1;
    } else if (a.idx > b.idx) {
      return 1;
    }
    return 0;
  };
}

function make_cmp_secondaries_by_evalue(allowed_bins) {
  function _best_allowed_spacing(secondary) {
    var i;
    if (allowed_bins == null) return secondary.spacings[0];
    for (i = 0; i < secondary.spacings.length; i++) {
      if (allowed_bins[secondary.spacings[i].bin]) return secondary.spacings[i];
    }
    return null;
  }
  return function(a, b) {
    var a_pv, b_pv;
    var cmp;
    if (a === b) return 0;
    if ((cmp = cmp_spacing_by_pvalue(_best_allowed_spacing(a), _best_allowed_spacing(b))) != 0) return cmp;
    if (a.idx < b.idx) {
      return -1;
    } else if (a.idx > b.idx) {
      return 1;
    }
    return 0;
  };
}

function show_motif_pair(pri_db, pri_id, sec_db, sec_id, sp_orient, sp_bin) {
  var i, primaryi, primary, secondary, spacing, block, manager;
  // first determine which primary
  for (i = 0; i < data.primaries.length; i++) {
    primary = data.primaries[i];
    if (primary.motif.db == pri_db && primary.motif.id == pri_id) {
      primaryi = i;
      break;
    }
  }
  if (i == data.primaries.length) return; // no such primary
  // determine secondary
  if (typeof sec_db == "number" && typeof sec_id == "string") {
    for (i = 0; i < primary.secondaries.length; i++) {
      secondary = primary.secondaries[i];
      if (secondary.motif.db == sec_db && secondary.motif.id == sec_id) break;
    }
    if (i == primary.secondaries.length) return; // no such secondary
  } else {
    secondary = primary.secondaries[0];
  }
  // determine spacing
  if (typeof sp_orient == "number" && typeof sp_gap == "number") {
    for (i = 0; i < secondary.spacings.length; i++) {
      spacing = secondary.spacings[i];
      if (spacing.orient == sp_orient && spacing.bin == sp_bin) break;
    }
    if (i == secondary.spacings.length) return; // no such spacing
  } else {
    spacing = secondary.spacings[0];
  }
  // get the block
  block = document.getElementById("primary_" + primaryi);
  if (!block) return;
  // get the manager
  manager = block.data_manager;
  manager.primary = primary;
  manager.secondary = secondary;
  manager.spacing = spacing;
  update_secondaries_list(manager);
}

/*
 * Process data
 *
 * Preprocess the data
 */
(function() {
  var primaryi, secondaryi, secondaryj, i;
  var primary, secondary, secondaries;
  for (primaryi = 0; primaryi < data.primaries.length; primaryi++) {
    primary = data.primaries[primaryi];
    secondaries = [];
    for (secondaryi = 0; secondaryi < primary.secondaries.length; secondaryi++) {
      data.secondary_motifs[primary.secondaries[secondaryi][0].idx].nonredundant = true;
      for (secondaryj = 0; secondaryj < primary.secondaries[secondaryi].length; secondaryj++) {
        secondary = primary.secondaries[secondaryi][secondaryj];
        secondary.motif = data.secondary_motifs[secondary.idx];
        secondary.cluster = primary.secondaries[secondaryi][0];
        secondaries.push(secondary);
      }
    }
    secondaries.sort(cmp_secondaries_by_evalue);
    primary.secondaries = secondaries;
  }
  for (i = 0; i < data.secondary_motifs.length; i++) {
    data.secondary_motifs[i].locked = false;
  }
})();
