var tomtom_alphabet = new Alphabet(data.alphabet, data.background);

var DelayLogoTask = function(logo, canvas) {
  this.logo = logo;
  this.canvas = canvas;
};

DelayLogoTask.prototype.run = function () {
  draw_logo_on_canvas(this.logo, this.canvas);
};

function make_logo(alphabet, pspm, rc, offset, className) {
  if (rc) pspm = pspm.copy().reverse_complement(alphabet);
  var logo = new Logo(alphabet, "");
  logo.add_pspm(pspm, offset);
  var canvas = document.createElement('canvas');
  canvas.height = 50;
  canvas.width = 0;
  canvas.className = className;
  size_logo_on_canvas(logo, canvas);
  add_draw_task(canvas, new DelayLogoTask(logo, canvas));
  return canvas;
}

function make_name(motif, should_link) {
  var name, link;
  name = motif.id + (typeof motif.alt === "string" ? " (" + motif.alt + ")" : "");
  if (should_link && typeof motif.url === "string") {
    link = document.createElement("a");
    link.href = motif.url;
    link.appendChild(document.createTextNode(name));
    return link;
  } else {
    return document.createTextNode(name);
  }
}

function make_logo_form() {
  $("logo_form").action = site_url + "/utilities/generate_logo";
  $("download_submit").addEventListener("click", function() {
    toggle_class(document.getElementById("download"), "hidden", true);
  }, false);
  $("download_cancel").addEventListener("click", function() {
    toggle_class(document.getElementById("download"), "hidden", true);
  }, false);
}

function make_alignment(query, target, rc, offset) {
  var logo, query_pspm, target_pspm;
  logo = new Logo(tomtom_alphabet, "");
  query_pspm = new Pspm(query.pwm, query.id, 0, 0, query.nsites);
  target_pspm = new Pspm(target.pwm, target.id, 0, 0, target.nsites);
  if (rc) target_pspm.reverse_complement(tomtom_alphabet);
  logo.add_pspm(target_pspm, (offset < 0 ? -offset : 0));
  logo.add_pspm(query_pspm, (offset < 0 ? 0 : offset));
  var canvas = document.createElement('canvas');
  canvas.height = 300;
  canvas.width = 0;
  canvas.title = "Click to show download options";
  size_logo_on_canvas(logo, canvas);
  add_draw_task(canvas, new DelayLogoTask(logo, canvas));
  return canvas;
}

function make_preview() {
  var i, j, qidx, tidx, query, target, match, matches, pview_table, pview_tbody, pview_row;
  var row, row_mlist, mlink;
  // scan the matches make a lookup from query index to matches index
  var lookup = {};
  for (i = 0; i < data.all_matches.length; i++) {
    lookup[data.all_matches[i].idx] = i;
  }

  pview_table = $("preview");
  pview_row = pview_table.querySelector(".pview_row");
  pview_tbody = pview_row.parentNode;
  pview_tbody.removeChild(pview_row);
  for (qidx = 0; qidx < data.queries.length; qidx++) {
    query = data.queries[qidx];
    matches = (lookup[qidx] != null ? data.all_matches[lookup[qidx]].matches : []);
    row = pview_row.cloneNode(true);
    row.querySelector(".pview_db").appendChild(document.createTextNode(data.query_dbs[query.db].name));
    row.querySelector(".pview_name").appendChild(document.createTextNode(query.id));
    if (typeof query.alt === "string") row.querySelector(".pview_alt").appendChild(document.createTextNode(query.alt));
    row.querySelector(".pview_logo").appendChild(make_logo(tomtom_alphabet, new Pspm(query.pwm)));
    row.querySelector(".pview_matches").appendChild(document.createTextNode(matches.length));
    row_mlist = row.querySelector(".pview_list");
    for (j = 0; j < matches.length; j++) {
      if (j > 0) row_mlist.appendChild(document.createTextNode(",\u2003 "));// EM Space
      match = matches[j];
      tidx = match.idx;
      target = data.targets[tidx];
      mlink = document.createElement("a");
      mlink.className = (qidx % 2 == 0 ? "ml1" : "ml2");
      mlink.appendChild(document.createTextNode(target.id));
      if (typeof target.alt === "string") mlink.appendChild(document.createTextNode("\u00A0(" + target.alt + ")"));
      mlink.href = "#match_" + qidx + "_" + tidx;
      row_mlist.appendChild(mlink);
    }
    pview_tbody.appendChild(row);
  }
}

function make_target_dbs() {
  var i, tdb_table, tdb_row, tdb_tbody, db;
  var row, match_count;
  tdb_table = $("tdbs");
  tdb_row = tdb_table.querySelector(".db_row");
  tdb_tbody = tdb_row.parentNode;
  tdb_tbody.removeChild(tdb_row);
  match_count = [];
  for (i = 0; i < data.target_dbs.length; i++) match_count[i] = 0;
  for (i = 0; i < data.targets.length; i++) {
    match_count[data.targets[i].db]++;
  }
  for (i = 0; i < data.target_dbs.length; i++) {
    db = data.target_dbs[i];
    row = tdb_row.cloneNode(true);
    row.querySelector(".db_name").appendChild(document.createTextNode(db.name));
    row.querySelector(".db_used").appendChild(document.createTextNode(db.loaded - db.excluded));
    row.querySelector(".db_matched").appendChild(document.createTextNode(match_count[i]));
    tdb_tbody.appendChild(row);
  }
  $("link_after_target_db").href = (data.all_matches.length > 0 ? "#query_" + data.all_matches[0].idx : "#settings");
}

function custom_logo(elem, qidx, tidx, rc, offset) {
  "use strict";
  function helper(qidx, tidx, rc, offset) {
    "use strict";
    return function() {
      if ($("dl_flip").value == "1") {
        $("dl_rc1").value = (rc ? "0" : "1");
        $("dl_rc2").value = "1";
        $("dl_shift").value = -(data.targets[tidx].len - (data.queries[qidx].len + offset));
      } else {
        $("dl_rc1").value = (rc ? "1" : "0");
        $("dl_rc2").value = "0";
        $("dl_shift").value = -offset;
      }
    };
  }
  var download_prompt, qidx, tidx, rc, offset, query, target, qmotif, tmotif;
  var flip_handler;
  // get the download prompt
  download_prompt = document.getElementById("download");
  // hide previous showing of the download prompt
  toggle_class(download_prompt, "hidden", true);
  if (typeof custom_logo.flip_handler !== "undefined") {
    $("dl_flip").removeEventListener("change", custom_logo.flip_handler, false);
  }
  // reset defaults on prompt
  $("dl_flip").value = "0";
  $("dl_width").value = "";
  $("dl_height").value = "";
  // configure download prompt to the selected logo
  query = data.queries[qidx];
  target = data.targets[tidx];
  tmotif = new Pspm(target.pwm, "1", 0, 0, target.nsites, target.evalue);
  qmotif = new Pspm(query.pwm, "2", 0, 0, query.nsites, query.evalue);
  $("dl_motifs").value = tmotif.as_meme({"with_header": true, "alphabet": tomtom_alphabet}) + qmotif.as_meme();
  $("dl_label1").value = target.id;
  $("dl_label2").value = query.id;
  flip_handler = helper(qidx, tidx, rc, offset);
  $("dl_flip").addEventListener("change", flip_handler, false);
  custom_logo.flip_handler = flip_handler;
  flip_handler();
  // position and display the download prompt
  position_popup(elem, download_prompt, 1);
  toggle_class(download_prompt, "hidden", false);
  // focus the download button
  $("download_submit").focus();
}

function make_custom_logo_handler(elem, qidx, tidx, rc, offset) {
  return function(evt) {
    custom_logo(elem, qidx, tidx, rc, offset);
    evt.preventDefault();
  };
}

function make_matches() {
  var i, j, matches_box, match_box, match_entry, box, links, entries, entry;
  var qidx, query, tidx, target, matches, overlap, rc, logo_container;
  var prev_query_link, next_query_link;
  matches_box = $("matches");
  match_box = matches_box.querySelector(".match_box");
  matches_box.removeChild(match_box);
  for (i = 0; i < data.all_matches.length; i++) {
    matches = data.all_matches[i].matches;
    qidx = data.all_matches[i].idx;
    query = data.queries[qidx];
    box = match_box.cloneNode(true);
    box.id = "query_" + qidx;
    box.querySelector(".query_name").appendChild(make_name(query, false));
    links = box.querySelector(".links");
    prev_query_link = (i == 0 ? "#target_dbs" : "#query_" + data.all_matches[i-1].idx);
    next_query_link = ((i + 1) >= data.all_matches.length ? "#settings" : "#query_" + data.all_matches[i+1].idx);
    links.querySelector(".prev").href = prev_query_link;
    links.querySelector(".next").href = next_query_link;
    match_entry = box.querySelector(".match_entry");
    entries = match_entry.parentNode;
    entries.removeChild(match_entry);
    for (j = 0; j < matches.length; j++) {
      tidx = matches[j].idx;
      target = data.targets[tidx];
      entry = match_entry.cloneNode(true);
      entry.id = "match_" + qidx + "_" + tidx;
      entry.querySelector(".match_name").appendChild(make_name(target, true));
      entry.querySelector(".match_db").appendChild(document.createTextNode(data.target_dbs[target.db].name));
      entry.querySelector(".match_pvalue").appendChild(document.createTextNode(matches[j].pv));
      entry.querySelector(".match_evalue").appendChild(document.createTextNode(matches[j].ev));
      entry.querySelector(".match_qvalue").appendChild(document.createTextNode(matches[j].qv));
      if (matches[j].off < 0) {
        overlap = Math.min(query.pwm.length + matches[j].off, target.pwm.length);
      } else if (matches[j].off == 0) {
        overlap = Math.min(query.pwm.length, target.pwm.length);
      } else {
        overlap = Math.min(query.pwm.length, target.pwm.length - matches[j].off);
      }
      entry.querySelector(".match_overlap").appendChild(document.createTextNode(overlap));
      entry.querySelector(".match_offset").appendChild(document.createTextNode(matches[j].off));
      rc = (typeof matches[j].rc === "boolean" && matches[j].rc);
      entry.querySelector(".match_orientation").appendChild(document.createTextNode(rc ? "Reverse Complement" : "Normal"));
      logo_container = entry.querySelector(".logo_container");
      logo_container.appendChild(make_alignment(query, target, rc, matches[j].off));
      entries.appendChild(entry);
      var custom_logo_handler = make_custom_logo_handler(logo_container.querySelector("canvas"), qidx, tidx, rc, matches[j].off);
      entry.querySelector(".download_btn").addEventListener("click", custom_logo_handler, false);
      entry.querySelector(".download_text").addEventListener("click", custom_logo_handler, false);
      entry.querySelector("canvas").addEventListener("click", custom_logo_handler, false);
      if (i != 0) entry.querySelector(".prev_query").href = prev_query_link;
      if ((i + 1) < data.all_matches.length) entry.querySelector(".next_query").href = next_query_link;
      if (j != 0) {
        entry.querySelector(".prev_target").href = "#match_" + qidx + "_" + matches[j-1].idx;
      } else if (i != 0) {
        var prev_query =  data.all_matches[i-1];
        var prev_qidx = prev_query.idx;
        var prev_matches = prev_query.matches;
        var prev_tidx = prev_matches[prev_matches.length - 1].idx;
        entry.querySelector(".prev_target").href = "#match_" + prev_qidx + "_" + prev_tidx;
      }
      if ((j + 1) < matches.length) {
        entry.querySelector(".next_target").href = "#match_" + qidx + "_" + matches[j+1].idx;
      } else if ((i + 1) < data.all_matches.length) {
        var next_query = data.all_matches[i+1];
        var next_qidx = next_query.idx;
        var next_matches = next_query.matches;
        var next_tidx = next_matches[0].idx;
        entry.querySelector(".next_target").href = "#match_" + next_qidx + "_" + next_tidx;
      }
    }

    matches_box.appendChild(box);
  }

}

function make_program() {
  $("link_before_program").href = (data.all_matches.length > 0 ? "#query_" + data.all_matches[data.all_matches.length - 1].idx : "#target_dbs");
  $("version").appendChild(document.createTextNode(data.version));
  $("release").appendChild(document.createTextNode(data.release));
  $("cmd").value = data.cmd.join(" ");
  $("runtime").appendChild(document.createTextNode(data.runtime.seconds));
}
