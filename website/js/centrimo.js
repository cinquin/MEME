
var alphabet = null;
var sequences = null;
var control_sequences = null;
var motifs = null;
var background = null;

function register_component(id, element, controler) {
  "use strict";
  if (id == "alphabet") {
    alphabet = controler;
    element.addEventListener("alphabet_changed", function (e) {
      if (sequences != null) sequences.set_custom_alphabet(e.detail.has_custom, e.detail.alphabet);
      if (control_sequences != null) control_sequences.set_custom_alphabet(e.detail.has_custom, e.detail.alphabet);
    }, false);
  } else if (id == "sequences") {
    sequences = controler;
    if (alphabet != null) {
      sequences.set_custom_alphabet(alphabet.has_custom_alphabet(), alphabet.get_custom_alphabet());
    }
  } else if (id == "control_sequences") {
    control_sequences = controler;
    if (alphabet != null) {
      control_sequences.set_custom_alphabet(alphabet.has_custom_alphabet(), alphabet.get_custom_alphabet());
    }
  } else if (id == "motifs") {
    motifs = controler;
  } else if (id == "background") {
    background = controler;
  }
}

function check() {
  "use strict";
  var alphs = null;
  if (alphabet != null) {
    alphs = alphabet.get_alphabets();
  }
  if (sequences != null) {
    if (!sequences.check(alphs)) return false;
    if (alphs == null) alphs = sequences.get_alphabets();
  }
  if ($("compare_on").checked && control_sequences != null) {
    if (!control_sequences.check(alphs)) return false;
    if (alphs == null) alphs = control_sequences.get_alphabets();
  }
  if (motifs != null && !motifs.check(alphs, alphabet != null && alphabet.get_alphabets() != null)) return false;
  if (!check_job_details()) return false;
  if (!check_num_value("minimum score threshold", "min_score", null, null, 5)) return false;
  if (!check_int_value("maximum region width", "max_region", 0, null, 200)) return false;
  if (!check_num_value("E-value threshold", "evalue_threshold", 0, null, 10)) return false;
  if (background != null && !background.check(alphs)) return false;
  return true;
}

function options_changed() {
  if (background != null && background.changed()) return true;
  if ($("strands").value !== "both") return true;
  if (!/^\s*5\s*$/.test($("min_score").value)) return true;
  if ($("opt_score").checked) return true;
  if ($("use_max_region").checked) return true;
  if (!/^\s*10\s*$/.test($("evalue_threshold").value)) return true;
  if (!$("store_ids").checked) return true;
  return false;
}

function options_reset(evt) {
  //reset options
  if (background != null) background.reset();
  $("strands").value = "both";
  $("min_score").value = 5;
  $("opt_score").checked = false;
  $("use_max_region").checked = false;
  $("max_region").disabled = true;
  $("max_region").value = 200;
  $("evalue_threshold").value = 10;
  $("store_ids").checked = true;
}

function update_compare() {
  $('compare_sequences_area').style.display = ($('compare_on').checked ? 'block' : 'none');
}

function fix_reset() {
  update_compare();
}

function on_form_submit(evt) {
  if (!check()) {
    evt.preventDefault();
  }
}

function on_form_reset(evt) {
  window.setTimeout(function(evt) {
    fix_reset();
  }, 50);
}


function on_load() {
  // add listeners for the enrichment mode
  $("compare_off").addEventListener("click", update_compare, false);
  $("compare_on").addEventListener("click", update_compare, false);
  // add listener to the form to check the fields before submit
  $("centrimo_form").addEventListener("submit", on_form_submit, false);
  $("centrimo_form").addEventListener("reset", on_form_reset, false);
}

// add a load
(function() {
  "use strict";
  window.addEventListener("load", function load(evt) {
    "use strict";
    window.removeEventListener("load", load, false);
    on_load();
  }, false);
})();
