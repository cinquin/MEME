
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
  if ($("generate_off").checked) {
    if (control_sequences != null) {
      if(!control_sequences.check(alphs)) return false;
      if (alphs == null) alphs = control_sequences.get_alphabets();
    }
  }
  if (motifs != null) {
    if(!motifs.check(alphs, alphabet != null && alphabet.get_alphabets() != null)) return false;
    if (alphs == null) alphs = motifs.get_alphabets();
  }
  if (!check_num_value("sequence match threshold", "pwm_threshold", null, null, 1)) return false;
  if (!check_num_value("motif match threshold", "pvalue_threshold", 0, 1, 0.0002)) return false;
  if (!check_num_value("p-value report threshold", "pvalue_report_threshold", 0, 1, 0.05)) return false;
  if (!check_job_details()) return false;
  if (background != null && !background.check(alphs)) return false;
  return true;
}

function options_changed() {
  if ($("scoring").value != "avg") return true;
  if ($("pvalue_threshold").value != 0.0002) return true;
  if ($("method").value != "ranksum") return true;
  if ($("pwm_threshold").value != 1) return true;
  if (!/^\s*0.05\s*$/.test($("pvalue_report_threshold").value)) return true;
  if (background != null && background.changed()) return true;
  return false;
}

function options_reset(evt) {
  //reset options
  $("scoring").value = "avg";
  $("pvalue_threshold").value = 0.0002;
  $('hits_area').style.display = 'none';
  $("method").value = "ranksum";
  $("pwm_threshold").value = 1;
  $('fisher_pwm_area').style.display = 'none';
  $("pvalue_report_threshold").value = 0.05;
  if (background != null) background.reset();
}

function fix_reset() {
  on_ch_generate();
  on_ch_scoring();
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

function on_ch_generate() {
  $('control').style.display = $('generate_off').checked ? 'block' : 'none';
}

function on_ch_method() {
  var method_value = $('method').value;
  //$('scoring_area').style.display = /(fisher)/.test(method_value) ? 'block' : 'none';
}

function on_ch_method() {
  var method_value = $('method').value;
  $('fisher_pwm_area').style.display = /(fisher)/.test(method_value) ? 'block' : 'none';
}

function on_ch_scoring() {
  var scoring_value = $('scoring').value;
  $('hits_area').style.display = /(totalhits)/.test(scoring_value) ? 'block' : 'none';
}

function on_load() {
  // add listeners for the control generation mode
  $("generate_off").addEventListener("click", on_ch_generate, false);
  $("generate_on").addEventListener("click", on_ch_generate, false);
  // add listener to the form to check the fields before submit
  $("ame_form").addEventListener("submit", on_form_submit, false);
  $("ame_form").addEventListener("reset", on_form_reset, false);
  $("method").addEventListener("change", on_ch_method, false);
  $("scoring").addEventListener("change", on_ch_scoring, false);
  on_ch_scoring();
  on_ch_method();
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
