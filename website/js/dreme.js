var alphabet = null;
var sequences = null;
var control_sequences = null;
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
  } else if (id == "background") {
    background = controler;
  }
}

function check() {
  "use strict";
  var alphs = null;
  if (alphabet != null) alphs = alphabet.get_alphabets();
  if (sequences != null) {
    if (!sequences.check(alphs)) return false;
    if (alphs == null) alphs = sequences.get_alphabets();
  }
  if ($("discr_on").checked && control_sequences != null) {
    if (!control_sequences.check(alphs)) return false;
  }
  if (!check_job_details()) return false;
  if (!check_num_value("E-value threshold", "ethresh", 0, 1, 0.05)) return false;
  if (!check_int_value("motif count", "nmotifs", 1, null, 10)) return false;
  return true;
}

function options_changed() {
  if (!/^\s*0\.05\s*$/.test($("ethresh").value)) return true;
  if ($("nmotifs_enable").checked) return true;
  if ($("norc").checked) return true;
  return false;
}

function options_reset(evt) {
  $("ethresh").value = 0.05;
  $("nmotifs_enable").checked = false;
  $("nmotifs").value = 10;
  $("nmotifs").disabled = true;
  $("norc").checked = false;
}

function fix_reset() {
  $('discr_sequences_area').style.display = ($('discr_on').checked ? 'block' : 'none');
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

function on_ch_discr() {
  $('discr_sequences_area').style.display = ($('discr_on').checked ? 'block' : 'none');
}

function on_load() {
  // add listeners for the motif discovery mode
  $("discr_off").addEventListener("click", on_ch_discr, false);
  $("discr_on").addEventListener("click", on_ch_discr, false);
  // add listener to the form to check the fields before submit
  $("dreme_form").addEventListener("submit", on_form_submit, false);
  $("dreme_form").addEventListener("reset", on_form_reset, false);
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
