
var motifs = null;
var sequences = null;

function register_component(id, element, controler) {
  "use strict";
  if (id == "alphabet") {
    alphabet = controler;
    element.addEventListener("alphabet_changed", function (e) {
      if (sequences != null) sequences.set_custom_alphabet(e.detail.has_custom, e.detail.alphabet);
    }, false);
  } else if (id == "motifs") {
    motifs = controler;
    element.addEventListener("motifs_loaded", function() {
      if (sequences != null) {
        sequences.set_expected_alphabet(motifs.get_alphabet());
      }
    }, false);
  } else if (id == "sequences") {
    sequences = controler;
    if (motifs != null) {
      sequences.set_expected_alphabet(motifs.get_alphabet());
    }
  }
}

function check() {
  "use strict";
  var alphs = null, seq_alphabet, translate_dna;
  if (alphabet != null) {
    alphs = alphabet.get_alphabets();
  }
  translate_dna = $("translate_dna").checked;
  if (motifs == null || sequences == null) {
    alert("A Javascript error occurred - please try clearing the cache and refreshing the page." +
        " If the problem persists then contact the support email for a solution.");
    return false;
  }
  if (motifs != null) {
    if (!motifs.check(alphs, alphabet != null && alphabet.get_alphabets() != null)) return false;
    if (alphs == null) alphs = motifs.get_alphabets();
  }
  if (translate_dna) {
    if (!AlphStd.PROTEIN.equals(motifs.get_alphabet())) {
      alert("Please provide protein motifs when using the translate DNA option.");
      return false;
    }
    var seq_alphabets = sequences.get_alphabets();
    if (seq_alphabets != null && !AlphabetUtil.any_equal(seq_alphabets, [AlphStd.DNA, AlphStd.RNA])) {
      alert("Please provide DNA or RNA sequences when using the translate DNA option.");
      return false;
    }
  }
  else if (sequences != null) {
    if (!sequences.check(alphs)) return false;
    if (alphs == null) alphs = sequences.get_alphabets();
  }
  if (!check_job_details()) return false;
  if (!check_num_value("sequence E-value threshold", "seq_evalue", 0, null, 10)) return false;
  if (!check_num_value("motif E-value threshold", "motif_evalue", 0, null, 0.05)) return false;
  return true;
}


function options_changed() {
  if ($("strands").value !== "combine") return true;
  if ($("translate_dna").checked) return true;
  if (!/^\s*10\s*$/.test($("seq_evalue").value)) return true;
  if ($("use_seq_comp").checked) return true;
  if ($("scale_m").checked) return true;
  if ($("motif_evalue_enable").checked) return true;
  if (! $("rem_corr").checked) return true;
  return false;
}

function options_reset(evt) {
  $("strands").value = "combine";
  $("translate_dna").checked = false;
  $("seq_evalue").value = "10";
  $("use_seq_comp").checked = false;
  $("scale_m").checked = false;
  $("motif_evalue_enable").checked = false;
  $("motif_evalue").disabled = true;
  $("motif_evalue").value = "0.05";
  $("rem_corr").checked = true;
}

function fix_reset() {
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
  // add listener to the form to check the fields before submit
  $("mast_form").addEventListener("submit", on_form_submit, false);
  $("mast_form").addEventListener("reset", on_form_reset, false);
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
