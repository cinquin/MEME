
var sequences = null;

function register_component(id, element, controler) {
  "use strict";
  if (id == "sequences") {
    sequences = controler;
  }
}

function on_alg_change() {
//  $("flanking_sequences").style.display = ($('alg_meme').checked ? 'none' : 'block');
  $("occurrences").style.display = ($('alg_simp').checked ? 'block' : 'none');
  $("motifx_thresholds").style.display = ($('alg_mtfx').checked ? 'block' : 'none');
//  $("num_motifs").style.display = ($('alg_meme').checked ? 'block' : 'none');
}

function on_filetype_change() {
  var selectFiletype = document.getElementsByName("fg_filetype")[0];
  $("sequence_column_name").style.display = (selectFiletype.value == "psm") ? 'block' : 'none';
  $("summary_table").style.display = (selectFiletype.value == "psm" ? 'block' : 'none');
  calculateSummary();
}

function add_psm() {
  "use strict";
  var tbody = $("psm_rows");
  var main_row = $("psm_main_row");
  var new_row = main_row.cloneNode(true);
  // reset inputs
  var inputs, options, i, j, groups;
  var num = tbody.rows.length + 1;
  var name_re = /^psm_1(.*)$/;
  new_row.id = "psm_row_" + num;
  inputs = new_row.querySelectorAll("input");
  for (i = 0; i < inputs.length; i++) {
    if (inputs[i].type == "checkbox") {
      inputs[i].checked = inputs[i].defaultChecked;
    } else {
      inputs[i].value = inputs[i].defaultValue;
    }
    if ((groups = name_re.exec(inputs[i].name)) != null) {
      inputs[i].name = "psm_" + num + groups[1];
    }
  }
  inputs = new_row.querySelectorAll("select");
  for (i = 0; i < inputs.length; i++) {
    if ((groups = name_re.exec(inputs[i].name)) != null) {
      inputs[i].name = "psm_" + num + groups[1];
    }
    options = inputs[i].querySelectorAll("select option");
    for (j = 0; j < options.length; j++) {
      options[j].selected = options[j].defaultSelected;
    }
  }
  
  // add the row
  tbody.appendChild(new_row);
  
  // add a anyFileChanged event to the the psm file selector psm_(num)
  var selectFile = document.getElementsByName("psm_" + num)[0];
  selectFile.addEventListener('change', calculateSummary, false);
}

function firstFileChanged(selectFile) {
  var filter_enable = document.getElementsByName("filter_enable")[0];
  var filter_field = document.getElementsByName("filter_field")[0];
  
  // Retrieve the first file from the FileList object
  var f = selectFile.files[0];
  
  if (f) {
    var r = new FileReader();
    r.onload = function(e) {
      var contents = e.target.result;
      var header = contents.split('\n')[0].split('\t');
      var firstline = contents.split('\n')[1].split('\t');
      
      filter_field.options.length=0;
      
      for (i = 0; i < header.length; i++) {
        if (!isNaN(firstline[i])) {
          filter_field.options.add(new Option(header[i], header[i], false, false));
        }
      }
      
      if (filter_field.options.length == 0) {
        filter_enable.checked = false;
      }
      filter_enable.disabled = (filter_field.options.length == 0);
      toggleFilter();
    }
    r.readAsText(f);
  } else {
    alert("Failed to load file");
    filter_field.options.length=0;
    filter_enable.checked = false;
    filter_enable.disabled = (filter_field.options.length == 0);
    toggleFilter();
  }
}

function updateSummary(summary) {
  var summaryTable = document.getElementById("summary_table");
  summaryTable.rows[0].cells[1].innerHTML = summary[0];
  summaryTable.rows[1].cells[1].innerHTML = summary[1];
}

function calculateSummary() {
  var summary = [0, 0];
  var modDict = {}
  
  var sequence_column_name = document.getElementsByName("psm_column_name")[0];
  var filter_enable = document.getElementsByName("filter_enable")[0];
  var filter_field = document.getElementsByName("filter_field")[0];

  var rows = document.getElementById("psm_table").rows;
  
  // for each row in the table of files, excluding the header
  for (i=1; i < rows.length; i++) {
    // get the file selector
    selectFile = document.getElementsByName("psm_" + i)[0];
    // Retrieve the first (and only!) File from the FileList object
    var files = selectFile.files;
    
    if (files) {
      for (var j=0, f; f=files[j]; j++) {
        var r = new FileReader();
        r.onload = (function(f) {
                    return function(e) {
                    var contents = e.target.result;
                    var numOfPtms = 0;
                    var fileLines = contents.split('\n');
                    var headerArray = fileLines[0].split('\t');
                    var indexOfSequence = headerArray.indexOf(sequence_column_name.value);
                    var modPattern = /.\[.*?\]/g;
                    // for each line in the file
                    for (k=1; k < fileLines.length; k++) {
                      unprocessed_sequence = fileLines[k].split('\t')[indexOfSequence];
                      sequence = unprocessed_sequence;
                      // if the peptide format is not in tide, we need to modify it...
                      if (unprocessed_sequence != undefined) {
                        if (unprocessed_sequence.length > 2 && unprocessed_sequence[1] == '.') {
                          sequence = "";
                          if (unprocessed_sequence[0] != '-') {
                            sequence += unprocessed_sequence[0];
                          }
                          sequence += unprocessed_sequence.substring(2, unprocessed_sequence.length-2);
                          if (unprocessed_sequence[unprocessed_sequence.length-1] != '-') {
                            sequence += unprocessed_sequence[unprocessed_sequence.length-1];
                          }
                        } else {
                          sequence = unprocessed_sequence;
                          plusIdx = sequence.indexOf('+');
                          minusIdx = sequence.indexOf('-');
                          if ((plusIdx != -1 && sequence[plusIdx-1] != '[') || (minusIdx != -1 && sequence[minusIdx-1] != '[')) {
                              modPattern = /.[\+|\-][0-9|\.]+/g;
                          }
                        }
                      }
                      while (match = modPattern.exec(sequence)) {
                        matchStr = match[0];
                        modDict[matchStr] = (modDict[matchStr] == undefined) ? 1 : modDict[matchStr]+1;
                      }
                    }
                    var totalMods = 0;
                    for (var key in modDict) {
                      totalMods += modDict[key];
                    }
                    summary[0] = totalMods;
                    summary[1] = Object.keys(modDict).length;
                    updateSummary(summary);
                    };
                    })(f);
        r.readAsText(f);
      }
    }
  }
}

function remove_psm() {
  "use strict";
  var tbody = $("psm_rows");
  if (tbody.rows.length > 1) {
    tbody.deleteRow(-1);
  }
  calculateSummary();
}

function toggleFilter() {
  var colheader = document.getElementById("filter_table").rows[0].cells;
  var colbody = document.getElementById("filter_table").rows[1].cells;
  var do_show = document.getElementsByName("filter_enable")[0].checked;
  colheader[1].style.display = do_show ? '' : 'none';
  colheader[2].style.display = do_show ? '' : 'none';
  colheader[3].style.display = do_show ? '' : 'none';
  colbody[1].style.display = do_show ? '' : 'none';
  colbody[2].style.display = do_show ? '' : 'none';
  colbody[3].style.display = do_show ? '' : 'none';

}

function check() {
  "use strict";
  if (sequences != null) {
    if (!sequences.check()) return false;
  }
  if (!check_job_details()) return false;
  if (!check_int_value("motif width", "width", 1, 51, 7)) return false;
  return true;
}


function options_changed() {
  if (!/^\s*5\s*$/.test($("occurs").value)) return true;
  if (!/^\s*7\s*$/.test($("width").value)) return true;
  if ($("single_per_mass").checked) return true;
  if (!$("eliminate_enable").checked) return true;
  if (!/^\s*7\s*$/.test($("eliminate_width").value)) return true;
  if (!/^\s*20\s*$/.test($("count_threshold").value)) return true;
  if (!/^\s*0.000001\s*$/.test($("score_threshold").value)) return true;
  return false;
}

function options_reset(evt) {
  $("occurs").value = 5;
  $("width").value = 7;
  $("single_per_mass").checked = false;
  $("eliminate_enable").checked = true;
  $("eliminate_width").value = 7;
  $("count_threshold").value = 20;
  $("score_threshold").value = 0.000001;
}

function fix_reset() {
}

function reloadpage() {
  location.reload();
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

