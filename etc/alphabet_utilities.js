function make_alpha_bg_table(alph, freqs) {
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
  table.className = "alpha_bg_table";
  // create header
  thead = document.createElement("thead");
  table.appendChild(thead);
  row = thead.insertRow(thead.rows.length);
  if (alph.has_complement()) {
    add_text_header_cell(row, "Name", "pop_alph_name");
    if (freqs != null) add_text_header_cell(row, "Freq.", "pop_alph_freq");
    if (alph.has_bg()) add_text_header_cell(row, "Bg.", "pop_alph_bg");
    add_text_header_cell(row, "");
    add_text_header_cell(row, "");
    add_text_header_cell(row, "");
    if (alph.has_bg()) add_text_header_cell(row, "Bg.", "pop_alph_bg");
    if (freqs != null) add_text_header_cell(row, "Freq.", "pop_alph_freq");
    add_text_header_cell(row, "Name", "pop_alph_name");
  } else {
    add_text_header_cell(row, "");
    add_text_header_cell(row, "Name", "pop_alph_name");
    if (freqs != null) add_text_header_cell(row, "Freq.", "pop_alph_freq");
    if (alph.has_bg()) add_text_header_cell(row, "Bg.", "pop_alph_bg");
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
      if (freqs != null) add_text_cell(row, "" + freqs[i].toFixed(3));
      if (alph.has_bg()) add_text_cell(row, "" + alph.get_bg_freq(i).toFixed(3));
      add_cell(row, colour_symbol(i)); 
      add_text_cell(row, "~");
      add_cell(row, colour_symbol(c)); 
      if (alph.has_bg()) add_text_cell(row, "" + alph.get_bg_freq(c).toFixed(3));
      if (freqs != null) add_text_cell(row, "" + freqs[c].toFixed(3));
      add_text_cell(row, alph.get_name(c));
    }
  } else {
    for (i = 0; i < alph.get_size_core(); i++) {
      row = tbody.insertRow(tbody.rows.length);
      add_cell(row, colour_symbol(i)); 
      add_text_cell(row, alph.get_name(i));
      if (freqs != null) add_text_cell(row, "" + freqs[i].toFixed(3));
      if (alph.has_bg()) add_text_cell(row, "" + alph.get_bg_freq(i).toFixed(3));
    }
  }
  return table;
}

