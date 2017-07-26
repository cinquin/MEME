//======================================================================
// start Alphabet object
//======================================================================
var Alphabet = function(alphabet, background) {
  "use strict";
  var i, j, sym, aliases, complement, comp_e_sym, ambigs, generate_background;
  generate_background = (background == null);
  if (generate_background) {
    background = [];
    for (i = 0; i < alphabet.ncore; i++) background[i] = 1.0 / alphabet.ncore;
  } else if (alphabet.ncore != background.length) {
    throw new Error("The background length does not match the alphabet length.");
  }
  this.name = alphabet.name;
  this.like = (alphabet.like != null ? alphabet.like.toUpperCase() : null);
  this.ncore = alphabet.ncore;
  this.symbols = alphabet.symbols;
  this.background = background;
  this.genbg = generate_background;
  this.encode = {};
  this.encode2core = {};
  this.complement = {};
  // check if all symbols are same case
  var seen_uc = false;
  var seen_lc = false;
  var check_case = function (syms) {
    var s, sym;
    if (typeof syms === "string") {
      for (s = 0; s < syms.length; s++) {
        sym = syms.charAt(s);
        if (sym >= 'a' && sym <= 'z') seen_lc = true;
        else if (sym >= 'A' && sym <= 'Z') seen_uc = true;
      }
    }
  };
  for (i = 0; i < this.symbols.length; i++) {
    check_case(this.symbols[i].symbol);
    check_case(this.symbols[i].aliases);
  }
  // now map symbols to indexes
  var update_array = function(array, syms, index) {
    var s, sym;
    if (typeof syms === "string") {
      for (s = 0; s < syms.length; s++) {
        sym = syms.charAt(s);
        array[sym] = index;
        // when only a single case is used, then encode as case insensitive
        if (seen_uc != seen_lc) {
          if (sym >= 'a' && sym <= 'z') {
            array[sym.toUpperCase()] = index;
          } else if (sym >= 'A' && sym <= 'Z') {
            array[sym.toLowerCase()] = index;
          }
        }
      }
    }
  }
  // map core symbols to index
  for (i = 0; i < this.ncore; i++) {
    update_array(this.encode2core, this.symbols[i].symbol, i);
    update_array(this.encode, this.symbols[i].symbol, i);
    update_array(this.encode2core, this.symbols[i].aliases, i);
    update_array(this.encode, this.symbols[i].aliases, i);
  }
  // map ambigous symbols to index
  ambigs = {};
  for (i = this.ncore; i < this.symbols.length; i++) {
    update_array(this.encode, this.symbols[i].symbol, i);
    update_array(this.encode, this.symbols[i].aliases, i);
    ambigs[this.symbols[i].equals] = i;
  }
  // determine complements
  for (i = 0; i < this.ncore; i++) {
    complement = this.symbols[i].complement;
    if (typeof complement === "string") {
      this.complement[i] = this.encode2core[complement];
    }
  }
  next_symbol:
  for (i = this.ncore; i < this.symbols.length; i++) {
    complement = "";
    for (j = 0; j < this.symbols[i].equals.length; j++) {
      comp_e_sym = this.complement[this.encode2core[this.symbols[i].equals.charAt(j)]];
      if (typeof comp_e_sym !== "number") continue next_symbol;
      complement += this.symbols[comp_e_sym].symbol;
    }
    complement = complement.split("").sort().join("");
    if (typeof ambigs[complement] === "number") {
      this.complement[i] = ambigs[complement];
    }
  }
  // determine case insensitivity
  this.case_insensitive = true;
  if (seen_uc == seen_lc) {
    // when there is a mixture of cases it probably won't
    // be case insensitive but we still need to check
    loop:
    for (i = 0; i < this.symbols.length; i++) {
      sym = this.symbols[i].symbol;
      if (sym >= 'A' && sym <= 'Z') {
        if (this.encode[sym.toLowerCase()] != i) {
          this.case_insensitive = false;
          break loop;
        }
      } else if (sym >= 'a' && sym <= 'z') {
        if (this.encode[sym.toUpperCase()] != i) {
          this.case_insensitive = false;
          break loop;
        }
      }
      aliases = this.symbols[i].aliases;
      if (aliases != null) {
        for (j = 0; j < aliases.length; j++) {
          sym = aliases.charAt(j);
          if (sym >= 'A' && sym <= 'Z') {
            if (this.encode[sym.toLowerCase()] != i) {
              this.case_insensitive = false;
              break loop;
            }
          } else if (sym >= 'a' && sym <= 'z') {
            if (this.encode[sym.toUpperCase()] != i) {
              this.case_insensitive = false;
              break loop;
            }
          }
        }
      }
    }
  }
  // normalise aliases to remove the prime symbol and eliminate
  // the alternate cases when the alphabet is case insensitive
  var seen, out;
  for (i = 0; i < this.symbols.length; i++) {
    sym = this.symbols[i].symbol;
    aliases = this.symbols[i].aliases;
    if (typeof aliases != "string") aliases = "";
    seen = {};
    out = [];
    if (this.case_insensitive) {
      sym = sym.toUpperCase();
      aliases = aliases.toUpperCase();
    }
    seen[sym] = true;
    for (j = 0; j < aliases.length; j++) {
      if (!seen[aliases.charAt(j)]) {
        seen[aliases.charAt(j)] = true;
        out.push(aliases.charAt(j));
      }
    }
    this.symbols[i].aliases = out.sort().join("");
  }
};
// return the name of the alphabet
Alphabet.prototype.get_alphabet_name = function() {
  return this.name;
};
// return if the alphabet can be complemented
Alphabet.prototype.has_complement = function() {
  return (typeof this.symbols[0].complement === "string");
};
// return true if an uppercase letter has the same meaning as the lowercase form
Alphabet.prototype.is_case_insensitive = function() {
  return this.case_insensitive;
};
// return the information content of an alphabet letter
Alphabet.prototype.get_ic = function() {
  return Math.log(this.ncore) / Math.LN2;
};
// return the count of the core alphabet symbols
Alphabet.prototype.get_size_core = function() {
  return this.ncore;
};
// return the count of all alphabet symbols
Alphabet.prototype.get_size_full = function() {
  return this.symbols.length;
};
// return the symbol for the given alphabet index
Alphabet.prototype.get_symbol = function(alph_index) {
  "use strict";
  if (alph_index < 0 || alph_index >= this.symbols.length) {
    throw new Error("Alphabet index out of bounds");
  }
  return this.symbols[alph_index].symbol;
};
// return the aliases for the given alphabet index
Alphabet.prototype.get_aliases = function(alph_index) {
  "use strict";
  if (alph_index < 0 || alph_index >= this.symbols.length) {
    throw new Error("Alphabet index out of bounds");
  }
  var sym_obj = this.symbols[alph_index];
  return (sym_obj.aliases != null ? sym_obj.aliases : "");
};
// return the name for the given alphabet index
Alphabet.prototype.get_name = function(alph_index) {
  "use strict";
  var sym;
  if (alph_index < 0 || alph_index >= this.symbols.length) {
    throw new Error("Alphabet index out of bounds");
  }
  sym = this.symbols[alph_index];
  return (typeof sym.name === "string" ? sym.name : sym.symbol);
};
// return the alphabet it is like or null
Alphabet.prototype.get_like = function() {
  "use strict";
  return this.like;
};
// return the index of the complement for the given alphabet index
Alphabet.prototype.get_complement = function(alph_index) {
  var comp_e_sym = this.complement[alph_index];
  if (typeof comp_e_sym === "number") {
    return comp_e_sym;
  } else {
    return -1;
  }
};
// return a string containing the core symbols
Alphabet.prototype.get_symbols = function() {
  "use strict";
  var i, core_symbols;
  core_symbols = "";
  for (i = 0; i < this.ncore; i++) {
    core_symbols += this.symbols[i].symbol;
  }
  return core_symbols;
};
// return if the background was not a uniform generated background
Alphabet.prototype.has_bg = function() {
  "use strict";
  return !this.genbg;
};
// get the background frequency for the index
Alphabet.prototype.get_bg_freq = function(alph_index) {
  "use strict";
  var freq, i, symbols;
  if (alph_index >= 0) {
    if (alph_index < this.ncore) {
      return this.background[alph_index];
    } else if (alph_index < this.symbols.length) {
      freq = 0;
      symbols = this.symbols[alph_index].equals;
      for (i = 0; i < symbols.length; i++) {
        freq += this.background[this.encode2core[symbols.charAt(i)]];
      }
      return freq;
    } 
  }
  throw new Error("The alphabet index is out of range.");
};
// get the colour of the index
Alphabet.prototype.get_colour = function(alph_index) {
  "use strict";
  if (alph_index < 0 || alph_index >= this.symbols.length) {
    throw new Error("BAD_ALPHABET_INDEX");
  }
  if (typeof this.symbols[alph_index].colour != "string") {
    return "black";
  }
  return "#" + this.symbols[alph_index].colour;
};
// get the rgb componets of the colour at the index
Alphabet.prototype.get_rgb = function(alph_index) {
  "use strict";
  if (alph_index < 0 || alph_index >= this.symbols.length) {
    throw new Error("BAD_ALPHABET_INDEX");
  }
  if (typeof this.symbols[alph_index].colour != "string") {
    return {"red": 0, "green": 0, "blue": 0};
  }
  var colour = this.symbols[alph_index].colour;
  var red = parseInt(colour.substr(0, 2), 16) / 255;
  var green = parseInt(colour.substr(2, 2), 16) / 255;
  var blue = parseInt(colour.substr(4, 2), 16) / 255;
  return {"red": red, "green": green, "blue": blue};
};
// convert a symbol into the index
Alphabet.prototype.get_index = function(letter) {
  "use strict";
  var alph_index;
  alph_index = this.encode[letter];
  if (typeof alph_index === "undefined") {
    return -1;
  }
  return alph_index;
};
// convert a symbol into the list of core indexes that it equals
Alphabet.prototype.get_indexes = function(letter) {
  "use strict";
  var alph_index, comprise_str, i, comprise_list;
  alph_index = this.encode[letter];
  if (typeof alph_index === "undefined") {
    throw new Error("Unknown letter");
  }
  comprise_str = this.symbols[alph_index].equals;
  comprise_list = [];
  if (typeof comprise_str == "string") {
    for (i = 0; i < comprise_str.length; i++) {
      comprise_list.push(this.encode2core[comprise_str.charAt(i)]);
    }
  } else {
    comprise_list.push(alph_index);
  }
  return comprise_list;
};
// check if a symbol is the primary way of representing the symbol in the alphabet
Alphabet.prototype.is_prime_symbol = function(letter) {
  var alph_index;
  alph_index = this.encode[letter];
  if (alph_index == null) return false;
  if (this.is_case_insensitive()) {
    return (this.symbols[alph_index].symbol.toUpperCase() == letter.toUpperCase());
  } else {
    return (this.symbols[alph_index].symbol == letter);
  }
};
// compare 2 alphabets
Alphabet.prototype.equals = function(other) {
  "use strict";
  var i, sym1, sym2;
  // first check that it's actually an alphabet object
  if (!(typeof other === "object" && other != null && other instanceof Alphabet)) {
    return false;
  }
  // second shortcircuit if it's the same object
  if (this === other) return true;
  // compare
  if (this.name !== other.name) return false;
  if (this.ncore !== other.ncore) return false;
  if (this.symbols.length !== other.symbols.length) return false;
  for (i = 0; i < this.symbols.length; i++) {
    sym1 = this.symbols[i];
    sym2 = other.symbols[i];
    if (sym1.symbol !== sym2.symbol) return false;
    if (sym1.aliases !== sym2.aliases) return false;
    if (sym1.name !== sym2.name) return false;
    if (typeof sym1.colour !== typeof sym2.colour || 
        (typeof sym1.colour === "string" && typeof sym2.colour === "string" &&
         parseInt(sym1.colour, 16) != parseInt(sym2.colour, 16))) {
      return false;
    }
    if (sym1.complement !== sym2.complement) return false;
    if (sym1.equals !== sym2.equals) return false;
  }
  return true;
};
Alphabet.prototype.check_core_subset = function(super_alph) {
  var complement_same = true;
  var seen_set = {};
  var sub_i, sub_symbol, super_i, super_symbol;
  for (sub_i = 0; sub_i < this.ncore; sub_i++) {
    sub_symbol = this.symbols[sub_i];
    super_i = super_alph.encode[sub_symbol.symbol]; 
    if (super_i == null) return 0;
    super_symbol = super_alph.symbols[super_i];
    if (seen_set[super_i]) return 0;
    seen_set[super_i] = true;
    // check complement
    if (sub_symbol.complement != null && super_symbol.complement != null) {
      if (super_alph.encode[sub_symbol.complement] != super_alph.encode[super_symbol.complement]) {
        complement_same = false;
      }
    } else if (sub_symbol.complement != null || super_symbol.complement != null) {
      complement_same = false;
    }
  }
  return (complement_same ? 1 : -1);
};
// convert a sequence to its reverse complement
Alphabet.prototype.invcomp_seq = function(seq) {
  "use strict";
  var syms, i, e_sym, comp_e_sym;
  if (!this.has_complement()) throw new Error("Alphabet must be complementable");
  syms = seq.split("");
  for (i = 0; i < syms.length; i++) {
    e_sym = this.encode[syms[i]];
    if (typeof e_sym === "undefined") {
      e_sym = this.ncore; // wildcard
    }
    comp_e_sym = this.complement[e_sym];
    if (typeof comp_e_sym === "undefined") {
      comp_e_sym = e_sym; // not complementable
    }
    syms[i] = this.symbols[comp_e_sym].symbol;
  }
  return syms.reverse().join("");
};
// convert the alphabet to the text version
Alphabet.prototype.as_text = function() {
  "use strict";
  function name_as_text(name) {
    var i, c, out;
    out = "\"";
    for (i = 0; i < name.length; i++) {
      c = name.charAt(i);
      if (c == "\"") {
        out += "\\\"";
      } else if (c == "/") {
        out += "\\/";
      } else if (c == "\\") {
        out += "\\\\";
      } else {
        out += c;
      }
    }
    out += "\"";
    return out;
  }
  function symbol_as_text(sym) {
    var out;
    out = sym.symbol;
    if (typeof sym.name === "string" && sym.name != sym.symbol) {
      out += " " + name_as_text(sym.name);
    }
    if (typeof sym.colour === "string") {
      out += " " + sym.colour;
    }
    return out;
  }
  var out, i, j, c, sym;
  out = "";
  // output core symbols with 2 way complements
  for (i = 0; i < this.ncore; i++) {
    c = this.complement[i];
    if (typeof c === "number" && i < c && this.complement[c] === i) {
      out += symbol_as_text(this.symbols[i]) + " ~ " + symbol_as_text(this.symbols[c]) + "\n";  
    }
  }
  // output core symbols with no complement
  for (i = 0; i < this.ncore; i++) {
    if (typeof this.complement[i] === "undefined") {
      out += symbol_as_text(this.symbols[i]) + "\n";
    }
  }
  // output ambiguous symbols that have comprising characters
  for (i = this.ncore; i < this.symbols.length; i++) {
    if (this.symbols[i].equals.length == 0) break;
    out += symbol_as_text(this.symbols[i]) + " = " + this.symbols[i].equals + "\n";
    if (typeof this.symbols[i].aliases === "string") {
      for (j = 0; j < this.symbols[i].aliases.length; j++) {
        if (this.symbols[i].aliases.charAt(j) == this.symbols[i].symbol) continue;
        out += this.symbols[i].aliases.charAt(j) + " = " + this.symbols[i].equals + "\n";
      }
    }
  }
  // output aliases of core symbols
  for (i = 0; i < this.ncore; i++) {
    if (typeof this.symbols[i].aliases === "string") {
      for (j = 0; j < this.symbols[i].aliases.length; j++) {
        if (this.symbols[i].aliases.charAt(j) == this.symbols[i].symbol) continue;
        out += this.symbols[i].aliases.charAt(j) + " = " + this.symbols[i].symbol + "\n";
      }
    }
  }
  // output gap symbols
  i = this.symbols.length - 1;
  if (this.symbols[i].equals.length == 0) {
    out += symbol_as_text(this.symbols[i]) + " =\n";
    if (typeof this.symbols[i].aliases === "string") {
      for (j = 0; j < this.symbols[i].aliases.length; j++) {
        if (this.symbols[i].aliases.charAt(j) == this.symbols[i].symbol) continue;
        out += this.symbols[i].aliases.charAt(j) + " =\n";
      }
    }
  }
  return out;
};
// output the alphabet as it appears in minimal MEME format
Alphabet.prototype.as_meme = function() {
  "use strict";
  function name_as_text(name) {
    var i, c, out;
    out = "\"";
    for (i = 0; i < name.length; i++) {
      c = name.charAt(i);
      if (c == "\"") {
        out += "\\\"";
      } else if (c == "/") {
        out += "\\/";
      } else if (c == "\\") {
        out += "\\\\";
      } else {
        out += c;
      }
    }
    out += "\"";
    return out;
  }
  if (this.equals(AlphStd.DNA)) {
    return "ALPHABET= ACGT\n";
  } else if (this.equals(AlphStd.PROTEIN)) {
    return "ALPHABET= ACDEFGHIKLMNPQRSTVWY\n";
  } else {
    return "ALPHABET" + 
      (this.name != null ? " " + name_as_text(this.name) : "") + 
      (this.like != null ? " " + this.like + "-LIKE" : "") + "\n" +
      this.as_text() + "END ALPHABET\n";
  }
};

// Returns a table showing all the letters in the alphabet
Alphabet.prototype.as_table = function() {
  "use strict";
  var i, j, row, th, td, aliases, equals, sym;
  var table = document.createElement("table");
  // create the core symbol header
  row = table.insertRow(table.rows.length);
  th = document.createElement("th");
  th.appendChild(document.createTextNode("Symbol(s)"));
  row.appendChild(th);
  th = document.createElement("th");
  th.appendChild(document.createTextNode("Name"));
  row.appendChild(th);
  th = document.createElement("th");
  if (this.has_complement()) {
    th.appendChild(document.createTextNode("Complement"));
  }
  row.appendChild(th);
  // list the core symbols
  for (i = 0; i < this.ncore; i++) {
    row = table.insertRow(table.rows.length);
    td = document.createElement("td");
    if (this.symbols[i].colour != null) {
      td.style.color = '#' + this.symbols[i].colour;
    }
    td.appendChild(document.createTextNode(this.symbols[i].symbol));
    aliases = this.get_aliases(i);
    if (aliases.length > 0) {
      td.appendChild(document.createTextNode(' ' + aliases.split('').join(' ')));
    }
    row.appendChild(td);
    td = document.createElement("td");
    if (this.symbols[i].name != null) {
      td.appendChild(document.createTextNode(this.symbols[i].name));
    }
    row.appendChild(td);
    td = document.createElement("td");
    if (this.symbols[i].complement != null) {
      td.style.color = this.get_colour(this.get_index(this.symbols[i].complement));
      td.appendChild(document.createTextNode(this.symbols[i].complement));
    }
    row.appendChild(td);
  }
  // create the ambiguous symbol header
  row = table.insertRow(table.rows.length);
  th = document.createElement("th");
  th.appendChild(document.createTextNode("Symbol(s)"));
  row.appendChild(th);
  th = document.createElement("th");
  th.appendChild(document.createTextNode("Name"));
  row.appendChild(th);
  th = document.createElement("th");
  th.appendChild(document.createTextNode("Matches"));
  row.appendChild(th);
  // list the ambiguous symbols
  for (i = this.ncore; i < this.symbols.length; i++) {
    row = table.insertRow(table.rows.length);
    td = document.createElement("td");
    if (this.symbols[i].colour != null) {
      td.style.color = '#' + this.symbols[i].colour;
    }
    td.appendChild(document.createTextNode(this.symbols[i].symbol));
    aliases = this.get_aliases(i);
    if (aliases.length > 0) {
      td.appendChild(document.createTextNode(' ' + aliases.split('').join(' ')));
    }
    row.appendChild(td);
    td = document.createElement("td");
    if (this.symbols[i].name != null) {
      td.appendChild(document.createTextNode(this.symbols[i].name));
    }
    row.appendChild(td);
    td = document.createElement("td");
    equals = this.symbols[i].equals.split('');
    for (j = 0; j < equals.length; j++) {
      if (j != 0) td.appendChild(document.createTextNode(' '));
      sym = document.createElement("span");
      sym.style.color = this.get_colour(this.get_index(equals[j]));
      sym.appendChild(document.createTextNode(equals[j]));
      td.appendChild(sym);
    }
    row.appendChild(td);
  }
  return table;
};

// returns a dictionary of the colours for EPS
Alphabet.prototype._as_eps_dict = function() {
  "use strict";
  var i, sym, rgb;
  var out = "/fullColourDict <<\n";
  for (i = 0; i < this.ncore; i++) {
    sym = this.get_symbol(i);
    sym = sym.replace(/\\/g, "\\\\");
    sym = sym.replace(/\(/g, "\\(");
    sym = sym.replace(/\)/g, "\\)");
    rgb = this.get_rgb(i);
    out += " (" + sym + ") [" + rgb.red.toFixed(4) + " " + rgb.green.toFixed(4) + " " + rgb.blue.toFixed(4) + "]\n";
  }
  out += ">> def\n";
  out += "/mutedColourDict <<\n";
  for (i = 0; i < this.ncore; i++) {
    sym = this.get_symbol(i);
    sym = sym.replace(/\\/g, "\\\\");
    sym = sym.replace(/\(/g, "\\(");
    sym = sym.replace(/\)/g, "\\)");
    rgb = Alphabet.lighten_colour(this.get_rgb(i));
    out += " (" + sym + ") [" + rgb.red.toFixed(4) + " " + rgb.green.toFixed(4) + " " + rgb.blue.toFixed(4) + "]\n";
  }
  out += ">> def\n";
  return out;
};

// return the alphabet name or a list of primary symbols
Alphabet.prototype.toString = function() {
  "use strict";
  if (this.name != null) {
    return this.name;
  } else {
    return this.get_symbols();
  }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Helper functions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// Convert a colour specified in RGB colourspace values into LAB colourspace
Alphabet.rgb2lab = function(rgb) {
  "use strict";
  var xyzHelper, labHelper;
  // XYZ helper
  xyzHelper = function(value) {
    if (value > 0.0445) {
      value = (value + 0.055) / 1.055;
      value = Math.pow(value, 2.4);
    } else {
      value /= 12.92;
    }
    value *= 100;
    return value;
  };
  // lab helper
  labHelper = function(value) {
    if (value > 0.008856) {
      value = Math.pow(value, 1.0 / 3.0);
    } else {
      value = (7.787 * value) + (16.0 / 116.0);
    }
    return value;
  };
  // convert into XYZ colourspace
  var c1, c2, c3;
  if (typeof rgb == "number") {
    c1 = xyzHelper(((rgb >> 16) & 0xFF) / 255.0);
    c2 = xyzHelper(((rgb >> 8) & 0xFF) / 255.0);
    c3 = xyzHelper((rgb & 0xFF) / 255.0);
  } else {
    c1 = xyzHelper(rgb.red);
    c2 = xyzHelper(rgb.green);
    c3 = xyzHelper(rgb.blue);
  }
  var x = (c1 * 0.4124) + (c2 * 0.3576) + (c3 * 0.1805);
  var y = (c1 * 0.2126) + (c2 * 0.7152) + (c3 * 0.0722);
  var z = (c1 * 0.0193) + (c2 * 0.1192) + (c3 * 0.9505);
  // convert into Lab colourspace
  c1 = labHelper(x / 95.047);
  c2 = labHelper(y / 100.0);
  c3 = labHelper(z / 108.883);
  var l = (116.0 * c2) - 16;
  var a = 500.0 * (c1 - c2);
  var b = 200.0 * (c2 - c3);
  return {"l": l, "a": a, "b": b};
};

// Convert a colour specified in HSV colourspace into RGB colourspace
Alphabet.hsv2rgb = function(hue, sat, value, output_object) {
  // achromatic (grey)
  var r = value;
  var g = value;
  var b = value;
  if (sat != 0) {
    var h = hue / 60.0;
    var i = Math.floor(h);
    var f = h - i;
    var p = value * (1.0 - sat);
    var q = value * (1.0 - (sat * f));
    var t = value * (1.0 - (sat * (1.0 - f)));
    if (i == 0) {
      r = value;
      g = t;
      b = p;
    } else if (i == 1) {
      r = q;
      g = value;
      b = p;
    } else if (i == 2) {
      r = p;
      g = value;
      b = t;
    } else if (i == 3) {
      r = p;
      g = q;
      b = value;
    } else if (i == 4) {
      r = t;
      g = p;
      b = value;
    } else {
      r = value;
      g = p;
      b = q;
    }
  }
  if (output_object) {
    return {"red": r, "green": g, "blue": b};
  } else {
    return (Math.floor(r * 255) << 15) | (Math.floor(g * 255) << 8) | (Math.floor(b * 255));
  }
};

// Calculate a distance score between two colours in LAB colourspace
Alphabet.lab_dist = function(lab1, lab2) {
  var c1 = Math.sqrt((lab1.l * lab1.l) + (lab1.a * lab1.a));
  var c2 = Math.sqrt((lab2.l * lab2.l) + (lab2.a * lab2.a));
  var dc = c1 - c2;
  var dl = lab1.l - lab2.l;
  var da = lab1.a - lab2.a;
  var db = lab1.b - lab2.b;
  // we don't want NaN due to rounding errors so fudge things a bit...
  var dh = 0;
  var dh_squared = (da * da) + (db * db) - (dc * dc);
  if (dh_squared > 0) {
    dh = Math.sqrt(dh_squared);
  }
  var first = dl;
  var second = dc / (1.0 + (0.045 * c1));
  var third = dh / (1.0 + (0.015 * c1));
  return Math.sqrt((first * first) + (second * second) + (third * third));
};

// convert an RGB value into a HSL value
Alphabet.rgb2hsl = function(rgb) {
  "use strict";
  var min, max, delta, h, s, l, r, g, b;
  if (typeof rgb == "number") {
    r = ((rgb >> 16) & 0xFF) / 255.0;
    g = ((rgb >> 8) & 0xFF) / 255.0;
    b = (rgb & 0xFF) / 255.0;
  } else {
    r = rgb.red;
    g = rgb.green;
    b = rgb.blue;
  }
  min = Math.min(r, g, b);
  max = Math.max(r, g, b);
  delta = max - min;
  l = min + (delta / 2);
  if (max == min) {
    h = 0; // achromatic (grayscale)
    s = 0;
  } else {
    if (l > 0.5) {
      s = delta / (2 - max - min);
    } else {
      s = delta / (max + min);
    }
    if (max == r) {
      h = (g - b) / delta;
      if (g < b) h += 6;
    } else if (max == g) {
      h = ((b - r) / delta) + 2;
    } else {
      h = ((r - g) / delta) + 4;
    }
    h /= 6;
  }
  return {"h": h, "s": s, "l": l};
};

// convert a HSL value into an RGB value
Alphabet.hsl2rgb = function(hsl, output_object) {
  "use strict";
  function _hue(p, q, t) {
    "use strict";
    if (t < 0) t += 1;
    else if (t > 1) t -= 1;
    if (t < (1.0 / 6.0)) {
      return p + ((q - p) * 6.0 * t);
    } else if (t < 0.5) {
      return q;
    } else if (t < (2.0 / 3.0)) {
      return p + ((q - p) * ((2.0 / 3.0) - t) * 6.0);
    } else {
      return p;
    }
  }
  var r, g, b, p, q;
  if (hsl.s == 0) {
    // achromatic (grayscale)
    r = hsl.l;
    g = hsl.l;
    b = hsl.l;
  } else {
    if (hsl.l < 0.5) {
      q = hsl.l * (1 + hsl.s);
    } else {
      q = hsl.l + hsl.s - (hsl.l * hsl.s);
    }
    p = (2 * hsl.l) - q;
    r = _hue(p, q, hsl.h + (1.0 / 3.0));
    g = _hue(p, q, hsl.h);
    b = _hue(p, q, hsl.h - (1.0 / 3.0));
  }
  if (output_object) {
    return {"red": r, "green": g, "blue": b};
  } else {
    return (Math.floor(r * 255) << 15) | (Math.floor(g * 255) << 8) | (Math.floor(b * 255));
  }
};

Alphabet.lighten_colour = function(rgb) {
  "use strict";
  var hsl = Alphabet.rgb2hsl(rgb);
  hsl.l += (1.0 - hsl.l) * 2 / 3;
  return Alphabet.hsl2rgb(hsl, typeof rgb != "number");
};

//======================================================================
// end Alphabet object
//======================================================================

//======================================================================
// start StandardAlphabet object
//======================================================================

// an extension of the alphabet object to support some additional fields 
// only present in standard alphabets.
var StandardAlphabet = function(enum_code, enum_name, alphabet_data) {
  Alphabet.apply(this, [alphabet_data]);
  this.enum_code = enum_code;
  this.enum_name = enum_name;
};
StandardAlphabet.prototype = Alphabet.prototype;
StandardAlphabet.prototype.constructor = StandardAlphabet;

// A unique code for this standard alphabet.
// This code will be a power of 2 to enable creation of bitsets for
// a selection of standard alphabets.
StandardAlphabet.prototype.get_code = function() {
  return this.enum_code;
};

// A unique name for this standard alphabet.
// this name will be all upper case and the same as the property that
// refers to this alphabet in the AlphStd collection.
StandardAlphabet.prototype.get_enum = function() {
  return this.enum_name;
};

//======================================================================
// end StandardAlphabet object
//======================================================================

// A collection of standard alphabets.
var AlphStd = {
  RNA: new StandardAlphabet(1, "RNA", {
    "name": "RNA",
    "like": "RNA",
    "ncore": 4,
    "symbols": [
      {"symbol": "A", "name": "Adenine", "colour": "CC0000"},
      {"symbol": "C", "name": "Cytosine", "colour": "0000CC"},
      {"symbol": "G", "name": "Guanine", "colour": "FFB300"},
      {"symbol": "U", "name": "Uracil", "colour": "008000",
        "aliases": "T"},
      {"symbol": "N", "name": "Any base", "equals": "ACGU", "aliases": "X."},
      {"symbol": "V", "name": "Not U", "equals": "ACG"},
      {"symbol": "H", "name": "Not G", "equals": "ACU"},
      {"symbol": "D", "name": "Not C", "equals": "AGU"},
      {"symbol": "B", "name": "Not A", "equals": "CGU"},
      {"symbol": "M", "name": "Amino", "equals": "AC"},
      {"symbol": "R", "name": "Purine", "equals": "AG"},
      {"symbol": "W", "name": "Weak", "equals": "AU"}, 
      {"symbol": "S", "name": "Strong", "equals": "CG"},
      {"symbol": "Y", "name": "Pyrimidine", "equals": "CU"},
      {"symbol": "K", "name": "Keto", "equals": "GU"}
    ]
  }), 
  DNA: new StandardAlphabet(2, "DNA", {
    "name": "DNA",
    "like": "DNA",
    "ncore": 4,
    "symbols": [
      {"symbol": "A", "name": "Adenine", "colour": "CC0000", "complement": "T"},
      {"symbol": "C", "name": "Cytosine", "colour": "0000CC", "complement": "G"},
      {"symbol": "G", "name": "Guanine", "colour": "FFB300", "complement": "C"},
      {"symbol": "T", "name": "Thymine", "colour": "008000", "complement": "A",
        "aliases": "U"},
      {"symbol": "N", "name": "Any base", "equals": "ACGT", "aliases": "X."},
      {"symbol": "V", "name": "Not T", "equals": "ACG"},
      {"symbol": "H", "name": "Not G", "equals": "ACT"},
      {"symbol": "D", "name": "Not C", "equals": "AGT"},
      {"symbol": "B", "name": "Not A", "equals": "CGT"},
      {"symbol": "M", "name": "Amino", "equals": "AC"},
      {"symbol": "R", "name": "Purine", "equals": "AG"},
      {"symbol": "W", "name": "Weak", "equals": "AT"}, 
      {"symbol": "S", "name": "Strong", "equals": "CG"},
      {"symbol": "Y", "name": "Pyrimidine", "equals": "CT"},
      {"symbol": "K", "name": "Keto", "equals": "GT"}
    ]
  }), 
  PROTEIN: new StandardAlphabet(4, "PROTEIN", {
    "name": "Protein",
    "like": "PROTEIN",
    "ncore": 20,
    "symbols": [
      {"symbol": "A", "name": "Alanine", "colour": "0000CC"},
      {"symbol": "C", "name": "Cysteine", "colour": "0000CC"},
      {"symbol": "D", "name": "Aspartic acid", "colour": "FF00FF"},
      {"symbol": "E", "name": "Glutamic acid", "colour": "FF00FF"},
      {"symbol": "F", "name": "Phenylalanine", "colour": "0000CC"},
      {"symbol": "G", "name": "Glycine", "colour": "FFB300"},
      {"symbol": "H", "name": "Histidine", "colour": "FFCCCC"},
      {"symbol": "I", "name": "Isoleucine", "colour": "0000CC"},
      {"symbol": "K", "name": "Lysine", "colour": "CC0000"},
      {"symbol": "L", "name": "Leucine", "colour": "0000CC"},
      {"symbol": "M", "name": "Methionine", "colour": "0000CC"},
      {"symbol": "N", "name": "Asparagine", "colour": "008000"},
      {"symbol": "P", "name": "Proline", "colour": "FFFF00"},
      {"symbol": "Q", "name": "Glutamine", "colour": "008000"},
      {"symbol": "R", "name": "Arginine", "colour": "CC0000"},
      {"symbol": "S", "name": "Serine", "colour": "008000"},
      {"symbol": "T", "name": "Threonine", "colour": "008000"},
      {"symbol": "V", "name": "Valine", "colour": "0000CC"},
      {"symbol": "W", "name": "Tryptophan", "colour": "0000CC"},
      {"symbol": "Y", "name": "Tyrosine", "colour": "33E6CC"},
      {"symbol": "X", "name": "Any amino acid", "equals": "ACDEFGHIKLMNPQRSTVWY", "aliases": "*."},
      {"symbol": "B", "name": "Asparagine or Aspartic acid", "equals": "DN"}, 
      {"symbol": "Z", "name": "Glutamine or Glutamic acid", "equals": "EQ"}, 
      {"symbol": "J", "name": "Leucine or Isoleucine", "equals": "IL"}
    ]
  })
};

//======================================================================
// start Symbol object
//======================================================================
var Symbol = function(alph_index, scale, alphabet) {
  "use strict";
  //variable prototype
  this.symbol = alphabet.get_symbol(alph_index);
  this.scale = scale;
  this.colour = alphabet.get_colour(alph_index);
};

Symbol.prototype.get_symbol = function() {
  "use strict";
  return this.symbol;
};

Symbol.prototype.get_scale = function() {
  "use strict";
  return this.scale;
};

Symbol.prototype.get_colour = function() {
  "use strict";
  return this.colour;
};

Symbol.prototype.toString = function() {
  "use strict";
  return this.symbol + " " + (Math.round(this.scale*1000)/10) + "%";
};

function compare_symbol(sym1, sym2) {
  "use strict";
  if (sym1.get_scale() < sym2.get_scale()) {
    return -1;
  } else if (sym1.get_scale() > sym2.get_scale()) {
    return 1;
  } else {
    return 0;
  }
}
//======================================================================
// end Symbol object
//======================================================================

//======================================================================
// start Pspm object
//======================================================================
var Pspm = function(matrix, name, ltrim, rtrim, nsites, evalue, pssm, alt) {
  "use strict";
  var row, col, data, row_sum, delta, evalue_re;
  if (typeof name !== "string") {
    name = "";
  }
  this.name = name;
  //construct
  if (matrix instanceof Pspm) {
    // copy constructor
    this.alph_length = matrix.alph_length;
    this.motif_length = matrix.motif_length;
    this.name = matrix.name;
    this.alt = matrix.alt;
    this.nsites = matrix.nsites;
    this.evalue = matrix.evalue;
    this.ltrim = matrix.ltrim;
    this.rtrim = matrix.rtrim;
    this.pspm = [];
    for (row = 0; row < matrix.motif_length; row++) {
      this.pspm[row] = [];
      for (col = 0; col < matrix.alph_length; col++) {
        this.pspm[row][col] = matrix.pspm[row][col];
      }
    }
    if (matrix.pssm != null) {
      this.pssm = [];
      for (row = 0; row < matrix.motif_length; row++) {
        this.pspm[row] = [];
        for (col = 0; col < matrix.alph_length; col++) {
          this.pssm[row][col] = matrix.pssm[row][col];
        }
      }
    }
  } else {
    // check parameters
    if (ltrim == null) {
      ltrim = 0;
    } else if (typeof ltrim !== "number" || ltrim % 1 !== 0 || ltrim < 0) {
      throw new Error("ltrim must be a non-negative integer, got: " + ltrim);
    }
    if (rtrim == null) {
      rtrim = 0;
    } else if (typeof rtrim !== "number" || rtrim % 1 !== 0 || rtrim < 0) {
      throw new Error("rtrim must be a non-negative integer, got: " + rtrim);
    }
    if (nsites != null) {
      if (typeof nsites !== "number" || nsites < 0) {
        throw new Error("nsites must be a positive number, got: " + nsites);
      } else if (nsites == 0) {
        nsites = null;
      }
    }
    if (evalue != null) {
      if (typeof evalue === "number") {
        if (evalue < 0) {
          throw new Error("evalue must be a non-negative number, got: " + evalue);
        }
      } else if (typeof evalue === "string") {
        evalue_re = /^((?:[+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|inf)$/;
        if (!evalue_re.test(evalue)) {
          throw new Error("evalue must be a non-negative number, got: " + evalue);
        }
      } else {
        throw new Error("evalue must be a non-negative number, got: " + evalue);
      }
    }
    // set properties
    this.name = name;
    this.alt = alt;
    this.nsites = nsites;
    this.evalue = evalue;
    this.ltrim = ltrim;
    this.rtrim = rtrim;
    if (typeof matrix === "string") {
      // string constructor
      data = parse_pspm_string(matrix);
      this.alph_length = data["alph_length"];
      this.motif_length = data["motif_length"];
      this.pspm = data["pspm"];
      if (this.evalue == null) {
        if (data["evalue"] != null) {
          this.evalue = data["evalue"];
        } else {
          this.evalue = 0;
        }
      }
      if (this.nsites == null) {
        if (typeof data["nsites"] === "number") {
          this.nsites = data["nsites"];
        } else {
          this.nsites = 20;
        }
      }
    } else {
      // assume pspm is a nested array
      this.motif_length = matrix.length;
      this.alph_length = (matrix.length > 0 ? matrix[0].length : 0);
      if (this.nsites == null) {
        this.nsites = 20;
      }
      if (this.evalue == null) {
        this.evalue = 0;
      }
      this.pspm = [];
      // copy pspm and check
      for (row = 0; row < this.motif_length; row++) {
        if (this.alph_length != matrix[row].length) {
          throw new Error("COLUMN_MISMATCH");
        }
        this.pspm[row] = [];
        row_sum = 0;
        for (col = 0; col < this.alph_length; col++) {
          this.pspm[row][col] = matrix[row][col];
          row_sum += this.pspm[row][col];
        }
        delta = 0.1;
        if (isNaN(row_sum) || (row_sum > 1 && (row_sum - 1) > delta) || 
            (row_sum < 1 && (1 - row_sum) > delta)) {
          throw new Error("INVALID_SUM");
        }
      }
      // copy pssm
      if (pssm != null) {
        this.pssm = [];
        for (row = 0; row < this.motif_length; row++) {
          this.pssm[row] = [];
          for (col = 0; col < this.alph_length; col++) {
            this.pssm[row][col] = pssm[row][col];
          }
        }
      }
    }
  }
};

Pspm.prototype.copy = function() {
  "use strict";
  return new Pspm(this);
};

Pspm.prototype.reverse = function() {
  "use strict";
  var x, y, temp, temp_trim;
  //reverse
  x = 0;
  y = this.motif_length-1;
  while (x < y) {
    temp = this.pspm[x];
    this.pspm[x] = this.pspm[y];
    this.pspm[y] = temp;
    x++;
    y--;
  }
  // reverse pssm (if defined)
  if (typeof this.pssm !== "undefined") {
    //reverse
    x = 0;
    y = this.motif_length-1;
    while (x < y) {
      temp = this.pssm[x];
      this.pspm[x] = this.pssm[y];
      this.pssm[y] = temp;
      x++;
      y--;
    }
  }
  //swap triming
  temp_trim = this.ltrim;
  this.ltrim = this.rtrim;
  this.rtrim = temp_trim;
  return this; //allow function chaining...
};

Pspm.prototype.reverse_complement = function(alphabet) {
  "use strict";
  var x, y, temp, i, row, c, temp_trim;
  if (this.alph_length != alphabet.get_size_core()) {
    throw new Error("The alphabet size does not match the size of the pspm.");
  }
  if (!alphabet.has_complement()) {
    throw new Error("The specified alphabet can not be complemented.");
  }
  // reverse motif
  this.reverse();
  //complement
  for (x = 0; x < this.motif_length; x++) {
    row = this.pspm[x];
    for (i = 0; i < row.length; i++) {
      c = alphabet.get_complement(i);
      if (c < i) continue;
      temp = row[i];
      row[i] = row[c];
      row[c] = temp;
    }
  }
  // complement pssm (if defined)
  if (typeof this.pssm !== "undefined") {
    //complement
    for (x = 0; x < this.motif_length; x++) {
      row = this.pssm[x];
      for (i = 0; i < row.length; i++) {
        c = alphabet.get_complement(i);
        if (c < i) continue;
        temp = row[i];
        row[i] = row[c];
        row[c] = temp;
      }
    }
  }
  return this; //allow function chaining...
};

Pspm.prototype.get_stack = function(position, alphabet, ssc) {
  "use strict";
  var row, stack_ic, alphabet_ic, stack, i, sym;
  if (this.alph_length != alphabet.get_size_core()) {
    throw new Error("The alphabet size does not match the size of the pspm.");
  }
  row = this.pspm[position];
  stack_ic = this.get_stack_ic(position, alphabet);
  if (ssc) stack_ic -= this.get_error(alphabet);
  alphabet_ic = alphabet.get_ic();
  stack = [];
  for (i = 0; i < this.alph_length; i++) {
    sym = new Symbol(i, row[i]*stack_ic/alphabet_ic, alphabet);
    if (sym.get_scale() <= 0) {
      continue;
    }
    stack.push(sym);
  }
  stack.sort(compare_symbol);
  return stack;
};

Pspm.prototype.get_stack_ic = function(position, alphabet) {
  "use strict";
  var row, H, i;
  if (this.alph_length != alphabet.get_size_core()) {
    throw new Error("The alphabet size does not match the size fo the pspm.");
  }
  row = this.pspm[position];
  H = 0;
  for (i = 0; i < this.alph_length; i++) {
    if (row[i] === 0) {
      continue;
    }
    H -= (row[i] * (Math.log(row[i]) / Math.LN2));
  }
  return alphabet.get_ic() - H;
};

Pspm.prototype.get_error = function(alphabet) {
  "use strict";
  if (this.nsites === 0) {
    return 0;
  }
  return (alphabet.get_size_core()-1) / (2 * Math.LN2 * this.nsites);
};

Pspm.prototype.get_motif_length = function() {
  "use strict";
  return this.motif_length;
};

Pspm.prototype.get_alph_length = function() {
  "use strict";
  return this.alph_length;
};

Pspm.prototype.get_left_trim = function() {
  "use strict";
  return this.ltrim;
};

Pspm.prototype.get_right_trim = function() {
  "use strict";
  return this.rtrim;
};

Pspm.prototype.as_best_match = function(alphabet) {
  "use strict";
  var match, odds, best_odds, best_index;
  var i, j;
  match = "";
  for (i = 0; i < this.motif_length; i++) {
    best_index = 0;
    best_odds = this.pspm[i][0] / alphabet.get_bg_freq(0);
    for (j = 1; j < this.alph_length; j++) {
      odds = this.pspm[i][j] / alphabet.get_bg_freq(j);
      if (odds > best_odds) {
        best_odds = odds;
        best_index = j;
      }
    }
    match += alphabet.get_symbol(best_index);
  }
  return match;
};

Pspm.prototype.as_count_matrix = function() {
  "use strict";
  var count, count_text, text;
  var i, j;
  text = "";
  for (i = 0; i < this.motif_length; i++) {
    if (i !== 0) {
      text += "\n";
    }
    for (j = 0; j < this.alph_length; j++) {
      if (j !== 0) {
        text += " ";
      }
      count = Math.round(this.nsites * this.pspm[i][j]);
      count_text = "" + count;
      // pad up to length of 4
      if (count_text.length < 4) {
        text += (new Array(5 - count_text.length)).join(" ") + count_text;
      } else {
        text += count_text;
      }
    }
  }
  return text; 
};

Pspm.prototype.as_probability_matrix = function() {
  "use strict";
  var text;
  var i, j;
  text = "";
  for (i = 0; i < this.motif_length; i++) {
    if (i !== 0) {
      text += "\n";
    }
    for (j = 0; j < this.alph_length; j++) {
      if (j !== 0) {
        text += " ";
      }
      text += this.pspm[i][j].toFixed(6);
    }
  }
  return text; 
};

Pspm.prototype.as_score_matrix = function(alphabet, pseudo) {
  "use strict";
  var me, score, out, row, col, score_text;
  me = this;
  if (typeof this.pssm === "undefined") {
    if (!(typeof alphabet === "object" && alphabet != null && alphabet instanceof Alphabet)) {
      throw new Error("The alphabet is required to generate the pssm.");
    }
    if (typeof pseudo === "undefined") {
      pseudo = 0.01;
    } else if (typeof pseudo !== "number" || pseudo < 0) {
      throw new Error("Expected positive number for pseudocount");
    }
    score = function(row, col) {
      "use strict";
      var p, bg, p2;
      p = me.pspm[row][col];
      bg = alphabet.get_bg_freq(col);
      p2 = (p * me.nsites + bg * pseudo) / (me.nsites + pseudo);
      return (p2 > 0 ? Math.round((Math.log(p2 / bg) / Math.LN2) * 100) : -10000);
    };
  } else {
    score = function(row, col) {
      "use strict";
      return me.pssm[row][col];
    };
  }
  out = "";
  for (row = 0; row < this.motif_length; row++) {
    for (col = 0; col < this.alph_length; col++) {
      if (col !== 0) {
        out += " ";
      }
      score_text = "" + score(row, col);
      // pad out to 6 characters
      if (score_text.length < 6) {
        out += (new Array(7 - score_text.length)).join(" ") + score_text;
      } else {
        out += score_text;
      }
    }
    out += "\n";
  }
  return out;
}

Pspm.prototype.as_pspm = function() {
  "use strict";
  return "letter-probability matrix: alength= " + this.alph_length + 
      " w= " + this.motif_length + " nsites= " + this.nsites + 
      " E= " + (typeof this.evalue === "number" ? 
          this.evalue.toExponential() : this.evalue) + "\n" +
      this.as_probability_matrix();
};

Pspm.prototype.as_pssm = function(alphabet, pseudo) {
  "use strict";
  return "log-odds matrix: alength= " + this.alph_length + 
      " w= " + this.motif_length + 
      " E= " + (typeof this.evalue == "number" ?
          this.evalue.toExponential() : this.evalue) + "\n" +
      this.as_score_matrix(alphabet, pseudo);
};

Pspm.prototype.as_meme = function(options) {
  var with_header, with_pspm, with_pssm, version, alphabet, bg_source, pseudocount, strands;
  var out, alen, i;
  // get the options
  if (typeof options !== "object" || options === null) {
    options = {};
  }
  with_header = (typeof options["with_header"] === "boolean" ? options["with_header"] : false);
  with_pspm = (typeof options["with_pspm"] === "boolean" ? options["with_pspm"] : false);
  with_pssm = (typeof options["with_pssm"] === "boolean" ? options["with_pssm"] : false);
  if (!with_pspm && !with_pssm) with_pspm = true;
  if (with_header) {
    if (typeof options["version"] === "string" && /^\d+(?:\.\d+){0,2}$/.test(options["version"])) {
      version = options["version"];
    } else if (typeof options["version"] === "number") {
      version = options["version"].toFixed(0);
    } else {
      version = "4";
    }
    if (typeof options["strands"] === "number" && options["strands"] === 1) {
      strands = 1;
    } else {
      strands = 2;
    }
    if (typeof options["bg_source"] === "string") {
      bg_source = options["bg_source"];
    } else {
      bg_source = "unknown source";
    }
    if (typeof options["alphabet"] === "object" && options["alphabet"] != null
        && options["alphabet"] instanceof Alphabet) {
      alphabet = options["alphabet"];
    } else {
      throw new Error("The alphabet is required to generate the header.");
    }
  }
  // now create the output
  out = "";
  if (with_header) {
    out = "MEME version " + version + "\n\n";
    out += alphabet.as_meme() + "\n";
    if (alphabet.has_complement()) { // assume DNA has both strands unless otherwise specified
      out += "strands: " + (strands === 1 ? "+" : "+ -") + "\n\n";
    }
    out += "Background letter frequencies (from " + bg_source + "):\n";
    alen = alphabet.get_size_core();
    for (i = 0; i < alen; i++) {
      if (i !== 0) {
        if (i % 9 === 0) { // maximum of nine entries per line
          out += "\n";
        } else {
          out += " ";
        }
      }
      out += alphabet.get_symbol(i) + " " + alphabet.get_bg_freq(i).toFixed(3);
    }
  }
  out += "\n\n";
  out += "MOTIF " + this.name + (this.alt == null ? "" : " " + this.alt);
  if (with_pssm) {
    out += "\n\n";
    out += this.as_pssm(options["alphabet"], options["pseudocount"]);
  }
  if (with_pspm) {
    out += "\n\n";
    out += this.as_pspm();
  }
  return out;
}

Pspm.prototype.toString = function() {
  "use strict";
  var str, i, row;
  str = "";
  for (i = 0; i < this.pspm.length; i++) {
    row = this.pspm[i];
    str += row.join("\t") + "\n";
  }
  return str;
};

function parse_pspm_properties(str) {
  "use strict";
  var parts, i, eqpos, before, after, properties, prop, num, num_re;
  num_re = /^((?:[+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|inf)$/;
  parts = trim(str).split(/\s+/);
  // split up words containing =
  for (i = 0; i < parts.length;) {
    eqpos = parts[i].indexOf("=");
    if (eqpos != -1) {
      before = parts[i].substr(0, eqpos);
      after = parts[i].substr(eqpos+1);
      if (before.length > 0 && after.length > 0) {
        parts.splice(i, 1, before, "=", after);
        i += 3;
      } else if (before.length > 0) {
        parts.splice(i, 1, before, "=");
        i += 2;
      } else if (after.length > 0) {
        parts.splice(i, 1, "=", after);
        i += 2;
      } else {
        parts.splice(i, 1, "=");
        i++;
      }
    } else {
      i++;
    }
  }
  properties = {};
  for (i = 0; i < parts.length; i += 3) {
    if (parts.length - i < 3) {
      throw new Error("Expected PSPM property was incomplete. "+
          "Remaing parts are: " + parts.slice(i).join(" "));
    }
    if (parts[i+1] !== "=") {
      throw new Error("Expected '=' in PSPM property between key and " +
          "value but got " + parts[i+1]); 
    }
    prop = parts[i].toLowerCase();
    num = parts[i+2];
    if (!num_re.test(num)) {
      throw new Error("Expected numeric value for PSPM property '" + 
          prop + "' but got '" + num + "'");
    }
    properties[prop] = num;
  }
  return properties;
}

function parse_pspm_string(pspm_string) {
  "use strict";
  var header_re, lines, first_line, line_num, col_num, alph_length, 
      motif_length, nsites, evalue, pspm, i, line, match, props, parts,
      j, prob;
  header_re = /^letter-probability\s+matrix:(.*)$/i;
  lines = pspm_string.split(/\n/);
  first_line = true;
  line_num = 0;
  col_num = 0;
  alph_length;
  motif_length;
  nsites;
  evalue;
  pspm = [];
  for (i = 0; i < lines.length; i++) {
    line = trim(lines[i]);
    if (line.length === 0) { 
      continue;
    }
    // check the first line for a header though allow matrices without it
    if (first_line) {
      first_line = false;
      match = header_re.exec(line);
      if (match !== null) {
        props = parse_pspm_properties(match[1]);
        if (props.hasOwnProperty("alength")) {
          alph_length = parseFloat(props["alength"]);
          if (alph_length != 4 && alph_length != 20) {
            throw new Error("PSPM property alength should be 4 or 20" +
                " but got " + alph_length);
          }
        }
        if (props.hasOwnProperty("w")) {
          motif_length = parseFloat(props["w"]);
          if (motif_length % 1 !== 0 || motif_length < 1) {
            throw new Error("PSPM property w should be an integer larger " +
                "than zero but got " + motif_length);
          }
        }
        if (props.hasOwnProperty("nsites")) {
          nsites = parseFloat(props["nsites"]);
          if (nsites <= 0) {
            throw new Error("PSPM property nsites should be larger than " +
                "zero but got " + nsites);
          }
        }
        if (props.hasOwnProperty("e")) {
          evalue = props["e"];
          if (evalue < 0) {
            throw new Error("PSPM property evalue should be " +
                "non-negative but got " + evalue);
          }
        }
        continue;
      }
    }
    pspm[line_num] = [];
    col_num = 0;
    parts = line.split(/\s+/);
    for (j = 0; j < parts.length; j++) {
      prob = parseFloat(parts[j]);
      if (prob != parts[j] || prob < 0 || prob > 1) {
        throw new Error("Expected probability but got '" + parts[j] + "'"); 
      }
      pspm[line_num][col_num] = prob;
      col_num++;
    }
    line_num++;
  }
  if (typeof motif_length === "number") {
    if (pspm.length != motif_length) {
      throw new Error("Expected PSPM to have a motif length of " + 
          motif_length + " but it was actually " + pspm.length);
    }
  } else {
    motif_length = pspm.length;
  }
  if (typeof alph_length !== "number") {
    alph_length = pspm[0].length;
    if (alph_length != 4 && alph_length != 20) {
      throw new Error("Expected length of first row in the PSPM to be " +
          "either 4 or 20 but got " + alph_length);
    }
  }
  for (i = 0; i < pspm.length; i++) {
    if (pspm[i].length != alph_length) {
      throw new Error("Expected PSPM row " + i + " to have a length of " + 
          alph_length + " but the length was " + pspm[i].length);
    }
  }
  return {"pspm": pspm, "motif_length": motif_length, 
    "alph_length": alph_length, "nsites": nsites, "evalue": evalue};
}
//======================================================================
// end Pspm object
//======================================================================

//======================================================================
// start Logo object
//======================================================================

var Logo = function(alphabet, options) {
  "use strict";
  this.alphabet = alphabet;
  this.fine_text = "";
  this.x_axis = 1;
  this.y_axis = true;
  this.xlate_nsyms = 1;
  this.xlate_start = null;
  this.xlate_end = null;
  this.pspm_list = [];
  this.pspm_column = [];
  this.rows = 0;
  this.columns = 0;
  if (typeof options === "string") {
    // the old method signature had fine_text here so we support that
    this.fine_text = options;
  } else if (typeof options === "object" && options != null) {
    this.fine_text = (typeof options.fine_text === "string" ? options.fine_text : "");
    this.x_axis = (typeof options.x_axis === "boolean" ? (options.x_axis ? 1 : 0) : 1);
    if (options.x_axis_hidden != null && options.x_axis_hidden) this.x_axis = -1;
    this.y_axis = (typeof options.y_axis === "boolean" ? options.y_axis : true);
    this.xlate_nsyms = (typeof options.xlate_nsyms === "number" ? options.xlate_nsyms : this.xlate_nsyms);
    this.xlate_start = (typeof options.xlate_start === "number" ? options.xlate_start : this.xlate_start);
    this.xlate_end = (typeof options.xlate_end === "number" ? options.xlate_end : this.xlate_end);
  }
};

Logo.prototype.add_pspm = function(pspm, column) {
  "use strict";
  var col;
  if (typeof column === "undefined") {
    column = 0;
  } else if (column < 0) {
    throw new Error("Column index out of bounds.");
  }
  this.pspm_list[this.rows] = pspm;
  this.pspm_column[this.rows] = column;
  this.rows++;
  col = column + pspm.get_motif_length();
  if (col > this.columns) {
    this.columns = col;
  }
};

Logo.prototype.get_columns = function() {
  "use strict";
  return this.columns;
};

Logo.prototype.get_xlate_nsyms = function() {
  "use strict";
  return this.xlate_nsyms;
};

Logo.prototype.get_xlate_start = function() {
  "use strict";
  return (this.xlate_start != null ? this.xlate_start : 0);
};

Logo.prototype.get_xlate_end = function() {
  "use strict";
  return (this.xlate_end != null ? this.xlate_end : this.columns * this.xlate_nsyms);
};

Logo.prototype.get_xlate_columns = function() {
  "use strict";
  return this.get_xlate_end() - this.get_xlate_start();
};

Logo.prototype.get_rows = function() {
  "use strict";
  return this.rows;
};

Logo.prototype.get_pspm = function(row_index) {
  "use strict";
  if (row_index < 0 || row_index >= this.rows) {
    throw new Error("INDEX_OUT_OF_BOUNDS");
  }
  return this.pspm_list[row_index];
};

Logo.prototype.get_offset = function(row_index) {
  "use strict";
  if (row_index < 0 || row_index >= this.rows) {
    throw new Error("INDEX_OUT_OF_BOUNDS");
  }
  return this.pspm_column[row_index];
};

Logo.prototype._as_eps_data = function(ssc, errbars) {
  var i, j, pos, stack_pos, pspm, stack, sym, out;
  out = "";
  for (i = 0; i < this.rows; i++) {
    out += "\nStartLine\n";
    // Indent
    for (j = 0; j < this.pspm_column[i]; j++) {
      out += "() startstack\nendstack\n\n";
    }
    pspm = this.pspm_list[i];
    if (pspm.get_left_trim() > 0) {
      out += "MuteColour\nDrawTrimEdge\n" + pspm.get_left_trim() + " DrawTrimBg\n";
    }
    for (pos = 0; pos < pspm.get_motif_length(); pos++) {
      if (pos != 0 && pos == pspm.get_left_trim()) { // enable full colour
        out += "DrawTrimEdge\nRestoreColour\n";
      } else if (pos == (pspm.get_motif_length() - pspm.get_right_trim())) {
        out += "MuteColour\n" + pspm.get_right_trim() + " DrawTrimBg\n";
      }
      out += "(" + (pos + 1) + ") startstack\n";
      stack = pspm.get_stack(pos, this.alphabet, ssc);
      for (stack_pos = 0; stack_pos < stack.length; stack_pos++) {
        sym = stack[stack_pos];
        out += " " + (sym.get_scale() * this.alphabet.get_ic()) + " (" + sym.get_symbol() + ") numchar\n";
      }
      if (errbars) {
        out += " " + pspm.get_error(this.alphabet) + " Ibeam\n";
      }
      out += "endstack\n\n";
    }
    if (pspm.get_right_trim() > 0 || pspm.get_left_trim() == pspm.get_motif_length()) {
      out += "RestoreColour\n";
    }
    out += "EndLine\n";
  }
  return out;
};

Logo.prototype.as_eps = function(options) {
  "use strict";
  if (this.xlate_nsyms != 1) throw new Error("Unsupported setting xlate_nsyms for EPS");
  if (this.xlate_start != null) throw new Error("Unsupported setting xlate_start for EPS");
  if (this.xlate_end != null) throw new Error("Unsupported setting xlate_end for EPS");

  var LOGOHEIGHT = 7.5; // default height of line in cm
  var cm2pts, height, width, now, ssc, errbars;
  if (typeof options === "undefined") {
    options = {};
  }
  cm2pts = 72 / 2.54;
  if (typeof options.logo_height == "number") {
    height = options.logo_height;
  } else {
    height = LOGOHEIGHT * this.rows;
  }
  if (typeof options.logo_width == "number") {
    width = options.logo_width;
  } else {
    width = this.columns + 2;
  }
  now = new Date();
  ssc = (typeof options.ssc == "boolean" ? options.ssc : false);
  errbars = (typeof options.show_error_bar == "boolean" ? options.show_error_bar : ssc);
  var values = {
    "LOGOHEIGHT": height,
    "LOGOWIDTH": width,
    "BOUNDINGHEIGHT": Math.round(height * cm2pts),
    "BOUNDINGWIDTH": Math.round(width * cm2pts),
    "LOGOLINEHEIGHT": (height / this.rows),
    "CHARSPERLINE": this.columns,
    "BARBITS": this.alphabet.get_ic(),
    "LOGOTYPE": (this.alphabet.has_complement() ? "NA" : "AA"),
    "CREATIONDATE": now.getDate() + "." + (now.getMonth() + 1) + "." + now.getFullYear() + " " + now.getHours() + ":" + now.getMinutes() + ":" + now.getSeconds(),
    "ERRORBARFRACTION": (typeof options.error_bar_fraction == "number" ? options.error_bar_fraction : 1.0),
    "TICBITS": (typeof options.ticbits == "number" ? options.ticbits : 1.0),
    "TITLE": (typeof options.title == "string" ? options.title : ""),
    "FINEPRINT": (typeof options.fineprint == "string" ? options.fineprint : this.fine_text),
    "XAXISLABEL": (typeof options.xaxislabel == "string" ? options.xaxislabel : ""),
    "YAXISLABEL": (typeof options.yaxislabel == "string" ? options.yaxislabel : "bits"),
    "SSC": ssc,
    "YAXIS": (typeof options.show_y_axis == "boolean" ? options.show_y_axis : this.y_axis),
    "SHOWENDS": (typeof options.show_ends == "boolean" ? options.show_ends : false),
    "ERRBAR": errbars,
    "OUTLINE": (typeof options.show_outline == "boolean" ? options.show_outline : false),
    "NUMBERING": (typeof options.show_numbering == "boolean" ? options.show_numbering : this.x_axis != 0),
    "SHOWINGBOX": (typeof options.show_box == "boolean" ? options.show_box : false),
    "CREATOR": (typeof options.creator == "string" ? options.creator : "motif_logo.js"),
    "FONTSIZE": (typeof options.label_font_size == "number" ? options.label_font_size : 12),
    "TITLEFONTSIZE": (typeof options.title_font_size == "number" ? options.title_font_size : 12),
    "SMALLFONTSIZE": (typeof options.small_font_size == "number" ? options.small_font_size : 6),
    "TOPMARGIN" : (typeof options.top_margin == "number" ? options.top_margin : 0.9),
    "BOTTOMMARGIN": (typeof options.bottom_margin == "number" ? options.bottom_margin : 0.9),
    "COLORDICT": this.alphabet._as_eps_dict(),
    "DATA": this._as_eps_data(ssc, errbars)
  };
  // now this requires that the script containing the template has been imported!
  return motif_logo_template(values);
};

//======================================================================
// end Logo object
//======================================================================

// calculate the exact size (in pixels) of an object drawn on the
// canvas assuming that the background of the canvas is transparent.
function canvas_bounds(ctx, cwidth, cheight) {
  "use strict";
  var data, r, c, top_line, bottom_line, left_line, right_line, 
      txt_width, txt_height;

  // extract the image data
  data = ctx.getImageData(0, 0, cwidth, cheight).data;

  // set initial values
  top_line = -1; bottom_line = -1; left_line = -1; right_line = -1;
  txt_width = 0; txt_height = 0;

  // Find the top-most line with a non-transparent pixel
  for (r = 0; r < cheight; r++) {
    for (c = 0; c < cwidth; c++) {
      if (data[r * cwidth * 4 + c * 4 + 3]) {
        top_line = r;
        break;
      }
    }
    if (top_line != -1) {
      break;
    }
  }
  
  // Only bother looking if we found at least one set pixel... 
  if (top_line != -1) {

    //find the last line with a non-transparent pixel
    for (r = cheight-1; r >= top_line; r--) {
      for(c = 0; c < cwidth; c++) {
        if(data[r * cwidth * 4 + c * 4 + 3]) {
          bottom_line = r;
          break;
        }
      }
      if (bottom_line != -1) {
        break;
      }
    }
    // calculate height
    txt_height = bottom_line - top_line + 1;

    // Find the left-most line with a non-transparent pixel
    for (c = 0; c < cwidth; c++) {
      for (r = top_line; r <= bottom_line; r++) {
        if (data[r * cwidth * 4 + c * 4 + 3]) {
          left_line = c;
          break;
        }
      }
      if (left_line != -1) {
        break;
      }
    }

    //find the right most line with a non-transparent pixel
    for (c = cwidth-1; c >= left_line; c--) {
      for(r = top_line; r <= bottom_line; r++) {
        if(data[r * cwidth * 4 + c * 4 + 3]) {
          right_line = c;
          break;
        }
      }
      if (right_line != -1) {
        break;
      }
    }
    txt_width = right_line - left_line + 1;
  }

  //return the bounds
  return {bound_top: top_line, bound_bottom: bottom_line, 
    bound_left: left_line, bound_right: right_line, width: txt_width, 
    height: txt_height};
}

//======================================================================
// start RasterizedAlphabet
//======================================================================

// Rasterize Alphabet
// 1) Measure width of text at default font for all symbols in alphabet
// 2) sort in width ascending
// 3) Drop the top and bottom 10% (designed to ignore outliers like 'W' and 'I')
// 4) Calculate the average as the maximum scaling factor (designed to stop I becoming a rectangular blob).
// 5) Assume scale of zero would result in width of zero, interpolate scale required to make perfect width font
// 6) Draw text onto temp canvas at calculated scale
// 7) Find bounds of drawn text
// 8) Paint on to another canvas at the desired height (but only scaling width to fit if larger).
var RasterizedAlphabet = function(alphabet, logo_scale, font, width) {
  "use strict";
  var default_size, safety_pad, canvas, ctx, middle, baseline, widths, sizes,
      i, sym, size, tenpercent, avg_width, scale, 
      target_width, target_height;
  //variable prototypes
  this.alphabet = alphabet;
  this.scale = logo_scale;
  this.sym_cache = {};
  this.stack_num_cache = [];
  this.scale_num_cache = [];
  // size of canvas
  default_size = 60; // size of measuring canvas
  safety_pad = 20; // pixels to pad around so we don't miss the edges
  // create a canvas to do our measuring
  canvas = document.createElement("canvas");
  if (!canvas.getContext) throw new Error("No canvas support");
  canvas.width = default_size + 2 * safety_pad;
  canvas.height = default_size + 2 * safety_pad;
  middle = Math.round(canvas.width / 2);
  baseline = Math.round(canvas.height - safety_pad);
  ctx = canvas.getContext('2d');
  if (!supports_text(ctx)) throw new Error("Canvas does not support text");
  ctx.font = font;
  ctx.textAlign = "center";
  ctx.translate(middle, baseline);
  // list of widths
  widths = [];
  sizes = [];
  //now measure each letter in the alphabet
  for (i = 0; i < alphabet.get_size_core(); ++i) {
    // reset the canvas
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = alphabet.get_colour(i);
    // draw the test text
    ctx.fillText(alphabet.get_symbol(i), 0, 0);
    //measure
    size = canvas_bounds(ctx, canvas.width, canvas.height);
    if (size.width === 0) throw new Error("Invisible symbol!");
    widths.push(size.width);
    sizes[i] = size;
  }
  //sort the widths
  widths.sort(function(a,b) {return a - b;});
  //drop 10% of the items off each end
  tenpercent = Math.floor(widths.length / 10);
  for (i = 0; i < tenpercent; ++i) {
    widths.pop();
    widths.shift();
  }
  //calculate average width
  avg_width = 0;
  for (i = 0; i < widths.length; ++i) {
    avg_width += widths[i];
  }
  avg_width /= widths.length;
  // calculate the target width
  target_width = width * this.scale * 2;
  // calculate scales
  for (i = 0; i < alphabet.get_size_core(); ++i) {
    sym = alphabet.get_symbol(i);
    size = sizes[i];
    // calculate scale
    scale = target_width / Math.max(avg_width, size.width);
    // estimate scaled height
    target_height = size.height * scale;
    // create an appropriately sized canvas
    canvas = document.createElement("canvas");
    canvas.width = target_width;
    canvas.height = target_height + safety_pad * 2;
    // calculate the middle
    middle = Math.round(canvas.width / 2);
    // calculate the baseline
    baseline = Math.round(canvas.height - safety_pad);
    // get the context and prepare to draw the rasterized text
    ctx = canvas.getContext('2d');
    ctx.font = font;
    ctx.fillStyle = alphabet.get_colour(i);
    ctx.textAlign = "center";
    ctx.translate(middle, baseline);
    ctx.save();
    ctx.scale(scale, scale);
    // draw the text
    ctx.fillText(sym, 0, 0);
    ctx.restore();
    this.sym_cache[sym] = {"image": canvas, "size": canvas_bounds(ctx, canvas.width, canvas.height)};
  }
};

RasterizedAlphabet.prototype.get_alphabet = function() {
  return this.alphabet;
};

RasterizedAlphabet.prototype.get_scale = function() {
  return this.scale;
};

RasterizedAlphabet.prototype.draw_stack_sym = function(ctx, letter, dx, dy, dWidth, dHeight) {
  "use strict";
  var entry, image, size;
  entry = this.sym_cache[letter];
  image = entry.image;
  size = entry.size;
  ctx.drawImage(image, 0, size.bound_top -1, image.width, size.height+1, dx, dy, dWidth, dHeight);
};

RasterizedAlphabet.prototype.draw_stack_num = function(ctx, font, stack_width, index) {
  var image, image_ctx, text_length;
  if (index >= this.stack_num_cache.length) {
    image = document.createElement("canvas");
    // measure the text
    image_ctx = image.getContext('2d');
    image_ctx.save();
    image_ctx.font = font;
    text_length = image_ctx.measureText("" + (index + 1)).width;
    image_ctx.restore();
    // resize the canvas to fit
    image.width = Math.ceil(stack_width);
    image.height = Math.ceil(text_length);
    // draw the text
    image_ctx = image.getContext('2d');
    image_ctx.translate(Math.round(stack_width / 2), 0);
    image_ctx.font = font;
    image_ctx.textBaseline = "middle";
    image_ctx.textAlign = "right";
    image_ctx.rotate(-(Math.PI / 2));
    image_ctx.fillText("" + (index + 1), 0, 0);
    this.stack_num_cache[index] = image;
  } else {
    image = this.stack_num_cache[index];
  }
  ctx.drawImage(image, 0, 0);
}

RasterizedAlphabet.prototype.draw_scale_num = function(ctx, font, num) {
  var image, image_ctx, text_size, m_length;
  if (num >= this.scale_num_cache.length) {
    image = document.createElement("canvas");
    // measure the text
    image_ctx = image.getContext('2d');
    image_ctx.font = font;
    text_size = image_ctx.measureText("" + num);
    if (text_size.actualBoundingBoxAscent && text_size.actualBoundingBoxDesent) {
      // resize the canvas to fit
      image.width = Math.ceil(text_size.width);
      image.height = Math.ceil(text_size.actualBoundingBoxAscent + text_size.actualBoundingBoxDesent);
      // draw the text
      image_ctx = image.getContext('2d');
      image_ctx.font = font;
      image_ctx.textAlign = "right";
      image_ctx.fillText("" + num, image.width, text_size.actualBoundingBoxAscent);
    } else {
      // measure width of 'm' to approximate height, we double it later anyway
      m_length = image_ctx.measureText("m").width;
      // resize the canvas to fit
      image.width = Math.ceil(text_size.width);
      image.height = Math.ceil(2 * m_length);
      // draw the text
      image_ctx = image.getContext('2d');
      image_ctx.font = font;
      image_ctx.textAlign = "right";
      image_ctx.textBaseline = "middle";
      image_ctx.fillText("" + num, image.width, m_length);
    }
    this.scale_num_cache[num] = image;
  } else {
    image = this.scale_num_cache[num];
  }
  ctx.drawImage(image, -image.width, -Math.round(image.height / 2))
}

//======================================================================
// end RasterizedAlphabet
//======================================================================

//======================================================================
// start LogoMetrics object
//======================================================================

var LogoMetrics = function(ctx, logo_columns, logo_rows, has_names, has_finetext, x_axis, y_axis) {
  "use strict";
  var i, row_height;
  //variable prototypes
  this.pad_top = (has_names ? 5 : 0);
  this.pad_left = (y_axis ? 10 : 0);
  this.pad_right = (has_finetext ? 15 : 0);
  this.pad_bottom = 0;
  this.pad_middle = 20;
  this.name_height = 14;
  this.name_font = "bold " + this.name_height + "px Times, sans-serif";
  this.name_spacer = 0;
  this.y_axis = y_axis;
  this.y_label = "bits";
  this.y_label_height = 12;
  this.y_label_font = "bold " + this.y_label_height + "px Helvetica, sans-serif";
  this.y_label_spacer = 3;
  this.y_num_height = 12;
  this.y_num_width = 0;
  this.y_num_font = "bold " + this.y_num_height + "px Helvetica, sans-serif";
  this.y_tic_width = 5;
  this.stack_pad_left = 0;
  this.stack_font = "bold 25px Helvetica, sans-serif";
  this.stack_height = 90;
  this.stack_width = 26;
  this.stacks_pad_right = 5;
  this.x_axis = x_axis;
  this.x_num_above = 2;
  this.x_num_height = 12;
  this.x_num_width = 0;
  this.x_num_font = "bold " + this.x_num_height + "px Helvetica, sans-serif";
  this.fine_txt_height = 6;
  this.fine_txt_above = 2;
  this.fine_txt_font = "normal " + this.fine_txt_height + "px Helvetica, sans-serif";
  this.letter_metrics = new Array();
  this.summed_width = 0;
  this.summed_height = 0;
  //calculate the width of the y axis numbers
  ctx.font = this.y_num_font;
  for (i = 0; i <= 2; i++) {
    this.y_num_width = Math.max(this.y_num_width, ctx.measureText("" + i).width);
  }
  //calculate the width of the x axis numbers (but they are rotated so it becomes height)
  if (x_axis == 1) {
    ctx.font = this.x_num_font;
    for (i = 1; i <= logo_columns; i++) {
      this.x_num_width = Math.max(this.x_num_width, ctx.measureText("" + i).width);
    }
  } else if (x_axis == 0) {
    this.x_num_height = 4;
    this.x_num_width = 4;
  } else {
    this.x_num_height = 0;
    this.x_num_width = 0;
  }
  
  //calculate how much vertical space we want to draw this
  //first we add the padding at the top and bottom since that's always there
  this.summed_height += this.pad_top + this.pad_bottom;
  //all except the last row have the same amount of space allocated to them
  if (logo_rows > 1) {
    row_height = this.stack_height + this.pad_middle;
    if (has_names) {
      row_height += this.name_height;
      //the label is allowed to overlap into the spacer
      row_height += Math.max(this.y_num_height/2, this.name_spacer); 
      //the label is allowed to overlap the space used by the other label
      row_height += Math.max(this.y_num_height/2, this.x_num_height + this.x_num_above); 
    } else {
      row_height += this.y_num_height/2; 
      //the label is allowed to overlap the space used by the other label
      row_height += Math.max(this.y_num_height/2, this.x_num_height + this.x_num_above); 
    }
    this.summed_height += row_height * (logo_rows - 1);
  }
  //the last row has the name and fine text below it but no padding
  this.summed_height += this.stack_height + (this.y_axis ? this.y_num_height/2 : 0);

  var fine_txt_total = (has_finetext ? this.fine_txt_height + this.fine_txt_above : 0);
  if (has_names) {
    this.summed_height += fine_txt_total + this.name_height;
    this.summed_height += Math.max((this.y_axis ? this.y_num_height/2 : 0), 
        this.x_num_height + this.x_num_above + this.name_spacer);
  } else {
    this.summed_height += Math.max((this.y_axis ? this.y_num_height/2 : 0), 
        this.x_num_height + this.x_num_above + fine_txt_total);
  }

  //calculate how much horizontal space we want to draw this
  //first add the padding at the left and right since that's always there
  this.summed_width += this.pad_left + this.pad_right;
  if (this.y_axis) {
    //add on the space for the y-axis label
    this.summed_width += this.y_label_height + this.y_label_spacer;
    //add on the space for the y-axis
    this.summed_width += this.y_num_width + this.y_tic_width;
  }
  //add on the space for the stacks
  this.summed_width += (this.stack_pad_left + this.stack_width) * logo_columns;
  //add on the padding after the stacks (an offset from the fine text)
  this.summed_width += this.stacks_pad_right;

};

//======================================================================
// end LogoMetrics object
//======================================================================

//found this trick at http://talideon.com/weblog/2005/02/detecting-broken-images-js.cfm
function image_ok(img) {
  "use strict";
  // During the onload event, IE correctly identifies any images that
  // weren't downloaded as not complete. Others should too. Gecko-based
  // browsers act like NS4 in that they report this incorrectly.
  if (!img.complete) {
    return false;
  }
  // However, they do have two very useful properties: naturalWidth and
  // naturalHeight. These give the true size of the image. If it failed
  // to load, either of these should be zero.
  if (typeof img.naturalWidth !== "undefined" && img.naturalWidth === 0) {
    return false;
  }
  // No other way of checking: assume it's ok.
  return true;
}
  
function supports_text(ctx) {
  "use strict";
  if (!ctx.fillText) {
    return false;
  }
  if (!ctx.measureText) {
    return false;
  }
  return true;
}

//draws the scale, returns the width
function draw_scale(ctx, metrics, alphabet_ic, raster) {
  "use strict";
  var tic_height, i;
  tic_height = metrics.stack_height / alphabet_ic;
  ctx.save();
  ctx.translate(metrics.y_label_height, metrics.y_num_height/2);
  //draw the axis label
  ctx.save();
  ctx.font = metrics.y_label_font;
  ctx.translate(0, metrics.stack_height/2);
  ctx.rotate(-(Math.PI / 2));
  ctx.textAlign = "center";
  ctx.fillText("bits", 0, 0);
  ctx.restore();

  ctx.translate(metrics.y_label_spacer + metrics.y_num_width, 0);

  //draw the axis tics
  ctx.save();
  ctx.translate(0, metrics.stack_height);
  for (i = 0; i <= alphabet_ic; i++) {
    //draw the number
    ctx.save();
    ctx.translate(-1, 0);
    raster.draw_scale_num(ctx, metrics.y_num_font, i);
    ctx.restore();
    //draw the tic
    ctx.fillRect(0, -1, metrics.y_tic_width, 2);
    //prepare for next tic
    ctx.translate(0, -tic_height);
  }
  ctx.restore();

  ctx.fillRect(metrics.y_tic_width - 2, 0, 2, metrics.stack_height)

  ctx.restore();
}

function draw_stack_num(ctx, metrics, row_index, raster) {
  "use strict";
  ctx.save();
  ctx.translate(0, Math.round(metrics.stack_height + metrics.x_num_above));
  if (metrics.x_axis == 1) {
    raster.draw_stack_num(ctx, metrics.x_num_font, metrics.stack_width, row_index);
  } else if (metrics.x_axis == 0) {
    // draw dots instead of the numbers (good for small logos)
    ctx.beginPath();
    var radius = Math.round(metrics.x_num_height / 2);
    ctx.arc(Math.round(metrics.stack_width / 2), radius, radius, 0, 2 * Math.PI, false);
    ctx.fill();
  }
  ctx.restore();
}

function draw_stack(ctx, metrics, symbols, raster) {
  "use strict";
  var preferred_pad, sym_min, i, sym, sym_height, pad;
  preferred_pad = 0;
  sym_min = 5;

  ctx.save();//1
  ctx.translate(0, metrics.stack_height);
  for (i = 0; i < symbols.length; i++) {
    sym = symbols[i];
    sym_height = metrics.stack_height * sym.get_scale();
    
    pad = preferred_pad;
    if (sym_height - pad < sym_min) {
      pad = Math.min(pad, Math.max(0, sym_height - sym_min));
    }
    sym_height -= pad;

    //translate to the correct position
    ctx.translate(0, -(pad/2 + sym_height));

    //draw
    raster.draw_stack_sym(ctx, sym.get_symbol(), 0, 0, metrics.stack_width, sym_height);
    //translate past the padding
    ctx.translate(0, -(pad/2));
  }
  ctx.restore();//1
}

function draw_dashed_line(ctx, pattern, start, x1, y1, x2, y2) {
  "use strict";
  var x, y, len, i, dx, dy, tlen, theta, mulx, muly, lx, ly;
  dx = x2 - x1;
  dy = y2 - y1;
  tlen = Math.pow(dx*dx + dy*dy, 0.5);
  theta = Math.atan2(dy,dx);
  mulx = Math.cos(theta);
  muly = Math.sin(theta);
  lx = [];
  ly = [];
  for (i = 0; i < pattern; ++i) {
    lx.push(pattern[i] * mulx);
    ly.push(pattern[i] * muly);
  }
  i = start;
  x = x1;
  y = y1;
  len = 0;
  ctx.beginPath();
  while (len + pattern[i] < tlen) {
    ctx.moveTo(x, y);
    x += lx[i];
    y += ly[i];
    ctx.lineTo(x, y);
    len += pattern[i];
    i = (i + 1) % pattern.length;
    x += lx[i];
    y += ly[i];
    len += pattern[i];
    i = (i + 1) % pattern.length;
  }
  if (len < tlen) {
    ctx.moveTo(x, y);
    x += mulx * (tlen - len);
    y += muly * (tlen - len);
    ctx.lineTo(x, y);
  }
  ctx.stroke();
}

function draw_trim_background(ctx, metrics, left_start, left_end, left_divider, right_start, right_end, right_divider) {
  "use strict";
  var left_size = left_end - left_start;
  var right_size = right_end - right_start;
  var line_x;

  ctx.save();//s8
  ctx.fillStyle = "rgb(240, 240, 240)";
  if (left_size > 0) {
    ctx.fillRect(left_start * metrics.stack_width, 0, left_size * metrics.stack_width, metrics.stack_height);
  }
  if (right_size > 0) {
    ctx.fillRect(right_start * metrics.stack_width, 0, right_size * metrics.stack_width, metrics.stack_height);
  }
  ctx.fillStyle = "rgb(51, 51, 51)";
  if (left_size > 0 && left_divider) {
    line_x = (left_end * metrics.stack_width) - 0.5;
    draw_dashed_line(ctx, [3], 0, line_x, 0, line_x, metrics.stack_height);
  }
  if (right_size > 0 && right_divider) {
    line_x = (right_start * metrics.stack_width) + 0.5;
    draw_dashed_line(ctx, [3], 0, line_x, 0, line_x, metrics.stack_height);
  }
  ctx.restore();//s8
}

function size_logo_on_canvas(logo, canvas, show_names, scale) {
  "use strict";
  var draw_name, draw_finetext, metrics;
  draw_name = (typeof show_names === "boolean" ? show_names : (logo.get_rows() > 1));
  draw_finetext = (logo.fine_text.length > 0);
  if (canvas.width !== 0 && canvas.height !== 0) {
    return;
  }
  metrics = new LogoMetrics(canvas.getContext('2d'), 
      logo.get_xlate_columns(), logo.get_rows(), draw_name, draw_finetext, logo.x_axis, logo.y_axis);
  if (typeof scale == "number") {
    //resize the canvas to fit the scaled logo
    canvas.width = metrics.summed_width * scale;
    canvas.height = metrics.summed_height * scale;
  } else {
    if (canvas.width === 0 && canvas.height === 0) {
      canvas.width = metrics.summed_width;
      canvas.height = metrics.summed_height;
    } else if (canvas.width === 0) {
      canvas.width = metrics.summed_width * (canvas.height / metrics.summed_height);
    } else if (canvas.height === 0) {
      canvas.height = metrics.summed_height * (canvas.width / metrics.summed_width);
    }
  }
}

function draw_logo_on_canvas(logo, canvas, show_names, scale) {
  "use strict";
  var i, draw_name, draw_finetext, ctx, metrics, raster, pspm_i, pspm, 
      offset, col_index, motif_position, ssc;
  ssc = false;
  draw_name = (typeof show_names === "boolean" ? show_names : (logo.get_rows() > 1));
  draw_finetext = (logo.fine_text.length > 0);
  ctx = canvas.getContext('2d');
  //assume that the user wants the canvas scaled equally so calculate what the best width for this image should be
  metrics = new LogoMetrics(ctx, logo.get_xlate_columns(), logo.get_rows(), draw_name, draw_finetext, logo.x_axis, logo.y_axis);
  if (typeof scale == "number") {
    //resize the canvas to fit the scaled logo
    canvas.width = metrics.summed_width * scale;
    canvas.height = metrics.summed_height * scale;
  } else {
    if (canvas.width === 0 && canvas.height === 0) {
      scale = 1;
      canvas.width = metrics.summed_width;
      canvas.height = metrics.summed_height;
    } else if (canvas.width === 0) {
      scale = canvas.height / metrics.summed_height;
      canvas.width = metrics.summed_width * scale;
    } else if (canvas.height === 0) {
      scale = canvas.width / metrics.summed_width;
      canvas.height = metrics.summed_height * scale;
    } else {
      scale = Math.min(canvas.width / metrics.summed_width, canvas.height / metrics.summed_height);
    }
  }
  // cache the raster based on the assumption that we will be drawing a lot
  // of logos the same size and alphabet
  if (typeof draw_logo_on_canvas.raster_cache === "undefined") {
    draw_logo_on_canvas.raster_cache = [];
  }
  for (i = 0; i < draw_logo_on_canvas.raster_cache.length; i++) {
    raster = draw_logo_on_canvas.raster_cache[i];
    if (raster.get_alphabet().equals(logo.alphabet) &&
        Math.abs(raster.get_scale() - scale) < 0.1) break;
    raster = null;
  }
  if (raster == null) {
    raster = new RasterizedAlphabet(logo.alphabet, scale, metrics.stack_font, metrics.stack_width);
    draw_logo_on_canvas.raster_cache.push(raster);
  }
  ctx = canvas.getContext('2d');
  ctx.save();//s1
  ctx.scale(scale, scale);
  ctx.save();//s2
  ctx.save();//s7
  //create margin
  ctx.translate(Math.round(metrics.pad_left), Math.round(metrics.pad_top));
  for (pspm_i = 0; pspm_i < logo.get_rows(); ++pspm_i) {
    pspm = logo.get_pspm(pspm_i);
    offset = logo.get_offset(pspm_i);
    //optionally draw name if this isn't the last row or is the only row 
    if (draw_name && (logo.get_rows() == 1 || pspm_i != (logo.get_rows()-1))) {
      ctx.save();//s4
      ctx.translate(Math.round(metrics.summed_width/2), Math.round(metrics.name_height));
      ctx.font = metrics.name_font;
      ctx.textAlign = "center";
      ctx.fillText(pspm.name, 0, 0);
      ctx.restore();//s4
      ctx.translate(0, Math.round(metrics.name_height + 
          Math.min(0, metrics.name_spacer - metrics.y_num_height/2)));
    }
    //draw scale
    if (logo.y_axis) draw_scale(ctx, metrics, logo.alphabet.get_ic(), raster);
    ctx.save();//s5
    //translate across past the scale
    if (logo.y_axis) {
      ctx.translate(Math.round(metrics.y_label_height + metrics.y_label_spacer + 
        metrics.y_num_width + metrics.y_tic_width), Math.round(metrics.y_num_height / 2));
    }
    //draw the trimming background
    if (pspm.get_left_trim() > 0 || pspm.get_right_trim() > 0) {
      var left_start = offset * logo.get_xlate_nsyms();
      var left_end = (offset + pspm.get_left_trim()) * logo.get_xlate_nsyms();
      var left_divider = true;
      if (left_end < logo.get_xlate_start() || left_start > logo.get_xlate_end()) {
        // no overlap
        left_start = 0;
        left_end = 0;
        left_divider = false;
      } else {
        if (left_start < logo.get_xlate_start()) {
          left_start = logo.get_xlate_start();
        }
        if (left_end > logo.get_xlate_end()) {
          left_end = logo.get_xlate_end();
          left_divider = false;
        }
        left_start -= logo.get_xlate_start();
        left_end -= logo.get_xlate_start();
        if (left_end < left_start) {
          left_start = 0;
          left_end = 0;
          left_divider = false;
        }
      }
      var right_end = (offset + pspm.get_motif_length()) * logo.get_xlate_nsyms();
      //var right_start = right_end - (pspm.get_left_trim() * logo.get_xlate_nsyms());
      var right_start = right_end - (pspm.get_right_trim() * logo.get_xlate_nsyms());
      var right_divider = true;
      if (right_end < logo.get_xlate_start() || right_start > logo.get_xlate_end()) {
        // no overlap
        right_start = 0;
        right_end = 0;
        right_divider = false;
      } else {
        if (right_start < logo.get_xlate_start()) {
          right_start = logo.get_xlate_start();
          right_divider = false;
        }
        if (right_end > logo.get_xlate_end()) {
          right_end = logo.get_xlate_end();
        }
        right_start -= logo.get_xlate_start();
        right_end -= logo.get_xlate_start();
        if (right_end < right_start) {
          right_start = 0;
          right_end = 0;
          right_divider = false;
        }
      }
      draw_trim_background(ctx, metrics, left_start, left_end, left_divider, right_start, right_end, right_divider);
    }
    //draw letters
    var xlate_col;
    for (xlate_col = logo.get_xlate_start(); xlate_col < logo.get_xlate_end(); xlate_col++) {
      ctx.translate(metrics.stack_pad_left,0);
      col_index = Math.floor(xlate_col / logo.get_xlate_nsyms());
      if (xlate_col % logo.get_xlate_nsyms() == 0) {
        if (col_index >= offset && col_index < (offset + pspm.get_motif_length())) {
          motif_position = col_index - offset;
          draw_stack_num(ctx, metrics, motif_position, raster);
          draw_stack(ctx, metrics, pspm.get_stack(motif_position, logo.alphabet, ssc), raster);
        }
      } else {
        if (col_index >= offset && col_index < (offset + pspm.get_motif_length())) {
          ctx.save();// s5.1
          ctx.translate(0, Math.round(metrics.stack_height));
          // TODO draw a dot or dash or something to indicate continuity of the motif
          ctx.restore(); //s5.1
        }
      }
      ctx.translate(Math.round(metrics.stack_width), 0);
    }
    ctx.restore();//s5
    ////optionally draw name if this is the last row but isn't the only row 
    if (draw_name && (logo.get_rows() != 1 && pspm_i == (logo.get_rows()-1))) {
      //translate vertically past the stack and axis's        
      ctx.translate(0, metrics.y_num_height/2 + metrics.stack_height + 
          Math.max(metrics.y_num_height/2, metrics.x_num_above + metrics.x_num_width + metrics.name_spacer));

      ctx.save();//s6
      ctx.translate(metrics.summed_width/2, metrics.name_height);
      ctx.font = metrics.name_font;
      ctx.textAlign = "center";
      ctx.fillText(pspm.name, 0, 0);
      ctx.restore();//s6
      ctx.translate(0, metrics.name_height);
    } else {
      //translate vertically past the stack and axis's        
      ctx.translate(0, metrics.y_num_height/2 + metrics.stack_height + 
          Math.max(metrics.y_num_height/2, metrics.x_num_above + metrics.x_num_width));
    }
    //if not the last row then add middle padding
    if (pspm_i != (logo.get_rows() -1)) {
      ctx.translate(0, metrics.pad_middle);
    }
  }
  ctx.restore();//s7
  if (logo.fine_text.length > 0) {
    ctx.translate(metrics.summed_width - metrics.pad_right, metrics.summed_height - metrics.pad_bottom);
    ctx.font = metrics.fine_txt_font;
    ctx.textAlign = "right";
    ctx.fillText(logo.fine_text, 0,0);
  }
  ctx.restore();//s2
  ctx.restore();//s1
}

function create_canvas(c_width, c_height, c_id, c_title, c_display) {
  "use strict";
  var canvas = document.createElement("canvas");
  //check for canvas support before attempting anything
  if (!canvas.getContext) {
    return null;
  }
  var ctx = canvas.getContext('2d');
  //check for html5 text drawing support
  if (!supports_text(ctx)) {
    return null;
  }
  //size the canvas
  canvas.width = c_width;
  canvas.height = c_height;
  canvas.id = c_id;
  canvas.title = c_title;
  canvas.style.display = c_display;
  return canvas;
}

function logo_1(alphabet, fine_text, pspm) {
  "use strict";
  var logo = new Logo(alphabet, fine_text);
  logo.add_pspm(pspm);
  return logo;
}

function logo_2(alphabet, fine_text, target, query, query_offset) {
  "use strict";
  var logo = new Logo(alphabet, fine_text);
  if (query_offset < 0) {
    logo.add_pspm(target, -query_offset);
    logo.add_pspm(query);
  } else {
    logo.add_pspm(target);
    logo.add_pspm(query, query_offset);
  }      
  return logo;
}

/*
 * Specifies an alternate source for an image.
 * If the image with the image_id specified has
 * not loaded then a generated logo will be used 
 * to replace it.
 *
 * Note that the image must either have dimensions
 * or a scale must be set.
 */
function alternate_logo(logo, image_id, scale) {
  "use strict";
  var image = document.getElementById(image_id);
  if (!image) {
    alert("Can't find specified image id (" +  image_id + ")");
    return;
  }
  //if the image has loaded then there is no reason to use the canvas
  if (image_ok(image)) {
    return;
  }
  //the image has failed to load so replace it with a canvas if we can.
  var canvas = create_canvas(image.width, image.height, image_id, image.title, image.style.display);
  if (canvas === null) {
    return;
  }
  //draw the logo on the canvas
  draw_logo_on_canvas(logo, canvas, null, scale);
  //replace the image with the canvas
  image.parentNode.replaceChild(canvas, image);
}

/*
 * Specifes that the element with the specified id
 * should be replaced with a generated logo.
 */
function replace_logo(logo, replace_id, scale, title_txt, display_style) {
  "use strict";
  var element = document.getElementById(replace_id);
  if (!replace_id) {
    alert("Can't find specified id (" + replace_id + ")");
    return;
  }
  //found the element!
  var canvas = create_canvas(50, 120, replace_id, title_txt, display_style);
  if (canvas === null) {
    return;
  }
  //draw the logo on the canvas
  draw_logo_on_canvas(logo, canvas, null, scale);
  //replace the element with the canvas
  element.parentNode.replaceChild(canvas, element);
}

/*
 * Fast string trimming implementation found at
 * http://blog.stevenlevithan.com/archives/faster-trim-javascript
 *
 * Note that regex is good at removing leading space but
 * bad at removing trailing space as it has to first go through
 * the whole string.
 */
function trim (str) {
  "use strict";
  var ws, i;
  str = str.replace(/^\s\s*/, '');
  ws = /\s/; i = str.length;
  while (ws.test(str.charAt(--i)));
  return str.slice(0, i + 1);
}
