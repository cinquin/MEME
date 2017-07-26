
// Custom error
var AlphabetError = function(message, reasons) {
  "use strict";
  this.message = message;
  this.stack = Error().stack;
  this.reasons = reasons;
};
AlphabetError.prototype = Object.create(Error.prototype);
AlphabetError.prototype.name = "AlphabetError";
AlphabetError.prototype.constructor = AlphabetError;
/*
 * The purpose of this parser is to aid in creating the Alphabet objects from
 * various file sources. The actual Alphabet object only supports taking
 * in a normalized and correctly ordered Javascript object definition.
 */
var AlphabetParser = function(handler) {
  this.handler = (handler != null ? handler : {});
  this.name = null;
  this.like = null;
  this.ncore = 0;
  this.error = false;
  this.lineNo = 0;
  this.seenHeader = false;
  this.seenAmbiguous = false;
  this.coreSymbols = {};
  this.allSymbols = {};
  this.symbolMap = {};
  // regular expressions
  var SYMBOL_RX = "[A-Za-z0-9\\?\\.\\*\\-]";
  var COLOUR_RX = "[A-Fa-f0-9]{6}";
  var NAME_RX = function(group) {
    // This strange combinatination of lookahead and back-reference is designed
    // to stop back-tracking: it emulates an "atomic group". Without it this
    // regex can be very slow! Note that the lookahead doesn't seem to capture
    // so there must be a capture group within.
    return "(?=(\"(?:[^\\\\\"]+|\\\\(?:[\"\\\\/bfnrt]|u[0-9A-Fa-f]{4}))*\"))\\" + group;
  };
  var CORE_RX = function(group) {
    return "(" + SYMBOL_RX + ")(?:\\s+" + NAME_RX(group + 1) + ")?(\\s+" + COLOUR_RX + ")?";
  };
  this.COLOUR_RE = new RegExp(COLOUR_RX);
  this.SYMBOL_RE = new RegExp(SYMBOL_RX);
  this.SYMBOLS_RE = new RegExp(SYMBOL_RX + "+");
  this.HEADER_RE = new RegExp("^\\s*ALPHABET(?:\\s+v1)?(?:\\s+" + NAME_RX(1) + ")?(?:\\s+(RNA|DNA|PROTEIN)-LIKE)?\\s*$");
  this.CORE_SINGLE_RE = new RegExp("^\\s*" + CORE_RX(1) + "\\s*$");
  this.CORE_PAIR_RE = new RegExp("^\\s*" + CORE_RX(1) + "\\s*~\\s*" + CORE_RX(4) + "\\s*$");
  this.AMBIG_RE = new RegExp("^\\s*" + CORE_RX(1) + "\\s*=\\s*(" + SYMBOL_RX + "+)\\s*$");
};

//
// throw an error
//
AlphabetParser.prototype._parse_error = function(message) {
  "use strict";
  this.error = true;
  throw new AlphabetError(message);
};

//
// Parse a header
//
AlphabetParser.prototype.parse_header = function(name, like) {
  "use strict";
  if (this.error) throw new Error("An error occurred previously");
  if (like != null && !like.match("^(?:RNA|DNA|PROTEIN)$"))
      throw new Error("The like parameter must be RNA, DNA or PROTEIN");
  if (!this.seenHeader) {
    this.seenHeader = true;
    this.name = name;
    this.like = like;
  } else {
    this._parse_error("Header already defined!");
  }
};

//
// Parse a symbol.
//
AlphabetParser.prototype.parse_symbol = function(sym, features) {
  "use strict";
  if (this.error) throw new Error("An error occurred previously");
  var i;
  if (features == null) features = {};
  if (!this.SYMBOL_RE.test(sym)) this._parse_error("Symbol invalid");
  var name = null;
  if (features.name != null && features.name != "") {
    name = features.name;
  }
  var colour = 0;
  if (features.colour != null && features.colour != "") {
    if (!this.COLOUR_RE.exec(features.colour)) this._parse_error("Colour invalid");
    colour = parseInt(features.colour, 16);
  }
  var complement = null;
  if (features.complement != null && features.complement != "") {
    if (!this.SYMBOL_RE.test(features.complement)) this._parse_error("Complement invalid");
    complement = features.complement;
  }
  var comprise = null;
  if (features.comprise != null && features.comprise != "") {
    if (!this.SYMBOLS_RE.test(features.comprise)) this._parse_error("Comprise invalid");
    // sort and remove duplicates
    comprise = features.comprise.split("").sort(AlphabetUtil.sym_compare);
    comprise = comprise.filter(function (item, pos, ary) { return pos == 0 || item != ary[pos - 1]; });
  }
  var aliases = [];
  if (features.aliases != null && features.aliases != "") {
    if (!this.SYMBOLS_RE.test(features.aliases)) this._parse_error("Aliases invalid");
    aliases = features.aliases.split("").sort(AlphabetUtil.sym_compare);
  }
  if (this.allSymbols[sym]) {
    this._parse_error("Symbol already defined");
    return;
  }
  for (i = 0; i < aliases.length; i++) {
    if (this.allSymbols[aliases[i]]) {
      this._parse_error("Symbol already defined");
      return;
    }
  }
  var symbol = {
    "sym": sym,
    "aliases": aliases,
    "name": name,
    "colour": colour,
    "complement": complement,
    "comprise": comprise
  };
  var key;
  if (comprise != null) {
    // ambiguous symbol!
    this.seenAmbiguous = true;
    for (i = 0; i < comprise.length; i++) {
      if (!this.coreSymbols[comprise[i]]) {
        this._parse_error("Unknown core symbol");
        return;
      }
    }
    key = comprise.join("");
  } else {
    // core symbol!
    if (this.seenAmbiguous) {
      this._parse_error("Core symbol after ambig");
      return;
    }
    key = sym;
    this.ncore++;
  }
  if (sym == '?' && (comprise == null || comprise.length != this.ncore)) {
    this._parse_error("Symbol '?' is reserved for use as a wildcard only.");
    return;
  }
  if (this.symbolMap[key] == null) {
    this.symbolMap[key] = [];
  }
  this.symbolMap[key].push(symbol);
  this.allSymbols[sym] = true;
  for (i = 0; i < aliases.length; i++) this.allSymbols[aliases[i]] = true;
  if (comprise == null) this.coreSymbols[sym] = true;
};

//
// Parse a line from an alphabet file.
//
AlphabetParser.prototype.parse_line = function(line) {
  "use strict";
  if (this.error) throw new Error("An error occurred previously");
  this.lineNo++;
  var line2 = AlphabetParser.strip_comments(line);
  if (line2 === "")  return;
  var m;
  if ((m = this.HEADER_RE.exec(line2)) != null) {
    this.parse_header(AlphabetParser.parse_name(m[1]), m[2]);
  } else if ((m = this.CORE_PAIR_RE.exec(line2)) != null) {
    this.parse_symbol(m[1], {name: AlphabetParser.parse_name(m[2]), colour: m[3], complement: m[4]});
    this.parse_symbol(m[4], {name: AlphabetParser.parse_name(m[5]), colour: m[6], complement: m[1]});
  } else if ((m = this.CORE_SINGLE_RE.exec(line2)) != null) {
    this.parse_symbol(m[1], {name: AlphabetParser.parse_name(m[2]), colour: m[3]});
  } else if ((m = this.AMBIG_RE.exec(line2)) != null) {
    this.parse_symbol(m[1], {name: AlphabetParser.parse_name(m[2]), colour: m[3], comprise: m[4]});
  } else {
    this._parse_error("Alphabet line does not match a known pattern!");
  }
};

//
// Private method
// Merge the list of symbols for each set of comprising characters
// and create an alias list that includes all the ways of referring
// to the symbol.
//
AlphabetParser.prototype._parse_done_merge_aliases = function() {
  "use strict";
  var key, i, symbol;
  var mergedSymbolMap = {};
  // merge aliases
  for (key in this.symbolMap) {
    if (!this.symbolMap.hasOwnProperty(key)) continue;
    // sort the list of symbols with the same comprising core syms
    // puts in order: letters, numbers, symbols.
    var symbolList = this.symbolMap[key].sort(AlphabetParser.symbol_compare);
    // compiles a merged list of aliases and picks the best name and colour
    var symbolName = null;
    var symbolColour = 0;
    var symbolAliases = [];
    for (i = 0; i < symbolList.length; i++) {
      symbol = symbolList[i];
      if (symbolName == null) symbolName = symbol.name;
      if (symbolColour == 0) symbolColour = symbol.colour;
      symbolAliases.push(symbol.sym);
      symbolAliases.push.apply(symbolAliases, symbol.aliases);
    }
    // ensure the aliases are sorted in a standard way
    symbolAliases.sort(AlphabetUtil.sym_compare);
    // create and store a new merged symbol object
    symbol = {
      "sym": symbolList[0].sym,
      "aliases": symbolAliases,
      "name": symbolName,
      "colour": symbolColour,
      "complement": symbolList[0].complement,
      "comprise": symbolList[0].comprise
    };
    mergedSymbolMap[key] = symbol;
  }
  // update the symbol map
  this.symbolMap = mergedSymbolMap;
};

AlphabetParser.prototype._parse_done_create_wildcard_if_missing = function() {
  "use strict";
  var sym;
  var comprise = [];
  for (sym in this.coreSymbols) {
    if (!this.coreSymbols.hasOwnProperty(sym)) continue;
    comprise.push(sym);
  }
  comprise.sort(AlphabetUtil.sym_compare);
  var key = comprise.join("");
  if (this.symbolMap[key] == null) {
    var symbol = {
      "sym": '?',
      "aliases": ['?'],
      "name": null,
      "colour": 0,
      "complement": null,
      "comprise": comprise
    };
    this.symbolMap[key] = symbol;
    this.allSymbols['?'] = true;
  }
};

//
// Private method
// Assigns colours to core symbols which don't have one assigned.
//
AlphabetParser.prototype._parse_done_set_missing_colours = function() {
  "use strict";
  var i, j, key, symbol, rgb;
  var ncolours = 0;
  var unique_colours = {};
  for (key in this.symbolMap) {
    // skip over inherited properties
    if (!this.symbolMap.hasOwnProperty(key)) continue;
    symbol = this.symbolMap[key];
    if (symbol.colour != 0) {
      if (!unique_colours[symbol.colour]) ncolours++;
      unique_colours[symbol.colour] = true;
    } else if (key.length == 1) {
      ncolours++;
    }
  }
  // generate Lab colourspace versions of each assigned colour
  var uniques = [];
  for (rgb in unique_colours) {
    // skip over inherited properties
    if (!this.symbolMap.hasOwnProperty(rgb)) continue;
    rgb = 0 + rgb; // convert back to number (the "in" operator converts it to a string)
    // unique colours
    uniques.push(Alphabet.rgb2lab(rgb));
  }
  // generate evenly distributed colours in HSV colourspace and convert to RGB and Lab
  var sat = 1.0;
  var value = 0.4;
  var step = 360.0 / ncolours;
  var colours = [];
  for (i = 0; i < ncolours; i++) {
    rgb = Alphabet.hsv2rgb(step * i, sat, value);
    colours.push({"rgb": rgb, "lab": Alphabet.rgb2lab(rgb)});
  }
  // for each of the unique colours find the closest colour in the generated colours and remove it
  while (uniques.length > 0) {
    var best_dist = null;
    var best_i = null;
    var best_j = null;
    for (i = 0; i < uniques.length; i++) {
      for (j = 0; j < colours.length; j++) {
        var dist = Alphabet.lab_dist(uniques[i], colours[j].lab);
        if (best_dist == null || dist < best_dist) {
          best_dist = dist;
          best_i = i;
          best_j = j;
        }
      }
    }
    if (best_dist == null) {
      throw new Error("Somehow we ran out of colours?!");
    }
    uniques.splice(best_i, 1);
    colours.splice(best_j, 1);
  }
  // assign the colours
  for (key in this.symbolMap) {
    // skip over inherited properties
    if (!this.symbolMap.hasOwnProperty(key)) continue;
    symbol = this.symbolMap[key];
    if (key.length == 1 && symbol.colour == 0) {
      symbol.colour = colours.pop().rgb;
    }
  }
};

AlphabetParser.prototype._parse_done_check_like_alphabet = function() {
  "use strict";
  if (this.like != null) {
    var core = null;
    var comp = null;
    if (this.like == "RNA") {
      core = "ACGU";
    } else if (this.like == "DNA") {
      core = "ACGT";
      comp = "TGCA";
    } else if (this.like == "PROTEIN") {
      core = "ACDEFGHIKLMNPQRSTVWY";
    }
    if (core != null) {
      var i, sym_obj, comp1, comp2;
      for (i = 0; i < core.length; i++) {
        sym_obj = this.symbolMap[core.charAt(i)];
        if (sym_obj == null) {
          this._parse_error("Not like " + this.like + " - missing symbol '" + core.charAt(i) + "'");
          return;
        }
        if (sym_obj.comprise != null) {
          this._parse_error("Not like " + this.like + " - expected core symbol '" + core.charAt(i) + "'");
          return;
        }
        comp1 = null;
        if (comp != null) comp1 = this.symbolMap[comp.charAt(i)];
        comp2 = null;
        if (sym_obj.complement != null) comp2 = this.symbolMap[sym_obj.complement];
        if (comp1 && comp1 != comp2) {
          this._parse_error("Not like " + this.like + " - complement incorrect for '" + core.charAt(i) + "'");
          return;
        }
      }
    }
  }
};

//
// Finish the parsing
//
AlphabetParser.prototype.parse_done = function() {
  "use strict";
  if (this.error) throw new Error("An error occurred previously");
  var i, key, symbol;
  this._parse_done_merge_aliases();
  this._parse_done_create_wildcard_if_missing();
  this._parse_done_set_missing_colours();
  this._parse_done_check_like_alphabet();
  // create the final variables which will actually be used
  var symbols = [];
  for (key in this.symbolMap) {
    if (!this.symbolMap.hasOwnProperty(key)) continue;
    symbols.push(this.symbolMap[key]);
  }
  symbols.sort(AlphabetParser.symbol_compare);
  // 
  var obj = {
    "name": this.name,
    "ncore": this.ncore,
    "symbols": []
  };
  if (this.like != null) obj["like"] = this.like.toUpperCase();
  for (i = 0; i < symbols.length; i++) {
    var symbol = symbols[i];
    var symObj = {"symbol": symbol.sym};
    if (symbol.aliases.length > 0) symObj.aliases = symbol.aliases.join("");
    if (symbol.name != null) symObj.name = symbol.name;
    if (symbol.colour != 0) symObj.colour = ("000000" + symbol.colour.toString(16)).substr(-6).toUpperCase();
    if (symbol.complement != null) symObj.complement = symbol.complement;
    if (symbol.comprise != null) symObj.equals = symbol.comprise.join("");
    obj.symbols.push(symObj);
  }
  return obj;
};

AlphabetParser.prototype.process_blob = function (blob) {
  "use strict";
  if (blob.size > 10000) {
    if (typeof this.handler.error == "function") {
      this.handler.error(true, "This file is too large to be an alphabet definition.");
    }
    return;
  }
  // setup the reader
  var me = this; // so we can access 'this' inside the closure
  var reader = new FileReader();
  reader.onload = function(evt) {
    "use strict";
    var chunk, format, message;
    // process the loaded chunk
    chunk = new Uint8Array(reader.result);
    // check magic numbers
    if ((format = unusable_format(chunk, 40, blob.name)) != null) {
      switch (format.type) {
        case FileType.ENCODING:
          message = "This is encoded as " + format.name + " which is not parsable as a alphabet.";
          break;
        case FileType.BINARY:
        case FileType.COMPRESSED:
        default:
          message = "This is a " + format.name + " file which is not parsable as a alphabet.";
          break;
      }
      if (typeof me.handler.error == "function") me.handler.error(true, message);
      return;
    }
    // convert to lines
    var decoder, line, code, letter;
    decoder = new Utf8Decoder
    decoder.set_source(chunk, true);
    line = "";
    while ((code = decoder.next()) != null) {
      letter = String.fromCharCode(code);
      if (letter == '\n') {
        // return the line - first remove any leading or trailing carrage return chars
        // as I'm assuming that they're part of the newline
        if (line.charAt(0) == '\r') {
          line = line.slice(1);
        }
        if (line.charAt(line.length - 1) == '\r') {
          line = line.slice(0, -1);
        }
        try {
          me.parse_line(line);
        } catch (e) {
          if (e instanceof AlphabetError) {
            if (typeof me.handler.error == "function") me.handler.error(true, e.message);
            return;
          } else {
            throw e;
          }
        }
        line = "";
        continue;
      }
      // add a letter to the line
      line += letter;
    }
    try {
      if (line.length > 0) me.parse_line(line);
      var alphabet_data = me.parse_done();
      if (typeof me.handler.data == "function") me.handler.data(alphabet_data);
    } catch (e) {
      if (e instanceof AlphabetError) {
        if (typeof me.handler.error == "function") me.handler.error(true, e.message, e.reasons);
        return;
      } else {
        throw e;
      }
    }
  };
  // read the next chunk
  reader.readAsArrayBuffer(blob);
};

// 
// Utility function
// Removes comments not within a quoted string
//
AlphabetParser.strip_comments = function(text) {
  "use strict";
  var instr = false;
  var eschr = false;
  var i, c;
  outer:
  for (i = 0; i < text.length; i++) {
    c = text.charAt(i);
    if (instr) {
      if (eschr) { // skip one
        eschr = false;
      } else {
        switch (c) {
          case '\\':
            eschr = true;
            break;
          case '"':
            instr = false;
            break;
        }
      }
    } else {
      switch (c) {
        case '"':
          instr = true;
          break;
        case '#':
          break outer;
      }
    }
  }
  return text.substring(0, i);
};

// 
// Utility function
//
AlphabetParser.parse_name = function(token) {
  "use strict";
  if (token == null || token === "") return null;
  var name = JSON.parse(token);
  if (typeof name !== "string") throw new Error("token is not just a string!");
  return name;
};

//
// Utility function
//
AlphabetParser.symbol_compare = function(symbol1, symbol2) {
  var i, cmp;
  // compare comprise
  if (symbol1.comprise != null && symbol2.comprise != null) {
    if (symbol1.comprise.length != symbol2.comprise.length) {
      return symbol2.comprise.length - symbol1.comprise.length;
    }
    for (i = 0; i < symbol1.comprise.length; i++) {
      cmp = AlphabetUtil.sym_compare(symbol1.comprise[i], symbol2.comprise[i]);
      if (cmp != 0) return cmp;
    }
  } else if (symbol1.comprise != null) {
    return 1;
  } else if (symbol2.comprise != null) {
    return -1;
  }
  // compare sym
  cmp = AlphabetUtil.sym_compare(symbol1.sym, symbol2.sym);
  if (cmp != 0) return cmp;
  // compare complement
  if (symbol1.complement != null && symbol2.complement != null) {
    cmp = AlphabetUtil.sym_compare(symbol1.complement, symbol2.complement);
    if (cmp != 0) return cmp;
  } else if (symbol1.complement != null) {
    return 1;
  } else if (symbol2.complement != null) {
    return -1;
  }
  // compare aliases
  if (symbol1.aliases.length != symbol2.aliases.length) {
    return symbol2.aliases.length - symbol1.aliases.length;
  }
  for (i = 0; i < symbol1.aliases.length; i++) {
    cmp = AlphabetUtil.sym_compare(symbol1.aliases[i], symbol2.aliases[i]);
    if (cmp != 0) return cmp;
  }
  // compare name
  if (symbol1.name != null && symbol2.name != null) {
    return (symbol1.name < symbol2.name ? -1 : (symbol1.name == symbol2.name ? 0 : 1));
  } else if (symbol1.name != null) {
    return 1;
  } else if (symbol2.name != null) {
    return -1;
  }
  // compare colour
  return symbol2.colour - symbol1.colour;
};


