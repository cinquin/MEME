var BgCharType = {
  JUNK:         0, 
  WHITESPACE:   1<<0, 
  COMMENT:      1<<1,
  CHAIN:        1<<2,
  NUMBER_START: 1<<3,
  NUMBER_MID:   1<<4
};

var BgErrorType = {
  FORMAT: 1,
  ENCODING: 2,
  ALPHABET_MISMATCH: 3,
  JUNK: 5,
  DUPLICATE: 7,
  MISSING_CHAIN: 8,
  MISSING_LENGTH: 9,
  INCORRECT_SUM: 11,
  BAD_PROBABILITY: 12,
  MISSING_PROBABILITY: 13
};

var BgAbortError = function(message) {
  "use strict";
  this.message = message;
  this.stack = Error().stack;
};
BgAbortError.prototype = Object.create(Error.prototype);
BgAbortError.prototype.name = "BgAbortError";
BgAbortError.prototype.constructor = BgAbortError;

var BgParserUtil = {};

BgParserUtil.add_map = function (map, letters, flags) {
  "use strict";
  var i;
  for (i = 0; i < letters.length; i++) {
    map[letters.charCodeAt(i)] |= flags;
  }
};

BgParserUtil.make_type_map = function () {
  "use strict";
  var i, map;
  map = new Int16Array(128);
  for (i = 0; i < map.length; i++) map[i] = BgCharType.JUNK;
  BgParserUtil.add_map(map, "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.*-?", BgCharType.CHAIN);
  BgParserUtil.add_map(map, "0123456789.", BgCharType.NUMBER_START);
  BgParserUtil.add_map(map, "0123456789.eE+-", BgCharType.NUMBER_MID);
  BgParserUtil.add_map(map, " \t\f\n\r", BgCharType.WHITESPACE);
  BgParserUtil.add_map(map, "#", BgCharType.COMMENT);
  return map;
};

BgParserUtil.compare_alphabet = function(syms, alph) {
  bad = [];
  missing = [];
  idxs = [];
  for (i = 0; i < syms.length; i++) {
    idxs[i] = alph.get_index(syms.charAt(i));
    if (idxs[i] == -1) bad.push(syms.charAt(i));
  }
  idxs.sort(function (a, b) { return a - b; });
  for (i = 0; i < idxs.length; i++) {
    if (idxs[i] != -1) break;
  }
  for (j = 0; i < idxs.length; i++, j++) {
    for (; idxs[i] > j; j++) {
      missing.push(alph.get_symbol(j));
    }
  }
  for (; j < alph.get_size_core(); j++) {
    missing.push(alph.get_symbol(j));
  }
  return {"bad": bad, "missing": missing};
};

BgParserUtil.TYPE_MAP = BgParserUtil.make_type_map();

var BgParser = function (handler) {
  "use strict";
  // Handler for callbacks.
  this.handler = handler;
  // Has the parser been told to stop parsing? We may be waiting on async processes so keep track.
  this.give_up = false;
  // Convert bytes into UTF-8 code points, track current line, column and byte offset.
  this.decoder = new LineDecoder();
  // The current parsing method, this changes as new input is seen.
  this.process = this._process_start;
  // Tracks the current known symbols.
  this.syms = "";
  // Converts a known symbol code into the number choosen to represent it. This number will equal the position in the syms string.
  this.sym_map = {};
  // Stores the text making up the symbol chain.
  this.chain = "";
  // Stores the text making up the probability of the symbol chain.
  this.prob = "";
  // Seen a chain of this length.
  this.chain_len = 1;
  // Store of all probabilities.
  this.background = new Float64Array(1);
  this.background[0] = -1;
};

BgParser.prototype._signal_stop = function () {
  if (typeof this.handler.progress == "function") this.handler.progress(1.0, this.syms);
  if (typeof this.handler.end == "function") this.handler.end(this.syms, this.background);
};

BgParser.prototype._error = function (offset, line, column, code, message, data) {
  if (window.console && console.log) {
    console.log("BG PARSER ERROR - offset: " + offset + "  line: " + line +
        "  column: " + column + "  code: " + code + "  message: " + message +
        "\ndata: " + JSON.stringify(data));
  }
  if (typeof this.handler.error == "function") {
    this.handler.error(offset, line, column, code, message, data);
  }
};

BgParser.prototype._error_here = function (code, message, data) {
  this._error(this.decoder.offset(), this.decoder.line(), this.decoder.column(),
      code, message, data);
};

BgParser.prototype._process_start = function (code, type) {
  "use strict";
  var error_code, error_message;
  // dispatch to the correct mode
  if ((type & BgCharType.COMMENT) != 0) {
    //nothing in line yet so we can just read the comment
    this.process = this._process_comment;
    return true;
  } else if ((type & BgCharType.WHITESPACE) != 0) {
    // whitespace at the start of the line can be safely ignored
    return true;
  } else if ((this.chain_len == 1 && (type & BgCharType.CHAIN) != 0) || this.sym_map[code] != null) {
    this.process = this._process_chain;
    this.chain = "";
    this.prob = "";
    return false;
  } else {
    // report error and treat as comment
    if ((type & BgCharType.CHAIN) != 0) { // not in alphabet
      this._error_here(BgErrorType.ALPHABET_MISMATCH,
          "Found a letter which does not match the previously established alphabet.",
          {"letter": String.fromCharCode(code)});
    } else { // junk
      this._error_here(BgErrorType.JUNK, 
          "Found a character that was not part of a state chain", {});
    }
    // skip until new line
    this.process = this._process_comment;
    return true;
  }
}

BgParser.prototype._process_comment = function (code, type) {
  // Blindly accept characters until we get to a new line
  if (this.decoder.column() == 0) {
    this.process = this._process_start;
    return false;
  }
  return true;
};

BgParser.prototype._process_chain = function (code, type) {
  "use strict";
  if ((type & BgCharType.WHITESPACE) != 0) {
    // white space, so end of the chain
    this.process = this._process_gap;
    return true;
  } else if ((type & BgCharType.CHAIN) != 0 && ((this.chain.length == 0 && this.chain_len == 1) || this.sym_map[code] != null)) {
    // extend the chain
    this.chain += String.fromCharCode(code);
    return true;
  } else {
    //error
    if ((type & BgCharType.CHAIN) != 0) {
      this._error_here(BgErrorType.ALPHABET_MISMATCH,
          "Found a symbol which does not exist in the previously established alphabet.",
          {"letter": String.fromCharCode(code)});
    } else {
      this._error_here(BgErrorType.JUNK, 
          "Bad character", {});
    }
    this.process = this._process_comment;
    return true;
  }
};

BgParser.prototype._process_gap = function (code, type) {
  "use strict";
  if (this.decoder.column() == 0) {
    this._error_here(BgErrorType.MISSING_PROBABILITY,
        "Did not find a probability to go with the chain", {});
    this.process = this._process_start;
    return false;
  } else if ((type & BgCharType.WHITESPACE) != 0) {
    return true;
  } else if ((type & BgCharType.NUMBER_START) != 0) {
    this.process = this._process_number;
    return false;
  } else {
    // error
    this._error_here(BgErrorType.JUNK, 
        "Found a character that was not the start of a number", {});
    this.process = this._process_comment;
    return true;
  }
};

BgParser.prototype._idx_to_chain = function(idx) {
  var value = idx + 1;
  var chain = "";
  while (value > 0) {
    var sym_index = ((value - 1) % this.syms.length);
    value -= (sym_index + 1);
    value /= sym_index;
    chain += this.syms.charAt(sym_index);
  }
  return chain;
};

BgParser.prototype._chain_to_idx = function(chain) {
  var i, idx;
  idx = this.sym_map[chain.charCodeAt(0)] + 1;
  for (i = 1; i < chain.length; i++) {
    idx = idx * this.syms.length + (this.sym_map[chain.charCodeAt(i)] + 1);
  }
  idx--; // index counts from 0 not 1;
  return idx; 
};


BgParser.prototype._validate_and_store_record = function() {
  "use strict";
  var i, index, prob, new_background;
  // validate the chain
  if (this.chain.length > this.chain_len) {
    // check that all previous chains of smaller lengths have been declared
    var start = this._chain_to_idx(string_mult(this.syms.charAt(0), this.chain_len));
    var end = this._chain_to_idx(string_mult(this.syms.slice(-1), this.chain_len));
    var missing_list = [];
    var sum = 0;
    for (i = start; i <= end; i++) {
      if (this.background[i] == -1) {
        // missing entry!
        missing_list.push(this._idx_to_chain(i));
      } else {
        sum += this.background[i];
      }
    }
    // report problems
    if (missing_list.length > 0) {
      // there are some missing entries for chain_len
      this._error_here(BgErrorType.MISSING_CHAIN,
          "Missing " + missing_list.length + " chain value(s) of length " +
          this.chain_len + ".", {missing: missing_list});
    } else if (Math.abs(sum - 1.0) > 0.1) {
      // there are no missing entries but they don't sum to 1
      this._error_here(BgErrorType.INCORRECT_SUM,
          "Summed probability of all length " + this.chain_len +
          " chains was " + sum + " when it should have been 1.0",
          {chain_length: this.chain_len, total: sum});
    }
    if (this.chain.length > (this.chain_len + 1)) {
      // we jumped by an increment of more than 1 so we are missing some
      this._error_here(BgErrorType.MISSING_LENGTH,
          "Missing chains for lengths " + (this.chain_len + 1) + " to " +
          (this.chain.length - 1) + ".",
          {start: (this.chain_len + 1), end: (this.chain.length - 1)});
    }
    // all core symbols should have been initilised by the first
    // chain of length 2 so we need to ensure they are sorted correctly
    if (this.chain_len == 1) {
      var bg = [];
      for (i = 0; i < this.syms.length; i++) {
        bg.push({sym: this.syms.charAt(i), prob: this.background[i]});
      }
      bg.sort(function (a, b) { return AlphabetUtil.sym_compare(a.sym, b.sym); });
      this.syms = "";
      for (i = 0; i < bg.length; i++) {
        this.syms += bg[i].sym;
        this.background[i] = bg[i].prob;
      }
    }
    // expand the background to fit
    var size = 0;
    for (i = 1; i <= this.chain.length; i++) size += Math.pow(this.syms.length, i);
    if (this.background.length < size) {
      new_background = new Float64Array(size);
      new_background.set(this.background);
      new_background.fill(-1, this.background.length);
      this.background = new_background;
    }
    this.chain_len = this.chain.length;
  }
  // Add unknown symbols to the alphabet.
  // Note that chains longer than 1 will not have got to this point if any of
  // the symbols were unknown.
  if (this.chain_len == 1 && this.chain.length == 1 && this.sym_map[this.chain.charCodeAt(0)] == null) {
    this.sym_map[this.chain.charCodeAt(0)] = this.syms.length;
    this.syms += this.chain;
    // expand the background to fit
    if (this.background.length < this.syms.length) {
      new_background = new Float64Array(this.background.length * 2);
      new_background.set(this.background);
      new_background.fill(-1, this.background.length);
      this.background = new_background;
    }
  }
  // convert chain to a position
  index = this._chain_to_idx(this.chain);
  // validate the probability
  if (this.prob == "") {
    this._error_here(BgErrorType.MISSING_PROBABILITY, "Did not find a probability to go with the chain", {});
  } else if (!/^(?:0(?:\.[0-9]+)?|[1-9][0-9]*(?:\.[0-9]+)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$/.test(this.prob)) {
    this._error_here(BgErrorType.BAD_PROBABILITY, "Number does not match expected format.", {});
    prob = 0;
  } else {
    prob = +this.prob;
    if (prob <= 0 || prob >= 1) {
      this._error_here(BgErrorType.BAD_PROBABILITY, "Number is not in expected range 0 < p < 1", {});
      prob = 0;
    }
  }
  // check for duplicates
  if (this.background[index] != -1) {
    this._error_here(BgErrorType.DUPLICATE,
        "Found a duplicated symbol chain " +
        this.chain + " = " + this.background[index], {});
  } else {
    // store the probability
    this.background[index] = prob;
  }
};

BgParser.prototype._process_number = function (code, type) {
  "use strict";
  var prob;
  if ((type & (BgCharType.NUMBER_START | BgCharType.NUMBER_MID)) != 0) {
    this.prob += String.fromCharCode(code);
    return true;
  } else if ((type & BgCharType.WHITESPACE) != 0) {
    this._validate_and_store_record();
    this.process = this._process_whitespace;
    return false;
  } else {
    this._error_here(BgErrorType.JUNK, "Expected a number.", {});
    this.process = this._process_comment;
    return true;
  }
};

BgParser.prototype._process_whitespace = function (code, type) {
  "use strict";
  if (this.decoder.column() == 0) {
    this.process = this._process_start;
    return false;
  } else if ((type & BgCharType.WHITESPACE) != 0) {
    return true;
  } else {
    this._error_here(BgErrorType.JUNK, "Expected whitespace", {});
    this.process = this._process_comment;
  }
};

BgParser.prototype._process_finish = function() {
  if (this.process == this._process_chain || 
      this.process == this._process_gap || 
      this.process == this._process_number) {
    this._validate_and_store_record();
  }
};

BgParser.prototype.process_blob = function (blob, offset, chunk_size) {
  "use strict";
  var reader, me;
  if (this.give_up) return;
  // set default parameter values
  if (typeof offset === "undefined") offset = 0;
  if (typeof chunk_size === "undefined") chunk_size = 1 << 12;
  // update the progress
  if (offset === 0 && typeof this.handler.begin == "function")
    this.handler.begin(blob.size);
  if (typeof this.handler.progress == "function")
    this.handler.progress(offset / blob.size, this.syms);
  // setup the reader
  me = this; // so we can access 'this' inside the closure
  reader = new FileReader();
  reader.onload = function(evt) {
    "use strict";
    var chunk, format, code, type, consumed, i;
    if (me.give_up) return;
    // process the loaded chunk
    chunk = new Uint8Array(reader.result);
    if (offset == 0) {
      if ((format = unusable_format(chunk, 40, blob.name)) != null) {
        me._error(0, 0, 0, BgErrorType.FORMAT,
            "The file format is not correct for a background file",
            {"format_type": format.type, "format_name": format.name});
        me._signal_stop();
        return;
      }
    }
    try {
      me.decoder.set_source(chunk, (offset + chunk_size) >= blob.size);
      while ((code = me.decoder.next()) != null) {
        type = (code < 128 ? BgParserUtil.TYPE_MAP[code] : BgCharType.JUNK);
        for (i = 0; i < 20; i++) {
          consumed = me.process(code, type);
          if (consumed) break;
        }
        if (i == 20) throw new Error("Infinite loop protection activated!");
      }
    } catch (e) {
      if (e instanceof EncodingError) {
        // report error and stop scan as something is wrong with the encoding
        me._error(e.offset, this.decoder.line(), this.decoder.column(),
            BgErrorType.ENCODING, e.message, {enc_code: e.code});
        me._signal_stop();
        return;
      } else if (e instanceof BgAbortError) {
        me._signal_stop();
        return;
      } else {
        throw e;
      }
    }
    if ((offset + chunk_size) >= blob.size) {
      me._process_finish();
      me._signal_stop();
    } else {
      // start loading the next chunk
      me.process_blob(blob, offset + chunk_size, chunk_size);
    }
  };
  // read the next chunk
  reader.readAsArrayBuffer(blob.slice(offset, offset + chunk_size));
};

BgParser.prototype.cancel = function () {
  this.give_up = true;
  this.handler = {};
};

var BgHandler = function(options) {
  this.reset();
  this.configure(options);
};

BgHandler.prototype.configure = function (options) {
  "use strict";
  if (typeof options != "object" || options == null) options = {};
  // configure file size
  if (typeof options.file_max == "number") {
    this.max_file_size = options.file_max;
  } else {
    this.max_file_size = 1 << 20;
  }
  // configure expected alphabet
  if (options.alphabet != null && options.alphabet instanceof Alphabet) {
    this.expected_alphabet = options.alphabet;
  } else {
    this.expected_alphabet = null;
  }
  this._calc_alph();
};

BgHandler.prototype.reset = function () {
  // have the file details changed?
  this.updated = false;
  // the part of the file processed
  this.fraction = 0;
  // bg details
  this.file_size = 0;
  this.syms = "";
  this.background = null;
  // these indicate when not UTF-8
  this.unusable_format_type = 0;
  this.unusable_format_name = null;
  this.encoding_error = null;
  // problems
  this.syntax_errors = new FileFaults();
  this.mismatches = new FileFaults();
  this.duplicates = new FileFaults();
  this.missing_chains = [];
  this.missing_ranges = [];
  this.bad_sums = [];
  // alphabet problems (updated when symbols change and expected_alphabet changes)
  this.bad_syms = [];
  this.missing_syms = [];
  this.guess_alphabet = null;
};

BgHandler.prototype._calc_alph = function() {
  "use strict";
  var info;
  if (this.expected_alphabet != null) {
    info = BgParserUtil.compare_alphabet(this.syms, this.expected_alphabet);
    this.bad_syms = info.bad;
    this.missing_syms = info.missing;
  } else {
    this.bad_syms = [];
    this.missing_syms = [];
  }
  if (this.expected_alphabet != null && this.bad_syms.length == 0 && this.missing_syms.length == 0) {
    this.guess_alphabet = this.expected_alphabet;
  } else if (this.syms == "ACGT") {
    this.guess_alphabet = AlphStd.DNA;
  } else if (this.syms == "ACGU") {
    this.guess_alphabet = AlphStd.RNA;
  } else if (this.syms == "ACDEFGHIKLMNPQRSTVWY") {
    this.guess_alphabet = AlphStd.PROTEIN;
  }
};

BgHandler.prototype.summary = function () {
  var error, warning, messages, reason, reasons, letters, add;
  var help;
  // setup
  error = false;
  warning = false;
  messages = [];
  // create closure to add messages
  add = function(is_error, message, reasons) {
    "use strict";
    messages.push({"is_error": is_error, "message": message, "reasons": reasons});
    if (is_error) error = true;
    else warning = true;
  };
  // file size warning
  if (this.max_file_size != null && this.file_size > this.max_file_size) {
    add(false, "Large file", ["File is " + (Math.round(this.file_size / (1<<20) )) + "MB"])
  }
  // encoding or format error
  help = " - re-save as plain text; either Unicode UTF-8 (no Byte Order Mark) or ASCII";
  if (this.unusable_format_name != null) {
    switch (this.unusable_format_type) {
      case FileType.ENCODING:
        add(true, "Bad encoding \"" + this.unusable_format_name + "\"" + help);
        break;
      case FileType.BINARY:
        add(true, "Bad format \"" + this.unusable_format_name + "\"" + help);
        break;
      case FileType.COMPRESSED:
        add(true, "Bad format \"" + this.unusable_format_name + "\" - must be decompressed first");
        break;
    }
    return {"error": error, "warning": warning, "messages": messages};
  } else if (this.encoding_error != null) {
    // report anything else that indicates it's not UTF-8
    add(true, "Bad encoding" + help, [this.encoding_error]);
    return {"error": error, "warning": warning, "messages": messages};
  }
  // syntax errors
  if (this.syntax_errors.faults() > 0) {
    reason = "Incorrect syntax on " +
      (this.syntax_errors.line_count == 1 ? "line " : "lines ") + this.syntax_errors.lines_str();
    add(true, "Syntax errors", [reason]);
  }
  // report non-alpha characters
  if (this.mismatches.faults() > 0) {
    letters = this.mismatches.letters_sep();
    if (this.mismatches.faults() == 1) {
      reason = "Found a character " + (letters.length > 0 ? letters + " " : "") +
        "not in the detected alphabet on ";
    } else {
      reason = "Found " + this.mismatches.faults() + " characters " + 
        (letters.length > 0 ? "(including " + letters + ") " : "") + 
        "not in the detected alphabet on ";
    }
    reason += (this.mismatches.line_count == 1 ? "line " : "lines ") +
      this.mismatches.lines_str();
    add(true, "Alphabet is inconsistant", [reason]);
  }
  // report duplicated chains
  if (this.duplicates.faults() > 0) {
    if (this.duplicates.faults() == 1) {
      reason = "Found a duplicated entry on line " + this.duplicates.lines_str();
    } else {
      reason = "Found " + this.duplicates.line_count +
        " duplicated entries on lines " + this.duplicates.lines_str();
    }
    add(true, "Duplicated entries", [reason]);
  }
  // report missing chains
  reasons = []
  if (this.missing_chains.length > 0) {
    if (this.missing_chains.length == 1) {
      reason = "Did not find the entry for " + this.missing_chains[0];
    } else {
      reason = "Did not find " + this.missing_chains.length + " entries for ";
      var end = Math.min(10, this.missing_chains.length - 1);
      for (i = 0; i < end; i++) {
        if (i != 0) reason += ", ";
        reason += this.missing_chains[i];
      }
      if (i == (this.missing_chains.length - 1)) {
        reason += " and " + this.missing_chains[i];
      } else {
        reason += ", ...";
      }
    }
    reasons.push(reason);
  }
  if (this.missing_ranges.length > 0) {
    if (this.missing_ranges.length == 1) {
      reason = "Did not find any entries of length ";
      if (this.missing_ranges[0].start == this.missing_ranges[0].end) {
        reason += this.missing_ranges[0].start;
      } else {
        reason += this.missing_ranges[0].start + " to " + this.missing_ranges[0].end;
      }
    } else {
      reason = "Did not find any entries of length ";
      end = Math.min(5, this.missing_ranges.length - 1);
      for (i = 0; i < end; i++) {
        if (i != 0) reason += ", ";
        if (this.missing_ranges[i].start == this.missing_ranges[i].end) {
          reason += this.missing_ranges[i].start;
        } else {
          reason += this.missing_ranges[i].start + " to " + this.missing_ranges[i].end;
        }
      }
      if (i == (this.missing_ranges.length - 1)) {
        reason += " or ";
        if (this.missing_ranges[i].start == this.missing_ranges[i].end) {
          reason += this.missing_ranges[i].start;
        } else {
          reason += this.missing_ranges[i].start + " to " + this.missing_ranges[i].end;
        }
      } else {
        reason += ", ...";
      }
    }
    reasons.push(reason);
  }
  if (reasons.length > 0) add(true, "Missing entries", reasons);
  // report incorrect sums
  if (this.bad_sums.length > 0) {
    if (this.bad_sums.length == 1) {
      reason = "The probabilities did not sum to 1 for the entries of length " + this.bad_sums[0];
    } else {
      reason = "The probabilities did not sum to 1 for multiple entry lengths ";
      end = Math.min(10, this.bad_sums.length - 1);
      for (i = 0; i < end; i++) {
        if (i != 0) reason += ", ";
        reason += this.bad_sums[i];
      }
      if (i == (this.bad_sums.length - 1)) {
        reason += " and " + this.bad_sums[i];
      } else {
        reason += ", ...";
      }
    }
    add(true, "Probabilities should sum to 1", [reason]);
  }

  if (this.bad_syms.length > 0 || (this.fraction == 1 && this.missing_syms > 0)) {
    reason = "The background alphabet " + 
      (this.guess_alphabet != null 
       ? "seems to be " + this.guess_alphabet : "could not be determined (" + this.syms + ")") +
      " but " + this.expected_alphabet + " was expected";
    add(true, "Background is wrong alphabet", [reason])
  }

  return {"error": error, "warning": warning, "messages": messages, 
    "syms": this.syms, "guess_alphabet": this.guess_alphabet};
};

BgHandler.prototype.begin = function (file_size) {
  this.reset();
  this.file_size = file_size;
  this.updated = true;
};

BgHandler.prototype.end = function (syms, background) {
  this.updated = true;
  this.syms = syms;
  this.background = background;
};

BgHandler.prototype.progress = function (fraction, syms) {
  this.fraction = fraction;
  if (this.syms != syms) {
    this.syms = syms;
    this._calc_alph();
  }
};

BgHandler.prototype.error = function (offset, line, column, code, message, data) {
  switch (code) {
    case BgErrorType.FORMAT:
      this.unusable_format_type = data["format_type"];
      this.unusable_format_name = data["format_name"];
      break;
    case BgErrorType.ENCODING:
      this.encoding_error = data["enc_code"];
      break;
    case BgErrorType.JUNK:
    case BgErrorType.BAD_PROBABILITY:
    case BgErrorType.MISSING_PROBABILITY:
      this.syntax_errors.add(line);
      break;
    case BgErrorType.ALPHABET_MISMATCH:
      this.mismatches.add(line, data["letter"]);
      break;
    case BgErrorType.DUPLICATE:
      this.duplicates.add(line);
      break;
    case BgErrorType.MISSING_CHAIN:
      this.missing_chains.push.apply(this.missing_chains, data["missing"]);
      break;
    case BgErrorType.MISSING_LENGTH:
      this.missing_ranges.push({start: data.start, end: data.end});
      break;
    case BgErrorType.INCORRECT_SUM:
      this.bad_sums.push(data["chain_length"]);
      break;
  }
  this.updated = true;
};

