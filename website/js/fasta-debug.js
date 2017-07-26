var FCharType = {
  JUNK:         0,    // anything that does not fit in another category
  DISALLOWED:   1<<0, // any character that is not allowed in a background file (ie NUL)
  WHITESPACE:   1<<1, // any space, tab or newline
  NEWLINE:      1<<2, // any newline
  COMMENT:      1<<3, // a comment (';')
  CHEVRON:      1<<4, // a cheveron ('>')
  SEQUENCE:     1<<5, // a character that might be in a alphabet
  ALPHA:        1<<6  // a character that is in an allowed alphabet
};

var FastaCheckerUtil = {};

FastaCheckerUtil.add_map = function (map, letters, flags) {
  "use strict";
  var i;
  for (i = 0; i < letters.length; i++) map[letters.charCodeAt(i)] |= flags;
};

FastaCheckerUtil.make_map = function (alphabets) {
  "use strict";
  var a, alphabet, i, j, aliases, map, sym;
  map = new Int16Array(128);
  for (i = 0; i < map.length; i++) map[i] = FCharType.JUNK;
  map[0x0] = FCharType.DISALLOWED;
  FastaCheckerUtil.add_map(map, " \t\f\n\r", FCharType.WHITESPACE);
  FastaCheckerUtil.add_map(map, "\n\r", FCharType.NEWLINE);
  FastaCheckerUtil.add_map(map, ";", FCharType.COMMENT);
  FastaCheckerUtil.add_map(map, ">", FCharType.CHEVRON);
  FastaCheckerUtil.add_map(map, "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.*-?", FCharType.SEQUENCE);
  for (a = 0; a < alphabets.length; a++) {
    alphabet = alphabets[a];
    for (i = 0; i < alphabet.get_size_full(); i++) {
      sym = alphabet.get_symbol(i);
      map[sym.charCodeAt(0)] |= FCharType.ALPHA;
      if (alphabet.is_case_insensitive()) {
        if (sym >= 'A' && sym <= 'Z') {
          map[sym.toLowerCase().charCodeAt(0)] |= FCharType.ALPHA;
        } else if (sym >= 'a' && sym <= 'z') {
          map[sym.toUpperCase().charCodeAt(0)] |= FCharType.ALPHA;
        }
      }
      aliases = alphabet.get_aliases(i);
      for (j = 0; j < aliases.length; j++) {
        sym = aliases.charAt(j);
        map[sym.charCodeAt(0)] |= FCharType.ALPHA;
        if (alphabet.is_case_insensitive()) {
          if (sym >= 'A' && sym <= 'Z') {
            map[sym.toLowerCase().charCodeAt(0)] |= FCharType.ALPHA;
          } else if (sym >= 'a' && sym <= 'z') {
            map[sym.toUpperCase().charCodeAt(0)] |= FCharType.ALPHA;
          }
        }
      }
    }
  }
  return map;
};


//******************************************************************************
// Fasta Checker
//******************************************************************************
var FastaChecker = function (handler, alphabets) {
  "use strict";
  if (alphabets == null) alphabets = [AlphStd.RNA, AlphStd.DNA, AlphStd.PROTEIN];
  var i;
  // store a reference to the handler
  this.handler = handler;
  // store the list of alphabets
  this.alphabets = alphabets;
  // enable the MEME extension "WEIGHTS"
  this.enable_weights = (typeof this.handler.enable_weights == "boolean" ? this.handler.enable_weights : false);
  // allow sequence characters which are a wildcard and are commonly used for masking
  this.allow_sequence_masking = (typeof this.allow_sequence_masking == "boolean" ? this.allow_sequence_masking : true);
  // allow ambiguous sequence characters
  this.allow_ambigs = (typeof this.handler.allow_ambigs == "boolean" ? this.handler.allow_ambigs : true);
  if (this.allow_ambigs) this.allow_sequence_masking = true;
  // only allow the sequence to be represented with uppercase letters
  // this is ignored for alphabets that are case-sensitive
  this.only_uppercase = (typeof this.handler.only_uppercase == "boolean" ? this.handler.only_uppercase : false);
  for (i = 0; i < alphabets.length; i++) {
    if (!alphabets[i].is_case_insensitive()) {
      this.only_uppercase = false;
      break;
    }
  }
  // setup a charCode mapping for the allowed characters in the sequence
  this.type_map = FastaCheckerUtil.make_map(alphabets);
  // UTF-8 reader
  this.decoder = new LineDecoder();
  // current parsing function
  this.process = this._process_start;
  // other stuff
  this.line_byte_offset = 0;
  this.file_symbols = "";
  this.file_counts = new Uint32Array(128);
  this.file_counts.fill(0);
  this.seq_byte_offset = 0;
  this.seq_line = 0;
  this.seq_name = "";
  this.seq_desc_len = 0;
  this.seq_offset = 0;
  this.after_comment = null;
  this.weight = "";
  this.weight_byte_offset = 0;
  this.weight_column = 0;
  this.weight_count = 0;
  // abort flag
  this.give_up = false;
};

// 
// Reports the sequence name, description length and sequence length to the handler.
//
FastaChecker.prototype._report_seq_info = function () {
  "use strict";
  // report
  if (typeof this.handler.info_seq == "function") {
    this.handler.info_seq(this.seq_byte_offset, this.seq_line, this.seq_name,
        this.seq_desc_len, this.seq_offset);
  }
};

// 
// Resets the sequence information to the defaults
//
FastaChecker.prototype._reset_seq_info = function() {
  this.seq_name = "";
  this.seq_desc_len = 0;
  this.seq_offset = 0;
  this.seq_byte_offset = this.decoder.offset();
  this.seq_line = this.decoder.line();
};

FastaChecker.prototype._process_start = function (code, type) {
  "use strict";
  if (this.decoder.column() == 0 && (type & FCharType.CHEVRON) != 0) {
    this._reset_seq_info();
    this.process = this._process_name;
  } else if (this.decoder.column() == 0 && (type & FCharType.COMMENT) != 0) {
    this.after_comment = this._process_start;
    this.process = this._process_comment;
    if (typeof this.handler.info_comment == "function") {
      this.handler.info_comment(this.decoder.offset(), this.decoder.line());
    }
  } else if (!((type & FCharType.WHITESPACE) != 0)) {
    // text before sequence
    if ((type & FCharType.CHEVRON) != 0) {
      // possible misplaced sequence start
      if (typeof this.handler.warn_chevron == "function") {
        this.handler.warn_chevron(this.decoder.offset(), this.decoder.line(),
            this.decoder.column());
      }
    }
    if (typeof this.handler.warn_junk == "function") {
      this.handler.warn_junk(this.decoder.offset(), this.decoder.line(),
          this.decoder.column(), String.fromCharCode(code));
    }
  }
  return true;
};

FastaChecker.prototype._process_comment = function (code, type) {
  if (this.decoder.column() == 0) {
    this.process = this.after_comment;
    return false;
  } else if ((type & FCharType.DISALLOWED) != 0) {
    if (typeof this.handler.warn_junk == "function") {
      this.handler.warn_junk(this.decoder.offset(), this.decoder.line(),
          this.decoder.column(), String.fromCharCode(code));
    }
  }
  return true;
};

FastaChecker.prototype._process_name = function (code, type) {
  "use strict";
  if ((type & FCharType.WHITESPACE) != 0) {
    // end of the sequence name
    if (this.enable_weights && this.seq_name == "WEIGHTS") {
      this.process = this._process_weights;
      this.weight = "";
      this.weight_byte_offset = this.decoder.offset() + 1;
      this.weight_column = this.decoder.column() + 1;
      this.weight_count = 0;
    } else {
      this.process = this._process_description;
      this.seq_desc_len = 0;
    }
    return false;
  } else if (!((type & FCharType.DISALLOWED) != 0)) {
    if ((type & FCharType.CHEVRON) != 0) {
      // possible misplaced sequence start
      if (typeof this.handler.warn_chevron == "function") {
        this.handler.warn_chevron(this.decoder.offset(), this.decoder.line(),
            this.decoder.column());
      }
    }
    this.seq_name += String.fromCharCode(code);
  } else {
    if (typeof this.handler.warn_junk == "function") {
      this.handler.warn_junk(this.decoder.offset(), this.decoder.line(),
          this.decoder.column(), String.fromCharCode(code));
    }
  }
  return true;
};

FastaChecker.prototype._process_description = function (code, type) {
  "use strict";
  if ((type & FCharType.NEWLINE) != 0) {
    // end of the description
    this.process = this._process_sequence;
  } else if (!((type & FCharType.DISALLOWED) != 0)) {
    this.seq_desc_len++;
    if ((type & FCharType.CHEVRON) != 0) {
      // possible misplaced sequence start
      if (typeof this.handler.warn_chevron == "function") {
        this.handler.warn_chevron(this.decoder.offset(), this.decoder.line(),
            this.decoder.column());
      }
    }
  } else {
    if (typeof this.handler.warn_junk == "function") {
      this.handler.warn_junk(this.decoder.offset(), this.decoder.line(),
          this.decoder.column(), String.fromCharCode(code));
    }
  }
  return true;
};

FastaChecker.prototype._process_weights = function (code, type) {
  "use strict";
  if ((type & FCharType.WHITESPACE) != 0) {
    if (this.weight.length > 0) {
      var is_weight = /^[10](?:\.\d+)?$/;
      var weight = parseFloat(this.weight);
      if (!is_weight.test(this.weight) || !isFinite(weight) || weight > 1 || weight <= 0) {
        if (typeof this.handler.warn_weight == "function") {
          this.handler.warn_weight(this.weight_byte_offset, this.decoder.line(), this.weight_column, this.weight);
        }
      }
      this.weight_count++;
    }
    this.weight = "";
    this.weight_byte_offset = this.byte_offst + 1;
    this.weight_column = this.decoder.column() + 1;
    if ((type & FCharType.NEWLINE) != 0) {
      if (typeof this.handler.info_weights == "function") {
        this.handler.info_weights(this.seq_byte_offset, this.seq_line, this.weight_count);
      }
      this.process = this._process_post_weights;
    }
  } else if (!((type & FCharType.DISALLOWED) != 0)) {
    this.weight += String.fromCharCode(code);
  } else {
    if (typeof this.handler.warn_junk == "function") {
      this.handler.warn_junk(this.decoder.offset(), this.decoder.line(),
          this.decoder.column(), String.fromCharCode(code));
    }
  }
  return true;
};

FastaChecker.prototype._process_post_weights = function (code, type) {
  "use strict";
  if (this.decoder.column() == 0 && (type & FCharType.CHEVRON) != 0) {
    // new sequence
    this._reset_seq_info();
    this.process = this._process_name;
  } else if (this.decoder.column() == 0 && (type & FCharType.COMMENT) != 0) {
    // comment
    this.after_comment = this._process_post_weights;
    this.process = this._process_comment;
    if (typeof this.handler.info_comment == "function") {
      this.handler.info_comment(this.decoder.offset(), this.decoder.line());
    }
  } else if (!((type & FCharType.WHITESPACE) != 0)) {
    // text before sequence
    if ((type & FCharType.CHEVRON) != 0) {
      // possible misplaced sequence start
      if (typeof this.handler.warn_chevron == "function") {
        this.handler.warn_chevron(this.decoder.offset(), this.decoder.line(),
            this.decoder.column());
      }
    }
    if (typeof this.handler.warn_junk == "function") {
      this.handler.warn_junk(this.decoder.offset(), this.decoder.line(),
          this.decoder.column(), String.fromCharCode(code));
    }
  }
  return true;
};

//
// This is used when reading sequence data and hence only valid sequence
// data or the start of a new sequence is allowed.
//
FastaChecker.prototype._process_sequence = function (code, type) {
  "use strict";
  var i;
  if (this.decoder.column() == 0) {
    if ((type & FCharType.CHEVRON) != 0) {
      // new sequence
      this._report_seq_info();
      this._reset_seq_info();
      this.process = this._process_name;
      return true;
    } else if ((type & FCharType.COMMENT) != 0) {
      // comment
      this.after_comment = this._process_sequence;
      this.process = this._process_comment;
      if (typeof this.handler.info_comment == "function") {
        this.handler.info_comment(this.decoder.offset(), this.decoder.line());
      }
      return true;
    }
  }

  // track potential alphabet symbols 
  if ((type & FCharType.SEQUENCE) != 0) {
    if ((this.file_counts[code] += 1) == 1) {
      this.file_symbols += String.fromCharCode(code);
    }
  }

  if ((type & FCharType.ALPHA) != 0) {
    this.seq_offset++;
  } else if (!((type & FCharType.WHITESPACE) != 0)) {
    // non whitespace
    if ((type & FCharType.CHEVRON) != 0) {
      // possible misplaced sequence start
      if (typeof this.handler.warn_chevron == "function") {
        this.handler.warn_chevron(this.decoder.offset(), this.decoder.line(),
            this.decoder.column());
      }
    }
    if (typeof this.handler.warn_seq == "function") {
      this.handler.warn_seq(this.decoder.offset(), this.decoder.line(),
          this.decoder.column(), this.seq_offset, this.seq_line, this.seq_name,
          type, String.fromCharCode(code));
    }
    this.seq_offset++;
  }
  return true;
};

// When we're done, call the approprate functions on the handler
FastaChecker.prototype._signal_stop = function() {
  if (typeof this.handler.progress == "function") this.handler.progress(1.0);
  if (typeof this.handler.end == "function") this.handler.end();
};

//******************************************************************************
// Public functions
//******************************************************************************

FastaChecker.prototype.process_blob = function (blob, offset, chunk_size) {
  "use strict";
  var reader, me, alphabet;
  if (this.give_up) return;
  me = this; // so we can access 'this' inside the closure
  // set default parameter values
  if (typeof offset === "undefined") offset = 0;
  if (typeof chunk_size === "undefined") chunk_size = 1 << 10;
  // update the progress
  if (offset === 0) {
    if (typeof this.handler.begin == "function") {
      this.handler.begin(
          blob.size, 
          function () { return me.guess_alphabets(); }
      );
    }
  }
  if (typeof this.handler.progress == "function") {
    this.handler.progress(offset / blob.size);
  }
  // setup the reader
  reader = new FileReader();
  reader.onload = function(evt) {
    "use strict";
    var i, chunk, format, code, type, consumed;
    if (me.give_up) return;
    // process the loaded chunk
    chunk = new Uint8Array(reader.result);
    if (offset == 0) { // check for common (but unreadable) file types
      if ((format = unusable_format(chunk, 400, blob.name)) != null) {
        // report error and stop scan as we don't have a chance of understanding this file
        if (typeof me.handler.error_format == "function") 
          me.handler.error_format(format.type, format.name);
        me._signal_stop();
        return;
      }
    }
    try {
      // validate the loaded chunk
      me.decoder.set_source(chunk, (offset + chunk_size) >= blob.size);
      while ((code = me.decoder.next()) != null) {
        type = code < 128 ? me.type_map[code] : FastaCheckerUtil.JUNK;
        for (i = 0; i < 20; i++) {
          consumed = me.process(code, type);
          if (consumed) break;
        }
        if (i == 20) throw new Error("Infinite loop protection activated!");
      }
    } catch (e) {
      if (e instanceof EncodingError) {
        // report error and stop scan as something is wrong with the encoding
        if (typeof me.handler.error_encoding == "function")
          me.handler.error_encoding(e.offset, e.code, e.message);
        me._signal_stop();
        return;
      } else {
        throw e;
      }
    }
    if ((offset + chunk_size) >= blob.size) {
      if (me.process === me._process_comment) me.process = me.after_comment;
      if (me.process === me._process_name) {
        if (me.enable_weights && me.seq_name == "WEIGHTS") {
          if (typeof me.handler.info_weights == "function") {
            me.handler.info_weights(me.seq_byte_offset, me.seq_line, this.weight_count);
          }
          me.process = me._process_post_weights;
        }
      }
      if (me.process === me._process_name || 
          me.process === me._process_description ||
          me.process === me._process_sequence) {
        me._report_seq_info();
      }
      me._signal_stop();
    } else {
      // start loading the next chunk
      me.process_blob(blob, offset + chunk_size, chunk_size);
    }
  };
  // read the next chunk
  reader.readAsArrayBuffer(blob.slice(offset, offset + chunk_size));
};

FastaChecker.prototype.cancel = function () {
  "use strict";
  this.give_up = true;
  this.handler = {};
};

FastaChecker.prototype.guess_alphabets = function() {
  "use strict";
  return AlphabetUtil.guess(this.alphabets, this.file_symbols, this.file_counts, true);
};


//******************************************************************************
// Fasta Handler
//******************************************************************************
var FastaHandler = function (options) {
  this.configure(options);
  this.reset();
};

FastaHandler.prototype.configure = function (options) {
  "use strict";
  var alphabets, alphabet_name;
  if (typeof options != "object" || options == null) options = {};
  // configure file size
  if (typeof options.file_max == "number") {
    this.max_file_size = options.file_max;
  } else {
    this.max_file_size = null; // 20 MB
  }
  // enable MEME weights extension?
  if (typeof options.weights == "boolean") {
    this.enable_weights = options.weights;
  } else {
    this.enable_weights = true;
  }
  // allow masking characters like N in DNA?
  if (typeof options.mask == "boolean") {
    this.allow_sequence_masking = options.masks;
  } else {
    this.allow_sequence_masking = true;
  }
  // allow other ambiguous characters?
  if (typeof options.ambigs == "boolean") {
    this.allow_ambigs = options.ambigs;
  } else {
    this.allow_ambigs = true;
  }
  // allow gap characters?
  if (typeof options.gaps == "boolean") {
    this.allow_gaps = options.gaps;
  } else {
    this.allow_gaps = false;
  }
  // require all sequence to be uppercase (can help with debugging)?
  if (typeof options.uppercase == "boolean") {
    this.only_uppercase = options.uppercase;
  } else {
    this.only_uppercase = false;
  }
  // specify a maximum name length?
  if (typeof options.max_name_len == "number" && options.max_name_len >= 1) {
    this.max_name_length = options.max_name_len;
  } else {
    this.max_name_length = null;
  }
  // specify a maximum description length?
  if (typeof options.max_desc_len == "number" && options.max_desc_len >= 1) {
    this.max_desc_length = options.max_desc_len;
  } else {
    this.max_desc_length = null;
  }
  // specify a minimum sequence length?
  if (typeof options.min_seq_len == "number" && options.min_seq_len >= 0) {
    this.min_seq_length = options.min_seq_len;
  } else {
    this.min_seq_length = null;
  }
  // specify a maximum sequence length?
  if (typeof options.max_seq_len == "number" && options.max_seq_len >= 1) {
    this.max_seq_length = options.max_seq_len;
  } else {
    this.max_seq_length = null;
  }
  // check that the minimum is not larger than the maximum
  if (this.min_seq_length != null && this.max_seq_length != null &&
      this.min_seq_length > this.max_seq_length) {
    throw new Error("Option min_seq must not be larger than option max_seq.");
  }
  // specifiy the maximum number of sequences
  if (typeof options.max_seq_count == "number" && options.max_seq_count >= 1) {
    this.max_seq_count = options.max_seq_count;
  } else {
    this.max_seq_count = null;
  }
  // specify the maximum total sequence length
  if (typeof options.max_seq_total == "number" && options.max_seq_total >= 1) {
    this.max_seq_total = options.max_seq_total;
  } else {
    this.max_seq_total = null;
  }
};

FastaHandler.prototype.reset = function () {
  // have the file details changed?
  this.updated = false;
  // what is the alphabet
  this.fn_alphabets = function() { return []; };
  // the part of the file processed
  this.fraction = 0;
  // fasta details
  this.file_size = 0;
  this.file_symbols = "";
  // these indicate when not UTF-8
  this.unusable_format_type = 0;
  this.unusable_format_name = null;
  this.encoding_error = null;
  // the lookup for sequence IDs
  this.name_lookup = {};
  // the count of sequences
  this.sequence_count = 0;
  this.sequence_total = 0;
  // the alphabet of the sequences
  this.alphabet = AlphType.UNKNOWN;
  // keep track of problems found
  this.missing_name = new FileFaults();
  this.long_name = new FileFaults();
  this.duplicate_name = new FileFaults();
  this.long_description = new FileFaults();
  this.short_sequence = new FileFaults();
  this.long_sequence = new FileFaults();
  this.comment = new FileFaults();
  this.junk = new FileFaults();
  this.bad_weight = new FileFaults();
  this.chevron = new FileFaults();
  // problems within the sequence data
  this.seq_errors = 0;// total count
  this.seq_gap = new FileFaults(); // gap characters
  this.seq_ambig = new FileFaults(); // ambiguous characters
  this.seq_mask = new FileFaults(); // masking character - Eg. N for DNA
  this.seq_lc = new FileFaults(); // lower case
  this.seq_non_alpha = new FileFaults();
  this.count_seqs_with_error = 0;
  this.last_seq_error_line = -1;
};

FastaHandler.prototype.summary = function () {
  "use strict";
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
  // alphabet warning
  if (this.sequence_count === 0) {
    add(true, "No sequences found");
  }
  // junk
  if (this.junk.faults() > 0) {
    reason = "Found some text that is not part of any sequence or weights on " +
      (this.junk.line_count == 1 ? "line " : "lines ") + this.junk.lines_str();
    add(true, "Junk text found", [reason]);
  }
  // weights
  if (this.bad_weight.faults() > 0) {
    if (this.bad_weight.faults() == 1) {
      reason = "Found a sequence weight that was not in the correct " +
        "range (0 < weight \u2264 1) on line " + this.bad_weight.lines_str();
    } else {
      reason = "Found " + this.bad_weight.faults() + " sequences weights that " +
        "were not in the correct range (0 < weight \u2264 1) on " +
        (this.bad_weight.line_count == 1 ? "line " : "lines ") +
        this.bad_weight.lines_str();
    }
    add(true, "Sequence weights out of range", [reason]);
  }
  // missing sequence identifiers
  if (this.missing_name.faults() > 0) {
    if (this.missing_name.faults() == 1) {
      reason = "Missed a sequence identifier on line ";
    } else {
      reason = "Missed " + this.missing_name.faults() +
        " sequence identifiers on lines ";
    }
    reason += this.missing_name.lines_str();
    add(true, "Sequence identifier missing - all sequences must have identifiers", [reason]);
  }
  // duplicated sequence identifiers
  if (this.duplicate_name.faults() > 0) {
    if (this.duplicate_name.faults() == 1) {
      reason = "Found a duplicate sequence identifier on line ";
    } else {
      reason = "Found " + this.duplicate_name.faults() + " duplicated sequence identifiers on lines ";
    }
    reason += this.duplicate_name.lines_str();
    add(true, "Sequence identifier duplicated - identifiers must be unique", [reason]);
  }
  // sequence length
  if (this.short_sequence.faults() > 0 || this.long_sequence.faults() > 0) {
    reasons = [];
    if (this.short_sequence.faults() > 0) {
      if (this.short_sequence.faults() == 1) {
        reason = "Found a sequence shorter than the minimum allowed length " + 
          "of " + this.min_seq_length + " on line ";
      } else {
        reason = "Found " + this.short_sequence.faults() + " sequences shorter " +
          "than the minimum allowed length of " + this.min_seq_length + " on lines ";
      }
      reason += this.short_sequence.lines_str();
      reasons.push(reason);
    }
    if (this.long_sequence.faults() > 0) {
      if (this.long_sequence.faults() == 1) {
        reason = "Found a sequence longer than the maximum allowed length " + 
          "of " + this.max_seq_length + " on line ";
      } else {
        reason = "Found " + this.long_sequence.faults() + " sequences longer " +
          "than the maximum allowed length of " + this.max_seq_length + " on lines ";
      }
      reason += this.long_sequence.lines_str();
      reasons.push(reason);
    }
    add(true, "Sequence length out of bounds", reasons);
  }
  // sequence data problems
  if (this.seq_errors > 0) {
    reasons = [];
    // report total errors
    reasons.push("Total of " + this.seq_errors + " bad characters in " +
        this.count_seqs_with_error + " sequences");
    // report gap characters
    if (this.seq_gap.faults() > 0) {
      if (this.seq_gap.faults() == 1) {
        reason = "Found a disallowed gap character " + 
          this.seq_gap.letters_sep() + " on ";
      } else {
        reason = "Found " + this.seq_gap.faults() + " disallowed gap " +
          "characters " + this.seq_gap.letters_sep() + " on ";
      }
      reason += (this.seq_gap.line_count == 1 ? "line " : "lines ") +
        this.seq_gap.lines_str();
      reasons.push(reason);
    }
    // report ambiguous characters
    if (this.seq_ambig.faults() > 0) {
      if (this.seq_ambig.faults() == 1) {
        reason = "Found a disallowed ambiguous character " + 
          this.seq_ambig.letters_sep() + " on ";
      } else {
        reason = "Found " + this.seq_ambig.faults() + " disallowed ambiguous " +
          "characters " + this.seq_ambig.letters_sep() + " on ";
      }
      reason += (this.seq_ambig.line_count == 1 ? "line " : "lines ") +
        this.seq_ambig.lines_str();
      reasons.push(reason);
    }
    // report masking characters
    if (this.seq_mask.faults() > 0) {
      if (this.seq_mask.faults() == 1) {
        reason = "Found a disallowed masking character " + 
          this.seq_mask.letters_sep() + " on ";
      } else {
        reason = "Found " + this.seq_mask.faults() + " disallowed masking " +
          "characters " + this.seq_mask.letters_sep() + " on ";
      }
      reason += (this.seq_mask.line_count == 1 ? "line " : "lines ") +
        this.seq_mask.lines_str();
      reasons.push(reason);
    }
    // report lowercase characters
    if (this.seq_lc.faults() > 0) {
      if (this.seq_lc.faults() == 1) {
        reason = "Found a disallowed lowercase character " + 
          this.seq_lc.letters_sep() + " on ";
      } else {
        reason = "Found " + this.seq_lc.faults() + " disallowed lowercase " +
          "characters " + this.seq_lc.letters_sep() + " on ";
      }
      reason += (this.seq_lc.line_count == 1 ? "line " : "lines ") +
        this.seq_lc.lines_str();
      reasons.push(reason);
    }
    // report non-alpha characters
    if (this.seq_non_alpha.faults() > 0) {
      letters = this.seq_non_alpha.letters_sep();
      if (this.seq_non_alpha.faults() == 1) {
        reason = "Found a character " + (letters.length > 0 ? letters + " " : "") +
          "not in the alphabet on ";
      } else {
        reason = "Found " + this.seq_non_alpha.faults() + " characters " + 
          (letters.length > 0 ? "(including " + letters + ") " : "") + 
          "not in the alphabet on ";
      }
      reason += (this.seq_non_alpha.line_count == 1 ? "line " : "lines ") +
        this.seq_non_alpha.lines_str();
      reasons.push(reason);
    }
    add(true, "Sequence contains bad characters", reasons);
  }
  // report comments
  if (this.comment.faults() > 0) {
    if (this.comment.faults() == 1) {
      reason = "Found a comment on line ";
    } else {
      reason = "Found " + this.comment.faults() + " comments on lines ";
    }
    reason += this.comment.lines_str();
    add(true, "Unsupported comments - comment lines must be removed", [reason]);
  }
  // report unusual chevrons
  if (this.chevron.faults() > 0) {
    if (this.chevron.faults() == 1) {
      reason = "Found a '>' character on line " + this.chevron.lines_str() +
        " which is not at the begining";
    } else {
      reason = "Found " + this.chevron.faults() + " '>' characters on lines " + 
        this.chevron.lines_str() + " where the '>' is not at the begining";
    }
    add(false, "Potential malformed sequence start", [reason]);
  }
  // report long sequence identifiers
  if (this.long_name.faults() > 0) {
    if (this.long_name.faults() == 1) {
      reason = "Found a sequence with a long identifier on line ";
    } else {
      reason = "Found " + this.long_name.faults() + 
        " sequences with long identifiers on lines ";
    }
    reason += this.long_name.lines_str();
    add(false, "Long sequence identifiers may cause problems", [reason]);

  }
  // report long sequence descriptions
  if (this.long_description.faults() > 0) {
    if (this.long_description.faults() == 1) {
      reason = "Found a sequence with a long description on line ";
    } else {
      reason = "Found " + this.long_description.faults() + 
        " sequences with long descriptions on lines ";
    }
    reason += this.long_description.lines_str();
    add(false, "Long sequence descriptions may cause problems", [reason]);
  }
  // report sequence count excess
  if (this.max_seq_count != null && this.sequence_count > this.max_seq_count) {
    add(true, "Too many sequences", ["Found " + this.sequence_count + 
        " sequences but this only accepts up to " + this.max_seq_count]);
  }
  // report sequence length total excess
  if (this.max_seq_total != null && this.sequence_total > this.max_seq_total) {
    add(true, "Combined sequence length exceeds maximum", 
        ["Found sequences with lengths totaling " + this.sequence_total +
        " but this only accepts a total length up to " + this.max_seq_total]);
  }
  // clear updated state
  this.updated = false;
  // return state
  return {"error": error, "warning": warning, "messages": messages};
};

FastaHandler.prototype.guess_alphabets = function () {
  return this.fn_alphabets();
}

// tracks the progress of reading the file
FastaHandler.prototype.progress = function (fraction) {
  "use strict";
  this.fraction = fraction;
};

// Reading of the file has begun
FastaHandler.prototype.begin = function (file_size, fn_alphabets) {
  "use strict";
  this.reset();
  this.file_size = file_size;
  this.fn_alphabets = fn_alphabets;
  this.updated = true;
};

// Reading of the file has finished (perhaps early due to an error)
FastaHandler.prototype.end = function () {
  "use strict";
  this.updated = true;
};

// Notes the existence of a sequence
FastaHandler.prototype.info_seq = function (offset, line, name, desc_len, len) {
  "use strict";
  var reason, below_max;
  // check if we're going to overstep the maximum sequence count
  this.sequence_count++;
  if (this.max_seq_count != null && this.sequence_count > this.max_seq_count) {
    this.updated = true;
  }
  // check if we're going to overstep the maximum total sequence length
  this.sequence_total += len;
  if (this.max_seq_total != null && this.sequence_total > this.max_seq_total) {
    this.updated = true;
  }
  // check for missing name
  if (name.length == 0) {
    this.missing_name.add(line);
    this.updated = true;
  } else {
    // check for long name
    if (this.max_name_length != null && name.length > this.max_name_length) {
      this.long_name.add(line);
      this.updated = true;
    }
    // check for duplicate name
    if (this.name_lookup[name]) {
      this.duplicate_name.add(line);
      this.updated = true;
    }
    // store name to check for duplicates later
    this.name_lookup[name] = true;
  }
  // check for long description
  if (this.max_desc_length != null && desc_len > this.max_desc_length) {
    this.long_description.add(line);
    this.updated = true;
  }
  // check sequence length
  if (this.min_seq_length != null && len < this.min_seq_length) {
    this.short_sequence.add(line);
    this.updated = true;
  } else if (this.max_seq_length != null && len > this.max_seq_length) {
    this.long_sequence.add(line);
    this.updated = true;
  }
};

// Notes the existence of a comment
FastaHandler.prototype.info_comment = function (offset, line) {
  "use strict";
  this.comment.add(line);
};

// Notes the existence of a set of weights (MEME extension to FASTA)
FastaHandler.prototype.info_weights = function (offset, line, num_weights) {
  "use strict";
};

// Warns about a character not associated with a sequence, comment or weights
FastaHandler.prototype.warn_junk = function (offset, line, column, character) {
  "use strict";
  this.junk.add(line);
  this.updated = true;
};

// Warns about a bad character in a sequence
FastaHandler.prototype.warn_seq = function (offset, line, column,
    seq_offset, seq_line, seq_name, type, character) {
  "use strict";
  var param, lines;
  this.seq_errors++;
  this.seq_non_alpha.add(line, character);
  // count the number of sequences with errors
  if (this.last_seq_error_line != seq_line) {
    this.count_seqs_with_error++;
    this.last_seq_error_line = seq_line;
  }
  this.updated = true;
};

// Warns about a bad weight (weights are a MEME extension to FASTA)
FastaHandler.prototype.warn_weight = function (offset, line, column, weight_text) {
  "use strict";
  this.bad_weight.add(line);
  this.updated = true;
};

// Warns about an oddly placed chevron character; 
// note that this same character may trigger a call to warn_junk, warn_seq or warn_weight
FastaHandler.prototype.warn_chevron = function (offset, line, column) {
  "use strict";
  this.chevron.add(line);
  this.updated = true;
};

// Parsing has stopped due to an unreadable file format
FastaHandler.prototype.error_format = function (type, name) {
  "use strict";
  this.unusable_format_type = type;
  this.unusable_format_name = name;
  this.updated = true;
};

// Parsing has stopped due to an error in the file encoding
FastaHandler.prototype.error_encoding = function (offset, type, reason) {
  "use strict";
  this.encoding_error = reason;
  this.updated = true;
};
