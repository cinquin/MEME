// Enumerations
var FileType = {ENCODING: 1, BINARY: 2, COMPRESSED: 3};
var EncErrorType = {GENERIC: 1, ILLEGAL_BYTE: 2, INCOMPLETE_CODE_UNIT: 3,
  OVERLONG_ENCODING: 4, INVALID_CODE_POINT: 5, NULL_BYTE: 6};

var AlphabetUtil = {};

//
// Utility function
//
AlphabetUtil.sym_compare = function(a, b) {
  "use strict";
  var a_letter = (a >= 'A' && a <= 'Z') || (a >= 'a' && a <= 'z');
  var b_letter = (b >= 'A' && b <= 'Z') || (b >= 'a' && b <= 'z');
  if (a_letter && b_letter) {
    return (a < b ? -1 : (a == b ? 0 : 1));
  } else if (a_letter) {
    return -1;
  } else if (b_letter) {
    return 1;
  }
  var a_num = (a >= '0' && a <= '9');
  var b_num = (b >= '0' && b <= '9');
  if (a_num && b_num) {
    return (a < b ? -1 : (a == b ? 0 : 1));
  } else if (a_num) {
    return -1;
  } else if (b_num) {
    return 1;
  }
  return (a < b ? -1 : (a == b ? 0 : 1));
};

AlphabetUtil.from_bitset = function(bitset) {
  var key, alphabet;
  var list = [];
  for (key in AlphStd) {
    if (!AlphStd.hasOwnProperty(key)) continue;
    alphabet = AlphStd[key];
    if (!(alphabet instanceof StandardAlphabet)) continue;
    if ((alphabet.get_code() & bitset) != 0) {
      list.push(alphabet);
    }
  }
  list.sort(function (alph_a, alph_b) { return alph_a.get_code() - alph_b.get_code(); });
  return list;
};

AlphabetUtil.to_bitset = function(list_of_alphabets) {
  var bitset = 0;
  var i, alphabet;
  for (i = 0; i < list_of_alphabets.length; i++) {
    alphabet = list_of_alphabets[i];
    if (!(alphabet instanceof StandardAlphabet)) continue;
    bitset |= alphabet.get_code();
  }
  return bitset;
};

AlphabetUtil.any_equal = function (list1, list2) {
  var i, j;
  if (list1 == null || list2 == null) return false;
  for (i = 0; i < list1.length; i++) {
    if (list1[i] == null) continue;
    for (j = 0; j < list2.length; j++) {
      if (list2[j] == null) continue;
      if (list1[i].equals(list2[j])) return true;
    }
  }
  return false;
};

AlphabetUtil.display = function (alphabets) {
  var i;
  display_name = "";
  for (i = 0; i < alphabets.length -1; i++) {
    if (i > 0) {
      display_name += ", ";
    }
    display_name += alphabets[i].get_alphabet_name();
  }
  if (alphabets.length > 1) {
    display_name += " or ";
  }
  if (alphabets.length >= 1) {
    display_name += alphabets[alphabets.length - 1].get_alphabet_name() + " alphabet";
  }
  return display_name;
};

AlphabetUtil.guess = function(alphabet_list, symbols, counts, only_best) {
  var i, j, sym, idx, count, alphabet, valid_alphabets, seen, seen_size, 
      prime_seen, prime_seen_size, prime_count, alt_count, other_count,
      total_count, score, prime_score;
  if (typeof only_best != "boolean") only_best = false;
  // filter the list so it only contains alphabets for with all the symbols are valid
  valid_alphabets = [];
  for (i = 0; i < alphabet_list.length; i++) {
    alphabet = alphabet_list[i];
    for (j = 0; j < symbols.length; j++) {
      if (alphabet.get_index(symbols.charAt(j)) == -1) break;
    }
    if (j == symbols.length) {
      valid_alphabets.push(alphabet);
    }
  }
  if (valid_alphabets.length == 0) return [];
  if (valid_alphabets.length == 1) return valid_alphabets;
  // calculate a score for the remaining alphabets
  var scored_alphabets = [];
  for (i = 0; i < valid_alphabets.length; i++) {
    alphabet = valid_alphabets[i];
    seen = {};
    seen_size = 0;
    prime_seen = {};
    prime_seen_size = 0;
    prime_count = 0;
    alt_count = 0;
    other_count = 0;
    unknown_count = 0;
    for (j = 0; j < symbols.length; j++) {
      sym = symbols.charAt(j);
      count = counts[symbols.charCodeAt(j)];
      idx = alphabet.get_index(sym);
      if (idx != -1) {
        if (idx < alphabet.get_size_core()) {
          if (!seen[idx]) {
            seen[idx] = true;
            seen_size++;
          }
          if (alphabet.is_prime_symbol(sym)) {
            if (!prime_seen[idx]) {
              prime_seen[idx] = true;
              prime_seen_size++;
            }
            prime_count += count;
          } else {
            alt_count += count;
          }
        } else if (idx == alphabet.get_size_core()) {
          if (alphabet.is_prime_symbol(sym)) {
            prime_count += count;
          } else {
            alt_count += count;
          }
        } else {
          other_count += count;
        }
      } else {
        unknown_count += count;
      }
    }
    total_count = prime_count + alt_count + other_count + unknown_count;
    score = ((prime_count + alt_count) / total_count) * (seen_size / alphabet.get_size_core());
    prime_score = (prime_count / total_count) * (prime_seen_size / alphabet.get_size_core());
    scored_alphabets.push({"alphabet": alphabet, "score": score, "prime_score": prime_score});
  }
  // sort the alphabets by score
  scored_alphabets.sort(
    function (a, b) { 
      if (a.score < b.score) {
        return 1;
      } else if (a.score > b.score) {
        return -1;
      } else if (a.prime_score < b.prime_score) {
        return 1;
      } else if (a.prime_score > b.prime_score) {
        return -1;
      } else {
        return 0;
      }
    }
  );
  var best_alphabets = [];
  best_alphabets[0] = scored_alphabets[0].alphabet;
  score = scored_alphabets[0].score;
  prime_score = scored_alphabets[0].prime_score;
  for (i = 1; i < scored_alphabets.length; i++) {
    if (only_best && (scored_alphabets[i].score != score || scored_alphabets[i].prime_score != prime_score)) break;
    best_alphabets[i] = scored_alphabets[i].alphabet;
  }
  return best_alphabets;
};

var AlphabetType = function(code, name, alphabets, custom, display_name, match_fn) {
  this.code = code;
  this.name = name;
  this.alphabets = alphabets;
  this.custom = (typeof custom == "boolean" ? custom : false);
  this.display_name = display_name;
  this.match_fn = match_fn;
};

AlphabetType.prototype.get_standard_alphabets = function() {
  return this.alphabets.slice();
};

AlphabetType.prototype.allow_custom_alphabets = function() {
  return this.custom;
};

AlphabetType.prototype.toString = function() {
  if (typeof this.display_name != "string") {
    this.display_name = AlphabetUtil.display(this.alphabets);
  }
  return this.display_name;
};

AlphabetType.prototype.matches = function(alphabet) {
  if (typeof this.match_fn == "function") {
    return this.match_fn(alphabet);
  }
  var i;
  for (i = 0; i < this.alphabets.length; i++) {
    if (this.alphabets[i].equals(alphabet)) return true;
  }
  return false;
};

AlphabetType.prototype.bitset = function() {
  return AlphabetUtil.to_bitset(this.alphabets);
};

var AlphType = {
  ANY_ALPHABET: 
    new AlphabetType(1, "ANY_ALPHABET", [AlphStd.RNA, AlphStd.DNA, AlphStd.PROTEIN],
        true,
        "any alphabet",
        function() {
          return true;
        }),
  COMPLEMENTABLE_ALPHABET: 
    new AlphabetType(2, "COMPLEMENTABLE_ALPHABET", [AlphStd.DNA], true,
        "complementable alphabet",
        function(alphabet) {
          return alphabet.has_complement();
        }),
  STANDARD_ALPHABET:
    new AlphabetType(3, "STANDARD_ALPHABET", [AlphStd.RNA, AlphStd.DNA, AlphStd.PROTEIN]),
  RNA_OR_DNA_ALPHABET:
    new AlphabetType(4, "RNA_OR_DNA_ALPHABET", [AlphStd.RNA, AlphStd.DNA]),
  DNA_OR_PROTEIN_ALPHABET:
    new AlphabetType(5, "DNA_OR_PROTEIN_ALPHABET", [AlphStd.DNA, AlphStd.PROTEIN]),
  RNA_ALPHABET:
    new AlphabetType(6, "RNA_ALPHABET", [AlphStd.RNA]),
  DNA_ALPHABET:
    new AlphabetType(7, "DNA_ALPHABET", [AlphStd.DNA]),
  PROTEIN_ALPHABET:
    new AlphabetType(8, "PROTEIN_ALPHABET", [AlphStd.PROTEIN])
};



/*
 * Allows searching a UTF-8 byte buffer (referred to as a haystack) with a
 * ASCII string (referred to as a needle) using Boyer-Moore string search.
 */
var BMSearchArray = function (needle, ignore_case) {
  "use strict";
  var i, code;
  if (typeof ignore_case !== "boolean") ignore_case = false;
  this.ignore_case = ignore_case;
  this.needle1 = new Uint8Array(needle.length);
  this.needle2 = null;
  if (ignore_case) this.needle2 = new Uint8Array(needle.length);
  for (i = 0; i < needle.length; i++) {
    code = needle.charCodeAt(i);
    if (code >= 0 && code < 128) {
      if (ignore_case) {
        if (code >= 65 && code <= 90) {
          this.needle1[i] = code;
          this.needle2[i] = code + 32;
        } else if (code >= 97 && code <= 122) {
          this.needle1[i] = code - 32;
          this.needle2[i] = code;
        } else {
          this.needle1[i] = code;
          this.needle2[i] = code;
        }
        this.needle1
      } else {
        this.needle1[i] = code;
      }
    } else {
      throw new Error("Search string contains a character outside the ASCII range.");
    }
  }
  this.char_table = this._make_char_table();
  this.offset_table = this._make_offset_table();
};

/*
 * Searches the haystack from the specified index. If the index is
 * negative it assumes that the needle was previously partially matched
 * by that amount prior to the current haystack.
 *
 * Will return an object containing the property "index" and property "complete"
 * where "index" is the match index and "complete" is true if a full match
 * was found (including previous partial match).
 */
BMSearchArray.prototype.indexIn = function (haystack, from_index) {
  "use strict";
  var i, j, h, no_overlap, overhang;
  if (typeof from_index !== "number") from_index = 0;
  if (from_index < -(this.needle1.length - 1) || from_index > haystack.length) 
    throw new Error("The from_index must be a index in the haystack");
  if (this.needle1.length == from_index) return from_index;
  if (this.ignore_case) {
    for (i = from_index + this.needle1.length - 1; i < haystack.length;) {
      j = this.needle1.length - 1;
      h = (i >= 0 ? haystack[i] : this.needle1[i - from_index]);
      while (this.needle1[j] == h || this.needle2[j] == h) {
        if (j == 0) return {"index": i, "complete": true};
        --i; --j;
        h = (i >= 0 ? haystack[i] : this.needle1[i - from_index]);
      }
      if (h >= 0 && h < 128) {
        // within ASCII range
        i += Math.max(this.offset_table[this.needle1.length - 1 - j], this.char_table[h]);
      } else {
        // definitely not in needle as outside ASCII range
        i += this.needle1.length;
      }
    }
    // try to find a partial match
    no_overlap = (haystack.length - 1) + this.needle1.length;
    while (i < no_overlap) {
      j = this.needle1.length - 1;
      // treat the overhanging portion as a match
      overhang = i - (haystack.length - 1);
      j -= overhang;
      i -= overhang;
      // continue testing the overlapping bit as normal
      h = (i >= 0 ? haystack[i] : this.needle1[i - from_index]);
      while (this.needle1[j] == h || this.needle2[j] == h) {
        if (j == 0) return {"index": i, "complete": false};
        --i; --j;
        h = (i >= 0 ? haystack[i] : this.needle1[i - from_index]);
      }
      // now with the mismatch we can only use the charTable
      // because we don't really know what's past the end
      // so using the suffix table doesn't make sense
      if (h >= 0 && h < 128) {
        // within ASCII range
        i += Math.max(this.needle1.length - j, this.char_table[h]);
      } else {
        // definitely not in needle as outside ASCII range
        i += this.needle1.length;
      }
    }
    return null;
  } else {
    for (i = from_index + this.needle1.length - 1; i < haystack.length;) {
      j = this.needle1.length - 1;
      h = (i >= 0 ? haystack[i] : this.needle1[i - from_index]);
      while (this.needle1[j] == h) {
        if (j == 0) return {"index": i, "complete": true};
        --i; --j;
        h = (i >= 0 ? haystack[i] : this.needle1[i - from_index]);
      }
      if (h >= 0 && h < 128) {
        // within ASCII range
        i += Math.max(this.offset_table[this.needle1.length - 1 -j], this.char_table[h]);
      } else {
        // definately not in needle as outside ASCII range
        i += this.needle1.length;
      }
    }
    // try to find a partial match
    no_overlap = (haystack.length - 1) + this.needle1.length;
    while (i < no_overlap) {
      j = this.needle1.length - 1;
      // treat the overhanging portion as a match
      overhang = i - (haystack.length - 1);
      i -= overhang;
      j -= overhang;
      // continue testing the overlapping bit as normal
      h = (i >= 0 ? haystack[i] : this.needle1[i - from_index]);
      while (this.needle1[i] == h) {
        if (j == 0) return {"index": i, "complete": false};
        --i; --j;
        h = (i >= 0 ? haystack[i] : this.needle1[i - from_index]);
      }
      // now with the mismatch we can only use the char_table
      // because we don't really know what's past the end
      // so using the suffix table doesn't make sense
      if (h >= 0 && h < 128) {
        i += Math.max(this.needle1.length - j, this.char_table[h]);
      } else {
        // definately not in needle as outside ASCII range
        i += this.needle1.length;
      }
    }
    return null;
  }
};

/*
 * Returns the length of the needle.
 */
BMSearchArray.prototype.length = function() {
  return this.needle1.length;
};

/*
 * Returns if the search is ignoring case.
 */
BMSearchArray.prototype.is_ignoring_case = function() {
  return this.ignore_case;
};

/*
 * Returns the needle.
 */
BMSearchArray.prototype.toString = function() {
  return String.fromCharCode.apply(null, this.needle1);
};

/*
 * Private function that builds the char table for the Boyer-Moore
 * string search.
 */
BMSearchArray.prototype._make_char_table = function () {
  "use strict";
  var table, i;
  table = new Int32Array(128);
  for (i = 0; i < table.length; i++) {
    table[i] = this.needle1.length;
  }
  for (i = 0; i < this.needle1.length - 1; ++i) {
    if (this.ignore_case) {
      table[this.needle1[i]] = this.needle1.length - 1 - i;
      table[this.needle2[i]] = this.needle1.length - 1 - i;
    } else {
      table[this.needle1[i]] = this.needle1.length - 1 - i;
    }
  }
  return table;
};

/*
 * Private function that builds the offset table for the Boyer-Moore
 * string search.
 */
BMSearchArray.prototype._make_offset_table = function () {
  "use strict";
  var table, i, last_prefix_position, slen;
  table = new Int32Array(this.needle1.length);
  last_prefix_position = this.needle1.length;
  for (i = this.needle1.length - 1; i >= 0; --i) {
    if (this._is_prefix(i + 1)) {
      last_prefix_position = i + 1;
    }
    table[this.needle1.length - 1 - i] = last_prefix_position - i + this.needle1.length - 1;
  }
  for (i = 0; i < this.needle1.length - 1; ++i) {
    slen = this._suffix_length(i);
    table[slen] = this.needle1.length - 1 - i + slen;
  }
  return table;
};

/*
 * Private function that determines if the substring begining at the
 * specified position and extending to the end is also a prefix.
 */
BMSearchArray.prototype._is_prefix = function(pos) {
  "use strict";
  var i, j;
  for (i = pos, j = 0; i < this.needle1.length; ++i, ++j) {
    if (this.needle1[i] != this.needle1[j]) {
      return false;
    }
  }
  return true;
};

/*
 * Private function that determines the longest substring
 * that matches the suffix ending at the specified position.
 */
BMSearchArray.prototype._suffix_length = function(pos) {
  "use strict";
  var len, i, j;
  len = 0;
  for (i = pos, j = this.needle1.length - 1; 
      i >= 0 && this.needle1[i] == this.needle1[j]; --i, --j) {
    len++;
  }
  return len;
};

// Encoding error
var EncodingError = function(message, offset) {
  "use strict";
  this.message = message;
  this.stack = Error().stack;
  this.offset = offset;
};
EncodingError.prototype = Object.create(Error.prototype);
EncodingError.prototype.name = "EncodingError";
EncodingError.prototype.code = EncErrorType.GENERIC;
EncodingError.prototype.constructor = EncodingError;

// Encoding error - Illegal 0xFE or 0xFF byte
var IllegalByteError = function(byte_value, byte_offset) {
  EncodingError.call(this, "Invalid UTF-8 start byte 0x" + 
      byte_value.toString(16) + " at offset "+ byte_offset + ".", byte_offset);
  this.byte_value = byte_value;
};
IllegalByteError.prototype = Object.create(EncodingError.prototype);
IllegalByteError.prototype.name = "IllegalByteError";
IllegalByteError.prototype.code = EncErrorType.ILLEGAL_BYTE;
IllegalByteError.constructor = IllegalByteError;
// Encoding error - Incomplete codepoint
var IncompleteCodeunitError = function (expected_length, found_length, codeunit_offset) {
  EncodingError.call(this, "Expected " + expected_length + " bytes in " +
      "the code-unit at offset " + codeunit_offset + " but only found " +
      found_length + " bytes.", codeunit_offset);
  this.expected_length = expected_length;
  this.found_length = found_length;
};
IncompleteCodeunitError.prototype = Object.create(EncodingError.prototype);
IncompleteCodeunitError.prototype.name = "IncompleteCodeunitError";
IncompleteCodeunitError.prototype.code = EncErrorType.INCOMPLETE_CODE_UNIT;
IncompleteCodeunitError.constructor = IncompleteCodeunitError;
// Encoding error - Overlong codepoint
var OverlongEncodingError = function (codepoint, encoded_length, codeunit_offset) {
  EncodingError.call(this, "Encoding of U+" + codepoint.toString(16) +
      " at offset " + codeunit_offset + " used " + encoded_length +
      " bytes which is more than allowed.", codeunit_offset);
  this.codepoint = codepoint;
  this.encoded_length = encoded_length;
};
OverlongEncodingError.prototype = Object.create(EncodingError.prototype);
OverlongEncodingError.prototype.name = "OverlongEncodingError";
OverlongEncodingError.prototype.code = EncErrorType.OVERLONG_ENCODING;
OverlongEncodingError.constructor = OverlongEncodingError;
// Encoding error - Invalid codepoints
var InvalidCodepointError = function (codepoint, codeunit_offset) {
  EncodingError.call(this, "The value U+" + codepoint.toString(16) +
      " at offset " + codeunit_offset + " is not a legal Unicode value.",
      codeunit_offset);
  this.codepoint = codepoint;
};
InvalidCodepointError.prototype = Object.create(EncodingError.prototype);
InvalidCodepointError.prototype.name = "InvalidCodepointError";
InvalidCodepointError.prototype.code = EncErrorType.INVALID_CODE_POINT;
InvalidCodepointError.constructor = InvalidCodepointError;
// Encoding error - Invalid codepoints
var NulByteError = function (byte_offset) {
  EncodingError.call(this, "NUL byte at offset " + byte_offset +
      " is not allowed in plain text.", byte_offset);
};
NulByteError.prototype = Object.create(EncodingError.prototype);
NulByteError.prototype.name = "NulByteError";
NulByteError.prototype.code = EncErrorType.NULL_BYTE;
NulByteError.constructor = NulByteError;

/*
 * Reads UTF-8 encoded bytes into codepoints.
 */
var Utf8Decoder = function () {
  this.source = null;
  this.last_source = false;
  this.pos = 0;
  this.code_len = 0;
  this.code_read = 0;
  this.code_data = 0;
  this.code_min = 0;
  this.source_offset = 0; // offset of the source array
};

/*
 * Reset the reader to the newly constructed state.
 */
Utf8Decoder.prototype.reset = function () {
  this.source = null;
  this.last_source = false;
  this.pos = 0;
  this.code_len = 0;
  this.code_read = 0;
  this.code_data = 0;
  this.code_min = 0;
  this.source_offset = 0;
};

/*
 * Return the position in the current byte array.
 */
Utf8Decoder.prototype.position = function () {
  return this.pos;
};

/*
 * Return the offset in all the sources.
 */
Utf8Decoder.prototype.blob_position = function () {
  return this.source_offset + this.pos;
};

/*
 * Set the byte array source. Throws an error if a non-exhausted source is
 * already present.
 */
Utf8Decoder.prototype.set_source = function (ab, last) {
  if (this.source != null && this.pos < this.source.length)
    throw new Error("Not finished with the current source yet!");
  if (this.last_source)
    throw new Error("Previous source was meant to be the last source!");
  this.source_offset += (this.source != null ? this.source.length : 0);
  this.source = ab;
  this.pos = 0;
  this.last_source = !!last;
};

/*
 * Get the next codepoint. Returns null if there are insufficient bytes. 
 * Throws an error when the encoding is incorrect.
 */
Utf8Decoder.prototype.next = function() {
  "use strict";
  var bite, i, temp;
  if (this.source == null || this.pos >= this.source.length) {
    this.code_len = 0;
    this.code_read = 0;
    return null;
  }
  if (this.code_len == this.code_read) {
    bite = this.source[this.pos++];
    if (bite <= 127) return bite;
    // outside of ASCII range
    if ((bite & 0xE0) == 0xC0) { // 110xxxxx = 2 bytes
      this.code_len = 2;
      this.code_data = bite & 0x1F;
      this.code_min = 0x7F; // must need more than 7 bits
    } else if ((bite & 0xF0) == 0xE0) { // 1110xxxx = 3 bytes
      this.code_len = 3;
      this.code_data = bite & 0x0F;
      this.code_min = 0x7FF; // must need more than 11 bits
    } else if ((bite & 0xF8) == 0xF0) { // 11110xxx = 4 bytes
      this.code_len = 4;
      this.code_data = bite & 0x07;
      this.code_min = 0xFFFF; // must need more than 16 bits
    } else if ((bite & 0xFC) == 0xF8) { // 111110xx = 5 bytes
      this.code_len = 5;
      this.code_data = bite & 0x03;
      this.code_min = 0x1FFFFF; // must need more than 21 bits
    } else if ((bite & 0xFE) == 0xFC) { // 1111110x = 6 bytes
      this.code_len = 6;
      this.code_data = bite & 0x01;
      this.code_min = 0x3FFFFFF; // must need more than 26 bits
    } else { // may be 10xxxxxx or the invalid bytes 0xFE or 0xFF
      this.code_len = 0;
      this.code_read = 1;
      throw new IllegalByteError(bite, this.source_offset + this.pos - 1);
    }
    this.code_read = 1;
  }
  // read each continuation byte
  for (; this.pos < this.source.length && this.code_read < this.code_len; this.pos++, this.code_read++) {
    bite = this.source[this.pos];
    if ((bite & 0xC0) == 0x80) {
      this.code_data = ((this.code_data << 6) | (bite & 0x3F));
    } else { // not 10xxxxxx
      temp = this.code_len;
      this.code_len = 0;
      throw new IncompleteCodeunitError(temp, this.code_read, 
          this.source_offset + this.pos - this.code_read);
    }
  }
  if (this.code_read == this.code_len) {
    // check that the encoding is in minimal form
    if (this.code_data <= this.code_min) {
      throw new OverlongEncodingError(this.code_data, this.code_len,
          this.source_offset + this.pos - this.code_len);
    } else if (this.code_data > 0x10FFFF) {
      this.code_len = 0;
      throw new InvalidCodepointError(this.code_data, this.source_offset + this.pos - this.code_read);
    }
    return this.code_data;
  } else if (this.last_source) {
    this.code_len = 0;
    throw new IncompleteCodeunitError(this.code_len, this.code_read,
        this.source_offset + this.pos - this.code_read);
  } else {
    return null;
  }
};

/*
 * Returns the byte offset of the codepoint returned by next.
 * When next returns null this returns the byte offset of the start
 * of the incomplete codepoint.
 */
Utf8Decoder.prototype.offset = function() {
  return this.source_offset + this.pos - this.code_read;
};

/*
 * Decode a UTF-8 file and keep track of the line number and column count.
 */
var LineDecoder = function() {
  Utf8Decoder.call(this);
  this.line_num = 0;
  this.col_num = -1;
  this.last_codepoint = 0;
  this.last_is_nl = false;
  this.last_is_nl_start = false;
};
LineDecoder.prototype = Object.create(Utf8Decoder.prototype);
LineDecoder.prototype.constructor = LineDecoder;

/*
 * Get the next character in a plain text file
 */
LineDecoder.prototype.next = function() {
  "use strict";
  var codepoint, is_nl, is_nl_start;
  // call superclass implementation to get the codepoint
  codepoint = Utf8Decoder.prototype.next.call(this);
  if (codepoint == null) return null;
  if (codepoint == 0) {
    throw new NulByteError(this.offset());
  }
  this.col_num++;
  is_nl = (codepoint == 10 || codepoint == 13);
  // see if the previous was a newline
  if (this.last_is_nl) {
    // see what the current char is
    if (!is_nl) {
      // not a newline, so line is incremented by previous newline
      this.line_num++;
      this.col_num = 0;
    } else {
      // current char is a newline, so line is only incremented if 
      // the previous is identical char or it is not a newline beginning
      if (this.last_codepoint == codepoint || !this.last_is_nl_start) {
        this.line_num++;
        this.col_num = 0;
        // so by deduction the current char is a newline beginning
        is_nl_start = true;
      }
    }
  } else {
    // previous char is not a newline
    if (is_nl) {
      // so by deduction the current char is a newline beginning
      is_nl_start = true;
    }
  }
  this.last_codepoint = codepoint;
  this.last_is_nl = is_nl;
  this.last_is_nl_start = is_nl_start;
  return codepoint;
};

/*
 * Get the line number of the previously returned character
 */
LineDecoder.prototype.line = function() {
  return this.line_num;
};

/*
 * Get the column number of the previously returned character
 */
LineDecoder.prototype.column = function() {
  return this.col_num;
};


/*
 * Does a byte comparison.
 */
function __is_bytes(chunk, byte_list) {
  "use strict";
  var i;
  if (chunk.length < byte_list.length) return false;
  for (i = 0; i < byte_list.length; i++) {
    if (chunk[i] != byte_list[i]) return false;
  }
  return true;
}


/*
 * Tries to return the file format if a file is NOT plain text UTF-8.
 *
 * Looks at the initial chunk of a file and tries to identify the file
 * encoding. Also looks for common file formats like rich text
 * format.
 */
function unusable_format(chunk, max_scan, file_name) {
  "use strict";
  var i, nulls, extension, sub;
  if (typeof max_scan != "number") max_scan = chunk.length;
  max_scan = Math.min(Math.max(4, max_scan), chunk.length);
  if (file_name != null) {
    extension = file_name.substring(file_name.lastIndexOf(".") + 1).toLowerCase();
  } else {
    extension = "";
  }
  if (chunk.length >= 2) {
    if (chunk[0] == 0xFE && chunk[1] == 0xFF) { 
      //UTF-16 big endian BOM
      return {name: "UTF-16", type: FileType.ENCODING};

    } else if (chunk[0] == 0xFF && chunk[1] == 0xFE && 
        (chunk.length < 4 || chunk[2] != 0x0 || chunk[3] != 0x0)) {
      //UTF-16 little endian BOM
      return {name: "UTF-16", type: FileType.ENCODING};
    }
  }
  if (__is_bytes(chunk, [0xEF, 0xBB, 0xBF])) {
    // evil!! UTF-8 BOM
    return {name: "UTF-8 with a byte order mark", type: FileType.ENCODING};
  }
  if (__is_bytes(chunk, [0xFF, 0xFE, 0x0, 0x0])) {
    // UTF-32 little endian BOM
    return {name: "UTF-32", type: FileType.ENCODING};
  }
  if (__is_bytes(chunk, [0x0, 0x0, 0xFE, 0xFF])) {
    // UTF-32 big endian BOM
    return {name: "UTF-32", type: FileType.ENCODING};
  }
  // test for rich text format
  if (__is_bytes(chunk, [0x7B, 0x5C, 0x72, 0x74, 0x66, 0x31])) {
    return {name: "rich text", type: FileType.BINARY};
  }
  // test for old MS Office, hoping to catch doc files here
  if (__is_bytes(chunk, [0xD0, 0xCF, 0x11, 0xE0, 0xA1, 0xB1, 0x1A, 0xE1])) {
    if (chunk.length > 520) { // check we can get 8 bytes at offset 512
      sub = chunk.subarray(512, 516);
      // check word doc sub header
      if (__is_bytes(sub, [0xEC, 0xA5, 0xC1, 0x00])) {
        return {name: "Microsoft Word document", type: FileType.BINARY};
      }
      // check powerpoint sub header
      if (extension.toUpperCase() == "PPT" ||
          __is_bytes(sub, [0x00, 0x6E, 0x1E, 0xF0]) ||
          __is_bytes(sub, [0x0F, 0x00, 0xE8, 0x03]) ||
          __is_bytes(sub, [0xA0, 0x46, 0x1D, 0xF0]) || 
          __is_bytes(sub, [0xFD, 0xFF, 0xFF, 0xFF, 0x0E, 0x00, 0x00, 0x00]) ||
          __is_bytes(sub, [0xFD, 0xFF, 0xFF, 0xFF, 0x1C, 0x00, 0x00, 0x00]) ||
          __is_bytes(sub, [0xFD, 0xFF, 0xFF, 0xFF, 0x43, 0x00, 0x00, 0x00])) {
        return {name: "Microsoft PowerPoint presentation", type: FileType.BINARY};
      }
      // check excel sub header
      if (extension.toUpperCase() == "XLS" ||
          __is_bytes(sub, [0x09, 0x08, 0x10, 0x00, 0x00, 0x06, 0x05, 0x00]) ||
          __is_bytes(sub, [0xFD, 0xFF, 0xFF, 0xFF, 0x10]) ||
          __is_bytes(sub, [0xFD, 0xFF, 0xFF, 0xFF, 0x1F]) ||
          __is_bytes(sub, [0xFD, 0xFF, 0xFF, 0xFF, 0x22]) ||
          __is_bytes(sub, [0xFD, 0xFF, 0xFF, 0xFF, 0x23]) ||
          __is_bytes(sub, [0xFD, 0xFF, 0xFF, 0xFF, 0x28]) ||
          __is_bytes(sub, [0xFD, 0xFF, 0xFF, 0xFF, 0x29])) {
        return {name: "Microsoft Excel spreadsheet", type: FileType.BINARY}
      }
    }
    return {name: "Microsoft Office document", type: FileType.BINARY};
  }
  // test for DOCX, PPTX, XLSX (any zip, but we test extension)
  if (__is_bytes(chunk, [0x50, 0x4B, 0x03, 0x04])) {
    if (extension.toUpperCase() == "DOCX") {
      return {name: "Microsoft Word document", type: FileType.BINARY};
    } else if (extension.toUpperCase() == "PPTX") {
      return {name: "Microsoft PowerPoint presentation", type: FileType.BINARY};
    } else if (extension.toUpperCase() == "XSLX") {
      return {name: "Microsoft Excel spreadsheet", type: FileType.BINARY};
    } else if (extension.toUpperCase() == "ZIP") {
      return {name: "ZIP archive", type: FileType.COMPRESSED};
    }
  }
  // test for GZIP
  if (__is_bytes(chunk, [0x1F, 0x8B, 0x08])) {
    return {name: "GZIP archive", type: FileType.COMPRESSED};
  }
  // test for PDF
  if (__is_bytes(chunk, [0x25, 0x50, 0x44, 0x46]) && extension.toUpperCase() == "PDF") {
    return {name: "Adobe PDF", type: FileType.BINARY};
  }

  // no magic number found, this is good news
  // assuming that this should be mostly ASCII (without U+0000)
  // do a scan for null bytes to see if there's a pattern
  // that indicates UTF-16 or UTF-32
  nulls = [0, 0, 0, 0];
  for (i = 0; i < max_scan; i++) {
    nulls[i % 4] += (chunk[i] === 0x0) ? 1 : 0;
  }
  if (nulls[0] + nulls[1] + nulls[2] + nulls[3] > 0) {
    if (nulls[0] != 0 && nulls[1] != 0 && nulls[2] != 0 && nulls[3] == 0) {
      // UTF-32 big endian
      return {name: "UTF-32", type: FileType.ENCODING};
    } else if (nulls[0] == 0 && nulls[1] != 0 && nulls[2] != 0 && nulls[3] != 0) {
      // UTF-32 little endian
      return {name: "UTF-32", type: FileType.ENCODING};
    } else if ((nulls[0] + nulls[2]) != 0 && (nulls[1] + nulls[3]) == 0) {
      // UTF-16 big endian
      return {name: "UTF-16", type: FileType.ENCODING};
    } else if ((nulls[0] + nulls[2]) == 0 && (nulls[1] + nulls[3]) != 0) {
      // UTF-16 little endian
      return {name: "UTF-16", type: FileType.ENCODING};
    } else {
      // binary format
      return {name: "unknown binary", type: FileType.BINARY};
    }
  }
  return null; // looks ok
};

//******************************************************************************
// File Faults
// keep track of the lines that faults occurred on.
//******************************************************************************
var FileFaults = function (line_max, letter_max) {
  this.count = 0; // count of faults
  // maximum number of line ranges to store
  this.line_max = (typeof line_max == "number" && line_max >= 1 ? line_max : 10);
  this.line_count = 0; // the total count of lines with faults
  this.lines = [];
  // maximum number of letters to store
  this.letter_max = (typeof letter_max == "number" && letter_max >= 0 ? letter_max : 5);
  this.letter_stored = 0; // tracks count stored (not the count seen)
  this.letters = {};
};

// {{{ FileFaults.control_char_re (very long regular expression)
FileFaults.control_char_re = /[\0-\x1F\x7F-\x9F\xAD\u0378\u0379\u037F-\u0383\u038B\u038D\u03A2\u0528-\u0530\u0557\u0558\u0560\u0588\u058B-\u058E\u0590\u05C8-\u05CF\u05EB-\u05EF\u05F5-\u0605\u061C\u061D\u06DD\u070E\u070F\u074B\u074C\u07B2-\u07BF\u07FB-\u07FF\u082E\u082F\u083F\u085C\u085D\u085F-\u089F\u08A1\u08AD-\u08E3\u08FF\u0978\u0980\u0984\u098D\u098E\u0991\u0992\u09A9\u09B1\u09B3-\u09B5\u09BA\u09BB\u09C5\u09C6\u09C9\u09CA\u09CF-\u09D6\u09D8-\u09DB\u09DE\u09E4\u09E5\u09FC-\u0A00\u0A04\u0A0B-\u0A0E\u0A11\u0A12\u0A29\u0A31\u0A34\u0A37\u0A3A\u0A3B\u0A3D\u0A43-\u0A46\u0A49\u0A4A\u0A4E-\u0A50\u0A52-\u0A58\u0A5D\u0A5F-\u0A65\u0A76-\u0A80\u0A84\u0A8E\u0A92\u0AA9\u0AB1\u0AB4\u0ABA\u0ABB\u0AC6\u0ACA\u0ACE\u0ACF\u0AD1-\u0ADF\u0AE4\u0AE5\u0AF2-\u0B00\u0B04\u0B0D\u0B0E\u0B11\u0B12\u0B29\u0B31\u0B34\u0B3A\u0B3B\u0B45\u0B46\u0B49\u0B4A\u0B4E-\u0B55\u0B58-\u0B5B\u0B5E\u0B64\u0B65\u0B78-\u0B81\u0B84\u0B8B-\u0B8D\u0B91\u0B96-\u0B98\u0B9B\u0B9D\u0BA0-\u0BA2\u0BA5-\u0BA7\u0BAB-\u0BAD\u0BBA-\u0BBD\u0BC3-\u0BC5\u0BC9\u0BCE\u0BCF\u0BD1-\u0BD6\u0BD8-\u0BE5\u0BFB-\u0C00\u0C04\u0C0D\u0C11\u0C29\u0C34\u0C3A-\u0C3C\u0C45\u0C49\u0C4E-\u0C54\u0C57\u0C5A-\u0C5F\u0C64\u0C65\u0C70-\u0C77\u0C80\u0C81\u0C84\u0C8D\u0C91\u0CA9\u0CB4\u0CBA\u0CBB\u0CC5\u0CC9\u0CCE-\u0CD4\u0CD7-\u0CDD\u0CDF\u0CE4\u0CE5\u0CF0\u0CF3-\u0D01\u0D04\u0D0D\u0D11\u0D3B\u0D3C\u0D45\u0D49\u0D4F-\u0D56\u0D58-\u0D5F\u0D64\u0D65\u0D76-\u0D78\u0D80\u0D81\u0D84\u0D97-\u0D99\u0DB2\u0DBC\u0DBE\u0DBF\u0DC7-\u0DC9\u0DCB-\u0DCE\u0DD5\u0DD7\u0DE0-\u0DF1\u0DF5-\u0E00\u0E3B-\u0E3E\u0E5C-\u0E80\u0E83\u0E85\u0E86\u0E89\u0E8B\u0E8C\u0E8E-\u0E93\u0E98\u0EA0\u0EA4\u0EA6\u0EA8\u0EA9\u0EAC\u0EBA\u0EBE\u0EBF\u0EC5\u0EC7\u0ECE\u0ECF\u0EDA\u0EDB\u0EE0-\u0EFF\u0F48\u0F6D-\u0F70\u0F98\u0FBD\u0FCD\u0FDB-\u0FFF\u10C6\u10C8-\u10CC\u10CE\u10CF\u1249\u124E\u124F\u1257\u1259\u125E\u125F\u1289\u128E\u128F\u12B1\u12B6\u12B7\u12BF\u12C1\u12C6\u12C7\u12D7\u1311\u1316\u1317\u135B\u135C\u137D-\u137F\u139A-\u139F\u13F5-\u13FF\u169D-\u169F\u16F1-\u16FF\u170D\u1715-\u171F\u1737-\u173F\u1754-\u175F\u176D\u1771\u1774-\u177F\u17DE\u17DF\u17EA-\u17EF\u17FA-\u17FF\u180F\u181A-\u181F\u1878-\u187F\u18AB-\u18AF\u18F6-\u18FF\u191D-\u191F\u192C-\u192F\u193C-\u193F\u1941-\u1943\u196E\u196F\u1975-\u197F\u19AC-\u19AF\u19CA-\u19CF\u19DB-\u19DD\u1A1C\u1A1D\u1A5F\u1A7D\u1A7E\u1A8A-\u1A8F\u1A9A-\u1A9F\u1AAE-\u1AFF\u1B4C-\u1B4F\u1B7D-\u1B7F\u1BF4-\u1BFB\u1C38-\u1C3A\u1C4A-\u1C4C\u1C80-\u1CBF\u1CC8-\u1CCF\u1CF7-\u1CFF\u1DE7-\u1DFB\u1F16\u1F17\u1F1E\u1F1F\u1F46\u1F47\u1F4E\u1F4F\u1F58\u1F5A\u1F5C\u1F5E\u1F7E\u1F7F\u1FB5\u1FC5\u1FD4\u1FD5\u1FDC\u1FF0\u1FF1\u1FF5\u1FFF\u200B-\u200F\u202A-\u202E\u2060-\u206F\u2072\u2073\u208F\u209D-\u209F\u20BB-\u20CF\u20F1-\u20FF\u218A-\u218F\u23F4-\u23FF\u2427-\u243F\u244B-\u245F\u2700\u2B4D-\u2B4F\u2B5A-\u2BFF\u2C2F\u2C5F\u2CF4-\u2CF8\u2D26\u2D28-\u2D2C\u2D2E\u2D2F\u2D68-\u2D6E\u2D71-\u2D7E\u2D97-\u2D9F\u2DA7\u2DAF\u2DB7\u2DBF\u2DC7\u2DCF\u2DD7\u2DDF\u2E3C-\u2E7F\u2E9A\u2EF4-\u2EFF\u2FD6-\u2FEF\u2FFC-\u2FFF\u3040\u3097\u3098\u3100-\u3104\u312E-\u3130\u318F\u31BB-\u31BF\u31E4-\u31EF\u321F\u32FF\u4DB6-\u4DBF\u9FCD-\u9FFF\uA48D-\uA48F\uA4C7-\uA4CF\uA62C-\uA63F\uA698-\uA69E\uA6F8-\uA6FF\uA78F\uA794-\uA79F\uA7AB-\uA7F7\uA82C-\uA82F\uA83A-\uA83F\uA878-\uA87F\uA8C5-\uA8CD\uA8DA-\uA8DF\uA8FC-\uA8FF\uA954-\uA95E\uA97D-\uA97F\uA9CE\uA9DA-\uA9DD\uA9E0-\uA9FF\uAA37-\uAA3F\uAA4E\uAA4F\uAA5A\uAA5B\uAA7C-\uAA7F\uAAC3-\uAADA\uAAF7-\uAB00\uAB07\uAB08\uAB0F\uAB10\uAB17-\uAB1F\uAB27\uAB2F-\uABBF\uABEE\uABEF\uABFA-\uABFF\uD7A4-\uD7AF\uD7C7-\uD7CA\uD7FC-\uF8FF\uFA6E\uFA6F\uFADA-\uFAFF\uFB07-\uFB12\uFB18-\uFB1C\uFB37\uFB3D\uFB3F\uFB42\uFB45\uFBC2-\uFBD2\uFD40-\uFD4F\uFD90\uFD91\uFDC8-\uFDEF\uFDFE\uFDFF\uFE1A-\uFE1F\uFE27-\uFE2F\uFE53\uFE67\uFE6C-\uFE6F\uFE75\uFEFD-\uFF00\uFFBF-\uFFC1\uFFC8\uFFC9\uFFD0\uFFD1\uFFD8\uFFD9\uFFDD-\uFFDF\uFFE7\uFFEF-\uFFFB\uFFFE\uFFFF]/g;
// }}}

FileFaults.prototype.add = function (line, letter) {
  "use strict";
  var last_line;
  this.count++;
  if (letter != null && this.letter_stored < this.letter_max && 
      !FileFaults.control_char_re.test(letter) && !this.letters[letter]) {
    this.letter_stored++;
    this.letters[letter] = true;
  }
  if (this.lines.length === 0) {
    this.line_count++;
    this.lines.push(line);
  } else {
    last_line = this.lines[this.lines.length - 1];
    if (typeof last_line == "number") {
      if (line == (last_line + 1)) {
        this.line_count++;
        this.lines.pop();
        this.lines.push({begin: last_line, end: line});
      } else if (line > last_line) {
        this.line_count++;
        if (this.lines.length < this.line_max) this.lines.push(line);
      }
    } else {
      if (line == (last_line.end + 1)) {
        this.line_count++;
        last_line.end = line;
      } else if (line > last_line.end) {
        this.line_count++;
        if (this.lines.length < this.line_max) this.lines.push(line);
      }
    }
  }
};

FileFaults.prototype.faults = function() {
  return this.count;
};

FileFaults.prototype.letters_sep = function() {
  "use strict";
  var letter, i, text;
  i = 0;
  text = "";
  for (letter in this.letters) {
    if (!this.letters.hasOwnProperty(letter)) continue;
    if (i == 0) { // first
    } else if (i == (this.letter_stored - 1)) { // last
      text += " and ";
    } else { // middle
      text += ", ";
    }
    text += "'" + letter + "'";
    i++;
  }
  return text;
};

FileFaults.prototype.letters_join = function() {
  "use strict";
  var letter, text;
  text = "";
  for (letter in this.letters) {
    if (!this.letters.hasOwnProperty(letter)) continue;
    text += letter;
  }
  if (text.length == 0) {
    return "";
  } else if (text.length == 1) {
    return "'" + text + "'";
  } else {
    return "\"" + text + "\"";
  }
};

FileFaults.prototype.lines_str = function() {
  "use strict";
  var text, i, range, count;
  text = "";
  count = 0;
  for (i = 0; i < this.lines.length; i++) {
    range = this.lines[i];
    count += (typeof range == "number" ? 1 : range.end - range.begin + 1);
    if (i == 0) { // no prefix on first
    } else if (i == (this.lines.length - 1) && count == this.line_count) {
      text += " and ";
    } else {
      text += ", ";
    }
    if (typeof range == "number") {
      text += (range + 1);
    } else {
      text += (range.begin + 1) + "-" + (range.end + 1);
    }
  }
  if (count < this.line_count) {
    text += ", ...";
  }
  return text;
};

/*
 * nop
 * Do nothing
 */
function nop() { }

/*
 * stop_evt
 *
 * Stop propagation of a event.
 */
function stop_evt(evt) {
  evt = evt || window.event;
  if (typeof evt.stopPropagation != "undefined") {
    evt.stopPropagation();
  } else {
    evt.cancelBubble = true;
  }
}



/*
 * has_class
 *
 * Checks for the existance of a class on a node.
 */
function has_class(node, cls) {
  var classes = node.className;
  var list = classes.replace(/^\s+/, '').replace(/\s+$/, '').split(/\s+/);
  for (var i = 0; i < list.length; i++) {
    if (list[i] == cls) {
      return true;
    }
  }
  return false;
}

/*
 * more_opts
 *
 * Shows more options by toggling an expansion section.
 * When hiding the expanded options it checks if they've
 * changed by calling a caller defined function and gives
 * the user a chance to abort. If the user continues the
 * options are hidden and reset to their defaults by calling
 * another caller defined function.
 */
function more_opts(node, fn_changed) {
  if (has_class(node, 'expanded')) {
    toggle_class(node, 'modified', fn_changed());
  }
  toggle_class(node, 'expanded');
}

/*
 * more_opts_key
 *
 * Calls more_opts when spacebar has been pressed.
 */
function more_opts_key(evt, node, fn_changed) {
  if (evt.which == 32 || evt.keyCode == 32) {
    more_opts(node, fn_changed);
  }
}

/*
 * more_opts_reset
 *
 * Prompts the user and resets the options if requested.
 */
function more_opts_reset(evt, node, title, fn_changed, fn_reset) {
  stop_evt(evt);
  if (!fn_changed()) return;
  if (!confirm("Reset the " + title + " to the defaults?")) return;
  fn_reset();
  toggle_class(node, 'modified', false);
}

/*
 * int_key
 *
 * Limits input to numbers and various control characters such as 
 * enter, backspace and delete.
 */
function int_key(evt) {
  evt = evt || window.event;
  var code = evt.keyCode;
  var keychar = String.fromCharCode(evt.which || evt.keyCode);
  var numre = /[0-9]/;
  // only allow 0-9 and various control characters
  if (code != 8 && // backspace
      code != 9 && // tab
      code != 13 && // enter
      code != 37 && // left arrow
      code != 38 && // up arrow
      code != 39 && // right arrow
      code != 40 && // down arrow
      code != 46 && // delete
      !numre.test(keychar)) {
    evt.preventDefault();
  }
}

/*
 * num_key
 *
 * Limits input to numbers, decimal, e and various control characters such as 
 * enter, backspace and delete.
 */
function num_key(evt) {
  evt = evt || window.event;
  var code = evt.keyCode;
  var keychar = String.fromCharCode(evt.which || evt.keyCode);
  var numre = /[0-9\-\.eE]/;
  // only allow 0-9, - and various control characters
  if (code != 8 && // backspace
      code != 9 && // tab
      code != 13 && // enter
      code != 37 && // left arrow
      code != 38 && // up arrow
      code != 39 && // right arrow
      code != 40 && // down arrow
      code != 46 && // delete
      !numre.test(keychar)) {
    evt.preventDefault();
  }
}

/*
 * Sums up the sizes of all the form fields in bytes.
 */
function calculate_form_size(form) {
  var size = 0;
  var i, j, field, blob;
  try {
    for (i = 0; i < form.length; i++) {
      field = form[i];
      if (field.disabled) continue;
      if (field.name == null || field.name == "") continue;
      if (field.tagName == "INPUT") {
        if (field.type.toUpperCase() == "FILE") {
          for (j = 0; j < field.files.length; j++) {
            size += field.files[j].size;
          }
        } else {
          blob = new Blob([field.value], {type: 'text/plain'});
          size += blob.size;
        }
      } else if (field.tagName == "TEXTAREA") {
        blob = new Blob([field.value], {type: 'text/plain'});
        size += blob.size;
      }
    }
  } catch (e) {
    if (window.console && console.log) {
      console.log("Failed to calculate form size: " + e);
    }
  }
  return size;
}


/*
 * Checks the value of the email field as defined in component_job_details.tmpl
 * I considered putting this in a separate script like the others 
 * (ie component_job_details.js) but it's so trivial it hardly seemed worth it.
 * Also checks the combined size of the form fields to hopefully avoid people
 * getting their large files reject with no explanation.
 */
function check_job_details(skip_email_check) {
  "use strict";
  var email, form_size;
  if (!skip_email_check) {
    email = $("email").value;
    // allow the email field to be empty
    if (!/^\s*$/.test(email)) {
      // if the field contains something then check that an @ is present
      if (email.indexOf("@") === -1) {
        alert("Please correct or remove the email used for recieving the job details.\n");
        return false;
      }
    }
  }
  form_size = calculate_form_size($("email").form) + 1024;
  if (form_size > (1024 * 1024 * 80)) {
    alert("The combined size of the form inputs exceeds 80MB and is too large to be accepted by the server.\n");
    return false;
  }
  return true;
};

/*
 * Checks the value of the integer field.
 */
function check_int_value(description, id, min, max, default_value) {
  "use strict";
  var elem, value, num;
  elem = $(id);
  if (elem.disabled) return true; // skip disabled elements
  value = elem.value;
  num = +(value);
  if (!/^\s*\d+\s*$/.test(value) || typeof num !== "number" || isNaN(num) || 
      (num % 1 != 0) ||
      (min != null && num < min) || 
      (max != null && num > max)) {
      alert("Please input a whole number" +
          (min != null ? " \u2265 " + min : "") + 
          (min != null && max != null ? " and" : "") + 
          (max != null ? " \u2264 " + max : "") +
          " for the " + description + "." + 
          (default_value != null ? " A typical value would be " +
           default_value + ".\n" : "")
      );

    return false;
  }
  return true;
}

function check_int_range(description1, id1, default1, description2, id2, default2, min, max) {
  "use strict";
  var elem1, elem2, num1, num2;
  if (!check_int_value(description1, id1, min, max, default1)) return false;
  if (!check_int_value(description2, id2, min, max, default2)) return false;
  elem1 = $(id1);
  elem2 = $(id2);
  if (elem1.disabled || elem2.disabled) return true;
  num1 = +(elem1.value);
  num2 = +(elem2.value);
  if (num1 > num2) {
    alert("Please choose values for the " + description1 + " and the " +
        description2 + " so the " + description1 + " is \u2264 the " +
        description2 + ".");
    return false;
  }
  return true;
}

/*
 * Checks the value of the numeric field.
 */
function check_num_value(description, id, min, max, default_value) {
  "use strict";
  var elem, num;
  elem = $(id);
  if (elem.disabled) return true; // skip disabled elements
  num = +(elem.value);
  if (typeof num !== "number" || isNaN(num) || 
      (min != null && num <= min) || 
      (max != null && num > max)) {
      alert("Please input a number" +
          (min != null ? " > " + min : "") + 
          (min != null && max != null ? " and" : "") + 
          (max != null ? " \u2264 " + max : "") +
          " for the " + description + "." + 
          (default_value != null ? " A typical value would be " +
           default_value + ".\n" : "")
      );
    return false;
  }
  return true;
}

function check_num_range(description1, id1, default1, description2, id2, default2, min, max) {
  "use strict";
  var elem1, elem2, num1, num2;
  if (!check_num_value(description1, id1, min, max, default1)) return false;
  if (!check_num_value(description2, id2, min, max, default2)) return false;
  elem1 = $(id1);
  elem2 = $(id2);
  if (elem1.disabled || elem2.disabled) return true;
  num1 = +(elem1.value);
  num2 = +(elem2.value);
  if (num1 > num2) {
    alert("Please choose values for the " + description1 + " and the " +
        description2 + " so the " + description1 + " is \u2264 the " +
        description2 + ".");
    return false;
  }
  return true;
}


/*
 * validate_background
 *
 * Takes a markov model background in text form, parses and validates it.
 * Returns an object with the boolean property ok. On success it is true
 * on error it is false. On error it has the additional properties: line and
 * reason. On success it has the property alphabet.
 */
function validate_background(bg_str, delta) {
  var lines, items, i, fields, prob, j, alphabet;
  var indexes, letters, prob_sum;
  if (typeof delta !== "number") delta = 0.01;
  lines = bg_str.split(/\r\n|\r|\n/g);
  items = [];
  for (i = 0; i < lines.length; i++) {
    //skip blank lines and comments
    if (/^\s*$/.test(lines[i]) || /^#/.test(lines[i])) {
      continue;
    }
    // split into two space separated items, allowing for leading and trailing space
    if(!(fields = /^\s*(\S+)\s+(\S+)\s*/.exec(lines[i]))) {
      return {ok: false, line: i+1, reason: "expected 2 space separated fields"};
    }
    prob = +(fields[2]);
    if (typeof prob !== "number" || isNaN(prob) || prob < 0 || prob > 1) {
      return {ok: false, line: i+1, 
        reason: "expected number between 0 and 1 inclusive for the probability field"};
    }
    items.push({key: fields[1], probability: prob, line: i+1});
  }
  // now sort entries
  items.sort( 
    function(a, b) {
      if (a.key.length != b.key.length) { 
        return a.key.length - b.key.length;
      } 
      return a.key.localeCompare(b.key); 
    } 
  );
  // determine the alphabet (all 1 character keys)
  alphabet = [];
  prob_sum = 0;
  for (i = 0; i < items.length; i++) {
    if (items[i].key.length > 1) {
      break;
    }
    alphabet.push(items[i].key);
    prob_sum += items[i].probability;
    if (i > 0 && alphabet[i] == alphabet[i-1]) {
      return {ok: false, line: Math.max(items[i-1].line, items[i].line), 
        reason: "alphabet letter seen before"}; // alphabet letter repeats!
    }
  }
  // test that the probabilities approximately sum to 1
  if (Math.abs(prob_sum - 1.0) > delta) {
    return {ok: false, line: items[alphabet.length -1].line, 
      reason: "order-0 probabilities do not sum to 1"};
  }
  //now iterate over all items checking the key is as expected
  indexes = [0, 0];
  order_complete = true;
  prob_sum = 0;
  for (j = alphabet.length; j < items.length; j++) {
    letters = items[j].key.split("");
    prob_sum += items[j].probability;
    // see if the length of the letters is expected
    if (letters.length != indexes.length) {
      return {ok: false, line: items[j].line, 
        reason: "missing or erroneous transition name"};
    }
    for (i = 0; i < letters.length; i++) {
      if (letters[i] != alphabet[indexes[i]]) {
        return {ok: false, line: items[j].line, 
          reason: "missing or erroneous transition name"};
      }
    }
    //increment the indexes starting from the right
    for (i = indexes.length - 1; i >= 0; i--) {
      indexes[i]++;
      if (indexes[i] >= alphabet.length) {
        indexes[i] = 0;
        continue;
      }
      break;
    }
    if (i == -1) {
      // all indexes overflowed so add another one. Technically we should 
      // add it at the start but as all existing entries are zero we
      // actually just add another 0 to the end.
      indexes.push(0);
      // test that the probabilities for the markov states approximately sum to 1
      if (Math.abs(prob_sum - 1.0) > delta) {
        return {ok: false, line: items[j].line, 
          reason: "order-" + (indexes.length - 2) + 
            " probabilities do not sum to 1"};
      }
      order_complete = true;
      prob_sum = 0;
    } else {
      order_complete = false;
    }
  }
  if (!order_complete) {
    // not enough states for a markov model of this order
    return {ok: false, line: lines.length, 
      reason: "not enough transitions for a complete order-" + 
      (indexes.length - 1) + " markov model background"};
  }
  return {ok: true, alphabet: alphabet.join("")};
}

function setup_enabopt() {
  "use strict";
  function _helper(target) {
    "use strict";
    var checkbox, label, input;
    // check if this has already been activated
    if (/\bactive\b/.test(target.className)) return;
    // get the components
    checkbox = target.querySelector("input[type=checkbox]");
    if (checkbox == null) return;
    label = target.querySelector("label");
    if (label == null) return;
    input = document.getElementById(label.htmlFor);
    if (input == null) return;
    // add listeners
    checkbox.addEventListener("click", function() {
      input.disabled = !checkbox.checked;
    }, false);
    label.addEventListener("click", function() {
      checkbox.click();
    }, false);
    // handle form resets
    checkbox.form.addEventListener("reset", function() {
      window.setTimeout(function() {
        input.disabled = !checkbox.checked;
      }, 50);
    }, false);
    // set to current state
    input.disabled = !checkbox.checked;
    // mark this as active
    target.className += " active";
  }
  var targets, i;
  targets = document.querySelectorAll(".enabopt");
  for (i = 0; i < targets.length; i++) {
    _helper(targets[i]);
  }
}

function setup_intonly() {
  "use strict";
  var targets, i;
  targets = document.querySelectorAll("input.intonly");
  for (i = 0; i < targets.length; i++) {
    targets[i].addEventListener("keypress", int_key, false); 
  }
}

function setup_numonly() {
  "use strict";
  var targets, i;
  targets = document.querySelectorAll("input.numonly");
  for (i = 0; i < targets.length; i++) {
    targets[i].addEventListener("keypress", num_key, false); 
  }
}

// anon function to avoid polluting global scope
(function() {
  "use strict";
  window.addEventListener("load", function load(evt) {
    window.removeEventListener("load", load, false);
    setup_help();
    setup_enabopt();
    setup_intonly();
    setup_numonly();
  }, false);

  // fix array.fill for typed arrays
  // adapted from https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Array/fill#Polyfill
  var array_fill = function(value) {

    // Steps 1-2.
    if (this == null) {
      throw new TypeError('this is null or not defined');
    }

    var O = Object(this);

    // Steps 3-5.
    var len = O.length >>> 0;

    // Steps 6-7.
    var start = arguments[1];
    var relativeStart = start >> 0;

    // Step 8.
    var k = relativeStart < 0 ?
      Math.max(len + relativeStart, 0) :
      Math.min(relativeStart, len);

    // Steps 9-10.
    var end = arguments[2];
    var relativeEnd = end === undefined ?
      len : end >> 0;

    // Step 11.
    var final = relativeEnd < 0 ?
      Math.max(len + relativeEnd, 0) :
      Math.min(relativeEnd, len);

    // Step 12.
    while (k < final) {
      O[k] = value;
      k++;
    }

    // Step 13.
    return O;
  };
  if (!Int8Array.prototype.fill) {
    Int8Array.prototype.fill = array_fill;
  }
  if (!Uint8Array.prototype.fill) {
    Uint8Array.prototype.fill = array_fill;
  }
  if (!Int16Array.prototype.fill) {
    Int16Array.prototype.fill = array_fill;
  }
  if (!Uint16Array.prototype.fill) {
    Uint16Array.prototype.fill = array_fill;
  }
  if (!Int32Array.prototype.fill) {
    Int32Array.prototype.fill = array_fill;
  }
  if (!Uint32Array.prototype.fill) {
    Uint32Array.prototype.fill = array_fill;
  }
  if (!Float32Array.prototype.fill) {
    Float32Array.prototype.fill = array_fill;
  }
  if (!Float64Array.prototype.fill) {
    Float64Array.prototype.fill = array_fill;
  }
  if (!Uint8ClampedArray.prototype.fill) {
    Uint8ClampedArray.prototype.fill = array_fill;
  }

})();
