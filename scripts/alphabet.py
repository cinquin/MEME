#!@WHICHPYTHON@

from __future__ import with_statement
import collections, itertools, json, math, re, string, sys
from xml.sax.saxutils import quoteattr

LabColour = collections.namedtuple('LabColour', ['l', 'a', 'b'])

def alphabetSymbolCompare(x, y):
    """Compares 2 strings of alphabet symbols"""
    if len(x) == len(y):
        for a, b in itertools.izip(x, y):
            # letters before numbers and symbols
            if a.isalpha():
                if b.isalpha():
                    if a < b:
                        return -1
                    elif b < a:
                        return 1
                    continue
                else:
                    return -1
            elif b.isalpha():
                return 1
            # numbers before symbols
            if a.isdigit():
                if b.isdigit():
                    if a < b:
                        return -1
                    elif b < a:
                        return 1
                    continue
                else:
                    return -1
            elif b.isdigit():
                return 1
            # both must be symbols, default to string comparison
            if a < b:
                return -1
            if b < a:
                return 1
        # no difference
        return 0
    else:
        return len(y) - len(x)

def _rgb2hsl(rgb):
    """Create HSL from an RGB integer"""
    r = ((rgb >> 16) & 0xFF) / 255
    g = ((rgb >> 8) & 0xFF) / 255
    b = (rgb & 0xFF) / 255
    min_c = min(r, g, b)
    max_c = max(r, g, b)
    delta_c = max_c - min_c
    l = min_c + (delta_c / 2)
    h = 0
    s = 0
    if max_c != min_c:
        if l > 0.5:
            s = delta_c / (2 - max_c - min_c)
        else:
            s = delta_c / (max_c + min_c)
        if max_c == r:
            h = (g - b) / delta_c
            if g < b:
                h += 6
        elif max_c == g:
            h = ((b - r) / delta_c) + 2
        else:
            h = ((r - g) / delta_c) + 4
        h /= 6
    return (h, s, l)

def _hsl2rgb(hsl):
    """Create an RGB integer from a HSL tuple"""
    def _hue(p, q, t):
        if t < 0:
            t += 1
        elif t > 1:
            t -= 1
        if t < (1.0 / 6.0):
            return p + ((q - p) * 6.0 * t)
        elif t < 0.5:
            return q
        elif t < (2.0 / 3.0):
            return p + ((q - p) * ((2.0 / 3.0) - t) * 6.0)
        else:
            return p
    h = hsl[0]
    s = hsl[1]
    l = hsl[2]
    if s == 0:
        # grayscale
        r = l
        g = l
        b = l
    else:
        if l < 0.5:
            q = l * (1 + s)
        else:
            q = l + s - (l * s)
        p = (2 * l) - q
        r = _hue(p, q, h + (1.0 / 3.0))
        g = _hue(p, q, h)
        b = _hue(p, q, h - (1.0 / 3.0))
    rgb = (int(round(r * 255)) << 16) | (int(round(g * 255)) << 8) | int(round(b * 255))
    return rgb

def _lighten(rgb):
    """Make a lighter version of the colour"""
    hsl = _rgb2hsl(rgb)
    hsl = (hsl[0], hsl[1], hsl[2] + ((1.0 - hsl[2]) * 2 / 3))
    return _hsl2rgb(hsl)

def _hsv2rgb(hue, sat, value):
    """Create an RGB integer from HSV"""
    # achromatic (grey)
    r = value
    g = value
    b = value
    if sat != 0:
        h = hue / 60
        i = math.floor(h)
        f = h - i
        p = value * (1.0 - sat)
        q = value * (1.0 - (sat * f))
        t = value * (1.0 - (sat * (1.0 - f)))
        if i == 0:
            r = value
            g = t
            b = p
        elif i == 1:
            r = q
            g = value
            b = p
        elif i == 2:
            r = p
            g = value
            b = t
        elif i == 3:
            r = p
            g = q
            b = value
        elif i == 4:
            r = t
            g = p
            b = value
        else:
            r = value
            g = p
            b = q
    red = int(round(r * 255))
    green = int(round(g * 255))
    blue = int(round(b * 255))
    rgb = (red << 16) | (green << 8) | blue
    return rgb

def _rgb2lab(rgb):
    """Convert an RGB integer to Lab tuple"""
    def xyzHelper(value):
        """Helper function for XYZ colourspace conversion"""
        c = value / 255
        if c > 0.0445:
            c = (c + 0.055) / 1.055
            c = math.pow(c, 2.4)
        else:
            c /= 12.92
        c *= 100
        return c
    def labHelper(value):
        """Helper function for Lab colourspace conversion"""
        c = value
        if c > 0.008856:
            c = math.pow(c, 1.0 / 3.0)
        else:
            c = (7.787 * c) + (16.0 / 116.0)
        return c
    # convert into XYZ colourspace
    c1 = xyzHelper((rgb >> 16) & 0xFF)
    c2 = xyzHelper((rgb >> 8) & 0xFF)
    c3 = xyzHelper(rgb & 0xFF)
    x = (c1 * 0.4124) + (c2 * 0.3576) + (c3 * 0.1805)
    y = (c1 * 0.2126) + (c2 * 0.7152) + (c3 * 0.0722)
    z = (c1 * 0.0193) + (c2 * 0.1192) + (c3 * 0.9505)
    # convert into Lab colourspace
    c1 = labHelper(x / 95.047)
    c2 = labHelper(y / 100.0)
    c3 = labHelper(z / 108.883)
    l = (116.0 * c2) - 16
    a = 500.0 * (c1 - c2)
    b = 200.0 * (c2 - c3)
    return LabColour(l, a, b)

def _lab_dist(lab1, lab2):
    """Calculate the distance between 2 colours in Lab colourspace"""
    c1 = math.sqrt((lab1.l * lab1.l) + (lab1.a * lab1.a))
    c2 = math.sqrt((lab2.l * lab2.l) + (lab2.a * lab2.a))
    dc = c1 - c2
    dl = lab1.l - lab2.l
    da = lab1.a - lab2.a
    db = lab1.b - lab2.b
    # we don't want NaN due to rounding errors so fudge things a bit...
    dh = 0
    dh_squared = (da * da) + (db * db) - (dc * dc)
    if dh_squared > 0:
        dh = math.sqrt(dh_squared)
    first = dl;
    second = dc / (1.0 + (0.045 * c1))
    third = dh / (1.0 + (0.015 * c1))
    return math.sqrt((first * first) + (second * second) + (third * third))

class AlphabetParseError(Exception):
    """An error in the alphabet being parsed"""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class AlphabetReaderSymbol(object):
    """Temporary information storage for an alphabet reader"""
    def __init__(self, symbol, name, colour, complement = None, comprise = None):
        self.symbol = symbol
        self.name = name
        self.colour = colour
        self.complement = complement
        self.comprise = comprise
        self.aliases = []
        self.alias = False

    def __cmp__(self, obj):
        if obj == None:
            return -1
        if not(isinstance(obj, AlphabetReaderSymbol)):
            return -1
        if self.comprise == None or len(self.comprise) == 1:
            if obj.comprise == None or len(obj.comprise) == 1:
                return alphabetSymbolCompare(self.symbol, obj.symbol)
            else:
                return 1
        else:
            if obj.comprise == None or len(obj.comprise) == 1:
                return -1
            else:
                return alphabetSymbolCompare(self.comprise, obj.comprise)


    def __str__(self):
        return repr(self.symbol)

    def isAmbiguous(self):
        return self.comprise != None and len(self.comprise) > 1;

class AlphabetReader(object):
    """Reads an alphabet and validates it before creating an alphabet object"""

    def __init__(self):
        # define some regular expression parts
        qr_sym = "[A-Za-z0-9?.*-]"
        qr_colour = "[0-9a-fA-F]{6}"
        qr_name = '("(?:[^\\\\"]+|\\\\["\\\\/bfnrt]|\\\\u[0-9A-Fa-f]{4})*")'
        qr_core = '(' + qr_sym + ')(?:\\s+' + qr_name + ')?(?:\\s+(' + qr_colour + '))?'
        # define some regular expressions
        self.header_re = re.compile('^\\s*ALPHABET(?:\\s+v1)?(?:\\s+' + qr_name + ')?(?:\\s+(DNA|RNA|PROTEIN)-LIKE)?\\s*$')
        self.core_single_re = re.compile('^\\s*' + qr_core + '\\s*$')
        self.core_pair_re = re.compile('^\\s*' + qr_core + '\\s*~\\s*' + qr_core + '\\s*$')
        self.ambig_re = re.compile('^\\s*' + qr_core + '\\s*=\\s*(' + qr_sym + '*)\\s*$')
        # track state
        self.parsed = False
        self.seen_header = False
        self.seen_symbol = False
        self.seen_ambig = False
        self.seen_lc = False
        self.seen_uc = False
        self.fully_complementable = False
        self.name = None
        self.like = None
        self.sym_lookup = {}
        self.core_syms = []
        self.ambig_syms = []

    def _checkComplements(self):
        """Ensure all referenced complement symbols actually exist."""
        fully_complementable = True
        for sym_obj in self.core_syms:
            if sym_obj.complement == None:
                fully_complementable = False
                continue
            complement_sym_obj = self.sym_lookup[sym_obj.complement]
            if complement_sym_obj == None:
                raise AlphabetParseError("core symbol " + sym_obj.symbol + " has complement " + sym_obj.complement + " which has not been defined")
            if complement_sym_obj.complement != sym_obj.symbol:
                raise RuntimeError("Alphabet symbol complement is not symmetric")
        self.fully_complementable = fully_complementable

    def _mergeAliases(self):
        """Check every ambiguity character to see if it has the same comprising
        characters as something else"""
        # make a list of symbols for each comprising character
        comprise_lookup = {}
        for sym in self.ambig_syms:
            if sym.comprise in comprise_lookup:
                comprise_lookup[sym.comprise].append(sym)
            else:
                comprise_lookup[sym.comprise] = [sym]
        for comprise, syms in comprise_lookup.iteritems():
            if len(comprise) == 1:
                # alias for a core symbol
                core_sym = self.sym_lookup[comprise]
                for sym in syms:
                    sym.alias = True
                    core_sym.aliases.append(sym.symbol)
            elif len(syms) > 1:
                # aliases for an ambiguous symbol
                prime_ambig = None
                for sym in syms:
                    if prime_ambig == None or alphabetSymbolCompare(sym.symbol, prime_ambig.symbol) < 0:
                        prime_ambig = sym
                for sym in syms:
                    if not(sym is prime_ambig):
                        sym.alias = True
                        prime_ambig.aliases.append(sym.symbol)
        # remove the aliases
        self.ambig_syms = [sym for sym in self.ambig_syms if not(sym.alias)]
        # sort the symbols
        self.core_syms.sort()
        self.ambig_syms.sort()

    def _createWildcardIfMissing(self):
        """Tries to find a wildcard but if none exists then creates one"""
        if len(self.ambig_syms) == 0 or len(self.ambig_syms[0].comprise) != len(self.core_syms):
            # print >> sys.stderr, ("Warning: wildcard symbol automatically generated")
            wildcard = AlphabetReaderSymbol('?', None, None, comprise = ''.join([sym.symbol for sym in self.core_syms]))
            self.ambig_syms.insert(0, wildcard)
            self.sym_lookup['?'] = wildcard

    def _setMissingColours(self):
        """Assigns colours to any core symbols that don't have one"""
        # count how many unique colours we need
        ncolours = 0
        unique_colours = set()
        for sym_obj in self.core_syms:
            if sym_obj.colour != None:
                unique_colours.add(sym_obj.colour)
            else:
                ncolours += 1
        for sym_obj in self.ambig_syms:
            if sym_obj.colour != None:
                unique_colours.add(sym_obj.colour)
        ncolours += len(unique_colours)
        # generate Lab colourspace versions of each assigned colour
        uniques = []
        for rgb in unique_colours:
            uniques.append(_rgb2lab(rgb))
        # generate evenly distributed colours in HSV colourspace and convert to RGB and Lab
        sat = 1.0
        value = 0.4
        step = 360 / ncolours
        colours = []
        for i in xrange(ncolours):
            rgb = _hsv2rgb(step * i, sat, value)
            colours.append((rgb, _rgb2lab(rgb)))
        # for each of the unique colours find the closest colour in the generated colours and remove it
        while uniques:
            best_dist = None
            best_i = None
            best_j = None
            for i in xrange(len(uniques)):
                for j in xrange(len(colours)):
                    dist = _lab_dist(uniques[i], colours[j][1])
                    if best_dist == None or dist < best_dist:
                        best_dist = dist
                        best_i = i
                        best_j = j
            if best_dist == None:
                raise RuntimeError("Somehow we ran out of colours?!")
            del uniques[best_i]
            del colours[best_j]
        # assign the colours
        for sym_obj in self.core_syms:
            if sym_obj.colour == None:
                sym_obj.colour = colours.pop(0)[0]
        for sym_obj in self.ambig_syms:
            if sym_obj.colour == None:
                sym_obj.colour = 0


    def _setMissingComplements(self):
        """Tries to find a complement for every ambiguous character"""
        comprise_lookup = {}
        for sym_obj in self.ambig_syms:
            comprise_lookup[sym_obj.comprise] = sym_obj
        for sym_obj in self.ambig_syms:
            comprise = sym_obj.comprise
            complement_list = []
            for sym in comprise:
                component_obj = self.sym_lookup[sym]
                if component_obj.complement == None:
                    break
                complement_list.append(component_obj.complement)
            else:
                # all the symbols had complements
                complement = ''.join(sorted(complement_list, cmp=alphabetSymbolCompare))
                if complement in comprise_lookup:
                    sym_obj.complement = comprise_lookup[complement].symbol

    def _checkLikeStandards(self):
        """Confirms that an alphabet is an extension of a specified standard alphabet"""
        if self.like != None:
            required_symbols = None
            required_complements = None
            if self.like == "RNA":
                required_symbols = "ACGU"
            elif self.like == "DNA":
                required_symbols = "ACGT"
                required_complements = "TGCA"
            elif self.like == "PROTEIN":
                required_symbols = "ACDEFGHIKLMNPQRSTVWY"
            if required_symbols != None:
                for i, sym in enumerate(required_symbols):
                    if not sym in self.sym_lookup:
                        raise AlphabetParseError("alphabet is not " + self.like + "-like; missing symbol " + sym)
                    sym_obj = self.sym_lookup[sym]
                    if sym_obj.comprise != None:
                        raise AlphabetParseError("alphabet is not " + self.like + "-like; symbol " + sym + " is ambiguous")
                    comp1_obj = None
                    if required_complements != None and required_complements[i] in self.sym_lookup:
                        comp1_obj = self.sym_lookup[required_complements[i]]
                    comp2_obj = None
                    if sym_obj.complement != None:
                        comp2_obj = self.sym_lookup[sym_obj.complement]
                    if not (comp1_obj is comp2_obj):
                        raise AlphabetParseError("alphabet is not " + self.like + "-like; symbol " + sym + " complement rules are incorrect")

    def parseHeader(self, name = None, like = None):
        if self.parsed:
            raise RuntimeError("Parsing is already done!")
        if self.seen_header:
            raise AlphabetParseError("repeated header")
        if self.seen_symbol:
            raise AlphabetParseError("header after symbol")
        if like != None and like != "DNA" and like != "RNA" and like != "PROTEIN":
            raise RuntimeError("If defined then \"like\" must be DNA, RNA or PROTEIN")
        self.name = name
        self.like = like
        self.seen_header = True

    def parseSymbol(self, symbol, name = None, colour = None, complement = None, comprise = None):
        # check assumptions
        if self.parsed:
            raise RuntimeError("Parsing is already done!")
        if not(isinstance(symbol, basestring) and len(symbol) == 1):
            raise RuntimeError("Expected symbol to be a single character")
        if name != None and not(isinstance(name, basestring)):
            raise RuntimeError("Expected name to be a string")
        if colour != None and not(isinstance(colour, (int, long))):
            raise RuntimeError("Expected colour to be a unsigned 24 bit number")
        if complement != None and not(isinstance(complement, basestring) and len(complement) == 1):
            raise RuntimeError("Expected complement symbol to be a single character")
        if comprise != None and not(isinstance(comprise, basestring)):
            raise RuntimeError("Expected name to be a string")
        if complement != None and comprise != None:
            raise RuntimeError("Expected only comprise or complement not both")
        # check state of reader
        if not(self.seen_header or self.seen_symbol):
            raise AlphabetParseError("expected header but found symbol")
        # check that the symbol is not yet defined
        if symbol in self.sym_lookup:
            raise AlphabetParseError("symbol " + str(symbol) + " is already used")
        # check any comprising symbols are defined, 
        # we can do this because core symbols must be defined before ambiguous symbols
        if comprise != None:
            # sort
            compriseList = sorted(comprise, cmp=alphabetSymbolCompare)
            # filter out duplicates
            comprise = ''.join([e for i, e in enumerate(compriseList) if i == 0 or compriseList[i-1] != e])
            for sym in comprise:
                comprise_obj = self.sym_lookup.get(sym)
                if comprise_obj == None:
                    raise AlphabetParseError("referenced symbol " + sym + " is unknown")
                elif comprise_obj.isAmbiguous():
                    raise AlphabetParseError("referenced symbol " + sym + " is ambiguous")
        elif self.seen_ambig:
            raise AlphabetParseError("unexpected core symbol (as ambiguous symbols seen)")
        # enforce the restriction that ? may only be a wildcard
        if symbol == "?":
            # check that this is the wildcard
            if comprise == None or comprise != ''.join(sorted([sym.symbol for sym in self.core_syms], cmp=alphabetSymbolCompare)):
                raise AlphabetParseError("symbol ? is reserved for wildcards only")
        # check for lower and upper case letters
        if symbol.islower():
            self.seen_lc = True
        elif symbol.isupper():
            self.seen_uc = True
        # create a symbol 
        symbol_obj = AlphabetReaderSymbol(symbol, name, colour, complement, comprise)
        if comprise == None:
            self.core_syms.append(symbol_obj)
        else:
            # first time we see an ambiguous symbol we check that all complements are defined
            if not(self.seen_ambig):
                self._checkComplements()
            self.ambig_syms.append(symbol_obj)
            self.seen_ambig = True
        self.sym_lookup[symbol] = symbol_obj
        self.seen_symbols = True

    def parseLine(self, line):
        def removeComments(line):
            i = 0
            # skip leading space
            while i < len(line):
                if not(line[i].isspace()):
                    break
                i += 1
            if i == len(line):
                return "" # line only contains space
            # test first non-space symbol (can not start string)
            if line[i] == "#":
                return "" # discard whole line
            i += 1
            # test following letters - have to check for string
            while i < len(line):
                c = line[i]
                if line[i] == "#":
                    return line[0:i] # found comment, trim line
                if line[i] == '"':
                    i += 1 # found string, search for end
                    while i < len(line):
                        if line[i] == '"' and line[i-1] != '\\':
                            break # end of string
                        i += 1
                i += 1
            return line # no comment in line
        def decodeName(name):
            if name != None:
                return json.loads(name)
            return None
        def decodeColour(colour):
            if colour != None:
                return int(colour, 16)
            return None
        # remove comments
        line = removeComments(line)
        # skip empty lines
        if len(line) == 0 or line.isspace():
            return
        # match the header
        match = self.header_re.match(line)
        if match != None:
            self.parseHeader(decodeName(match.group(1)), match.group(2))
            return
        # match a complementary pair of core symbols
        match = self.core_pair_re.match(line)
        if match != None:
            self.parseSymbol(match.group(1), decodeName(match.group(2)), decodeColour(match.group(3)), complement=match.group(4))
            self.parseSymbol(match.group(4), decodeName(match.group(5)), decodeColour(match.group(6)), complement=match.group(1))
            return
        # match a single core symbol
        match = self.core_single_re.match(line)
        if match != None:
            self.parseSymbol(match.group(1), decodeName(match.group(2)), decodeColour(match.group(3)))
            return
        # match an ambiguous or alias symbol
        match = self.ambig_re.match(line)
        if match != None:
            self.parseSymbol(match.group(1), decodeName(match.group(2)), decodeColour(match.group(3)), comprise=match.group(4))
            return
        # nothing matched
        raise AlphabetParseError("unrecognised pattern")
        
    def parseDone(self):
        if not(self.seen_ambig):
            self._checkComplements()
        self._mergeAliases()
        self._createWildcardIfMissing()
        self._setMissingComplements()
        self._setMissingColours()
        self._checkLikeStandards()
        self.parsed = True

    def parseFile(self, filename):
        with open(filename) as fh:
            for line in fh:
                self.parseLine(line.strip())
        self.parseDone()
        return Alphabet(self)

class AlphabetSymbol(object):
    """A symbol of the alphabet"""
    def __init__(self, index, symbol, aliases, name, colour):
        self.index = index
        self.symbol = symbol
        self.aliases = list(aliases)
        self.name = name
        self.colour = colour
        # the index of the complement symbol
        self.complement = None
        # a frozen set of indexes of core comprising symbols (core symbols use their own index)
        self.comprise = None
        # a tuple of 2 indexes of smaller symbols (possibly ambiguous) that when combined make this one
        self.pair = None

    def __eq__(self, other):
        if isinstance(other, AlphabetSymbol):
            if self.index != other.index:
                return False
            if self.symbol != other.symbol:
                return False
            if self.aliases != other.aliases:
                return False
            if self.complement != other.complement:
                return False
            if self.comprise != other.comprise:
                return False
            return True
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def xmlid(self):
        """Provide an XML compatible ID for a symbol"""
        if (self.symbol >= 'A' and self.symbol <= 'Z') or (self.symbol >= 'a' and self.symbol <= 'z'):
            return self.symbol
        elif self.symbol >= '0' and self.symbol <= '9':
            return 'n' + self.symbol
        else:
            return 'x{:02X}'.format(ord(self.symbol))

    def asText(self):
        out = [self.symbol]
        if self.name != None:
            out.append(" ")
            out.append(json.dumps(self.name))
        if self.colour != 0:
            out.append(" {:06X}".format(self.colour))
        return "".join(out)

                

class Alphabet(object):
    """Biological alphabet class.
    This defines the set of symbols from which various objects can be built, e.g. sequences and motifs.
    """

    def __init__(self, reader):
        if not(reader.parsed):
            raise RuntimeError("The reader must finish parsing before an alphabet can be created from it")
        # create a symbol object for each symbol in the reader
        self.ncore = len(reader.core_syms)
        self.ignore_case = reader.seen_lc != reader.seen_uc
        self.symbols = []
        self.lookup = {}
        self.find = {}
        self.comp_symbols = {}
        self.name = reader.name
        self.like = reader.like
        self.fully_complementable = reader.fully_complementable
        index = 0
        for sym_obj in itertools.chain(reader.core_syms, reader.ambig_syms):
            alph_sym = AlphabetSymbol(index, sym_obj.symbol, sym_obj.aliases, sym_obj.name, sym_obj.colour)
            self.symbols.append(alph_sym)
            for sym in itertools.chain(sym_obj.symbol, sym_obj.aliases):
                self.lookup[sym] = alph_sym
                if self.ignore_case:
                    if (sym >= 'A' and sym <= 'Z') or (sym >= 'a' and sym <= 'z'):
                        self.lookup[sym.swapcase()] = alph_sym
            index += 1
        # link complementary symbols
        for sym_obj in itertools.chain(reader.core_syms, reader.ambig_syms):
            if sym_obj.complement != None:
                alph_sym_1 = self.lookup[sym_obj.symbol]
                alph_sym_2 = self.lookup[sym_obj.complement]
                alph_sym_1.complement = alph_sym_2.index
        # link comprising symbols
        # core symbols have themselves in the comprise list
        for i in xrange(self.ncore):
            self.symbols[i].comprise = frozenset([i])
        # ambiguous symbols link to the core symbols
        for sym_obj in reader.ambig_syms:
            alph_sym = self.lookup[sym_obj.symbol]
            comprise_list = []
            for comprise_symbol in sym_obj.comprise:
                comprise_list.append(self.lookup[comprise_symbol].index)
            alph_sym.comprise = frozenset(comprise_list)
        # make a lookup for the comprising symbol sets
        for alph_sym in self.symbols:
            self.find[alph_sym.comprise] = alph_sym
        # make a tuple containing all core symbols
        self.core = tuple([sym.symbol for sym in self.symbols if sym.index < self.ncore])
        # make a tuple containing indexes of all complementary symbols
        if self.fully_complementable:
            self.complements = tuple([sym.complement for sym in self.symbols if sym.index < self.ncore])
        else:
            self.complements = None
        # make a lookup comprising all complements (for SPEED)
        if self.fully_complementable:
            for alph_sym in self.symbols:
              self.comp_symbols[alph_sym.symbol] = self.getComplement(alph_sym.symbol)
        # try to determine a pair of symbols that make up each ambiguous symbol
        symbols_by_count = [[] for _ in xrange(self.ncore + 1)]
        for alph_sym in self.symbols:
            symbols_by_count[len(alph_sym.comprise)].append(alph_sym)
        for alph_sym in self.symbols[self.ncore:]:
            for j in xrange(1, len(alph_sym.comprise)):
                for part in symbols_by_count[j]:
                    if not(part.comprise < alph_sym.comprise):
                        continue
                    remain = alph_sym.comprise - part.comprise
                    if remain in self.find:
                        alph_sym.pair = (part.index, self.find[remain].index)
                        break
                if alph_sym.pair != None:
                    break

    def __repr__(self):
        if self.name != None and len(self.name) > 0:
            return self.name
        return "".join(self.core)

    def __eq__(self, other):
        if isinstance(other, Alphabet):
            if self.ncore != other.ncore:
                return False
            if self.ignore_case != other.ignore_case:
                return False
            if self.symbols != other.symbols:
                return False
            return True
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def getColour(self, sym):
        """Get the colour for this symbol"""
        if isinstance(sym, (int, long)):
            if sym < len(self.symbols):
                return self.symbols[sym].colour
            raise RuntimeError("Symbol index " + sym + " does not exist in alphabet")
        else:
            if sym in self.lookup:
                return self.lookup[sym].colour
            raise RuntimeError("Symbol " + sym + " does not exist in alphabet")

    def getMutedColour(self, sym):
        """Get the muted colour for this symbol"""
        if isinstance(sym, (int, long)):
            if sym < len(self.symbols):
                return _lighten(self.symbols[sym].colour)
            raise RuntimeError("Symbol index " + str(sym) + " does not exist in alphabet")
        else:
            if sym in self.lookup:
                return _lighten(self.lookup[sym].colour)
            raise RuntimeError("Symbol " + sym + " does not exist in alphabet")


    def getWildcard(self):
        """Get the wildcard symbol"""
        if len(self.symbols) > self.ncore:
            return self.symbols[self.ncore].symbol
        else:
            return None

    def getSymbols(self):
        """Retrieve a tuple with all core symbols, immutable membership and order"""
        return self.core

    def getComplements(self):
        """Retrieve a tuple with all core symbol complement indicies, immutable"""
        return self.complements

    def isValidSymbol(self, sym):
        """Check if the symbol is a member of alphabet"""
        return sym in self.lookup

    def isCoreSymbol(self, sym):
        """Check if the symbol is a core member of alphabet"""
        if sym in self.lookup:
            return self.lookup[sym].index < self.ncore
        return False

    def isWildcardSymbol(self, sym):
        """Check if the symbol is the wildcard of the alphabet"""
        if sym in self.lookup:
            return self.lookup[sym].index == self.ncore
        return False

    def isPrimeSymbol(self, sym):
        """Check if the symbol is a member of the alphabet and the normal way to identify the letter"""
        if sym in self.lookup:
            if self.lookup[sym].symbol.lower() == sym.lower():
                return True
        return False

    def getIndex(self, sym):
        """Retrieve the index of the symbol (immutable)"""
        if sym in self.lookup:
            return self.lookup[sym].index
        return None

    def getComprisingIndexes(self, sym):
        """Retrieve the indexes of the comprising core symbols (immutable)"""
        if isinstance(sym, (int, long)):
            if sym < len(self.symbols):
                return self.symbols[sym].comprise
        else:
            if sym in self.lookup:
                return self.lookup[sym].comprise
        return None

    def getPairIndexes(self, sym):
        """Retrieve the indexes of the comprising core symbols (immutable)"""
        if isinstance(sym, (int, long)):
            if sym < len(self.symbols):
                return self.symbols[sym].pair
        else:
            if sym in self.lookup:
                return self.lookup[sym].pair
        return None


    def getSymbol(self, index):
        """Return the symbol at the given index"""
        if index < len(self.symbols):
            return self.symbols[index].symbol
        return None

    def getXmlId(self, index):
        if index < len(self.symbols):
            return self.symbols[index].xmlid()
        return None

    def getAliases(self, index):
        """Returns all symbols variants at the given index"""
        alph_sym = self.symbols[index]
        syms = []
        for sym in itertools.chain(alph_sym.symbol, alph_sym.aliases):
            syms.append(sym)
            if self.ignore_case and sym.isalpha():
                syms.append(sym.swapcase())
        return "".join(syms)


    def findSymbol(self, comprise):
        """Find a symbol that matches the string of comprising symbols.
        Duplicates and aliases are allowed in the string but all
        included symbols must be known."""
        symi_set = set()
        for sym in comprise:
            if not sym in self.lookup:
                return None
            symi_set |= self.lookup[sym].comprise
        symi_set = frozenset(symi_set)
        if symi_set in self.find:
            return self.find[symi_set].symbol
        return None

    def isComplementable(self):
        """Do all core symbols have complements"""
        return self.fully_complementable

    def getComplementIndex(self, sym):
        """Retrieve the index of the complement of the symbol"""
        if isinstance(sym, (int, long)):
            if sym < len(self.symbols):
                return self.symbols[sym].complement
        else:
            if sym in self.lookup:
                return self.lookup[sym].complement
        return None


    def getComplement(self, sym):
        """Retrieve the complement of the symbol"""
        complementi = self.getComplementIndex(sym)
        if complementi != None:
            return self.symbols[complementi].symbol
        return None

    def getComplementFast(self, sym):
        """Retrieve the complement of the symbol faster"""
        return self.comp_symbols[sym]

    def isValidString(self, symstr):
        """Check if the string contains only symbols that belong to the alphabet"""
        for sym in symstr:
            if not(self.isValidSymbol(sym)):
                return False
        return True

    def getLen(self):
        """Retrieve the number of core symbols in the alphabet"""
        return self.ncore

    def getWildLen(self):
        """Retrieve the number of core symbols plus the wildcard in the alphabet"""
        return min(len(self.symbols), self.ncore + 1)

    def getFullLen(self):
        """Retrieve the full count of core and ambiguous symbols in the alphabet"""
        return len(self.symbols)

    def encodeString(self, symstr):
        """Encode the symstr as indexes of the alphabet"""
        return [self.getIndex(sym) for sym in symstr]

    def getName(self):
        """Get the alphabet name"""
        if self.name != None and len(self.name) > 0:
            return self.name
        else:
            return "".join(self.core)

    def translator(self, allow_ambig = True):
        """Create a translation table that will convert sequences to the prime symbol.
        Optionally convert ambiguous symbols to the wildcard symbol."""
        src = []
        dest = []
        for alph_sym in self.symbols:
            dsym = alph_sym.symbol if alph_sym.index <= self.ncore or allow_ambig else self.getWildcard()
            if self.ignore_case and alph_sym.symbol.isalpha():
                src.append(alph_sym.symbol.swapcase())
                dest.append(dsym)
            for sym in alph_sym.aliases:
                src.append(sym)
                dest.append(dsym)
                if self.ignore_case and sym.isalpha():
                    src.append(sym.swapcase())
                    dest.append(dsym)
        return string.maketrans("".join(src), "".join(dest))


    def asText(self):
        """Create the text representation of the alphabet"""
        out = []
        # output core symbols with complements
        for i in xrange(self.ncore):
            sym1 = self.symbols[i]
            # check symbol has a complement (which we assume is symetric)
            if sym1.complement == None:
                continue
            sym2 = self.symbols[sym1.complement]
            if sym2.index < i:
                continue
            out.append(sym1.asText())
            out.append(" ~ ")
            out.append(sym2.asText())
            out.append("\n")
        # output core symbols without complements
        for i in xrange(self.ncore):
            sym = self.symbols[i]
            if sym.complement != None:
                continue
            out.append(sym.asText())
            out.append("\n")
        # output ambiguous symbols
        for i in xrange(self.ncore, len(self.symbols)):
            sym = self.symbols[i]
            comprise_str = "".join(sorted([self.symbols[symi].symbol for symi in sym.comprise], cmp=alphabetSymbolCompare))
            out.append(sym.asText())
            out.append(" = ")
            out.append(comprise_str)
            out.append("\n")
            if sym.aliases != None:
                for alias in sym.aliases:
                    out.append(alias)
                    out.append(" = ")
                    out.append(comprise_str)
                    out.append("\n")
        # output aliases for core symbols
        for i in xrange(self.ncore):
            sym = self.symbols[i]
            if sym.aliases != None:
                for alias in sym.aliases:
                    out.append(alias)
                    out.append(" = ")
                    out.append(sym.symbol)
                    out.append("\n")
        return "".join(out)
            
    def asXML(self, pad="    ", indent="  "):
        """Create the XML representation of the alphabet"""
        out = []
        out.append("{:s}<alphabet name={:s}".format(pad, quoteattr(self.getName())))
        if self.like != None:
            out.append(" like=\"{:s}\"".format(self.like.lower()))
        out.append(">\n")
        for sym in self.symbols:
            out.append("{:s}{:s}<letter id=\"{:s}\" symbol={:s}".format(pad, indent, sym.xmlid(), quoteattr(sym.symbol)))
            if len(sym.aliases) > 0:
                out.append(" aliases={:s}".format(quoteattr("".join(sym.aliases))))
            if len(sym.comprise) == 1:
                if sym.complement != None:
                    csym = self.symbols[sym.complement];
                    out.append(" complement={:s}".format(quoteattr(csym.symbol)))
            else:
                out.append(" equals={:s}".format(quoteattr("".join(self.symbols[i].symbol for i in sym.comprise))))
            if sym.name != None and len(sym.name) > 0:
                out.append(" name={:s}".format(quoteattr(sym.name)))
            if sym.colour != 0:
                out.append(" colour=\"{:06X}\"".format(sym.colour))
            out.append("/>\n")
        out.append("{:s}</alphabet>\n".format(pad))
        return "".join(out)

def dna():
    """Create a DNA alphabet"""
    factory = AlphabetReader()
    factory.parseHeader("DNA", like="DNA")
    factory.parseSymbol('A', "Adenine", 0xCC0000, complement='T')
    factory.parseSymbol('C', "Cytosine", 0x0000CC, complement='G')
    factory.parseSymbol('G', "Guanine", 0xFFB300, complement='C')
    factory.parseSymbol('T', "Thymine", 0x008000, complement='A')
    factory.parseSymbol('U', "Uracil", comprise="T")
    factory.parseSymbol('W', "Weak", comprise="AT")
    factory.parseSymbol('S', "Strong", comprise="CG")
    factory.parseSymbol('M', "Amino", comprise="AC")
    factory.parseSymbol('K', "Keto", comprise="GT")
    factory.parseSymbol('R', "Purine", comprise="AG")
    factory.parseSymbol('Y', "Pyrimidine", comprise="CT")
    factory.parseSymbol('B', "Not A", comprise="CGT")
    factory.parseSymbol('D', "Not C", comprise="AGT")
    factory.parseSymbol('H', "Not G", comprise="ACT")
    factory.parseSymbol('V', "Not T", comprise="ACG")
    factory.parseSymbol('N', "Any base", comprise="ACGT")
    factory.parseSymbol('X', "Any base", comprise="ACGT")
    factory.parseSymbol('.', "Any base", comprise="ACGT")
    factory.parseDone()
    return Alphabet(factory)

def rna():
    """Create a RNA alphabet"""
    factory = AlphabetReader()
    factory.parseHeader("RNA", like="RNA")
    factory.parseSymbol('A', "Adenine", 0xCC0000)
    factory.parseSymbol('C', "Cytosine", 0x0000CC)
    factory.parseSymbol('G', "Guanine", 0xFFB300)
    factory.parseSymbol('U', "Uracil", 0x008000)
    factory.parseSymbol('T', "Thymine", comprise="U")
    factory.parseSymbol('W', "Weak", comprise="AU")
    factory.parseSymbol('S', "Strong", comprise="CG")
    factory.parseSymbol('M', "Amino", comprise="AC")
    factory.parseSymbol('K', "Keto", comprise="GU")
    factory.parseSymbol('R', "Purine", comprise="AG")
    factory.parseSymbol('Y', "Pyrimidine", comprise="CU")
    factory.parseSymbol('B', "Not A", comprise="CGU")
    factory.parseSymbol('D', "Not C", comprise="AGU")
    factory.parseSymbol('H', "Not G", comprise="ACU")
    factory.parseSymbol('V', "Not U", comprise="ACG")
    factory.parseSymbol('N', "Any base", comprise="ACGU")
    factory.parseSymbol('X', "Any base", comprise="ACGU")
    factory.parseSymbol('.', "Any base", comprise="ACGU")
    factory.parseDone()
    return Alphabet(factory)

def protein():
    """Create a protein alphabet"""
    factory = AlphabetReader()
    factory.parseHeader("Protein", like="PROTEIN")
    factory.parseSymbol('A', "Alanine", 0x0000CC)
    factory.parseSymbol('R', "Arginine", 0xCC0000)
    factory.parseSymbol('N', "Asparagine", 0x008000)
    factory.parseSymbol('D', "Aspartic acid", 0xFF00FF)
    factory.parseSymbol('C', "Cysteine", 0x0000CC)
    factory.parseSymbol('E', "Glutamic acid", 0xFF00FF)
    factory.parseSymbol('Q', "Glutamine", 0x008000)
    factory.parseSymbol('G', "Glycine", 0xFFB300)
    factory.parseSymbol('H', "Histidine", 0xFFCCCC)
    factory.parseSymbol('I', "Isoleucine", 0x0000CC)
    factory.parseSymbol('L', "Leucine", 0x0000CC)
    factory.parseSymbol('K', "Lysine", 0xCC0000)
    factory.parseSymbol('M', "Methionine", 0x0000CC)
    factory.parseSymbol('F', "Phenylalanine", 0x0000CC)
    factory.parseSymbol('P', "Proline", 0xFFFF00)
    factory.parseSymbol('S', "Serine", 0x008000)
    factory.parseSymbol('T', "Threonine", 0x008000)
    factory.parseSymbol('W', "Tryptophan", 0x0000CC)
    factory.parseSymbol('Y', "Tyrosine", 0x33E6CC)
    factory.parseSymbol('V', "Valine", 0x0000CC)
    factory.parseSymbol('B', "Asparagine or Aspartic acid", comprise="ND")
    factory.parseSymbol('Z', "Glutamine or Glutamic acid", comprise="QE")
    factory.parseSymbol('J', "Leucine or Isoleucine", comprise="LI")
    factory.parseSymbol('X', "Any amino acid", comprise="ARNDCEQGHILKMFPSTWYV")
    factory.parseSymbol('*', "Any amino acid", comprise="ARNDCEQGHILKMFPSTWYV")
    factory.parseSymbol('.', "Any amino acid", comprise="ARNDCEQGHILKMFPSTWYV")
    factory.parseDone()
    return Alphabet(factory)

predef = (dna(), rna(), protein())

def getByName(name, alphabets = predef):
    """Retrieve a pre-defined alphabet by name.
    Currently, "Protein", "DNA" and "RNA" are available.
    Example:
    >>> alpha = sequence.getAlphabet('Protein')
    >>> alpha.getSymbols()
    will retrieve the 20 amino acid alphabet and output the tuple:
    ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
    """
    for alphabet in alphabets:
        if alphabet.name == name:
            return alphabet
    return None

def getBySeq(seq, alphabets = predef):
    """Retrieve a pre-defined alphabet by looking at symbol frequency."""
    counts = {}
    for sym in seq:
        if sym in counts:
            counts[sym] += 1
        else:
            counts[sym] = 1
    valid_alphabets = []
    # first filter the alphabets by invalid symbols
    for alphabet in alphabets:
        for sym in counts:
            if not(alphabet.isValidSymbol(sym)):
                break
        else:
            valid_alphabets.append(alphabet)
    if len(valid_alphabets) == 0:
        return None
    elif len(valid_alphabets) == 1:
        return valid_alphabets[0]
    # determine the fraction of the letters 
    best_normal_score = 0
    best_prime_score = 0
    best_alphabets = []
    for alphabet in valid_alphabets:
        seen = set()
        prime_seen = set()
        prime = 0
        alt = 0
        other = 0
        for sym, count in counts.iteritems():
            if alphabet.isCoreSymbol(sym):
                seen.add(alphabet.getIndex(sym))
                if alphabet.isPrimeSymbol(sym):
                    prime_seen.add(alphabet.getIndex(sym))
                    prime += count
                else:
                    alt += count
            elif alphabet.isWildcardSymbol(sym):
                if alphabet.isPrimeSymbol(sym):
                    prime += count
                else:
                    alt += count
            else:
                other += count
        total = prime + alt + other
        normal_score = ((prime + alt) / total) * (len(seen) / alphabet.getLen())
        prime_score = (prime / total) * (len(prime_seen) / alphabet.getLen())
        if normal_score > best_normal_score or (normal_score == best_normal_score and prime_score > best_prime_score):
            best_normal_score = normal_score
            best_prime_score = prime_score
            best_alphabets = [alphabet]
        elif normal_score == best_normal_score and prime_score == best_prime_score:
            best_alphabets.append(alphabet)
    if len(best_alphabets) == 0:
        return valid_alphabets[0]
    return best_alphabets[0]

def loadFromFile(filename):
    reader = AlphabetReader()
    return reader.parseFile(filename)

def main():
    if len(sys.argv) > 1:
        alphabet = loadFromFile(sys.argv[1])

if __name__ == '__main__': main()
