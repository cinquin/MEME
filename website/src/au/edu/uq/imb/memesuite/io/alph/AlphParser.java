package au.edu.uq.imb.memesuite.io.alph;

import au.edu.uq.imb.memesuite.data.Alph;

import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.ParseException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A class for parsing an alphabet.
 * If the alphabet is in text format then it can be parsed line by line with parseLine,
 * otherwise the more specific calls to parseHeader, parseSymbol and parseDone will be needed.
 */
public class AlphParser {
  private static final String SYMBOL_RX = "[A-Za-z0-9\\?\\.\\*\\-]";
  private static final String COLOUR_RX = "[A-Fa-f0-9]{6}";
  private static final String NAME_RX = "\"((?:[^\\\\\"]+|\\\\[\"\\/bfnrt]|\\\\u[0-9A-Fa-f]{4})*)\"";
  private static final String CORE_RX = "(" + SYMBOL_RX + ")(?:\\s+" + NAME_RX + ")?(?:\\s+(" + COLOUR_RX + "))?";
  private static final Pattern COLOUR_RE = Pattern.compile(COLOUR_RX);
  private static final Pattern SYMBOL_RE = Pattern.compile(SYMBOL_RX);
  private static final Pattern SYMBOLS_RE = Pattern.compile(SYMBOL_RX + "+");
  private static final Pattern STARS_RE = Pattern.compile("\\s*\\**\\s*");
  private static final Pattern HEADER_RE = Pattern.compile("^\\s*ALPHABET(?:\\s+v1)?(?:\\s+" + NAME_RX + ")?(?:\\s+(RNA|DNA|PROTEIN)-LIKE)?\\s*$");
  private static final Pattern CORE_SINGLE_RE = Pattern.compile("^\\s*" + CORE_RX + "\\s*$");
  private static final Pattern CORE_PAIR_RE = Pattern.compile("^\\s*" + CORE_RX + "\\s*~\\s*" + CORE_RX + "\\s*$");
  private static final Pattern AMBIG_RE = Pattern.compile("^\\s*" + CORE_RX + "\\s*=\\s*(" + SYMBOL_RX + "*)\\s*$");
  private int lineNo = 0;
  private boolean error = false;
  private boolean seenHeader = false;
  private boolean seenAmbiguous = false;
  private String name;
  private Alph.Like like = null;
  private boolean seenUpperCase = false;
  private boolean seenLowerCase = false;
  private Set<Character> allSymbols = new TreeSet<Character>();
  private SortedSet<Character> coreSymbols = new TreeSet<Character>(Alph.SYM_COMPARATOR);
  private Map<SortedSet<Character>,List<Symbol>> symbolMap = new TreeMap<SortedSet<Character>, List<Symbol>>(new SymSetComparator());

  /**
   * Removes any text after a # as long as it is not in a piece of quoted text
   */
  private static String stripComments(String text) {
    boolean instr = false;
    boolean escape = false;
    int i;
    outer:
    for (i = 0; i < text.length(); i++) {
      char c = text.charAt(i);
      if (instr) {
        if (escape) { // skip one
          escape = false;
        } else {
          switch (c) {
            case '\\':
              escape = true;
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
  }

  public static String parseName(String token) throws ParseException {
    StringBuilder builder = new StringBuilder(token.length());
    int i, c;
    int hex = 0;
    boolean escape = false;
    for (i = 0; i < token.length(); i += Character.charCount(c)) {
      c = token.codePointAt(i);
      if (hex > 0) {
        if (!((c >= 'A' && c <= 'F') || (c >= 'a' && c <= 'f') || (c >= '0' && c <= '9'))) {
          throw new ParseException("Missing hexadecimal number after unicode escape sequence", i);
        }
        hex--;
        if (hex == 0) {
          c = Integer.parseInt(token.substring(i - 3, i + 1), 16);
          builder.appendCodePoint(c);
        }
      } else if (escape) {
        switch (c) {
          case '"':
            builder.append('"');
            break;
          case '\\':
            builder.append('\\');
            break;
          case '/':
            builder.append('/');
            break;
          case 'b':
            builder.append('\b');
            break;
          case 'f':
            builder.append('\f');
            break;
          case 'n':
            builder.append('\n');
            break;
          case 'r':
            builder.append('\r');
            break;
          case 't':
            builder.append('\t');
            break;
          case 'u':
            hex = 4;
            break;
          default:
            throw new ParseException("Unknown escape sequence", i - 1);
        }
        escape = false;
      } else {
        if (c == '\\') { // test for backslash
          escape = true;
        } else {
          builder.appendCodePoint(c);
        }
      }
    }
    if (escape) throw new ParseException("Incomplete escape sequence", i - 1);
    if (hex > 0) throw new ParseException("Incomplete escape sequence", i - 1);
    return builder.toString();
  }

  private static String parseNameOrError(String token) {
    if (token == null) return null;
    if (token.isEmpty()) return null;
    try {
      return parseName(token);
    } catch (ParseException e) {
      throw new Error(e);
    }
  }

  protected void parseError(AlphParseExceptionType type, Character symbol) throws AlphParseException {
    error = true;
    throw new AlphParseException(lineNo, type, symbol);
  }


  public void parseHeader(String name, String like) throws AlphParseException {
    if (error) throw new IllegalStateException("An error occurred previously");
    if (like != null && !like.matches("^(?:RNA|DNA|PROTEIN)$")) {
      throw new IllegalArgumentException("Argument \"like\" must be RNA, DNA or PROTEIN");
    }
    if (!seenHeader) {
      seenHeader = true;
      this.name = name;
      this.like = (like != null ? Alph.Like.valueOf(like) : null);
    } else {
      parseError(AlphParseExceptionType.HEADER_REDEFINED, null);
    }
  }

  public void parseSymbol(String symStr, String nameStr, String colourStr, String complementStr, String compriseStr, String aliasesStr) throws AlphParseException {
    if (error) throw new IllegalStateException("An error occurred previously");
    // validate values
    if (!SYMBOL_RE.matcher(symStr).matches()) throw new IllegalArgumentException("Symbol invalid");
    if (colourStr != null && !colourStr.isEmpty() && !COLOUR_RE.matcher(colourStr).matches())
      throw new IllegalArgumentException("Colour string must be 6 digit hexadecimal number");
    if (complementStr != null && !complementStr.isEmpty() && !SYMBOL_RE.matcher(complementStr).matches())
      throw new IllegalArgumentException("Complement symbol invalid");
    if (compriseStr != null && !compriseStr.isEmpty() && !SYMBOLS_RE.matcher(compriseStr).matches())
      throw new IllegalArgumentException("Comprise symbols invalid");
    if (aliasesStr != null && !aliasesStr.isEmpty() && !SYMBOLS_RE.matcher(aliasesStr).matches())
      throw new IllegalArgumentException("Aliases symbols invalid");
    // convert the values
    char sym = symStr.charAt(0);
    SortedSet<Character> aliases = new TreeSet<Character>(Alph.SYM_COMPARATOR);
    if (aliasesStr != null) for (int i = 0; i < aliasesStr.length(); i++) aliases.add(aliasesStr.charAt(i));
    String name = (nameStr == null || nameStr.isEmpty() ? null : nameStr.replaceAll("\\p{C}", "")); // remove all control characters
    int colour = (colourStr == null || colourStr.isEmpty() ? 0 : Integer.parseInt(colourStr, 16));
    Character complement = (complementStr == null || complementStr.isEmpty() ? null : complementStr.charAt(0));
    SortedSet<Character> comprise = null;
    if (compriseStr != null && !compriseStr.isEmpty()) {
      comprise = new TreeSet<Character>(Alph.SYM_COMPARATOR);
      for (int i = 0; i < compriseStr.length(); i++) comprise.add(compriseStr.charAt(i));
    }
    seenUpperCase  |= Character.isUpperCase(sym);
    seenLowerCase |= Character.isLowerCase(sym);
    for (char alias : aliases) {
      seenUpperCase  |= Character.isUpperCase(alias);
      seenLowerCase |= Character.isLowerCase(alias);
    }
    if (allSymbols.contains(sym)) {
      parseError(AlphParseExceptionType.REDEFINED_SYMBOL, sym);
      return;
    }
    for (char alias : aliases) {
      if (allSymbols.contains(alias)) {
        parseError(AlphParseExceptionType.REDEFINED_SYMBOL, alias);
        return;
      }
    }
    Symbol symbol = new Symbol(sym, aliases, name, colour, complement, comprise);
    SortedSet<Character> key;
    if (comprise != null) {
      seenAmbiguous = true;
      for (char c : comprise) {
        if (!coreSymbols.contains(c)) {
          parseError(AlphParseExceptionType.UNKNOWN_SYMBOL, c);
          return;
        }
      }
      key = comprise;
    } else {
      if (seenAmbiguous) {
        parseError(AlphParseExceptionType.CORE_AFTER_AMBIG, sym);
        return;
      }
      key = new TreeSet<Character>(Alph.SYM_COMPARATOR);
      key.add(sym);
      key = Collections.unmodifiableSortedSet(key);
    }
    if (sym == '?' && (comprise == null || comprise.size() != coreSymbols.size())) {
      parseError(AlphParseExceptionType.RESERVED_SYMBOL, sym);
      return;
    }
    if (symbolMap.containsKey(key)) {
      List<Symbol> symbols = symbolMap.get(key);
      symbols.add(symbol);
    } else {
      List<Symbol> symbols = new ArrayList<Symbol>(1);
      symbols.add(symbol);
      symbolMap.put(key, symbols);
    }
    allSymbols.add(sym);
    for (char alias : aliases) allSymbols.add(alias);
    if (comprise == null) coreSymbols.add(sym);
  }

  public void parseCore(String symStr, String nameStr, String colourStr, String complementStr, String aliasesStr) throws AlphParseException {
    parseSymbol(symStr, nameStr, colourStr, complementStr, null, aliasesStr);
  }

  public void parseCore(String symStr, String nameStr, String colourStr, String complementStr) throws AlphParseException {
    parseSymbol(symStr, nameStr, colourStr, complementStr, null, null);
  }

  public void parseCore(String symStr, String nameStr, String colourStr) throws AlphParseException {
    parseSymbol(symStr, nameStr, colourStr, null, null, null);
  }
  
  public void parseAmbig(String symStr, String nameStr, String colourStr, String compriseStr, String aliasesStr) throws AlphParseException {
    parseSymbol(symStr, nameStr, colourStr, null, compriseStr, aliasesStr);
  }

  public void parseAmbig(String symStr, String nameStr, String compriseStr, String aliasesStr) throws AlphParseException {
    parseSymbol(symStr, nameStr, null, null, compriseStr, aliasesStr);
  }

  public void parseAmbig(String symStr, String nameStr, String compriseStr) throws AlphParseException {
    parseSymbol(symStr, nameStr, null, null, compriseStr, null);
  }

  public void parseLine(String line) throws AlphParseException {
    if (error) throw new IllegalStateException("An error occurred previously");
    Matcher m;
    lineNo++;
    line = stripComments(line);
    if ("".equals(line)) return;
    if ((m = STARS_RE.matcher(line)).matches()) {
      return;
    } else if ((m = HEADER_RE.matcher(line)).matches()) {
      parseHeader(parseNameOrError(m.group(1)), m.group(2));
    } else if ((m = CORE_PAIR_RE.matcher(line)).matches()) {
      parseSymbol(m.group(1), parseNameOrError(m.group(2)), m.group(3), m.group(4), null, null);
      parseSymbol(m.group(4), parseNameOrError(m.group(5)), m.group(6), m.group(1), null, null);
    } else if ((m = CORE_SINGLE_RE.matcher(line)).matches()) {
      parseSymbol(m.group(1), parseNameOrError(m.group(2)), m.group(3), null, null, null);
    } else if ((m = AMBIG_RE.matcher(line)).matches()) {
      parseSymbol(m.group(1), parseNameOrError(m.group(2)), m.group(3), null, m.group(4), null);
    } else {
      parseError(AlphParseExceptionType.UNKNOWN_PATTERN, null);
    }
  }

  public Alph parseDone() throws AlphParseException {
    if (error) throw new IllegalStateException("An error occurred previously");
    // merge all the aliased symbols into a single symbol
    Map<SortedSet<Character>,Symbol> compriseMap = new TreeMap<SortedSet<Character>, Symbol>(new SymSetComparator());
    for (Map.Entry<SortedSet<Character>,List<Symbol>> entry : symbolMap.entrySet()) {
      // Order the symbols in order of priority
      // This will put the core symbol first, followed by aliases in order letters, numbers, symbols.
      List<Symbol> symbols = entry.getValue();
      Collections.sort(symbols);
      // find values for the name and colour and merge all the aliases
      SortedSet<Character> aliases = new TreeSet<Character>(Alph.SYM_COMPARATOR);
      String name = null;
      Integer colour = null;
      for (Symbol symbol : symbols) {
        if (name == null) name = symbol.name;
        if (colour == null) colour = symbol.colour;
        aliases.add(symbol.sym);
        aliases.addAll(symbol.aliases);
      }
      // when all the symbols are in either lower or upper case then add the other letter case
      if (seenLowerCase != seenUpperCase) {
        SortedSet<Character> aliasesBothCase = new TreeSet<Character>(Alph.SYM_COMPARATOR);
        aliasesBothCase.addAll(aliases);
        if (seenUpperCase) {
          for (char c : aliases) if (c >= 'A' && c <= 'Z') aliasesBothCase.add(Character.toLowerCase(c));
        } else {
          for (char c : aliases) if (c >= 'a' && c <= 'z') aliasesBothCase.add(Character.toUpperCase(c));
        }
        aliases = aliasesBothCase;
      }
      // the first symbol will define the core sym, the complement and for an ambiguous symbol the comprising symbols
      compriseMap.put(entry.getKey(), new Symbol(symbols.get(0), aliases, name, colour));
    }
    // add a wildcard if one is missing
    if (!compriseMap.containsKey(coreSymbols)) {
      compriseMap.put(coreSymbols, new Symbol('?', Collections.singleton('?'), null, null, null, coreSymbols));
    }
    // once we have merged all the aliases then complements need to be extended to the ambiguous symbols
    outer:
    for (Map.Entry<SortedSet<Character>,Symbol> entry : compriseMap.entrySet()) {
      Set<Character> syms = entry.getKey();
      if (syms.size() != 1) {
        SortedSet<Character> complements = new TreeSet<Character>();
        for (char sym : syms) {
          Symbol symbol = compriseMap.get(new TreeSet<Character>(Collections.singleton(sym)));
          if (symbol == null) throw new Error("Could not find referenced core symbol"); // error should have been reported already
          if (symbol.complement == null) continue outer;
          complements.add(symbol.complement);
        }
        Symbol complement = compriseMap.get(complements);
        if (complement != null) {
          entry.setValue(new Symbol(entry.getValue(), complement.sym));
        }
      }
    }
    // confirm that we are deriving the claimed alphabet (if any)
    if (like != null) {
      String core = null;
      String comp = null;
      switch (like) {
        case RNA: core = "ACGU"; break;
        case DNA: core = "ACGT"; comp = "TGCA"; break;
        case PROTEIN: core = "ACDEFGHIKLMNPQRSTVWY"; break;
      }
      if (core != null) {
        for (int i = 0; i < core.length(); i++) {
          Symbol symbol = compriseMap.get(new TreeSet<Character>(Collections.singleton(core.charAt(i))));
          Symbol symbolComplement = symbol != null && symbol.complement != null ?
              compriseMap.get(new TreeSet<Character>(Collections.singleton(symbol.complement))) : null;
          Symbol expectedComplement = comp != null ?
              compriseMap.get(new TreeSet<Character>(Collections.singleton(comp.charAt(i)))) : null;
          if (symbol == null || symbol.comprise != null || ((expectedComplement != null)  && (symbolComplement != expectedComplement))) {
            parseError(AlphParseExceptionType.NOT_LIKE, null);
          }
        }
      }
    }
    // finally we create the alphabet object
    Alph alph = new Alph(this.name, this.like, false, false);
    // Add all the symbols
    for (Symbol symbol : compriseMap.values()) {
      alph.addSymbol(symbol.sym, symbol.aliases, symbol.name, symbol.colour);
    }
    // Add all the complements and comprising
    for (Symbol symbol : compriseMap.values()) {
      if (symbol.complement != null && Alph.symCompare(symbol.sym, symbol.complement) <= 0) {
        alph.setComplementaryPair(symbol.sym, symbol.complement);
      }
      if (symbol.comprise != null) {
        alph.setComprisingSet(symbol.sym, symbol.comprise);
      }
    }
    alph.lock();
    return alph;
  }

  public static Alph parseFile(Path file) throws IOException, AlphParseException {  BufferedReader reader = null;
    try {
      reader = Files.newBufferedReader(file, Charset.forName("UTF-8"));
      AlphParser parser = new AlphParser();
      String line;
      while ((line = reader.readLine()) != null) {
        parser.parseLine(line);
      }
      Alph alph = parser.parseDone();
      reader.close(); reader = null;
      return alph;
    } finally {
      if (reader != null) {
        try {
          reader.close();
        } catch (IOException e) {/* ignore */}
      }
    }
  }

  public static Alph standardDNA() {
    try {
      AlphParser parser = new AlphParser();
      parser.parseHeader("DNA", "DNA");
      // core symbols
      parser.parseCore("A", "Adenine", "CC0000", "T");
      parser.parseCore("C", "Cytosine", "0000CC", "G");
      parser.parseCore("G", "Guanine", "FFB300", "C");
      parser.parseCore("T", "Thymine", "008000", "A", "U");
      // ambiguous symbols
      parser.parseAmbig("W", "Weak", "AT");
      parser.parseAmbig("S", "Strong", "CG");
      parser.parseAmbig("M", "Amino", "AC");
      parser.parseAmbig("K", "Keto", "GT");
      parser.parseAmbig("R", "Purine", "AG");
      parser.parseAmbig("Y", "Pyrimidine", "CT");
      parser.parseAmbig("B", "Not A", "CGT");
      parser.parseAmbig("D", "Not C", "AGT");
      parser.parseAmbig("H", "Not G", "ACT");
      parser.parseAmbig("V", "Not T", "ACG");
      parser.parseAmbig("N", "Any base", "ACGT", "X.");
      // process
      return parser.parseDone();
    } catch (AlphParseException e) {
      throw new Error("Parsing errors should not happen for correct input!", e);
    }
  }

  public static Alph standardRNA() {
    try {
      AlphParser parser = new AlphParser();
      parser.parseHeader("RNA", "RNA");
      // core symbols
      parser.parseCore("A", "Adenine", "CC0000");
      parser.parseCore("C", "Cytosine", "0000CC");
      parser.parseCore("G", "Guanine", "FFB300");
      parser.parseCore("U", "Uracil", "008000", null, "T");
      // ambiguous symbols
      parser.parseAmbig("W", "Weak", "AU");
      parser.parseAmbig("S", "Strong", "CG");
      parser.parseAmbig("M", "Amino", "AC");
      parser.parseAmbig("K", "Keto", "GU");
      parser.parseAmbig("R", "Purine", "AG");
      parser.parseAmbig("Y", "Pyrimidine", "CU");
      parser.parseAmbig("B", "Not A", "CGU");
      parser.parseAmbig("D", "Not C", "AGU");
      parser.parseAmbig("H", "Not G", "ACU");
      parser.parseAmbig("V", "Not U", "ACG");
      parser.parseAmbig("N", "Any base", "ACGU", "X.");
      // process
      return parser.parseDone();
    } catch (AlphParseException e) {
      throw new Error("Parsing errors should not happen for correct input!", e);
    }
  }

  public static Alph standardProtein() {
    try {
      AlphParser parser = new AlphParser();
      parser.parseHeader("Protein", "PROTEIN");
      // core symbols
      parser.parseCore("A", "Alanine", "0000CC");
      parser.parseCore("R", "Arginine", "CC0000");
      parser.parseCore("N", "Asparagine", "008000");
      parser.parseCore("D", "Aspartic acid", "FF00FF");
      parser.parseCore("C", "Cysteine", "0000CC");
      parser.parseCore("E", "Glutamic acid", "FF00FF");
      parser.parseCore("Q", "Glutamine", "008000");
      parser.parseCore("G", "Glycine", "FFB300");
      parser.parseCore("H", "Histidine", "FFCCCC");
      parser.parseCore("I", "Isoleucine", "0000CC");
      parser.parseCore("L", "Leucine", "0000CC");
      parser.parseCore("K", "Lysine", "CC0000");
      parser.parseCore("M", "Methionine", "0000CC");
      parser.parseCore("F", "Phenylalanine", "0000CC");
      parser.parseCore("P", "Proline", "FFFF00");
      parser.parseCore("S", "Serine", "008000");
      parser.parseCore("T", "Threonine", "008000");
      parser.parseCore("W", "Tryptophan", "0000CC");
      parser.parseCore("Y", "Tyrosine", "33E6CC");
      parser.parseCore("V", "Valine", "0000CC");
      // ambiguous symbols
      parser.parseAmbig("B", "Asparagine or Aspartic acid", "ND");
      parser.parseAmbig("Z", "Glutamine or Glutamic acid", "QE");
      parser.parseAmbig("J", "Leucine or Isoleucine", "LI");
      parser.parseAmbig("X", "Any amino acid", "ARNDCEQGHILKMFPSTWYV", "*.");
      // process
      return parser.parseDone();
    } catch (AlphParseException e) {
      throw new Error("Parsing errors should not happen for correct input!", e);
    }
  }

  private static final class Symbol implements Comparable<Symbol> {
    public final char sym;
    public final SortedSet<Character> aliases;
    public final String name;
    public final Integer colour;
    public final Character complement;
    public final SortedSet<Character> comprise;

    public Symbol(char sym, Set<Character> aliases, String name, Integer colour, Character complement, Set<Character> comprise) {
      this.sym = sym;
      this.aliases = Collections.unmodifiableSortedSet(new TreeSet<Character>(aliases));
      this.name = name;
      this.colour = colour;
      this.complement = complement;
      this.comprise = (comprise != null ? Collections.unmodifiableSortedSet(new TreeSet<Character>(comprise)) : null);
    }

    public Symbol(Symbol base, Set<Character> aliases, String name, Integer colour) {
      this.sym = base.sym;
      this.aliases = Collections.unmodifiableSortedSet(new TreeSet<Character>(aliases));
      this.name = name;
      this.colour = colour;
      this.complement = base.complement;
      this.comprise = base.comprise;
    }

    public Symbol(Symbol base, char complement) {
      this.sym = base.sym;
      this.aliases = base.aliases;
      this.name = base.name;
      this.colour = base.colour;
      this.complement = complement;
      this.comprise = base.comprise;
    }

    public Set<Character> compriseSet() {
      if (comprise != null) {
        return comprise;
      } else {
        return Collections.singleton(sym);
      }
    }

    @Override
    public int compareTo(Symbol other) {
      if (other == null) return -1;
      // compare comprise
      if (this.comprise != null && other.comprise != null) {
        if (this.comprise.size() != other.comprise.size()) return other.comprise.size() - this.comprise.size();
        Iterator<Character> localComprise = this.comprise.iterator();
        Iterator<Character> otherComprise = other.comprise.iterator();
        while (localComprise.hasNext()) {
          int cmp = Alph.symCompare(localComprise.next(), otherComprise.next());
          if (cmp != 0) return cmp;
        }
      } else if (this.comprise != null) {
        return 1;
      } else if (other.comprise != null) {
        return -1;
      }
      // compare sym
      {
        int cmp = Alph.symCompare(this.sym, other.sym);
        if (cmp != 0) return cmp;
      }
      // compare complement
      if (this.complement != null && other.complement != null) {
        int cmp = Alph.symCompare(this.complement, other.complement);
        if (cmp != 0) return cmp;
      } else if (this.complement != null) {
        return 1;
      } else if (other.complement != null) {
        return -1;
      }
      // compare aliases
      {
        if (this.aliases.size() != other.aliases.size()) return other.aliases.size() - this.aliases.size();
        Iterator<Character> localAliases = this.aliases.iterator();
        Iterator<Character> otherAliases = other.aliases.iterator();
        while (localAliases.hasNext()) {
          int cmp = Alph.symCompare(localAliases.next(), otherAliases.next());
          if (cmp != 0) return cmp;
        }
      }
      // compare name
      if (this.name != null && other.name != null) {
        int cmp = this.name.compareTo(other.name);
        if (cmp != 0) return cmp;
      } else if (this.name != null) {
        return 1;
      } else if (other.name != null) {
        return -1;
      }
      // compare colour
      if (this.colour != null && other.colour != null) {
        return this.colour - other.colour;
      } else if (this.colour != null) {
        return 1;
      } else if (other.colour != null) {
        return -1;
      }
      // everything equal
      return 0;
    }
  }

  public static class AlphParseException extends Exception {
    public final int line;
    public final AlphParseExceptionType type;
    public final Character symbol;
    public AlphParseException(int line, AlphParseExceptionType type, Character symbol) {
      this.line = line;
      this.type = type;
      this.symbol = symbol;
    }
    public String toString() {
      return "Alphabet parse exception \"" + type + "\" on line " + line + " due to symbol " + symbol;
    }
  }

  public enum AlphParseExceptionType {
    UNKNOWN_PATTERN("the text does not match a known pattern"),
    HEADER_REDEFINED("the header has been listed more than once"),
    CORE_AFTER_AMBIG("the core symbol must be listed before any ambiguous symbols"),
    UNKNOWN_SYMBOL("the core symbol is unknown, all core symbols must be listed before ambiguous symbols"),
    REDEFINED_SYMBOL("the symbol has already been defined"),
    RESERVED_SYMBOL("the question mark '?' has been reserved for use as a wildcard only"),
    NOT_LIKE("the alphabet is not like the suggested standard alphabet");

    final String description;

    private AlphParseExceptionType(String description) {
      this.description = description;
    }

    public String toString() {
      return description;
    }
  }

  public class SymSetComparator implements Comparator<SortedSet<Character>> {
    @Override
    public int compare(SortedSet<Character> syms1, SortedSet<Character> syms2) {
      if (syms1.size() == syms2.size()) {
        Iterator<Character> iter1, iter2;
        iter1 = syms1.iterator();
        iter2 = syms2.iterator();
        while (iter1.hasNext()) {
          int cmp = Alph.symCompare(iter1.next(), iter2.next());
          if (cmp != 0) return cmp;
        }
        return 0;
      } else if (syms1.size() < syms2.size()) {
        return -1;
      } else {
        return 1;
      }
    }
  }
}
