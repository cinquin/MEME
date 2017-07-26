package au.edu.uq.imb.memesuite.data;

import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.IOException;
import java.util.*;

/**
 * An immutable alphabet representation
 */
public final class Alph implements JsonWr.JsonValue, Comparable<Alph> {
  public static final Comparator<Character> SYM_COMPARATOR = new SymComparator();
  private boolean locked;
  private Boolean caseInsensitive;
  private List<Symbol> symbols;
  private Map<Character,Symbol> lookup;
  private String name;
  private Like like;
  private boolean seenUppercase;
  private boolean seenLowercase;
  private int ncore;
  private boolean complementable;

  /**
   * Used for construction of the Alphabet object by the AlphParser!
   * Construct an empty unlocked Alphabet object.
   * @param name the name of the alphabet or null if not given.
   * @param like the base alphabet or null if no base.
   * @param seenUppercase does the alphabet use upper case?
   * @param seenLowercase does the alphabet use lower case?
   */
  public Alph(String name, Like like, boolean seenUppercase, boolean seenLowercase) {
    this.locked = false;
    this.symbols = new ArrayList<Symbol>();
    this.lookup = new TreeMap<Character, Symbol>();
    this.name = name;
    this.like = like;
    this.seenUppercase = seenUppercase;
    this.seenLowercase = seenLowercase;
    this.ncore = 0;
  }

  /**
   * Used for construction of the Alphabet object by the AlphParser!
   * This adds a symbol to the alphabet.
   * The Symbols MUST be added in the order needed by the alphabet.
   * @param sym the primary character used to refer to this Symbol.
   * @param aliases all characters that can be used to refer to this Symbol including the primary character.
   * @param name the name of the symbol or null if not given.
   * @param colour the colour of the symbol in 24bit RGB or null if not given.
   */
  public void addSymbol(char sym, Set<Character> aliases, String name, Integer colour) {
    assertUnlocked();
    Symbol symbol = new Symbol(symbols.size(), sym, aliases, name, colour);
    symbols.add(symbol);
    for (char alias : aliases) {
      lookup.put(alias, symbol);
    }
  }

  /**
   * Used for construction of the Alphabet object by the AlphParser!
   * This pairs two symbols as complements of each other. The symbols may be the same.
   * The symbols must have already been added by a call to addSymbol before calling.
   * This method does not check the values given beyond that they exist as Symbol objects.
   * @param sym1 the first symbol.
   * @param sym2 the second symbol.
   */
  public void setComplementaryPair(char sym1, char sym2) {
    assertUnlocked();
    Symbol symbol1 = lookupAndCheck(sym1);
    Symbol symbol2 = lookupAndCheck(sym2);
    symbol1.complement = symbol2;
    symbol2.complement = symbol1;
  }

  /**
   * Used for construction of the Alphabet object by the AlphParser!
   * This sets the comprising symbols of an ambiguous symbol.
   * @param sym the ambiguous symbol.
   * @param comprisingSyms the comprising symbols of the ambiguous symbol.
   */
  public void setComprisingSet(char sym, Set<Character> comprisingSyms) {
    assertUnlocked();
    Symbol symbol = lookupAndCheck(sym);
    TreeSet<Symbol> comprisingSymbols = new TreeSet<Symbol>();
    for (char comprisingSym : comprisingSyms) {
      comprisingSymbols.add(lookupAndCheck(comprisingSym));
    }
    symbol.comprise = Collections.unmodifiableSortedSet(comprisingSymbols);
  }

  /**
   * Used for construction of the Alphabet object by the AlphParser!
   * This locks out future changes so the object can be passed around as if immutable.
   */
  public void lock() {
    assertUnlocked();
    int i;
    this.complementable = true;
    for (i = 0; i < symbols.size(); i++) {
      Symbol symbol = symbols.get(i);
      if (symbol.comprise != null) {
        ncore = i;
        break;
      } else if (symbol.complement == null) {
        complementable = false;
      }
    }
    this.locked = true;
  }

  /**
   * Returns the number of core symbols in the alphabet.
   * @return the size of the alphabet core.
   */
  public int sizeCore() {
    assert(locked);
    return ncore;
  }

  public boolean isCaseInsensitive() {
    assert(locked);
    // cached result
    if (caseInsensitive != null) return caseInsensitive;
    // assume true and disprove
    caseInsensitive = true;
    for (Symbol symbol : symbols) {
      if (symbol.sym >= 'A' || symbol.sym <= 'Z') {
        if (lookup.get(Character.toLowerCase(symbol.sym)) != symbol) {
          caseInsensitive = false;
          break;
        }
      } else if (symbol.sym >= 'a' || symbol.sym <= 'z') {
        if (lookup.get(Character.toUpperCase(symbol.sym)) != symbol) {
          caseInsensitive = false;
          break;
        }
      }
      for (Character alias : symbol.aliases) {

        if (alias >= 'A' || alias <= 'Z') {
          if (lookup.get(Character.toLowerCase(alias)) != symbol) {
            caseInsensitive = false;
            break;
          }
        } else if (alias >= 'a' || alias <= 'z') {
          if (lookup.get(Character.toUpperCase(alias)) != symbol) {
            caseInsensitive = false;
            break;
          }
        }
      }
    }
    return caseInsensitive;
  }

  public boolean isComplementable() {
    assert(locked);
    return complementable;
  }

  /**
   * Returns the name of the alphabet.
   * @return the alphabet name.
   */
  public String getName() {
    assert(locked);
    return name;
  }

  /**
   * Returns the classification of the alphabet.
   * @return the alphabet name.
   */
  public Like getLikeAlphabet() {
    assert(locked);
    return like;
  }

  /**
   * Returns the name of the alphabet or if a name is not given then the string of core symbols.
   * @return a user readable representation of the alphabet.
   */
  public String toString() {
    if (locked) {
      if (name != null) {
        return name;
      } else {
        StringBuilder builder = new StringBuilder();
        for (Symbol symbol : symbols) {
          if (symbol.isCore()) builder.append(symbol.sym);
        }
        return builder.toString();
      }
    } else {
      return super.toString() + " (construction incomplete!)";
    }
  }

  /**
   * Returns the symbol with the given index.
   * @return the symbol object.
   */
  public Symbol index(int idx) {
    assert(locked);
    return symbols.get(idx);
  }

  /**
   * Lookup the symbol with the given sym.
   * @param sym the character used to refer to a symbol.
   * @return the symbol object.
   */
  public Symbol find(char sym) {
    assert(locked);
    return lookup.get(sym);
  }

  /**
   * Lookup the symbol with the given sym.
   * @param sym the one character string used to refer to a symbol.
   * @return the symbol object.
   */
  public Symbol find(String sym) {
    assert(locked);
    if (sym.length() != 1) return null;
    return lookup.get(sym.charAt(0));
  }

  /**
   * Check if the symbol exists in the alphabet.
   * @param sym the symbol to check.
   * @return true if the symbol exists.
   */
  public boolean isValidSym(char sym) {
    assert(locked);
    return lookup.containsKey(sym);
  }

  /**
   * Check if a set of symbols exists in the alphabet.
   * @param syms the symbols to check.
   * @return true if all the symbols exist.
   */
  public boolean isValidSymSet(Set<Character> syms) {
    assert(locked);
    for (char sym : syms) if (!isValidSym(sym)) return false;
    return true;
  }

  /**
   * Check if counts of ASCII symbols would be accepted by this alphabet.
   * Same as calling with the set of characters that have a non-zero count.
   * @param counts the counts of all ASCII symbols (indexes 0-127) and a extra count for all the non-ASCII symbols (index 128).
   * @return true if all the ASCII symbols with non-zero counts would be accepted and no non-ASCII symbols were encountered.
   */
  public boolean isValidSymSet(long[] counts) {
    if (counts.length != 129) throw new IllegalArgumentException("Counts array must have exactly 129 entries");
    if (counts[128] != 0) return false;
    for (int i = 0; i < 128; i++) {
      if (counts[i] > 0 && !isValidSym((char) i)) return false;
    }
    return true;
  }

  /**
   * Check if this Alph object is unlocked. Used for ensuring the Alph object can only be changed before being locked.
   * @throws IllegalStateException if this Alph object is locked.
   */
  private void assertUnlocked() throws IllegalStateException {
    if (locked) throw new IllegalStateException("Alphabet has been locked. No more changes are allowed.");
  }

  /**
   * Check if the symbol exists and return the related object if it does.
   * @param sym the symbol to lookup.
   * @return the object that stores information about the sym.
   * @throws IllegalStateException if the symbol is unknown.
   */
  private Symbol lookupAndCheck(char sym) throws IllegalStateException {
    Symbol symbol = lookup.get(sym);
    if (symbol == null) throw new IllegalStateException("Symbol " + sym + " is unknown!");
    return symbol;
  }

  public static int checkCoreSubset(Alph subAlph, Alph superAlph) {
    boolean complementSame = true;
    BitSet bitset = new BitSet(128);
    for (int i = 0; i < subAlph.sizeCore(); i++) {
      Symbol subSymbol = subAlph.index(i);
      Symbol superSymbol = superAlph.find(subSymbol.sym);
      if (superSymbol == null) return 0; // missing symbol so not a superset!
      assert(superSymbol.index < 128); // there's only 127 ASCII letters
      if (bitset.get(superSymbol.index)) return 0;
      bitset.set(superSymbol.index);
      // check complement
      if (subSymbol.complement != null && superSymbol.complement != null) {
        if (superAlph.find(subSymbol.complement.sym) != superSymbol.complement) {
          complementSame = false;
        }
      } else if (subSymbol.complement != null || superSymbol.complement != null) {
        complementSame = false;
      }
    }
    return (complementSame ? 1 : -1);
  }

  /**
   * Return the index of the best match to the alphabet by a frequency analysis of a map of ASCII characters to counts.
   * @param counts the characters and their counts.
   * @return the closest matching standard alphabet.
   */
  public static int guess(List<Alph> options, Map<Character,Long> counts) {
    List<Integer> validAlphabets = new ArrayList<Integer>(options.size());
    for (int i = 0; i < options.size(); i++) validAlphabets.add(i);
    // filter out alphabets that don't support all the symbols
    Iterator<Integer> it = validAlphabets.iterator();
    while (it.hasNext()) if (!options.get(it.next()).isValidSymSet(counts.keySet())) it.remove();
    // if that removed all the alphabets then reset and try another approach
    if (validAlphabets.isEmpty()) for (int i = 0; i < options.size(); i++) validAlphabets.add(i);
    // when there is only one alphabet left that is our guess
    if (validAlphabets.size() == 1) return validAlphabets.get(0);
    // otherwise try to rank the alphabets
    double bestScore = 0;
    double bestPrimeScore = 0;
    List<Integer> bestAlphabets = new ArrayList<Integer>();
    for (int index : validAlphabets) {
      Alph alph = options.get(index);
      Set<Integer> seen = new TreeSet<Integer>();
      Set<Integer> primeSeen = new TreeSet<Integer>();
      long primeCount = 1; // set count to 1 to avoid divide by zero
      long altCount = 0;
      long otherCount = 0;
      long unknownCount = 0;
      for (Map.Entry<Character,Long> entry : counts.entrySet()) {
        char sym = entry.getKey();
        long count = entry.getValue();
        Alph.Symbol symbol = alph.find(sym);
        if (symbol != null) {
          if (symbol.isCore()) {
            seen.add(symbol.index);
            if (symbol.isPrimeSym(sym)) {
              primeSeen.add(symbol.index);
              primeCount += count;
            } else {
              altCount += count;
            }
          } else if (symbol.isWildcard()) {
            if (symbol.isPrimeSym(sym)) {
              primeCount += count;
            } else {
              altCount += count;
            }
          } else {
            otherCount += count;
          }
        } else {
          unknownCount += count;
        }
      }
      long totalCount = primeCount + altCount + otherCount + unknownCount;
      double score = ((double)(primeCount + altCount) / totalCount) * ((double)seen.size() / alph.sizeCore());
      double primeScore = ((double)primeCount / totalCount) * ((double)primeSeen.size() / alph.sizeCore());
      if (score > bestScore || (score == bestScore && primeScore > bestPrimeScore)) {
        bestAlphabets.clear();
        bestAlphabets.add(index);
        bestScore = score;
        bestPrimeScore = primeScore;
      } else if (score == bestScore && primeScore == bestPrimeScore) {
        bestAlphabets.add(index);
      }
    }
    return bestAlphabets.get(0);
  }

  public static int symCompare(char sym1, char sym2) {
    if (Character.isLetter(sym1)) {
      if (Character.isLetter(sym2)) {
        return sym1 - sym2;
      } else {
        return -1;
      }
    } else if (Character.isLetter(sym2)) {
      return 1;
    }
    if (Character.isDigit(sym1)) {
      if (Character.isDigit(sym2)) {
        return sym1 - sym2;
      } else {
        return -1;
      }
    } else if (Character.isDigit(sym2)) {
      return 1;
    }
    return sym1 - sym2;
  }

  @Override
  public int compareTo(Alph other) {
    if (this == other) return 0;
    if (other == null) return -1;
    if (this.like == other.like) {
      if (symbols.size() != other.symbols.size()) return symbols.size() - other.symbols.size();
      for (int i = 0; i < symbols.size(); i++) {
        int cmp = symbols.get(i).compareTo(other.symbols.get(i));
        if (cmp != 0) return cmp;
      }
      return name.compareTo(other.name);
    } else if (this.like == null) {
      return 1; // put new alphabets after ones based on existing alphabets
    } else if (other.like == null) {
      return -1; // put new alphabets after ones based on existing alphabets
    } else {
      return (this.like.ordinal() < other.like.ordinal() ? -1 : 1);
    }
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;
    Alph other = (Alph) o;
    return like == other.like && !(name != null ? !name.equals(other.name) : other.name != null) && symbols.equals(other.symbols);
  }

  @Override
  public int hashCode() {
    int result = symbols.hashCode();
    result = 31 * result + (name != null ? name.hashCode() : 0);
    result = 31 * result + (like != null ? like.hashCode() : 0);
    return result;
  }

  @Override
  public void outputJson(JsonWr out) throws IOException {
    out.startObject();
    if (name != null) out.property("name", name);
    if (like != null) out.property("like", like.name().toLowerCase());
    out.property("ncore", ncore);
    out.property("symbols");
    out.startArray();
    for (Symbol symbol : symbols) out.value(symbol);
    out.endArray();
    out.endObject();
  }

  public enum Like { RNA, DNA, PROTEIN }

  private static class SymComparator implements Comparator<Character> {
    @Override
    public int compare(Character c1, Character c2) {
      if (c1 != null && c2 != null) {
        return symCompare(c1, c2);
      } else if (c1 != null) {
        return -1;
      } else if (c2 != null) {
        return 1;
      } else {
        return 0;
      }
    }
  }

  public final class Symbol implements Comparable<Symbol>, JsonWr.JsonValue {
    public final int index;
    public final char sym;
    public final SortedSet<Character> aliases;
    public final String name;
    public final Integer colour;
    private Symbol complement;
    private SortedSet<Symbol> comprise;

    private Symbol(int index, char sym, Set<Character> aliases, String name, Integer colour) {
      this.index = index;
      this.sym = sym;
      this.aliases = Collections.unmodifiableSortedSet(new TreeSet<Character>(aliases));
      this.name = name;
      this.colour = colour;
      this.complement = null;
      this.comprise = null;
    }

    public boolean isCore() {
      assert(locked);
      return comprise == null;
    }

    public boolean isWildcard() {
      assert(locked);
      return comprise != null && comprise.size() == ncore;
    }

    public Symbol getComplement() {
      assert(locked);
      return complement;
    }

    public Set<Symbol> getComprising() {
      assert(locked);
      if (comprise != null) return comprise;
      return Collections.singleton(this);
    }

    public boolean isPrimeSym(char sym) {
      assert(locked);
      return isCaseInsensitive() ? Character.toLowerCase(this.sym) == Character.toLowerCase(sym) : this.sym == sym;
    }

    @Override
    public boolean equals(Object obj) {
      return obj instanceof Symbol && compareTo((Symbol) obj) == 0;
    }

    protected int nonRecursiveHashCode() {
      int result = (int) sym;
      result = 31 * result + aliases.hashCode();
      result = 31 * result + (name != null ? name.hashCode() : 0);
      result = 31 * result + (colour != null ? colour.hashCode() : 0);
      result = 31 * result + (comprise != null ? comprise.hashCode() : 0);
      return result;

    }

    @Override
    public int hashCode() {
      int result = nonRecursiveHashCode();
      result = 31 * result + (complement != null ? complement.nonRecursiveHashCode() : 0);
      return result;
    }

    @Override
    public int compareTo(Symbol other) {
      if (this == other) return 0;
      if (other == null) return -1;
      // compare comprise
      if (this.comprise != null && other.comprise != null) {
        if (this.comprise.size() != other.comprise.size()) return other.comprise.size() - this.comprise.size();
        Iterator<Symbol> localComprise = this.comprise.iterator();
        Iterator<Symbol> otherComprise = other.comprise.iterator();
        while (localComprise.hasNext()) {
          int cmp = symCompare(localComprise.next().sym, otherComprise.next().sym);
          if (cmp != 0) return cmp;
        }
      } else if (this.comprise != null) {
        return 1; // other must be a core symbol
      } else if (other.comprise != null) {
        return -1; // this must be a core symbol
      }
      // compare sym
      {
        int cmp = symCompare(this.sym, other.sym);
        if (cmp != 0) return cmp;
      }
      // compare complement
      if (this.complement != null && other.complement != null) {
        int cmp = symCompare(this.complement.sym, other.complement.sym);
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

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      out.property("symbol", Character.toString(sym));
      {
        StringBuilder syms = new StringBuilder();
        for (char asym : aliases) {
          if (asym == sym) continue;
          if (seenLowercase != seenUppercase) {
            if (asym >= 'A' && asym <= 'Z') {
              if (seenLowercase) continue;
            } else if (asym >= 'a' && asym <= 'z') {
              if (seenUppercase) continue;
            }
          }
          syms.append(asym);
        }
        if (syms.length() > 0) out.property("aliases", syms.toString());
      }
      if (name != null) out.property("name", name);
      if (colour != null) out.property("colour", String.format("%06X", colour));
      if (isCore()) {
        if (getComplement() != null) out.property("complement", Character.toString(getComplement().sym));
      } else {
        StringBuilder syms = new StringBuilder(comprise.size());
        for (Symbol symbol : comprise) syms.append(symbol.sym);
        out.property("equals", syms);
      }
      out.endObject();
    }
  }
}
