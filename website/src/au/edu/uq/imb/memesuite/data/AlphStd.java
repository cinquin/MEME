package au.edu.uq.imb.memesuite.data;

import au.edu.uq.imb.memesuite.io.alph.AlphParser;

import java.util.*;

/**
 * An enumeration of standard alphabets
 */
public enum AlphStd {
  RNA(AlphParser.standardRNA()), DNA(AlphParser.standardDNA()), PROTEIN(AlphParser.standardProtein());

  /**
   * The alphabet representation of this standard alphabet.
   */
  private final Alph alph;

  /**
   * Create a standard alphabet.
   * @param alph the actual representation of the alphabet.
   */
  private AlphStd(Alph alph) {
    this.alph = alph;
  }

  /**
   * Get the name of the alphabet.
   * @return the name of the alphabet.
   */
  public String toString() {
    return this.alph.getName();
  }

  /**
   * Get the actual alphabet.
   * @return the actual alphabet.
   */
  public Alph getAlph() {
    return alph;
  }

  /**
   * Get the bitset representation of this alphabet.
   * @return the bitset representation of this alphabet.
   */
  public int toInt() {
    return (1 << this.ordinal());
  }

  /**
   * Guesses the alphabet from an array of counts of ASCII characters.
   * The array is not modified.
   * @param counts counts of all 128 ASCII letters and a 129th entry to count all non-ascii.
   * @return the closest matching standard alphabet.
   */
  public static AlphStd guess(long[] counts) {
    if (counts.length != 129) throw new IllegalArgumentException("Counts array must have exactly 129 entries");
    Map<Character,Long> charCounts = new TreeMap<Character, Long>();
    for (int i = 0; i < 128; i++) {
      if (counts[i] > 0) {
        charCounts.put((char)i, counts[i]);
      }
    }
    return guess(charCounts);
  }

  /**
   * Guess the alphabet from a map of ASCII characters to counts.
   * @param counts the characters and their counts.
   * @return the closest matching standard alphabet.
   */
  public static AlphStd guess(Map<Character,Long> counts) {
    // create a list of alphabets from the enumeration
    List<Alph> alphabets = new ArrayList<Alph>(AlphStd.values().length);
    for (int i = 0; i < AlphStd.values().length; i++) alphabets.add(AlphStd.values()[i].alph);
    // find the best match
    int bestIndex = Alph.guess(alphabets, counts);
    // return the best matching enum
    return AlphStd.values()[bestIndex];
  }

  /**
   * Convert from a bit representation to the Enum value.
   * @param alphaBit an integer with a single bit set.
   * @return the standard alphabet represented by that bit.
   */
  public static AlphStd fromInt(int alphaBit) {
    int pos = Integer.numberOfTrailingZeros(alphaBit);
    return AlphStd.values()[pos];
  }

  /**
   * Attempt to find a matching alphabet enum for the alphabet.
   * @param alph the alphabet.
   * @return the matching alphabet enum or null on failure.
   */
  public static AlphStd fromAlph(Alph alph) {
    Alph.Like like = alph.getLikeAlphabet();
    if (like == null) return null;
    // this assumes we name everything in the Like enum the same as in the AlphStd enum
    AlphStd alphStd = AlphStd.valueOf(like.name());
    assert(alphStd != null);
    if (alphStd.alph.equals(alph)) return alphStd;
    return null;
  }

  /**
   * Attempt to find a matching alphabet enum for the alphabet core symbols.
   * @param core the alphabet core symbols in order.
   * @return the matching alphabet enum or null on failure.
   */
  public static AlphStd fromCore(String core) {
    outer:
    for (AlphStd alphStd : AlphStd.values()) {
      Alph alph = alphStd.getAlph();
      if (alph.sizeCore() != core.length()) continue;
      for (int i = 0; i < core.length(); i++) {
        if (!alph.index(i).isPrimeSym(core.charAt(i))) continue outer;
      }
      return alphStd;
    }
    return null;
  }
}
