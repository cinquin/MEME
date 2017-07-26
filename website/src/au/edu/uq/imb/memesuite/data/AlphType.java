package au.edu.uq.imb.memesuite.data;

import java.util.Collections;
import java.util.EnumSet;
import java.util.Set;

/**
 * The different types of alphabet that can be accepted by a input.
 */
public enum AlphType {
  ANY_ALPHABET(EnumSet.allOf(AlphStd.class)),
  COMPLEMENTABLE_ALPHABET(EnumSet.of(AlphStd.DNA)),
  STANDARD_ALPHABET(EnumSet.allOf(AlphStd.class)),
  RNA_OR_DNA_ALPHBET(EnumSet.of(AlphStd.RNA, AlphStd.DNA)),
  DNA_OR_PROTEIN_ALPHABET(EnumSet.of(AlphStd.DNA, AlphStd.PROTEIN)),
  RNA_ALPHABET(EnumSet.of(AlphStd.RNA)),
  DNA_ALPHABET(EnumSet.of(AlphStd.DNA)),
  PROTEIN_ALPHABET(EnumSet.of(AlphStd.PROTEIN)),
  ;

  /** The unmodifiable set of standard alphabets that this type matches. */
  private Set<AlphStd> allowedStandardAlphabets;

  /**
   * Construct a alphabet type taking the set of standard alphabets that this type matches.
   * @param alphStds the set of standard alphabets.
   */
  private AlphType(EnumSet<AlphStd> alphStds) {
    allowedStandardAlphabets = Collections.unmodifiableSet(alphStds);
  }

  /**
   * Get the set of standard alphabets that are permissible for this alphabet type.
   * @return the set of permissible standard alphabets.
   */
  public Set<AlphStd> getStandardAlphabets() {
    return allowedStandardAlphabets;
  }

  /**
   * Check if the given alphabet is covered by this alphabet type.
   * @param alph the alphabet.
   * @return true if the alphabet is covered by this alphabet type.
   */
  public boolean matches(Alph alph) {
    return matches(this, alph);
  }

  /**
   * Check if the given standard alphabet is covered by this alphabet type.
   * @param alphStd the standard alphabet.
   * @return true if the standard alphabet is covered by this alphabet type.
   */
  public boolean matches(AlphStd alphStd) {
    return matches(this, alphStd);
  }

  /**
   * Check if an alphabet matches the definition given by the type.
   * @param type a type of alphabet.
   * @param alph a specific alphabet.
   * @return true if the alphabet meets the definition given by the type.
   */
  public static boolean matches(AlphType type, Alph alph) {
    switch (type) {
      case ANY_ALPHABET:
        return true;
      case COMPLEMENTABLE_ALPHABET:
        return alph.isComplementable();
      default:
        for (AlphStd alphStd : type.getStandardAlphabets()) {
          if (alphStd.getAlph().equals(alph)) return true;
        }
        return false;
    }
  }

  /**
   * Check if a standard alphabet matches the definition given by the type.
   * @param type a type of alphabet.
   * @param alphStd a specific standard alphabet.
   * @return true if the alphabet meets the definition given by the type.
   */
  public static boolean matches(AlphType type, AlphStd alphStd) {
    switch (type) {
      case ANY_ALPHABET:
        return true;
      case COMPLEMENTABLE_ALPHABET:
        return alphStd.getAlph().isComplementable();
      case STANDARD_ALPHABET:
        return true;
      default:
        return type.getStandardAlphabets().contains(alphStd);
    }
  }
}
