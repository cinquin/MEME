package au.edu.uq.imb.memesuite.data;

import au.edu.uq.imb.memesuite.util.JsonWr;
import au.edu.uq.imb.memesuite.util.SampleStats;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;

/**
 * Calculate statistics on the sequences
 */
public final class SequenceStats implements SequenceInfo {
  private boolean frozen;
  private SampleStats stats;
  private long seqLength;
  private long[] asciiCounts;
  private AlphStd alphabetGuess;

  /**
   * Class constructor
   */
  public SequenceStats() {
    this.stats = new SampleStats(true);
    this.seqLength = 0;
    this.asciiCounts = new long[129];
    Arrays.fill(this.asciiCounts, 0);
  }

  /**
   * Process the sequence and update the statistics
   * @param sequence the complete sequence
   * @see #addSeqPart(CharSequence)
   * @see #endSeq()
   */
  public void addSeq(CharSequence sequence) {
    this.addSeqPart(sequence);
    this.endSeq();
  }

  /**
   * Process the sequence and update the statistics
   * @param data a raw UTF-8 encoded buffer containing the complete sequence
   * @param start the start index of the sequence in the buffer
   * @param length the length of the sequence in the buffer
   * @see #addSeqPart(byte[], int, int)
   * @see #endSeq()
   */
  public void addSeq(byte[] data, int start, int length) {
    this.addSeqPart(data, start, length);
    this.endSeq();
  }

  /**
   * Add part of a sequence to be processed.
   * When the full sequence has been process with this method endSeq should be called.
   * @param sequencePart The characters that comprise part of a sequence
   * @see #endSeq()
   */
  public void addSeqPart(CharSequence sequencePart) {
    if (this.frozen) throw new IllegalStateException("SequenceStats is frozen.");
    this.seqLength += sequencePart.length();
    int codepoint;
    for (int i = 0; i < sequencePart.length(); i += Character.charCount(codepoint)) {
      codepoint = Character.codePointAt(sequencePart, i);
      if (codepoint >= 0 && codepoint <= 127) {
        this.asciiCounts[codepoint]++;
      } else {
        this.asciiCounts[128]++;
      }
    }
  }

  /**
   * Add part of a sequence to be processed.
   * When the full sequence has been process with this method endSeq should be called.
   * @param data a raw UTF-8 encoded buffer containing part of a sequence
   * @param start the start index of the part of the sequence in the buffer
   * @param length the length of the sequence part in the buffer
   * @see #endSeq()
   */
  public void addSeqPart(byte[] data, int start, int length) {
    if (this.frozen) throw new IllegalStateException("SequenceStats is frozen.");
    this.seqLength += length;
    int end = start + length;
    for (int i = start; i < end; i++) {
      byte codepoint = data[i];
      if (codepoint >= 0 && codepoint <= 127) {
        this.asciiCounts[codepoint]++;
      } else if ((codepoint & 0xC0) != 0x80) { // skip continuation bytes
        this.asciiCounts[128]++;
      }
    }
  }

  /**
   * Complete processing of a sequence and update the statistics to include the sequence.
   */
  public void endSeq() {
    if (this.frozen) throw new IllegalStateException("SequenceStats is frozen.");
    this.stats.update(this.seqLength);
    this.seqLength = 0;
  }

  /**
   * Stop further changes.
   * This allows the sequenceStat object to be returned without making
   * a copy as it becomes immutable.
   */
  public void freeze() {
    if (!this.frozen) {
      this.frozen = true;
      this.alphabetGuess = AlphStd.guess(this.asciiCounts);
    }
  }

  /**
   * Return the total count of sequences
   * @return The total count of sequences
   */
  public long getSequenceCount() {
    return this.stats.getCount();
  }

  /**
   * Return the total length of the sequence data
   * @return The total length of the sequence data
   */
  public long getTotalLength() {
    return this.stats.getTotal().longValueExact();
  }

  /**
   * Return the length of the smallest sequence seen or 0 when none have been seen.
   * @return The length of the smallest sequence
   */
  public long getMinLength() {
    return (long)this.stats.getSmallest();
  }

  /**
   * Return the length of the largest sequence seen or 0 when none have been seen.
   * @return The length of the largest sequence
   */
  public long getMaxLength() {
    return (long)this.stats.getLargest();
  }

  /**
   * Returns the average of all the sequence lengths seen or 0 when none have
   * been seen.
   * @return The average sequence length
   */
  public double getAverageLength() {
    return this.stats.getMean();
  }

  /**
   * Returns the standard deviation of all the sequence lengths seen or 0 when
   * less than 2 sequence lengths have been seen.
   * @return The standard deviation of all sequence lengths
   */
  public double getStandardDeviationLength() {
    return this.stats.getStandardDeviation();
  }

  /**
   * Try to guess the alphabet of the sequences by looking at the character frequencies.
   *
   * @return The most likely alphabet of the sequences
   */
  public AlphStd guessAlphabet() {
    if (this.frozen) return this.alphabetGuess;
    return AlphStd.guess(this.asciiCounts);
  }

  /**
   * Check that the alphabet would be acceptable by checking the symbols used in the sequence.
   *
   * @param alphabet the alphabet to test.
   * @return true if all symbols in the sequence are valid for the given alphabet.
   */
  @Override
  public boolean checkAlphabet(Alph alphabet) {
    return alphabet.isValidSymSet(this.asciiCounts);
  }

  /**
   * Check that one of the alphabets would be acceptable by checking the symbols used in the sequence.
   * @param alphabets the alphabets to test.
   * @return true if any of alphabets are acceptable.
   */
  public boolean checkAlphabets(Collection<Alph> alphabets) {
    for (Alph alphabet : alphabets) if (checkAlphabet(alphabet)) return true;
    return false;
  }
  /**
   * Check that one of the alphabets would be acceptable by checking the symbols used in the sequence.
   * @param alphabets the alphabets to test.
   * @return true if any of alphabets are acceptable.
   */
  public boolean checkAlphabetsEn(Collection<AlphStd> alphabets) {
    for (AlphStd alphabet : alphabets) if (checkAlphabet(alphabet.getAlph())) return true;
    return false;
  }

  /**
   *
   * @param type the type of alphabet that should be accepted.
   * @param restriction a specific restriction on the type of alphabet if specified.
   * @return true if a standard alphabet specified by the type matches or if an alphabet specified by the restriction matches.
   */
  public boolean checkAlphabets(AlphType type, Collection<Alph> restriction) {
    if (restriction != null && restriction.size() > 0) {
      for (Alph candidate : restriction) {
        if (AlphType.matches(type, candidate)) {
          if (checkAlphabet(candidate)) return true;
        }
      }
    } else {
      return checkAlphabetsEn(type.getStandardAlphabets());
    }
    return false;
  }

  @Override
  public void outputJson(JsonWr out) throws IOException {
    out.startObject();
    out.property("type", "sequence");
    out.property("alphabet", guessAlphabet().name());
    out.property("count", getSequenceCount());
    out.property("min", getMinLength());
    out.property("max", getMaxLength());
    out.property("avg", getAverageLength());
    out.property("total", getTotalLength());
    out.endObject();
  }
}
