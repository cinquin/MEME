#File: MotifUtils.pm
#Project: generic
#Author: James Johnson (though used parts of exising scripts)
#Created: June 2010 (though the repackaged scripts are much older)

package MotifUtils;

use strict;
use warnings;


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(seq_to_intern matrix_to_intern intern_to_iupac intern_to_meme load_alphabet read_background_file meme_header motif_id motif_width parse_double round);

use Alphabet qw(dna protein);
use Fcntl qw(O_RDONLY);
use POSIX qw(strtod);

##############################################################################
# Internal format description
##############################################################################
#
# The internal format for the motif is as follows:
#
# hash {
#   bg => hash {
#     <residue 1> => scalar (background prior probability of residue 1)
#     <residue 2> => scalar (background prior probability of residue 2)
#      ...
#     <residue N> => scalar (background prior probability of residue N)
#     dna => scalar (true if the bg is dna)
#     source => scalar (source description of bg)
#   }
#   strands => scalar (number of strands used in a dna motif)
#   id => scalar (name of motif)
#   alt => scalar (alternate name of motif)
#   url => scalar (url of motif)
#   width => scalar (width of the motif)
#   sites => scalar (number of sites that made up motif)
#   pseudo => scalar (total pseudocount added)
#   evalue => scalar (evalue of motif)
#   pspm => hash {
#     <residue 1> => array [
#       scalar (probability of residue 1 at position 1)
#       scalar (probability of residue 1 at position 2)
#       ...
#       scalar (probability of residue 1 at position M)
#     ]
#     <residue 2> => array [
#       scalar (probability of residue 2 at position 1)
#       scalar (probability of residue 2 at position 2)
#       ...
#       scalar (probability of residue 2 at position M)
#     ]
#   }
#   pssm (optional for meme files that have a better pssm) => hash {
#     <residue 1> => array [
#       scalar (probability of residue 1 at position 1)
#       scalar (probability of residue 1 at position 2)
#       ...
#       scalar (probability of residue 1 at position M)
#     ]
#     <residue 2> => array [
#       scalar (probability of residue 2 at position 1)
#       scalar (probability of residue 2 at position 2)
#       ...
#       scalar (probability of residue 2 at position M)
#     ]
#   }
# }
#
#
#
#

##############################################################################
# public functions
##############################################################################

#
# Takes one or more iupac sequences with posssible bracket expressions 
# and constructs a probability matrix in this modules internal format. 
#
# Usage:
#  seq_to_intern(
#   \%background,
#   $sequence_samples,
#   $sites_per_sample,
#   $pseudo_total,
#   id => $optional_id,
#   alt => $optional_alt,
#   strands => $optional_strands,
#   url => $optional_url
#   );
#
# Parameters
#  bg                 - a reference to the background model also defines the alphabet
#  sequence samples   - one or more iupac sequences seperated by newlines.
#  sites per sample   - the number of sites attributed to each sample;
#  pseudo total       - a pseudocount that is distributed according to the background model
#
# Options
#  id                 - use value for motif name
#  alt                - use value for motif alternate name 
#  strands            - only relevant to dna; the strands used to create the motif
#  url                - use value for motif information website
#  also               - other characters that should be accepted, will be treated as ANY char
#
sub seq_to_intern {
  my ($bg, $seq_lines, $sites_per_sample, $pseudo_total, %opt) = @_;
  my ($id, $alt, $strands, $url, $also) = ($opt{id}, $opt{alt}, $opt{strands}, $opt{url}, $opt{also});
  my $alph = $bg->{alph};
  $strands = ($alph->has_complement() ? 2 : 1) unless defined $strands;
  $alt = '' unless defined($alt);
  # split the sequence lines into many sequences
  my @seqs = split(/\n/, $seq_lines);

  #keep track of errors related to parsing seq
  my @errors = ();

  #keep track of counts from multple samples
  my @matrix = ();
  my $name = '';

  # get the dimensions
  my $width = 0;

  my $sample_count = 0;
  for (my $line = 0; $line < scalar(@seqs); $line++) {
    my ($sets, $convert_errors) = seq_to_sets($seqs[$line], $alph, $also, $line + 1);
    next if (scalar(@{$sets}) == 0 && !@{$convert_errors}); #skip empty lines
    push(@errors, @{$convert_errors});
    if ($width == 0) {
      $width = scalar(@{$sets});
      if ($width > 0) { #usable line found
        #initilize matrix to all zeros
        my $matrix_size = $alph->size_core() * $width;
        for (my $i = 0; $i < $matrix_size; $i++) {
          $matrix[$i] = 0;
        }
        #create a name (more accurate than intern_to_iupac for protein)
        $name = sets_to_seq($alph, $sets) unless $id;
      }
    }
    if ($width == 0 || scalar(@{$sets}) != $width) { #skip bad lines
      push(@errors, errmsg($line + 1, $seqs[$line], "Bad width. Sample skipped.", 0));
      next;
    }
    for (my $pos = 0; $pos < $width; $pos++) {
      my $set = $sets->[$pos];
      my $fract = $sites_per_sample / scalar(keys %{$set});
      foreach my $sym ( keys %{$set} ) {
        $matrix[$pos * $alph->size_core() + $alph->index($sym)] += $fract;
      }
    }
    $sample_count++;
  } # parse line

  if ($width == 0) {
    push(@errors, "Motif has no width.\n");
    return (undef, \@errors);
  }

  my $sites = $sample_count * $sites_per_sample;
  #initilize the internal motif datastructure
  my %motif = (id => $id , alt => $alt, url => $url, strands => $strands, width => $width, 
      sites => $sites, pseudo => $pseudo_total, bg => $bg, evalue => 0, pspm => {});
  my $motif_pspm = $motif{pspm};

  my $total = $sites + $pseudo_total;
  for (my $a = 0; $a < $alph->size_core(); $a++) {
    my $residue = $alph->char($a);
    my $pseudo_count = $bg->{$residue} * $pseudo_total;
    $motif_pspm->{$residue} = [];
    for (my $pos = 0; $pos < $width; $pos++) {
      $motif_pspm->{$residue}->[$pos] = ($pseudo_count + $matrix[$pos * $alph->size_core() + $a]) / $total;
    }
  }

  unless ($id) {
    if ($sample_count > 1) { # if multiple samples were used then convert to iupac
      $name = intern_to_iupac(\%motif);
    }
    $motif{id} = $name;
  }

  return (\%motif, \@errors);
}


#
# Converts a motif in matrix form into the internal format.
#
# Usage:
#  matrix_to_intern(
#   \%bg,
#   $matrix,
#   $orientation,
#   $site_count,
#   $pseudo_total,
#   alt => $alt
#   );
#
# Parameters
#  bg            - a reference to the background model also defines if motif is dna or protein
#  matrix        - a matrix; columns are space separated and rows are newline separated
#  orientation   - the matrix orientation; allowed values are 'auto', 'col' and 'row'
#  site count    - the number of sites in the matrix; ignored for count matrices
#  pseudo total  - a pseudocount that is distributed according to the background model
#
# Options
#  id            - use value for motif name
#  alt           - use value for motif alternate name 
#  strands       - only relevant to dna; the strands used to create the motif
#  url           - use value for motif information website
#  rescale       - force the use of the site count, even when the matrix is already counts
#
sub matrix_to_intern {
  my ($bg, $matrix, $orientation, $site_count, $pseudo_total, %opt) = @_;
  my ($id, $alt, $strands, $url, $rescale) = ($opt{id}, $opt{alt}, $opt{strands}, $opt{url}, $opt{rescale});
  my $alph = $bg->{alph};
  $strands = ($alph->has_complement() ? 2 : 1) unless defined $strands;
  $alt = '' unless defined($alt);
  $rescale = 0 unless defined($rescale);

  #initilize the internal motif datastructure
  #NB width refers to the motif width not the matrix width
  my %motif = (id => $id, alt => $alt, url => $url, strands => $strands, width => 0, 
      sites => 0, pseudo => $pseudo_total, bg => $bg, evalue => 0, pspm => {});
  my $motif_pspm = $motif{pspm};
  
  # lookups
  my $asize = $alph->size_core();

  # dimensions (set defaults)
  my $height = 0; 
  my $width = ($orientation eq 'row' ? $asize : -1); 

  # parse matrix into an array
  my $total = 0;
  my @matrix_array = ();
  my @errors = ();
  my @lines = split(/\n/, $matrix);
  for (my $line_i = 0; $line_i < scalar(@lines); $line_i++) {
    my $line = $lines[$line_i];
    $line =~ s/^\s+//;#trim left
    $line =~ s/\s+$//;#trim right
    # skip empty lines
    next if $line eq '';
    my @nums = split(/\s+/, $line);
    if ($width == -1) {
      $width = scalar(@nums); #expected width
    } elsif (scalar(@nums) != $width) {
      push(@errors, "Error: Expected $width elements on line $line_i but got " . scalar(@nums) . ".");
      return (undef, \@errors);
    }
    # count this row
    $height++;
    for (my $num_i = 0; $num_i < scalar(@nums); $num_i++) {
      my $num = parse_double($nums[$num_i]);
      if (not defined $num) {
        push(@errors, "Warning: Element $num_i is not a number.");
        $num = 0; #specify default
      }
      push(@matrix_array, $num);
      $total += $num;
    }
  }
  # check for an empty matrix
  if ($height == 0) {
    push(@errors, "Error: Empty matrix.");
    return (undef, \@errors);
  }
  # check if the dimensions are plausible
  if ($orientation eq 'col' && $height != $asize) {
    push(@errors, "Error: Expected $asize rows but got $height.");
    return (undef, \@errors);
  }
  if ($orientation eq 'auto' && $width != $asize && $height != $asize) {
    push(@errors, "Error: Expected either the row or column count to be the alphabet size $asize but got $width by $height.");
    return (undef, \@errors);
  }
  # now determine the orientation and transform the matrix to a row matrix
  if ($orientation eq 'auto') {
    if ($width == $asize && $height == $asize) {
      # need to guess orientation so use variance
      my $avg = $total / $asize;
      my $row_variance = 0;
      my $col_variance = 0;
      for (my $i = 0; $i < $asize; $i++) {
        my $row_sum = 0;
        my $col_sum = 0;
        for (my $j = 0; $j < $asize; $j++) {
          $row_sum += $matrix_array[$i * $width + $j];
          $col_sum += $matrix_array[$j * $width + $i];
        }
        $row_variance += ($row_sum - $avg) ** 2; # row sum to the power of 2
        $col_variance += ($col_sum - $avg) ** 2; # col sum to the power of 2;
      }
      $row_variance /= $asize;
      $col_variance /= $asize;
      if ($row_variance > $col_variance) { # probably a column matrix
        @matrix_array = transpose_matrix(\@matrix_array, $width, $height);
      }
    } elsif ($height == $asize) { # dimensions suggest column matrix
      @matrix_array = transpose_matrix(\@matrix_array, $width, $height);
      $height = $width;
      $width = $asize;
    }
  } elsif ($orientation eq 'col') {
    @matrix_array = transpose_matrix(\@matrix_array, $width, $height);
    $height = $width;
    $width = $asize;
  }
  #should now have a row matrix
  if ($rescale) {
    if (not defined($site_count)) {
      push(@errors, "Error: Rescale option requires a site count. Can't continue without a site count.");
      return (undef, \@errors);
    }
  } else {
    #detect probability matrix
    my $row_avg = $total / $height;
    if (int($row_avg + 0.5) > 1) { # probably a count matrix
      $site_count = int($row_avg + 0.5); #int truncates toward zero so this rounds for positive numbers
    } elsif (not defined($site_count)) {
      push(@errors, "Error: Expected count matrix but got probability matrix. Can't continue without site count.");
      return (undef, \@errors);
    }
  }

  #set sites and width
  $motif{sites} = $site_count;
  $motif{width} = $height;

  # normalise each row individually
  my $row_total = $site_count + $pseudo_total;
  for (my $y = 0; $y < $height; ++$y) {
    my $sum = 0; #calculate this row's sum
    for (my $x = 0; $x < $width; ++$x) {
      $sum += $matrix_array[$y*$width + $x];
    }
    if ($sum == 0) {
      push(@errors, "Error: Motif position $y summed to zero.");
      return (undef, \@errors);
    }
    # normalise
    for (my $x = 0; $x < $width; ++$x) { 
      $matrix_array[$y*$width + $x] /= $sum;
    }
    # copy to motif
    for (my $x = 0; $x < $width; ++$x) {
      my $residue = $alph->char($x);
      # calculate probabilities
      my $count = $matrix_array[$y*$width + $x] * $site_count;
      my $pseudocount = $bg->{$residue} * $pseudo_total;
      $motif_pspm->{$residue}->[$y] = ($count + $pseudocount) / $row_total;
    }
  }

  unless (defined $id) {
    $id = intern_to_iupac(\%motif);
    $motif{id} = $id;
  }

  # now make motif
  return (\%motif, \@errors);
}

#
# intern_to_iupac(
#   $motif
#   );
#
# Measures the eucledian distance between each motif position and each alphabet letter
# and selects the closest letter.
#
sub intern_to_iupac {
  my ($motif) = @_;
  my $alph = $motif->{bg}->{alph};
  my $out = '';
  # create a list of frequencies for each possible alphabet letter
  my @sym_options = ();
  for (my $i = 0; $i < $alph->size_full(); $i++) {
    my $sym = $alph->char($i);
    my $set = $alph->comprise($sym);
    my $freq = (1.0 / scalar(keys %{$set}));
    my %freqs = ();
    for (my $j = 0; $j < $alph->size_core(); $j++) {
      my $a = $alph->char($j);
      $freqs{$a} = ($set->{$a} ? $freq : 0);
    }
    push(@sym_options, {sym => $sym, freqs => \%freqs});
  }

  # determine the best match to each position
  for (my $i = 0; $i < $motif->{width}; $i++) {
    my $best_dist = 'inf';
    my $best_sym = '?';
    for (my $j = 0; $j < scalar(@sym_options); $j++) {
      my $sym_opt = $sym_options[$j];
      my $freqs = $sym_opt->{freqs};
      my $dist = 0;
      for (my $k = 0; $k < $alph->size_core(); $k++) {
        my $a = $alph->char($k);
        $dist += ($freqs->{$a} - $motif->{pspm}->{$a}->[$i]) ** 2;
      }
      $dist = $dist ** 0.5;
      if ($dist < $best_dist) {
        $best_dist = $dist;
        $best_sym = $sym_opt->{sym};
      }
    }
    $out .= $best_sym;
  }
  return $out;
}

#
# intern_to_meme(
#   $motif,
#   $add_pssm,
#   $add_pspm,
#   $add_header
#   );
#
#   Returns a motif in minimal meme format with the header optimal so multiple motifs can be concatenated.
#   Warning: MEME does special stuff with pseudo counts on proteins to generate the log odds matrix which 
#            this approach does not do!
#
sub intern_to_meme {
  my ($motif, $add_pssm, $add_pspm, $add_header) = @_;

  my %bg = %{$motif->{bg}};
  my $alph = $bg{alph};

  my $output = "";

  #
  # get the text for the PSPM (and PSSM)
  #
  my $log_odds = "log-odds matrix: alength= ".$alph->size_core().
      " w= ".$motif->{width}." E= ".$motif->{evalue}."\n";
  my $letter_prob = "letter-probability matrix: alength= ".$alph->size_core().
      " w= ".$motif->{width}." nsites= ".$motif->{sites}." E= ".$motif->{evalue}."\n";
  for (my $i = 0; $i < $motif->{width}; $i++) {
    for (my $j = 0; $j < $alph->size_core(); $j++) {
      my $a = $alph->char($j);
      #get the probability
      my $prob = $motif->{pspm}->{$a}->[$i];
      #append to the letter probability matrix
      $letter_prob .= sprintf("%10.6f\t", $prob);
      my $score;
      if ($motif->{pssm}) {
        $score = $motif->{pssm}->{$a}->[$i];
      } else {
        #calculate the log odds as log of zero is undefined give it a value of -10000 (TODO check does this make sense)
        $score = $prob > 0 ? round((log($prob/$bg{$a}) / log(2.0)) * 100) : -10000;
      }
      #append to the log odds matrix
      $log_odds .= sprintf("%6d\t", $score);
    }
    $letter_prob .= "\n";
    $log_odds .= "\n";
  }
  $letter_prob .= "\n";
  $log_odds .= "\n";

  # Print the motif.
  $output .= meme_header($motif->{bg}, $motif->{strands}) if $add_header;
  $output .= "MOTIF ".$motif->{id}.(defined($motif->{alt}) ? " ".$motif->{alt} : "")."\n\n";
  $output .= $log_odds if ($add_pssm);
  $output .= $letter_prob if ($add_pspm);
  $output .= "URL ".$motif->{url}."\n\n" if ($motif->{url});

  return($output);
}


#
# load_alphabet()
#
# Loads an alphabet using the default to pick DNA or Protein when
# no alphabet file is specified.
#
sub load_alphabet {
  my ($default_dna, $file) = @_;
  if (defined($file)) {
    return new Alphabet($file);
  } elsif ($default_dna) {
    return dna();
  } else {
    return protein();
  }
}

#
# read_background_file()
#
# Read a background file in fasta-get-markov format.
# If $bg_file is not defined, sets background to uniform.
# Contains the key/value 'source' that has the background 
# source description.
#
# Returns background as a hash.
#
sub read_background_file {
  my ($alph, $bg_file) = @_;

  # get the letters in the alphabet 
  my (%bg);

  $bg{alph} = $alph;

  # initialize the background to uniform if no file given
  if (! defined $bg_file) {
    $bg{source} = "uniform background";
    my $freq = 1.0 / $alph->size_core();
    for (my $i = 0; $i < $alph->size_core(); $i++) {
      $bg{$alph->char($i)} = $freq;
    }
    return(%bg);
  }
  
  # read the background file
  $bg{source} = "file `$bg_file'";
  open(BG_FILE, $bg_file) || die("Can't open $bg_file.\n");
  my $total_bg = 0;
  while (<BG_FILE>) {
    next if (/^#/); # skip comments
    my ($a, $f) = split;
    next unless (length($a) == 1); # skip higher order model
    die("Symbol '$a' is not a core alphabet symbol") unless($alph->is_core($a));
    my $real_a = $alph->encode($a);
    die("Symbol '$a' equates to '$real_a' in this alphabet and another equalivent symbol has already been seen") if exists $bg{$real_a};
    $bg{$real_a} = $f;
    $total_bg += $f;
  }
  close BG_FILE;

  # make sure they sum to 1 by normalizing
  for (my $i = 0; $i < $alph->size_core(); $i++) {
    my $a = $alph->char($i);
    die("The letter '$a' was not given a value in the background file: $bg_file\n") unless exists $bg{$a};
    $bg{$a} /= $total_bg;
  }

  return(%bg);
}  # background file

#
# meme_header(
#   \%background,
#   $strands
#   );
#
#   Returns the header part of minimal meme format.
#   Strands is only used if the alphabet is complementable and defaults to 2.
#
#
sub meme_header {
  my ($bg, $strands) = @_;
  my $alph = $bg->{alph};
  $strands = 2 unless defined $strands;
  # Print the file MEME header.
  my $output;
  $output = "MEME version 4\n\n";
  if ($alph->is_dna() || $alph->is_rna() || $alph->is_protein()) {
    #$output .= $alph->to_text() . "\n";	// not necessary to use full alphabet defn
    my $core = $alph->get_core();
    $output .= "ALPHABET= $core\n\n";
  } else {
    $output .= $alph->to_text() . "\n";
  }
  if ($alph->has_complement()) {
    if ($strands == 2) {
      $output .= "strands: + -\n\n";
    } else {
      $output .= "strands: +\n\n";
    }
  }
  $output .= "Background letter frequencies (from " . ($bg->{source} ? $bg->{source} : "unknown source") . "):\n";
  for (my $i = 0; $i < $alph->size_core(); $i++) {
    my $a = $alph->char($i);
    $output .= sprintf("%s %.5f ", $a, $bg->{$a});
  }
  $output .= "\n\n";
  return($output);
}

#
# motif_id(
#   $motif
#   );
#
# Returns the motif id
#
sub motif_id {
  my $motif = shift;
  return $motif->{id};
}

#
# motif_width(
#   $motif
#   );
#
# Returns the width of the motif
#
sub motif_width {
  my $motif = shift;
  return $motif->{width};
}

#
# parse_double
#
# Parses a number using strtod
#
sub parse_double {
  my $str = shift;
  $str =~ s/^\s+//;#trim left
  $str =~ s/\s+$//;#trim right
  $! = 0; #clear errno
  my($num, $unparsed) = strtod($str);
  if (($str eq '') || ($unparsed != 0) || $!) {
      return undef;
  } else {
      return $num;
  } 
} 

sub round {
  # int truncates which is the same as always rounding toward zero
  if ($_[0] > 0) {
    return int($_[0] + 0.5);
  } else {
    return int($_[0] - 0.5);
  }
}

#
# transpose_matrix
#
# returns a transposed copy of the matrix.
#
sub transpose_matrix {
  my ($matrix_ref, $width, $height) = @_;
  my @matrix = @$matrix_ref;
  die("Dimensions don't match matrix size.") unless scalar(@matrix) == $width * $height;
  my @copy = @matrix;

  my $copy_width = $height;
  my $copy_height = $width;
  for (my $x = 0; $x < $width; ++$x) {
    for (my $y = 0; $y < $height; ++$y) {
      my $copy_x = $y;
      my $copy_y = $x;
      $copy[$copy_x + $copy_y * $copy_width] = $matrix[$x + $y * $width];
    }
  }
  return @copy;
}

#
# seq_to_sets
#
# returns an array of sets representing the sequence.
#
sub seq_to_sets {
  my ($seq, $alph, $also, $lineno) = @_;

  #error positions for making error msg
  my @unrecognised = ();
  my @unexpected_open = ();
  my @unexpected_close = ();
  my @nocontent = ();

  #parse the iupac codes and bracket expressions into an array of sets
  my @sets = ();

  my @chars = split(//, $seq);
  my $inside_bracket = 0;
  my %bracket_value = ();
  my $bracket_start;
  my $bracket_negate = 0;
  my %empty_set = ();
  for (my $i = 0; $i < scalar(@chars); $i++) {
    my $char = $chars[$i];
    if (!$inside_bracket) { #not inside a bracket expression
      if ($char eq '[') { #found the start of a bracket expression
        $inside_bracket = 1;
        %bracket_value = ();
        $bracket_start = $i;
        $bracket_negate = 0;
      } elsif (defined($also) && $char =~ m/^[\Q$also\E]$/) { # found an additional wildcard
        push(@sets, $alph->negate_set(\%empty_set));
      } elsif ($alph->is_known($char)) { #found a normal character
        push(@sets, $alph->comprise($char)); #add it to the sets
      } elsif ($char =~ m/^\s$/i) { #found whitespace
         #ignore
      } elsif ($char eq ']') { #misplaced bracket 
        push(@unexpected_close, $i); #add error
      } else { #unrecognised character
        push(@unrecognised, $i); #add error
      }
    } else { #inside a bracket expression
      if ($bracket_start == ($i - 1) && $char eq '^') { #negate bracket expression
        $bracket_negate = 1;
      } elsif ($char eq ']') { #found the end of a bracket expression
        $inside_bracket = 0;
        %bracket_value = %{$alph->negate_set(\%bracket_value)} if ($bracket_negate);
        if (%bracket_value) { #bracket had content so add it to the sets
          my %copy = %bracket_value;
          push(@sets, \%copy);
        } else { #no valid characters within the bracket expression
          push(@nocontent, $bracket_start);
        }
      } elsif (defined($also) && $char =~ m/^[\Q$also\E]$/) { # found an additional wildcard
        my $set = $alph->negate_set(\%empty_set);
        @bracket_value{keys %{$set}} = values %{$set};
      } elsif ($alph->is_known($char)) { #found a normal character
        #add to the set of options that the bracket expression could equal
        my $set = $alph->comprise($char);
        @bracket_value{keys %{$set}} = values %{$set};
      } elsif ($char =~ m/^\s$/i) { #found whitespace
        #ignore
      } elsif ($char eq '[') { #misplaced bracket
        push(@unexpected_open, $i); #add error
      } else { #unrecognised character
        push(@unrecognised, $i); #add error
      }
    }
  }
  if ($inside_bracket) { #missing close bracket
    if (%bracket_value) { #bracket had content so add it to the sets
      push(@sets, \%bracket_value);
    } else { #no valid characters within the bracket expression
      push(@nocontent, $bracket_start); #add error
    }
  }
  my @errors = ();
  push(@errors, errmsg($lineno, $seq, "Unrecognised", @unrecognised));
  push(@errors, errmsg($lineno, $seq, "Unexpected ']'", @unexpected_close));
  push(@errors, errmsg($lineno, $seq, "Unexpected '['", @unexpected_open));
  push(@errors, errmsg($lineno, $seq, "No content in '[...]'", @nocontent));
  push(@errors, errmsg($lineno, $seq, "Missing ']'", $bracket_start)) if $inside_bracket;
  return (\@sets, \@errors);
}

#
# errmsg
#
# outputs two lines, the first showing the line where the problem occured,
# the second showing the positions of the error with ^ and the error message.
#
sub errmsg {
  my ($lineno, $line, $msg, @errors) = @_;
  return () unless @errors;
  my $pos = 0;
  my $indent = ' ' x length("$lineno: ");
  my $str = "$lineno: $line\n$indent";
  foreach my $i (@errors) {
    next if $i < $pos;
    $str .= ' ' x ($i - $pos);
    $str .= '^';
    $pos = $i + 1;
  }
  $str .= "  $msg";
  return ($str);
}

#
# sets_to_seq
#
# Converts an array of sets into a equalivent sequence making use of IUPAC codes and regular expression bracket expressions.
#
sub sets_to_seq {
  my ($alph, $sets) = @_;
  # generate a seq
  my $seq = '';
  for (my $i = 0; $i < scalar(@{$sets}); $i++) {
    $seq .= $alph->regex($sets->[$i]);
  }
  return $seq;
}
