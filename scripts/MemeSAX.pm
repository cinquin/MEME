package MemeSAX;

use base 'CheckingSAX';

use strict;
use warnings;

use CheckingSAX;
use Alphabet;

my $any_re = qr/^(.*)$/;
my $sci_re = qr/[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?/;
my $sci_or_inf_re = qr/(?:$sci_re|inf)/;
my $num_re = qr/^([+-]?$sci_or_inf_re)$/;
my $num_trim_re = qr/^\s*([+-]?$sci_or_inf_re)\s*$/;
my $pos_num_re = qr/^([+]?$sci_or_inf_re)$/;
my $pos_num_trim_re = qr/^\s*([+]?$sci_or_inf_re)\s*$/;
my $int_re = qr/^(\d+)$/;
my $int_trim_re = qr/^\s*(\d+)\s*$/;
my $float_re = qr/^(\d+(?:\.\d+)?)$/;
my $float_trim_re = qr/^\s*(\d+(?:\.\d+)?)\s*$/;
my $trim_re = qr/^\s*(.*?)\s*$/;
my $text_re = qr/^(.*\S.*)$/;
my $nospace_re = qr/^(\S+)$/;
my $sym_re = qr/[A-Za-z0-9\?\.\*\-]/;
my $sym_id_re = qr/^([a-zA-Z]|n[0-9]|x[a-zA-Z0-9]{2})$/;

my $ST_IN_SSITES = [{ELE => 'scanned_site', RPT => $RPT_ANY, ATRS => [
      {ATR => 'motif_id', VAL => qr/^motif_(\d+)$/},
      {ATR => 'strand', VAL => qr/^(plus|minus|none)$/},
      {ATR => 'position', VAL => $int_re},
      {ATR => 'pvalue', VAL => \&_value_pvalue}]}];
my $ST_IN_SSS = [{ELE => 'scanned_sites', RPT => $RPT_ALO, TO => $ST_IN_SSITES, ATRS => [
      {ATR => 'sequence_id', VAL => qr/^sequence_(\d+)$/},
      {ATR => 'pvalue', VAL => \&_value_pvalue},
      {ATR => 'num_sites', VAL => $int_re}]}];
my $ST_IN_MODEL = [
  &_struc_setting('command_line', $trim_re),
  &_struc_setting('host', $trim_re),
  &_struc_setting('type', \&_value_modtype),
  &_struc_setting('nmotifs', $int_trim_re),
  &_struc_setting('evalue_threshold', \&_value_evalue),
  &_struc_setting('object_function', $trim_re), #TODO stricter RE
  &_struc_setting('use_llr', $trim_re),
  &_struc_setting('min_width', $int_trim_re),
  &_struc_setting('max_width', $int_trim_re),
  # These optional values were in Tim's experimental "new" MEME branch...
  &_struc_setting('mask_type', $any_re, 1),
  &_struc_setting('min_key', $any_re, 1),
  &_struc_setting('max_key', $any_re, 1),
  &_struc_setting('negfile', $any_re, 1),
  &_struc_setting('kmer', $any_re, 1),
  &_struc_setting('ncopies', $any_re, 1),
  # End of "new" (currently unused) values
  &_struc_setting('wg', $num_trim_re, 1),
  &_struc_setting('ws', $num_trim_re, 1),
  &_struc_setting('endgaps', \&_value_yesno, 1),
  &_struc_setting('substring', \&_value_yesno),
  &_struc_setting('minsites', $int_trim_re),
  &_struc_setting('maxsites', $int_trim_re),
  &_struc_setting('wnsites', $num_trim_re),
  &_struc_setting('spmap', qr/^\s*(pam|uni)\s*$/),
  &_struc_setting('spfuzz', $num_trim_re),
  &_struc_setting('prior', qr/^\s*(addone|dirichlet|dmix|megadmix|megap)\s*$/),
  &_struc_setting('beta', $num_trim_re),
  &_struc_setting('maxiter', $int_trim_re),
  &_struc_setting('distance', $num_trim_re),
  &_struc_setting('num_sequences', $int_trim_re),
  &_struc_setting('num_positions', $int_trim_re),
  &_struc_setting('seed', qr/^\s*(\-?\d+)\s*$/),
  &_struc_setting('ctfrac', $num_trim_re),
  &_struc_setting('maxwords', qr/^\s*(\-?\d+)\s*$/),
  &_struc_setting('strands', qr/^\s*(both|forward|none)\s*$/),
  &_struc_setting('priors_file', $trim_re),
  &_struc_setting('reason_for_stopping', $trim_re),
  &_struc_alph_freqs('background_frequencies', [{ATR => 'source', VAL => $text_re}])
];
my $ST_IN_TRSET = [ # training set
  &_struc_alphabet(ID2SYM_AS => 'id_to_sym', SYM2ID_AS => 'sym_to_id', ALPH_AS => 'alph'),
  {ELE => 'sequence', RPT => $RPT_ALO, ATRS => [
      {ATR => 'id', VAL => qr/^sequence_(\d+)$/},
      {ATR => 'name', VAL => $text_re},
      {ATR => 'length', VAL => $int_re},
      {ATR => 'weight', VAL => $pos_num_re}]},
  &_struc_alph_freqs('letter_frequencies')
];
my $ST_IN_MEME = [
  {ELE => 'training_set', TO => $ST_IN_TRSET, ATRS => [
      {ATR => 'datafile', VAL => $text_re},
      {ATR => 'length', VAL => $int_re}]},
  {ELE => 'model', TO => $ST_IN_MODEL, SILENT => 1, END => \&_call_handle_settings},
  {ELE => 'motifs', TO => [&_struc_motif()]},
  {ELE => 'scanned_sites_summary', RPT => $RPT_OPT, TO => $ST_IN_SSS, ATRS => [
      {ATR => 'p_thresh', VAL => $num_re}]}];
my $ST_START = [
  {ELE => 'MEME', NAME => 'meme', RPT => $RPT_ONE, TO => $ST_IN_MEME, ATRS => [
      {ATR => 'version', VAL => qr/^(\d+)\.(\d+)\.(\d+)$/, EXPAND => 3},
      {ATR => 'release', VAL => $text_re}
    ]}];

#
# Return the defininition of a MEME XML file.
#
sub _xml_def {
  return $ST_START;
}

#
# Check value is either anr, oops, tcm or zoops.
# If anr is specified then convert it to tcm before returning.
#
sub _value_modtype {
  my ($self, $value) = @_;
  if ($value =~ qr/^\s*(anr|oops|tcm|zoops)\s*$/) {
    my $type = $1;
    $type = 'tcm' if ($type eq 'anr');
    return $type;
  }
  die("expected anr, oops, tcm or zoops\n");
}

sub _struc_setting {
  my ($element, $check, $opt) = @_;
  my $fn_store = sub {
    my ($self, $value) = @_;
    $self->{settings}->{$element} = $value;
  };
  return (
    {
      ELE => $element, 
      RPT => ($opt ? $RPT_OPT : $RPT_ONE),
      VAL => $check,
      SILENT => 1,
      END => $fn_store
    }
  );
}

sub _struc_contributing_site {
  my ($lflank, $seq, $rflank);
  my $fn_init = sub {
    my ($self, $attrs) = @_;
    $lflank = '';
    $seq = '';
    $rflank = '';
  };
  my $fn_left_flank = sub {
    my ($self, $left_flank) = @_;
    $lflank = $left_flank;
  };
  my $fn_letter_ref = sub {
    my ($self, $letter_id) = @_;
    $seq .= $self->{id_to_sym}->{$letter_id};
  };
  my $fn_right_flank = sub {
    my ($self, $right_flank) = @_;
    $rflank = $right_flank;
  };
  my $fn_contributing_site = sub {
    my ($self, $seq_num, $pos, $strand, $pvalue) = @_;
    return ($seq_num, $pos, $strand eq 'minus', $pvalue, $lflank, $seq, $rflank);
  };
  return (
    {
      ELE => 'contributing_site',
      NAME => 'site',
      HANDLE => $fn_contributing_site,
      WATCH => $fn_init,
      RPT => $RPT_ALO,
      ATRS => [
        {ATR => 'sequence_id', VAL => qr/^sequence_(\d+)$/},
        {ATR => 'position', VAL => $int_re},
        {ATR => 'strand', VAL => qr/^(plus|minus|none)$/},
        {ATR => 'pvalue', VAL => \&_value_pvalue}
      ],
      TO => [
        {ELE => 'left_flank', VAL => qr/^\s*($sym_re*)\s*$/, END => $fn_left_flank},
        {
          ELE => 'site',
          TO => [
            {
              ELE => 'letter_ref', 
              RPT => $RPT_ALO,
              END => $fn_letter_ref,
              ATRS => [
                {ATR => 'letter_id', VAL => $sym_id_re}
              ]
            }
          ]
        },
        {ELE => 'right_flank', VAL => qr/^\s*($sym_re*)\s*$/, END => $fn_right_flank}
      ]
    }
  );
}

sub _struc_alph_array {
  my ($callback, %options) = @_;
  my %values;
  my $fn_init = sub {
    %values = ();
  };
  my $fn_letter_value = sub {
    my ($self, $id, $value) = @_;
    if (defined($values{$id})) {
      $self->_error("Value repeated for letter ID \"$id\".");
    }
    $values{$id} = $value;
  };
  my $fn_result = sub {
    my ($self) = @_;
    my @alph_array = ();
    for (my $i = 0; $i < $self->{alph}->size_core(); $i++) {
      my $sym = $self->{alph}->char($i);
      my $id = $self->{sym_to_id}->{$sym};
      my $value = $values{$id};
      unless (defined($value)) {
        $self->_error("Value missing for letter ID \"$id\".");
        $value = 0;
      }
      $alph_array[$i] = $value;
    }
    $callback->(@alph_array);
  };
  return (
    {
      ELE => 'alphabet_array',
      RPT => $RPT_ALO,
      SILENT => 1,
      START => $fn_init,
      END => $fn_result,
      TO => [
        {
          ELE => 'value',
          RPT => $RPT_ALO,
          SILENT => 1,
          END => $fn_letter_value,
          ATRS => [
            {ATR => 'letter_id', VAL => $sym_id_re}
          ],
          VAL => $num_trim_re
        }
      ]
    }
  );
}
sub _struc_alph_freqs {
  my ($element, $atrs) = @_;
  $atrs = [] unless defined $atrs;
  my @freqs;
  my $fn_array = sub {
    @freqs = @_;
  };
  my $fn_result = sub {
    my ($self, @atrs) = @_;
    return (@atrs, \@freqs);
  };
  return (
    {
      ELE => $element,
      HANDLE => $fn_result,
      TO => [&_struc_alph_array($fn_array)], 
      ATRS => $atrs
    }
  );
}

sub _struc_alph_matrix {
  my ($callback, %options) = @_;
  my @alph_matrix;
  my $fn_init = sub {
    @alph_matrix = ();
  };
  my $fn_alph_array = sub {
    my @alph_array = @_;
    push(@alph_matrix, \@alph_array);
  };
  my $fn_result = sub {
    $callback->(@alph_matrix);
  };
  return (
    {
      ELE => 'alphabet_matrix',
      SILENT => 1,
      START => $fn_init,
      END => $fn_result,
      TO => [
        &_struc_alph_array($fn_alph_array)
      ]
    }
  );
}
sub _struc_motif {
  my ($motif_num, $name, $alt, $width, $sites, $llr, $ic, $re, $bayes_threshold, $evalue, $elapsed_time, $url, @scores, @probs, $regexp);
  my $fn_motif_init = sub {
    my ($self, $atrs);
    ($self, $atrs, $motif_num, $name, $alt, $width, $sites, $llr, $ic, $re, $bayes_threshold, $evalue, $elapsed_time, $url) = @_;
  };
  my $fn_scores = sub {
    @scores = @_;
  };
  my $fn_probs = sub {
    @probs = @_;
  };
  my $fn_regexp = sub {
    my ($self, $re) = @_;
    $regexp = $re;
  };
  my $fn_start_motif = sub {
    my ($self, $attrs) = @_;
    return (IDX => $motif_num - 1, NAME => $name, ALT => $alt, WIDTH => $width, SITES => $sites,
      LLR => $llr, IC => $ic, RE => $re, BAYES => $bayes_threshold,
      EVALUE => $evalue, ELAPSED => $elapsed_time, URL => $url,
      REGEXP => $regexp, PSM => \@scores, PWM => \@probs);
  };
  return (
    {
      ELE => 'motif',
      RPT => $RPT_ANY,
      SILENT => 1,
      START => $fn_motif_init,
      ATRS => [
        {ATR => 'id', VAL => qr/^motif_(\d+)$/},
        {ATR => 'name', VAL => $text_re},
        {ATR => 'alt', VAL => $text_re},
        {ATR => 'width', VAL => $int_re},
        {ATR => 'sites', VAL => $int_re},
        {ATR => 'llr', VAL => $num_re},
        {ATR => 'ic', VAL => $num_re},
        {ATR => 're', VAL => $num_re},
        {ATR => 'bayes_threshold', VAL => $num_re},
        {ATR => 'e_value', VAL => $num_re},
        {ATR => 'elapsed_time', VAL => $num_re},
        {ATR => 'url', OPT => qr/^(\S*)$/}
      ],
      TO => [
        {
          ELE => 'scores',
          SILENT => 1,
          TO => [&_struc_alph_matrix($fn_scores)]
        },
        {
          ELE => 'probabilities',
          SILENT => 1,
          TO => [&_struc_alph_matrix($fn_probs)]
        },
        {
          ELE => 'regular_expression',
          RPT => $RPT_OPT,
          SILENT => 1,
          END => $fn_regexp,
          VAL => qr/^\s*((?:$sym_re|\[$sym_re+\])+)\s*$/
        },
        {
          ELE => 'contributing_sites',
          NAME => 'motif', # by a bit of trickery we emit calls to start_motif, handle_site, end_motif
          START => $fn_start_motif,
          TO => [&_struc_contributing_site]
        }
      ]
    }
  );
}

sub _call_handle_settings {
  my ($self) = @_;
  my %settings = %{$self->{settings}};
  $self->_call('handle_settings', %settings);
}
