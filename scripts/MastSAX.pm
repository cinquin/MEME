package MastSAX;

use base 'CheckingSAX';

use strict;
use warnings;

use CheckingSAX;
use Data::Dumper;
use Alphabet;

my $any_re = qr/^(.*)$/;
my $any_multiline_re = qr/^([\S\s]*)$/;
my $sci_re = qr/[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?/;
my $sci_or_inf_re = qr/(?:$sci_re|inf)/;
my $num_re = qr/^([+-]?$sci_or_inf_re)$/;
my $int_re = qr/^(\d+)$/;
my $float_re = qr/^(\d+(?:\.\d+)?)$/;
my $text_re = qr/^(.*\S.*)$/;
my $nospace_re = qr/^(\S+)$/;

my $ST_IN_SEG = [
  {ELE => 'data', VAL => \&_value_sequence_seq},
  {ELE => 'hit', RPT => $RPT_ALO, ATRS => [
      {ATR => 'pos', VAL => $int_re},
      {ATR => 'idx', VAL => $int_re},
      {ATR => 'rc', OPT => \&_value_yn},
      {ATR => 'pvalue', VAL => \&_value_pvalue},
      {ATR => 'match', VAL => qr/^([\+ ]+)$/},
      {ATR => 'translation', OPT => \&_value_motif_seq}
    ]}
];
my $ST_IN_SEQUENCE = [
  # score appears once normaly and twice when strands are scanned separately
  {ELE => 'score', RPT => $RPT_ALO, ATRS => [
      {ATR => 'strand', VAL => qr/^(both|forward|reverse)$/},
      {ATR => 'combined_pvalue', VAL => \&_value_pvalue},
      {ATR => 'evalue', VAL => \&_value_evalue},
      {ATR => 'frame', OPT => qr/^([a-z])$/}
    ]},
  {ELE => 'seg', RPT => $RPT_ANY, TO => $ST_IN_SEG, ATRS => [
      {ATR => 'start', VAL => $int_re}
    ]}
];
my $ST_IN_MAST = [
  &_struc_command_line(),
  {ELE => 'settings', ATRS => [
      {ATR => 'strand_handling', VAL => qr/(combine|separate|norc|unstranded)/},
      {ATR => 'max_correlation', VAL => $num_re},
      {ATR => 'remove_correlated', VAL => \&_value_yn},
      {ATR => 'max_seq_evalue', VAL => \&_value_evalue},
      {ATR => 'adjust_hit', VAL => \&_value_yn},
      {ATR => 'max_hit_pvalue', VAL => \&_value_pvalue},
      {ATR => 'max_weak_pvalue', VAL => \&_value_pvalue}
    ]},
  &_struc_alphabet(ALPH_AS => 'alph', ATRS_AS => 'alph_atrs'),
  &_struc_alphabet(ELE => 'sequence_alphabet', RPT => $RPT_OPT, ALPH_AS => 'seq_alph'),
  {ELE => 'translate', RPT => $RPT_OPT, ATRS => [
      {ATR => 'num_sequence', VAL => $int_re},
      {ATR => 'num_motif', VAL => $int_re}
    ]},
  {ELE => 'background', ATRS => [
      #{ATR => 'from', VAL => qr/^(sequence|file|preset)$/},
      #{ATR => 'file', OPT => $text_re},
      {ATR => 'source', OPT => $text_re},
      \&_get_alphabet_atrs
    ], WATCH => \&_watch_ele_background},
  {ELE => 'motif_dbs', TO => [
    {ELE => 'motif_db', RPT => $RPT_ALO, HANDLE => \&_handle_ele_motif_db,
      ATRS => [
        {ATR => 'source', VAL => $text_re},
        {ATR => 'name', OPT => $text_re},
        {ATR => 'last_mod_date', VAL => $text_re}
      ], TO => [
        {ELE => 'background', ATRS => [
            \&_get_alphabet_atrs
          ], WATCH => \&_watch_ele_motif_db_background},
      ]}
    ]},
  {ELE => 'motifs', TO => [
    {ELE => 'motif', RPT => $RPT_ANY, HANDLE => \&_handle_ele_motif,
      ATRS => [
        {ATR => 'db', VAL => $int_re},
        {ATR => 'id', VAL => $nospace_re},
        {ATR => 'alt', OPT => $nospace_re},
        {ATR => 'length', VAL => $int_re},
        {ATR => 'nsites', OPT => $int_re},
        {ATR => 'evalue', OPT => \&_value_evalue},
        {ATR => 'bad', OPT => \&_value_yn},
        {ATR => 'url', OPT => $nospace_re}
      ], WATCH => \&_watch_ele_motif,
      TO => [
        {ELE => 'pos', RPT => $RPT_ALO, ATRS => [\&_get_alphabet_atrs], WATCH => \&_watch_ele_pos}, 
      ]
    }
  ]},
  {ELE => 'correlations', TO => [
    {ELE => 'correlation', RPT => $RPT_ANY, ATRS => [
        {ATR => 'idx_a', VAL => $int_re},
        {ATR => 'idx_b', VAL => $int_re},
        {ATR => 'value', VAL => $num_re}
      ]},
    ]},
  {ELE => 'nos', RPT => $RPT_OPT, TO => [
    {ELE => 'expect', RPT => $RPT_ALO, ATRS => [
        {ATR => 'gap', OPT => $int_re},
        {ATR => 'idx', VAL => $int_re}
      ]}
    ]},
  {ELE => 'sequence_dbs', TO => [
    {ELE => 'sequence_db', RPT => $RPT_ALO, ATRS => [
        {ATR => 'source', VAL => $text_re},
        {ATR => 'name', OPT => $text_re},
        {ATR => 'last_mod_date', VAL => $text_re},
        {ATR => 'seq_count', VAL => $int_re},
        {ATR => 'residue_count', VAL => $int_re},
        {ATR => 'link', OPT => $text_re},
      ]}
    ]},
  {ELE => 'sequences', TO => [
    {ELE => 'sequence', RPT => $RPT_ANY, TO => $ST_IN_SEQUENCE, ATRS => [
        {ATR => 'db', VAL => $int_re},
        {ATR => 'name', VAL => $text_re},
        {ATR => 'comment', VAL => $any_re},
        {ATR => 'length', VAL => $int_re}
      ]}
    ]},
  {ELE => 'runtime', ATRS => [
      {ATR => 'host', VAL => $text_re},
      {ATR => 'when', VAL => $text_re},
      {ATR => 'cycles', VAL => $int_re},
      {ATR => 'seconds', VAL => $float_re}
    ]}
];

my $ST_START = [
  {ELE => 'mast', RPT => $RPT_ONE, TO => $ST_IN_MAST, ATRS => [
    {ATR => 'version', VAL => qr/^(\d+\.\d+\.\d+)$/},
    {ATR => 'release', VAL => $text_re}
    ]}];

sub _xml_def {
  return $ST_START;
}

sub _init {
  my $self = shift;
  $self->SUPER::_init(@_);
}

sub _get_alphabet_atrs {
  my ($self) = @_;
  # return the ordered list of alphabet attributes
  return @{$self->{alph_atrs}};
}

sub _watch_ele_background {
  my ($self, $attrs, $from, $file, @probs) = @_;
  # enforce the availabilty of the file attribute when the background is source from a file
  if ($from eq 'file') {
    $self->_error('background/@file missing') unless defined $file;
  }
}

sub _watch_ele_motif_db_background {
  my ($self, $attrs, @probs) = @_;
  $self->{motif_bg} = \@probs;
}

sub _handle_ele_motif_db {
  my ($self, $source, $name, $last_mod_date) = @_;
  return ($source, $name, $last_mod_date, $self->{motif_bg});
}

sub _watch_ele_motif {
  my ($self) = @_;
  $self->{matrix} = [];
}

sub _watch_ele_pos {
  my ($self, $attrs, @probs) = @_;
  push(@{$self->{matrix}}, [@probs]);
}

sub _handle_ele_motif {
  my ($self, $db, $id, $alt, $length, $nsites, $evalue, $bad, $url) = @_;
  return ($db, $id, $alt, $length, $nsites, $evalue, $bad, $url, $self->{matrix});
}

#
# Check value is valid sequence data according to the motif alphabet.
# Whitespace is removed.
#
sub _value_motif_seq {
  my ($self, $value) = @_;
  my $seq = $value;
  $seq =~ s/\s//g;
  my ($pos, $sym) = $self->{alph}->find_unknown($seq);
  die("contains non-alphabet symbol $sym\n") if $sym;
  return $seq;
}

#
# Check value is valid sequence data according to the sequence alphabet (which may also be the motif alphabet).
# Whitespace is removed.
#
sub _value_sequence_seq {
  my ($self, $value) = @_;
  my $seq = $value;
  $seq =~ s/\s//g;
  my $alph = (defined($self->{seq_alph}) ? $self->{seq_alph} : $self->{alph});
  my ($pos, $sym) = $alph->find_unknown($seq);
  die("contains non-alphabet symbol $sym\n") if $sym;
  return $seq;
}
