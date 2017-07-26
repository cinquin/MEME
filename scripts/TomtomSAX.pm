package TomtomSAX;

use base 'CheckingSAX';

use strict;
use warnings;

use CheckingSAX;

my $num_re = qr/^((?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|inf)$/;
my $num_trim_re = qr/^\s*((?:[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)|inf)\s*$/;
my $int_re = qr/^(\d+)$/;
my $int_trim_re = qr/^\s*(\d+)\s*$/;
my $float_re = qr/^(\d+(?:\.\d+)?)$/;
my $float_trim_re = qr/^\s*(\d+(?:\.\d+)?)\s*$/;
my $trim_re = qr/^\s*(.*?)\s*$/;
my $text_re = qr/^(.*\S.*)$/;
my $nospace_re = qr/^(\S+)$/;

my $ST_IN_QUERY = [
  {ELE => 'target', RPT => $RPT_ANY, ATRS => [
      {ATR => 'idx', VAL => $int_re},
      {ATR => 'rc', OPT => \&_value_yn},
      {ATR => 'off', VAL => qr/^(-?\d+)$/},
      {ATR => 'pv', VAL => \&_value_pvalue},
      {ATR => 'ev', VAL => \&_value_evalue},
      {ATR => 'qv', VAL => \&_value_pvalue}]}
];
my $ST_IN_MATCHES = [
  {ELE => 'query', RPT => $RPT_ANY, TO => $ST_IN_QUERY, ATRS => [
      {ATR => 'idx', VAL => $int_re}]}
];
my $ST_IN_MOTIF = [
  {ELE => 'pos', RPT => $RPT_ALO, ATRS => [\&_get_alphabet_atrs]}, 
];
my $ST_IN_MOTIFS = [
  {ELE => 'motif', RPT => $RPT_ANY, TO => $ST_IN_MOTIF, ATRS => [
      {ATR => 'db', VAL => $int_re},
      {ATR => 'id', VAL => $nospace_re},
      {ATR => 'alt', OPT => $nospace_re},
      {ATR => 'length', VAL => $int_re},
      {ATR => 'nsites', OPT => $int_re},
      {ATR => 'evalue', OPT => \&_value_evalue},
      {ATR => 'url', OPT => $nospace_re}]}
];
my $ST_IN_DBS = [
  {ELE => 'db', RPT => $RPT_ALO, ATRS => [
      {ATR => 'source', VAL => $text_re},
      {ATR => 'name', VAL => $text_re},
      {ATR => 'loaded', VAL => $int_re},
      {ATR => 'excluded', VAL => $int_re},
      {ATR => 'last_mod_date', VAL => $text_re}]}
];
my $ST_IN_MODEL = [
  {ELE => 'command_line', VAL => $text_re},
  {ELE => 'distance_measure', ATRS => [
      {ATR => 'value', VAL => qr/^(pearson|allr|ed|kullback|sandelin|blic1|blic5|llr1|llr5)$/}]},
  {ELE => 'threshold', VAL => \&_value_evalue, ATRS => [
      {ATR => 'type', VAL => qr/^(qvalue|evalue)$/}]},
  &_struc_alphabet(ATRS_AS => 'alphabet_atrs'),
  {ELE => 'strands', VAL => qr/^(both|forward|none)$/},
  {ELE => 'background', ATRS => [
      {ATR => 'from', VAL => qr/^(first_target|file)$/},
      {ATR => 'file', OPT => $text_re},
      \&_get_alphabet_atrs
    ], WATCH => \&_watch_ele_background},
  {ELE => 'host', VAL => $trim_re},
  {ELE => 'when', VAL => $trim_re},
  {ELE => 'description', RPT => $RPT_OPT, VAL => $trim_re}
];
my $ST_IN_TOMTOM = [
  {ELE => 'model', TO => $ST_IN_MODEL}, 
  {ELE => 'query_dbs', TO => $ST_IN_DBS},
  {ELE => 'target_dbs', TO => $ST_IN_DBS},
  {ELE => 'queries', TO => $ST_IN_MOTIFS},
  {ELE => 'targets', TO => $ST_IN_MOTIFS},
  {ELE => 'matches', TO => $ST_IN_MATCHES},
  {ELE => 'runtime', ATRS => [
      {ATR => 'cycles', VAL => $int_re}, 
      {ATR => 'seconds', VAL => $float_re}]}
];
my $ST_START = [
  {ELE => 'tomtom', TO => $ST_IN_TOMTOM, ATRS => [
      {ATR => 'version', VAL => qr/^(\d+)\.(\d+)\.(\d+)$/, EXPAND => 3},
      {ATR => 'release', VAL => $text_re}
    ]}];

sub _xml_def {
  return $ST_START;
}

sub _init {
  my $self = shift;
  $self->{alphabet} = [];
  $self->SUPER::_init(@_);
}

sub _get_alphabet_atrs {
  return @{$_[0]->{alphabet_atrs}};
}

sub _watch_ele_background {
  my ($self, $attrs, $from, $file, @probs) = @_;
  if ($from eq 'file') {
    $self->_error('background/@file missing') unless defined $file;
  }
}
