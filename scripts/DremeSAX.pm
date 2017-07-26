package DremeSAX;

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
my $allowed_sym_re = qr/[!-~]/; # Allow all visible ASCII
my $sym_re = qr/^($allowed_sym_re)$/;
my $syms_re = qr/^($allowed_sym_re+)$/;
my $sym_id_re = qr/^([a-zA-Z]|n[0-9]|x[a-zA-Z0-9]{2})$/;

my $ST_IN_MOTIF = [
  {ELE => 'pos', RPT => $RPT_ALO, ATRS => [\&_get_alphabet_atrs]}, 
  #{ELE => 'match', RPT => $RPT_ALO, ATRS => [
  {ELE => 'match', RPT => $RPT_ANY, ATRS => [
      {ATR => 'seq', VAL => $syms_re},
      {ATR => 'p', VAL => $int_re},
      {ATR => 'n', VAL => $int_re},
      {ATR => 'pvalue', VAL => \&_value_pvalue},
      {ATR => 'evalue', VAL => \&_value_evalue}
    ]}
];
my $ST_IN_MOTIFS = [
  {ELE => 'motif', RPT => $RPT_ANY, TO => $ST_IN_MOTIF, ATRS => [
      {ATR => 'id', VAL => $nospace_re},
      {ATR => 'alt', VAL => $nospace_re},
      {ATR => 'seq', VAL => $syms_re},
      {ATR => 'length', VAL => $int_re},
      {ATR => 'nsites', VAL => $int_re},
      {ATR => 'p', VAL => $int_re},
      {ATR => 'n', VAL => $int_re},
      {ATR => 'pvalue', VAL => \&_value_pvalue},
      {ATR => 'evalue', VAL => \&_value_evalue},
      {ATR => 'unerased_evalue', VAL => \&_value_evalue}
    ]}];
my $ST_IN_MODEL = [
  {ELE => 'command_line', VAL => $text_re},
  {ELE => 'positives', ATRS => [
      {ATR => 'name'},
      {ATR => 'count', VAL => $int_re},
      {ATR => 'file'},
      {ATR => 'last_mod_date'}
    ]},
  {ELE => 'negatives', ATRS => [
      {ATR => 'name', VAL => $text_re},
      {ATR => 'count', VAL => $int_re},
      {ATR => 'from', VAL => qr/^(file|shuffled)$/},
      {ATR => 'file', OPT => $text_re},
      {ATR => 'last_mod_date', OPT => $text_re},
    ], WATCH => \&_watch_ele_negatives},
  &_struc_alphabet(ATRS_AS => 'alphabet_atrs'),
  {ELE => 'strands', VAL => qr/^(both|given|none)$/},
  {ELE => 'background', ATRS => [
      {ATR => 'from', VAL => qr/^(dataset)$/, SILENT => 1},
      \&_get_alphabet_atrs
    ]},
  {ELE => 'stop', ATRS => [
      {ATR => 'evalue', OPT => \&_value_evalue}, 
      {ATR => 'count', OPT => $int_re}, 
      {ATR => 'time', OPT => $int_re}]},
  {ELE => 'ngen', VAL => $int_trim_re},
  {ELE => 'add_pv_thresh', VAL => $num_trim_re},
  {ELE => 'seed', VAL => $int_trim_re},
  {ELE => 'host', VAL => $trim_re},
  {ELE => 'when', VAL => $trim_re},
  {ELE => 'description', RPT => $RPT_OPT, VAL => qr/^([\S\s]*)$/}
];
my $ST_IN_DREME = [
  {ELE => 'model', TO => $ST_IN_MODEL}, 
  {ELE => 'motifs', TO => $ST_IN_MOTIFS},
  {ELE => 'run_time', ATRS => [
      {ATR => 'cpu', VAL => $trim_re}, 
      {ATR => 'real', VAL => $float_re}, 
      {ATR => 'stop', VAL => qr/^(evalue|count|time)$/}]}
];
my $ST_START = [
  {ELE => 'dreme', TO => $ST_IN_DREME, ATRS => [
      {ATR => 'version', VAL => qr/^(\d+)\.(\d+)\.(\d+)$/, EXPAND => 3},
      {ATR => 'release', VAL => $text_re}
    ]}];

sub _xml_def {
  return $ST_START;
}

sub _init {
  my $self = shift;
  $self->{alphabet_atrs} = [];
  $self->SUPER::_init(@_);
}

sub _get_alphabet_atrs {
  return @{$_[0]->{alphabet_atrs}};
}

sub _watch_ele_negatives {
  my ($self, $attrs, $name, $count, $from, $file, $last_mod) = @_;
  if ($from eq 'file') {
    $self->_error('negatives/@file missing') unless defined $file;
    $self->_error('negatives/@last_mod missing') unless defined $last_mod;
  }
}
