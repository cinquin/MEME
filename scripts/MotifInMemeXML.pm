package MotifInMemeXML;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

use MemeSAX;

sub new {
  my $classname = shift;
  my $self = {};
  bless($self, $classname);
  $self->_init(@_);
  return $self;
}

sub parse_more {
  my $self = shift;
  my ($buffer) = @_;
  $self->{parser}->parse_more($buffer);
}

sub parse_done {
  my $self = shift;
  $self->{parser}->parse_done();
}

sub has_errors {
  my $self = shift;
  return ($self->{parser}->has_errors() or $self->{alph}->has_errors());
}

sub get_errors {
  my $self = shift;
  return ($self->{parser}->get_errors(), $self->{alph}->get_errors());
}

sub get_motifs {
  my $self = shift;
  my @motifs = @{$self->{queue}};
  $self->{queue} = [];
  return @motifs;
}

sub has_motif {
  my $self = shift;
  return scalar(@{$self->{queue}}) > 0;
}

sub next_motif {
  my $self = shift;
  return pop(@{$self->{queue}});
}

sub _init {
  my $self = shift;

  my $sax = new MemeSAX($self, 
    handle_alphabet => \&_alphabet,
    handle_sequence => \&_sequence,
    handle_background_frequencies => \&_background,
    handle_settings => \&_settings,
    start_motif => \&_start_motif,
    handle_site => \&_site,
    end_motif => \&_end_motif,
  );
  $self->{parser} = $sax;
  $self->{queue} = [];
  $self->{strands} = undef;
  $self->{motif} = undef;
  $self->{seq_names} = []; # maps sequence indexes to names
  $self->{bg} = undef;
  $self->{alph} = undef;
  $self->{counter} = undef;
}

sub _alphabet {
  my $self = shift;
  my ($alphabet) = @_;
  $self->{alph} = $alphabet;
}

sub _sequence {
  my $self = shift;
  my ($idx, $name, $width, $weight) = @_;
  $self->{seq_names}->[$idx] = $name;
}

sub _settings {
  my $self = shift;
  my (%settings) = @_;
  $self->{strands} = ($settings{strands} eq 'both' ? 2 : ($settings{strands} eq 'forward' ? 1 : 0));
}

sub _background {
  my $self = shift;
  my ($source, $freqs) = @_;
  my %bg = (alph => $self->{alph}, source => $source);
  for (my $i = 0; $i < $self->{alph}->size_core(); $i++) {
    $bg{$self->{alph}->char($i)} = $freqs->[$i];
  }
  $self->{bg} = \%bg;
}

sub _remap_matrix {
  my ($alph, $matrix) = @_;
  my $len = scalar(@{$matrix});
  my %out = ();
  for (my $a = 0; $a < $alph->size_core(); $a++) {
    my @values = ();
    $out{$alph->char($a)} = \@values;
    for (my $c = 0; $c < $len; $c++) {
      $values[$c] = $matrix->[$c]->[$a];
    }
  }
  return \%out;
}

sub _start_motif {
  my $self = shift;
  my (%values) = @_;

  my $alph = $self->{alph};

  my $motif = {
    bg => $self->{bg},
    strands => $self->{strands},
    id => $values{NAME}, 
    alt => $values{ALT},
    url => $values{URL}, 
    width => $values{WIDTH}, 
    sites => $values{SITES},
    pseudo => 0,
    evalue => $values{EVALUE},
    pspm => &_remap_matrix($alph, $values{PWM}),
    pssm => &_remap_matrix($alph, $values{PSM}),
    contributing_sites => []
  };

  $self->{motif} = $motif;
}

sub _site {
  my $self = shift;
  my ($seq_idx, $pos, $is_rc, $pvalue, $lflank, $site, $rflank) = @_;
  push(@{$self->{motif}->{contributing_sites}}, {
    sequence => $self->{seq_names}->[$seq_idx],
    pos => $pos,
    strand => ($self->{alph}->has_complement() ? ($is_rc ? -1 : 1) : 0),
    pvalue => $pvalue,
    lflank => $lflank,
    site => $site,
    rflank => $rflank
  });
}

sub _end_motif {
  my $self = shift;
  push(@{$self->{queue}}, $self->{motif});
  $self->{motif} = undef;
}


1;
