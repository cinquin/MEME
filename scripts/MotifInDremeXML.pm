package MotifInDremeXML;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();

use DremeSAX;

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

  my $sax = new DremeSAX($self, 
    'handle_alphabet' => \&_handle_alphabet,
    'handle_strands' => \&_handle_strands,
    'handle_background' => \&_handle_background,
    'start_motif' => \&_start_motif,
    'end_motif' => \&_end_motif,
    'handle_pos' => \&_handle_pos
  );
  $self->{parser} = $sax;
  $self->{queue} = [];
  $self->{motif} = undef;
  $self->{alph} = undef;
  $self->{strands} = undef;
  $self->{bg} = {};

}

sub _handle_alphabet {
  my ($self, $alphabet) = @_;
  $self->{alph} = $alphabet;
  $self->{bg}->{alph} = $self->{alph};
}

sub _handle_strands {
  my $self = shift;
  my ($strands) = @_;
  $self->{strands} = ($strands eq 'both' ? 2 : ($strands eq 'forward' ? 1 : 0));
}

sub _handle_background {
  my $self = shift;
  my (@probs) = @_;
  my $bg = $self->{bg};
  my $alph = $self->{alph};
  for (my $i = 0; $i < $alph->size_core(); $i++) {
    my $sym = $alph->char($i);
    $bg->{$sym} = $probs[$i];
  }
}

sub _start_motif {
  my $self = shift;
  my ($id, $alt, $seq, $length, $nsites, $p, $n, $pvalue, $evalue, $unerased_evalue) = @_;
  
  my %pspm = ();
  my $alph = $self->{alph};
  for (my $i = 0; $i < $alph->size_core(); $i++) {
    my $sym = $alph->char($i);
    $pspm{$sym} = [];
  }

  my $motif = {
    bg => $self->{bg},
    strands => $self->{strands},
    id => $seq, 
    alt => $alt,
    width => $length, 
    sites => $nsites,
    pseudo => 0,
    evalue => $evalue,
    pspm => \%pspm,
  };

  $self->{motif} = $motif;
}

sub _end_motif {
  my $self = shift;
  push(@{$self->{queue}}, $self->{motif});
  $self->{motif} = undef;
}

sub _handle_pos {
  my $self = shift;
  my (@probs) = @_;
  my $alph = $self->{alph};
  my $pspm = $self->{motif}->{pspm};
  for (my $i = 0; $i < $alph->size_core(); $i++) {
    my $sym = $alph->char($i);
    push(@{$pspm->{$sym}}, $probs[$i]);
  }
}

1;
