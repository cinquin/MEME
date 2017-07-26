package Alphabet;

use strict;
use warnings;
use feature 'unicode_strings';

use Carp;
use Fcntl qw(O_RDONLY);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(rna dna protein);

my $SYM_RE = qr/[A-Za-z0-9\?\.\*\-]/;
my $COLOUR_RE = qr/[0-9A-Fa-f]{6}/;
my $NAME_RE = qr/"((?:[^\\"]+|\\["\/bfnrt]|\\u[0-9A-Fa-f]{4})*)"/;
my $CORE_RE = qr/($SYM_RE)(?:\s+$NAME_RE)?(?:\s+($COLOUR_RE))?/;
my $COMMENT_RE = qr/^(\s*(?:[^#](?:[^"#]*|"(?:\\.|[^"\\]*)*")*)?)#.*$/;
my $HEADER_RE = qr/^\s*ALPHABET(?:\s+v1)?(?:\s+$NAME_RE)?(?:\s+(RNA|DNA|PROTEIN)-LIKE)?\s*$/;
my $CORE_SINGLE_RE = qr/^\s*$CORE_RE\s*$/;
my $CORE_PAIR_RE = qr/^\s*$CORE_RE\s*~\s*$CORE_RE\s*$/;
my $AMBIG_RE = qr/^\s*$CORE_RE\s*=\s*($SYM_RE*)\s*$/;

sub new {
  my $classname = shift;
  my $self = {};
  bless($self, $classname);
  $self->_init(@_);
  return $self;
}

sub _init {
  my $self = shift;
  my ($file) = @_;
  $self->{errors} = [];
  $self->{file_name} = undef;
  $self->{line_no} = 0;
  $self->{seen_header} = 0; #FALSE
  $self->{seen_symbol} = 0; #FALSE
  $self->{seen_ambig} = 0; #FALSE
  $self->{seen_lc} = 0; #FALSE
  $self->{seen_uc} = 0; #FALSE
  $self->{name} = undef;
  $self->{sym_lookup} = {};
  $self->{core_syms} = [];
  $self->{ambig_syms} = [];
  $self->{parsed} = 0; #FALSE
  if (defined($file)) {
    $self->parse_file($file);
    if (@{$self->{errors}}) {
      foreach my $error ( @{$self->{errors}} ) {
        warn($error);
      }
      croak("Invalid alphabet file \"$file\"\n");
    }
  }
}

sub _error {
  my $self = shift;
  my ($reason) = @_;
  my $line_info = '';
  my $file_info = '';
  $line_info = 'on line ' . $self->{line_no} . ' ' if ($self->{line_no});
  $file_info = 'file "' . $self->{file_name} . '"' if ($self->{file_name});
  push(@{$self->{errors}}, "Invalid alphabet $file_info- $line_info$reason\n");
}

sub parse_header {
  my $self = shift;
  my ($name, $like) = @_;
  die("Already finished parsing") if $self->{parsed};
  $self->_error('repeated header') if ($self->{seen_header});
  $self->_error('header after symbol') if ($self->{seen_symbol});
  $self->_error('like must be "DNA", "RNA" or "PROTEIN" (ignoring case)') if ($like && $like !~ m/^(RNA|DNA|PROTEIN)$/i);
  $self->{name} = $name;
  $self->{like} = ($like ? uc($like) : '');
  $self->{seen_header} = 1; #TRUE
}

sub uniq {
  my %seen;
  grep !$seen{$_}++, @_;
}

sub parse_symbol {
  my $self = shift;
  my ($sym, $name, $colour, $complement, $comprise, $aliases) = @_;
  die("Already finished parsing") if $self->{parsed};
  $self->_error('expected header but found symbol') unless ($self->{seen_header} || $self->{seen_symbol});
  if (defined($self->{sym_lookup}->{$sym})) {
    $self->_error("symbol $sym is already used");
    return; # can't really proceed with a duplicate
  }
  $self->{seen_lc} = 1 if ($sym ge 'a' && $sym le 'z');
  $self->{seen_uc} = 1 if ($sym ge 'A' && $sym le 'Z');
  my $symbol = {sym => $sym, aliases => [], name => $name, colour => $colour};
  if (defined($aliases)) {
    push(@{$symbol->{aliases}}, split(//, $aliases));
  }
  if (defined($comprise)) { # ambiguous symbol
    my @equals = split //, $comprise;
    @equals = sort _symbol_cmp @equals;
    @equals = &uniq(@equals);
    for (my $i = 0; $i < scalar(@equals); $i++) {
      my $csym = $self->{sym_lookup}->{$equals[$i]};
      if (!defined($csym)) {
        $self->_error('referenced symbol ' . $equals[$i] . ' is unknown');
      } elsif (defined($csym->{comprise})) {
        $self->_error('referenced symbol ' . $equals[$i] . ' is ambiguous');
      }
    }
    $symbol->{comprise} = join('', @equals);
    push(@{$self->{ambig_syms}}, $symbol);
    $self->{seen_ambig} = 1;
    if ($sym eq '?' && length($symbol->{comprise}) < scalar(@{$self->{core_syms}})) {
      $self->_error('symbol ? is reserved as a wildcard');
    }
  } else { # core symbol
    $self->_error('symbol ? is reserved as a wildcard') if ($sym eq '?');
    $self->_error('unexpected core symbol (as ambiguous symbols seen)') if ($self->{seen_ambig});
    $symbol->{complement} = $complement;
    push(@{$self->{core_syms}}, $symbol);
  }
  $self->{sym_lookup}->{$sym} = $symbol;
  $self->{seen_symbol} = 1;
}

sub parse_line {
  my $self = shift;
  die("Already finished parsing") if $self->{parsed};
  my ($line) = @_;
  $self->{line_no} += 1;
  $line =~ s/$COMMENT_RE/$1/;
  return unless ($line =~ m/\S/); # skip empty lines
  if ($line =~ m/$HEADER_RE/) {
    $self->parse_header($1, $2);
  } elsif ($line =~ m/$CORE_PAIR_RE/) {
    $self->parse_symbol($1, &_decode_JSON_string($2), &_decode_colour($3), $4);
    $self->parse_symbol($4, &_decode_JSON_string($5), &_decode_colour($6), $1);
  } elsif ($line =~ m/$CORE_SINGLE_RE/) {
    $self->parse_symbol($1, &_decode_JSON_string($2), &_decode_colour($3));
  } elsif ($line =~ m/$AMBIG_RE/) {
    $self->parse_symbol($1, &_decode_JSON_string($2), &_decode_colour($3), undef, $4);
  } else {
    $self->_error('unrecognised pattern');
  }
}

sub _add_lookup {
  my ($lookup, $target, $character, $both_case) = @_;
  $lookup->{$character} = $target;
  if ($both_case) {
    if ($character ge 'A' && $character le 'Z') {
      $lookup->{chr(ord($character) + 32)} = $target;
    } elsif ($character ge 'a' && $character le 'z') {
      $lookup->{chr(ord($character) - 32)} = $target;
    }
  }
}

sub parse_done {
  my $self = shift;
  die("Already finished parsing") if $self->{parsed};
  # now sort the symbols
  @{$self->{core_syms}} = sort _symobj_cmp @{$self->{core_syms}};
  @{$self->{ambig_syms}} = sort _symobj_cmp @{$self->{ambig_syms}};
  # check if there is a wildcard
  if (scalar(@{$self->{ambig_syms}}) == 0 || length($self->{ambig_syms}->[0]->{comprise}) != scalar(@{$self->{core_syms}})) {
    #warn("Warning: wildcard symbol automatically generated\n");
    my $wildcard_comprise = '';
    for (my $i = 0; $i < scalar(@{$self->{core_syms}}); $i++) {
      $wildcard_comprise .= $self->{core_syms}->[$i]->{sym};
    }
    unshift(@{$self->{ambig_syms}}, {sym => '?', aliases => [], comprise => $wildcard_comprise});
  }
  # create the final list of symbols and determine aliases
  my %comprise_lookup = ();
  $self->{ncore} = scalar(@{$self->{core_syms}});
  $self->{syms} = [];
  for (my $i = 0; $i < $self->{ncore}; $i++) {
    my $sym = $self->{core_syms}->[$i];
    $comprise_lookup{$sym->{sym}} = $sym;
    $sym->{index} = scalar(@{$self->{syms}});
    push(@{$self->{syms}}, $sym);
  }
  for (my $i = 0; $i < scalar(@{$self->{ambig_syms}}); $i++) {
    my $sym = $self->{ambig_syms}->[$i];
    my $real = $comprise_lookup{$sym->{comprise}};
    if (defined($real)) {
      push(@{$real->{aliases}}, $sym->{sym});
    } else {
      $comprise_lookup{$sym->{comprise}} = $sym;
      $sym->{index} = scalar(@{$self->{syms}});
      push(@{$self->{syms}}, $sym);
    }
  }
  # cleanup some entries we don't need
  delete $self->{sym_lookup};
  delete $self->{core_syms};
  delete $self->{ambig_syms};
  # map comprising characters to objects
  for (my $i = $self->{ncore}; $i < scalar(@{$self->{syms}}); $i++) {
    my $sym = $self->{syms}->[$i];
    my @list = map { $comprise_lookup{$_} } split(//, $sym->{comprise});
    $sym->{comprise} = \@list;
  }
  #determine complements
  my $fully_complementable = 1;
  for (my $i = 0; $i < $self->{ncore}; $i++) {
    my $sym = $self->{syms}->[$i];
    unless (defined($sym->{complement})) {
      $fully_complementable = 0;
      next;
    }
    $sym->{complement} = $comprise_lookup{$sym->{complement}};
  }
  AMBIG: for (my $i = $self->{ncore}; $i < scalar(@{$self->{syms}}); $i++) {
    my $sym = $self->{syms}->[$i];
    my @complements = ();
    for (my $j = 0; $j < scalar(@{$sym->{comprise}}); $j++) {
      my $target = $sym->{comprise}->[$j]->{complement};
      next AMBIG unless defined $target;
      push(@complements, $target->{sym});
    }
    @complements = sort _symbol_cmp @complements;
    my $complement = $comprise_lookup{join('', @complements)};
    next AMBIG unless defined $complement;
    $sym->{complement} = $complement;
  }
  # create lookup
  my $both_case = ($self->{seen_lc} != $self->{seen_uc});
  $self->{lookup} = {};
  for (my $i = 0; $i < scalar(@{$self->{syms}}); $i++) {
    my $sym = $self->{syms}->[$i];
    &_add_lookup($self->{lookup}, $sym, $sym->{sym}, $both_case);
    for (my $j = 0; $j < scalar(@{$sym->{aliases}}); $j++) {
      &_add_lookup($self->{lookup}, $sym, $sym->{aliases}->[$j], $both_case);
    }
  }
  $self->{comprise_lookup} = \%comprise_lookup;
  $self->{fully_complementable} = $fully_complementable;
  # create regular expressions to match something that isn't a valid symbol
  my $all_symbols = quotemeta(join('', (keys %{$self->{lookup}})));
  $self->{re_not_valid} = qr/[^$all_symbols]/;
  my $ambig_symbols = '';
  for (my $i = $self->{ncore}; $i < scalar(@{$self->{syms}}); $i++) {
    my $entry = $self->{syms}->[$i];
    my $prime_sym = $entry->{sym};
    $ambig_symbols .= $prime_sym;
    if ($both_case && $prime_sym =~ m/^[A-Za-z]$/) {
      my $prime_sym_copy = $prime_sym; # copy
      $prime_sym_copy =~ tr/A-Za-z/a-zA-Z/; # swap case
      $ambig_symbols .= $prime_sym_copy;
    }
    foreach my $alias_sym (@{$entry->{aliases}}) {
      $ambig_symbols .= $alias_sym;
      if ($both_case && $alias_sym =~ m/^[A-Za-z]$/) {
        my $alias_sym_copy = $alias_sym; # copy
        $alias_sym_copy =~ tr/A-Za-z/a-zA-Z/; # swap case
        $ambig_symbols .= $alias_sym_copy;
      }
    }
  }
  $self->{re_ambig} = qr/[\Q$ambig_symbols\E]/;
  $self->{re_ambig_only} = qr/^[\Q$ambig_symbols\E]*$/;
  #$ambig_symbols = quotemeta($ambig_symbols);
  #$self->{re_ambig} = qr/[$ambig_symbols]/;
  #$self->{re_ambig_only} = qr/^[$ambig_symbols]*$/;

  # create translatation strings
  # Strings to translate aliases of core symbols to prime representation
  my ($tr_src_core_alias_to_prime, $tr_dst_core_alias_to_prime) = ('', '');
  for (my $i = 0; $i < $self->{ncore}; $i++) {
    my $entry = $self->{syms}->[$i];
    my $prime_sym = $entry->{sym};
    # alternate case of primary symbol?
    if ($both_case && $prime_sym =~ m/^[A-Za-z]$/) {
      my $prime_sym_copy = $prime_sym; # copy
      $prime_sym_copy =~ tr/A-Za-z/a-zA-Z/; # swap case
      $tr_src_core_alias_to_prime .= $prime_sym_copy;
      $tr_dst_core_alias_to_prime .= $prime_sym;
    }
    # alias symbols
    foreach my $alias_sym (@{$entry->{aliases}}) {
      $tr_src_core_alias_to_prime .= $alias_sym;
      $tr_dst_core_alias_to_prime .= $prime_sym;
      # alternate case of alias symbol?
      if ($both_case && $alias_sym =~ m/^[A-Za-z]$/) {
        my $alias_sym_copy = $alias_sym; # copy
        $alias_sym_copy =~ tr/A-Za-z/a-zA-Z/; # swap case
        $tr_src_core_alias_to_prime .= $alias_sym_copy;
        $tr_dst_core_alias_to_prime .= $prime_sym;
      }
    }
  }
  $tr_src_core_alias_to_prime =~ s/([\/\-])/\\$1/g; # escape forward slash and dash
  $tr_dst_core_alias_to_prime =~ s/([\/\-])/\\$1/g; # escape forward slash and dash
  $self->{tr_src_core_alias_to_prime} = $tr_src_core_alias_to_prime;
  $self->{tr_dst_core_alias_to_prime} = $tr_dst_core_alias_to_prime;
  # Strings to translate aliases of ambiguous symbols to prime representation
  my ($tr_src_ambig_alias_to_prime, $tr_dst_ambig_alias_to_prime) = ('', '');
  for (my $i = $self->{ncore}; $i < scalar(@{$self->{syms}}); $i++) {
    my $entry = $self->{syms}->[$i];
    my $prime_sym = $entry->{sym};
    # alternate case of primary symbol?
    if ($both_case && $prime_sym =~ m/^[A-Za-z]$/) {
      my $prime_sym_copy = $prime_sym; # copy
      $prime_sym_copy =~ tr/A-Za-z/a-zA-Z/; # swap case
      $tr_src_ambig_alias_to_prime .= $prime_sym_copy;
      $tr_dst_ambig_alias_to_prime .= $prime_sym;
    }
    # alias symbols
    foreach my $alias_sym (@{$entry->{aliases}}) {
      $tr_src_ambig_alias_to_prime .= $alias_sym;
      $tr_dst_ambig_alias_to_prime .= $prime_sym;
      # alternate case of alias symbol?
      if ($both_case && $alias_sym =~ m/^[A-Za-z]$/) {
        my $alias_sym_copy = $alias_sym;
        $alias_sym_copy =~ tr/A-Za-z/a-zA-Z/;
        $tr_src_ambig_alias_to_prime .= $alias_sym_copy;
        $tr_dst_ambig_alias_to_prime .= $prime_sym;
      }
    }
  }
  $tr_src_ambig_alias_to_prime =~ s/([\/\-])/\\$1/g; # escape forward slash and dash
  $tr_dst_ambig_alias_to_prime =~ s/([\/\-])/\\$1/g; # escape forward slash and dash
  $self->{tr_src_ambig_alias_to_prime} = $tr_src_ambig_alias_to_prime;
  $self->{tr_dst_ambig_alias_to_prime} = $tr_dst_ambig_alias_to_prime;
  # Strings to translate ambiguous symbols and their aliases to the wildcard
  my ($tr_src_ambig_to_wild, $tr_dst_ambig_to_wild) = ('', '');
  my $wildcard = $self->{syms}->[$self->{ncore}]->{sym};
  for (my $i = $self->{ncore}; $i < scalar(@{$self->{syms}}); $i++) {
    my $entry = $self->{syms}->[$i];
    my $prime_sym = $entry->{sym};
    # primary symbol to wildcard (unless it is the wildcard)
    if ($i != $self->{ncore}) {
      $tr_src_ambig_to_wild .= $prime_sym;
      $tr_dst_ambig_to_wild .= $wildcard;
    }
    # alternate case of primary symbol?
    if ($both_case && $prime_sym =~ m/^[A-Za-z]$/) {
      my $prime_sym_copy = $prime_sym;
      $prime_sym_copy =~ tr/A-Za-z/a-zA-Z/;
      $tr_src_ambig_to_wild .= $prime_sym_copy;
      $tr_dst_ambig_to_wild .= $wildcard;
    }
    # alias symbols
    foreach my $alias_sym (@{$entry->{aliases}}) {
      $tr_src_ambig_to_wild .= $alias_sym;
      $tr_dst_ambig_to_wild .= $wildcard;
      # alternate case of alias symbol?
      if ($both_case && $alias_sym =~ m/^[A-Za-z]$/) {
        my $alias_sym_copy = $alias_sym;
        $alias_sym_copy =~ tr/A-Za-z/a-zA-Z/;
        $tr_src_ambig_to_wild .= $alias_sym_copy;
        $tr_dst_ambig_to_wild .= $wildcard;
      }
    }
  }
  $tr_src_ambig_to_wild =~ s/([\/\-])/\\$1/g; # escape forward slash and dash
  $tr_dst_ambig_to_wild =~ s/([\/\-])/\\$1/g; # escape forward slash and dash
  $self->{tr_src_ambig_to_wild} = $tr_src_ambig_to_wild;
  $self->{tr_dst_ambig_to_wild} = $tr_dst_ambig_to_wild;
  # Strings to translate unknown symbols to the wildcard
  my ($tr_src_unknown_to_wild, $tr_dst_unknown_to_wild) = ('', '');
  my $range_start = 0;
  for (my $i = 33; $i <= 126; $i++) {
    if (defined($self->{lookup}->{chr($i)})) {
      if ($i > ($range_start + 1)) {
        $tr_src_unknown_to_wild .= sprintf("\\%03o-\\%03o", $range_start, $i - 1);
      } elsif ($i == ($range_start + 1)) {
        $tr_src_unknown_to_wild .= sprintf("\\%03o", $range_start);
      }
      $range_start = $i + 1;
    }
  }
  $tr_src_unknown_to_wild .= sprintf("\\%03o-\\377", $range_start);
  $tr_dst_unknown_to_wild = ($wildcard eq '-' || $wildcard eq '/' ? "\\" : "") .$wildcard; # tr will extend using the last character
  $self->{tr_src_unknown_to_wild} = $tr_src_unknown_to_wild;
  $self->{tr_dst_unknown_to_wild} = $tr_dst_unknown_to_wild;
  # check 'like' a standard alphabet
  if ($self->{like}) {
    my ($req_core, $req_comp, $name);
    if ($self->{like} eq 'RNA') {
      $req_core = 'ACGU';
      $name = 'RNA';
    } elsif ($self->{like} eq 'DNA') {
      $req_core = 'ACGT';
      $req_comp = 'TGCA';
      $name = 'DNA';
    } elsif ($self->{like} eq 'PROTEIN') {
      $req_core = 'ACDEFGHIKLMNPQRSTVWY';
      $name = 'protein';
    }
    if (defined $req_core) {
      for (my $i = 0; $i < length($req_core); $i++) {
        my $sym = substr($req_core, $i, 1);
        my $symobj = $self->{lookup}->{$sym};
        unless (defined $symobj) {
          $self->_error("not $name like - expected symbol '$sym'");
        }
        unless ($symobj->{index} < $self->{ncore}) {
          $self->_error("not $name like - symbol '$sym' must be a core symbol");
        }
        if (defined $req_comp) {
          my $comp = substr($req_comp, $i, 1);
          my $compobj = $self->{lookup}->{$comp};
          if (not defined $compobj || not defined $symobj->{complement} || $symobj->{complement} != $compobj) {
            $self->_error("not $name like - symbol '$sym' complement must be '$comp'");
          }
        } else {
          #if (defined $symobj->{complement}) {
          #  $self->_error("not $name like - symbol '$sym' must not have complement");
          #}
        }
      }
    }
  }
  # done
  $self->{parsed} = 1; # TRUE
}

sub parse_file {
  my $self = shift;
  die("Already finished parsing") if $self->{parsed};
  my ($file) = @_;
  $self->{file_name} = $file;
  my ($fh, $line);
  sysopen($fh, $file, O_RDONLY);
  while($line = <$fh>) {
    chomp($line);
    $self->parse_line($line);
  }
  close($fh);
  $self->parse_done();
}

sub _sym_to_text {
  my ($sym) = @_;
  my $name_info = '';
  if (defined($sym->{name})) {
    $name_info = ' "' . &_encode_JSON_string($sym->{name}) . '"';
  }
  my $colour_info = '';
  if (defined($sym->{colour})) {
    $colour_info = sprintf(" %.06X", $sym->{colour});
  }
  return $sym->{sym} . $name_info . $colour_info;
}

sub to_text {
  my $self = shift;
  my %options = @_;
  my $out = '';
  my $gap = undef;
  # set option defaults
  $options{print_header} = 1 unless defined $options{print_header};
  $options{print_footer} = 1 unless defined $options{print_footer};
  # output header
  if ($options{print_header}) {
    $out .= "ALPHABET " . $self->name(1) . 
        ($self->{like} ? ' ' . $self->{like} . '-LIKE' : '') . "\n";
  }
  # output all paired core symbols
  for (my $i = 0; $i < $self->{ncore}; $i++) {
    my $sym = $self->{syms}->[$i];
    my $csym = $sym->{complement};
    if (defined($csym) && $sym->{index} < $csym->{index}) {
      $out .= &_sym_to_text($sym) . ' ~ ' . &_sym_to_text($csym) . "\n";
    }
  }
  # output all unpared core symbols
  for (my $i = 0; $i < $self->{ncore}; $i++) {
    my $sym = $self->{syms}->[$i];
    if (!defined($sym->{complement})) {
      $out .=  &_sym_to_text($sym) . "\n";
    }
  }
  # output ambiguous symbols
  for (my $i = $self->{ncore}; $i < scalar(@{$self->{syms}}); $i++) {
    my $sym = $self->{syms}->[$i];
    if (scalar(@{$sym->{comprise}}) == 0) {
      $gap = $sym;
      next;
    }
    my $comprise = '';
    for (my $j = 0; $j < scalar(@{$sym->{comprise}}); $j++) {
      $comprise .= $sym->{comprise}->[$j]->{sym};
    }
    $out .= &_sym_to_text($sym) . ' = ' . $comprise . "\n";
    for (my $j = 0; $j < scalar(@{$sym->{aliases}}); $j++) {
      $out .= $sym->{aliases}->[$j] . ' = ' . $comprise . "\n";
    }
  }
  # output core symbol aliases
  for (my $i = 0; $i < $self->{ncore}; $i++) {
    my $sym = $self->{syms}->[$i];
    for (my $j = 0; $j < scalar(@{$sym->{aliases}}); $j++) {
      $out .= $sym->{aliases}->[$j] . ' = ' . $sym->{sym} . "\n";
    }
  }
  # output gap symbol and aliases
  if ($gap) {
    $out .= &_sym_to_text($gap) . " =\n"; 
    for (my $j = 0; $j < scalar(@{$gap->{aliases}}); $j++) {
      $out .= $gap->{aliases}->[$j] . " =\n";
    }
  }
  # output footer
  if ($options{print_footer}) {
    $out .= "END ALPHABET\n";
  }
  return $out;
}

sub to_json() {
  my $self = shift;
  my ($jw) = @_;
  $jw->start_object_value();
  $jw->str_prop("name", $self->{name});
  $jw->str_prop("like", lc($self->{like})) if ($self->{like});
  $jw->num_prop("ncore", $self->{ncore});
  $jw->property("symbols");
  $jw->start_array_value();
  for (my $i = 0; $i < scalar(@{$self->{syms}}); $i++) {
    my $sym = $self->{syms}->[$i];
    $jw->start_object_value();
    $jw->str_prop("symbol", $sym->{sym});
    $jw->str_prop("aliases", join("", @{$sym->{aliases}})) if (scalar(@{$sym->{aliases}}) > 0);
    $jw->str_prop("name", $sym->{name}) if (defined($sym->{name}));
    $jw->str_prop("colour", sprintf("%06X", $sym->{colour})) if (defined($sym->{colour}) && $sym->{colour} > 0);
    if ($i < $self->{ncore}) {
      $jw->str_prop("complement", $sym->{complement}->{sym}) if (defined($sym->{complement}));
    } else {
      my $comprise = '';
      for (my $j = 0; $j < scalar(@{$sym->{comprise}}); $j++) {
        $comprise .= $sym->{comprise}->[$j]->{sym};
      }
      $jw->str_prop("equals", $comprise);
    }
    $jw->end_object_value();
  }
  $jw->end_array_value();
  $jw->end_object_value();
}

sub has_errors {
  my $self = shift;
  return scalar(@{$self->{errors}});
}

sub get_errors {
  my $self = shift;
  return @{$self->{errors}};
}

sub name {
  my $self = shift;
  my ($should_encode) = @_;
  if ($should_encode) {
    return '"' . &_encode_JSON_string($self->{name}) . '"';
  } else {
    return $self->{name};
  }
}

sub has_complement {
  my $self = shift;
  die("Not finished parsing") unless $self->{parsed};
  return $self->{fully_complementable};
}

sub size_core {
  my $self = shift;
  return $self->{ncore};
}

sub size_full {
  my $self = shift;
  return scalar(@{$self->{syms}});
}

#
# Returns a scalar, which is the concatenation of all the
# possible ambiguity codes of the alphabet.
#
sub get_ambig_syms {
    my $self = shift;
    my $ambig_syms;
    foreach my $sym (split(//, $self->get_syms())) {
        $ambig_syms .= $sym unless ($self->is_core($sym));
    }
    return $ambig_syms;
}

#
# Returns a scalar, which is the concatenation of all the
# possible symbols of the alphabet.
#
sub get_syms {
    my $self = shift;
    my $res = '';
    foreach my $symbol (@{$self->{syms}}) {
        $res .= $symbol->{sym};
    }
    return $res;
}

#
# Returns a scalar, which is the concatenation of all the
# core symbols of the alphabet.
#
sub get_core {
    my $self = shift;
    my $core_syms;
    foreach my $sym (split(//, $self->get_syms())) {
        $core_syms .= $sym if ($self->is_core($sym));
    }
    return $core_syms;
}

#
# Returns a dictionary containing the core symbols that
# each ambiguity code of the alphabet is composed of.
#
sub get_alternates_dict {
    my $self = shift;
    my %alternates;
    foreach my $ambig_base (split(//, $self->get_ambig_syms())) {
        %alternates = (%alternates, %{$self->comprise($ambig_base)});
    }
    return %alternates;
}

#
# Is the symbol part of the core alphabet.
#
sub is_core {
  my $self = shift;
  my ($letter) = @_;
  my $sym = $self->{lookup}->{$letter};
  return defined($sym) && $sym->{index} < $self->{ncore};
}

#
# Is the symbol part of the alphabet.
#
sub is_known {
  my $self = shift;
  my ($letter) = @_;
  my $sym = $self->{lookup}->{$letter};
  return defined($sym);
}

#
# Does the sequence only contain ambiguous symbols.
# This might be useful if you want to throw away something containing only
# ambiguous symbols as that has no real useful information.
#
sub is_ambig_seq {
  my $self = shift;
  my ($seq) = @_;
  my $ambig_only_re = $self->{re_ambig_only};
  return scalar($seq =~ m/$ambig_only_re/);
}

#
# Does the sequence contain ambiguous symbols.
#
sub is_part_ambig_seq {
  my $self = shift;
  my ($seq) = @_;
  my $ambig_re = $self->{re_ambig};
  return scalar($seq =~ m/$ambig_re/);
}

sub char {
  my $self = shift;
  my ($index) = @_;
  return undef if ($index < 0 || $index > scalar(@{$self->{syms}}));
  return $self->{syms}->[$index]->{sym};
}

#
# Converts a symbol into the index in the alphabet.
#
sub index {
  my $self = shift;
  my ($letter) = @_;
  my $sym = $self->{lookup}->{$letter};
  return undef unless defined $sym;
  return $sym->{index};
}

#
# Converts a symbol that could be an alias into the primary symbol.
#
sub encode {
  my $self = shift;
  my ($letter) = @_;
  my $sym = $self->{lookup}->{$letter};
  return undef unless defined $sym;
  return $sym->{sym};
}

#
# Takes a letter and returns the set of comprising core symbols.
#
sub comprise {
  my $self = shift;
  my ($letter) = @_;
  my $sym = $self->{lookup}->{$letter};
  my %set = ();
  if (defined($sym->{comprise})) {
    for (my $i = 0; $i < scalar(@{$sym->{comprise}}); $i++) {
      $set{$sym->{comprise}->[$i]->{sym}} = 1;
    }
  } else {
    $set{$sym->{sym}} = 1;
  }
  return \%set;
}

#
# Converts a set of core symbols into either a matching ambigous symbol or
# a group of core symbols.
#
sub regex {
  my $self = shift;
  my ($set) = @_;
  my $key = join('', sort _symbol_cmp (keys %{$set}));
  my $sym = $self->{comprise_lookup}->{$key};
  if (defined($sym)) {
    return $sym->{sym};
  } else {
    return '['.$key.']';
  }
}

#
# Takes a set of core symbols and creates a set containing the other core symbols.
#
sub negate_set {
  my $self = shift;
  my ($set) = @_;
  my %negated_set = ();
  for (my $i = 0; $i < $self->{ncore}; $i++) {
    my $a = $self->{syms}->[$i]->{sym};
    if (!($set->{$a})) {
      $negated_set{$a} = 1;
    }
  }
  return \%negated_set;
}

#
# Finds the first unknown symbol.
#
sub find_unknown {
  my $self = shift;
  my ($seq) = @_;
  my $not_valid = $self->{re_not_valid};
  if ($seq =~ m/($not_valid)/) {
    return ($-[1], $1);
  } else {
    return undef;
  }
}

#
# Reverse complements a sequence using the alphabet.
#
sub rc_seq {
  my $self = shift;
  croak("Alphabet does not allow complementing.") unless $self->has_complement();
  my ($sequence) = @_;
  my @seq = split //, $sequence;
  @seq = reverse @seq;
  @seq = map { $self->{lookup}->{$_}->{complement}->{sym} } @seq;
  return join('', @seq);
}

#
# Translate a sequence using the alphabet.
#
# NO_ALIASES    Convert any alias symbols into the main symbol.
# NO_AMBIGS     Convert any ambiguous symbols into the main wildcard symbol.
# NO_UNKNOWN    Convert any unknown symbols into the main wildcard symbol.
#
sub translate_seq {
  my $self = shift;
  my ($sequence, %options) = @_;
  my $no_aliases = (defined($options{NO_ALIASES}) ? $options{NO_ALIASES} : 0);
  my $no_ambigs = (defined($options{NO_AMBIGS}) ? $options{NO_AMBIGS} : 0);
  my $no_unknown = (defined($options{NO_UNKNOWN}) ? $options{NO_UNKNOWN} : 0);
  my ($src, $dst) = ('', '');
  if ($no_aliases) {
    $src .= $self->{tr_src_core_alias_to_prime};
    $dst .= $self->{tr_dst_core_alias_to_prime};
    if ($no_ambigs) {
      $src .= $self->{tr_src_ambig_to_wild};
      $dst .= $self->{tr_dst_ambig_to_wild};
    } else {
      $src .= $self->{tr_src_ambig_alias_to_prime};
      $dst .= $self->{tr_dst_ambig_alias_to_prime};
    }
  } elsif ($no_ambigs) {
    $src .= $self->{tr_src_ambig_to_wild};
    $dst .= $self->{tr_dst_ambig_to_wild};
  }
  if ($no_unknown) {
    $src .= $self->{tr_src_unknown_to_wild};
    $dst .= $self->{tr_dst_unknown_to_wild};
  }
  eval "\$sequence =~ tr/$src/$dst/";
  return $sequence;
}

#
# Compare 2 different alphabet objects for equality
#
sub equals {
  my $self = shift;
  my ($other) = @_;
  return 0 if ($self->size_core() != $other->size_core());
  return 0 if ($self->size_full() != $other->size_full());
  return 0 if ($self->has_complement() != $other->has_complement());
  return 0 if (defined($self->{name}) != defined($other->{name}));
  if (defined($self->{name}) && defined($other->{name})) {
    return 0 if ($self->{name} ne $other->{name});
  }
  my $nsyms = $self->size_full();
  for (my $i = 0; $i < $nsyms; $i++) {
    my $self_sym = $self->{syms}->[$i];
    my $other_sym = $other->{syms}->[$i];
    return 0 if ($self_sym->{sym} ne $other_sym->{sym});
    return 0 if (defined($self_sym->{name}) != defined($other_sym->{name}));
    if (defined($self_sym->{name}) && defined($other_sym->{name})) {
      return 0 if ($self_sym->{name} ne $other_sym->{name});
    }
    return 0 if (defined($self_sym->{colour}) != defined($other_sym->{colour}));
    if (defined($self_sym->{colour}) && defined($other_sym->{colour})) {
      return 0 if ($self_sym->{colour} != $other_sym->{colour});
    }
    return 0 if (scalar(@{$self_sym->{aliases}}) != scalar(@{$other_sym->{aliases}}));
    my $naliases = scalar(@{$self_sym->{aliases}});
    for (my $j = 0; $j < $naliases; $j++) {
      return 0 if ($self_sym->{aliases}->[$j] ne $other_sym->{aliases}->[$j]);
    }
  }
  $nsyms = $self->size_core();
  for (my $i = 0; $i < $nsyms; $i++) {
    my $self_sym = $self->{syms}->[$i];
    my $other_sym = $other->{syms}->[$i];
    return 0 if (defined($self_sym->{complement}) != defined($other_sym->{complement}));
    if (defined($self_sym->{complement}) && defined($other_sym->{complement})) {
      return 0 if ($self_sym->{complement}->{sym} ne $other_sym->{complement}->{sym});
    }
  }
  $nsyms = $self->size_full();
  for (my $i = $self->size_core(); $i < $nsyms; $i++) {
    my $self_sym = $self->{syms}->[$i];
    my $other_sym = $other->{syms}->[$i];
    return 0 if (scalar(@{$self_sym->{comprise}}) != scalar(@{$other_sym->{comprise}}));
    for (my $j = 0; $j < scalar(@{$self_sym->{comprise}}); $j++) {
      return 0 if ($self_sym->{comprise}->[$j]->{sym} ne $other_sym->{comprise}->[$j]->{sym});
    }
  }
  return 1;
}

sub is_dna {
  my $self = shift;
  $self->{cached_is_dna} = $self->equals(&dna()) unless defined $self->{cached_is_dna};
  return $self->{cached_is_dna};
}

sub is_rna {
  my $self = shift;
  $self->{cached_is_rna} = $self->equals(&rna()) unless defined $self->{cached_is_rna};
  return $self->{cached_is_rna};
}

sub is_protein {
  my $self = shift;
  $self->{cached_is_protein} = $self->equals(&protein()) unless defined $self->{cached_is_protein};
  return $self->{cached_is_protein};
}

#
# Helper function to decode a string encoded using the method defined in the JSON specification.
#
sub _decode_JSON_string {
  my ($text) = @_;
  return undef unless defined $text;
  $text =~ s/\\(["\\\/])/$1/g;
  $text =~ s/\\b/\x{8}/g;
  $text =~ s/\\f/\x{C}/g;
  $text =~ s/\\n/\x{A}/g;
  $text =~ s/\\r/\x{D}/g;
  $text =~ s/\\t/\x{9}/g;
  $text =~ s/\\u([0-9A-Fa-f]{4})/chr(hex($1))/ge;
  return $text;
}

#
# Helper function to encode a string using the method defined in the JSON specification.
#
sub _encode_JSON_string {
  my ($text) = @_;
  $text =~ s/(["\\\/])/\\$1/g;
  $text =~ s/\x{8}/\\b/g;
  $text =~ s/\x{C}/\\f/g;
  $text =~ s/\x{A}/\\n/g;
  $text =~ s/\x{D}/\\r/g;
  $text =~ s/\x{9}/\\t/g;
  $text =~ s/([\x{0}-\x{1F}\x{7F}-\x{9F}\x{2028}\x{2029}])/sprintf("\\u%.04X",ord($1))/ge;
  return $text;
}

#
# Helper function to convert the colour into a number.
#
sub _decode_colour {
  my ($text) = @_;
  return undef unless defined $text;
  return hex($text);
}

# _symbol_cmp
# Sorts so letters are before numbers which are before symbols.
# Otherwise sorts using normal string comparison.
sub _symbol_cmp($$) {
  my ($a, $b) = @_;
  if (length($a) == length($b)) {
    for (my $i = 0; $i < length($a); $i++) {
      my ($sym_a, $sym_b);
      $sym_a = ord(substr($a, $i, 1));
      $sym_b = ord(substr($b, $i, 1));
      # check to see if either is a letter
      my ($is_letter_a, $is_letter_b);
      $is_letter_a = (($sym_a >= ord('A') && $sym_a <= ord('Z')) || ($sym_a >= ord('a') && $sym_a <= ord('z')));
      $is_letter_b = (($sym_b >= ord('A') && $sym_b <= ord('Z')) || ($sym_b >= ord('a') && $sym_b <= ord('z')));
      if ($is_letter_a) {
        if ($is_letter_b) {
          if ($sym_a < $sym_b) {
            return -1;
          } elsif ($sym_b < $sym_a) {
            return 1;
          }
          next;
        } else {
          return -1;
        }
      } elsif ($is_letter_b) {
        return 1;
      }
      # check to see if either is a number
      my ($is_num_a, $is_num_b);
      $is_num_a = ($sym_a ge '0' && $sym_a le '9');
      $is_num_b = ($sym_b ge '0' && $sym_b le '9');
      if ($is_num_a) {
        if ($is_num_b) {
          if ($sym_a < $sym_b) {
            return -1;
          } elsif ($sym_b < $sym_a) {
            return 1;
          }
          next;
        } else {
          return -1;
        }
      } elsif ($is_num_b) {
        return 1;
      }
      # both must be symbols
      if ($sym_a < $sym_b) {
        return -1;
      } elsif ($sym_b < $sym_a) {
        return 1;
      }
    }
    return 0;
  } elsif (length($a) > length($b)) {
    return -1;
  } else {
    return 1;
  }
}

sub _symobj_cmp($$) {
  my ($sym1, $sym2) = @_;
  my ($ret);
  # compare comprise
  if (defined($sym1->{comprise}) && defined($sym2->{comprise})) {
    $ret = &_symbol_cmp($sym1->{comprise}, $sym2->{comprise});
    return $ret if $ret != 0;
  } elsif (defined($sym1->{comprise})) {
    return 1; # sym2 first because it is a core symbol
  } elsif (defined($sym2->{comprise})) {
    return -1; # sym1 first because it is a core symbol
  }
  # compare symbol
  $ret = &_symbol_cmp($sym1->{sym}, $sym2->{sym});
  return $ret if $ret != 0;
  # compare complement
  if (defined($sym1->{complement}) && defined($sym2->{complement})) {
    $ret = &_symbol_cmp($sym1->{complement}, $sym2->{complement});
    return $ret if $ret != 0;
  } elsif (defined($sym1->{complement})) {
    return 1;
  } elsif (defined($sym2->{complement})) {
    return -1;
  }
  # compare aliases
  if (scalar(@{$sym1->{aliases}}) != scalar(@{$sym2->{aliases}})) {
    # sort length decending
    return scalar(@{$sym2->{aliases}}) <=> scalar(@{$sym1->{aliases}});
  }
  for (my $i = 0; $i < scalar(@{$sym1->{aliases}}); $i++) {
    $ret = &_symbol_cmp($sym1->{aliases}->[$i], $sym2->{aliases}->[$i]);
    return $ret if $ret != 0;
  }
  # compare name
  if (defined($sym1->{name}) && defined($sym2->{name})) {
    $ret = $sym1->{name} cmp $sym2->{name};
    return $ret if $ret != 0;
  } elsif (defined($sym1->{name})) {
    return 1;
  } elsif (defined($sym2->{name})) {
    return -1;
  }
  # compare colour
  if (defined($sym1->{colour}) && defined($sym2->{colour})) {
    return $sym1->{colour} <=> $sym2->{colour};
  } elsif (defined($sym1->{colour})) {
    return 1;
  } elsif (defined($sym2->{colour})) {
    return -1;
  }
  return 0;
}
#
# Return the RNA alphabet.
#
sub rna {
  my $alph = new Alphabet();
  $alph->parse_header('RNA', 'RNA');
  # core symbols
  $alph->parse_symbol('A', 'Adenine', 0xCC0000);
  $alph->parse_symbol('C', 'Cytosine', 0x0000CC);
  $alph->parse_symbol('G', 'Guanine', 0xFFB300);
  $alph->parse_symbol('U', 'Uracil', 0x008000, undef, undef, 'T');
  # ambiguous symbols
  $alph->parse_symbol('W', 'Weak', undef, undef, 'AU');
  $alph->parse_symbol('S', 'Strong', undef, undef, 'CG');
  $alph->parse_symbol('M', 'Amino', undef, undef, 'AC');
  $alph->parse_symbol('K', 'Keto', undef, undef, 'GU');
  $alph->parse_symbol('R', 'Purine', undef, undef, 'AG');
  $alph->parse_symbol('Y', 'Pyrimidine', undef, undef, 'CU');
  $alph->parse_symbol('B', 'Not A', undef, undef, 'CGU');
  $alph->parse_symbol('D', 'Not C', undef, undef, 'AGU');
  $alph->parse_symbol('H', 'Not G', undef, undef, 'ACU');
  $alph->parse_symbol('V', 'Not U', undef, undef, 'ACG');
  $alph->parse_symbol('N', 'Any base', undef, undef, 'ACGU', 'X.');
  # for now treat '.' (gap) as a wildcard
  # process
  $alph->parse_done();
  # error check
  die("Uhoh, this shouldn't happen:\n" + join("\n", $alph->get_errors())) if ($alph->has_errors());
  # return
  return $alph;
}

#
# Return the DNA alphabet.
#
sub dna {
  my $alph = new Alphabet();
  $alph->parse_header('DNA', 'DNA');
  # core symbols
  $alph->parse_symbol('A', 'Adenine', 0xCC0000, 'T');
  $alph->parse_symbol('C', 'Cytosine', 0x0000CC, 'G');
  $alph->parse_symbol('G', 'Guanine', 0xFFB300, 'C');
  $alph->parse_symbol('T', 'Thymine', 0x008000, 'A', undef, 'U');
  # ambiguous symbols
  $alph->parse_symbol('W', 'Weak', undef, undef, 'AT');
  $alph->parse_symbol('S', 'Strong', undef, undef, 'CG');
  $alph->parse_symbol('M', 'Amino', undef, undef, 'AC');
  $alph->parse_symbol('K', 'Keto', undef, undef, 'GT');
  $alph->parse_symbol('R', 'Purine', undef, undef, 'AG');
  $alph->parse_symbol('Y', 'Pyrimidine', undef, undef, 'CT');
  $alph->parse_symbol('B', 'Not A', undef, undef, 'CGT');
  $alph->parse_symbol('D', 'Not C', undef, undef, 'AGT');
  $alph->parse_symbol('H', 'Not G', undef, undef, 'ACT');
  $alph->parse_symbol('V', 'Not T', undef, undef, 'ACG');
  $alph->parse_symbol('N', 'Any base', undef, undef, 'ACGT', 'X.');
  # for now treat '.' (gap) as a wildcard
  # process
  $alph->parse_done();
  die("Uhoh, this shouldn't happen:\n" + join("\n", $alph->get_errors())) if ($alph->has_errors());
  # return
  return $alph;
}

#
# Return the Protein alphabet.
#
sub protein {
  my $alph = new Alphabet();
  $alph->parse_header('Protein', 'PROTEIN');
  # core symbols
  $alph->parse_symbol('A', 'Alanine', 0x0000CC);
  $alph->parse_symbol('R', 'Arginine', 0xCC0000);
  $alph->parse_symbol('N', 'Asparagine', 0x008000);
  $alph->parse_symbol('D', 'Aspartic acid', 0xFF00FF);
  $alph->parse_symbol('C', 'Cysteine', 0x0000CC);
  $alph->parse_symbol('E', 'Glutamic acid', 0xFF00FF);
  $alph->parse_symbol('Q', 'Glutamine', 0x008000);
  $alph->parse_symbol('G', 'Glycine', 0xFFB300);
  $alph->parse_symbol('H', 'Histidine', 0xFFCCCC);
  $alph->parse_symbol('I', 'Isoleucine', 0x0000CC);
  $alph->parse_symbol('L', 'Leucine', 0x0000CC);
  $alph->parse_symbol('K', 'Lysine', 0xCC0000);
  $alph->parse_symbol('M', 'Methionine', 0x0000CC);
  $alph->parse_symbol('F', 'Phenylalanine', 0x0000CC);
  $alph->parse_symbol('P', 'Proline', 0xFFFF00);
  $alph->parse_symbol('S', 'Serine', 0x008000);
  $alph->parse_symbol('T', 'Threonine', 0x008000);
  $alph->parse_symbol('W', 'Tryptophan', 0x0000CC);
  $alph->parse_symbol('Y', 'Tyrosine', 0x33E6CC);
  $alph->parse_symbol('V', 'Valine', 0x0000CC);
  # ambiguous symbols
  $alph->parse_symbol('B', 'Asparagine or Aspartic acid', undef, undef, 'ND');
  $alph->parse_symbol('Z', 'Glutamine or Glutamic acid', undef, undef, 'QE');
  $alph->parse_symbol('J', 'Leucine or Isoleucine', undef, undef, 'LI');
  $alph->parse_symbol('X', 'Any amino acid', undef, undef, 'ARNDCEQGHILKMFPSTWYV', '*.');
  # for now treat '*' (stop codon) and '.' (gap) as a wildcard
  # process
  $alph->parse_done();
  # error check
  die("Uhoh, this shouldn't happen:\n" + join("\n", $alph->get_errors())) if ($alph->has_errors());
  # return
  return $alph;
}
