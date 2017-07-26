package CheckingSAX;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw($RPT_ONE $RPT_OPT $RPT_ALO $RPT_ANY _value_pvalue _value_evalue _value_yn _value_yesno _value_colour _struc_command_line _struc_alphabet);
our @EXPORT_OK = qw();

use Data::Dumper;
use XML::Parser::Expat;
use Carp;
use Alphabet;

our $RPT_ONE = 0; # occur once only (default)
our $RPT_OPT = 1; # optionaly occur once
our $RPT_ALO = 2; # repeat one or more times
our $RPT_ANY = 3; # repeat any number of times

my $SAX_DEBUG = !!($ENV{'SAX_DEBUG'});

sub new {
  my $classname = shift;
  my $self = {};
  bless($self, $classname);
  $self->_init(@_);
  return $self;
}

sub set_handlers {
  my $self = shift;
  my %handlers = @_;
  my ($key, $value);
  while (($key, $value) = each(%handlers)) {
    $self->{handlers}->{$key} = $value;
  }
}

sub parse_more {
  my $self = shift;
  my ($buffer) = @_;
  if (@{$self->{errors}}) {
    die("Unable to parse more as there are unhandled errors:\n" . join("\n", @{$self->{errors}}));
  }
  $self->{parser}->parse_more($buffer);
}

sub parse_done {
  my $self = shift;
  if (@{$self->{errors}}) {
    die("Unable to finish parsing as there are unhandled errors:\n" . join("\n", @{$self->{errors}}));
  }
  $self->{parser}->parse_done();
}

sub has_errors {
  my $self = shift;
  return scalar(@{$self->{errors}});
}

sub get_errors {
  my $self = shift;
  return @{$self->{errors}};
}

# sub classes should override this method
sub _xml_def {
  return [];
}

sub _init {
  my $self = shift;
  my ($userdata, %handlers) = @_;

  my $parser = new XML::Parser::ExpatNB; # non blocking parser
  $parser->setHandlers(Start => sub {$self->_handle_start(@_);});
  $parser->setHandlers(End => sub {$self->_handle_end(@_);});
  $parser->setHandlers(Char => sub {$self->_handle_char(@_);});
  $self->{parser} = $parser;
  $self->{userdata} = $userdata;
  $self->{handlers} = {};
  $self->{expected} = [];
  $self->{errors} = [];
  $self->{text} = undef;
  $self->{quiet} = 0;
  $self->_expand($self->_xml_def());
  $self->set_handlers(%handlers);
}

sub _expand {
  my $self = shift;
  my ($state) = @_;
  my @substates = @{$state};
  for (my $i = scalar(@substates) - 1; $i >= 0; $i--) {
    my $expect = {%{$substates[$i]}, SEEN => 0};
    $expect->{RPT} = 0 unless(defined($expect->{RPT}));
    push(@{$self->{expected}}, $expect);
  }
}

sub _find {
  my $self = shift;
  my ($element) = @_;
  my $top;
  # make sure the element is expected
  while (1) {
    $top = ${$self->{expected}}[-1];
    die("Unexpected element " . $element . "\n") unless (defined($top));
    # break out of loop if we've found the element
    last if ($top->{ELE} eq $element);
    # make sure this element can be discarded
    unless ($top->{SEEN} > 0 || $top->{RPT} == $RPT_OPT || $top->{RPT} == $RPT_ANY) {
      die("Missing expected element " . $top->{ELE} . ". Found " . $element . " instead\n"); 
    }
    pop(@{$self->{expected}});
  }
  return $top;
}

sub _handle_start {
  my $self = shift;
  my ($expat, $element, %attributes) = @_;

  my $top = $self->_find($element);
  # check that we can have more of this element
  unless ($top->{SEEN} == 0 || $top->{RPT} == $RPT_ALO || $top->{RPT} == $RPT_ANY) {
    die("Too many of element " . $top->{ELE} . "\n");
  }
  # increment the sighting count
  $top->{SEEN} += 1;
  # process the attributes
  my @params = ();
  if (defined($top->{ATRS})) {
    # process the attributes
    my @expected = @{$top->{ATRS}};
    PROC_ATTRS: for (my $i = 0; $i < scalar(@expected); $i++) {
      my $expect = $expected[$i];
      if (ref($expect) eq 'CODE') {
        splice(@expected, $i, 1, $expect->($self));
        redo PROC_ATTRS;
      }
      my $attribute = $expect->{ATR};
      my $count = (defined($expect->{TUPLE}) ? $expect->{TUPLE} : $expect->{EXPAND});
      $count = 1 unless (defined $count);
      my $value = $attributes{$attribute};
      if (defined($value)) {
        my $handler = (defined($expect->{VAL}) ? $expect->{VAL} : $expect->{OPT});
        if (defined($handler)) {
          if (ref($handler) eq 'Regexp') {
            unless ($value =~ $handler) {
              $self->_error("$element\@$attribute has invalid value \"$value\"");
            } elsif ($count == 1) {
              $value = $1;    
            } elsif ($count == 2) {
              $value = [$1, $2];
            } elsif ($count == 3) {
              $value = [$1, $2, $3];
            } elsif ($count == 4) {
              $value = [$1, $2, $3, $4];
            } elsif ($count == 5) {
              $value = [$1, $2, $3, $4, $5];
            } elsif ($count == 6) {
              $value = [$1, $2, $3, $4, $5, $6];
            } elsif ($count == 7) {
              $value = [$1, $2, $3, $4, $5, $6, $7];
            } elsif ($count == 8) {
              $value = [$1, $2, $3, $4, $5, $6, $7, $8];
            } elsif ($count == 9) {
              $value = [$1, $2, $3, $4, $5, $6, $7, $8, $9];
            } else {
              die("Value count larger than 9 can not be handled by a regular expression");
            }
          } elsif (ref($handler) eq 'CODE') {
            eval { $value = $handler->($self, $value, $count); };
            $self->_error("$element\@$attribute has invalid value \"$value\" - $@") if ($@);
            if ($count > 1) {
              die("Expected array reference, got " .ref($value)) unless (ref($value) eq 'ARRAY');
              die("Incorrect number of values returned") unless (scalar(@{$value}) == $count);
            }
          } else {
            warn "Unable to use handler type " . ref($handler) . "\n";
          }
        } elsif ($count != 1) {
          die("Expected a handler for a non-unary count");
        }
      } else {
        $self->_error("$element\@$attribute missing") if (defined($expect->{VAL}));
        $value = [(undef) x $count] if ($count > 1);
      }
      unless ($expect->{SILENT}) {
        if (defined($expect->{EXPAND})) {
          push(@params, @{$value});
        } else {
          push(@params, $value);
        }
      }
    }
  }
  if (defined($top->{START})) {
    @params = $top->{START}($self, \%attributes, @params);
  } elsif (defined($top->{WATCH})) {
    $top->{WATCH}($self, \%attributes, @params);
  }
  $top->{PARAMS} = \@params;
  if ($top->{HANDLE} || $self->{quiet}) {
    $self->{quiet} += 1;
  }
  unless ($top->{SILENT} || $self->{quiet}) {
    if (defined($top->{TO})) {
      #has sub states so need to call 'start_'
      my $call_name = 'start_' . (defined($top->{NAME}) ? $top->{NAME} : $element);
      $self->_call($call_name, @params);
    } else {
      # do the callback on end prefixed with the name 'handle_'
    }
  }
  # expand states
  $self->_expand($top->{TO}) if (defined($top->{TO}));
  $self->{text} = (defined($top->{VAL}) ? '' : undef);
}

sub _handle_end {
  my $self = shift;
  my ($expat, $element) = @_;
  my $top = $self->_find($element);
  # process value
  my $value = $self->{text};
  $self->{text} = undef;
  my $handler = $top->{VAL};
  if (defined($handler)) {
    if (ref($handler) eq 'Regexp') {
      unless ($value =~ $handler) {
        $self->_error("$element has invalid value \"$value\"");
      } else {
        $value = $1;    
      }
    } elsif (ref($handler) eq 'CODE') {
      eval { $value = $handler->($self, $value); };
      $self->_error("$element has invalid value \"$value\" - $@") if ($@);
    } else {
      warn "Unable to use handler type " . ref($handler) . "\n";
    }
  }
  my @value_param = (defined($value) ? ($value) : ());

  my @params;
  my $processor = (defined($top->{HANDLE}) ? $top->{HANDLE} : $top->{END});
  if (defined($processor)) {
    @params = $processor->($self, @{$top->{PARAMS}}, @value_param);
  } elsif (defined($top->{TO})) {
    @params = @value_param;
  } else {
    @params = (@{$top->{PARAMS}}, @value_param)
  }

  $self->{quiet} -= 1 if ($self->{quiet});

  unless ($top->{SILENT} || $self->{quiet}) {
    my $name = (defined($top->{NAME}) ? $top->{NAME} : $element);
    if (defined($top->{TO}) && !defined($top->{HANDLE})) {
      #has sub states so need to call 'end_'
      $self->_call('end_' . $name, @params);
    } else {
      # call handle
      $self->_call('handle_' . $name, @params);
    }
  }

}

sub _handle_char {
  my $self = shift;
  my ($expat, $string) = @_;
  $self->{text} .= $string if defined $self->{text};
}

sub _error {
  my $self = shift;
  my ($message) = @_;
  push(@{$self->{errors}}, $message);
}

sub _call {
  my $self = shift;
  my ($call, @args) = @_;
  if ($SAX_DEBUG) {
    print "CALL $_[0](";
    {
      local $Data::Dumper::Terse = 1;
      local $Data::Dumper::Indent = 0;
      print Dumper(@_[1 .. $#_]);
    }
    print ")\n";
  }
  if (defined($self->{handlers}->{$call}) && !$self->has_errors()) {
    $self->{handlers}->{$call}($self->{userdata}, @args);
  }
}

my $attribute_re = qr/^([a-zA-Z_:][-a-zA-Z0-9_:.]*)$/;
my $any_multiline_re = qr/^([\S\s]*)$/;
my $text_re = qr/^(.*\S.*)$/;
my $sci_re = qr/[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?/;
my $sci_or_inf_re = qr/(?:$sci_re|inf)/;
my $num_re = qr/^([+-]?$sci_or_inf_re)$/;
my $allowed_sym_re = qr/[A-Za-z0-9\.\*\?\-]/;
my $sym_re = qr/^($allowed_sym_re)$/;
my $syms_re = qr/^($allowed_sym_re+)$/;

#
# Check value is a p-value.
#
sub _value_pvalue {
  my ($self, $value) = @_;
  confess("No parameters provided! Check that you used a function reference and didn't call it by mistake!") unless @_;
  confess("Value is undefined") unless defined $value;
  if ($value =~ $num_re) {
    my $num = $1 + 0;
    if ($num < 0 || $num > 1) {
      die("not in range 0 to 1\n");
    }
    return $1;
  }
  die("not a number\n");
}

#
# Check value is an E-value.
#
sub _value_evalue {
  my ($self, $value) = @_;
  confess("No parameters provided! Check that you used a function reference and didn't call it by mistake!") unless @_;
  confess("Value is undefined") unless defined $value;
  if ($value =~ $num_re) {
    my $num = $1 + 0;
    if ($num < 0) {
      die("should not be negative\n");
    }
    return $1;
  }
  die("not a number\n");
}

#
# Check value is either y or n.
# Convert y to boolean true (1) and n to boolean false (0) before returning.
#
sub _value_yn {
  my ($self, $value) = @_;
  confess("No parameters provided! Check that you used a function reference and didn't call it by mistake!") unless @_;
  confess("Value is undefined") unless defined $value;
  if ($value eq 'y') {
    return 1;
  } elsif ($value eq 'n') {
    return 0;
  }
  die("expected y or n\n");
}

#
# Check value is either yes or no.
# Convert yes to boolean true (1) and no to boolean false (0) before returning.
#
sub _value_yesno {
  my ($self, $value) = @_;
  confess("No parameters provided! Check that you used a function reference and didn't call it by mistake!") unless @_;
  confess("Value is undefined") unless defined $value;
  if ($value eq 'yes') {
    return 1;
  } elsif ($value eq 'no') {
    return 0;
  }
  die("expected yes or no\n");
}

#
# Check value is a hexadecimal number with exactly 6 digits.
# Convert into a decimal number before returning.
#
sub _value_colour {
  my ($self, $value) = @_;
  confess("No parameters provided! Check that you used a function reference and didn't call it by mistake!") unless @_;
  confess("Value is undefined") unless defined $value;
  return hex($value) if ($value =~ m/^[0-9A-Fa-f]{6}$/);
  die("expected hexadecimal colour value\n");
}

sub _struc_command_line {
  my %options = @_;
  my $tag = (defined $options{ELE} ? $options{ELE} : 'command_line');
  my $repeat = (defined $options{RPT} ? $options{RPT} : $RPT_ONE);
  my @args; 
  my $start_cmd = sub {
    my ($self) = @_;
    @args = ();
  };
  my $end_arg = sub {
    my ($self, $arg) = @_;
    push(@args, $arg);
  };
  my $end_cmd = sub {
    my ($self) = @_;
    return @args;
  };
  return (
    {
      ELE => $tag,
      RPT => $repeat,
      HANDLE => $end_cmd,
      WATCH => $start_cmd, 
      TO => [
        {ELE => 'arg', RPT => $RPT_ALO, VAL => $any_multiline_re, END => $end_arg}
      ]
    }
  );
}

sub _struc_alphabet {
  my %options = @_;
  my $tag = (defined $options{ELE} ? $options{ELE} : 'alphabet');
  my $repeat = (defined $options{RPT} ? $options{RPT} : $RPT_ONE);
  my $store_alph_as = (defined $options{ALPH_AS} ? $options{ALPH_AS} : '');
  my $store_atrs_as = (defined $options{ATRS_AS} ? $options{ATRS_AS} : '');
  my $store_id2sym_as = (defined $options{ID2SYM_AS} ? $options{ID2SYM_AS} : '');
  my $store_sym2id_as = (defined $options{SYM2ID_AS} ? $options{SYM2ID_AS} : '');
  my $alph;
  my %sym_to_id;
  my %id_to_sym;
  # handle start tag of alphabet
  my $fn_alphabet_start = sub {
    my ($self, $attrs, $name, $like) = @_;
    $alph = new Alphabet();
    %sym_to_id = ();
    %id_to_sym = ();
    $alph->parse_header($name, $like);
  };
  # handle alphabet letter
  my $fn_alphabet_letter = sub {
    my ($self, $attrs, $id, $symbol, $aliases, $complement, $equals, $name, $colour) = @_;
    $alph->parse_symbol($symbol, $name, $colour, $complement, $equals, $aliases);
    $id_to_sym{$id} = $symbol;
    $sym_to_id{$symbol} = $id;
  };
  # handle end tag of alphabet
  my $fn_alphabet_end = sub {
    my ($self) = @_;
    $alph->parse_done();
    if ($store_alph_as) {
      $self->{$store_alph_as} = $alph;
    }
    if ($store_atrs_as) {
      my @atrs = ();
      for (my $i = 0; $i < $alph->size_core(); $i++) {
        my $symbol = $alph->char($i);
        my $id = $sym_to_id{$symbol};
        push(@atrs, {ATR => $id, VAL => $num_re});
      }
      $self->{$store_atrs_as} = \@atrs;
    }
    if ($store_id2sym_as) {
      $self->{$store_id2sym_as} = \%id_to_sym;
    }
    if ($store_sym2id_as) {
      $self->{$store_sym2id_as} = \%sym_to_id;
    }
    # return the alphabet
    return ($alph);
  };
  # return the structure to parse an alphabet
  return (
    {
      ELE => $tag,
      RPT => $repeat,
      HANDLE => $fn_alphabet_end,
      ATRS => [
        {ATR => 'name', VAL => $text_re},
        {ATR => 'like', OPT => qr/^(rna|dna|protein)$/}
      ],
      WATCH => $fn_alphabet_start,
      TO => [
        {
          ELE => 'letter',
          RPT => $RPT_ALO,
          ATRS => [
            {ATR => 'id', VAL => $attribute_re},
            {ATR => 'symbol', VAL => $sym_re},
            {ATR => 'aliases', OPT => $syms_re},
            {ATR => 'complement', OPT => $sym_re},
            {ATR => 'equals', OPT => $syms_re},
            {ATR => 'name', OPT => $text_re},
            {ATR => 'colour', OPT => \&_value_colour}
          ],
          WATCH => $fn_alphabet_letter
        }
      ]
    }
  );
}
