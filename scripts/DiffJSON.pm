package DiffJSON;

use warnings;
use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(diff_json);

use Fcntl qw(O_RDONLY);
use Data::Dumper;

use JsonRdr;
use JsonWr qw(convert_string);

my $TYPE_START_DATA = 1;
my $TYPE_END_DATA = 2;
my $TYPE_PROPERTY = 3;
my $TYPE_ATOM_NULL = 4;
my $TYPE_ATOM_BOOL = 5;
my $TYPE_ATOM_STRING = 6;
my $TYPE_ATOM_NUMBER = 7;
my $TYPE_START_OBJECT = 8;
my $TYPE_END_OBJECT = 9;
my $TYPE_START_ARRAY = 10;
my $TYPE_END_ARRAY = 11;

#
# Take 2 json files (json1, json2) and a patterns to ignore.
#
# If json1 or json2 is invalid json it should return a sentence describing
# the first error it finds including where in the file it is.
#
# If json1 is identical to json2 it should return an empty string.
#
# If json1 is different to json2 is should return a sentence describing the 
# first difference it finds including where in the file it is.
#
# Properties should be ordered identically and have identical values.
#
# Every part of the json file is given a name and users can pass regular
# expressions which are matched against these names to determine if an
# element and all its sub elements should be ignored.
#
sub diff_json {
  my ($jsonfilename1, $jsonfilename2, @patterns) = @_;
  my ($fh1, $fh2);

  # make regular expressions from the patterns to ignore
  my @ignore_regex = ();
  for my $pattern (@patterns) {
    push(@ignore_regex, qr/$pattern/);
  }

  sysopen($fh1, $jsonfilename1, O_RDONLY) or die("Failed to open $jsonfilename1");
  sysopen($fh2, $jsonfilename2, O_RDONLY) or die("Failed to open $jsonfilename2");

  my ($s1, $p1) = &make_parser(@ignore_regex);
  my ($s2, $p2) = &make_parser(@ignore_regex);

  my $message = '';
  while (!$message && &has_more($s1) && &has_more($s2)) {
    # get some more events
    $message = &parse_more($s1, $p1, $fh1) unless &has_event($s1);
    last if $message;
    $message = &parse_more($s2, $p2, $fh2) unless &has_event($s2);
    last if $message;
    # compare all the events we have
    my ($e1, $e2);
    while (&has_event($s1) && &has_event($s2)) {
      $e1 = &take_event($s1);
      $e2 = &take_event($s2);
      if (!&event_equal($e1, $e2)) {
        $message = "A difference was detected between the reference " . 
            "and tested json files.\n".
            "Got " . $e1->{name} . " but expected " . $e2->{name} . ".\n";
        last;
      }
    }
  }

  if (!$message && (&has_more($s1) || &has_more($s2))) {
    # I don't think this state can occur because the end tag
    # would not match
    $message = "Unmatched events\n";
  }

  close($fh1);
  close($fh2);

  return $message;
}

sub event_equal {
  my ($e1, $e2) = @_;

  if ($e1->{type} != $e2->{type}) {
    return 0;
  } else {
    my $ty = $e1->{type};
    if ($ty == $TYPE_START_DATA || $ty == $TYPE_PROPERTY || $ty == $TYPE_ATOM_STRING) {
      return ($e1->{data} eq $e2->{data});
    } elsif ($ty == $TYPE_ATOM_NUMBER) {
      return ($e1->{data} == $e2->{data});
    } elsif ($ty == $TYPE_ATOM_BOOL) {
      return ((!!$e1->{data}) == (!!$e2->{data}));
    }
    return 1;
  }
}

sub make_parser {
  my $state = {
    ignore => \@_,
    finished => 0,
    error => "",
    ignoring => 0,
    chain => [],
    events => []  # queue of {type => TYPE, name => FULL_NAME, data => TEXT_OR_UNDEF}
  };
  my $parser = new JsonRdr($state, 
    start_data => \&_start_data,
    end_data => \&_end_data,
    start_property => \&_start_property,
    end_property => \&_end_property,
    atom => \&_atom,
    start_object => \&_start_object,
    end_object => \&_end_object,
    start_array => \&_start_array,
    end_array => \&_end_array
  );

  return ($state, $parser);
}

sub parse_more {
  my ($state, $parser, $fh) = @_;
  return if ($state->{finished});
  eval {
    my $buffer;
    # read until the state contains some data to compare
    while (scalar(@{$state->{events}}) == 0 && !$state->{error}) {
      my $count = read($fh, $buffer, 100);
      die("Error reading file") unless defined($count);
      if ($count == 0) {
        $state->{finished} = 1;
        $parser->parse_done();
        last;
      }
      $parser->parse_more($buffer);
    }
  };
  if ($@) {
    return $@;
  } else {
    return '';
  }
}

sub has_event {
  my ($state) = @_;
  return scalar(@{$state->{events}}) > 0;
}

sub take_event {
  my ($state) = @_;
  my $event = shift(@{$state->{events}});
  return $event;
}

sub has_more {
  my ($state) = @_;
  return !$state->{finished} || scalar(@{$state->{events}}) > 0
}

sub _chain_to_name {
  my ($state) = @_;

  my $name = '';
  foreach my $link (@{$state->{chain}}) {
    my $type = $link->{type};
    if ($type == $TYPE_START_DATA) {
      $name .= $link->{data};
    } elsif ($type == $TYPE_START_OBJECT) {
      $name .= ':';
    } elsif ($type == $TYPE_PROPERTY) {
      $name .= $link->{data};
    } elsif ($type == $TYPE_START_ARRAY) {
      $name .= '@' . ($link->{index} > 0 ? $link->{index} : '');
    } else {
      $name .= '=';
      if ($type == $TYPE_ATOM_NULL) {
        $name .= 'null';
      } elsif ($type == $TYPE_ATOM_BOOL) {
        $name .= ($link->{data} ? 'true' : 'false');
      } elsif ($type == $TYPE_ATOM_NUMBER) {
        $name .= $link->{data};
      } else {
        $name .= &JsonWr::convert_string($link->{data});
      }
    }
  }
  return $name;
}

sub _should_ignore {
  my ($state) = @_;
  my $name = &_chain_to_name($state);
  for my $ignore_re (@{$state->{ignore}}) {
    return 1 if ($name =~ m/$ignore_re/);
  }
  return 0;
}

sub _queue_event {
  my ($state, $type, $data) = @_;
  my $name = &_chain_to_name($state);
  push(@{$state->{events}}, {type => $type, data => $data, name => $name});
}

sub _start_event {
  my ($state, $type, $data) = @_;
  unless ($state->{ignoring}) {
    if (&_should_ignore($state)) {
      $state->{ignoring} = 1;
    } else {
      &_queue_event($state, $type, $data);
    }
  } else {
    $state->{ignoring} += 1;
  }
}

sub _end_event {
  my ($state, $type, $data) = @_;
  unless ($state->{ignoring}) {
    &_queue_event($state, $type, $data);
  } else {
    $state->{ignoring} -= 1;
  }
}

sub _start_data {
  my ($state, $var_name) = @_;
  die("var_name is not defined") unless defined $var_name;
  push(@{$state->{chain}}, {type => $TYPE_START_DATA, data => $var_name});
  &_start_event($state, $TYPE_START_DATA, $var_name);
}

sub _end_data {
  my ($state) = @_;
  pop(@{$state->{chain}});
  &_end_event($state, $TYPE_END_DATA);
}

sub _start_property {
  my ($state, $name) = @_;
  die("property name is not defined") unless defined $name;
  push(@{$state->{chain}}, {type => $TYPE_PROPERTY, data => $name});
  &_start_event($state, $TYPE_PROPERTY, $name);
}

sub _end_property {
  my ($state) = @_;
  pop(@{$state->{chain}});
  $state->{ignoring} -= 1 if ($state->{ignoring});
}

sub _atom {
  my ($state, $value, $is_null, $is_bool, $is_string, $is_num, $is_int, $is_sci) = @_;
  die("value is not defined") unless ($is_null || defined $value);
  my $type;
  if ($is_null) {
    $type = $TYPE_ATOM_NULL;
  } elsif ($is_bool) {
    $type = $TYPE_ATOM_BOOL;
  } elsif ($is_string) {
    $type = $TYPE_ATOM_STRING;
  } elsif ($is_num) {
    $type = $TYPE_ATOM_NUMBER;
  } else {
    die("Value has indeterminate type");
  }
  my $last = $state->{chain}->[-1];
  if ($last->{type} == $TYPE_START_ARRAY) {
    $last->{index} += 1;
  }

  push(@{$state->{chain}}, {type => $type, data => $value});
  &_queue_event($state, $type, $value) unless ($state->{ignoring});
  pop(@{$state->{chain}});
}

sub _start_object {
  my ($state) = @_;
  my $last = $state->{chain}->[-1];
  if ($last->{type} == $TYPE_START_ARRAY) {
    $last->{index} += 1;
  }
  push(@{$state->{chain}}, {type => $TYPE_START_OBJECT});
  &_start_event($state, $TYPE_START_OBJECT);
}

sub _end_object {
  my ($state) = @_;
  pop(@{$state->{chain}});
  &_end_event($state, $TYPE_END_OBJECT);
}

sub _start_array {
  my ($state) = @_;
  my $last = $state->{chain}->[-1];
  if ($last->{type} == $TYPE_START_ARRAY) {
    $last->{index} += 1;
  }
  push(@{$state->{chain}}, {type => $TYPE_START_ARRAY, index => 0});
  &_start_event($state, $TYPE_START_ARRAY);
}

sub _end_array {
  my ($state) = @_;
  pop(@{$state->{chain}});
  &_end_event($state, $TYPE_END_ARRAY);
}
