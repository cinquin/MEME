#!/usr/bin/perl
use strict;
use warnings;

use Carp;
use Fcntl qw(O_RDONLY);

sub _convert_string {
  my ($string) = @_;
  croak("string undefined") unless defined $string;
  my $storage = '"';
  for (my $i = 0; $i < length($string); $i++) {
    my $letter = substr($string, $i, 1);
    my $num = ord($letter);
    if ($letter eq '"') {
      $storage .= "\\\"";
    } elsif ($letter eq "\\") {
      $storage .= "\\\\";
    } elsif ($letter eq '/') {
      $storage .= "\\/";
    } elsif ($num >= 0 && $num <= 0x1F || $num >= 0x7F && $num <= 0x9F ||
      $num == 0x2028 || $num == 0x2029) {
      # control character, so encode
      # Note that techinically U+2028 and U+2029 are valid JSON but they are
      # not valid Javascript as they encode 'LINE SEPARATOR' (U+2028) and
      # 'PARAGRAPH SEPARATOR' (U+2029) which are treated by Javascript as
      # newline characters and are hence not valid in strings.
      if ($num == 8) { #backspace
        $storage .= "\\b";
      } elsif ($num == 9) { #tab
        $storage .= "\\t";
      } elsif ($num == 10) { # line feed (newline)
        $storage .= "\\n";
      } elsif ($num == 12) { # form feed
        $storage .= "\\f";
      } elsif ($num == 13) { # carriage return
        $storage .= "\\r";
      } else {
        $storage .= sprintf("\\u%.04x", $num);
      }
    } else {
      $storage .= $letter;
    }
  }
  $storage .= '"';
  return $storage;
}

# Read in EPS template
my ($fn_name, $template_path) = @ARGV;
die("Function name not defined") unless (defined($fn_name));
die("Function name does not match allowed pattern.") unless ($fn_name =~ m/^[a-z_][a-z0-9_]*$/);
die("Template file not specified") unless (defined($template_path) && -e $template_path);

my $fh;
sysopen($fh, $template_path, O_RDONLY);
print "function $fn_name(inputs) {\n";
print "  function _input(name) {\n";
print "    if (typeof inputs[name] === \"undefined\") {\n";
print "      throw new Error(\"Missing template variable: \" + name);\n";
print "    }\n";
print "    return inputs[name];\n";
print "  }\n";
my $line;
my $first = 1;
my $var_re = qr/\{\$([A-Za-z_][A-Z-a-z0-9_]*)\}/;
print "  return (\n";
while ($line = <$fh>) {
  print " +\n" unless $first;
  my $pos = 0;
  while ($line =~ m/$var_re/g) {
    print " + " if $pos > 0;
    if ($pos < $-[0]) {
      print _convert_string(substr($line, $pos, $-[0] - $pos));
      print " + ";
    }
    print "_input(\"" . $1 . "\")";
    $pos = $+[0];
  }
  print " + " if $pos > 0;
  print _convert_string(substr($line, $pos));
  $first = 0;
}
print "" if $first;
print "\n  );\n}";

