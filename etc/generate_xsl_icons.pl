#!/usr/bin/perl
use strict;
use warnings;

use Fcntl qw(O_RDONLY);
use File::Basename;
use MIME::Base64;

print "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
print "<xsl:stylesheet version=\"1.0\" xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\">\n";

# Read in list of png files, output XLS variables
foreach my $image_path (@ARGV) {
  my $name = fileparse($image_path);
  my $mime;
  if ($name =~ m/\.png$/) {
    $mime = "image/png";
  } elsif ($name =~ m/\.gif$/) {
    $mime = "image/gif";
  } elsif ($name =~ m/\.jpg$/) {
    $mime = "image/jpeg";
  } else {
    die("unknown image type");
  }
  my $fh;
  sysopen($fh, $image_path, O_RDONLY);
  my $raw = do{ local $/ = undef; <$fh>; };
  close($fh);
  print "  <xsl:variable name=\"", $name, "\" select=\"'data:", $mime,
      ";base64,", encode_base64($raw, ""), "'\"/>\n";
}

print "</xsl:stylesheet>";
