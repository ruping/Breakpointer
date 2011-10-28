package CIGAR;

use strict;
use Carp;

sub cigarlen {

  my ($class, $cigar) = @_;
  my $ciglen = 0;

  while ( $cigar =~ /(\d+)[MDNS]/g ) {
    $ciglen += $1;
  }

  return $ciglen;
}

1;
