package Chrstart;

use strict;
use Carp;

sub getchrpos{

  my ($class, $samfiles) = @_;

  my $chr_old;

  my %chr_start;

  foreach my $samfile (@{$samfiles}) {

    open SAMFILE, "$samfile" or croak "The sam file read error!";

    while ( <SAMFILE> ) {

      next if /^@/;
      chomp;

      my @cols = split /\t/;
      my $chr  = $cols[2];

      $chr = 'chr'.$chr unless ($chr =~ /^chr/);
      $chr .= 'T' if ($chr eq 'chrM');

      if ($chr ne $chr_old) {
        $chr_start{$chr}{$samfile} = tell SAMFILE;
      }

      $chr_old = $chr;

    }

    close SAMFILE;
    print STDERR "$samfile chr_start loaded\n";

  }

  return %chr_start;
}

1;
