package Fof;

use strict;

#
# file or filenames?
#

sub fileorfilename {

  my ($class, $input) = @_;

  my @fof = ();

  my @files = map( <${_}>, $input );

  open IP, "$files[0]";

  my $nr=0;

  my $line;

  while ( ($line = <IP> ) && ( $nr < 20 || scalar(@fof) != 0)) {

      next if ($line =~ /^#/);   #skip comment

      if ($line =~ /^(\S+)/) {

         if (-r $1) {
            push @fof,$1;
         }
      }
      $nr++;
  }
  close IP;

  if (scalar(@fof) != 0) {
    @files = @fof;
  }

  return @files;

}

1;


