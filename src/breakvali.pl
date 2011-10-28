use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use Fof;

my $umr_file;
my $readlen;
my $ermis_file;
my $reads_file;
my $verbose;

GetOptions (
             "umr=s"      => \$umr_file,
             "ermis=s"    => \$ermis_file,
             "readlen=i"  => \$readlen,
             "reads=s"    => \$reads_file,
             "verbose|v"  => \$verbose,
             "help|h"     => sub{
                              print "\nusage: $0 [options]\n\n";
                              print "\t--help\t\tprint this help message\n";
                              print "\t--umr\t\tthe unmapped reads file or a file listing the names of the umr-files (fasta|fastq)\n";
                              print "\t--readlen\tthe length of the read \n";
                              print "\t--reads\t\tthe reads file or file of filenames (just for gff)\n";
                              print "\t--ermis\t\tthe mismatch file generated by breakpointer-breakmis \n\n\n";
                              exit 0;
                             },
           );


my @umr_files = Fof->fileorfilename($umr_file);
print STDERR "Input unmapped reads file is:\n".(join ("\n", @umr_files))."\n" if $verbose;

# reads files to read
my %reads;

if ($reads_file) {

  my @reads_files = Fof->fileorfilename($reads_file);
  print STDERR "Input reads file is:\n".(join ("\n", @reads_files))."\n" if $verbose;

  foreach my $reads (@reads_files) {
    my $type;
    my $rname;
    my $line;
    open READS, "$reads";
    while ( <READS> ) {
      chomp;
      $type = 'q' if (/^@/ && !defined $type);
      $type = 'a' if (/^>/ && !defined $type);
      if ($type eq 'q') {
        if ($_ =~ /^@(.+)$/) {
          $rname = $1;
          $_ = <READS>;
          chomp;
          $reads{$rname} = $_;
          $_ = <READS>;
          $_ = <READS>;
        }
      }
      if ($type eq 'a') {
        if ($_ =~ /^>(.+?)\|/) {
          $rname = $1;
          $_ = <READS>;
          chomp;
          $reads{$rname} = $_;
        }
      }
    }
    close READS;
    print STDERR "$reads loaded\n" if $verbose;
  }
}


# store Ends-Skewed data
# build hash for check
my %hash;     # for exists check
my %info;     # remember the id

my %ER;
open ER, "$ermis_file";
while ( <ER> ) {
  next if /^#/;
  chomp;
  my @cols = split /\t/;
  $cols[8] =~ /^ID=(.+?);.+BinomialScore=(.+?);.+seedseq=(.+)$/;
  my $ID = $1;
  my $S_B = $2;
  my $seedseq = $3;
  my $S_M = $cols[5];
  next if $seedseq eq 'RME';

  if ($reads_file) { #gff, need reads
    $seedseq =~ /^(.+?)\|.+\[(.+)\]/;
    my $rname  = $1;
    my $cutype = $2;
    if ($cutype eq '+p'){
      $seedseq = substr($reads{$rname}, 0, 25);
    }
    elsif ($cutype eq '+s'){
      $seedseq = substr($reads{$rname}, $readlen-25, 25);
    }
    elsif ($cutype eq '-p'){
      my $readseq = reverse($reads{$rname});
      $readseq =~ tr/ACGT/TGCA/;
      $seedseq = substr($readseq, 0, 25);
    }
    else {
      my $readseq = reverse($reads{$rname});
      $readseq =~ tr/ACGT/TGCA/;
      $seedseq = substr($readseq, $readlen-25, 25);
    }
    print STDERR "$rname\t$cutype\t$seedseq\n" if $verbose;
  }


  # strand + of seedseq
  $hash{$seedseq} = 0;
  push (@{$info{$seedseq}{'p'}}, $ID);

  # strand - of seedseq
  my $RC = reverse($seedseq);
  $RC =~ tr/ACGT/TGCA/;
  $hash{$RC} = 0;
  push (@{$info{$RC}{'n'}}, $ID);


  $ER{$ID}{'SB'}  = $S_B;
  $ER{$ID}{'SM'}  = $S_M;
  $ER{$ID}{'DD1'} = join("\t", $cols[0],$cols[1],$cols[2],$cols[3],$cols[4]);
  $ER{$ID}{'DD2'} = join("\t", $cols[6],$cols[7],$cols[8]);
}
close ER;
print STDERR "$ermis_file loaded\n" if $verbose;
my $numberreport = scalar(keys %hash);
print STDERR "$numberreport\n" if $verbose;


# check each unmappable read

my $count;
foreach my $umr (@umr_files) {
  my $umr_type;
  open UMR, "$umr";
  while ( <UMR> ) {

    chomp;
    $umr_type = 'q' if (/^@/ && !defined $umr_type);
    $umr_type = 'a' if (/^>/ && !defined $umr_type);

    if ($umr_type eq 'a') {
      if ( $_ =~ /^>/ ) {
        $_ = <UMR>;
        my $string = $_;
        for (my $i=0; $i<=$readlen-25; $i++) {
          $string = substr($_, $i, 25);
          if (exists $hash{$string}) {
            if ($info{$string}{$i} eq '') {
              $hash{$string}++;
              $info{$string}{$i} = 1;
            }
            last;
          }
        }
        $count++;
      }
    }
    if ($umr_type eq 'q') {
      if ( $_ =~ /^@/ ) {
        $_ = <UMR>;
        my $string = $_;
        for (my $i=0; $i<=$readlen-25; $i++) {
          $string = substr($_, $i, 25);
          if (exists $hash{$string}) {
            if ($info{$string}{$i} eq '') {
              $hash{$string}++;
              $info{$string}{$i} = 1;
            }
            last;
          }
        }
        $count++;
        $_ = <UMR>;
        $_ = <UMR>;
      }
    }
    print STDERR "$count unmapped reads processed\n" if ($verbose && $count%1000000 == 0);
  }
  close UMR;
  print STDERR "$umr finished\n";
}


# count the number of supporting reads
my %SU;
foreach my $key (keys %hash) {
  my $times = $hash{$key};
  foreach my $id (@{$info{$key}{'p'}}){
    $SU{$id}{'times'}  += $times;
  }
  foreach my $id (@{$info{$key}{'n'}}){
    $SU{$id}{'times'}  += $times;
  }
}

#clear mem
%hash = ();
%info = ();


# sorting

my $count = 1;
foreach my $ID ( sort {$ER{$a}{'SB'} <=> $ER{$b}{'SB'} } keys %ER) {
  next if $SU{$ID}{'times'} == 0;
  $ER{$ID}{'rank_SB'} = $count;
  $count++;
}

$count = 1;
foreach my $ID ( sort {$ER{$a}{'SM'} <=> $ER{$b}{'SM'} } keys %ER) {
  next if $SU{$ID}{'times'} == 0;
  $ER{$ID}{'rank_SM'} = $count;
  $count++;
}

#combining
foreach my $ID (keys %ER) {
  next if $SU{$ID}{'times'} == 0;
  my $SB = $ER{$ID}{'SB'};
  my $SM = $ER{$ID}{'SM'};
  my $SU = $SU{$ID}{'times'};
  my $rank_SB   = sprintf("%.3f", $ER{$ID}{'rank_SB'}*100/($count-1));
  my $rank_SM   = sprintf("%.3f", $ER{$ID}{'rank_SM'}*100/($count-1));
  my $RankScore = sprintf("%.3f", ($rank_SB+$rank_SM)/2);
  $ER{$ID}{'DD2'} .= ";MismatchScore=".$SM.";SU=".$SU.";rank_SB=".$rank_SB.";rank_SM=".$rank_SM;
  printf ("%s\n", join("\t", $ER{$ID}{'DD1'}, $RankScore, $ER{$ID}{'DD2'}));
}

exit 0;
