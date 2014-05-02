#!/usr/bin/perl

#******************************************************************************#
#  breakpointer_run.pl @ Breakpointer (version 0.1)
#  The pipeline script of Breakpointer.
#
#  (c) 2011 - Sun Ruping
#  Dept. Vingron (Computational Mol. Bio.)
#  Max-Planck-Institute for Molecular Genetics
#  Ihnestr. 73, D-14195, Berlin, Germany
#  EMAIL: ruping@molgen.mpg.de
#
#  Breakpointer is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License
#******************************************************************************#

use strict;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use Fof;
use Chrstart;

my $noexecute = 0;
my $runlevels = 0;
my $winsize;
my $mapfile = "";
my $readlen = 0;
my $out_dir = "./";
my $unique = 0;
my $qual_clip = 0;
my $mistag = "MD";
my $unmapped = "";
my $basename = "";
my $help;
my $BP = "$RealBin/";

#if no argument is given, print help
if (@ARGV == 0) {
  helpm();
} else {
  printf STDERR "\n# $0 %s\n",join(" ",@ARGV);
}

GetOptions (
            "runlevel=s"   => \$runlevels,
            "noexecute"    => \$noexecute,
            "windowsize=i" => \$winsize,
            "mapping=s"    => \$mapfile,
            "readlen=i"    => \$readlen,
            "outdir=s"     => \$out_dir,
            "unique=i"     => \$unique,
            "mistag=s"     => \$mistag,
            "qualclip"     => \$qual_clip,
            "unmap=s"      => \$unmapped,
            "basename=s"   => \$basename,
            "help|h"       => \$help,
	   );


#if help, print help
helpm() if ($help);

my %runlevel;
if ($runlevels != 0) {
  foreach my $r (split /,/,$runlevels) {
    my $from=1;
    my $to=20;
    if ($r=~/^(\d+)/) {
      $from=$1;
    }
    if ($r=~/\-(\d+)$/) {
      $to=$1;
    } elsif ($r!~/\-/) {
      $to=$from;
    }
    for (my $i=$from;$i<=$to;$i++) {
      $runlevel{$i}=1;
    }
  }
}
else{
  for (1..3){
   $runlevel{$_} = 1;
  }
}


###
###runlevel0.5: preparation
###

printtime();
printf STDERR "Checking options:\n";

if ($runlevels == 0){
  print STDERR "no runlevel is given, will take default (1-3).\n";
}
else{
  printf STDERR "runlevel ".join(",", sort keys %runlevel)." will be run.\n";
}

if (!$winsize){
  printf STDERR "warning: no window size is given, will take default.\n";
}
else{
  printf STDERR "window size is $winsize.\n";
}

if ($mapfile eq ""){
  printf STDERR "warning: no mapping file is given, stop running.\n";
  helpm();
}
else {
  $basename = basename($mapfile);
  unless ( -r $mapfile ) {
    printf STDERR "warning: $mapfile is not readable or not present.\n";
    exit(0);
  }
  my @mapping_files = Fof->fileorfilename($mapfile);
  printf STDERR "mapping files are:\n".join("\n",@mapping_files)."\n";
  foreach my $mapping ( @mapping_files ) {
    unless ( -r $mapping) {
      printf STDERR "warning: $mapping is not readable or not existing, please check the file mode.\n";
      exit(0);
    }
  }
}

if (-w $out_dir) {
  printf STDERR "output directory is $out_dir.\n";
}
else {
  printf STDERR "$out_dir is not writable\n";
  exit(0);
}
printf STDERR "output basename is $basename.\n";

if ($unique == 0){
  printf STDERR "take all mapped reads.\n";
}
else{
  printf STDERR "take only uniquely mapped reads.\n";
}

printf STDERR "unique tag in the mapping file is $tag_uniq.\n";
printf STDERR "value for the unique tag is $val_uniq.\n";

if ($readlen == 0){
  printf STDERR "warning: no readlen is given, will take default 36bp.\n";
}
else{
  printf STDERR "read length is $readlen.\n";
}

if ($unmapped eq ""){
  printf STDERR "warning: no unmapped reads file is given, runlevel3 will be skipped\n";
  delete $runlevel{3};
}
else {
  unless ( -r $unmapped ) {
    printf STDERR "warning: $unmapped is not readable or not present.\n";
    exit(0);
  }
  my @umr_files = Fof->fileorfilename($unmapped);
  printf STDERR "unmapped reads files is:\n".join("\n",@umr_files)."\n";
  foreach my $umr (@umr_files) {
    unless (-r $umr and $umr ne "") {
       printf STDERR "warning: $umr is not readable, please check the file mode.\n";
       exit(0);
    }
  }
}

#print Dumper(\%runlevel);

###
###runlevel1:  scan the read alignment, searching for depth skewed regions
###
$runlevels=1;
if (exists $runlevel{$runlevels}) {
  printtime();
  printf STDERR "RUNLEVEL 1: scan the read alignment, searching for depth skewed regions\n";

  my $endskew = $basename."\.endskew";
  my $op_mapfile = "--mapping $mapfile";
  my $op_winsize = "";
  if ($winsize){$op_winsize = "--windowsize $winsize";}
  my $op_unique = "";
  if ($unique != 0) {$op_unique = "--unique $unique";}
  my $op_readlen = "";
  if ($readlen != 0) {
    $op_readlen = "--readlen $readlen";
  }

  my $cmd = "$BP/breakpointer $op_mapfile $op_winsize $op_readlen $op_unique >$out_dir/$endskew";
  if (-e "$out_dir/$endskew") {
    printf STDERR "$out_dir/$endskew exists, skip running RUNLEVEL 1\n";
  } else {
    RunCommand($cmd,$noexecute);
  }
  printf "RUNLEVEL 1 done\n";
}



###
###runlevel2: mismatch screening for each depth skewed region
###
$runlevels=2;
if (exists $runlevel{$runlevels}) {
  printtime();
  printf STDERR "RUNLEVEL 2: mismatch screening for each depth skewed region\n";

  my $region_f = "$out_dir/$basename.endskew";
  unless (-r $region_f){
    printf STDERR "warning: $region_f is not readable!\n";
  }

  my $endskewmis = $basename."\.endskew\.mis\.gff";
  my $op_regionf = "--region $region_f";
  my $op_mapping = "--mapping $mapfile";
  my $op_readlen = "";
  if ($readlen != 0) {
    $op_readlen = "--readlen $readlen";
  }
  my $op_unique = "";
  if ($unique != 0) {
    $op_unique = "--unique $unique";
  }
  my $op_qualclip = "--qualclip no";
  if ($qual_clip){
    $op_qualclip = "--qualclip phred33";
  }
  my $op_mistag = "";
  if ($mistag ne "MD"){
    $op_mistag = "--mistag $mistag";
  }

  my $cmd = "$BP/breakmis $op_regionf $op_mapping $op_readlen $op_unique $op_qualclip $op_mistag >$out_dir/$endskewmis";
  if (-e "$out_dir/$endskewmis") {
    printf STDERR "$out_dir/$endskewmis exists, skip running RUNLEVEL 2\n";
  } else {
    RunCommand($cmd,$noexecute);
  }
  printf STDERR "RUNLEVEL 2 done\n";

}

###
###runlevel3: validation using unmapped reads
###

$runlevels=3;
if ( exists $runlevel{$runlevels} ) {
  printtime();
  printf STDERR "RUNLEVEL 3: validation using unmapped reads\n";

  my $vali = $basename."\.endskew\.mis\.vali\.gff";
  my $op_umr = "--umr $unmapped";
  my $op_readlen = "--readlen 36";
  if ($readlen != 0) {
    $op_readlen = "--readlen $readlen";
  }
  my $ermis = $basename."\.endskew\.mis\.gff";
  unless (-e "$out_dir/$ermis") {
    printf "warning: $out_dir/$ermis is not readable.\n";
    exit(0);
  }
  my $op_ermis = "--ermis $out_dir/$ermis";

  my $cmd = "perl $BP/breakvali.pl $op_umr $op_readlen $op_ermis --verbose >$out_dir/$vali";
  if (-e "$out_dir/$vali") {
    printf STDERR "$out_dir/$vali exists, skip running RUNLEVEL 3\n";
  } else {
    RunCommand($cmd,$noexecute);
  }
  printf STDERR "RUNLEVEL 3 done\n";
}

printtime();
printf STDERR "DONE.\n";


###
###subroutine region----------------------------
###

sub RunCommand {
  my ($command,$noexecute) = @_ ;
  print STDERR "$command\n";
  unless ($noexecute) {
    system($command);
  }
}


sub printtime {
  my @time = localtime(time);
  printf STDERR "\n[".($time[5]+1900)."\/".($time[4]+1)."\/".$time[3]." ".$time[2].":".$time[1].":".$time[0]."]\t";
}

sub helpm {
  print "\nBreakpointer v0.1 usage: $0 [options]\n\nOptions:\n";
  print "\t--windowsize\t<int>\t\tthe window size, default is 10 for < 50bp reads, 20 for longer reads/variable length reads.\n";
  print "\t--readlen\t<int>\t\tthe length of the read (default: using variable read length)\n";
  print "\t--mapping\t<string>\tthe mapping file in BAM format. It could be an individual BAM file or a file listing the filenames of multiple BAM files (line seperated).\n\t\t\t\t\tAll the BAM files must be sorted SAMELY according to chromosomes and coordinates. They should contain header tag \"\@HD\tVN:1.0\tSO:coordinate\".\n";
  print "\t--outdir\t<string>\tthe output directory (default: current directory)\n";
  print "\t--unique\t<0/1>\t\t0: take all the alignments (default), 1: take only unique alinged reads.\n\t\t\t\t\tIf your BAM files only contain uniquely mapped reads or only a few non-unique reads, we recommand to leave it as default (0).\n\t\t\t\t\tIf the BAM files contain many multi-location alignments, it is better to set it to 1.\n\t\t\t\t\tHowever, since different mappers generate different tags for uniqueness, if 1 is set, user shoule provide unique tag info (see tag/val_uniq).\n";
  print "\t--mistag\t<string>\tthe bam tag for mismatch string, usually it is MD, but user can define it by using this option.\n";
  print "\t--qualclip\t\t\twhether to do the quality clipping for mismatch screening, default no.\n";
  print "\t--noexecute\t\t\tRunning pipeline without executing the program, for testing purpose only.\n";
  print "\t--runlevel\t<int>\t\tThe stages of runlevel, 3 in total, either set with individual level \"1\" or multi levels like \"1-3\" (default) \n";
  print "\t\t\t\t\trunlevel 1: scan the read alignment, searching for depth skewed regions;\n";
  print "\t\t\t\t\trunlevel 2: mismatch screeing for each depth skewed region;\n";
  print "\t\t\t\t\trunlevel 3: validate each candidate region by looking for support from unmappable reads.\n";
  print "\t--unmap\t\t<string>\tFile containing unmapped reads, either one file or a file listing the names of multiple files. must be fasta/fastq format.\n";
  print "\t--basename\t<string>\tthe basename of the output files (default: take the basename of the mapping files)\n";
  print "\t--help\t\t\t\tprint this help message.\n\n\n";
  exit 0;
}


exit 22;
