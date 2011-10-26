Breakpointer is a fast tool for locating breakpoints of structural variants (SV) from the alignment of single end reads (SE) produced by next generation sequencing (NGS). It adopts a heuristic method in searching for local mapping signature created by structural variants. Breakpointer is a command line tool that runs under linux system. You can redistribute it and/or modify it under the terms of the GNU General Public License.

Learn More
---

update soon
    

Installation
---

First, download breakpointer, unzip.

	unzip ruping-Breakpointer-XXX.zip
	cd ruping-Breakpointer-XXX

Second, make sure you have installed Bamtools (https://github.com/pezmaster31/bamtools). Then write down the bamtools_directory where ./lib/ and ./include/ sub-directories are located.

Third, run make in following way

	make BAMTOOLS_ROOT=/bamtools_directory/

You will see a directory called "breakpointer", within which you will find the pipeline script and binaries.


Usage
---

Run perl script breakpointer_run.pl:

	perl breakpointer_run.pl [options]
    

Or you could run step by step:

	breakpointer [options]
	breakmis [options]


Options
---
	--winsize <int> the window size, default is 10 for < 50bp reads, 20 for longer reads.
	--readlen <int> the length of the read (now only support fixed length)
	--mapping <string> the mapping file in BAM format. It could be an individual BAM file or a file listing the filenames of multiple BAM files (line seperate).All the BAM files must be sorted SAMELY according to chromosomes and coordinates. They should contain header tag "@HD   VN:1.0  SO:coordinate".
	--outdir <string> the output directory (default: current directory)
	--unique <0/1> 0: take all the alignments (default), 1: take only unique alinged reads. If your BAM files only contain uniquely mapped reads or only a few non-unique reads, we recommand to leave it as default (0). If the BAM files contain many multi-location alignments, it is better to set it to 1. However, since different mappers generate different tags for uniqueness, if 1 is set, user shoule provide unique tag info (see tag/val_uniq).
	--tag_uniq <string> the tag in the BAM file denotating whether a read is uniquely mapped (default "XT" is taken as output from BWA).
	--val_uniq <int> the value for the above tag of uniquely mapped reads (default "85" is taken as from the output from BWA).
	--noexecute Running pipeline without executing the program, for testing purpose only.
	--runlevel <int> The stages of runlevel, 3 in total, either set with individual level "1" or multi levels like "1-3" (default). runlevel 1: scan the read alignment, searching for depth skewed regions; runlevel 2: mismatch screeing for each depth skewed region; runlevel 3: validate each candidate region by looking for support from unmappable reads.
	--unmap <string> File containing unmapped reads, either one file or a file listing the names of multiple files. must be fasta/fastq format.
	--help print this help message.



Contact
---
Sun Ruping
Dept. Vingron (Computational Molecular Biology)
Max Planck Institute for Molecular Genetics
Ihnestr. 63-73, D-14195 Berlin, Germany

Email: ruping@molgen.mpg.de
Project Website: https://github.com/ruping/Breakpointer
