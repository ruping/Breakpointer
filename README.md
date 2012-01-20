Breakpointer is a fast tool for locating sequence breakpoints from the alignment of single end reads (SE) produced by next generation sequencing (NGS). It adopts a heuristic method in searching for local mapping signatures created by insertion/deletions (indels) or more complex structural variants(SVs). With current NGS single-end sequencing data, the output regions by Breakpoint mainly contain the approximate breakpoints of indels and a limited number of large SVs. Notably, Breakpointer can uncover breakpoints of insertions which are longer than the read length. Breakpointer also can find breakpoints of many variants located in repetitive regions. The regions can be used not only as a extra support for SV predictions by other tools (such as by split-read method), but also can serve as a database for searching variants which might be missed by other tools. Breakpointer is a command line tool that runs under linux system. You can redistribute it and/or modify it under the terms of the GNU General Public License.


Learn More
---

Breakpointer takes advanage of two local mapping features of single-end reads as a consequence of indel/SVs: 1) non-uniform read distribution (depth skewness) and 2) misalignments at the boundaries of indel/SVs. We summarize these features as "breakpoint signature". Breakpointer proceeds in three stages in capturing this signature. It is implemented in C++ and perl. Input is the file or files containing alignments of single-end reads against a reference genome (in .BAM format). Output is the predicted regions containing potential breakpoints of SVs (in .GFF format). To be able to read in .BAM files, Breakpointer requires bamtools API, which users should install beforehand.
    

Installation
---

First, download breakpointer, unzip.

	unzip ruping-Breakpointer-XXX.zip
	cd ruping-Breakpointer-XXX

Second, make sure you have installed Bamtools (https://github.com/pezmaster31/bamtools). Then write down the bamtools_directory where ./lib/ and ./include/ sub-directories are located. Currently bamtools 2.0 has been tested.

Third, run make in following way to install

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
	--winsize   <int>   the window size, default is 10 for < 50bp reads, 20 for longer/variable length reads.
	--readlen   <int>   the length of the read (default: using variable read length).
	--mapping   <string>  the mapping file in BAM format. It could be an individual BAM file or a file listing the filenames of multiple BAM files (line seperate).All the BAM files must be sorted SAMELY according to chromosomes and coordinates. They should contain header tag "@HD   VN:1.0  SO:coordinate".
	--outdir   <string>  the output directory (default: current directory)
	--unique   <0/1>   0: take all the alignments (default), 1: take only unique alinged reads. If your BAM files only contain uniquely mapped reads or only a few non-unique reads, we recommand to leave it as default (0). If the BAM files contain many multi-location alignments, it is better to set it to 1. However, since different mappers generate different tags for uniqueness, if 1 is set, user shoule provide unique tag info (see tag/val_uniq).
	--tag_uniq  <string>  the tag in the BAM file denotating whether a read is uniquely mapped (default "XT" is taken as output from BWA).
	--val_uniq  <int>   the value for the above tag of uniquely mapped reads (default "85" is taken as from the output from BWA).
	--mistag   <string>  the bam tag for mismatch string, usually it is MD, but user can define it by using this option.
	--qualclip        whether to do the quality clipping for mismatch screening, default no.
	--noexecute        Running pipeline without executing the program, for testing purpose only.
	--runlevel  <int>   The stages of runlevel, 3 in total, either set with individual level "1" or multi levels like "1-3" (default). runlevel 1: scan the read alignment, searching for depth skewed regions; runlevel 2: mismatch screeing for each depth skewed region; runlevel 3: validate each candidate region by looking for support from unmappable reads.
	--unmap   <string>   File containing unmapped reads, either one file or a file listing the names of multiple files. must be fasta/fastq format.
	--help             print this help message.



Contact
---
Sun Ruping

Dept. Vingron (Computational Molecular Biology)
Max Planck Institute for Molecular Genetics. Ihnestr. 63-73, D-14195 Berlin, Germany

Email: ruping@molgen.mpg.de

Project Website: https://github.com/ruping/Breakpointer
