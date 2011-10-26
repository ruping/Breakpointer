Breakpointer is a fast tool for locating breakpoints of structural variants (SV) from the alignment of single end reads (SE) produced by next generation sequencing (NGS). It adopts a heuristic method in searching for local mapping signature created by structural variants. Breakpointer is a command line tool that runs under linux system. You can redistribute it and/or modify it under the terms of the GNU General Public License.

Learn More
---

update soon
    

Installation
---

First, download breakpointer, unzip.

	unzip ruping-Breakpointer-XXX.zip
	cd ruping-Breakpointer-XXX

Second, make sure you have installed Bamtools (https://github.com/pezmaster31/bamtools).
Third, write down the bamtools_directory where ./lib/ and ./include/ sub-directories are located.
Last, make in following way

	make BAMTOOLS_ROOT=/bamtools_directory/

5. You will see a directory called "breakpointer", within which you will find the pipeline script.

(There are also built-in Linux X86_64 binaries inside ./prebuilt/)

Usage
---

Run perl script breakpointer_run.pl:

	perl breakpointer_run.pl [options]
    

Or you could run step by step:

	breakpointer [options]
	breakmis [options]


Important
---

to be added



Contact
---
Sun Ruping
Dept. Vingron (Computational Molecular Biology)
Max Planck Institute for Molecular Genetics
Ihnestr. 63-73, D-14195 Berlin, Germany

Email: ruping@molgen.mpg.de
Project Website: https://github.com/ruping/Breakpointer
