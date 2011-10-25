Breakpointer is a novel tool for locating breakpoints of structural variants (SV) from the alignment of single end reads (SE) produced by next generation sequencing (NGS). It adopts a heuristic method in searching for local mapping signature as a consequence of structural variants. You can redistribute it and/or modify it under the terms of the GNU General Public License.

Learn More
---

update soon
    

Installation
---

1. Make sure you have installed Bamtools (https://github.com/pezmaster31/bamtools).
2. Change the directories in Makefile to wherebam/lib/ and wherebam/include/ .
3. make

There are also built-in Linux X86_64 binaries of breakpointer and breakmis

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
Ihnestra<DF>e 63-73, D-14195 Berlin, Germany

Email: ruping@molgen.mpg.de
Project Website: https://github.com/ruping/Breakpointer
