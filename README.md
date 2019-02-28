# SLR
SLR is a scaffolding tool based on long reads and contig classification.
=========
License
=========

Copyright (C) 2014 Junwei Luo(luojunwei@hpu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Junwei Luo(luojunwei@hpu.edu.cn),
College of Computer Science and Technology,
Henan Polytechnic University,
Jiaozuo,
454000,
China


FLRO
=================

1) Introduction
```
    FLRO is a tool which aims to detect overlap regions among long reads..
    The input long read data of FLRO is the long reads (fasta format).
```
2) Before installing and running
```
    Please build and install DSK which can count kmer frequency of long reads
```
3) Installing.
```
    FLRO should run on Linux operating sysetm with gcc. We test FLRO using gcc4.6.3 on Ubuntu.
    Create a main directory (eg:FLRO). Copy all source code to this directory.
	cd FLRO
	make all
```
4) Running.
```
    Step 1: Use DSK to crete the kmer frequency file.
    Step 2: FLRO <long-read-file> <kmer-frequency-file> <kmer-length> <step> <low-kmer-frequency> <high-kmer-frequency> <output-file-name> <thread-number>
    Example:
	FLRO sra.fasta kmer 13 1 2 100 result.fa 10
```
5) Output.
```
    The output file "output-file-name" is the overlap result.
```
