# FLRO
FLRO can detect overlap regions among long reads.
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
    Note:
    	Each line in the kmer-frequency-file should be "kmer kmer-frequency". 
	For example:
	AGTCCAGGCCGGG 3
	GAAATCCAGCCGC 6
	AACCGGCGAATCG 3
	TATTTTAACATTC 2
	TATGGCCGATGAA 4
	AAAGCCGAAGCCT 3
	CATCTTCACATCA 2
	ATAAGTGATAGCT 4
	TCGGCCATATTAC 4
	ATTATTGCAATAC 6
     An example of command line is shown below.
	FLRO sra.fasta kmer 13 1 2 100 result.fa 10
```
5) Output.
```
    The output file "output-file-name" is the overlap result.
     The first column is the first read number.
    The second column is the second read number.
    The third column is aligning orientation. 0 represents forward aligning. 1 represents reverse aligning.
    The fourth column is starting position in the first read.
    The fifth column is ending position in the first read.
    The sixth column is the length of the first read.
    The seventh column is starting position in the second read.
    The eighth column is ending position in the second read.
    The ninth column is the length of the second read.
    
    One example in the result file is shown below.
    36423,1,0,5326,9923,9923,1,4364,10479
    It means the region [5326,9923] in the first read is overlapped with the region [1,4364] in the second read.
```
