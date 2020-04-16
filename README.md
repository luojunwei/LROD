# LROD (Version 1.0)
LROD can detect overlap regions among long reads.
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


LROD
=================

1) Introduction
```
    LROD is a tool which aims to detect overlap regions among long reads..
    The input long read data of LROD is the long reads (fasta format).
```
2) Before installing and running
```
    Please build and install DSK which can count kmer frequency of long reads
```
3) Installing.
```
    LROD should run on Linux operating sysetm with gcc. We test LROD using gcc4.6.3 on Ubuntu.
    Create a main directory (eg:LROD). Copy all source code to this directory.
	cd LROD
	make all
```
4) Running.
```
    Step 1: Use DSK to crete the kmer frequency file.
        dsk -file <long-read-file.fa>  -kmer-size 15
	dsk2ascii -file <long-read-file.h5> -out kmer-frequency-file.txt
    Step 2: LROD -r <long-read-file> -c <kmer-frequency-file> -o result-file [options]
    	-r long-read-file: input file with fasta format;
	-c kmer-frequency-file: each line in the kmer-frequency-file should be "kmer kmer-frequency";
	-o result-file: result file;
	-t count: thread count (default 1);
	-k kmerLength: kmer length (default 15);
	-q smallKmerLength: small kmer length (default 9);
	-f minimumKmerFrequency: minimum kmer frequency (default 2);
	-m maxKmerFrequencyRatio: maximum kmer frequency ratio (default 0.9);
	-s kmerStep: kmer step (default 1);
	-d distance: the small distance used for determining whether two common kmers are consistant (default 400);
	-e distance: the large distance used for determining whether two common kmers are consistant (default 1500);
	-a min-overlap-length: the minimum overlap length between two long reads (default 500);
	-b length-ratio: the maximum length ratio between two aligned regions (default 0.3); 
	-h -help  Show rules of use for LROD.
	
    Note:
    	Each line in the kmer-frequency-file should be "kmer kmer-frequency". 
	For example:
	AGTCCAGGCCGGGAA 3
	GAAATCCAGCCGCCG 6
	AACCGGCGAATCGGA 3
	TATTTTAACATTCTC 2
	TATGGCCGATGAATT 4
	AAAGCCGAAGCCTAG 3
	CATCTTCACATCAGA 2
	ATAAGTGATAGCTTC 4
	TCGGCCATATTACCA 4
	ATTATTGCAATACTT 6
     An example of command line is shown below. The output file is result.csv.
	LROD -r sra.fasta -c kmer-cout -o result
```
5) Output.
```
    The output file "output-file-name.csv" is the overlap result.
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
    It means the region [5326,9923] in the first read is forward overlapped with the region [1,4364] in the second read.
```
6) Test.
```
    In the directory Test, there are long read file and its k-mer frequency file. Users can use this dataset to run LROD.
```
