#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<iostream>
#include<ctype.h>
#include <time.h>
#include <malloc.h>

#include <getopt.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>


#include "kmer.h"
#include "read.h"
#include "aligning.h"
#include "bitarray.h"

using namespace std;


int main(int argc,char** argv)
{
	
	long int maxSize = 1000000;
	char * StrLine = (char *)malloc(sizeof(char)*maxSize);

	char * readFile = NULL;
	char * binaryKmerFile = NULL;
	char * outputKmerFile = NULL;
	long int step = 1;
	long int kmerLength = 13;
	long int smallKmerLength = 9;
	long int threadCount = 1;
	
	long int smallIntervalDistance = 400;
	long int largeIntervalDistance = 1500;
	long int overlapLengthCutOff = 500;
	float lengthRatio = 0.3;
	int frequencyCutOff = 3;
	
	
	int ch = 0;
	while ((ch = getopt(argc, argv, "c:r:o:m:n:d:k:e:a:s:t:f:q:b:")) != -1) {
		switch (ch) {
			case 'r': readFile = (char *)(optarg); break;
			case 'c': binaryKmerFile = (char *)(optarg); break;
			case 'o': outputKmerFile = (char *)optarg; break;
			case 'k': kmerLength = atoi(optarg); break;
			case 'q': smallKmerLength = atoi(optarg); break;
			case 's': step = atoi(optarg); break;
			
			case 'd': smallIntervalDistance = atoi(optarg); break;
			case 'e': largeIntervalDistance = atoi(optarg); break;
			case 'a': overlapLengthCutOff = atoi(optarg); break;
			case 't': threadCount = atoi(optarg); break;
			case 'b': lengthRatio = atof(optarg); break;

			default: break; 
		}
	}
	
	ReadSetHead * readSetHead = GetReadSetHead(readFile,StrLine,maxSize);
	
	KmerHashTableHead * kmerHashTableHead = (KmerHashTableHead *)malloc(sizeof(KmerHashTableHead));

	KmerReadNodeHead * kmerReadNodeHead = InitKmerReadNodeHead(binaryKmerFile,readSetHead,kmerLength,step,kmerHashTableHead,frequencyCutOff);
	
	GetCommonKmerHeadAllThread(kmerHashTableHead, kmerReadNodeHead, readSetHead, kmerLength, readFile, outputKmerFile, step, threadCount, smallKmerLength, smallIntervalDistance, largeIntervalDistance, overlapLengthCutOff, lengthRatio);

	return 0;

}
