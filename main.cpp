#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<iostream>
#include<ctype.h>
#include <time.h>
#include <malloc.h>
#include <string.h>


#include "kmer.h"
#include "read.h"
#include "aligning.h"
#include "bitarray.h"

using namespace std;


int main(int argc,char** argv)
{
	
	int maxSize = 100000;
	char * StrLine = (char *)malloc(sizeof(char)*maxSize);
	int kmerLength = atoi(argv[3]);
	int i,j;
	int readcount;
	
	int fileNameLen = strlen(argv[1]);
	char * readFile = (char *)malloc(sizeof(char)*fileNameLen + 10);
	strcpy(readFile, argv[1]);
	fileNameLen = strlen(argv[2]);
	char * binaryKmerFile = (char *)malloc(sizeof(char)*fileNameLen + 10);
	strcpy(binaryKmerFile, argv[2]);
	
	int step =  atoi(argv[4]);
	int Min = atoi(argv[5]);
	int Max = atoi(argv[6]);
	fileNameLen = strlen(argv[7]); 
	char * outputKmerFile = (char *)malloc(sizeof(char)*fileNameLen + 10);
	strcpy(outputKmerFile, argv[7]);
	
	ReadSetHead * readSetHead = GetReadSetHead(readFile,StrLine,maxSize);
	free(StrLine);
	
	readcount = readSetHead->readCount;
	int indexNumCut =  atoi(argv[3]);
	if(indexNumCut == 0 || indexNumCut>readcount){
		indexNumCut = readcount;
	}
	
	KmerHashTableHead * kmerHashTableHead = (KmerHashTableHead *)malloc(sizeof(KmerHashTableHead));
	
	KmerReadNodeHead * kmerReadNodeHead = InitKmerReadNodeHead(binaryKmerFile,readFile,kmerLength,step,Min,Max,kmerHashTableHead);
	
	GetCommonKmerHeadAllThread(kmerHashTableHead, kmerReadNodeHead, readSetHead, kmerLength, readFile, outputKmerFile, step, atoi(argv[8]));
		
	return 0;

}
