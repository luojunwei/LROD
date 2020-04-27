#ifndef Kmer_H_INCLUDED 
#define Kmer_H_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
//#include <stdint.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "read.h"
#include "bitarray.h"
using namespace std;

typedef struct KmerHashNode{
    unsigned int kmer;
	int startPositionInArray;
	//unsigned int fre;
}KmerHashNode;

typedef struct KmerHashTableHead{
    KmerHashNode * kmerHashNode;
	unsigned long int allocationCount;
	long int min;
	long int max;
}KmerHashTableHead;

typedef struct KmerReadNode{
    unsigned int kmer;
	unsigned int readIndex; // Indicates the index containing the read of the kmer
	unsigned int position;  //Indicates the position of the kmer in the read
	bool orientation;  //Indicates the direction of the kmer in the read
}KmerReadNode;

typedef struct KmerReadNodeHead{
    KmerReadNode * kmerReadNode;
	long int realCount;
    long int allocationCount;
	long int kmerLength;
	long int startReadIndex;  
	long int endReadIndex;
}KmerReadNodeHead;


char *  strup(char * str);

void ReverseComplementKmer(char * kmer, long int kmerLength);   //Inverse complementary function of kmer

unsigned int hash32shift(unsigned int key);

unsigned int Hash(unsigned int kmer, unsigned int max);
//This function is used to check whether a certain kmer exists in the kmer hash table
long int SearchKmerHashTable(KmerHashTableHead * kmerHashTableHead, unsigned int kmer);  

void sort(KmerReadNodeHead * a, long int left, long int right);

bool DetectSameKmer(char * kmer, long int kmerLength);  

KmerReadNodeHead * InitKmerReadNodeHead(char * address, ReadSetHead * readSetHead, long int kmerLength, long int step, KmerHashTableHead * kmerHashTableHead, int frequencyCutOff);

KmerHashTableHead * GetKmerHashTableHead(char * address, ReadSetHead * readSetHead, long int kmerLength, long int step, long int min, float maxRatio);

int GetKmerHashTableHead_UnitTest(KmerHashTableHead * kmerHashTableHead);

KmerReadNodeHead * GetKmerReadNodeHeadSub(ReadSetHead * readSetHead, long int kmerLength, long int step, long int intervalCount);

int GetKmerReadNodeHeadSub_UnitTest(KmerReadNodeHead * kmerReadNodeHead);

void InitKmerReadNodeHeadSub(ReadSetHead * readSetHead, KmerReadNodeHead * kmerReadNodeHead, KmerHashTableHead * kmerHashTableHead, long int kmerLength, long int step, long int startReadIndex, long int endReadIndex);

#endif