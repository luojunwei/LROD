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
	unsigned int startPositionInArray;
	//unsigned int fre;
}KmerHashNode;

typedef struct KmerHashTableHead{
    KmerHashNode * kmerHashNode;
	unsigned long int allocationCount;
}KmerHashTableHead;

typedef struct KmerReadNode{
    unsigned int kmer;
	unsigned int readIndex;
	unsigned int position;
	bool orientation;
}KmerReadNode;

typedef struct KmerReadNodeHead{
    KmerReadNode * kmerReadNode;
	long int realCount;
    long int allocationCount;
	long int kmerLength;
}KmerReadNodeHead;


char *  strup(char * str);

void ReverseComplementKmer(char * kmer, long int kmerLength);

unsigned int hash32shift(unsigned int key);

unsigned int Hash(unsigned int kmer, unsigned int max);

long int SearchKmerHashTable(KmerHashTableHead * kmerHashTableHead, unsigned int kmer);

void sort(KmerReadNodeHead * a, long int left, long int right);

bool DetectSameKmer(char * kmer, long int kmerLength);

KmerReadNodeHead * InitKmerReadNodeHead(char * address, ReadSetHead * readSetHead, long int kmerLength, long int step, KmerHashTableHead * kmerHashTableHead, int frequencyCutOff);



#endif