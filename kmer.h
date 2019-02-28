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
//#include "FindNumber.h"
#include "bitarray.h"
using namespace std;

typedef struct KmerHashNode{
    unsigned int kmer;
	unsigned long int startPositionInArray;
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
	int kmerLength;
}KmerReadNodeHead;


char *  strup(char * str);

void ReverseComplementKmer(char * kmer, int kmerLength);

unsigned int hash32shift(unsigned int key);

unsigned long int Hash(unsigned int kmer, unsigned long int max);

long int SearchKmerHashTable(KmerHashTableHead * kmerHashTableHead, unsigned int kmer);

void sort(KmerReadNodeHead * a, long int left, long int right);

KmerReadNodeHead * InitKmerReadNodeHead(char * address, char * file, int kmerLength, int step, int min, int max, KmerHashTableHead * kmerHashTableHead);







#endif