#ifndef Aligning_H_INCLUDED 
#define Aligning_H_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "read.h"
#include "kmer.h"
#include "bitarray.h"

using namespace std;

typedef struct CommonKmer{
	unsigned long  int readIndex;
	unsigned long  int leftPosition;
	unsigned long  int rightPosition;
	bool orientation;
}CommonKmer;

typedef struct CommonKmerHead{
	CommonKmer * commonKmer;
	unsigned long int readIndex;
	unsigned long int realCount;
	unsigned long int allocationCount;
}CommonKmerHead;

typedef struct ArcIndex{
	unsigned long  int startIndex;
	unsigned long  int endIndex;
}ArcIndex;
 
typedef struct AdjGraph{
	unsigned long  int dataLeft;
	unsigned long  int dataRight;
	unsigned long  int eadgCount;
	bool visit;
}AdjGraph;
 
typedef struct AdjGraphHead{
	AdjGraph * graph;
	unsigned long  int realCountGraph;
	unsigned long  int allocationCountGraph;
	AdjGraph * reverseGraph;
	unsigned long  int reverseRealCountGraph;
	unsigned long  int reverseAllocationCountGraph;
	ArcIndex * arcIndex;
	unsigned long  int realCountArc;
	unsigned long  int allocationCountArc;
	int * distanceToSource;
	int * edgeToSource;
	unsigned long  int nodeCount;
	long int leftStart;
	long int rightStart;
	long int leftEnd;
	long int rightEnd;
}AdjGraphHead;

typedef struct GetCommonKmerHeadP{
	KmerHashTableHead * kmerHashTableHead;
	KmerReadNodeHead * kmerReadNodeHead;
	ReadSetHead * readSetHead;
	int kmerLength;
	char * readFile;
	char * outputFile;
	unsigned long  int step;
	int threadIndex;
	int totalThreadNumber;
}GetCommonKmerHeadP;

void ReAllocateCommonKmer(CommonKmerHead * commonKmerHead);

void InsertCommonToTwoReadAligningHead(CommonKmerHead * commonKmerHead, KmerReadNodeHead * kmerReadNodeHead, KmerHashTableHead * kmerHashTableHead, long int hashIndex, unsigned long int readIndex);

CommonKmerHead * GetCommonKmerHead(KmerHashTableHead * kmerHashTableHead, KmerReadNodeHead * kmerReadNodeHead, ReadSetHead * readSetHead, int kmerLength, char * readFile,  char * outputFile, unsigned long  int step);

void sort(CommonKmer * a, int left, int right);

void sortGraph(AdjGraph * graph, int left, int right);

void DestroyGraph(AdjGraphHead * G);

int Overlap_Display(AdjGraphHead * G,int leftIndex,int rightIndex,bool orien,int leftLen,int rightLen,FILE * fp);

double GetLongestPathInGraph(AdjGraphHead * G, bool orien);

int CreatGraph(AdjGraphHead * G, CommonKmer * commonKmer, unsigned long  int startIndex, unsigned long  int endIndex, int a);

void GetOverlapResult(AdjGraphHead * G, CommonKmerHead * commonKmerHead, ReadSetHead * readSetHead,FILE * fp);

void ReAllocateAdjGraph(AdjGraphHead * G, int a);

void ReAllocateArcIndex(AdjGraphHead * G);

void * GetCommonKmerHeadThread(void * arg);

CommonKmerHead * GetCommonKmerHeadAllThread(KmerHashTableHead * kmerHashTableHead, KmerReadNodeHead * kmerReadNodeHead, ReadSetHead * readSetHead, int kmerLength, char * readFile, char * outputFile, unsigned long  int step, int totalThreadNumber);

#endif