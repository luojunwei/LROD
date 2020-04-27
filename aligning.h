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
	long int readIndex;
	long int leftPosition;
	long int rightPosition;
	bool orientation;
}CommonKmer;

typedef struct CommonKmerHead{
	CommonKmer * commonKmer;
	long int readIndex;
	long int realCount;
	long int allocationCount;
}CommonKmerHead;

typedef struct CommonKmerHeadArray{
	CommonKmerHead * commonKmerHead;
	long int count;
}CommonKmerHeadArray;

typedef struct ArcIndex{
	long int startIndex;
	long int endIndex;
	float weight;
}ArcIndex;
 
typedef struct AdjGraph{
	long int dataLeft;
	long int dataRight;
	bool visit;
}AdjGraph;
 
typedef struct AdjGraphHead{
	AdjGraph * graph;
	long int realCountGraph;
	long int allocationCountGraph;
	AdjGraph * reverseGraph;
	long int reverseRealCountGraph;
	long int reverseAllocationCountGraph;
	ArcIndex * arcIndex;
	long int realCountArc;
	long int allocationCountArc;
	float * distanceToSource;
	long int * edgeToSource;
	bool * visited;
	long int largestIntervalDistance;
	long int nodeCount;
	long int leftStart;
	long int rightStart;
	long int leftEnd;
	long int rightEnd;
	char * localLeftRead;
	char * localRightRead;
	long int smallKmerLength;
	long int kmerLength;
	float lengthRatio;
	long int overlapLengthCutOff;
}AdjGraphHead;

typedef struct GetCommonKmerHeadP{
	KmerHashTableHead * kmerHashTableHead;
	KmerReadNodeHead * kmerReadNodeHead;
	ReadSetHead * readSetHead;
	long int kmerLength;
	char * readFile;
	char * outputFile;
	long int step;
	long int threadIndex;
	long int totalThreadNumber;
	long int smallIntervalDistance;
	long int largeIntervalDistance;
	long int overlapLengthCutOff;
	long int smallKmerLength;
	float lengthRatio;
	long int startReadIndex;
}GetCommonKmerHeadP;

void ReAllocateCommonKmer(CommonKmerHead * commonKmerHead);

void InsertCommonToTwoReadAligningHead(CommonKmerHead * commonKmerHead, KmerReadNodeHead * kmerReadNodeHead, KmerHashTableHead * kmerHashTableHead, long int hashIndex, unsigned long int readIndex, unsigned long  long int position, bool orien);

int GetCommonKmerHeadAllThreadNew_UnitTest(char * result, long int readCount);

CommonKmerHead * GetCommonKmerHead(KmerHashTableHead * kmerHashTableHead, KmerReadNodeHead * kmerReadNodeHead, ReadSetHead * readSetHead, long int kmerLength, char * readFile,  char * outputFile, unsigned long  long int step);

void sort(CommonKmer * a, long int left, long int right);

void sortGraph(AdjGraph * graph, long int left, long int right);

void DestroyGraph(AdjGraphHead * G);

long int Overlap_Display(AdjGraphHead * G,long int leftIndex,long int rightIndex,bool orien,long int leftLen,long int rightLen,FILE * fp, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, char * leftRead, char * rightRead);

long int Overlap_DisplayLocalRegion(AdjGraphHead * G,long int leftLen,long int rightLen);

long int AddEdgeInGraph(AdjGraphHead * G, bool orientation, long int leftIndex , long int rightIndex, long int largestIntervalDistance, long int maxIntervalDistance, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, char * leftRead, char * rightRead);

long int CreatGraphSinglePath(AdjGraphHead * G, CommonKmer * commonKmer, unsigned long int startIndex, unsigned long int endIndex, long int a, long int leftIndex, long int rightIndex, long int leftLen, long int rightLen, FILE * fp, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, char * leftRead, char * rightRead);

long int CreatGraphLocalRegion(AdjGraphHead * G, long int distance);

void GetOverlapResult(AdjGraphHead * G, CommonKmerHead * commonKmerHead, ReadSetHead * readSetHead, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, FILE * fp);

void ReAllocateAdjGraph(AdjGraphHead * G, long int a);

void ReAllocateArcIndex(AdjGraphHead * G);

void * GetCommonKmerHeadThread(void * arg);

CommonKmerHead * GetCommonKmerHeadAllThread(KmerHashTableHead * kmerHashTableHead, KmerReadNodeHead * kmerReadNodeHead, ReadSetHead * readSetHead, long int kmerLength, char * readFile, char * outputFile, unsigned long  long int step, long int totalThreadNumber, long int smallKmerLength, long int smallIntervalDistance, long int largeIntervalDistance, long int overlapLengthCutOff, float lengthRatio, long int startReadIndex);


CommonKmerHead * GetCommonKmerHeadAllThreadNew(KmerHashTableHead * kmerHashTableHead, KmerReadNodeHead * kmerReadNodeHead, ReadSetHead * readSetHead, long int kmerLength, char * readFile, char * outputFile, unsigned long  long int step, long int totalThreadNumber, long int smallKmerLength, long int smallIntervalDistance, long int largeIntervalDistance, long int overlapLengthCutOff, float lengthRatio, long int subReadCount);


void DetectCommon(CommonKmerHead * commonKmerHead, long int position, char * kmer, char * read, long int readLength, long int kmerLength, long int distance);

long int GetCommonShorterKmer(AdjGraphHead * G, CommonKmerHead * commonKmerHead, char * leftRead, char * rightRead, long int leftStartPosition, long int leftEndPosition, long int rightStartPosition, long int rightEndPosition, long int kmerLength, bool orientation);

void swapCommonKmer(CommonKmer *a, long int left, long int right);

void downToMaxHeap(CommonKmer *a, long int bgn, long int end);

void BubbleSort(CommonKmer *a, long int left, long int right);

void buildMaxHeap(CommonKmer *a, long int bgn, long int end);

void heapSort(CommonKmer *a, long int bgn, long int end);

void SubRemoveMultipleSameKmer(CommonKmerHead * commonKmerHead, long int startIndex, long int endIndex);

void RemoveMultipleSameKmer(CommonKmerHead * commonKmerHead);

void RemoveLowNumberKmer(CommonKmerHead * commonKmerHead, long int * forwardKmerCount, long int * reverseKmerCount, long int readCount);

#endif