#ifndef Read_H_INCLUDED 
#define Read_H_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

using namespace std;

typedef struct ReadSet{
	char * read;
	long int readLength;
//	long int allReadLen;
	
}ReadSet;

typedef struct ReadSetHead{
	ReadSet * readSet;
	long int readCount;
}ReadSetHead;


typedef struct ReadToKmerIndex{
	long int * index;
	long int indexCount;
}ReadToKmerIndex;

typedef struct ReadToKmerIndexSet{
	ReadToKmerIndex * readToKmerIndex;
	long int readCount;
}ReadToKmerIndexSet;


ReadSetHead * GetReadSetHead(char *filename,char *StrLine, long int maxSize);

long int max (long int x, long int y);

long int min (long int x, long int y);

#endif