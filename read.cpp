#ifndef Read_CPP_INCLUDED 
#define Read_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "read.h"

using namespace std;

long int max (long int x, long int y)
{
	long int z;
	if (x > y) z = x;
	else z = y;
	return (z);
}

long int min (long int x, long int y)
{
	long int z;
	if (x < y) z = x;
	else z = y;
	return (z);
}

ReadSetHead * GetReadSetHead(char *filename,char *StrLine, long int maxSize){
	
	ReadSetHead * readSetHead = (ReadSetHead * )malloc(sizeof(ReadSetHead));
	
	FILE *fp; 
	long int m=1;
    if((fp = fopen(filename,"r")) == NULL){ 
        printf("error!"); 
        return NULL; 
    } 
	
    readSetHead->readCount = 0;
	while((fgets(StrLine, maxSize, fp)) != NULL){	
		if(StrLine[0]=='>'){
			readSetHead->readCount++;
		}
	}
	fclose(fp);
	
	printf("Number of long reads: %ld;\n",readSetHead->readCount);
	
	readSetHead->readSet = (ReadSet *)malloc(sizeof(ReadSet)*readSetHead->readCount);
	for(long int i = 0; i <readSetHead->readCount; i++){
		readSetHead->readSet[i].readLength = 0;
	}
	
	FILE *fp1; 
	if((fp1 = fopen(filename,"r")) == NULL){ 
        printf("error!"); 
        return NULL; 
    } 
	long int readIndex = -1;
	long int j = 0;
	long int allReadLen = 0;
	while((fgets(StrLine, maxSize, fp1)) != NULL){

		j++;

		if(StrLine[0]=='>'){
			readIndex++;
			continue;
		}
		
		readSetHead->readSet[readIndex].readLength = strlen(StrLine);
		if(StrLine[readSetHead->readSet[readIndex].readLength - 1]=='\n' || StrLine[readSetHead->readSet[readIndex].readLength - 1]=='\r'){
			readSetHead->readSet[readIndex].readLength--;
		}
		
		readSetHead->readSet[readIndex].read = (char *)malloc(sizeof(char)*(readSetHead->readSet[readIndex].readLength + 1));
		strncpy(readSetHead->readSet[readIndex].read, StrLine, readSetHead->readSet[readIndex].readLength);
		readSetHead->readSet[readIndex].read[readSetHead->readSet[readIndex].readLength] = '\0';
		
	}
	
	return readSetHead;
}


#endif