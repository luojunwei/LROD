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

ReadSetHead * GetReadSetHead(char *filename,char *StrLine, int maxSize){
	
	ReadSetHead * readSetHead = (ReadSetHead * )malloc(sizeof(ReadSetHead));
	
	FILE *fp; 
	int m=1;
    if((fp = fopen(filename,"r")) == NULL) 
    { 
        printf("error!"); 
        return NULL; 
    } 
	
    readSetHead->readCount = 0;
	while(fgets(StrLine,maxSize,fp)){
		if(StrLine[0]=='\n' || StrLine[0]=='\r'){break;}
		if(feof(fp)){break;}
		
		if(StrLine[0]=='>'){
			readSetHead->readCount++;
		}
	}
	fclose(fp);
	
	readSetHead->readSet = (ReadSet *)malloc(sizeof(ReadSet)*readSetHead->readCount);
	for(int i = 0; i <readSetHead->readCount; i++){
		readSetHead->readSet[i].readLength = 0;
	}
	
	FILE *fp1; 
	if((fp1 = fopen(filename,"r")) == NULL) 
    { 
        printf("error!"); 
        return NULL; 
    } 
	int readIndex = -1;
	int j = 0;
	long int allReadLen = 0;

	while(fgets(StrLine,maxSize,fp1)){
		if(StrLine[0]=='\n' || StrLine[0]=='\r'){break;}
		if(feof(fp1)){break;}
		j++;
		if(j%2 == 1){
			readIndex++;
			continue;
		}
		readSetHead->readSet[readIndex].readLength = strlen(StrLine);
	}

	return readSetHead;
}


#endif