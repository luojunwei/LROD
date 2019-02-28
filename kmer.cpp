#ifndef Kmer_CPP_INCLUDED 
#define Kmer_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
//#include <stdint.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include<assert.h> 

#include "kmer.h"
#include "bitarray.h"

char *  strup(char * str)     
{  
    assert(str);                  
    char *ret = str;              
    while(*str != '\0')           
    {     
           if((*str >= 'a')&&(*str <= 'z'))
        {  
            *str = *str -32;          
            str++;  
        }  
        else  
            str++;  
    }  
    return ret;             
}  

unsigned int hash32shift(unsigned int key) 
{ 
  key = ~key + (key << 15); 
  key = key ^ (key >> 12); 
  key = key + (key << 2); 
  key = key ^ (key >> 4); 
  key = key * 2057; 
  key = key ^ (key >> 16); 
  return key; 
}


unsigned long int Hash(unsigned int kmer, unsigned long int max)  
{  
   return (hash32shift(kmer)*kmer) % max;  
} 

void ReverseComplementKmer(char * kmer, int kmerLength){
	for(int k = 0; k < kmerLength/2; k++){
		char temp = kmer[k];
		kmer[k] = kmer[kmerLength -1 - k];
		kmer[kmerLength -1 - k] = temp;
	}
	
	for(int i = 0; i < kmerLength; i++){
		if(kmer[i] == 'A' || kmer[i] == 'a'){
			kmer[i] = 'T';
		}else if(kmer[i] == 'T' || kmer[i] == 't'){
			kmer[i] = 'A';
		}else if(kmer[i] == 'G' || kmer[i] == 'g'){
			kmer[i] = 'C';
		}else if(kmer[i] == 'C' || kmer[i] == 'c'){
			kmer[i] = 'G';
		}else if(kmer[i] == 'N' || kmer[i] == 'n'){
			kmer[i] = 'N';
		}
	}
}


long int SearchKmerHashTable(KmerHashTableHead * kmerHashTableHead, unsigned int kmer){

	unsigned long int hashIndex = Hash(kmer, kmerHashTableHead->allocationCount);
	while(true){
		if(kmerHashTableHead->kmerHashNode[hashIndex].kmer == 0){
			return -1;;
		}
		if(kmerHashTableHead->kmerHashNode[hashIndex].kmer == kmer + 1){
			return hashIndex;
		}else{
			hashIndex = (hashIndex + 1)%kmerHashTableHead->allocationCount;
		}
	}
	return -1;
}


void sort(KmerReadNode * a, long int left, long int right)
{
    if(left >= right){
        return ;
    }
    long int i = left;
    long int j = right;
    unsigned int key = a[left].kmer;
	unsigned int position = a[left].position;
	unsigned int readIndex = a[left].readIndex;
	bool orientation = a[left].orientation;
	
    while(i < j)                       
    {
        while(i < j && key <= a[j].kmer){
            j--;
        }
         
        a[i].kmer = a[j].kmer;
		a[i].readIndex = a[j].readIndex;
		a[i].position = a[j].position;
		a[i].orientation = a[j].orientation;
        
         
        while(i < j && key >= a[i].kmer){
            i++;
        }
         
        a[j].kmer = a[i].kmer;
		a[j].readIndex = a[i].readIndex;
		a[j].position = a[i].position;
		a[j].orientation = a[i].orientation;
    }
     
    a[i].kmer = key;
	a[i].readIndex = readIndex;
	a[i].position = position;
	a[i].orientation = orientation;

    sort(a, left, i - 1);
    sort(a, i + 1, right);
                   
}


KmerReadNodeHead * InitKmerReadNodeHead(char * address, char * file,int kmerLength,int step ,int min, int max, KmerHashTableHead * kmerHashTableHead){

	FILE * fp; 
    if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	long int hashIndex = 0;
	int maxSize = 1000;
	char * line = (char *)malloc(sizeof(char)*maxSize);
	char * kmer = (char *)malloc(sizeof(char)*(kmerLength + 1));
	char * kmerC = (char *)malloc(sizeof(char)*(10 + 1));
	long int kmerCount = 0;
	long int length = 0;
	int frequency = 0;
	long int averragefrequency=0;
	long int allKmerFrequency = 0;
	while((fgets(line, maxSize, fp)) != NULL){
		length = strlen(line);
		strncpy(kmer, line + kmerLength + 1, length - kmerLength - 1);
		kmer[length-kmerLength-1] = '\0';
		frequency = atoi(kmer);	
		if(frequency <= max && frequency >= min){
			kmerCount++;
			allKmerFrequency = allKmerFrequency + frequency;
		}
	}
	
	fclose(fp);

	kmerHashTableHead->allocationCount = kmerCount*1.5;
	kmerHashTableHead->kmerHashNode = (KmerHashNode *)malloc(sizeof(KmerHashNode)*kmerHashTableHead->allocationCount);
	
	for(long int i = 0; i < kmerHashTableHead->allocationCount; i++){
		kmerHashTableHead->kmerHashNode[i].kmer = 0;
		kmerHashTableHead->kmerHashNode[i].startPositionInArray = 0;
	}
	
	if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	
	unsigned int kmerInteger = 0;
	
	while((fgets(line, maxSize, fp)) != NULL){
		
		strncpy(kmer, line, kmerLength);
		kmer[kmerLength]='\0';
		length = strlen(line);
		strncpy(kmerC, line + kmerLength + 1, length - kmerLength - 1);
		kmerC[length-kmerLength-1] = '\0';
		frequency = atoi(kmerC);	
		if(!(frequency <= max && frequency >= min)){
			continue;
		}
		
		SetBitKmer(&kmerInteger, kmerLength, kmer);
		
		hashIndex = Hash(kmerInteger, kmerHashTableHead->allocationCount);
		while(true){
			if(kmerHashTableHead->kmerHashNode[hashIndex].kmer == 0){
				kmerHashTableHead->kmerHashNode[hashIndex].kmer = kmerInteger + 1;
				break;
			}else{
				hashIndex = (hashIndex + 1) % kmerHashTableHead->allocationCount;
			}
		}
	}
	fclose(fp);

	
	
	KmerReadNodeHead * kmerReadNodeHead = (KmerReadNodeHead *)malloc(sizeof(KmerReadNodeHead));
	kmerReadNodeHead->realCount = 0;
	kmerReadNodeHead->allocationCount = allKmerFrequency*1.1;
	kmerReadNodeHead->kmerReadNode = (KmerReadNode *)malloc(sizeof(KmerReadNode)*kmerReadNodeHead->allocationCount);
	
	for(long int i = 0; i < kmerReadNodeHead->allocationCount; i++){
		kmerReadNodeHead->kmerReadNode[i].kmer = 0;
	}
	
	if((fp = fopen(file, "r")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }
	
	unsigned int readIndex = 0;
	long int Size = 100000;
	long int readLength = 0;
	char str[Size];
	char * kmer1 = (char *)malloc(sizeof(char)*(kmerLength + 1));
	char * kmer2 = (char *)malloc(sizeof(char)*(kmerLength + 1));
	while(!feof(fp)){
		fgets(str, Size, fp);
		if(str[0]=='>'){
			readIndex++;
			continue;
		}
		
		readLength = strlen(str);
		for(int j = 0; j < readLength - kmerLength + 1;j = j+step){
			strncpy(kmer1,str + j, kmerLength);
			kmer1[kmerLength] = '\0';
			SetBitKmer(&kmerInteger, kmerLength, kmer1);
			hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
			if(hashIndex!=-1){
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].kmer = kmerInteger;
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].readIndex = readIndex; 
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].position= j; 
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].orientation = true;
				kmerReadNodeHead->realCount++;
			}else{
				ReverseComplementKmer(kmer1, kmerLength);
				SetBitKmer(&kmerInteger, kmerLength, kmer1);
				hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
				if(hashIndex!=-1){
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].kmer = kmerInteger;
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].readIndex = readIndex; 
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].position= j; 
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].orientation = false;
					kmerReadNodeHead->realCount++;
				}
			}	
		}
	}
	
	sort(kmerReadNodeHead->kmerReadNode, 0, kmerReadNodeHead->realCount - 1);
	kmerInteger = kmerReadNodeHead->kmerReadNode[0].kmer + 1;
	for(long int i = 0; i < kmerReadNodeHead->realCount; i++){
		if(kmerReadNodeHead->kmerReadNode[i].kmer != kmerInteger){
			kmerInteger = kmerReadNodeHead->kmerReadNode[i].kmer;
			hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
			kmerHashTableHead->kmerHashNode[hashIndex].startPositionInArray = i;
		}
	}
	
	fclose(fp);
	free(line);
	free(kmer);
	free(kmerC);
	free(kmer1);

	return kmerReadNodeHead;
}


#endif