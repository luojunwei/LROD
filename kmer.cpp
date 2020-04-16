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
#include <time.h>

#include "kmer.h"
#include "bitarray.h"

char * strup(char * str)     
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


unsigned int Hash(unsigned int kmer, unsigned int max)  
{  
   return (hash32shift(kmer)*kmer) % max;  
} 




void ReverseComplementKmer(char * kmer, long int kmerLength){
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

	long int hashIndex = Hash(kmer, kmerHashTableHead->allocationCount);
	while(true){
		if(kmerHashTableHead->kmerHashNode[hashIndex].kmer == 0){
			return -1;;
		}
		if(kmerHashTableHead->kmerHashNode[hashIndex].kmer == kmer + 1){
			return hashIndex;
		}else{
			hashIndex = (hashIndex + 5)%kmerHashTableHead->allocationCount;
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
		
		if(i < j){
			a[i].kmer = a[j].kmer;
			a[i].readIndex = a[j].readIndex;
			a[i].position = a[j].position;
			a[i].orientation = a[j].orientation;
			i++;
		}
         
        while(i < j && key >= a[i].kmer){
            i++;
        }
		
		if(i < j){
			a[j].kmer = a[i].kmer;
			a[j].readIndex = a[i].readIndex;
			a[j].position = a[i].position;
			a[j].orientation = a[i].orientation;
			j--;
		}
         
    }
     
    a[i].kmer = key;
	a[i].readIndex = readIndex;
	a[i].position = position;
	a[i].orientation = orientation;

    sort(a, left, i - 1);
    sort(a, i + 1, right);
                   
}

bool DetectSameKmer(char * kmer, long int kmerLength){
	long int i = 1;
	for(; i < kmerLength; i++){
		if(kmer[0] != kmer[i]){
			break;
		}
	}
	if(i < kmerLength){
		return true;
	}
	return false;

}

KmerReadNodeHead * InitKmerReadNodeHead(char * address, ReadSetHead * readSetHead, long int kmerLength, long int step , KmerHashTableHead * kmerHashTableHead, int frequencyCutOff){

	FILE * fp; 
    if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	long int hashIndex = 0;
	long int maxSize = 100000;
	char * line = (char *)malloc(sizeof(char)*maxSize);
	char * kmer = (char *)malloc(sizeof(char)*(kmerLength + 1));
	char * kmerC = (char *)malloc(sizeof(char)*(10 + 1));
	long int kmerCount = 0;
	long int length = 0;
	int frequency = 0;
	
	long int arrayCount = 1000000;
	int * freArray = (int *)malloc(sizeof(int)*arrayCount);
	for(long int i = 0; i < arrayCount; i++){
		freArray[i] = 0;
	}
	cout<<"bb"<<endl;
	long int allKmerFrequency = 0;
	while((fgets(line, maxSize, fp)) != NULL){
		
		length = strlen(line);
		
		strncpy(kmer, line + kmerLength + 1, length - kmerLength - 1);
		kmer[length-kmerLength-1] = '\0';
		if(DetectSameKmer(kmer, kmerLength) != true){
			continue;
		}
		
		frequency = atoi(kmer);	
		if(frequency > arrayCount - 10){
			continue;
		}
		freArray[frequency]++;
		kmerCount++;
		allKmerFrequency = allKmerFrequency + frequency;
	}
	
	printf("The number of kmer types: %ld;\n",kmerCount);
	printf("Sum of frequencies for different kmer: %ld;\n",allKmerFrequency);
	
	float acc = 0;
	long int max = 0;
	long int min = 0;
	
	for(long int i = 0; i < arrayCount; i++){
		if(freArray[i] != 0){
			acc = acc + (float)(freArray[i]*i)/allKmerFrequency;
			if(acc > 0.9 && max == 0){
				max = i;
				break;
			}
		}
	}
	
	if(min < 2){
		min = 2;
	}
	min = 2;
	if(min > max){
		cout<<"min is larger than max!"<<endl;
		exit(0);
	}
	
	printf("The range of kmer frequency is:[%ld,%ld];\n",min,max);
	
	fclose(fp);

	kmerCount = 0;
	allKmerFrequency = 0;
	

    if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	while((fgets(line, maxSize, fp)) != NULL){
		strncpy(kmer, line, kmerLength);
		kmer[kmerLength]='\0';
		if(DetectSameKmer(kmer, kmerLength) != true){
			continue;
		}
		
		length = strlen(line);
		strncpy(kmerC, line + kmerLength + 1, length - kmerLength - 1);
		kmerC[length-kmerLength-1] = '\0';
		frequency = atoi(kmerC);
		
		if(frequency > arrayCount - 10){
			continue;
		}
		if(frequency <= max && frequency >= min){
			kmerCount++;
			allKmerFrequency = allKmerFrequency + frequency;
		}
	}
	
	printf("There are %ld kmer for overlap detection;\n",kmerCount);
	
	fclose(fp);
	
	kmerHashTableHead->allocationCount = kmerCount*1.2;
	kmerHashTableHead->kmerHashNode = (KmerHashNode *)malloc(sizeof(KmerHashNode)*kmerHashTableHead->allocationCount);
	
	for(unsigned long int i = 0; i < kmerHashTableHead->allocationCount; i++){
		kmerHashTableHead->kmerHashNode[i].kmer = 0;
		kmerHashTableHead->kmerHashNode[i].startPositionInArray = -1;
	}
	
	if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	
	unsigned long int kmerInteger = 0;

	while((fgets(line, maxSize, fp)) != NULL){

		strncpy(kmer, line, kmerLength);
		kmer[kmerLength]='\0';
		if(DetectSameKmer(kmer, kmerLength) != true){
			continue;
		}
		length = strlen(line);
		strncpy(kmerC, line + kmerLength + 1, length - kmerLength - 1);
		kmerC[length-kmerLength-1] = '\0';
		frequency = atoi(kmerC);	
		if(frequency > arrayCount - 10){
			continue;
		}
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
	cout<<"bb3333:"<<kmerReadNodeHead->allocationCount<<endl;

	long int readLength = 0;

	char * kmer1 = (char *)malloc(sizeof(char)*(kmerLength + 1));
	

	
	for(long int i = 0; i < readSetHead->readCount; i++){
		readLength = readSetHead->readSet[i].readLength;
		for(int j = 0; j < readLength - kmerLength + 1-step ; j = j+step){
			
			strncpy(kmer1,readSetHead->readSet[i].read + j, kmerLength);
			kmer1[kmerLength] = '\0';
			if(DetectSameKmer(kmer1, kmerLength) != true){
				continue;
			}
			SetBitKmer(&kmerInteger, kmerLength, kmer1);
			hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
			if(hashIndex!=-1){
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].kmer = kmerInteger;
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].readIndex = i+1; 
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].position= j; 
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].orientation = true;
				kmerReadNodeHead->realCount++;
			}else{
				ReverseComplementKmer(kmer1, kmerLength);
				SetBitKmer(&kmerInteger, kmerLength, kmer1);
				hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
				if(hashIndex!=-1){
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].kmer = kmerInteger;
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].readIndex = i+1; 
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
	free(line);
	free(kmer);
	free(kmerC);
	free(kmer1);
	kmerReadNodeHead->kmerLength = kmerLength;
	
	return kmerReadNodeHead;
}

KmerHashTableHead * GetKmerHashTableHead(char * address, ReadSetHead * readSetHead, long int kmerLength, long int step, long int min, float maxRatio){
	KmerHashTableHead * kmerHashTableHead = (KmerHashTableHead *)malloc(sizeof(KmerHashTableHead));
	FILE * fp; 
    if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	long int hashIndex = 0;
	long int maxSize = 100000;
	char * line = (char *)malloc(sizeof(char)*maxSize);
	char * kmer = (char *)malloc(sizeof(char)*(kmerLength + 1));
	char * kmerC = (char *)malloc(sizeof(char)*(10 + 1));
	long int kmerCount = 0;
	long int length = 0;
	int frequency = 0;
	
	long int arrayCount = 1000000;
	int * freArray = (int *)malloc(sizeof(int)*arrayCount);
	for(long int i = 0; i < arrayCount; i++){
		freArray[i] = 0;
	}
	
	long int allKmerFrequency = 0;
	while((fgets(line, maxSize, fp)) != NULL){
		
		length = strlen(line);
		
		strncpy(kmer, line + kmerLength + 1, length - kmerLength - 1);
		kmer[length-kmerLength-1] = '\0';
		if(DetectSameKmer(kmer, kmerLength) != true){
			continue;
		}
		
		frequency = atoi(kmer);	
		if(frequency > arrayCount - 10){
			continue;
		}
		freArray[frequency]++;
		kmerCount++;
		allKmerFrequency = allKmerFrequency + frequency;
	}
	
	printf("The number of kmer types: %ld;\n",kmerCount);
	printf("Sum of frequencies for different kmer: %ld;\n",allKmerFrequency);
	
	float acc = 0;
	long int max = 0;
	
	if(maxRatio > 1){
		maxRatio = 1;
	}
	
	for(long int i = 0; i < arrayCount; i++){
		if(freArray[i] != 0){
			acc = acc + (float)(freArray[i]*i)/allKmerFrequency;
			if(acc > maxRatio && max == 0){
				max = i;
				break;
			}
		}
	}
	
	if(min < 2){
		printf("minimumKmerFrequency < 2, LROD sets minimumKmerFrequency = 2!");
		min = 2;
	}

	if(min > max){	
		cout<<"The paramter minimumKmerFrequency is larger than maxKmerFrequencyRatio. Please increas the value of maxKmerFrequencyRatio or decrease the value of minimumKmerFrequency!"<<endl;
		exit(0);
	}
	
	printf("The range of kmer frequency is:[%ld,%ld];\n",min,max);
	
	fclose(fp);

	kmerCount = 0;
	allKmerFrequency = 0;
	

    if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	while((fgets(line, maxSize, fp)) != NULL){
		strncpy(kmer, line, kmerLength);
		kmer[kmerLength]='\0';
		if(DetectSameKmer(kmer, kmerLength) != true){
			continue;
		}
		
		length = strlen(line);
		strncpy(kmerC, line + kmerLength + 1, length - kmerLength - 1);
		kmerC[length-kmerLength-1] = '\0';
		frequency = atoi(kmerC);
		
		if(frequency > arrayCount - 10){
			continue;
		}
		if(frequency <= max && frequency >= min){
			kmerCount++;
			allKmerFrequency = allKmerFrequency + frequency;
		}
	}
	
	printf("There are %ld kmer for overlap detection;\n",kmerCount);
	
	fclose(fp);
	
	kmerHashTableHead->allocationCount = kmerCount*2;
	kmerHashTableHead->kmerHashNode = (KmerHashNode *)malloc(sizeof(KmerHashNode)*kmerHashTableHead->allocationCount);
	
	for(unsigned long int i = 0; i < kmerHashTableHead->allocationCount; i++){
		kmerHashTableHead->kmerHashNode[i].kmer = 0;
		kmerHashTableHead->kmerHashNode[i].startPositionInArray = -1;
	}
	
	if((fp = fopen(address, "r")) == NULL){
        printf("%s, does not exist!", address);
        exit(0);
    }
	
	unsigned long int kmerInteger = 0;

	while((fgets(line, maxSize, fp)) != NULL){

		strncpy(kmer, line, kmerLength);
		kmer[kmerLength]='\0';
		if(DetectSameKmer(kmer, kmerLength) != true){
			continue;
		}
		length = strlen(line);
		strncpy(kmerC, line + kmerLength + 1, length - kmerLength - 1);
		kmerC[length-kmerLength-1] = '\0';
		frequency = atoi(kmerC);	
		if(frequency > arrayCount - 10){
			continue;
		}
		if(!(frequency <= max && frequency >= min)){
			continue;
		}
		//cout<<line<<endl;
		
		SetBitKmer(&kmerInteger, kmerLength, kmer);
		//cout<<kmerInteger<<endl;
		hashIndex = Hash(kmerInteger, kmerHashTableHead->allocationCount);
		while(true){
			if(kmerHashTableHead->kmerHashNode[hashIndex].kmer == 0){
				kmerHashTableHead->kmerHashNode[hashIndex].kmer = kmerInteger + 1;
				break;
			}else{
				hashIndex = (hashIndex + 5) % kmerHashTableHead->allocationCount;
			}
		}
	}
	fclose(fp);
	
	free(line);
	free(kmer);
	free(kmerC);

	
	return kmerHashTableHead;
}




KmerReadNodeHead * GetKmerReadNodeHeadSub(ReadSetHead * readSetHead, long int kmerLength, long int step, long int intervalCount){
	
	
	long int max = 0;
	long int allKmerFrequency = 0;
	long int startReadIndex = 0;
	long int endReadIndex = startReadIndex + intervalCount;
	while(true){
		if(startReadIndex >= readSetHead->readCount){
			break;
		}
		if(endReadIndex >= readSetHead->readCount){
			endReadIndex = readSetHead->readCount - 1;
		}
		allKmerFrequency = 0;
		for(long int i = startReadIndex; i <= endReadIndex; i++){
			allKmerFrequency = allKmerFrequency + readSetHead->readSet[i].readLength - kmerLength;
		}
		if(allKmerFrequency > max){
			max = allKmerFrequency;
		}
		startReadIndex = endReadIndex + 1;
		endReadIndex = endReadIndex + intervalCount;
	}
	max = max/step;
	
	
	KmerReadNodeHead * kmerReadNodeHead = (KmerReadNodeHead *)malloc(sizeof(KmerReadNodeHead));
	kmerReadNodeHead->realCount = 0;
	kmerReadNodeHead->allocationCount = max*1.1;
	kmerReadNodeHead->kmerReadNode = (KmerReadNode *)malloc(sizeof(KmerReadNode)*kmerReadNodeHead->allocationCount);
	
	for(long int i = 0; i < kmerReadNodeHead->allocationCount; i++){
		kmerReadNodeHead->kmerReadNode[i].kmer = 0;
	}
	
	return kmerReadNodeHead;
}

void InitKmerReadNodeHeadSub(ReadSetHead * readSetHead, KmerReadNodeHead * kmerReadNodeHead, KmerHashTableHead * kmerHashTableHead, long int kmerLength, long int step, long int startReadIndex, long int endReadIndex){
	kmerReadNodeHead->realCount = 0;
	for(unsigned long int i = 0; i < kmerHashTableHead->allocationCount; i++){
		kmerHashTableHead->kmerHashNode[i].startPositionInArray = -1;
	}
	
	for(long int i = 0; i < kmerReadNodeHead->allocationCount; i++){
		kmerReadNodeHead->kmerReadNode[i].kmer = 0;
	}
	
	long int readLength = 0;

	char * kmer1 = (char *)malloc(sizeof(char)*(kmerLength + 1));
	
	long int hashIndex = 0;
	unsigned long int kmerInteger = 0;
	
	for(long int i = startReadIndex; i <= endReadIndex; i++){
		readLength = readSetHead->readSet[i].readLength;
		for(int j = 0; j < readLength - kmerLength + 1-step ; j = j+step){
			
			strncpy(kmer1,readSetHead->readSet[i].read + j, kmerLength);
			kmer1[kmerLength] = '\0';
			
			SetBitKmer(&kmerInteger, kmerLength, kmer1);
			hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
			if(hashIndex!=-1){
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].kmer = kmerInteger;
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].readIndex = i+1; 
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].position= j; 
				kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].orientation = true;
				kmerReadNodeHead->realCount++;
			}else{
				ReverseComplementKmer(kmer1, kmerLength);
				SetBitKmer(&kmerInteger, kmerLength, kmer1);
				hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
				if(hashIndex!=-1){
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].kmer = kmerInteger;
					kmerReadNodeHead->kmerReadNode[kmerReadNodeHead->realCount].readIndex = i+1; 
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

	kmerReadNodeHead->kmerLength = kmerLength;
	kmerReadNodeHead->startReadIndex = startReadIndex;
	kmerReadNodeHead->endReadIndex = endReadIndex;
	
}


#endif