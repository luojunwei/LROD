#ifndef Aligning_CPP_INCLUDED 
#define Aligning_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <ctype.h>
#include <time.h>

#include "aligning.h"
#include "read.h"

using namespace std;

void ReAllocateCommonKmer(CommonKmerHead * commonKmerHead){
	unsigned long int maxCount = commonKmerHead->allocationCount*1.5;
	CommonKmer * commonKmer = (CommonKmer *)malloc(sizeof(CommonKmer)*maxCount);
	for(long int i = 0; i < commonKmerHead->allocationCount; i++){
		commonKmer[i].readIndex = commonKmerHead->commonKmer[i].readIndex;
		commonKmer[i].leftPosition = commonKmerHead->commonKmer[i].leftPosition;
		commonKmer[i].rightPosition = commonKmerHead->commonKmer[i].rightPosition;
		commonKmer[i].orientation = commonKmerHead->commonKmer[i].orientation;
	}
	free(commonKmerHead->commonKmer);
	commonKmerHead->commonKmer = commonKmer;
	commonKmerHead->allocationCount = maxCount;
	
}


void InsertCommonToTwoReadAligningHead(CommonKmerHead * commonKmerHead, KmerReadNodeHead * kmerReadNodeHead, KmerHashTableHead * kmerHashTableHead, long int hashIndex, unsigned long int readIndex, unsigned long  long int position, bool orien){
	if(kmerHashTableHead->kmerHashNode[hashIndex].startPositionInArray == -1){
		return;
	}
	unsigned long int i = kmerHashTableHead->kmerHashNode[hashIndex].startPositionInArray;
	unsigned long int kmer = kmerHashTableHead->kmerHashNode[hashIndex].kmer - 1;


	for(; i < kmerReadNodeHead->realCount; i++){
		if(kmerReadNodeHead->kmerReadNode[i].kmer != kmer){
			break;
		}
		if(kmerReadNodeHead->kmerReadNode[i].readIndex <= readIndex && readIndex <= kmerReadNodeHead->endReadIndex){
			continue;
		}
		commonKmerHead->commonKmer[commonKmerHead->realCount].readIndex = kmerReadNodeHead->kmerReadNode[i].readIndex;
		
		commonKmerHead->commonKmer[commonKmerHead->realCount].rightPosition = kmerReadNodeHead->kmerReadNode[i].position;
		commonKmerHead->commonKmer[commonKmerHead->realCount].leftPosition = position;
		if(orien == kmerReadNodeHead->kmerReadNode[i].orientation){
			commonKmerHead->commonKmer[commonKmerHead->realCount].orientation = 0;
		}else{
			commonKmerHead->commonKmer[commonKmerHead->realCount].orientation = 1;
		}
		commonKmerHead->realCount++;
		if(commonKmerHead->realCount >= commonKmerHead->allocationCount){
			ReAllocateCommonKmer(commonKmerHead);
		}
	}
	
}

CommonKmerHead * GetCommonKmerHeadAllThreadNew(KmerHashTableHead * kmerHashTableHead, KmerReadNodeHead * kmerReadNodeHead, ReadSetHead * readSetHead, long int kmerLength, char * readFile, char * outputFile, unsigned long  long int step, long int totalThreadNumber, long int smallKmerLength, long int smallIntervalDistance, long int largeIntervalDistance, long int overlapLengthCutOff, float lengthRatio, long int subReadCount){
	
	long int startReadIndex = 0;
	long int endReadIndex = subReadCount;
	FILE * fp = NULL;
	if((fp = fopen(outputFile,"w")) == NULL){
        printf("%s, does not exist!", outputFile);
        exit(0);
   	}
	fclose(fp);
	
	while(true){
		if(startReadIndex >= readSetHead->readCount){
			break;
		}
		if(endReadIndex >= readSetHead->readCount){
			endReadIndex = readSetHead->readCount - 1;
		}
		
		InitKmerReadNodeHeadSub(readSetHead, kmerReadNodeHead, kmerHashTableHead, kmerLength, step, startReadIndex, endReadIndex);
		
		GetCommonKmerHeadAllThread(kmerHashTableHead, kmerReadNodeHead, readSetHead, kmerLength, readFile, outputFile, step, totalThreadNumber, smallKmerLength, smallIntervalDistance, largeIntervalDistance, overlapLengthCutOff, lengthRatio,startReadIndex);
		
		startReadIndex = endReadIndex + 1;
		endReadIndex = endReadIndex + subReadCount;

	}
	
}

CommonKmerHead * GetCommonKmerHeadAllThread(KmerHashTableHead * kmerHashTableHead, KmerReadNodeHead * kmerReadNodeHead, ReadSetHead * readSetHead, long int kmerLength, char * readFile, char * outputFile, unsigned long  long int step, long int totalThreadNumber, long int smallKmerLength, long int smallIntervalDistance, long int largeIntervalDistance, long int overlapLengthCutOff, float lengthRatio, long int startReadIndex){
	
	pthread_t tid[totalThreadNumber];
    
    long int i = 0;
    GetCommonKmerHeadP * getCommonKmerHeadP = new GetCommonKmerHeadP[totalThreadNumber];
    long int * threadSignal = new long int[totalThreadNumber]; 
    
	printf("The number of enabled threads is: %ld;\n",totalThreadNumber);
	printf("The length of kmer is: %ld;\n",kmerLength);
	printf("The length of smallkmer is: %ld;\n",smallKmerLength);
	printf("The minimum overlap length of each pair reads is: %ld;\n",overlapLengthCutOff);
	
		
    for(i = 0; i< totalThreadNumber; i++){
        getCommonKmerHeadP[i].kmerHashTableHead = kmerHashTableHead;
        getCommonKmerHeadP[i].kmerReadNodeHead = kmerReadNodeHead;
		getCommonKmerHeadP[i].readSetHead = readSetHead;
		getCommonKmerHeadP[i].kmerLength = kmerLength;
		getCommonKmerHeadP[i].readFile = readFile;
		getCommonKmerHeadP[i].smallIntervalDistance = smallIntervalDistance;
		getCommonKmerHeadP[i].largeIntervalDistance = largeIntervalDistance;
		getCommonKmerHeadP[i].overlapLengthCutOff = overlapLengthCutOff;
		getCommonKmerHeadP[i].smallKmerLength = smallKmerLength;
		getCommonKmerHeadP[i].lengthRatio = lengthRatio;
		getCommonKmerHeadP[i].startReadIndex = startReadIndex;
		
		char * outputFileTemp = (char *)malloc(sizeof(char)*(strlen(outputFile)+10));
		sprintf(outputFileTemp, "%s-%ld", outputFile, i);
		getCommonKmerHeadP[i].outputFile = outputFileTemp;
		
		getCommonKmerHeadP[i].step = step;
        getCommonKmerHeadP[i].threadIndex = i;
        getCommonKmerHeadP[i].totalThreadNumber = totalThreadNumber;
             
        if(pthread_create(&tid[i], NULL, GetCommonKmerHeadThread, (void *)&getCommonKmerHeadP[i])){
            cout<<"create thread wrong!"<<endl;
            exit(1);
       }      
    }
    
    for(i = 0; i < totalThreadNumber; i++){
        pthread_join(tid[i], NULL);
    }
	
	
	FILE * fp = NULL;

	long int Size = 10000;
	char * str = (char *)malloc(sizeof(char)*Size);
	for(i = 0; i < totalThreadNumber; i++){
        FILE * fpTemp = NULL;
		
		if((fpTemp = fopen(getCommonKmerHeadP[i].outputFile,"r")) == NULL){
        	printf("%s, does not exist!", getCommonKmerHeadP[i].outputFile);
        	exit(0);
    	}
		
		if((fp = fopen(outputFile,"a")) == NULL){
        	printf("%s, does not exist!", outputFile);
        	exit(0);
   		}
		setbuf(fp,NULL);
		
		while((fgets(str, Size, fpTemp)) != NULL){
			setbuf(fp,NULL);
			long int len = strlen(str);
			if(str[len - 1]=='\n' || str[len - 1]=='\r'){
				str[len - 1] = '\0';
			}
			fprintf(fp, "%s\n", str);
			fflush(fp);

		}
		fflush(fp);
		fclose(fp);
		fclose(fpTemp);
		remove(getCommonKmerHeadP[i].outputFile);
	}
	
}

void * GetCommonKmerHeadThread(void * arg){
	GetCommonKmerHeadP * getCommonKmerHeadP = (GetCommonKmerHeadP *)arg;
	 
	KmerHashTableHead * kmerHashTableHead = getCommonKmerHeadP->kmerHashTableHead;
    KmerReadNodeHead * kmerReadNodeHead = getCommonKmerHeadP->kmerReadNodeHead;
	ReadSetHead * readSetHead = getCommonKmerHeadP->readSetHead;
    long int kmerLength = getCommonKmerHeadP->kmerLength;
	char * readFile = getCommonKmerHeadP->readFile;
    char * outputFile = getCommonKmerHeadP->outputFile;
	long int step = getCommonKmerHeadP->step;
    long int threadIndex = getCommonKmerHeadP->threadIndex;
	long int totalThreadNumber = getCommonKmerHeadP->totalThreadNumber;
	long int smallIntervalDistance = getCommonKmerHeadP->smallIntervalDistance;
	long int largeIntervalDistance = getCommonKmerHeadP->largeIntervalDistance;
	long int overlapLengthCutOff = getCommonKmerHeadP->overlapLengthCutOff;
	long int startReadIndex = getCommonKmerHeadP->startReadIndex;

	AdjGraphHead * G = (AdjGraphHead *)malloc(sizeof(AdjGraphHead));
	G->allocationCountGraph = 20000;
	G->graph = (AdjGraph*)malloc(sizeof(AdjGraph)* G->allocationCountGraph);
	G->realCountGraph = 0;
	G->reverseAllocationCountGraph = 20000;
	G->reverseGraph = (AdjGraph*)malloc(sizeof(AdjGraph)* G->reverseAllocationCountGraph);
	G->reverseRealCountGraph = 0;
	G->allocationCountArc = 20000;
	G->arcIndex = (ArcIndex*)malloc(sizeof(ArcIndex)* G->allocationCountArc);
	G->realCountArc = 0;
	G->nodeCount = 0;
	G->kmerLength = getCommonKmerHeadP->kmerLength;
	G->smallKmerLength = getCommonKmerHeadP->smallKmerLength;
	G->lengthRatio = getCommonKmerHeadP->lengthRatio;
	G->overlapLengthCutOff = getCommonKmerHeadP->overlapLengthCutOff;
	
	G->largestIntervalDistance = smallIntervalDistance;
	
	for(long int n=0;n<G->allocationCountGraph;n++){
		G->graph[n].dataLeft = -1;
		G->graph[n].dataRight = -1;
		
		G->graph[n].visit = 0;
		G->reverseGraph[n].dataLeft = -1;
		G->reverseGraph[n].dataRight = -1;

		G->reverseGraph[n].visit = 0;
		G->arcIndex[n].startIndex = 0;
		G->arcIndex[n].endIndex = 0;
	}
	
	
	
	long int maxCount = 100000;
	
	CommonKmerHead * commonKmerHead = (CommonKmerHead *)malloc(sizeof(CommonKmerHead));
	commonKmerHead->commonKmer = (CommonKmer *)malloc(sizeof(CommonKmer)*maxCount);
	commonKmerHead->realCount = 0;
	commonKmerHead->allocationCount = maxCount;
	
	
	
	AdjGraphHead * localG = (AdjGraphHead *)malloc(sizeof(AdjGraphHead));
	localG->largestIntervalDistance = largeIntervalDistance;
	localG->allocationCountGraph = 2*largeIntervalDistance;
	localG->graph = (AdjGraph*)malloc(sizeof(AdjGraph)* localG->allocationCountGraph);
	localG->realCountGraph = 0;
	localG->allocationCountArc = 20000;
	localG->arcIndex = (ArcIndex*)malloc(sizeof(ArcIndex)* localG->allocationCountArc);
	localG->realCountArc = 0;
	localG->nodeCount = 0;
	localG->localLeftRead = (char *)malloc(sizeof(char)*(localG->allocationCountGraph + 2*kmerLength));
	localG->localRightRead = (char *)malloc(sizeof(char)*(localG->allocationCountGraph + 2*kmerLength));
	localG->lengthRatio = getCommonKmerHeadP->lengthRatio;
	localG->kmerLength = getCommonKmerHeadP->kmerLength;
	localG->smallKmerLength = getCommonKmerHeadP->smallKmerLength;
	CommonKmerHead * localCommonKmerHead = (CommonKmerHead *)malloc(sizeof(CommonKmerHead));
	localCommonKmerHead->allocationCount = 2*largeIntervalDistance;
	localCommonKmerHead->commonKmer = (CommonKmer *)malloc(sizeof(CommonKmer)*localCommonKmerHead->allocationCount*2);
	localCommonKmerHead->realCount = 0;
	
	
	
	unsigned long int kmerInteger = 0;
	long int hashIndex;
	unsigned long int readIndex = 0;
	long int Size = 200000;
	long int readLength = 0;

	char * kmer1 = (char *)malloc(sizeof(char)*(kmerLength + 1));
	char * kmer2 = (char *)malloc(sizeof(char)*(kmerLength + 1));
	FILE * fp = fopen(readFile,"r");
	FILE * fp1 = fopen(outputFile,"w");
	
	long int * forwardKmerCount = (long int *)malloc(sizeof(long int)*readSetHead->readCount);
	long int * reverseKmerCount = (long int *)malloc(sizeof(long int)*readSetHead->readCount);
	
	char buf[1024];
	time_t timep;
	double sencond;
	int pi = 0;
	float per = 0.0;
	
	pi = readSetHead->readCount/totalThreadNumber;
	
	for(long int i = startReadIndex; i < readSetHead->readCount; i++){
		
		commonKmerHead->readIndex = i + 1;
		readIndex = i + 1;
		if(readIndex % totalThreadNumber != threadIndex){
			continue;
		}
		
		if(readIndex % pi == 0 && totalThreadNumber > 1){
			per = (float(readIndex)/readSetHead->readCount);
			printf("Thread %ld starts to detect overlap!\n", threadIndex);
		}
		
		readLength = readSetHead->readSet[i].readLength;
		for(long int j = 0; j < readLength - kmerLength + 1;j = j + step){
			strncpy(kmer1,readSetHead->readSet[i].read + j, kmerLength);
			kmer1[kmerLength] = '\0';
			
			if(DetectSameKmer(kmer1, kmerLength) != true){
				continue;
			}
			
			SetBitKmer(&kmerInteger, kmerLength, kmer1);
			hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
			if(hashIndex!=-1){
				InsertCommonToTwoReadAligningHead(commonKmerHead, kmerReadNodeHead, kmerHashTableHead, hashIndex, readIndex, j, 1);
			}else{
				ReverseComplementKmer(kmer1, kmerLength);
				SetBitKmer(&kmerInteger, kmerLength, kmer1);
				hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
					InsertCommonToTwoReadAligningHead(commonKmerHead, kmerReadNodeHead, kmerHashTableHead, hashIndex, readIndex, j, 0);
			}	
		}
		
		RemoveLowNumberKmer(commonKmerHead, forwardKmerCount, reverseKmerCount, readSetHead->readCount);

		if(commonKmerHead->realCount > 100000){
			heapSort(commonKmerHead->commonKmer, 0, commonKmerHead->realCount - 1);
		}else{
			sort(commonKmerHead->commonKmer, 0, commonKmerHead->realCount - 1);
		}
		
		RemoveMultipleSameKmer(commonKmerHead);
		
		GetOverlapResult(G, commonKmerHead, readSetHead, localG, localCommonKmerHead, fp1);

		commonKmerHead->realCount = 0;
				
	}
	
	fclose(fp);
	fflush(fp1);
	fclose(fp1);
	
}


void SubRemoveMultipleSameKmer(CommonKmerHead * commonKmerHead, long int startIndex, long int endIndex){
	bool token = false;
	bool token1 = false;
	for(long int j = startIndex; j < endIndex; j++){
		token = false;
		for(long int m = j + 1; m <= endIndex; m++){
			if(commonKmerHead->commonKmer[m].leftPosition == -1){
				continue;
			}
			
			if(commonKmerHead->commonKmer[j].leftPosition == commonKmerHead->commonKmer[m].leftPosition){
				token = true;
				commonKmerHead->commonKmer[m].leftPosition = -1;
				continue;
			}
			
			if(abs(commonKmerHead->commonKmer[j].leftPosition - commonKmerHead->commonKmer[m].leftPosition) == abs(commonKmerHead->commonKmer[j].rightPosition - commonKmerHead->commonKmer[m].rightPosition) 
			  && abs(commonKmerHead->commonKmer[j].leftPosition - commonKmerHead->commonKmer[m].leftPosition) < 10){
				token1 = true;
				commonKmerHead->commonKmer[m].leftPosition = -1;
			}
		}
		if(token != false){
			commonKmerHead->commonKmer[j].leftPosition = -1;
		}
	}
}

void RemoveMultipleSameKmer(CommonKmerHead * commonKmerHead){
	long int startIndex = 0;
	long int readIndex = commonKmerHead->commonKmer[0].readIndex;
	long int endIndex = -1;
	for(long int i = 0; i < commonKmerHead->realCount; i++){
		
		if(commonKmerHead->commonKmer[i].readIndex != readIndex){
			endIndex = i - 1;
			
			SubRemoveMultipleSameKmer(commonKmerHead, startIndex, endIndex);
			
			startIndex = i;
			readIndex = commonKmerHead->commonKmer[i].readIndex;
		}
	}
	endIndex = commonKmerHead->realCount - 1;
	
	SubRemoveMultipleSameKmer(commonKmerHead, startIndex, endIndex);
	
	long int shiftCount = 0;
	
	for(long int i = 0; i < commonKmerHead->realCount; i++){
		if(commonKmerHead->commonKmer[i].leftPosition == -1){
			shiftCount++;
		}else if(shiftCount != 0){
			commonKmerHead->commonKmer[i - shiftCount].leftPosition = commonKmerHead->commonKmer[i].leftPosition;
			commonKmerHead->commonKmer[i - shiftCount].rightPosition = commonKmerHead->commonKmer[i].rightPosition;
			commonKmerHead->commonKmer[i - shiftCount].orientation = commonKmerHead->commonKmer[i].orientation;
			commonKmerHead->commonKmer[i - shiftCount].readIndex = commonKmerHead->commonKmer[i].readIndex;
		}
	}

	commonKmerHead->realCount = commonKmerHead->realCount - shiftCount;
}



void RemoveLowNumberKmer(CommonKmerHead * commonKmerHead, long int * forwardKmerCount, long int * reverseKmerCount, long int readCount){
	
	for(long int i = 0; i < readCount; i++){
		forwardKmerCount[i] = 0;
		reverseKmerCount[i] = 0;
		
		
	}
	for(long int i = 0; i < commonKmerHead->realCount; i++){
		if(commonKmerHead->commonKmer[i].orientation == 0){
			forwardKmerCount[commonKmerHead->commonKmer[i].readIndex - 1]++;
		}else{
			reverseKmerCount[commonKmerHead->commonKmer[i].readIndex - 1]++;
		}
	}

	for(long int i = 0; i < commonKmerHead->realCount; i++){
		if(forwardKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] < 15 && commonKmerHead->commonKmer[i].orientation == 0){
			commonKmerHead->commonKmer[i].leftPosition = -1;
		}
		if(forwardKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] >= 15 && commonKmerHead->commonKmer[i].orientation == 0 && forwardKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] <= reverseKmerCount[commonKmerHead->commonKmer[i].readIndex - 1]){
			commonKmerHead->commonKmer[i].leftPosition = -1;
		}
		
		if(reverseKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] < 15 && commonKmerHead->commonKmer[i].orientation == 1){
			commonKmerHead->commonKmer[i].leftPosition = -1;
		}
		
		if(reverseKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] >= 15 && commonKmerHead->commonKmer[i].orientation == 1 && forwardKmerCount[commonKmerHead->commonKmer[i].readIndex - 1] >= reverseKmerCount[commonKmerHead->commonKmer[i].readIndex - 1]){
			commonKmerHead->commonKmer[i].leftPosition = -1;
		}
		
	}
	
	long int shiftCount = 0;
	
	for(long int i = 0; i < commonKmerHead->realCount; i++){
		if(commonKmerHead->commonKmer[i].leftPosition == -1){
			shiftCount++;
		}else if(shiftCount != 0){
			commonKmerHead->commonKmer[i - shiftCount].leftPosition = commonKmerHead->commonKmer[i].leftPosition;
			commonKmerHead->commonKmer[i - shiftCount].rightPosition = commonKmerHead->commonKmer[i].rightPosition;
			commonKmerHead->commonKmer[i - shiftCount].orientation = commonKmerHead->commonKmer[i].orientation;
			commonKmerHead->commonKmer[i - shiftCount].readIndex = commonKmerHead->commonKmer[i].readIndex;
		}
	}
	
	commonKmerHead->realCount = commonKmerHead->realCount - shiftCount;
}


void sortGraph(AdjGraph * graph, long int left, long int right){
	if(left >= right){
        return ;
    }
    long int i = left;
    long int j = right;
    long int key = graph[left].dataLeft;
	long int dataRight = graph[left].dataRight;
	
    while(i < j){
        while(i < j && key <= graph[j].dataLeft){
            j--;
        }
		
		if(i < j){
			graph[i].dataLeft = graph[j].dataLeft;
			graph[i].dataRight = graph[j].dataRight;
			i++;
		}

         
        while(i < j && key > graph[i].dataLeft){
            i++;
        }
		
		if(i < j){
			graph[j].dataLeft = graph[i].dataLeft;
			graph[j].dataRight = graph[i].dataRight;
			j--;
		}
        
    }
    
    graph[i].dataLeft = key;
	graph[i].dataRight = dataRight;
	
    sortGraph(graph, left, i - 1);
    sortGraph(graph, i + 1, right);
}

void swapCommonKmer(CommonKmer *a, long int left, long int right){
	long int temp = a[left].readIndex;
 	a[left].readIndex = a[right].readIndex;
	a[right].readIndex = temp;
	
	temp = a[left].leftPosition;
 	a[left].leftPosition = a[right].leftPosition;
	a[right].leftPosition = temp;
	
	temp = a[left].rightPosition;
 	a[left].rightPosition = a[right].rightPosition;
	a[right].rightPosition = temp;
	
	bool temp1 = a[left].orientation;
 	a[left].orientation = a[right].orientation;
	a[right].orientation = temp1;
}

void downToMaxHeap(CommonKmer *a, long int bgn, long int end){
    long int child;
    long int parent = bgn;
    while ((child = parent * 2 + 1) < end)
    {
        if ((child < end - 1) && ((a[child].readIndex < a[child + 1].readIndex) || (a[child].readIndex == a[child + 1].readIndex && a[child].leftPosition < a[child + 1].leftPosition)))
            ++child;   
        if ((a[child].readIndex > a[parent].readIndex) || (a[child].readIndex == a[parent].readIndex && a[child].leftPosition > a[parent].leftPosition))
            swapCommonKmer(a, child, parent);
        else
            break;
        parent = child;
    }
}

void buildMaxHeap(CommonKmer *a, long int bgn, long int end){
    if (bgn >= end - 1)
        return;

    int parent = end / 2 - 1;
    while (parent >= 0)
    {
        downToMaxHeap(a, parent, end);
        --parent;
    }
}

void heapSort(CommonKmer *a, long int bgn, long int end){
    buildMaxHeap(a, bgn, end);

    while (end > 1)
    {
        swapCommonKmer(a, 0, --end);
        downToMaxHeap(a, 0, end);
    }
}


void sort(CommonKmer *a, long int left, long int right)
{

	if(left >= right){
        return ;
    }

    long int i = left;
    long int j = right;
    long int key = a[left].readIndex;
	long int leftPosition = a[left].leftPosition;
	long int rightPosition = a[left].rightPosition;
    bool orientation = a[left].orientation;
	
    while(i < j){
        while(i < j && (key < a[j].readIndex || (key == a[j].readIndex && leftPosition <= a[j].leftPosition))){
            j--;
        }
		
		if(i < j){
			a[i].readIndex = a[j].readIndex;
			a[i].leftPosition = a[j].leftPosition;
			a[i].rightPosition = a[j].rightPosition;
			a[i].orientation = a[j].orientation;
			i++;
		}
      
        while(i < j && (key > a[i].readIndex || (key == a[i].readIndex && leftPosition > a[i].leftPosition))){
            i++;
        }
		
		if(i < j){
			a[j].readIndex = a[i].readIndex;
			a[j].leftPosition = a[i].leftPosition;
			a[j].rightPosition = a[i].rightPosition;
			a[j].orientation = a[i].orientation;
			j--;
		}
  
    }

    a[i].readIndex = key;
	a[i].leftPosition = leftPosition;
	a[i].rightPosition = rightPosition;
	a[i].orientation = orientation;

    sort(a, left, i - 1);
    sort(a, i + 1, right); 
}


void DestroyGraph(AdjGraphHead * G){
	G->realCountGraph = 0;
	G->realCountArc = 0;
}


long int Overlap_DisplayLocalRegion(AdjGraphHead * G,long int leftLen,long int rightLen){
	if(G->realCountArc <= 0){
		return 0;
	}
	long int leftStartpos=-1,leftEndpos=-1,rightStartpos=-1,rightEndpos=-1;
	
	leftStartpos = G->graph[G->arcIndex[0].startIndex].dataLeft;
	leftEndpos = G->graph[G->arcIndex[G->realCountArc - 1].endIndex].dataLeft;
	rightStartpos = G->graph[G->arcIndex[0].startIndex].dataRight;
	rightEndpos = G->graph[G->arcIndex[G->realCountArc - 1].endIndex].dataRight;

	long int distance = abs(leftLen - rightLen);
	
	distance = 300;

	if((leftStartpos > distance || rightStartpos > distance) || (leftLen - leftEndpos > distance || rightLen - rightEndpos > distance)){
		return 0;
	}else{
		return 1;
	}
	
	

}


long int Overlap_Display_Graph(AdjGraphHead * G, long int leftIndex,long int rightIndex,bool orien,long int leftLen,long int rightLen,FILE * fp, long int leftStartpos, long int leftEndpos, long int rightStartpos, long int rightEndpos, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, char * leftRead, char * rightRead){
	if(leftIndex == rightIndex){
		return 0;
	}
	long int kmerLength = G->kmerLength;
	
	long int i;
	long int t = 0;

	long int overlenleft,overlenright;
	long int MinOverLen,MaxOverLen;

	long int tempLeftStart = leftStartpos;
	long int tempLeftEnd = leftEndpos;
	long int tempRightStart = rightStartpos;
	long int tempRightEnd = rightEndpos;
	
	MinOverLen = min(leftEndpos -leftStartpos,rightEndpos - rightStartpos);
	MaxOverLen = max(leftEndpos -leftStartpos,rightEndpos - rightStartpos);

	if(MinOverLen < 300){return 0;} 
	if((float)(MaxOverLen - MinOverLen)/MaxOverLen > G->lengthRatio ){return 0;} 

	t=G->largestIntervalDistance;
	long int maxIntervalDistance = localG->largestIntervalDistance;
	
	long int localLeftStart = 0;
	long int localLeftEnd = 0;
	long int localRightStart = 0;
	long int localRightEnd = 0;

	long int dd = 0;
	
	long int alignLen = 30;
	float alignRtio = 0.6;

	if(orien==0){
		if((rightStartpos > leftStartpos && leftStartpos > maxIntervalDistance) || (rightStartpos <= leftStartpos && rightStartpos > maxIntervalDistance)){ 
			return 0;
		}
		if((leftLen-leftEndpos > rightLen-rightEndpos && rightLen-rightEndpos > maxIntervalDistance) || (leftLen-leftEndpos <= rightLen-rightEndpos && leftLen-leftEndpos > maxIntervalDistance)){ 
			return 0;
		}
		
		if(rightStartpos > leftStartpos){
			if(leftStartpos > maxIntervalDistance){
				return 0;
			}else if(leftStartpos > t){
				localLeftStart = 0;
				localLeftEnd = leftStartpos - 1;
				localRightStart = rightStartpos - leftStartpos;
				localRightEnd = rightStartpos - 1;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 0);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			rightStartpos = rightStartpos - leftStartpos;
			leftStartpos = 1;
		}else{
			if(rightStartpos > maxIntervalDistance){
				return 0;
			}else if(rightStartpos > t){
				localLeftStart = leftStartpos - rightStartpos;
				localLeftEnd = leftStartpos - 1;
				localRightStart = 0;
				localRightEnd = rightStartpos - 1;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 0);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			leftStartpos = leftStartpos-rightStartpos;
			rightStartpos = 1;
		}
		
		if(leftLen-leftEndpos > rightLen-rightEndpos){
			if(rightLen-rightEndpos > maxIntervalDistance){
				return 0;
			}else if(rightLen-rightEndpos > t){
				localLeftStart = leftEndpos + kmerLength;
				localLeftEnd = leftEndpos + rightLen - rightEndpos - 1;
				localRightStart = rightEndpos + kmerLength;
				localRightEnd = rightLen - 1;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 0);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			leftEndpos = rightLen - rightEndpos + leftEndpos;
			rightEndpos = rightLen;
			
		}else{
			if(leftLen-leftEndpos > maxIntervalDistance){
				return 0;
			}else if(leftLen-leftEndpos > t){
				localLeftStart = leftEndpos + kmerLength;
				localLeftEnd = leftLen - 1;
				localRightStart = rightEndpos + kmerLength;
				localRightEnd = rightEndpos + leftLen - leftEndpos - 1;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 0);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			rightEndpos = leftLen - leftEndpos + rightEndpos;
			leftEndpos = leftLen;
		}
		
	}

	if(orien==1){
		
		if((rightLen-rightEndpos > leftStartpos && leftStartpos > maxIntervalDistance) || (rightLen-rightEndpos <= leftStartpos && rightLen-rightEndpos > maxIntervalDistance)){
			return 0;
		}
		if((leftLen-leftEndpos > rightStartpos && rightStartpos > maxIntervalDistance) || (leftLen-leftEndpos <= rightStartpos && leftLen-leftEndpos > maxIntervalDistance)){
			return 0;
		}
		
		if(rightLen-rightEndpos > leftStartpos){
			if(leftStartpos > maxIntervalDistance){
				return 0;
			}else if(leftStartpos > t){

				localLeftStart = 0;
				localLeftEnd = leftStartpos - 1;
				localRightStart = rightEndpos + kmerLength;
				localRightEnd = rightEndpos + kmerLength + leftStartpos;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 1);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			rightEndpos = rightEndpos + leftStartpos;
			leftStartpos = 1;
		}else{
			if(rightLen-rightEndpos > maxIntervalDistance){
				return 0;
			}else if(rightLen-rightEndpos > t){

				localLeftStart = leftStartpos - rightLen + rightEndpos + kmerLength;
				localLeftEnd = leftStartpos - 1;
				localRightStart = rightEndpos + kmerLength;
				localRightEnd = rightLen - 1;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 1);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			leftStartpos = leftStartpos - rightLen + rightEndpos;
			rightEndpos = rightLen;
		}
		
		if(leftLen-leftEndpos > rightStartpos){
			if(rightStartpos > maxIntervalDistance){ 
				return 0;
			}else if(rightStartpos > t){

				localLeftStart = leftEndpos + kmerLength;
				localLeftEnd = leftEndpos + kmerLength + rightStartpos - 1;
				localRightStart = 0;
				localRightEnd = rightStartpos - 1;

				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 1);
				
				if(ss < alignRtio){
					return 0;
				}
				
			}
			leftEndpos = leftEndpos + rightStartpos;
			rightStartpos = 1;
		}else{
			if(leftLen-leftEndpos > maxIntervalDistance){
				return 0;
			}else if(leftLen-leftEndpos > t){

				localLeftStart = leftEndpos + kmerLength;
				localLeftEnd = leftLen - 1;
				localRightStart = rightStartpos - leftLen + leftEndpos + kmerLength;
				localRightEnd = rightStartpos - 1;
				
				long int ss = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, localLeftStart, localLeftEnd, localRightStart, localRightEnd, G->smallKmerLength, 1);

				if(ss < alignRtio){
					return 0;
				}
				
			}
			
			rightStartpos = rightStartpos - leftLen + leftEndpos;
			leftEndpos = leftLen;
			
		}

	}

	overlenleft = leftEndpos - leftStartpos;
	overlenright = rightEndpos - rightStartpos;
	
	long int tempMaxLen = MaxOverLen;
	long int tempMinLen = MinOverLen;
	
	MinOverLen = min(overlenleft,overlenright);
	MaxOverLen = max(overlenleft,overlenright);

	if(MaxOverLen<G->overlapLengthCutOff){ return 0;}
	if((float)(MaxOverLen - MinOverLen)/MaxOverLen > G->lengthRatio ){ return 0;} 
	
	if(leftIndex > rightIndex){	
		fprintf(fp,"%ld,%ld,%d,%ld,%ld,%ld,%ld,%ld,%ld\n",leftIndex,rightIndex,orien,
												leftStartpos,leftEndpos,leftLen,rightStartpos,rightEndpos,rightLen);
	}else{
		fprintf(fp,"%ld,%ld,%d,%ld,%ld,%ld,%ld,%ld,%ld\n",rightIndex,leftIndex,orien,
												rightStartpos,rightEndpos,rightLen,leftStartpos,leftEndpos,leftLen);
	}

	return 1;

}

void ReAllocateAdjGraph(AdjGraphHead * G, long int a){
	
	if(a == 0){
		G->allocationCountGraph = G->allocationCountGraph*1.5;
		AdjGraph * graph = (AdjGraph *)malloc(sizeof(AdjGraph)*G->allocationCountGraph);
		for(long int i = 0; i < G->realCountGraph; i++){
			graph[i].dataLeft = G->graph[i].dataLeft;
			graph[i].dataRight = G->graph[i].dataRight;

			graph[i].visit = G->graph[i].visit;
		}
		free(G->graph);
		G->graph = graph;
	}else{
		G->reverseAllocationCountGraph = G->reverseAllocationCountGraph*1.5;
		AdjGraph * graph = (AdjGraph *)malloc(sizeof(AdjGraph)*G->reverseAllocationCountGraph);
		for(long int i = 0; i < G->reverseRealCountGraph; i++){
			graph[i].dataLeft = G->reverseGraph[i].dataLeft;
			graph[i].dataRight = G->reverseGraph[i].dataRight;

			graph[i].visit = G->reverseGraph[i].visit;
		}
		free(G->reverseGraph);
		G->reverseGraph = graph;
	}
		
}


void ReAllocateArcIndex(AdjGraphHead * G){
	
	G->allocationCountArc = G->allocationCountArc*2;
	
	ArcIndex * arcIndex = (ArcIndex *)malloc(sizeof(ArcIndex)*G->allocationCountArc);
	
	for(long int i = 0; i < G->realCountArc; i++){
		arcIndex[i].startIndex = G->arcIndex[i].startIndex;
		arcIndex[i].endIndex = G->arcIndex[i].endIndex;
	}
	
	free(G->arcIndex);
	
	G->arcIndex = arcIndex;
}


long int AddEdgeInGraph(AdjGraphHead * G, bool orientation, long int leftIndex , long int rightIndex, long int largestIntervalDistance, long int maxIntervalDistance, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, char * leftRead, char * rightRead){
	AdjGraph * graph = NULL;
	if(orientation == 0){
		if(G->graph[rightIndex].visit == true){
			return 0;
		}
		graph = G->graph;
	}else{
		if(G->reverseGraph[rightIndex].visit == true){
			return 0;
		}
		graph = G->reverseGraph;
	}
	
	long int kmerLength = 13;
	
	long int n1 = graph[leftIndex].dataLeft;
	long int n2 = graph[leftIndex].dataRight;
	long int m1 = graph[rightIndex].dataLeft;
	long int m2 = graph[rightIndex].dataRight;
			
		
	if(abs(m1-n1) > maxIntervalDistance || abs(m2-n2) > maxIntervalDistance){
		return 0;
	}
	
	long int minvalue = min(abs(m1-n1), abs(m2-n2));
	long int maxvalue = max(abs(m1-n1), abs(m2-n2));
	
	long int temp = 0;
	float ss1 = 0;
	long int alignLength = 30;
	
	
	if(orientation == 0 && float(maxvalue - minvalue)/maxvalue < G->lengthRatio && ((n1 < m1 && n2 < m2) || (n1>m1 && n2>m2))){
		
		if(abs(m1 - n1) <= largestIntervalDistance && abs(m2 - n2) <= largestIntervalDistance){
			temp = 1;
		}else if(abs(m1 - n1) <= maxIntervalDistance && abs(m2 - n2) <= maxIntervalDistance){
			if(float(maxvalue - minvalue)/maxvalue < 0.1){
				temp = 1;
			}else if(n1 < m1 && n2 < m2){
				temp = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, n1, m1, n2, m2, G->smallKmerLength, 0);
			}else if(n1>m1 && n2>m2){
				temp = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, m1, n1, m2, n2, G->smallKmerLength, 0);
			}
		}
		if(temp == 1 || ss1 > 0.5){
			G->arcIndex[G->realCountArc].startIndex = leftIndex;
			G->arcIndex[G->realCountArc].endIndex = rightIndex;
			G->arcIndex[G->realCountArc].weight = ((((float)minvalue/maxvalue) + (1 - (float)maxvalue/maxIntervalDistance))/2) * (maxvalue);
			G->realCountArc++;
			if(G->realCountArc >= G->allocationCountArc){
				ReAllocateArcIndex(G);
			}
			return 1;
		}else{
			return 2;
		}
		
		
	}
	
	if(orientation == 1 && float(maxvalue - minvalue)/maxvalue < G->lengthRatio && ((n1<m1 && n2>m2) || (n1>m1 && n2<m2))){
		
		if(abs(m1 - n1) <= largestIntervalDistance && abs(m2 - n2) <= largestIntervalDistance){
			temp = 1;
		}else if(abs(m1 - n1) <= maxIntervalDistance && abs(m2 - n2) <= maxIntervalDistance){
			if(float(maxvalue - minvalue)/maxvalue < 0.1){
				temp = 1;
			}else if(n1<m1 && n2>m2){
				temp = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, n1, m1, m2, n2, G->smallKmerLength, 1);
			}else if(n1>m1 && n2<m2){
				temp = GetCommonShorterKmer(localG, localCommonKmerHead, leftRead, rightRead, m1, n1, n2, m2, G->smallKmerLength, 1);
			}
		}
		if(temp == 1 || ss1 > 0.5){
			G->arcIndex[G->realCountArc].startIndex = leftIndex;
			G->arcIndex[G->realCountArc].endIndex = rightIndex;
			G->arcIndex[G->realCountArc].weight = ((((float)minvalue/maxvalue) + (1 - (float)maxvalue/maxIntervalDistance))/2) * (maxvalue);
			G->realCountArc++;
			if(G->realCountArc >= G->allocationCountArc){
				ReAllocateArcIndex(G);
			}
			return 1;
		}else{
			return 2;
		}
		
	}
	
	return 0;
}


long int CreatGraphSinglePath(AdjGraphHead * G, CommonKmer * commonKmer, unsigned long int startIndex, unsigned long int endIndex, long int a, long int leftIndex, long int rightIndex, long int leftLen, long int rightLen, FILE * fp, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, char * leftRead, char * rightRead) {
	
	AdjGraph * graph = NULL;
	long int realCountGraph = 0;
	if(a == 0){
		graph = G->graph;
		realCountGraph = G->realCountGraph;
	}else{
		graph = G->reverseGraph;
		realCountGraph = G->reverseRealCountGraph;
	}
	
	long int largestIntervalDistance = G->largestIntervalDistance;
	long int maxIntervalDistance = localG->largestIntervalDistance;
	
	for(long int n = 0; n < realCountGraph; n++){
		graph[n].visit = 0;
	}

	startIndex = 0;
	endIndex = 0;
	long int rightReadLength = strlen(rightRead);
	
	for(long int j = 0; j < realCountGraph - 1; j++){
		if(graph[j].visit == 1){
			continue;
		}
		
		if(graph[j].visit == 0 && a == 0){
			if(graph[j].dataLeft > graph[j].dataRight){
				if(graph[j].dataRight > maxIntervalDistance){
					continue;
				}
			}else{
				if(graph[j].dataLeft > maxIntervalDistance){
					continue;
				}
			}
		}
		
		if(graph[j].visit == 0 && a == 1){
			if(graph[j].dataLeft > rightReadLength - graph[j].dataRight){
				if(rightReadLength - graph[j].dataRight > maxIntervalDistance){
					continue;
				}
			}else{
				if(graph[j].dataLeft > maxIntervalDistance){
					continue;
				}
			}
		}
		startIndex = j;
		if(realCountGraph - startIndex <= 2){
			return 0;
		}
		endIndex = startIndex + 1;
		long int lastEndIndex = endIndex;
		long int iniStartIndex = startIndex;
		bool token = false;
		
		int edgeCount = 0;
		
		while(endIndex < realCountGraph && startIndex < realCountGraph){
			
			long int edge = AddEdgeInGraph(G, a, startIndex, endIndex, largestIntervalDistance, maxIntervalDistance, localG, localCommonKmerHead, leftRead, rightRead);
			
			if(edge == 1){
				graph[startIndex].visit = 1;
				graph[endIndex].visit = 1;
				lastEndIndex = endIndex;
				startIndex = endIndex;
				endIndex = startIndex + 1;
				token = true;
				edgeCount++;
			}else if(edge == 0){
				endIndex++;
				continue;
			}else{
				break;
			}
		}
		if(token == true && edgeCount > 2){
			long int result = 0;
			
			if(a == 0){
				result = Overlap_Display_Graph(G, leftIndex,rightIndex,a,leftLen,rightLen,fp,graph[iniStartIndex].dataLeft, graph[lastEndIndex].dataLeft, graph[iniStartIndex].dataRight, graph[lastEndIndex].dataRight, localG,localCommonKmerHead, leftRead, rightRead);
			}else{
				result = Overlap_Display_Graph(G, leftIndex,rightIndex,a,leftLen,rightLen,fp,graph[iniStartIndex].dataLeft, graph[lastEndIndex].dataLeft, graph[lastEndIndex].dataRight, graph[iniStartIndex].dataRight, localG,localCommonKmerHead, leftRead, rightRead);
			}
			
			if(result == 1){
				return 1;
			}
		}
	}
	
	return 0;	
	
}


long int CreatGraphLocalRegion(AdjGraphHead * G, long int distance) {
	long int count,throld;
	unsigned long  long int i,j,k,m,n;
	long int n1,n2,m1,m2;
	unsigned long int f1,f2;
	count = 0;
	long int minvalue;
	long int maxvalue;
	AdjGraph * graph = G->graph;
	long int realCountGraph = G->realCountGraph;

	long int firstIntervalDistance = distance;
	long int largestIntervalDistance = 2*distance;
	long int edgeCount = 0;
	
	for(j=0; j<realCountGraph; ++j){
		edgeCount = 0;
		for(k=j+1; k<realCountGraph && k < j + 5; ++k){
			
			n1 = graph[j].dataLeft;
			n2 = graph[j].dataRight;
			m1 = graph[k].dataLeft;
			m2 = graph[k].dataRight;

			if(abs(m1-n1) > largestIntervalDistance && abs(m2-n2) > largestIntervalDistance){
				break;
			}
			minvalue = min(abs(m1-n1),abs(m2-n2));
			maxvalue = max(abs(m1-n1),abs(m2-n2));
			if(float(maxvalue - minvalue)/maxvalue >= G->lengthRatio){
				continue;
			}
			if(n1<m1 && n2<m2){
				if(abs(m1-n1) < firstIntervalDistance && abs(m2-n2) < firstIntervalDistance || float(maxvalue - minvalue)/maxvalue < G->lengthRatio){	
					G->arcIndex[G->realCountArc].startIndex = j;
					G->arcIndex[G->realCountArc].endIndex = k;
					G->arcIndex[G->realCountArc].weight =1;
					G->realCountArc++;
					if(G->realCountArc >= G->allocationCountArc){
						ReAllocateArcIndex(G);
					}
				}
			}else if(n1>m1 && n2>m2){
				if(abs(m1-n1) < firstIntervalDistance && abs(m2-n2) < firstIntervalDistance || float(maxvalue - minvalue)/maxvalue < G->lengthRatio){
					G->arcIndex[G->realCountArc].startIndex = j;
					G->arcIndex[G->realCountArc].endIndex = k;
					G->arcIndex[G->realCountArc].weight =1;
					G->realCountArc++;  
					if(G->realCountArc >= G->allocationCountArc){
						ReAllocateArcIndex(G);
					}
				}
			}else{continue;}
				
		}
	}
	
	
}


void GetOverlapResult(AdjGraphHead * G, CommonKmerHead * commonKmerHead, ReadSetHead * readSetHead, AdjGraphHead * localG, CommonKmerHead * localCommonKmerHead, FILE * fp){
	long int leftIndex, rightIndex, orien;
	long int leftLen,rightLen;
	long int realCommonKmerCount;
	long int s = 0 ;
	unsigned long  long int readIndex = commonKmerHead->commonKmer[0].readIndex;
	unsigned long  long int startIndex = 0;
	unsigned long  long int endIndex = 0;
	bool orientation;
	
	for(long int m = 0; m < commonKmerHead->realCount; m++){
		if(commonKmerHead->commonKmer[m].readIndex != readIndex){
			endIndex = m - 1;
			
			leftIndex = commonKmerHead->readIndex;
			rightIndex = readIndex;
			
			leftLen = readSetHead->readSet[leftIndex - 1].readLength;
			rightLen = readSetHead->readSet[rightIndex - 1].readLength;
			
			if(G->realCountGraph < 5 && G->reverseRealCountGraph < 5){
			
				G->realCountGraph = 0;
				G->realCountArc = 0;
				G->reverseRealCountGraph = 0;
				
				if(commonKmerHead->commonKmer[m].orientation == 0){
					G->graph[G->realCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition;
					G->graph[G->realCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition;
					G->graph[G->realCountGraph].visit = 1;
					G->realCountGraph++;
					if(G->realCountGraph >= G->allocationCountGraph){
						ReAllocateAdjGraph(G, 0);
					}
				}else{
					G->reverseGraph[G->reverseRealCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition;
					G->reverseGraph[G->reverseRealCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition;
					G->reverseGraph[G->reverseRealCountGraph].visit = 1;
					G->reverseRealCountGraph++;
					if(G->reverseRealCountGraph >= G->reverseAllocationCountGraph){
						ReAllocateAdjGraph(G, 1);
					}
				}
				
				startIndex = m;
				readIndex = commonKmerHead->commonKmer[m].readIndex;
				continue;
			}
			if(G->realCountGraph > G->reverseRealCountGraph){
				orientation = false;
			}else{
				orientation = true;
			}

			char * leftRead = readSetHead->readSet[leftIndex - 1].read;
			char * rightRead = readSetHead->readSet[rightIndex - 1].read;
			
			long int result = CreatGraphSinglePath(G,commonKmerHead->commonKmer, startIndex, endIndex, orientation,leftIndex, rightIndex, leftLen, rightLen, fp, localG,localCommonKmerHead, leftRead, rightRead);
			
			G->realCountGraph = 0;
			G->realCountArc = 0;
			G->reverseRealCountGraph = 0;
			
			startIndex = m;
			readIndex = commonKmerHead->commonKmer[m].readIndex;
		}
			
		if(commonKmerHead->commonKmer[m].orientation == 0){
			G->graph[G->realCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition;
			G->graph[G->realCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition;

			G->graph[G->realCountGraph].visit = 1;
			G->realCountGraph++;
			if(G->realCountGraph >= G->allocationCountGraph){
				ReAllocateAdjGraph(G, 0);
			}
		}else{
			G->reverseGraph[G->reverseRealCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition;
			G->reverseGraph[G->reverseRealCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition;

			G->reverseGraph[G->reverseRealCountGraph].visit = 1;
			G->reverseRealCountGraph++;
			if(G->reverseRealCountGraph >= G->reverseAllocationCountGraph){
				ReAllocateAdjGraph(G, 1);
			}
		}
		
	}
	
	

	endIndex = commonKmerHead->realCount - 1;
	if(G->realCountGraph >= 5 || G->reverseRealCountGraph >= 5){
		if(G->realCountGraph > G->reverseRealCountGraph){
			orientation = false;
		}else{
			orientation = true;
		}
		leftIndex = commonKmerHead->readIndex;
		rightIndex = readIndex;
		leftLen = readSetHead->readSet[leftIndex - 1].readLength;
		rightLen = readSetHead->readSet[rightIndex - 1].readLength;
		char * leftRead = readSetHead->readSet[leftIndex - 1].read;
		char * rightRead = readSetHead->readSet[rightIndex - 1].read;	
		
		long int result = CreatGraphSinglePath(G,commonKmerHead->commonKmer, startIndex, endIndex, orientation, leftIndex, rightIndex, leftLen, rightLen, fp, localG,localCommonKmerHead, leftRead, rightRead);
	
		G->realCountGraph = 0;
		G->realCountArc = 0;
		G->reverseRealCountGraph = 0;
		
	}
	
}


void DetectCommon(CommonKmerHead * commonKmerHead, long int position, char * kmer, char * read, long int readLength, long int kmerLength, long int distance){
	bool t = true;
	distance = max(distance, 100) + 100;
	long int i = max(0, position - distance);
	for(; i < readLength - kmerLength + 1 && i < position + distance; i++){
		t = true;
		for(long int j = 0; j < kmerLength; j++){
			if(kmer[j] != read[i + j]){
				t = false;
				break;
			}
		}
		if(t == true){
			commonKmerHead->commonKmer[commonKmerHead->realCount].rightPosition = i;
			commonKmerHead->commonKmer[commonKmerHead->realCount].leftPosition = position;		
			commonKmerHead->realCount++;
			if(commonKmerHead->realCount >= commonKmerHead->allocationCount){
				ReAllocateCommonKmer(commonKmerHead);
			}
		}
	}
}

long int GetCommonShorterKmer(AdjGraphHead * G, CommonKmerHead * commonKmerHead, char * leftRead, char * rightRead, long int leftStartPosition, long int leftEndPosition, long int rightStartPosition, long int rightEndPosition, long int kmerLength, bool orientation)
{
	
	
	leftStartPosition = leftStartPosition + kmerLength;
	rightStartPosition = rightStartPosition + kmerLength;
	commonKmerHead->realCount = 1;
	G->realCountGraph = 0;
	G->realCountArc = 0;
	commonKmerHead->commonKmer[0].rightPosition = 0;
	commonKmerHead->commonKmer[0].leftPosition = 0;	
	
	long int tempLeftReadLength = leftEndPosition - leftStartPosition + 1;
	long int tempRightReadLength = rightEndPosition - rightStartPosition + 1;
	
	strncpy(G->localLeftRead, leftRead + leftStartPosition, tempLeftReadLength);

	strncpy(G->localRightRead, rightRead + rightStartPosition, tempRightReadLength);
	
	G->localLeftRead[tempLeftReadLength] = '\0';
	G->localRightRead[tempRightReadLength] = '\0';

	char * tempLeftRead = G->localLeftRead;
	char * tempRightRead = G->localRightRead;
	
	
	
	if(orientation == 1){		
		ReverseComplementKmer(tempRightRead, tempRightReadLength);
	}
	
	char kmer[kmerLength + 1];
	long int kmerCount = tempLeftReadLength - kmerLength + 1;
	
	for(long int i = 0; i < kmerCount; i++){
		strncpy(kmer, tempLeftRead + i, kmerLength);
		kmer[kmerLength] = '\0';
		if(DetectSameKmer(kmer,kmerLength) != true){
			continue;
		}
		DetectCommon(commonKmerHead, i, kmer, tempRightRead, tempRightReadLength, kmerLength, abs(tempLeftReadLength - tempRightReadLength));
	}
	
	commonKmerHead->commonKmer[commonKmerHead->realCount].rightPosition = tempRightReadLength - kmerLength;
	commonKmerHead->commonKmer[commonKmerHead->realCount].leftPosition = tempLeftReadLength - kmerLength;	
	commonKmerHead->realCount++;
	
	RemoveMultipleSameKmer(commonKmerHead);
	
	long int num = 2 + tempLeftReadLength/300;
	if(commonKmerHead->realCount < num){
		return 0;
	}
	
	for(long int i = 0; i < commonKmerHead->realCount; i++){
		G->graph[G->realCountGraph].dataLeft = commonKmerHead->commonKmer[i].leftPosition;
		G->graph[G->realCountGraph].dataRight = commonKmerHead->commonKmer[i].rightPosition;
		G->realCountGraph++;
		if(G->realCountGraph >= G->allocationCountGraph){
			ReAllocateAdjGraph(G, 0);
		}
	}
	
	CreatGraphLocalRegion(G,300);
	
	long int result = Overlap_DisplayLocalRegion(G,tempLeftReadLength,tempRightReadLength);
	
	return result;

}



#endif
