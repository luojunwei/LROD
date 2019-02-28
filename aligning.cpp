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

#include "aligning.h"
#include "read.h"

using namespace std;

unsigned long int yyCount0 = 0;
unsigned long int all = 0;
unsigned long int simCount = 0;

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


void InsertCommonToTwoReadAligningHead(CommonKmerHead * commonKmerHead, KmerReadNodeHead * kmerReadNodeHead, KmerHashTableHead * kmerHashTableHead, long int hashIndex, unsigned long int readIndex){
	unsigned long int i = kmerHashTableHead->kmerHashNode[hashIndex].startPositionInArray;
	unsigned long  int kmer = kmerHashTableHead->kmerHashNode[hashIndex].kmer - 1;
	unsigned long  int position = 0;
	bool orien = false;
	for(; i < kmerReadNodeHead->allocationCount; i++){
		if(kmerReadNodeHead->kmerReadNode[i].readIndex == readIndex){
			position = kmerReadNodeHead->kmerReadNode[i].position;
			orien = kmerReadNodeHead->kmerReadNode[i].orientation;
			break;
		}
	}
	
	i = kmerHashTableHead->kmerHashNode[hashIndex].startPositionInArray;
	
	for(; i < kmerReadNodeHead->allocationCount; i++){
		if(kmerReadNodeHead->kmerReadNode[i].kmer != kmer){
			break;
		}
		if(kmerReadNodeHead->kmerReadNode[i].readIndex <= readIndex){
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

CommonKmerHead * GetCommonKmerHeadAllThread(KmerHashTableHead * kmerHashTableHead, KmerReadNodeHead * kmerReadNodeHead, ReadSetHead * readSetHead, int kmerLength, char * readFile, char * outputFile, unsigned long  int step, int totalThreadNumber){
	
	pthread_t tid[totalThreadNumber];
    
    int i = 0;
    GetCommonKmerHeadP * getCommonKmerHeadP = new GetCommonKmerHeadP[totalThreadNumber];
    int * threadSignal = new int[totalThreadNumber]; 
    
    for(i = 0; i< totalThreadNumber; i++){
        getCommonKmerHeadP[i].kmerHashTableHead = kmerHashTableHead;
        getCommonKmerHeadP[i].kmerReadNodeHead = kmerReadNodeHead;
		getCommonKmerHeadP[i].readSetHead = readSetHead;
		getCommonKmerHeadP[i].kmerLength = kmerLength;
		getCommonKmerHeadP[i].readFile = readFile;
		
		char * outputFileTemp = (char *)malloc(sizeof(char)*(strlen(outputFile)+10));
		sprintf(outputFileTemp, "%s-%d", outputFile, i);
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
	
	
	FILE * fp = fopen(outputFile,"w");
	int Size = 10000;
	char * str = (char *)malloc(sizeof(char)*Size);
	for(i = 0; i < totalThreadNumber; i++){
        FILE * fpTemp = fopen(getCommonKmerHeadP[i].outputFile,"r");
		while(!feof(fpTemp)){
			fgets(str, Size, fpTemp);
			fprintf(fp, "%s", str);
		}
		fclose(fpTemp);
		remove(getCommonKmerHeadP[i].outputFile);
    }
	fflush(fp);
	fclose(fp);
	
	
	
}

void * GetCommonKmerHeadThread(void * arg){
	GetCommonKmerHeadP * getCommonKmerHeadP = (GetCommonKmerHeadP *)arg;
	 
	KmerHashTableHead * kmerHashTableHead = getCommonKmerHeadP->kmerHashTableHead;
    KmerReadNodeHead * kmerReadNodeHead = getCommonKmerHeadP->kmerReadNodeHead;
	ReadSetHead * readSetHead = getCommonKmerHeadP->readSetHead;
    int kmerLength = getCommonKmerHeadP->kmerLength;
	char * readFile = getCommonKmerHeadP->readFile;
    char * outputFile = getCommonKmerHeadP->outputFile;
	int step = getCommonKmerHeadP->step;
    int threadIndex = getCommonKmerHeadP->threadIndex;
	int totalThreadNumber = getCommonKmerHeadP->totalThreadNumber;

	
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
	for(int n=0;n<G->allocationCountGraph;n++){
		G->graph[n].dataLeft = -1;
		G->graph[n].dataRight = -1;
		G->graph[n].eadgCount = 0;
		G->graph[n].visit = 0;
		G->reverseGraph[n].dataLeft = -1;
		G->reverseGraph[n].dataRight = -1;
		G->reverseGraph[n].eadgCount = 0;
		G->reverseGraph[n].visit = 0;
		G->arcIndex[n].startIndex = 0;
		G->arcIndex[n].endIndex = 0;
	}
	
	
	
	int maxCount = 100000;
	
	CommonKmerHead * commonKmerHead = (CommonKmerHead *)malloc(sizeof(CommonKmerHead));
	commonKmerHead->commonKmer = (CommonKmer *)malloc(sizeof(CommonKmer)*maxCount);
	commonKmerHead->realCount = 0;
	commonKmerHead->allocationCount = maxCount;
	
	unsigned int kmerInteger = 0;
	long int hashIndex;
	unsigned long  int readIndex = 0;
	long int Size = 200000;
	long int readLength = 0;
	char str[Size];
	char * kmer1 = (char *)malloc(sizeof(char)*(kmerLength + 1));
	char * kmer2 = (char *)malloc(sizeof(char)*(kmerLength + 1));
	FILE * fp = fopen(readFile,"r");
	FILE * fp1 = fopen(outputFile,"w");
	while(!feof(fp)){
		fgets(str, Size, fp);
		if(str[0]=='>'){
			readIndex++;
			commonKmerHead->readIndex = readIndex;
			continue;
		}
		if(readIndex % totalThreadNumber != threadIndex){
			continue;
		}
		
		readLength = strlen(str);
		for(int j = 0; j < readLength - kmerLength;j = j+step){
			strncpy(kmer1,str + j, kmerLength);
			kmer1[kmerLength] = '\0';
			SetBitKmer(&kmerInteger, kmerLength, kmer1);
			hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
			if(hashIndex!=-1){
				InsertCommonToTwoReadAligningHead(commonKmerHead, kmerReadNodeHead, kmerHashTableHead, hashIndex, readIndex);
			}else{
				ReverseComplementKmer(kmer1, kmerLength);
				SetBitKmer(&kmerInteger, kmerLength, kmer1);
				hashIndex = SearchKmerHashTable(kmerHashTableHead, kmerInteger);
					InsertCommonToTwoReadAligningHead(commonKmerHead, kmerReadNodeHead, kmerHashTableHead, hashIndex, readIndex);
			}	
		}
		sort(commonKmerHead->commonKmer, 0, commonKmerHead->realCount - 1);
		
		GetOverlapResult(G, commonKmerHead, readSetHead, fp1);
		
		commonKmerHead->realCount = 0;
	}
	fclose(fp);
	fflush(fp1);
	fclose(fp1);
	
}


void sortGraph(AdjGraph * graph, int left, int right){
	if(left >= right){
        return ;
    }
    unsigned long  int i = left;
    unsigned long  int j = right;
    unsigned long  int key = graph[left].dataLeft;
	unsigned long  int dataRight = graph[left].dataRight;
	unsigned long  int eadgCount = graph[left].eadgCount;
	
    while(i < j){
        while(i < j && key <= graph[j].dataLeft){
            j--;
        }
		
        graph[i].dataLeft = graph[j].dataLeft;
		graph[i].dataRight = graph[j].dataRight;
		graph[i].eadgCount = graph[j].eadgCount;
        
         
        while(i < j && key >= graph[i].dataLeft){
            i++;
        }
       
        graph[j].dataLeft = graph[i].dataLeft;
		graph[j].dataRight = graph[i].dataRight;
		graph[j].eadgCount = graph[i].eadgCount;
    }
    
    graph[i].dataLeft = key;
	graph[i].dataRight = dataRight;
	graph[i].eadgCount = eadgCount;
	
    sortGraph(graph, left, i - 1);
    sortGraph(graph, i + 1, right);

}


void sort(CommonKmer *a, int left, int right)
{
    if(left >= right){
        return ;
    }
	
    unsigned long  int i = left;
    unsigned long  int j = right;
    unsigned long  int key = a[left].readIndex;
	unsigned long  int leftPosition = a[left].leftPosition;
	unsigned long  int rightPosition = a[left].rightPosition;
    bool orientation = a[left].orientation;
	
    while(i < j){
        while(i < j && key <= a[j].readIndex){
            j--;
        }
		
        a[i].readIndex = a[j].readIndex;
		a[i].leftPosition = a[j].leftPosition;
		a[i].rightPosition = a[j].rightPosition;
		a[i].orientation = a[j].orientation;
      
        while(i < j && key >= a[i].readIndex){
            i++;
        }
       
        a[j].readIndex = a[i].readIndex;
		a[j].leftPosition = a[i].leftPosition;
		a[j].rightPosition = a[i].rightPosition;
		a[j].orientation = a[i].orientation;
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

double GetLongestPathInGraph(AdjGraphHead * G, bool orien){
	
	AdjGraph * graph = NULL;
	int realCountGraph = 0;
	if(orien == 0){
		graph = G->graph;
		realCountGraph = G->realCountGraph;
	}else{
		graph = G->reverseGraph;
		realCountGraph = G->reverseRealCountGraph;
	}
	
	if(G->nodeCount < realCountGraph){
		G->distanceToSource = (int*)malloc(sizeof(int)*realCountGraph);
		G->edgeToSource = (int*)malloc(sizeof(int)*realCountGraph);
		G->nodeCount = realCountGraph;
	}
	for(long int i = 0; i < realCountGraph; i++){
		G->distanceToSource[i] = -100;
		G->edgeToSource[i] = -100;
	}
	long int maxValue = -1;
	G->leftStart = 0;
	G->rightStart = 0;
	G->leftEnd = 0;
	G->rightEnd = 0;
	
	long int maxIndex = 0;
	long int previousIndex = -1;
	long int maxOverlapLen = 0;
	long int minOverlapLen = 0;
	double ratio = 100;
	int cc = 0;
	for(long int i = 0; i< realCountGraph && cc < 4; i++){	
		G->distanceToSource[i] = 0;
		if(i != 0 && graph[i].dataLeft - graph[i - 1].dataLeft != 1){
			cc++;
		}
		for(long int j = 0; j < G->realCountArc; j++){	
			if(G->arcIndex[j].startIndex < i){
				continue;
			}
			if(G->distanceToSource[G->arcIndex[j].endIndex] < G->distanceToSource[G->arcIndex[j].startIndex] + 1){
				G->distanceToSource[G->arcIndex[j].endIndex] = G->distanceToSource[G->arcIndex[j].startIndex] + 1;
				G->edgeToSource[G->arcIndex[j].endIndex] = G->arcIndex[j].startIndex;
			}
		}
		bool token = false;
		for(long int p = i; p < realCountGraph; p++){
			if(G->distanceToSource[p] > maxValue){
				maxValue = G->distanceToSource[p];
				maxIndex = p;
				if(orien == 0){
					G->leftStart = graph[i].dataLeft;
					G->rightStart = graph[i].dataRight;
					G->leftEnd = graph[p].dataLeft;
					G->rightEnd = graph[p].dataRight;
				}else{
					G->leftStart = graph[i].dataLeft;
					G->rightStart = graph[p].dataRight;
					G->leftEnd = graph[p].dataLeft;
					G->rightEnd = graph[i].dataRight;
				}
				token = true;
			}
		}
		if(token != false){
			ratio = 1;
		}

		
		for(long int i = 0; i < realCountGraph; i++){
			G->distanceToSource[i] = -100;
			G->edgeToSource[i] = -100;
		}
	}
	return ratio;
}


int Overlap_Display(AdjGraphHead * G,int leftIndex,int rightIndex,bool orien,int leftLen,int rightLen,FILE * fp){
	
	AdjGraph * graph = NULL;
	int realCountGraph = 0;
	if(orien == 0){
		graph = G->graph;
		realCountGraph = G->realCountGraph;
	}else{
		graph = G->reverseGraph;
		realCountGraph = G->reverseRealCountGraph;
	}
	
	int i;
	int t = 0;
	int temp1;
	int count,temp=5,F=0,leftemp,rightemp;
	int leftStartpos=-1,leftEndpos=-1,rightStartpos=-1,rightEndpos=-1;
	int stempleftStartpos,stempleftEndpos,stemprightStartpos,stemprightEndpos;
	
	int overlenleft,overlenright;
	int MinOverLen,MaxOverLen;
	int max_i = 0,keep_i=-1;
	int temp2=0;
	int left_x,right_y;
	double SumXY=1,SumX=1,SumY=1;
	double Similarity;
	
	double ratio = GetLongestPathInGraph(G, orien);
	
	leftStartpos = G->leftStart;
	leftEndpos = G->leftEnd;
	rightStartpos = G->rightStart;
	rightEndpos = G->rightEnd;
	
	MinOverLen = min(leftEndpos -leftStartpos,rightEndpos - rightStartpos);
	MaxOverLen = max(leftEndpos -leftStartpos,rightEndpos - rightStartpos);
	
	if(MinOverLen < 500){return 0;} 
	
	all++;
	t=800;
	if(orien==0){
		
		if(rightStartpos > leftStartpos){
			if(leftStartpos > t){
				yyCount0++; 
				return 0;
			}
			rightStartpos = rightStartpos - leftStartpos;
			leftStartpos = 1;
		}else{
			if(rightStartpos > t){
				yyCount0++; 
				return 0;
			}
			leftStartpos = leftStartpos-rightStartpos;
			rightStartpos = 1;
		}
		
		if(leftLen-leftEndpos > rightLen-rightEndpos){
			if(rightLen-rightEndpos > t){
				yyCount0++; 
				return 0;
			}
			leftEndpos = rightLen - rightEndpos + leftEndpos;
			rightEndpos = rightLen;
			
		}else{
			if(leftLen-leftEndpos > t){
				yyCount0++; 
				return 0;
			}
			rightEndpos = leftLen - leftEndpos + rightEndpos;
			leftEndpos = leftLen;
		}
		
	}
	if(orien==1){
		
		if(rightLen-rightEndpos > leftStartpos){
			if(leftStartpos > t){
				yyCount0++; 
				return 0;
			}
			rightEndpos = rightEndpos + leftStartpos;
			leftStartpos = 1;
		}else{
			if(rightLen-rightEndpos > t){
				yyCount0++; 
				return 0;
			}
			leftStartpos = leftStartpos - rightLen + rightEndpos;
			rightEndpos = rightLen;
		}
		
		if(leftLen-leftEndpos > rightStartpos){
			if(rightStartpos > t){
				yyCount0++; 
				return 0;
			}
			leftEndpos = leftEndpos + rightStartpos;
			rightStartpos = 1;
		}else{
			if(leftLen-leftEndpos > t){
				yyCount0++; 
				return 0;
			}
			leftEndpos = leftLen;
			rightStartpos = rightStartpos - leftLen + leftEndpos;
			
		}

	}
	
	overlenleft = leftEndpos - leftStartpos;
	overlenright = rightEndpos - rightStartpos;
	MinOverLen = min(overlenleft,overlenright);
	MaxOverLen = max(overlenleft,overlenright);
	
	if(MinOverLen<500){ return 0;}
	if(MinOverLen*3 <= MaxOverLen *2){ return 0;} 
	
	if(leftIndex > rightIndex){
		fprintf(fp,"%d,%d,%d,%d,%d,%d,%d,%d,%d\n",leftIndex,rightIndex,orien,
												leftStartpos,leftEndpos,leftLen,rightStartpos,rightEndpos,rightLen);
	}else{
		fprintf(fp,"%d,%d,%d,%d,%d,%d,%d,%d,%d\n",rightIndex,leftIndex,orien,
												rightStartpos,rightEndpos,rightLen,leftStartpos,leftEndpos,leftLen);
	}
	
	return 1;

}

void ReAllocateAdjGraph(AdjGraphHead * G, int a){
	
	if(a == 0){
		G->allocationCountGraph = G->allocationCountGraph*1.5;
		AdjGraph * graph = (AdjGraph *)malloc(sizeof(AdjGraph)*G->allocationCountGraph);
		for(long int i = 0; i < G->realCountGraph; i++){
			graph[i].dataLeft = G->graph[i].dataLeft;
			graph[i].dataRight = G->graph[i].dataRight;
			graph[i].eadgCount = G->graph[i].eadgCount;
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
			graph[i].eadgCount = G->reverseGraph[i].eadgCount;
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

int CreatGraph(AdjGraphHead * G, CommonKmer * commonKmer, unsigned long  int startIndex, unsigned long  int endIndex, int a) {
	int count,throld;
	unsigned long  int i,j,k,m,n;
	int n1,n2,m1,m2;
	count = 0;
	int minvalue;
	AdjGraph * graph = NULL;
	int realCountGraph = 0;
	if(a == 0){
		graph = G->graph;
		realCountGraph = G->realCountGraph;
	}else{
		graph = G->reverseGraph;
		realCountGraph = G->reverseRealCountGraph;
	}
	
	sortGraph(graph, 0, realCountGraph - 1);

	for(j=0; j<realCountGraph; ++j){
		
		for(k=j+1; k<realCountGraph && k < j + 10 ; ++k){
			
			n1 = graph[j].dataLeft;
			n2 = graph[j].dataRight;
			m1 = graph[k].dataLeft;
			m2 = graph[k].dataRight;
			if(a == 0){
				if(n1<m1 && n2<m2){
					minvalue = max((m1-n1),(m2-n2));
					throld = 150 + minvalue*0.15;
					if(throld>500){throld = 500;}
					if(abs((m1-n1)-(m2-n2))< throld){
						G->arcIndex[G->realCountArc].startIndex = j;
						G->arcIndex[G->realCountArc].endIndex = k;
						G->realCountArc++;
						graph[j].eadgCount++;
						graph[k].visit = 0;
						if(G->realCountArc >= G->allocationCountArc){
							ReAllocateArcIndex(G);
						}
					}
				}else if(n1>m1 && n2>m2){
					minvalue = max((n1-m1),(n2-m2));
					throld = 150 + minvalue*0.15;
					if(throld>500){throld = 500;}
					if(abs((n1-m1)-(n2-m2))< throld){
						G->arcIndex[G->realCountArc].startIndex = j;
						G->arcIndex[G->realCountArc].endIndex = k;
						G->realCountArc++;
						graph[j].eadgCount++;
						graph[k].visit = 0;    
						if(G->realCountArc >= G->allocationCountArc){
							ReAllocateArcIndex(G);
						}
					}
				}else{continue;}
			}
			if(a == 1){
				if(n1<m1 && n2>m2){
					minvalue = max((m1-n1),(n2-m2));
					throld = 150 + minvalue*0.15;
					if(throld>500){throld = 500;}
					if(abs((m1-n1)-(n2-m2))<throld){
						G->arcIndex[G->realCountArc].startIndex = j;
						G->arcIndex[G->realCountArc].endIndex = k;
						G->realCountArc++;
						graph[j].eadgCount++;
						graph[k].visit = 0;
						if(G->realCountArc >= G->allocationCountArc){
							ReAllocateArcIndex(G);
						}
					}
				}else if(n1>m1 && n2<m2){
					minvalue = max((n1-m1),(m2-n2));
					throld = 150 + minvalue*0.15;
					if(throld>500){throld = 500;}
					if(abs((n1-m1)-(m2-n2))<throld){
						G->arcIndex[G->realCountArc].startIndex = j;
						G->arcIndex[G->realCountArc].endIndex = k;
						G->realCountArc++;
						graph[j].eadgCount++;
						graph[k].visit = 0;
						if(G->realCountArc >= G->allocationCountArc){
							ReAllocateArcIndex(G);
						}
					}
				}else{continue;}
			}					
		}
	}
	
	
}

void GetOverlapResult(AdjGraphHead * G, CommonKmerHead * commonKmerHead, ReadSetHead * readSetHead,FILE * fp){
	int leftIndex, rightIndex, orien;
	long int leftLen,rightLen;
	long int realCommonKmerCount;
	int s = 0 ;
	unsigned long  int readIndex = commonKmerHead->commonKmer[0].readIndex;
	unsigned long  int startIndex = 0;
	unsigned long  int endIndex = 0;
	bool orientation;
	
	for(long int m = 0; m < commonKmerHead->realCount; m++){
		if(commonKmerHead->commonKmer[m].readIndex != readIndex){
			endIndex = m - 1;
			if(G->realCountGraph < 10 && G->reverseRealCountGraph < 10){
				G->realCountGraph = 0;
				G->realCountArc = 0;
				G->reverseRealCountGraph = 0;
				
				if(commonKmerHead->commonKmer[m].orientation == 0){
					G->graph[G->realCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition;
					G->graph[G->realCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition;
					G->graph[G->realCountGraph].eadgCount = 0;
					G->graph[G->realCountGraph].visit = 1;
					G->realCountGraph++;
					if(G->realCountGraph >= G->allocationCountGraph){
						ReAllocateAdjGraph(G, 0);
					}
				}else{
					G->reverseGraph[G->reverseRealCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition;
					G->reverseGraph[G->reverseRealCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition;
					G->reverseGraph[G->reverseRealCountGraph].eadgCount = 0;
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
			leftIndex = commonKmerHead->readIndex;
			rightIndex = readIndex;
			leftLen = readSetHead->readSet[leftIndex - 1].readLength;
			rightLen = readSetHead->readSet[rightIndex - 1].readLength;
			
			CreatGraph(G,commonKmerHead->commonKmer,startIndex, endIndex, orientation);
			
			int result = Overlap_Display(G,leftIndex,rightIndex,orientation,leftLen,rightLen,fp);
			
			G->realCountGraph = 0;
			G->realCountArc = 0;
			G->reverseRealCountGraph = 0;
			if(result != 1){
				CreatGraph(G,commonKmerHead->commonKmer,startIndex, endIndex, !orientation);
				Overlap_Display(G,leftIndex,rightIndex,!orientation,leftLen,rightLen,fp);
				G->realCountGraph = 0;
				G->realCountArc = 0;
				G->reverseRealCountGraph = 0;
			}
			startIndex = m;
			readIndex = commonKmerHead->commonKmer[m].readIndex;
		}
			
		if(commonKmerHead->commonKmer[m].orientation == 0){
			G->graph[G->realCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition;
			G->graph[G->realCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition;
			G->graph[G->realCountGraph].eadgCount = 0;
			G->graph[G->realCountGraph].visit = 1;
			G->realCountGraph++;
			if(G->realCountGraph >= G->allocationCountGraph){
				ReAllocateAdjGraph(G, 0);
			}
		}else{
			G->reverseGraph[G->reverseRealCountGraph].dataLeft = commonKmerHead->commonKmer[m].leftPosition;
			G->reverseGraph[G->reverseRealCountGraph].dataRight = commonKmerHead->commonKmer[m].rightPosition;
			G->reverseGraph[G->reverseRealCountGraph].eadgCount = 0;
			G->reverseGraph[G->reverseRealCountGraph].visit = 1;
			G->reverseRealCountGraph++;
			if(G->reverseRealCountGraph >= G->reverseAllocationCountGraph){
				ReAllocateAdjGraph(G, 1);
			}
		}
		
	}
}




#endif
