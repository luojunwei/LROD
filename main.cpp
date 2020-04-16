
 #include<unistd.h>
 
 
 #include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<iostream>
#include<ctype.h>
#include <time.h>
#include <malloc.h>

#include <getopt.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
 
 
#include "kmer.h"
#include "read.h"
#include "aligning.h"
#include "bitarray.h"
 
 using namespace std;
 
 void print_usage()
 {
     printf("\nLROD can detect overlap regions among long reads.\n");
	 printf("\nThe command format is as follows: ");
     printf("\nLROD -r <long-read-file> -c <kmer-frequency-file> -o result-file [options]\n");
     printf("\nUsage:\n");
     printf("\t-r long-read-file: input file with fasta format\n");
     printf("\t-c kmer-frequency-file: each line in the kmer-frequency-file should be \"kmer kmer-frequency\"\n");
     printf("\t-o result-file: result file\n");
     printf("\t-t count: thread count (default 1)\n");
	 printf("\t-k kmerLength: kmer length (default 15)\n");
	 printf("\t-q smallKmerLength: small kmer length (default 9)\n");
	 printf("\t-f minimumKmerFrequency: minimum kmer frequency (default 2)\n");
	 printf("\t-m maxKmerFrequencyRatio: maximum kmer frequency ratio (it should be smaller than 1, by default 0.9)\n");
	 printf("\t-s maximumStep: kmer step (default 1)\n");
	 printf("\t-k kmerLength: kmer length (default 13)\n");
	 printf("\t-d distance: the small distance used for determining whether two common kmers are consistant (default 400)\n");
	 printf("\t-e distance: the large distance used for determining whether two common kmers are consistant (default 1500)\n");
	 printf("\t-a min-overlap-length: the minimum overlap length between two long reads (default 500)\n");
	 printf("\t-b length-ratio: the maximum length ratio between two aligned regions (default 0.3)\n");
     printf("\t-h, -help\n");
 }
 
 int main(int argc, char* argv[])
 {
    if(argc < 2){
		print_usage();
		exit(0);
	}
	long int maxSize = 1000000;
	char * StrLine = (char *)malloc(sizeof(char)*maxSize);

	  
	 char * readFile = NULL;
	char * binaryKmerFile = NULL;
	//char * outputKmerFile = NULL;
	char * outputKmerFile = (char *)malloc(sizeof(char)*500);
	FILE * fp4;
	strcpy(outputKmerFile, "output");
	
	long int step = 1;
	long int kmerLength = 15;
	long int smallKmerLength = 9;
	long int threadCount = 1;
	
	long int smallIntervalDistance = 400;
	long int largeIntervalDistance = 1500;
	long int overlapLengthCutOff = 500;
	float lengthRatio = 0.3;
	int frequencyCutOff = 3;
	long int minimumKmerFrequency = 2;
	float maxKmerFrequencyRatio = 0.9;
    
 
     struct option long_options[] = { 
         {"readFile", required_argument, NULL, 'r'},
         {"binaryKmerFile", required_argument, NULL, 'c'},
         {"outputKmerFile", optional_argument, NULL, 'o'},
		 
         {"step",  optional_argument, NULL, 's'},
         {"kmerLength", optional_argument, NULL, 'k'},
         {"smallKmerLength", optional_argument, NULL, 'q'},
		 {"minimumKmerFrequency", optional_argument, NULL, 'f'},
		 {"maxKmerFrequencyRatio", optional_argument, NULL, 'm'},
         {"threadCount", optional_argument, NULL, 't'},
         {"smallIntervalDistance", optional_argument, NULL, 'd'},
         {"largeIntervalDistance", optional_argument, NULL, 'e'},
         {"overlapLengthCutOff", optional_argument, NULL, 'a'},
         {"lengthRatio", optional_argument, NULL, 'b'},
         {"help", no_argument, NULL, 'h'},
         {0, 0, 0, 0}
     };  
     
  int ch = 0;
  while((ch = getopt_long(argc, argv, "c:r:o:m:n:d:k:e:a:s:t:f:q:b:h", long_options, NULL)) != -1) {
        
        switch (ch) {
            case 'r': readFile = (char *)(optarg); break;
			case 'c': binaryKmerFile = (char *)(optarg); break;
			case 'o': outputKmerFile = (char *)optarg; break;
			case 'k': kmerLength = atoi(optarg); break;
			case 'q': smallKmerLength = atoi(optarg); break;
			case 'f': minimumKmerFrequency = atoi(optarg); break;
			case 'm': maxKmerFrequencyRatio = atof(optarg); break;
			case 's': step = atoi(optarg); break;
			
			case 'd': smallIntervalDistance = atoi(optarg); break;
			case 'e': largeIntervalDistance = atoi(optarg); break;
			case 'a': overlapLengthCutOff = atoi(optarg); break;
			case 't': threadCount = atoi(optarg); break;
			case 'b': lengthRatio = atof(optarg); break;
            case 'h':
                print_usage();
                return 0;
            default:
                return -1;
        }
    }
    
	 
	
	if((fp4 = fopen(readFile, "r")) == NULL){
		printf("Missing long read file!\n");
		printf("\nFor detailed usage of LROD, please command: -h or -help!\n");
        exit(0);
    }
	fclose(fp4);
	
	if((fp4 = fopen(binaryKmerFile, "r")) == NULL){
		printf("Missing kmer frequency file!\n");
		printf("\nFor detailed usage of LROD, please command: -h or -help!\n");
		exit(0);
    }
	fclose(fp4);
	
	strcat(outputKmerFile, ".csv");
	if((fp4 = fopen(outputKmerFile, "w")) == NULL){
        printf("%s, does not exist!\n", outputKmerFile);
        exit(0);
    }
	fclose(fp4);
	
    
	printf("\nStart to load long reads!\n");
 
	ReadSetHead * readSetHead = GetReadSetHead(readFile,StrLine,maxSize);
	
	if(readSetHead->readCount <= 1){
		printf("The number of reads is smaller than one!\n");
        exit(0);
	}
	printf("Start to construct kmer hash table!\n");
	KmerHashTableHead * kmerHashTableHead = GetKmerHashTableHead(binaryKmerFile, readSetHead, kmerLength, step, minimumKmerFrequency, maxKmerFrequencyRatio);
	
	long int subReadCount = 50000;
	KmerReadNodeHead * kmerReadNodeHead = GetKmerReadNodeHeadSub(readSetHead, kmerLength, step, subReadCount);
	
	printf("Prepare to detect overlap among long reads!\n");
	GetCommonKmerHeadAllThreadNew(kmerHashTableHead, kmerReadNodeHead, readSetHead, kmerLength, readFile, outputKmerFile, step, threadCount, smallKmerLength, smallIntervalDistance, largeIntervalDistance, overlapLengthCutOff, lengthRatio, subReadCount);
	
	printf("done!\n");
	printf("Result file name is: %s!\n",outputKmerFile);
	 
    return 0;
 }