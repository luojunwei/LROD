#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
//#define CATCH_CONFIG_RUNNER
#include "/usr/local/include/catch2/catch.hpp"

#include <unistd.h>

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
 
 
#include "../kmer.h"
#include "../read.h"
#include "../aligning.h"
#include "../bitarray.h"
 
using namespace std;

void print_usage(){
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
 
int OverlapDetect(int argc, char * argv[]){
    if(argc < 3){
		printf("Missing long read file or k-mer frequency file!\n");
		printf("For detailed usage of LROD, please command: -h or -help!\n");
		return 1;
	}
	long int maxSize = 1000000;
	char * StrLine = (char *)malloc(sizeof(char)*maxSize);

	  
	char * readFile = argv[1];
	char * binaryKmerFile = argv[2];
	char * outputKmerFile = (char *)malloc(sizeof(char)*500);
	FILE * fp4;
	strcpy(outputKmerFile, "output");
	
	long int step = 1;  //Step length when extracting kmer from reads
	long int kmerLength = atoi(argv[3]);  //Extract the length of kmer from the reads
	long int smallKmerLength = 9;  //Kmer length in the second stage of detection
	long int threadCount = 1;  //Number of threads
	
	long int smallIntervalDistance = 400;
	long int largeIntervalDistance = 1500;
	long int overlapLengthCutOff = 500;
	float lengthRatio = 0.3;
	int frequencyCutOff = 3;
	long int minimumKmerFrequency = 2;  
	float maxKmerFrequencyRatio = 0.9;
 
    
    
	
	if((fp4 = fopen(readFile, "r")) == NULL){
		printf("Missing long read file!\n");
		printf("For detailed usage of LROD, please command: -h or -help!\n");
        return 2;
    }
	fclose(fp4);
	
	if((fp4 = fopen(binaryKmerFile, "r")) == NULL){
		printf("Missing kmer frequency file!\n");
		printf("For detailed usage of LROD, please command: -h or -help!\n");
		return 3;
    }
	fclose(fp4);
	
	strcat(outputKmerFile, ".csv");
	if((fp4 = fopen(outputKmerFile, "w")) == NULL){
        printf("%s, does not exist!\n", outputKmerFile);
        return 4;
    }
	fclose(fp4);
	
    
	printf("\nStart to load long reads!\n");
	//Store long reads in an array
	ReadSetHead * readSetHead = GetReadSetHead(readFile,StrLine,maxSize);
	
	if(readSetHead->readCount <= 1){
		printf("The number of reads is smaller than one!\n");
        return 5;
	}
	printf("Start to construct kmer hash table!\n");
	
	KmerHashTableHead * kmerHashTableHead = GetKmerHashTableHead(binaryKmerFile, readSetHead, kmerLength, step, minimumKmerFrequency, maxKmerFrequencyRatio);
	if(GetKmerHashTableHead_UnitTest(kmerHashTableHead) == 0){
		return 6;
	}
	long int subReadCount = 50000;
	//Create a kmer hash table space based on the reads length
	KmerReadNodeHead * kmerReadNodeHead = GetKmerReadNodeHeadSub(readSetHead, kmerLength, step, subReadCount);
	
	printf("Prepare to detect overlap among long reads!\n");
	//Start multi-thread sequence alignment
	GetCommonKmerHeadAllThreadNew(kmerHashTableHead, kmerReadNodeHead, readSetHead, kmerLength, readFile, outputKmerFile, step, threadCount, smallKmerLength, smallIntervalDistance, largeIntervalDistance, overlapLengthCutOff, lengthRatio, subReadCount);
	
	if(GetCommonKmerHeadAllThreadNew_UnitTest(outputKmerFile, readSetHead->readCount) == 0){
		return 7;
	}
	
	printf("done!\n");
	printf("Result file name is: %s!\n",outputKmerFile);
	 
    return 0;
}


TEST_CASE( "Normal Test", "[OverlapDetect]") {
	
	char * argv1[4];
	argv1[0] = (char *)malloc(sizeof(char)*500);
	argv1[1] = (char *)malloc(sizeof(char)*500);
	argv1[2] = (char *)malloc(sizeof(char)*500);
	argv1[3] = (char *)malloc(sizeof(char)*500);

	strcpy(argv1[0], "ss");
	strcpy(argv1[1], "../test/long_read.fa");
	strcpy(argv1[2], "../test/kmer_file.txt");
	strcpy(argv1[3], "15");
	
	REQUIRE( OverlapDetect(9, argv1) == 0 );
	
	free(argv1[0]);
	free(argv1[1]);
	free(argv1[2]);
	free(argv1[3]);

}


TEST_CASE( "Read file test", "[OverlapDetect2]") {
	
	char * argv1[4];
	argv1[0] = (char *)malloc(sizeof(char)*500);
	argv1[1] = (char *)malloc(sizeof(char)*500);
	argv1[2] = (char *)malloc(sizeof(char)*500);
	argv1[3] = (char *)malloc(sizeof(char)*500);

	strcpy(argv1[0], "ss");
	strcpy(argv1[1], "./data/read_1.fa");
	strcpy(argv1[2], "./data/kmer_frequency_1.txt");
	strcpy(argv1[3], "15");
	
	REQUIRE( OverlapDetect(9, argv1) == 5 );
	
	free(argv1[0]);
	free(argv1[1]);
	free(argv1[2]);
	free(argv1[3]);

}


TEST_CASE( "Normal test 1", "[OverlapDetect2]") {
	
	char * argv1[4];
	argv1[0] = (char *)malloc(sizeof(char)*500);
	argv1[1] = (char *)malloc(sizeof(char)*500);
	argv1[2] = (char *)malloc(sizeof(char)*500);
	argv1[3] = (char *)malloc(sizeof(char)*500);

	strcpy(argv1[0], "ss");
	strcpy(argv1[1], "./data/read_100.fa");
	strcpy(argv1[2], "./data/kmer_frequency_100.txt");
	strcpy(argv1[3], "15");
	
	REQUIRE( OverlapDetect(9, argv1) == 0 );
	
	free(argv1[0]);
	free(argv1[1]);
	free(argv1[2]);
	free(argv1[3]);

}

TEST_CASE( "Kmer file NULL", "[OverlapDetect2]") {
	
	char * argv1[4];
	argv1[0] = (char *)malloc(sizeof(char)*500);
	argv1[1] = (char *)malloc(sizeof(char)*500);
	argv1[2] = (char *)malloc(sizeof(char)*500);
	argv1[3] = (char *)malloc(sizeof(char)*500);

	strcpy(argv1[0], "ss");
	strcpy(argv1[1], "./data/read_100.fa");
	strcpy(argv1[2], "./data/kmer_frequency_xx.txt");
	strcpy(argv1[3], "15");
	
	REQUIRE( OverlapDetect(9, argv1) == 3 );
	
	free(argv1[0]);
	free(argv1[1]);
	free(argv1[2]);
	free(argv1[3]);

}

TEST_CASE( "Read file NULL", "[OverlapDetect2]") {
	
	char * argv1[4];
	argv1[0] = (char *)malloc(sizeof(char)*500);
	argv1[1] = (char *)malloc(sizeof(char)*500);
	argv1[2] = (char *)malloc(sizeof(char)*500);
	argv1[3] = (char *)malloc(sizeof(char)*500);

	strcpy(argv1[0], "ss");
	strcpy(argv1[1], "./data/read_xx.fa");
	strcpy(argv1[2], "./data/kmer_frequency_100.txt");
	strcpy(argv1[3], "15");
	
	REQUIRE( OverlapDetect(9, argv1) == 2 );
	
	free(argv1[0]);
	free(argv1[1]);
	free(argv1[2]);
	free(argv1[3]);

}












