CC=g++

CPPFLAGS = -g -Wall -O3

LROD: main.o read.o kmer.o aligning.o bitarray.o
	$(CC) -o $@ $^ -lpthread
	
main.o: main.cpp read.h kmer.h aligning.h bitarray.h
	$(CC) -c main.cpp

bitarray.o: bitarray.cpp 
	$(CC) -c bitarray.cpp 

read.o: read.cpp read.h
	$(CC) -c read.cpp
	
sortContigSet.o: kmer.cpp read.h bitarray.h
	$(CC) -c kmer.cpp
	
aligning.o: aligning.cpp kmer.h read.h bitarray.h
	$(CC) -c aligning.cpp


	
all: LROD
	rm -f *.o

clean:
	rm -f *.o
	rm LROD
