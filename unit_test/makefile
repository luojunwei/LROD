CC=g++

CPPFLAGS = -g -Wall -O3

LROD_Unit_Test: unit_test.o read.o kmer.o aligning.o bitarray.o
	$(CC) -o $@ $^ -lpthread
	
unit_test.o: unit_test.cpp ../read.h ../kmer.h ../aligning.h ../bitarray.h
	$(CC) -c unit_test.cpp

bitarray.o: ../bitarray.cpp 
	$(CC) -c ../bitarray.cpp 

read.o: ../read.cpp ../read.h
	$(CC) -c ../read.cpp

kmer.o: ../kmer.cpp ../read.h ../bitarray.h
	$(CC) -c ../kmer.cpp
	
aligning.o: ../aligning.cpp ../kmer.h ../read.h ../bitarray.h
	$(CC) -c ../aligning.cpp


	
all: LROD_Unit_Test

	rm -f *.o

clean:
	rm -f *.o
	rm LROD_Unit_Test
