#ifndef BITS_H_INCLUDED 
#define BITS_H_INCLUDED 

#define BASE 2 
#define SHIFT 2 
#define MASK 3 
  
void SetBitTemp(unsigned long int * bit_array, unsigned long int bit_number, char value);

void SetBitKmer(unsigned long int * bit_array, long int len, char * value);

char GetBit(unsigned long long int * bit_array, unsigned int bit_number);

void GetBit(unsigned long long int * bit_array, unsigned int bit_number, int len, char * value);

void SetBit(char bit_array[], unsigned int bit_number, char value);

char GetBit(char bit_array[], unsigned int bit_number);

void SetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length, char * value) ;

char * GetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length);

void GetBit(char bit_array[], unsigned int bit_start_number, unsigned int bit_length, char * value);

#endif 
