#include<stdio.h>

/////////////////////////// Global variables /////////////////////////////////
#define ARRAY_SIZE 9
#define L 9             // L = 9
#define MAX_SIZE 20     // Maximum number of "digits" for x and m arrays
#define BASE (1 << 29)  // Base b = 2^29
#define MASK ((1U << 29) - 1)

/////////////////////////// Function declaration ////////////////////////////

// Parses the input from a file or string into a 29-bit representation.
void ParsingArray(long long bitArray[], char bitString[], FILE *file, int l);
// Converts the result of arithmetic operations into a 29-bit representation.
void ConversionTo29bit(long long C[], int l);
// Multiplies two large integers using the schoolbook algorithm.
void MultSchoolBook(long long A[], int len_A, long long B[], int len_B, long long C[]);
// Subtracts two large integers represented as arrays.
void subtract(long long *A, long long *B, int len);
void subtract_any(long long *A, long long *B);
// Performs Barrett reduction for modular arithmetic.
void barrett_reduction(long long *x, int x_len, long long *m, int m_len, long long *result, long long *mu);
// function for computing power of 2
unsigned long long power_of_two(int n);
// function for print an array
void printArray(long long arr[], int size);
// function for divide an array with base^exp
void divide_by_base_exp(long long *num, int num_len, long long *result, int exp);
// function for addition of two 256-bit integer
void Addition_256(long long *A, long long *B, long long *C);
// function for reseting an array to 0
void InitializeToZero(long long *arr);
// function for converting an integer to binary
void IntToBin(int *bin, int num);
// function for comparing two array of integer
int compare(long long *a, long long *b, int len);
// function for mudulo base^exp
void modulo_by_base_exp(long long *num, int num_len, long long *result, int exp);