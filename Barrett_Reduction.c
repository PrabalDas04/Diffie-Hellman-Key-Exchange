#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Barrett_Reduction.h"

////////////////////////////////// Function Definition ///////////////////////////////////

// calculate power of 2
unsigned long long power_of_two(int n)
{
    return 1ULL << n;  // Left shift 1 by n positions
}

// parse the string into array of cell size 29-bit
// output : bitArray[0,1,...,8]
void ParsingArray(long long bitArray[], char bitString[], FILE *file, int l)
{
	int i, k = 0, j = 28;
	
    // Ensure the bit string is no longer than l bits
    // (l+1) is sizeof(bitString)
    if (fgets(bitString, l+1, file) == NULL || strlen(bitString) != l) {
        fprintf(stderr, "Invalid bit string format in file.\n");
        fclose(file);
        exit(0);
    }
	
	// iteratively take 29 bits and store the corresponding integer value
    for (i = 255; i > 23 ; i--)
    {
    	// get the integer rep of the 29 bits in the k-cell
		bitArray[k] = bitArray[k] + power_of_two(28-j) * (bitString[i] - '0');
		j--;
		if(j == -1)
		{
			k++;
			j = 28;
		}
     }

     // last 24 bits are kept in 8-th cell
    for (i = 23; i >= 0; i--)
    	bitArray[8] = bitArray[8] + power_of_two(23-i) * (bitString[i] - '0');   
}

// school book multiplication
// output : C = A * B
void MultSchoolBook(long long A[], int len_A, long long B[], int len_B, long long C[])
{
    for (int i = 0; i < len_A; i++)
    {
        for (int j = 0; j < len_B; j++)
        {
            C[i + j] += A[i] * B[j]; // Add product of coefficients
        }
    }
}

// function for kepping the string into 29-bit size cells
// l is the degree of C
// output : C[0,1,...,l+1]
void ConversionTo29bit(long long C[], int l)
{
	long long temp = C[0];
    int i;
    
    for(i = 0; i < l; i++)
    {
    	// extract the last 29 bits
    	C[i] = temp & 0x1FFFFFFF;
    	// shift temp to 29 bits and store in temp
    	temp >>= 29;
    	// addtion with temp
    	temp = temp + C[i+1];
    }
    C[l] = temp & 0x1FFFFFFF;
    temp >>= 29;
    C[l+1] = temp;
}

// Function to compare two arrays
// output : 1 if a > b, -1 if a < b and 0 if a = b
int compare(long long *a, long long *b, int len)
{
    for (int i = len - 1; i >= 0; i--) 
    {
        if (a[i] > b[i]) return 1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

// function for printing array arr
void printArray(long long arr[], int size)
{
    for (int i = 0; i < size; i++)
    {
        printf("%d-th block : %lld\n", i, arr[i]);
    }
    printf("\n");
}


// Function to divide an array by b^exp
void divide_by_base_exp(long long *num, int num_len, long long *result, int exp)
{
    // Shift right by 29 * exp bits (equivalent to dividing by b^exp)
    for (int i = exp; i < num_len; i++)
        result[i - exp] = num[i];
}

// Function to subtract two arrays, where A > B
// output : A <-- (A - B)
void subtract(long long *A, long long *B, int len)
{
    int carry = 0;
    
    for(int i = 0; i < len; i++)
    {
    	if(A[i] < B[i] + carry)
    	{
    		A[i] = (A[i] + (MASK + 1)) - B[i] - carry;
    		carry = 1;
    	}
    	else{
    		A[i] -= B[i] + carry;
    		carry = 0;
    	}
    }
}

// subtraction function where A and B arbitrary
// output : A <-- (A - B)
void subtract_any(long long *A, long long *B)
{
	// q is the 256-bit prime, which is same as "prime"
	long long q[MAX_SIZE] = {535425013,174332635,444665496,192778653,388389189,518147849,304619691,363717891,15281728,0,0,0,0,0,0,0,0,0,0,0};
	
	// If A < B, then enter 'If'
	if(compare(A, B, MAX_SIZE) < 0)
	{
		// do (B - A) --> B
		subtract(B, A, MAX_SIZE);
		// do (q - (B - A)) --> q
		subtract(q, B, MAX_SIZE);
		// copy the result in array A
		for(int i = 0; i < MAX_SIZE; i++)
			A[i] = q[i];
	}
	else
		subtract(A, B, MAX_SIZE);
}

// addition function for 256-bit int
// output : C <-- (A + B)
void Addition_256(long long *A, long long *B, long long *C)
{
	long long TString[MAX_SIZE] = 
 	{450887704,490307913,387807083,403879883,291135210,307268612,110539282,24605042,70628772,35,0,0,0,0,0,0,0,0,0,0};
 	
 	long long prime[MAX_SIZE] = 
 	{535425013,174332635,444665496,192778653,388389189,518147849,304619691,363717891,15281728,0,0,0,0,0,0,0,0,0,0,0};
	long long result[MAX_SIZE] = {0};
	// component wise addition
	for(int i = 0; i < MAX_SIZE; i++)
		result[i] = A[i] + B[i];
		
	// conversion to base
	ConversionTo29bit(result,8);
	//initialization to 0
	InitializeToZero(C);
	// barrett reduction
	barrett_reduction(result, 10, prime, MAX_SIZE, C, TString);
}

// Function to compute num mod b^exp
void modulo_by_base_exp(long long *num, int num_len, long long *result, int exp)
{
    // Only keep the lowest `exp` "digits"
    for (int i = 0; i < exp; i++)
        result[i] = num[i];
}

// function for barrett reduction
// output : result = x (mod m)
void barrett_reduction(long long *x, int x_len, long long *m, int m_len, long long *result, long long *mu)
{
    // initialization to 0
    long long q1[MAX_SIZE] = {0};
    long long q2[MAX_SIZE] = {0};
    long long q3[MAX_SIZE] = {0};
    long long r1[MAX_SIZE] = {0};

    // Step 1: Compute q1 = floor(x / b^(L-1))
    divide_by_base_exp(x, MAX_SIZE, q1, L - 1);

    // Step 2: Compute q2 = q1 * mu
    // length of q1 = 10 = mu(or T)
    MultSchoolBook(q1, 10, mu, 10, q2);
    // conversion into 29-bit size cell
     // degree of q2 is 18
     ConversionTo29bit(q2, 18);

    // Step 3: Compute q3 = floor(q2 / b^(L+1))
    divide_by_base_exp(q2, MAX_SIZE, q3, L + 1);

    // Step 5: Compute r1 = (q3 * m)
    MultSchoolBook(q3, 10, m, m_len, r1);
    // conversion into 29-bit size cell
     // degree of r1 is 17
     ConversionTo29bit(r1, 17);

    // Step 6: Compute r = x - r1
    subtract_any(x, r1);

    // Step 8: While r >= m, reduce r by m
    while (compare(x, m, MAX_SIZE) >= 0)
        subtract_any(x, m);

    // Copy the result
    for (int i = 0; i < MAX_SIZE; i++) {
        result[i] = x[i];
    }  
}

// function for converting integer to binary
// output: array 'bin', binary rep of 'num'
void IntToBin(int *bin, int num)
{
    for(int i = 0; i < 29; i++)
    	bin[i] = 0;
    for (int i = 28; i >= 0; i--)
        bin[i] = (num >> i) & 1;  // Extract the bit at position i
}

// function for initializing the array arr to 0
void InitializeToZero(long long *arr)
{
	for(int j = 0; j < MAX_SIZE; j++)
 		arr[j] = 0;
}
