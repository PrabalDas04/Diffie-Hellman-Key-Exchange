#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "Barrett_Reduction.h"

long long* SelectArray(long long *A, long long *B, int b);
long long* SelectArrayWithInit(long long *A, long long *B, int b);
void L2R_Exponentiation(long long *prod, long long *prime, long long *TString, long long *num, long long *exp);
void R2L_Exponentiation(long long *prod, long long *prime, long long *TString, long long *num, long long *exp);
void Montogomery_Exponentiation(long long *prime, long long *TString, long long *num, char *exp);
void Montogomery_Exponentiation_WO_Branching(long long *prime, long long *TString, long long *num, char *exp);



// function for selecting an array based on bit b
long long* SelectArray(long long *A, long long *B, int b)
{
    // select Aor B according to b
    return (long long *)((uintptr_t)A * (1 - b) + (uintptr_t)B * b);
}

// function for selecting an array based on 'b' and then initialize it with 0
long long* SelectArrayWithInit(long long *A, long long *B, int b)
{
    // Select array pointer based on b
    long long *selected_array = (long long *)((uintptr_t)A * (1 - b) + (uintptr_t)B * b);

    // Initialize the selected array to zero
    for (int i = 0; i < MAX_SIZE; i++)
        selected_array[i] = 0;

    return selected_array;
}


// function to compute left to right exponentiation: num ^ exp
// output : prod <-- num ^ exp
void L2R_Exponentiation(long long *prod, long long *prime, long long *TString, long long *num, long long *exp)
{
	long long prod_temp[MAX_SIZE] = {0};  // temporary memory
	int bin[29] = {0};	// array for storing binary rep of integer
	int j = 8, count = 0, k = 28;
	
	// this loop will iterate 261 times(29*9), since during int to binary conversion each number is represented in 29 bit string
	for(int i = 0; i < 261; i++)
 	{
 		// doing prod * prod
 		MultSchoolBook(prod, 9, prod, 9, prod_temp);
 		ConversionTo29bit(prod_temp, 16);
 		// prod memory is set to 0
 		InitializeToZero(prod);
 		barrett_reduction(prod_temp, 18, prime, MAX_SIZE, prod, TString);
 		
 		// prod_temp memory is set to 0
 		InitializeToZero(prod_temp);
 		
 		// count for updating 'bin' array
 		if(count == 0)
 		{
 			count = 29;
 			// Integer to Binary conversion
 			IntToBin(bin, exp[j]);
 			// index for msb of exp[j]
 			k = 28;
 			j--;
 		}
 		
 		// checking the k-th bit of exp[j] is 1 or not
 		if(bin[k] == 1)
 		{
 			// doing prod * num
 			MultSchoolBook(prod, 9, num, 9, prod_temp);
 			ConversionTo29bit(prod_temp, 16);
 			barrett_reduction(prod_temp, 18, prime, MAX_SIZE, prod, TString);
 			
 			// prod_temp memory is set to 0
 			InitializeToZero(prod_temp);
 		}
 		
 		// decrement of counters
 		count--;
 		k--;
 	}
}


// function to compute right to left exponentiation: num ^ exp
// output : prod <-- num ^ exp
void R2L_Exponentiation(long long *prod, long long *prime, long long *TString, long long *num, long long *exp)
{
	long long prod_temp[MAX_SIZE] = {0};  // temporary memory
	int bin[29] = {0};	// array for storing binary rep of integer
	int j = 0, count = 0, k = 0;
	
	// this loop will iterate 261 times(29*9), since during int to binary conversion each number is represented in 29 bit string
	for(int i = 0; i < 261; i++)
 	{
 		// count for updating 'bin' array
 		if(count == 0)
 		{
 			count = 29;
 			// Integer to Binary conversion
 			IntToBin(bin, exp[j]);
 			// index for lsb of exp[j]
 			k = 0;
 			j++;
 		}
 		
 		// checking the k-th bit of exp[j] is 1 or not
 		if(bin[k] == 1)
 		{
 			// doing prod * num
 			MultSchoolBook(prod, 9, num, 9, prod_temp);
 			ConversionTo29bit(prod_temp, 16);
 			barrett_reduction(prod_temp, 18, prime, MAX_SIZE, prod, TString);
 			
 			// prod_temp memory is set to 0
 			InitializeToZero(prod_temp);
 		}
 		
 		// doing num * num
 		MultSchoolBook(num, 9, num, 9, prod_temp);
 		ConversionTo29bit(prod_temp, 16);
 		
 		// num memory is set to 0
 		InitializeToZero(num);
 		barrett_reduction(prod_temp, 18, prime, MAX_SIZE, num, TString);
 		
 		// prod_temp memory is set to 0
 		InitializeToZero(prod_temp);
 		
 		// decrese the value
 		count--;
 		k++;
 	}
}


// Montogomery exponentiation for num ^ exp
// output : num ^ exp
void Montogomery_Exponentiation(long long *prime, long long *TString, long long *num, char *exp)
{
	long long prod_temp[MAX_SIZE] = {0};  // temporary memory
	long long R[MAX_SIZE] = {0};	// initialization
	
	// R = num * num
	MultSchoolBook(num, 9, num, 9, prod_temp);
 	ConversionTo29bit(prod_temp, 16);
 	barrett_reduction(prod_temp, 18, prime, MAX_SIZE, R, TString);
 		
	// prod_temp memory is set to 0
 	InitializeToZero(prod_temp);
	
	// from 2nd msb to lsb, exp[0] is msb
	for(int i = 1; i < 256; i++)
	{
		if(exp[i] == '0')
 		{
 			// doing num * R
 			MultSchoolBook(num, 9, R, 9, prod_temp);
 			ConversionTo29bit(prod_temp, 16);
 			// R memory is set to 0, for safety o/w some cell may not be overwrite, keeping previous content
 			InitializeToZero(R);
 			barrett_reduction(prod_temp, 18, prime, MAX_SIZE, R, TString);
 			
 			// prod_temp memory is set to 0
 			InitializeToZero(prod_temp);
 			
 			// doing num * num
 			MultSchoolBook(num, 9, num, 9, prod_temp);
 			ConversionTo29bit(prod_temp, 16);
 			// num memory is set to 0, for safety o/w some cell may not be overwrite, keeping previous content
 			InitializeToZero(num);
 			barrett_reduction(prod_temp, 18, prime, MAX_SIZE, num, TString);
 			
 			// prod_temp memory is set to 0
 			InitializeToZero(prod_temp);
 		}
		else
		{
			// doing num * R
			MultSchoolBook(num, 9, R, 9, prod_temp);
 			ConversionTo29bit(prod_temp, 16);
 			// num memory is set to 0, for safety o/w some cell may not be overwrite, keeping previous content
 			InitializeToZero(num);
 			barrett_reduction(prod_temp, 18, prime, MAX_SIZE, num, TString);
 			
 			// prod_temp memory is set to 0
 			InitializeToZero(prod_temp);
 			
 			// doing R * R
 			MultSchoolBook(R, 9, R, 9, prod_temp);
 			ConversionTo29bit(prod_temp, 16);
 			// R memory is set to 0, for safety o/w some cell may not be overwrite, keeping previous content
 			InitializeToZero(R);
 			barrett_reduction(prod_temp, 18, prime, MAX_SIZE, R, TString);
		
			// prod_temp memory is set to 0
 			InitializeToZero(prod_temp);
		}
	}
}


// Montogomery without branching if-else condition
// output : num <-- num ^ exp
void Montogomery_Exponentiation_WO_Branching(long long *prime, long long *TString, long long *num, char *exp)
{
	long long prod_temp[MAX_SIZE] = {0};  // temporary memory
	long long R[MAX_SIZE] = {0};	// initialization
	long long Sq_num[MAX_SIZE] = {0};
	long long Sq_R[MAX_SIZE] = {0};
	long long *ptr;
	
	// R = num * num
	MultSchoolBook(num, 9, num, 9, prod_temp);
 	ConversionTo29bit(prod_temp, 16);
 	barrett_reduction(prod_temp, 18, prime, MAX_SIZE, R, TString);
 		
	// prod_temp memory is set to 0
 	InitializeToZero(prod_temp);
	
	// from 2nd msb to lsb, exp[0] is msb
	for(int i = 1; i < 256; i++)
	{
		// doing num * R
 		MultSchoolBook(num, 9, R, 9, prod_temp);
 		ConversionTo29bit(prod_temp, 16);
 		//First select the array and initialize it with 0, then do the barrett reduction
 		barrett_reduction(prod_temp, 18, prime, MAX_SIZE, SelectArrayWithInit(R, num, exp[i] - '0'), TString);
		InitializeToZero(prod_temp);
		
		// doing Sq_num = num * num
 		MultSchoolBook(num, 9, num, 9, Sq_num);
 		
 		// doing Sq_R = R * R
 		MultSchoolBook(R, 9, R, 9, Sq_R);
 		
 		//selecting array according to exp[i]
		ptr = SelectArray(Sq_num, Sq_R, exp[i] - '0');
		
		// copying content of ptr to prod_temp
		for(int j = 0; j < MAX_SIZE; j++)
			prod_temp[j] = ptr[j];
	
		ConversionTo29bit(prod_temp, 16);
		barrett_reduction(prod_temp, 18, prime, MAX_SIZE, SelectArrayWithInit(num, R, exp[i] - '0'), TString);
		
		// initializing the arrays to avoid memory replace error
		InitializeToZero(prod_temp);
		InitializeToZero(Sq_num);
		InitializeToZero(Sq_R);
	}
}



int main()
{
    // to compute prod = num ^ exp, num is a generator
 	long long num[MAX_SIZE] = {2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

 	// final result will be stored here, initialized to 1
 	long long prod[MAX_SIZE] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    long long TString[MAX_SIZE] = 
 	{450887704,490307913,387807083,403879883,291135210,307268612,110539282,24605042,70628772,35,0,0,0,0,0,0,0,0,0,0};
 	long long prime[MAX_SIZE] = 
 	{535425013,174332635,444665496,192778653,388389189,518147849,304619691,363717891,15281728,0,0,0,0,0,0,0,0,0,0,0};

////////////////////////////////////////// L2R & R2L Multiplication ////////////////////////
    char expString[257] = {0};
 	long long exp[9] = {0};
 	
	// "exp_256.txt" stores the 256-bit exponent
 	FILE *file5 = fopen("exp_256.txt", "r");
    // Read the bit string from the file
    if (!file5) {
        perror("Unable to open file");
        return 1;
    }
 	
    // parsing the string into 29-bit string
 	ParsingArray(exp, expString, file5, 256);
 	fclose(file5);

 	// Left to Right exponentiation
 	L2R_Exponentiation(prod, prime, TString, num, exp);
 	printf("after left to right Expo: \n");
 	printArray(prod, MAX_SIZE);
 	
 	// Right to Left exponentiation
 	//R2L_Exponentiation(prod, prime, TString, num, exp);
 	//printf("after right to left Expo: \n");
 	//printArray(prod, MAX_SIZE);


//////////////////////////////// Montgomery Multiplication ///////////////////////
/* 	char exp[257] = {0};
 	
 	FILE *file5 = fopen("exp_256.txt", "r");
    // Read the bit string from the file
    if (!file5) {
        perror("Unable to open file");
        return 1;
    }
 	
 	if (fgets(exp, 257, file5) == NULL || strlen(exp) != 256) {
        fprintf(stderr, "Invalid bit string format in file.\n");
        fclose(file5);
        return 0;
    }

 	// Montogomery multiplication
 	//Montogomery_Exponentiation(prime, TString, num, exp);
 	//printf("Montogomery: \n");
 	//printArray(num, MAX_SIZE);
 	
 	//Montogomery_Exponentiation_WO_Branching(prime, TString, num, exp);
 	//printf("Montogomery without Branching: \n");
 	//printArray(num, MAX_SIZE);
*/
    return 0;
}