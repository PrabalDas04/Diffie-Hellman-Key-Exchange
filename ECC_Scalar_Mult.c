#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "Barrett_Reduction.h"

///////////////////////// function declaration ////////////////////////////
void L2R_Exponentiation(long long *prod, long long *prime, long long *TString, long long *num, long long *exp);
void elliptic_curve_add(long long *x1, long long *y1, long long *x2, long long *y2, long long *x3, long long *y3);
void L2R_Elliptic_Multiple(long long *x3, long long *y3, long long *x4, long long *y4, long long *exp);
void Inverse_Using_Fermat(long long *x, long long *x_inv);
int is_zero(long long *array);



// function for computing inverse of x by Fermat's theorem
void Inverse_Using_Fermat(long long *x, long long *x_inv)
{
	long long exp_for_inv[MAX_SIZE] = 
 	{535425011,174332635,444665496,192778653,388389189,518147849,304619691,363717891,15281728,0,0,0,0,0,0,0,0,0,0,0};
 	// final result will be stored here, initialized to 1
 	long long prod[MAX_SIZE] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
 	// T-value
 	long long TString[MAX_SIZE] = 
 	{450887704,490307913,387807083,403879883,291135210,307268612,110539282,24605042,70628772,35,0,0,0,0,0,0,0,0,0,0};
 	// prime value
 	long long prime[MAX_SIZE] = 
 	{535425013,174332635,444665496,192778653,388389189,518147849,304619691,363717891,15281728,0,0,0,0,0,0,0,0,0,0,0};
 	
 	// prod = x ^ exp_for_inv i.e x ^ (p-2)
	L2R_Exponentiation(prod, prime, TString, x, exp_for_inv);
	
	//copying the result in x_inv
	for(int i = 0; i < MAX_SIZE; i++)
		x_inv[i] = prod[i];
}

// function for checking array is 0 or not
// output 1 if array = 0 o/w 1
int is_zero(long long *array)
{
    for (int i = 0; i < MAX_SIZE; i++) 
    {
        if (array[i] != 0) 
            return 0; // Not zero
    }
    return 1; // All elements are zero
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



// Function for elliptic curve point addition
// this function also includes point doubling
// points: (x1, y1), (x2, y2); output: (x3, y3) <-- (x1, y1) + (x2, y2)
void elliptic_curve_add(long long *x1, long long *y1, long long *x2, long long *y2, long long *x3, long long *y3)
{
	//long long a[MAX_SIZE] =
	//{535425004,174332635,444665496,192778653,388389189,518147849,304619691,363717891,15281728,0,0,0,0,0,0,0,0,0,0,0};
	long long a[MAX_SIZE] = {0}; 
    long long TString[MAX_SIZE] = 
 	{450887704,490307913,387807083,403879883,291135210,307268612,110539282,24605042,70628772,35,0,0,0,0,0,0,0,0,0,0};
 	
 	long long prime[MAX_SIZE] = 
 	{535425013,174332635,444665496,192778653,388389189,518147849,304619691,363717891,15281728,0,0,0,0,0,0,0,0,0,0,0};
    
    long long m[MAX_SIZE] = {0};   // Slope
    // Temporary variables
    long long temp1[MAX_SIZE] = {0};      
    long long temp2[MAX_SIZE] = {0};
    long long temp3[MAX_SIZE] = {0};
    // array rep. of 3 and 2
    long long array_three[MAX_SIZE] = {0};
    array_three[0] = 3;
    long long array_two[MAX_SIZE] = {0};
    array_two[0] = 2;

    // If P is the point at infinity, return Q
    if (is_zero(x1) && is_zero(y1))
    {
        for (int i = 0; i < MAX_SIZE; i++)
        {
            x3[i] = x2[i];
            y3[i] = y2[i];
        }
        return;
    }

    // If Q is the point at infinity, return P
    if (is_zero(x2) && is_zero(y2))
    {
        for (int i = 0; i < MAX_SIZE; i++)
        {
            x3[i] = x1[i];
            y3[i] = y1[i];
        }
        return;
    }

    // If P == -Q, return point at infinity
    Addition_256(y1, y2, temp1); // temp1 = y1 + y2
    if (memcmp(x1, x2, MAX_SIZE * sizeof(long long)) == 0 && memcmp(temp1, temp2, MAX_SIZE * sizeof(long long)) == 0)
    {
        for (int i = 0; i < MAX_SIZE; i++)
        {
            x3[i] = 0;
            y3[i] = 0;
        }
        return;
    }
	
	// if x1 = x2 and y1 = y2
    if (memcmp(x1, x2, MAX_SIZE * sizeof(long long)) == 0 && memcmp(y1, y2, MAX_SIZE * sizeof(long long)) == 0)
    {
        // Point doubling case
        MultSchoolBook(x1, MAX_SIZE, x1, MAX_SIZE, temp2); 
        ConversionTo29bit(temp2, 16);
 		barrett_reduction(temp2, 18, prime, MAX_SIZE, temp1, TString);  // temp1 = x1^2
        InitializeToZero(temp2);
        
        MultSchoolBook(array_three, MAX_SIZE, temp1, MAX_SIZE, temp2);  
        ConversionTo29bit(temp2, 8);
        InitializeToZero(temp1);
 		barrett_reduction(temp2, 10, prime, MAX_SIZE, temp1, TString);  // temp1 = 3 * x1^2
        InitializeToZero(temp2);
        
        Addition_256(temp1, a, temp1);      // temp1 = 3 * x1^2 + a
       
        MultSchoolBook(array_two, MAX_SIZE, y1, MAX_SIZE, temp2);  
        ConversionTo29bit(temp2, 8);
 		barrett_reduction(temp2, 18, prime, MAX_SIZE, temp3, TString);  // temp3 = 2 * y1
        InitializeToZero(temp2);   
        
        Inverse_Using_Fermat(temp3, temp2);   // temp2 = (2 * y1)^(-1)
        InitializeToZero(temp3);
        
        MultSchoolBook(temp1, MAX_SIZE, temp2, MAX_SIZE, temp3);  
        ConversionTo29bit(temp3, 16);
 		barrett_reduction(temp3, 18, prime, MAX_SIZE, m, TString);  // m = (3 * x1^2 + a) / (2 * y1)
    	
 		// initialization to 0
        InitializeToZero(temp1); 
        InitializeToZero(temp2);
        InitializeToZero(temp3);
    } 
    else
    {
        // memory copy: temp2 <-- y1 and temp1 <--y2
        memcpy(temp2, y1, MAX_SIZE * sizeof(long long));
 		memcpy(temp1, y2, MAX_SIZE * sizeof(long long));
        
        subtract_any(temp1, temp2);        // temp1 = y2 - y1
        InitializeToZero(temp2);
        
        // memory copy: temp3 <-- x1 and temp2 <--x2
        memcpy(temp3, x1, MAX_SIZE * sizeof(long long));
        memcpy(temp2, x2, MAX_SIZE * sizeof(long long));
        
        subtract_any(temp2, temp3);        // temp2 = x2 - x1
        InitializeToZero(temp3);
        Inverse_Using_Fermat(temp2, temp3);  // temp3 = (x2 - x1)^(-1)
        InitializeToZero(temp2);
        
        MultSchoolBook(temp1, MAX_SIZE, temp3, MAX_SIZE, temp2);  
        ConversionTo29bit(temp2, 16);
 		barrett_reduction(temp2, 18, prime, MAX_SIZE, m, TString);  // m = (y2 - y1) / (x2 - x1)
        
        // initialization to 0
        InitializeToZero(temp1); 
        InitializeToZero(temp2);
        InitializeToZero(temp3);
    }
	// temp2 = m * m
    MultSchoolBook(m, MAX_SIZE, m, MAX_SIZE, temp2);   
    ConversionTo29bit(temp2, 16);
 	barrett_reduction(temp2, 18, prime, MAX_SIZE, temp1, TString);  // temp1 = m^2
 	InitializeToZero(temp2);
 	
 	// copying the data of x1 and x2 to tem2, temp3 resp
 	// this is needed as after subtraction both the function argument value will be changed
 	memcpy(temp2, x1, MAX_SIZE * sizeof(long long));
 	memcpy(temp3, x2, MAX_SIZE * sizeof(long long));
    subtract_any(temp1, temp2);	// temp1 = m^2 - x1
    subtract_any(temp1, temp3);	// temp1 = m^2 - x1 - x2
            
    memcpy(x3, temp1, MAX_SIZE * sizeof(long long)); // x3 = temp1 mod p
	
	// memory copy: temp2 <-- x1 and temp3 <--x2
	memcpy(temp2, x1, MAX_SIZE * sizeof(long long));
 	memcpy(temp3, x2, MAX_SIZE * sizeof(long long));
	
	subtract_any(temp2, temp1);		// temp2 = x1 - x3
	InitializeToZero(temp1);
	
    MultSchoolBook(m, MAX_SIZE, temp2, MAX_SIZE, temp1);   
    ConversionTo29bit(temp1, 16);
    InitializeToZero(temp2);
 	barrett_reduction(temp1, 18, prime, MAX_SIZE, temp2, TString);	// temp2 = m * (x1 - x3)
    
    // memory copy: temp1 <-- y1
    memcpy(temp1, y1, MAX_SIZE * sizeof(long long));
    subtract_any(temp2, temp1);		// temp2 = m * (x1 - x3) - y1
    
    // memory copy: y3 <-- temp2  
    memcpy(y3, temp2, MAX_SIZE * sizeof(long long));
}

// elliptic curve point multiplication using left-to-right method
// output : (x4, y4) <-- exp * (x3, y3)
void L2R_Elliptic_Multiple(long long *x3, long long *y3, long long *x4, long long *y4, long long *exp)
{
	int bin[29] = {0};	// array for storing binary rep of integer
	int j = 8, count = 0, k = 28;
	
	// temporary variable
	long long x5[MAX_SIZE] = {0};
    long long y5[MAX_SIZE] = {0}; 
	
	// this loop will iterate 261 times(29*9), since during int to binary conversion each number is represented in 29 bit string
	for(int i = 0; i < 261; i++)
 	{
 		// doing point doubling of p4 = (x4, y4) and stored in (x5, y5)
 		elliptic_curve_add(x4, y4, x4, y4, x5, y5);
 		// (x4, y4) <-- (x5, y5)
 		for(int i = 0; i < MAX_SIZE; i++)
 		{
 			x4[i] = x5[i]; 
 			y4[i] = y5[i];
 		}
 		// initialization to 0
 		InitializeToZero(x5);
 		InitializeToZero(y5);
 		
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
 		if(bin[k] == 1)
 		{
 			// doing (x5,y5) = (x4, y4) + (x3, y3)
 			elliptic_curve_add(x4, y4, x3, y3, x5, y5);
 			// (x4, y4) <-- (x5, y5)
 			for(int i = 0; i < MAX_SIZE; i++)
 			{
 				x4[i] = x5[i]; 
 				y4[i] = y5[i];
 			}
 			// initialization to 0
 			InitializeToZero(x5);
 			InitializeToZero(y5);
 		}
 		count--;
 		k--;
 	}
}


int main()
{
    long long x1[MAX_SIZE] = {1, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    long long y1[MAX_SIZE] = {1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    long long x2[MAX_SIZE] = {1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    long long y2[MAX_SIZE] = {1, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    long long x3[MAX_SIZE] = {0};
    long long y3[MAX_SIZE] = {0};
    
    elliptic_curve_add(x1, y1, x2, y2, x3, y3);
    
    printf("x3 :\n");
    printArray(x3, MAX_SIZE);
    printf("y3 :\n");
    printArray(y3, MAX_SIZE);
    
    long long x4[MAX_SIZE] = {0};
    long long y4[MAX_SIZE] = {0};
   
    long long exp[MAX_SIZE] = {1005178, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    L2R_Elliptic_Multiple(x3, y3, x4, y4, exp);
    
    printf("x4 :\n");
    printArray(x4, MAX_SIZE);
    printf("y4 :\n");
    printArray(y4, MAX_SIZE);
    
    return 0;

}


