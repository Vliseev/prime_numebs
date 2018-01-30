/**************************************************************
 *
 * tutor7.c
 *
   Probable-primes test example.
   The Fermat test and (more reliable) Miller-Rabin test
   are invoked for test of probable primality.

   Compile with:

   % cc -O tutor7.c tools.c giants.c -lm -o tutor7

   and run as:

   % tutor7

   whence various probable primes are reported, including
   some well-known, faulty reports for some of the tests.
   As a final option, you may manually enter a number for
   probable-primality testing.


 *	c. 1998 Perfectly Scientific, Inc.
 *	All Rights Reserved.
 *
 *
 *************************************************************/

/* include files */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#ifdef _WIN32 
#include <process.h>
#endif
#include <string.h>
#include "giants.h"
#include "tools.h"

/**************************************************************
 *
 *	Main Function
 *
 **************************************************************/
#define WORDS 4000   /* 16*WORDS being the maximum bit length
					   of an operand (prior to reduction (mod N)). */
main(
	int	argc,
	char 	*argv[]
)
{   giant p = newgiant(WORDS), x = newgiant(WORDS);
    int j, d, test;

	init_tools(WORDS);
    /* Let us first test p = F_7 = 2^128 + 1. */
	itog(1, p); gshiftleft(128, p); iaddg(1, p);
	printf("Analyzing character of ");
	gout(p); fflush(stdout);
	if(pseudointq(2, p)) {
		printf("Base-2 pseudoprime. \n");
		if(prime_probable(p)) {
			printf("Probable prime (MR test). \n");
		} else printf("Composite.\n");
	} else {
		printf("Composite.\n");
	}
	printf("\n");

    /* Next, let us test a peculiar number:
       p = 3125031751 = 151 * 751 * 28351.
     */
	itog(151, p); smulg(751, p); smulg(28351, p);
	printf("Analyzing character of ");
	gout(p); fflush(stdout);
	if(pseudointq(7, p)) {
		printf("Base-7 pseudoprime. \n");
		if(prime_probable(p)) {
			printf("Probable prime (MR test). \n");
		} else printf("Composite.\n");
	} else {
		printf("Composite.\n");
	}
	printf("\n");

    /* Next, let us test the repunit number:
       p = (10^23-1)/9, a known prime.
     */
	powerg(10, 23, p);
	itog(1, x); subg(x, p); idivg(9, p);
	printf("Analyzing character of ");
	gout(p); fflush(stdout);
	if(pseudointq(5, p)) {
		printf("Base-5 pseudoprime. \n");
		if(prime_probable(p)) {
			printf("Probable prime (MR test). \n");
		} else printf("Composite.\n");
	} else {
		printf("Composite.\n");
	}
	printf("\n");

    /* Finally, manual test. */
	printf("Enter a p for testing:\n"); fflush(stdout);
	gin(p);
	printf("Analyzing character of ");
	gout(p); fflush(stdout);
	if(pseudointq(3, p)) {
		printf("Base-3 pseudoprime. \n");
		if(prime_probable(p)) {
			printf("Probable prime (MR test). \n");
		} else printf("Composite.\n");
	} else {
		printf("Composite.\n");
	}	
}

