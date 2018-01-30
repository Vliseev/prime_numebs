/**************************************************************
 *
 * tutor5.c
 *
   Example program for trying to yank a considerable factor,
   very quickly, from an integer N.

   Compile with

   % cc -O tutor5.c fmodule.c tools.c giants.c ellmont.c -lm -o tutor5

   and run as:

   % tutor5

   whence you enter a number N to be factored.


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
#include "fmodule.h"

/**************************************************************
 *
 *	Main Function
 *
 **************************************************************/
#define WORDS 100   /* 16*WORDS being the maximum bit length
					   of an operand (prior to reduction (mod N)). */
main(
	int	argc,
	char 	*argv[]
)
{   giant N = newgiant(WORDS), x = newgiant(WORDS);
    int j, ct;
	unsigned int curve, B, C, S;

	init_fmodule(WORDS);  /* Start up factoring engine. */
	verbose(0);  /* Let's run silent. */
	gin(N);

	ct = sieve(N, 65536);
	for(j=0; j<ct; j++) {
		printf("%d^%d\n", (int)(prime_list()[j]), (int)(exponent_list()[j]));
		fflush(stdout);
	}
	
	if(isone(N)) exit(0);  /* Exit if sieve saturates N. */
	if(pseudointq(5,N)) goto PP;  /* Exit if N is pseudo, base 5. */

verbose(1);

	if(pollard_rho(N, x, 100000)) {
		gout(x);
	    if(pseudointq(5,N)) goto PP; 
								  /* Exit if N is pseudo, base 5. */
	}


/* Composite output. */
	printf("Composite:\n");
	gout(N);
	exit(0);

PP:  /* Probable prime output. */
	printf("Probable prime:\n");
	gout(N);
	exit(0);

}

