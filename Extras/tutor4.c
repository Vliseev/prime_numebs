/**************************************************************
 *
 *	tutor4.c
 *
   Example program for testing factoring module fmodule.[ch].

   Compile with

   % cc -O tutor4.c fmodule.c tools.c giants.c ellmont.c -lm -o tutor4

   and run as:

   % tutor4

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

int			
psi_rand(
	void
)
{
	unsigned short	hi;
	unsigned short	low;
	time_t			tp;
	int				result;

	time(&tp);
	low = (unsigned short)rand();
	hi = (unsigned short)rand();
	result = ((hi << 16) | low) ^ ((int)tp);

	return (result & 0x7fffffff);
}

main(
	int	argc,
	char 	*argv[]
)
{   giant N = newgiant(WORDS), x = newgiant(WORDS);
    int j, ct;
	unsigned int curve, B, C, S;

	init_fmodule(WORDS);  /* Start up factoring engine. */
	verbose(1);
	gin(N);

printf("Commencing sieve...\n"); fflush(stdout);

	ct = sieve(N, 65536);
	for(j=0; j<ct; j++) {
		printf("%d^%d\n", (int)(prime_list()[j]), (int)(exponent_list()[j]));
		fflush(stdout);
	}
	
	if(isone(N)) exit(0);  /* Exit if sieve saturates N. */
	if(pseudointq(5,N)) {  /* Exit if N is pseudo, base 5. */
				gout(N);
				exit(0);
	}
printf("Cofactor is composite...commencing Pollard rho method...\n");
fflush(stdout);
	if(pollard_rho(N, x, 15000, 1)) {
		gout(x);
	    if(pseudointq(5,N)) { /* Exit if N is pseudo, base 5. */
				gout(N);
				exit(0);
		}
	}
printf("Cofactor is composite...commencing Pollard-(p-1) method...\n");
fflush(stdout);

	if(pollard_pminus1(N, x, 10000, 1)) {
		gout(x);
	    if(pseudointq(5,N)) { /* Exit if N is pseudo, base 5. */
				gout(N);
				exit(0);
		}
	}

printf("Cofactor is composite...commencing elliptic curve method...\n");
fflush(stdout);

	curve = 0;
	while(++curve) {
    	if(curve < 4) {
			B = 1000;
		} else if(curve < 10) {
			B = 10000;
		} else if(curve < 33) {
			B = 100000;
		} else B = 1000000;
    	C = 50*B;
    	S = psi_rand();
    	printf("Curve %u: ", curve);
		if(ecm(N, x, S, B, C)) {
		gout(x);
	    if(pseudointq(5,N)) { /* Exit if N is pseudo, base 5. */
				gout(N);
				exit(0);
		}
		printf("\n"); fflush(stdout);
	}
}
}

