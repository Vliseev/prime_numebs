/**************************************************************
 *
 *	tutor2.c
 *
   Example program for Montgomery-format elliptic functions found in ellmont.[ch].
   The primary feature of Montgomery format is that x-coordinates alone are
   calculated, so the elliptic arithmetic is efficient.  A generalized
   Montgomery scenario arises from:

		y^2 = x^3 + c x^2 + a x + b,

   which is an overspecified  parameterization
   (see comments in ellmont.c) but is neverthless
   useful in actual implementation.

   Furthermore, using what we call Montgomery-Brent parameterization, 
   one can find a random curve parameter pair (An/Ad), and also
   initial point {x,z} on the curve

     y^2 = x^3 + (An/Ad-2) x^2 + x  (mod p)

   such that the curve order is guaranteed to be divisible by 12 (a property
   deemed advantageous for factoring problems), and no explicit inversions of any
   kind need occur during an entire factoring trial.

   Note that for standard cryptographic work, the projective format of
   ellproj.c is preferred, for the Montgomery format ignores y coordinates.
   Thus, though Montgomery format has certain applications to cryptography,
   it is indicated most strongly in the related but different
   domains of factoring and primality proving.

   Compile with

   % cc -O tutor2.c ellmont.c giants.c -lm -o tutor2

   and run as:

   % tutor2

   whence a curve is chosen, a point P = {X,Z} on said curve established, and
   a few multiples of P are produced.  Note that these should agree (in their
   x-coordinates, after normalization) with the results from tutor.c

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
#include "ellmont.h"

static giant t1 = NULL;
point_mont pt1, pt2;

void
ensure(int sh)
{
	if(!t1) t1 = newgiant(sh);
	pt1 = new_point_mont(sh);
	pt2 = new_point_mont(sh);
}

void
ma_out_mont(point_mont pt) 
/* Mathematica list format, output. */
{
	printf("{");
	gout(pt->x); printf(",");
	gout(pt->z); printf("}\n");
}

/**************************************************************
 *
 *	Main Function
 *
 **************************************************************/
#define WORDS 100

main(
	int	argc,
	char 	*argv[]
)
{   giant a = newgiant(WORDS), b = newgiant(WORDS), c = newgiant(WORDS),
		  p = newgiant(WORDS);
	unsigned int j, k;

	ensure(WORDS);
	itog(1, p); gshiftleft(127, p); itog(1, a); subg(a, p); /* p = 2^127-1. */
	itog(3, a);
	itog(17, b);  
	itog(0, c);  /* The elliptic curve will be y^2 = x^3 + a x + b. */
    itog(3, t1);  /* Start with seed = 3 for random-point search. */
	find_point_mont(pt1,t1,a,b,c,p); /* {x, z = 1} is now a normalized point. */
	printf("Starting point P:\n");
    ma_out_mont(pt1);
    printf("The first few multiples of P are:\n");
	for(j=0; j <= 4; j++) {
		itog(j, t1);
        ell_mul_mont(pt1, pt2, t1, a, b, c, p);
		printf("%d * P = \n",j);
		ma_out_mont(pt2);
		normalize_mont(pt2,p);
		printf("==\n");
		ma_out_mont(pt2);
	}

	printf("Give a multiplier k:\n"); fflush(stdout);
	gin(t1);
    ell_mul_mont(pt1, pt2, t1, a, b, c, p);
	printf("k * P = \n");
    ma_out_mont(pt2);
	normalize_mont(pt2,p);
	printf("==\n");
	ma_out_mont(pt2);

}

