/**************************************************************
 *
 *	tutor1.c
 *
   Example program for projective coordinate elliptic algebra.as
   in ellproj.[ch].
   The elliptic curve is over F_p with p = 2^127-1,
   
   y^2 = x^3 + a x + b

   with a = 3, b = 17 fixed for the example. 
   An initial point P is generated, using find_point_proj().  Then
   a few multiples of P are reported, followed by an option to input
   an arbitrary multiple of P.

   Compile with

   % cc -O tutor1.c ellproj.c tools.c giants.c -lm -o tutor1

   and run as:

   % tutor1

   with an eventual prompt to type (or paste) in an actual 
   multiplier k, whence a report of the elliptic multiple 
   k*P emanates.
   
   The curve with a = 3, b = 17 has order

   |E| = 170141183460469231731611440180631708799
       = 3 * 7 * 7 * 13 * 28027 * 
		 3176670344634393565948814741467

   and the point P reported has the big prime order 317667....
   All of this can be verified in run-time fashion;
   i.e. one may verify that

     |E| * P = {1,1,0} = point at infinity

   but that

     (3 * 7 * 7 * 13 * 28027) * P = 53559597 * P   !=   {1,1,0}

   so that P has the big prime order stated.	

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
#include "ellproj.h"
#include "tools.h"

static giant t1 = NULL;
point_proj pt1, pt2;


void
ensure(int sh)
{
	if(!t1) t1 = newgiant(sh);
	pt1 = new_point_proj(sh);
	pt2 = new_point_proj(sh);
}

void
ma_out_proj(point_proj pt) 
/* Mathematica list format, output. */
{
	printf("{");
	gout(pt->x); printf(",");
	gout(pt->y); printf(",");
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
{   giant a = newgiant(WORDS), b = newgiant(WORDS), p = newgiant(WORDS);
	int j;

    init_tools(WORDS);  /* Initialize algorithms. */
	ensure(WORDS);
	itog(1, p); gshiftleft(127, p); itog(1, a); subg(a, p); /* p = 2^127-1. */
	itog(3, a);
	itog(17, b);  /* The elliptic curve will be y^2 = x^3 + a x + b. */
    itog(3, t1);  /* Start with seed = 3 for random-point search. */
	find_point_proj(pt1,t1,a,b,p); /* {x, y, z = 1} is now a normalized point. */
	printf("Starting point P:\n");
    ma_out_proj(pt1);
    printf("The first few multiples of P are:\n");
	for(j=0; j <= 4; j++) {
		itog(j, t1);
        ell_mul_proj(pt1, pt2, t1, a, p);
		printf("%d * P = \n",j);
		ma_out_proj(pt2);
		normalize_proj(pt2,p);
		printf("==\n");
		ma_out_proj(pt2);
	}

	printf("Give a multiplier k:\n"); fflush(stdout);
	gin(t1);
    ell_mul_proj(pt1, pt2, t1, a, p);
	printf("k * P = \n");
    ma_out_proj(pt2);
	normalize_proj(pt2,p);
	printf("==\n");
	ma_out_proj(pt2);

}

