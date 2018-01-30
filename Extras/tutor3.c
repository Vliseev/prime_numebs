/**************************************************************
 *
 *	tutor3.c
 *
   Example program for affine-format elliptic functions found 
   in ellaffi.[ch].

   Compile with

   % cc -O tutor3.c ellaffi.c tools.c giants.c -lm -o tutor3

   and run as:

   % tutor3

   whence a curve is chosen, a point P = {x,y} on said curve 
   established, and a few multiples of P are produced.  
   The results should agree with those of tutor.c (projective)
   and tutor2.c (Montgomery) scearios.

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
#include "ellaffi.h"
#include "tools.h"

static giant t1 = NULL;
point_affi pt1, pt2;


void
ensure(int sh)
{
	if(!t1) t1 = newgiant(sh);
	pt1 = new_point_affi(sh);
	pt2 = new_point_affi(sh);
}

void
ma_out_affi(point_affi pt) 
/* Mathematica list format, output. */
{
	if(pt->z == 0) {
		printf("POINT-AT-INFINITY\n");
		return;
	}
	printf("{");
	gout(pt->x); printf(",");
	gout(pt->y); 
    printf("}\n");
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
{   giant a = newgiant(WORDS), b = newgiant(WORDS), 
		  p = newgiant(WORDS);
	int j;

    init_tools(WORDS);  /* Initialize algorithms. */
	ensure(WORDS);
	itog(1, p); gshiftleft(127, p); itog(1, a); subg(a, p); /* p = 2^127-1. */
	itog(3, a);
	itog(17, b);  /* The elliptic curve will be y^2 = x^3 + a x + b. */
    itog(3, t1);  /* Start with seed = 3 for point search. */
	find_point_affi(pt1,t1,a,b,p); /* {x, y, z = (bool)1} is now a normalized point. */
	printf("Starting point P:\n");
    ma_out_affi(pt1);
    printf("The first few multiples of P are:\n");
	for(j=0; j <= 4; j++) {
		itog(j, t1);
        ell_mul_affi(pt1, pt2, t1, a, p);
		printf("%d * P = \n",j);
		ma_out_affi(pt2);
	}

	printf("Give a multiplier k:\n"); fflush(stdout);
	gin(t1);
    ell_mul_affi(pt1, pt2, t1, a, p);
	printf("k * P = \n");
    ma_out_affi(pt2);
}

