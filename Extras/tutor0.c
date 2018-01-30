/**************************************************************
 *
 *	tutor0.c
 *
   Example program for projective coordinate elliptic algebra.as
   in ellproj.[ch].
   The elliptic curve is over F_p with p = 2^127-1,
   
   y^2 = x^3 + a x + b

   with a = 3, b = 17 fixed for the example. 
   The program simply generates random (X, Y, 1) on the curve.

   Compile with

   % cc -O tutor0.c ellproj.c tools.c giants.c -lm -o tutor0

   and run as:

   % tutor0

   and points will be generated.

 *	Updates:
 *    8 Jan 02    REC - Klinger's WIN32 random seed.
 *		3 Apr 98    REC - Creation
 *
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

#define WORDS 100

void
ma_out_proj(point_proj pt) 
/* Mathematica list format, output. */
{
	printf("{");
	gout(pt->x); printf(",");
	gout(pt->y); printf(",");
	gout(pt->z); printf("}\n");
}

main(
	int	argc,
	char 	*argv[]
)
{   giant a = newgiant(WORDS), b = newgiant(WORDS), 
		  p = newgiant(WORDS), x = newgiant(WORDS);
    point_proj pt = new_point_proj(WORDS);
	int j, seed;

	init_tools(WORDS);  /* Initialize algorithms. */
	itog(1, p); gshiftleft(192, p);
    itog(399, a); subg(a,p); /* Field is F_p, p = 2^192-399. */
	itog(3, a);
	itog(17, b);  /* Curve will be y^2 = x^3 + a x + b (mod p). */
    printf("p = "); gout(p);
    printf("a = "); gout(a);
    printf("b = "); gout(b);
    printf("Now we generate random points on the curve:\n");
    for(j=0; j<10; j++) {
  
      #ifdef _WIN32
        srand((unsigned)time(NULL));
        itog(rand(), x);  /* Seed the point search. */
      #endif
      #ifndef _WIN32
        itog(random(), x); /* Seed the point search. */
      #endif
      find_point_proj(pt,x,a,b,p); 
			/* {x, y, z = 1} is now a normalized point. */
    	ma_out_proj(pt);
	}
}

