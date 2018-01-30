/**************************************************************
 *
 * tutor6.c
 *
   Complex-multiplication (CM) representation finder.
   This program accepts an input number (perhaps a probable prime)
   p and attempts to find representations:

       4p = u^2 + |D| v^2

   where D is a fundamental (negative) discriminant from
   disc.h.  A primary application of such representations is,
   that if we can factor p + 1 +- u (or one of a few other
   similar expressions when D = -3, -4) into f * r, where
   r is again a smaller probable prime, then there
   is good chance that elliptic curve primality proving (ECPP)
   can be invoked, whereby if r is indeed prime (and certain curve
   calculations go through) then p is prime. 

   Compile with:

   % cc -O tutor6.c tools.c giants.c -lm -o tutor6

   and run as:

   % tutor6

   and input the prime.


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
#include "disc.h"

/**************************************************************
 *
 *	Main Function
 *
 **************************************************************/
#define WORDS 400   /* 16*WORDS being the maximum bit length
					   of an operand (prior to modular reduction). */
main(
	int	argc,
	char 	*argv[]
)
{   giant p = newgiant(WORDS), u = newgiant(WORDS),
		  v = newgiant(WORDS);
    int j, d, test;

	init_tools(WORDS);
    gin(p);  /* input prime such as 618970019642690137449562111 */
	for(j = 0; j < 20; j++) {
		d = disc[j];
		printf("%d: ", d); fflush(stdout);
		test = cornacchia4(p, d, u, v);
		printf("%d\n", test);
    if(test) {
		   printf("u = "); gout(u);
		   printf("v = "); gout(v);
    }
	}
}

