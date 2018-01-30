/**************************************************************
 *
 * ellaffi.c
 *
   Fast algorithms for fundamental elliptic curve arithmetic,
   standard (affine) format.  Such algorithms apply in domains 
   such as:
    -- factoring
    -- primality studies (e.g. rigorous primality proofs)
    -- elliptic curve cryptography (ECC) 
  
   AFFINE FORMAT

   Functions are supplied herein for affine format
   of points.  Alternative formats differ in their
   range of applicability, efficiency, and so on.
   Primary advantages of the standard, affine format herein are:
    -- Simply operation structure, relative to the
       alternative formats (as in ellproj.[ch], ellmont.[ch])
    -- Low operation count: for general addition per se, 
	   2 muls, 1 square, and 1 inverse in the field 
	   (except note that the typical
       inefficiency of the inverse for odd characteristic
       fields is the main motivation for alternative
       formats)
	-- Low memory requirement, compared to, say, the
       projective format of ellproj.[ch]

   The elliptic curve is over F_p, with p > 3 prime, and Weierstrass
   parameterization:

      y^2 = x^3 + a x + b

   The affine-format coordinates are actually stored in
   the form {x, y, z}, with the z-component signifying, upon
   vanishing, the point-at-infinity.  (And thus z need not be
   kept in multiprecision fashion -- we keep it as a Bollean int.) 
   
   The basic point-multiplication function is

      ell_mul_affi()

   which obtains the result k * P for given point P and integer
   multiplier k.  
 
   REFERENCES

   Crandall R and Pomerance C 1998, "Prime numbers: a computational
		perspective," Springer-Verlag, manuscript

   Solinas J 1998, IEEE P1363 Annex A (draft standard)

   LEGAL AND PATENT NOTICE

   This and related PSI library source code is intended solely for 
   educational and research applications, and should not be used
   for commercial purposes without explicit permission from PSI
   (not to mention proper clearance of legal restrictions beyond
   the purview of PSI).  
   The code herein will not perform cryptography per se,
   although the number-theoretical algorithms herein -- all of which 
   are in the public domain -- can be used in principle to effect 
   what is known as elliptic curve cryptography (ECC).  Note that 
   there are strict constraints on how cryptography may be used, 
   especially in regard to exportability.
   Therefore one should avoid any casual distribution of actual 
   cryptographic software.  Along similar lines, because of various 
   patents, proprietary to Apple Computer, Inc., and perhaps to other 
   organizations, one should not tacitly assume that an ECC scheme is 
   unconstrained.  For example,the commercial use of certain fields 
   F_p^k (i.e., fixation of certain primes p) is covered in Apple 
   patents.

 *	Updates:
 *    8 Jan 02    REC - rearranged functions in .[ch]
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
#include "ellaffi.h"
#include "tools.h"


/* global variables */

static giant t0 = NULL, t1 = NULL, t2 = NULL, t3 = NULL;

/**************************************************************
 *
 *	Maintenance functions
 *
 **************************************************************/

point_affi
new_point_affi(int shorts)
{
	point_affi pt;

	if(t0 == NULL) init_ell_affi(shorts);
	pt = (point_affi) malloc(sizeof(point_struct_affi));
	pt->x = newgiant(shorts);
	pt->y = newgiant(shorts);
	pt->z = 1;
	return(pt);
}

void
free_point_affi(point_affi pt)
{
	free(pt->x); free(pt->y);
	free(pt);
}

void
ptop_affi(point_affi pt1, point_affi pt2)
{
	gtog(pt1->x, pt2->x);
	gtog(pt1->y, pt2->y);
}

int
init_ell_affi(int shorts) 
/* Called by new_point_affi(), to set up giant registers. */
{	
	t0 = newgiant(shorts);
	t1 = newgiant(shorts);
	t2 = newgiant(shorts);
	t3 = newgiant(shorts);
}


/**************************************************************
 *
 *	Elliptic curve operations
 *
 **************************************************************/

/* Begin affine-format functions for 
 
   y^2 = x^3 + a x + b.

   A point is kept as a triple {x, y, z}, with z a Boolean
   int and point-at-infinity = {1,1,0}.

 */

void
ell_double_affi(point_affi pt, giant a, giant p)
/* pt := 2 pt on the curve. */
{	
	giant x = pt->x, y = pt->y;

	if(isZero(y) || (pt->z == 0)) {
		itog(1,x); itog(1,y); pt->z = 0;
		return;
	}	
	gtog(x, t1); squareg(t1); modg(p, t1);
	smulg(3, t1);
	addg(a, t1); /* t1 := 3x^2 + a. */
	gtog(y, t2); gshiftleft(1, t2); 
	binvg(p, t2); /* t1 := 1/(2y). */
	mulg(t2, t1); modg(p, t1);  /* t1 = slope m. */
	gtog(t1, t2); squareg(t2); modg(p, t2);
	gtog(x, t3);
	subg(x, t2); subg(x, t2); gtog(t2, x); modg(p, x);
		/* x := m^2 - 2x. */
	subg(x, t3); mulg(t1, t3); subg(y, t3);
	gtog(t3, y); modg(p, y);
}
 
void
ell_add_affi(point_affi pt0, point_affi pt1, giant a, giant p)
/* pt0 := pt0 - pt1 on the curve. */
{   
	giant x0 = pt0->x, y0 = pt0->y,
		  x1 = pt1->x, y1 = pt1->y;
 
	if(pt0->z == 0) {
		gtog(x1,x0); gtog(y1,y0); pt0->z = pt1->z;
		return;
	}
	if(pt1->z == 0) return;
	if(gcompg(x0, x1) == 0) {
		if(gcompg(y0, y1) == 0) {
			ell_double_affi(pt0, a, p);
			return;
		}
		itog(1,x0); itog(1,y0); pt0->z = 0;
	}
	gtog(x1, t2); subg(x0, t2); modg(p, t2);
	binvg(p, t2);  /* t2 := 1/(x1-x0). */
	gtog(y1, t1); subg(y0, t1); mulg(t2, t1);
	modg(p, t1); /* t1 := slope m. */
	gtog(t1, t2); squareg(t2); modg(p, t2);
	subg(x0, t2); subg(x1, t2); modg(p, t2);
		/* t2 := m^2 - x0 - x1. */
	gtog(x1, t3); subg(t2, t3); mulg(t1, t3);
	subg(y1, t3); modg(p, t3);
	gtog(t2, x0);
	gtog(t3, y0); 
}


void
ell_neg_affi(point_affi pt, giant p)
/* pt := -pt on the curve. */
{
	negg(pt->y); modg(p, pt->y);
}

void
ell_sub_affi(point_affi pt0, point_affi pt1, giant a, giant p)
/* pt0 := pt0 - pt1 on the curve. */
{
	ell_neg_affi(pt1, p);
	ell_add_affi(pt0, pt1, a, p);
	ell_neg_affi(pt1, p);
}

void
ell_mul_affi(point_affi pt0, point_affi pt1, giant k, giant a, giant p)
/* General elliptic multiplication;
   pt1 := k*pt0 on the curve, 
   with k an arbitrary integer. 
 */
{	
	giant x = pt0->x, y = pt0->y,
		  xx = pt1->x, yy = pt1->y;
	int ksign, hlen, klen, b, hb, kb;
    
	if(isZero(k)) {
		itog(1, xx);
		itog(1, yy);
		pt1->z = 0;
		return;
	}
    ksign = k->sign;
	if(ksign < 0) negg(k);
	gtog(x,xx); gtog(y,yy); pt1->z = pt0->z;
	gtog(k, t0); addg(t0, t0); addg(k, t0); /* t0 := 3k. */
	hlen = bitlen(t0);
	klen = bitlen(k);
	for(b = hlen-2; b > 0; b--) {
		ell_double_affi(pt1,a,p);
		hb = bitval(t0, b);
		if(b < klen) kb = bitval(k, b); else kb = 0;
		if((hb != 0) && (kb == 0))
			ell_add_affi(pt1, pt0, a, p);
		else if((hb == 0) && (kb !=0))
			ell_sub_affi(pt1, pt0, a, p);
	}
	if(ksign < 0) {
		ell_neg_affi(pt1, p);
		k->sign = -k->sign;
	}
}


void
find_point_affi(point_affi pt, giant seed, giant a, giant b, giant p)
/* Starting with seed, finds a point {x,y, (int)1} on curve.  */
{	giant x = pt->x, y = pt->y;

    modg(p, seed);
	while(1) {
		gtog(seed, x);
		squareg(x); modg(p, x);
		addg(a, x);
		mulg(seed,x); addg(b, x);
		modg(p, x); /* x := seed^3 + a seed + b. */
		if(sqrtmod(p, x)) break;  /* test cubic form for having root. */
		iaddg(1, seed);
	}
  gtog(x,y);
	gtog(seed,x);
	pt->z = 1;
}


