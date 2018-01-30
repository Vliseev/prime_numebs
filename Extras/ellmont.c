/**************************************************************
 *
 *	ellmont.c
 *
   Fast algorithms for fundamental elliptic curve arithmetic,
   projective format.  Such algorithms apply in domains such as:
    -- factoring
    -- primality studies (e.g. rigorous primality proofs)
    -- elliptic curve cryptography (ECC) 
  
   MONTGOMERY FORMAT

   Functions are supplied herein for Montgomery format
   of points.  Alternative formats differ in their
   range of applicability, efficiency, and so on.
   Primary advantages of the Montgomery format herein are:
    -- No explicit inversions (until perhaps one such at the end of
       an elliptic multiply operation)
    -- Furthermore, no inversions whatsoever for factoring
       studies
    -- Fairly low operation count (~6 muls for point doubling,
       ~8 muls for point addition)
    -- No reference is made to y-coordinates at any time (of
       course this could also be construed as a disadvantage)  

   Here the curve is again over F_p, or sometimes over F_N where N is
   a composite number to be factored (and so F_N is a 
   "pseudo-field").The Montgomery parameterization is taken to be
   the generalization:

      y^2 = x^3 + c x^2 + a x + b

   for some of the functions, while the Montgomery-Brent
   specialization: c = An/Ad-2, a = 1, b = 0 is assumed for
   yet other functions.  Note that the (a,b,c) parameterization
   is actually over-specified, since any such can be brought down
   to Weierstrass form y^2 = x^3 + a' x + b'.  However, using
   the (a,b,c) format allows for many special cases to be
   handled efficietly in one code set, as is done herein.  Note that
   the Montgomery-Brent curves

     	y^2 = x^3 + (An/Ad-2) x^2 + x

   are intended primarily for factoring.  Indeed, herein are
   functions for guaranteeing that the order of such a curve is
   divisible by 12.

   The Montgomery-format points are of the form {X, Z}, with 
   y-coordinates ignored, the true x-coordinate being x = X/Z.
   The point multiplication functions are:

      	ell_mul_mont()
	  	ell_mul_brent()

   which obtain the {X/Z} pair for a chosen multiple k * P.  Again,
   a call to normalize_mont() after the multiply will give the
   explicit x-coordinate of the multiple.  Methods of
   choosing curves/points are embodied in the functions 

   		find_point_mont();  -- finds x of a point on curve
   		find_curve_point_brent(),  -- find a curve and x of a point

   Montgomery format is indicated for factoring
   and/or primality studies (where one usually may, with impunity, 
   ignore y-coordinates), yet sometimes also for cryptography 
   (such as x-embedding and El Gamal signature schemes).  
   The Montgomery operations are actually
   faster than the projective ones, but of course the projective 
   format involves more information (y coordinates).
  
   REFERENCES

   Brent R, Crandall R, Dilcher K, and Van Halewyn C 1997,
        "New factors of Fermat numbers," manuscript

   Crandall R, U.S. Patents #5159632 (1992), #5271061 (1993),
        #5463690 (1994), "Method and apparatus for public key 
        exchange in a cryptographic system."

   Crandall R, 1996 U. S. Patent #5581616, "Method and apparatus
        for Digital Signature Authentication."

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
#include "ellmont.h"



/* global variables */

static giant t0 = NULL, t1 = NULL, t2 = NULL, t3 = NULL, t4 = NULL,
	         t5 = NULL, t6 = NULL, t7 = NULL;

static point_mont ptmp0, ptmp1, ptmp2;

/**************************************************************
 *
 *	Maintenance functions
 *
 **************************************************************/

point_mont
new_point_mont(int shorts)
{
	point_mont pt;

	if(t0 == NULL) init_ell_mont(shorts);
	pt = (point_mont) malloc(sizeof(point_struct_mont));
	pt->x = newgiant(shorts);
	pt->z = newgiant(shorts);
	return(pt);
}

void
free_point_mont(point_mont pt)
{
	free(pt->x); free(pt->z);
	free(pt);
}

void
ptop_mont(point_mont pt1, point_mont pt2)
{
	gtog(pt1->x, pt2->x);
	gtog(pt1->z, pt2->z);
}

int
init_ell_mont(int shorts) 
/* Called by new_point_mont(), to set up giant registers. */
{	
	t0 = newgiant(shorts);
	t1 = newgiant(shorts);
	t2 = newgiant(shorts);
	t3 = newgiant(shorts);
	t4 = newgiant(shorts);
	t5 = newgiant(shorts);
	t6 = newgiant(shorts);
	t7 = newgiant(shorts);
    ptmp0 = new_point_mont(shorts);
    ptmp1 = new_point_mont(shorts);
    ptmp2 = new_point_mont(shorts);
}

/**************************************************************
 *
 *	Elliptic curve operations
 *
 **************************************************************/

/* Begin Montgomery-format functions for 
 
   y^2 = x^3 + c x^2 + a x + b.

   These are useful in cryptography, assuming y-coordinates
   are not needed.  The basic elliptic multiply is inversion-free.
   The functions are also useful to test other formats.  For
   example, x-coordinates herein should naturally agree with 
   Weierstrass-affine or Weierstrass-projective
   x-coordinate, when c = 0.  Alternatively, the Montgomery-Brent
   format, so useful for factoring, can be achieved with a = 1,
   b = 0, and c = (An/Ad-2), where Ad is a denominator term
   that eliminates all inversions throughout factoring.
   Similarly, other special curves, such as the CM curves
   for discriminants -3, -4, can be tested by such settings
   as (a,b,c) = (0,b,0) or (a,0,0), and so on.  For these
   reasons the routines herein have special case branching
   in case one or more parameters vanishes.

   The generalization of the Montgomery multiplication ladder
   is achieved through the relations (given here in TeX format): 

$$ X_+/Z_+ = elladd(X_1/Z_1, X_2/Z_2, X_-/Z_-)$$
$$ = {{Z_-} \over {X_-}}{{(X_1 X_2 - A Z_1 Z_2)^2 - 4 B(X_1 Z_2 + X_2 Z_1 + C Z_1 Z_2) Z_1 Z_2} \over {(X_1 Z_2 - X_2 Z_1)^2}}.$$

$$ X_+/Z_+ = elldouble(X_1/Z_1)$$
$$ = {{(X_1^2 - AZ_1^2)^2 - 4 B (2 X_1 + CZ_1) Z_1^3} \over {4 Z_1 (X_1^3 + CX_1^2 Z_1 + A X_1 Z_1^2 + BZ_1^3)}}.$$

   The elliptic arithmetic has bene optimized for the special cases,
   with the general case (i.e. general parameter triple (a,b,c))
   performed via numerator/denominator functions, between which
   some redundancy is expected.

 */

static void numer_double(giant x, giant z, giant res, giant a, giant b, giant c, giant N)
/* Numerator algebra, internal to ell_mul_mont() supporting routines:
   res := (x^2 - a z^2)^2 - 4 b (2 x + c z) z^3.
 */
{
    gtog(x, t1); squareg(t1); modg(N, t1);
    gtog(z, res); squareg(res); modg(N, res);
    gtog(res, t2);
    if(!isZero(a) ) {
        if(!isone(a)) { 
	    mulg(a, res); modg(N, res);
        }
        subg(res, t1); modg(N, t1);
    }
    squareg(t1); modg(N, t1);
    /* t1 := (x^2 - a z^2)^2. */
    if(isZero(b))  { 
	    gtog(t1, res);
		return;
    }
    if(!isZero(c)) { 
	gtog(z, res); mulg(c, res); modg(N, res);
    } else {
        itog(0, res);
    }
    addg(x, res); addg(x, res); mulg(b, res);
    modg(N, res);
    gshiftleft(2, res); mulg(z, res); modg(N, res);
    mulg(t2, res); modg(N, res);
    negg(res); addg(t1, res);
    modg(N, res);
}

static void denom_double(giant x, giant z, giant res, giant a, giant b, giant c, giant N)
/* Denominator algebra, internal to ell_mul_mont() supporting routines:
    res = 4 z (x^3 + c x^2 z + a x z^2 + b z^3). */
{
    gtog(x, res); gtog(z, t1);
    if(!isZero(c)) {
	gtog(c, t2); mulg(t1, t2); modg(N, t2);
	addg(t2, res);
    }
    mulg(x, res); modg(N, res);
    squareg(t1); modg(N, t1);
    if(!isZero(a)) {
	gtog(t1, t2);
    	mulg(a, t2); modg(N, t2);
    	addg(t2, res);
    }
    mulg(x, res); modg(N, res);
    if(!isZero(b)) {
	mulg(z, t1); modg(N, t1);
    	mulg(b, t1); modg(N, t1);
    	addg(t1, res);
    }
    mulg(z, res); gshiftleft(2, res);
    modg(N, res);
}

static void numer_add(giant x1, giant z1, giant x2, giant z2, giant res,
	giant a, giant b, giant c, giant N)
/* Numerator algebra, internal to ell_mul_mont() supporting routines:
    res := (x1 x2 - a z1 z2)^2 -
  	          4 b(x1 z2 + x2 z1 + c z1 z2) z1 z2
 */
{
    gtog(x1, t1); mulg(x2, t1); modg(N, t1);
    gtog(z1, t2); mulg(z2, t2); modg(N, t2);
    gtog(t1, res);
    if(!isZero(a)) {
	gtog(a, t3);
      	mulg(t2, t3); modg(N, t3);
      	subg(t3, res);
    }
    squareg(res); modg(N, res);
    if(isZero(b)) return;
    if(!isZero(c)) {
        gtog(c, t3);
    	mulg(t2, t3); modg(N, t3);
    } else itog(0, t3);
    gtog(z1, t4); mulg(x2, t4); modg(N, t4);
    addg(t4, t3);
    gtog(x1, t4); mulg(z2, t4); modg(N, t4);
    addg(t4, t3); mulg(b, t3); modg(N, t3);
    mulg(t2, t3); gshiftleft(2, t3); modg(N, t3);
    subg(t3, res);
    modg(N, res);
}

static void denom_add(giant x1, giant z1, giant x2, giant z2, giant res, giant N)
/* Denominator algebra, internal to ell_mul_mont() supporting routines:
    res := (x1 z2 - x2 z1)^2
 */
{
    gtog(x1, res); mulg(z2, res); modg(N, res);
    gtog(z1, t1); mulg(x2, t1); modg(N, t1);
    subg(t1, res); squareg(res); modg(N, res);
}


void
ell_even_mont(point_mont pt1, point_mont pt2, giant a, giant b, giant c, giant N)
/* Doubling routine, internal to Montgomery elliptic mul. */
{   
	giant x1 = pt1->x, z1 = pt1->z, x2 = pt2->x, z2 = pt2->z;
    
    if(isZero(b) && isone(a)) {  /* Montgomery-Brent case. */
	gtog(x1, t1); squareg(t1); modg(N, t1); /* t1 := x1^2. */
	gtog(z1, t2); squareg(t2); modg(N, t2); /* t2 := z1^2. */

    gtog(x1, t3); mulg(z1, t3); modg(N, t3);
	gtog(t3, z2); mulg(c, z2); modg(N, z2);
	addg(t1, z2); addg(t2, z2); mulg(t3, z2); gshiftleft(2, z2);
        modg(N, z2);  /* z2 := 4 x1 z1 (x1^2 + c x1 z1 + z1^2). */
        gtog(t1, x2); subg(t2, x2); squareg(x2); modg(N, x2);
						/* x2 := (x1^2 - z1^2)^2. */
     }
     else if(isZero(a) && isZero(c)) { /* CM curve case. */
	 gtog(x1, t1);
	 squareg(t1); modg(N, t1);
     mulg(x1, t1); modg(N, t1);   	/* t1 := x^3. */
     gtog(z1, t2);
	 squareg(t2); modg(N, t2);
     mulg(z1, t2); modg(N, t2);		/* t2 := z1^3 */
     mulg(b, t2); modg(N, t2); 	/* t2 := b z1^3. */
     gtog(t1, t3); addg(t2, t3);		/* t3 := x^3 + b z1^3 */	
     mulg(z1, t3); modg(N, t3);		/* t3 = z1 ( x^3 + b z1^3 ) */
	 gshiftleft(2, t3); modg(N, t3);	/* t3 = 4 z1 (x1^3 + b z1^3) */
     gshiftleft(3, t2);			/* t2 = 8 b z1^3 */
     subg(t2, t1);				/* t1 = x^3 - 8 b z1^3 */
	 mulg(x1, t1); modg(N, t1);		/* t1 = x1 (x1^3 - 8 b z1^3) */
	 gtog(t3, z2);
	 gtog(t1, x2);
    }
    else {
	numer_double(x1, z1, t0, a, b, c, N);
	denom_double(x1, z1, t7, a, b, c, N);
	gtog(t0, x2); gtog(t7, z2);
    }
}

void
ell_odd_mont(point_mont pt1, point_mont pt2, point_mont ptorg, giant a, giant b, giant c, giant N)
/* Adding routine, internal to Montgomery-Brent elliptic mul. */
{
	giant x1 = pt1->x, z1 = pt1->z, x2 = pt2->x, z2 = pt2->z,
		  xor = ptorg->x, zor = ptorg->z;

   if(isZero(b) && isone(a)) {  /* Montgomery-Brent case. */

	gtog(x1, t1); addg(z1, t1);  		/* t1 := x1 + z1. */
	gtog(x2, t2); subg(z2, t2);  		/* t2 := x2 - z2. */
	gtog(x1, t3); subg(z1, t3);  		/* t3 := x1 - z1. */
	gtog(x2, t4); addg(z2, t4);  		/* t4 := x2 + z2. */
	mulg(t2, t1); modg(N, t1);	   /* t1 := (x1 + z1)(x2 - z2) */
	mulg(t4, t3); modg(N, t3);	   /* t4 := (x2 + z2)(x1 - z1) */
	gtog(t1, z2); subg(t3, z2); /* z2 := ((x1 + z1)(x2 - z2) - x2) */
	squareg(z2); modg(N,  z2);
	mulg(xor, z2); modg(N, z2);
	gtog(t1, x2); addg(t3, x2);
	squareg(x2); modg(N, x2);
	mulg(zor, x2); modg(N, x2);
    }
    else if(isZero(a) && isZero(c)) {  /* CM curve case. */
	gtog(x1, t1); mulg(x2, t1);  modg(N, t1);  /* t1 := x1 x2. */
	gtog(z1, t2); mulg(z2, t2);  modg(N, t2);  /* t2 := z1 z2. */
	gtog(x1, t3); mulg(z2, t3);  modg(N, t3);  /* t3 := x1 z2. */
	gtog(z1, t4); mulg(x2, t4);  modg(N, t4);  /* t4 := x2 z1. */
	gtog(t3, z2); subg(t4, z2); squareg(z2); modg(N, z2);
	mulg(xor, z2); modg(N, z2);
	gtog(t1, x2); squareg(x2); modg(N, x2);
	addg(t4, t3); mulg(t2, t3); modg(N, t3);
	mulg(b, t3); modg(N, t3);
	addg(t3, t3); addg(t3, t3);
	subg(t3, x2); mulg(zor, x2); modg(N, x2);
    }
    else {
	numer_add(x1, z1, x2, z2, t0, a, b, c, N);
	mulg(zor, t0); modg(N, t0);
	denom_add(x1, z1, x2, z2, t7, N);
	mulg(xor, t7); modg(N, t7);
	gtog(t0, x2); gtog(t7, z2);
    }
}


void
ell_mul_mont(point_mont pt0, point_mont pt1, giant n, giant a, giant b, giant c, giant N)
/* General elliptic multiply, Montgomery format.
   pt1 := n * pt0,
   where n is a nonnegative integer. 
 */
{
	giant x0 = pt0->x, z0 = pt0->z, x1 = pt1->x, z1 = pt1->z;
    int bits = bitlen(n)-2;

    if(isZero(n)) {
		itog(0, z1);
		return;
    }
    gtog(x0, x1); gtog(z0, z1);
	if (isone(n)) return;
	if ((n->n[0] == 2) && (n->sign == 1))
	{
		ell_even_mont(pt1, pt1, a, b, c, N);
		return;
	}
	gtog(x0, ptmp0->x);
	gtog(z0, ptmp0->z);
    
	ell_even_mont(pt1, ptmp1, a, b, c, N);

	do {
	   if(bitval(n, bits--))
		{
			ell_odd_mont(ptmp1, pt1, ptmp0, a, b, c, N);
			ell_even_mont(ptmp1, ptmp1, a, b, c, N);
		}
		else
		{
			ell_odd_mont(pt1, ptmp1, ptmp0, a, b, c, N);
			ell_even_mont(pt1, pt1, a, b, c, N);
		}
     } while (bits >= 0);
}

void
normalize_mont(point_mont pt, giant p)
/* Obtain actual x,z values via normalization:
   {x,z} := {x/z, 1}.
 */

{	giant x = pt->x, z = pt->z;

	if(isZero(z)) {
		return;
	}
	binvaux(p, z); 
	mulg(z, x); modg(p, x);
	itog(1, z);
}


void
find_point_mont(point_mont pt, giant seed, giant a, giant b, giant c, giant p)
/* Starting with seed, finds a random (Montgomery) point {x,1} on curve y^2 = x^3 + c x^2 + a x + b.  This method is valid for any prime characteristic > 3, since we do not actually require a square root (y-coordinate) in the field.
 */
{	giant x = pt->x, z = pt->z;

	gtog(p, t1); iaddg(1, t1); gshiftright(1, t1); /* t1 := (p+1)/2. */

    while(1) {
		gtog(seed, x);
	    addg(c, x);
		mulg(seed,x); modg(p, x);
		addg(a, x);
		mulg(seed,x); addg(b, x);
		modg(p, x); /* x := seed^3 + c seed^2 + a seed + b. */
		gtog(x, z);
		powermodg(x, t1, p);
	    if(gcompg(z, x) == 0) break;
		iaddg(1, seed);
	}
	gtog(seed,x);
	itog(1, z);
}

/* Begin Montgomery-Brent-format functions for 
 
   y^2 = x^3 + (An/Ad-2) x^2 + x.

   These are especially useful in factoring studies.
   The parameter (An/Ad-2) has an explicit denominator Ad
   which one continually "carries around"
   in order to avoid comp;etely any inversions throughout a 
   factoring process.
 */

void
ell_even_brent(point_mont pt1, point_mont pt2, giant An, giant Ad, giant N)
/* Doubling routine, internal to Montgomery elliptic mul. */
{   
	giant x1 = pt1->x, z1 = pt1->z, x2 = pt2->x, z2 = pt2->z;
    
	gtog(x1, t1);
	addg(z1, t1);
	squareg(t1);
	modg(N, t1);
	gtog(x1, t2);
	subg(z1, t2);
	squareg(t2);
	modg(N, t2);
	gtog(t1, t3);
	subg(t2, t3);
	gtog(t2, x2);
	mulg(t1, x2);
	gshiftleft(2, x2);
	modg(N, x2);
	mulg(Ad, x2);
	modg(N, x2);
	mulg(Ad, t2);
	gshiftleft(2, t2);
	modg(N, t2);
	gtog(t3, t1);
	mulg(An, t1);
	modg(N, t1);
	addg(t1, t2);
	mulg(t3, t2);
	modg(N, t2);
	gtog(t2,z2);
}


void
ell_odd_brent(point_mont pt1, point_mont pt2, point_mont ptorg, giant N)
/* Adding routine, internal to Montgomery elliptic mul. */
{
	giant x1 = pt1->x, z1 = pt1->z, x2 = pt2->x, z2 = pt2->z,
		  xorg = ptorg->x, zorg = ptorg->z;

	gtog(x1, t1);
	subg(z1, t1);
	gtog(x2, t2);
	addg(z2, t2);
	mulg(t1, t2);
	modg(N, t2);
	gtog(x1, t1);
	addg(z1, t1);
	gtog(x2, t3);
	subg(z2, t3);
	mulg(t3, t1);
	modg(N, t1);
	gtog(t2, x2);
	addg(t1, x2);
	squareg(x2);
	modg(N, x2);
	gtog(t2, z2);
	subg(t1, z2);
	squareg(z2);
	modg(N, z2);
	mulg(zorg, x2);
	modg(N, x2);
	mulg(xorg, z2);
	modg(N, z2);
}


void
ell_mul_brent(point_mont pt0, point_mont pt1, giant n, giant An, giant Ad, giant N)
/* General elliptic multiply, Montgomery-Brent format.
   pt1 := n * pt0,
   where n is nonnegative. 
 */
{
	giant x0 = pt0->x, z0 = pt0->z, x1 = pt1->x, z1 = pt1->z;
    int b = bitlen(n)-2;

    if(isZero(n)) {
		itog(0, z1);
		return;
    }
    gtog(x0, x1); gtog(z0, z1);
	if (isone(n)) return;
	if ((n->n[0] == 2) && (n->sign == 1))
	{
		ell_even_brent(pt1, pt1, An, Ad, N);
		return;
	}
	gtog(x0, ptmp0->x);
	gtog(z0, ptmp0->z);
    
	ell_even_brent(pt1, ptmp1, An, Ad, N);

	do {
	   if(bitval(n, b--))
		{
			ell_odd_brent(ptmp1, pt1, ptmp0, N);
			ell_even_brent(ptmp1, ptmp1, An, Ad, N);
		}
		else
		{
			ell_odd_brent(pt1, ptmp1, ptmp0, N);
			ell_even_brent(pt1, pt1, An, Ad, N);
		}
     } while (b >= 0);
}

void
ell_mul_int_brent(point_mont pt, unsigned int n, giant An, giant Ad, giant N)
/* Simple elliptic multiply, for (unsigned int) multiples 
   n, Montgomery-Brent format.
   pt := n * pt,
   where n is a nonnegative integer. 
 */
{
    itog(n, t4);
	ell_mul_brent(pt, ptmp2, t4, An, Ad, N);
	ptop_mont(ptmp2, pt);
}

void
normalize_brent(point_mont pt, giant p)
/* Obtain actual x,z values via normalization:
   {x,z} := {x/z, 1}.
 */
{	normalize_mont(pt, p);
}

/* From R. P. Brent, priv. comm. 1996:
Let s > 5 be a pseudo-random seed,

	u/v = (s^2 - 5)/(4s)

Then starting point is {x_1, z_1} where

	x_1 = (u/v)^3
and
	a = (v-u)^3(3u+v)/(4u^3 v) - 2
*/

void
find_curve_point_brent(point_mont pt, int seed, giant An, giant Ad, giant N)
/* Find a curve and point living on found curve.
 */
{
	giant x = pt->x, z = pt->z;

	itog(seed, t5);
	gtog(t5, t4);
	squareg(t4);
	itog(5, t2);
	subg(t2, t4);
	modg(N, t4);
	addg(t5, t5);
	addg(t5, t5);
	modg(N, t5);
	gtog(t4, x);
	squareg(x);
	modg(N, x);
	mulg(t4, x);
	modg(N, x);
	gtog(t5, z);
	squareg(z);
	modg(N, z);
	mulg(t5, z);
	modg(N, z);

	/* Now for A. */
	gtog(t5, t2);
	subg(t4, t2);
	gtog(t2, t3);
	squareg(t2);
	modg(N, t2);
	mulg(t3, t2);
	modg(N, t2);  /* (v-u)^3. */
	gtog(t4, t3);
	addg(t3, t3);
	addg(t4, t3);
	addg(t5, t3);
	modg(N, t3);
	mulg(t3, t2);
	modg(N, t2);  /* (v-u)^3 (3u+v). */
	gtog(t5, t3);
	mulg(t4, t3);
	modg(N, t3);
	squareg(t4);
	modg(N, t4);
	mulg(t4, t3);
	modg(N, t3);
	addg(t3, t3);
	addg(t3, t3);
	modg(N, t3);
	gtog(t3, Ad);
	gtog(t2, An);  /* An/Ad is now c + 2. */
}



