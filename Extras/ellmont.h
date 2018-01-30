/**************************************************************
 *
 *	ellmont.h
 *
 *	Header file for ellmont.c
 *
 *	Updates:
 *		3 Apr 98  REC - Creation
 *
 *	c. 1998 Perfectly Scientific, Inc.
 *	All Rights Reserved.
 *
 *
 *************************************************************/

/* definitions */

typedef struct
{
	 giant x;
	 giant z;
} point_struct_mont;

typedef point_struct_mont *point_mont;

point_mont
new_point_mont(int shorts);

void
free_point_mont(point_mont pt);

void
ptop_mont(point_mont pt1, point_mont pt2);

int
init_ell_mont(int shorts);

static void numer_double(giant x, giant z, giant res, giant a, giant b, giant c, giant N);

static void denom_double(giant x, giant z, giant res, giant a, giant b, giant c, giant N);

static void numer_add(giant x1, giant z1, giant x2, giant z2, giant res,
	giant a, giant b, giant c, giant N);

static void denom_add(giant x1, giant z1, giant x2, giant z2, giant res, giant N);

void
ell_even_mont(point_mont pt1, point_mont pt2, giant a, giant b, giant c, giant N);

void
ell_odd_mont(point_mont pt1, point_mont pt2, point_mont ptorg, giant a, giant b, giant c, giant N);

void
ell_mul_mont(point_mont pt0, point_mont pt1, giant n, giant a, giant b, giant c, giant N);

void
normalize_mont(point_mont pt, giant p);

void
find_point_mont(point_mont pt, giant seed, giant a, giant b, giant c, giant p);

void
ell_even_brent(point_mont pt1, point_mont pt2, giant An, giant Ad, giant N);

void
ell_odd_brent(point_mont pt1, point_mont pt2, point_mont ptorg, giant N);

void
ell_mul_brent(point_mont pt0, point_mont pt1, giant n, giant An, giant Ad, giant N);

void
normalize_brent(point_mont pt, giant p);

void
find_curve_point_brent(point_mont pt, int seed, giant An, giant Ad, giant N);

