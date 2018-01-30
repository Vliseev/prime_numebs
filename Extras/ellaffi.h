/**************************************************************
 *
 *	ellaffi.h
 *
 *	Header file for ellaffi.c
 *
 *	Updates:
 *		3 Apr 98    REC - Creation
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
	 giant y; 
	 int  z;
} point_struct_affi;

typedef point_struct_affi *point_affi;

point_affi
new_point_affi(int shorts);

void
free_point_affi(point_affi pt);

void
ptop_affi(point_affi pt1, point_affi pt2);

int
init_ell_affi(int shorts);


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
ell_double_affi(point_affi pt, giant a, giant p);
 
void
ell_add_affi(point_affi pt0, point_affi pt1, giant a, giant p);

void
ell_neg_affi(point_affi pt, giant p);

void
ell_sub_affi(point_affi pt0, point_affi pt1, giant a, giant p);

void
ell_mul_affi(point_affi pt0, point_affi pt1, giant k, giant a, giant p);

void
find_point_affi(point_affi pt, giant seed, giant a, giant b, giant p);

