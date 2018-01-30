/* AKS.c
   Implementation of the AKS primality-proving algorithm.
   
   This code implements an unpublished variant due to
   H. W. Lenstra, Jr.  In said variant, one establishes
   that p is not a prime power a^(b>1), then defines:
   
   v > Round[Log[2,p]]^2.
   
   Then if r is such that ord(p, r) > v,
   we check the relation
   
   (x-a)^p = x^p-a (mod (x^r-1, p))
   
   for a = 1 through phi(r)-1.
   In the code, we force r prime (so that phi(r) = r-1, simply)
   and perform fast polynomial arithmetic, as in
   Algorithm 9.6.1, and using other algorithms from
   Chapter 9, of
   R. Crandall and C. Pomerance, Prime numbers: a computational
   perspective, Springer 2001/2.
   
   Usage: (link files taken from "PrimeKit," at www.perfsci.com)
   compile:
   $ cc -O AKS.c giants.c tools.c
   run:
   $ a.out
   Enter a candidate which is not a prime power:
   
   This is not "production" code; rather, it is a quick experiment
   to assess the actual runtime of the Lenstra variant.
   Lenstra proves O~(log^8 p) (the original AKS paper was O~(log^12 p)
   and it is believed that the basic ideas are really yielding
   O~(log^6 p) for the runtime.
   Some example primes p and times T to prove, on a decent
   workstation, are as follows:
   
   p        T (seconds, very roughly, on a "decent workstation")
   
   7001                                 1
   70001                                3
   700001                              15
   2147483647                         200
   1125899906842679                  4000  ("about an hour")         
   618970019642690137449562111     100000  ("about a day")
   
   Note that the timing behavior follows pretty well the rule
   
   T ~ C log^6 p,
   
   where C ~ 10^(-6) seconds.
   
   A note is in order on effective complexity.  It is tantalizing
   that the tough loop (below: "for(a = 1; ...") is
   "embarrallelm" meaning embarrassingly parallelizable.
   That is, one can take say phi(r)-1 machines, and
   just do a single a value on each machine.  In this way,
   given the above timings, a 1000-digit prime might be doable
   in "about a year."
   
   Thanks to C. Pomerance, J. Buhler, H. Lenstra for ideas 
   relating to this experiment.
   
   R. E. Crandall
   Center for Advanced Computation
   Sep 2002
   
*/

#include <stdio.h>
#include<assert.h>
#include <math.h>
#include <stdlib.h>
#include "giants.h"

#define P_BREAK 32
#define MAX_DIGS 24
#define MAX_COEFFS (1<<19)
typedef struct
         {
         int deg;    
         giant *coe;
         } polystruct;   
typedef polystruct *poly;

static poly p1, pbuff, pbuff2;
static giant p, coe, tmp, *mcand, globx, globy;
static unsigned int r;

void mulp(poly x, poly y);
giant *newa(int n);

poly
newpoly(int coeffs) {
        poly pol;
        pol = (poly) malloc(sizeof(polystruct));
        pol->coe = (giant *)newa(coeffs);
        return(pol);
}

int
log_2(int n) {
        int c = 1, d = -1;
        while(c <= n) {
                c <<= 1;
                ++d;
        }
        return)d);   
}

void
justifyp(poly x) {
        int j;
        for(j = x->deg; j >= 0; j--) {
                if(!isZero(x->coe[j])) break;
        }
        x->deg = (j>0)? j : 0;
}


void 
atoa(giant *a, giant *b, int n) {
        int j;
        for(j=0; j<n; j++) gtog(a[j], b[j]);
}

void
ptop(poly x, poly y)
/* y := x. */
{
        y->deg = x->deg;
        atoa(x->coe, y->coe, y->deg+1);
}

void quickmodg(giant g, giant x) 
{       int sgn = x->sign;

        if(sgn == 0) return;
        if(sgn > 0) {
                if(gcompg(x, g) >= 0) subg(g, x);
                return;
        }
        addg(g,x);
        return;
}


void
addp(poly x, poly y)
/* y += x. */
{
        int d = x->deg, j;

        if(y->deg > d) d = y->deg;
        for(j = 0; j <= d; j++) {
                if((j <= x->deg) && (j <= y->deg)) {
                        addg(x->coe[j], y->coe[j]);
                        quickmodg(p, y->coe[j]);
                        continue;
                }
                if((j <= x->deg) && (j > y->deg)) {
                        gtog(x->coe[j], y->coe[j]);
                        quickmodg(p, y->coe[j]);
                        continue;
                }
        }
        y->deg = d;
        justifyp(y);
}

void
grammarmulp(poly a, poly b) 
/* b *= a. */
{
        int dega = a->deg, degb = b->deg, deg = dega + degb;
        register int d, kk, first, diffa;

        for(d=deg; d>=0; d--) {
                diffa = d-dega;
                itog(0, coe);
                for(kk=0; kk<=d; kk++) {
                        if((kk>degb)||(kk<diffa)) continue;
                        gtog(b->coe[kk], tmp);
                        mulg(a->coe[d-kk], tmp);
                        modg(p, tmp);
                        addg(tmp, coe);
                        quickmodg(p, coe);
                }
                gtog(coe, mcand[d]);
        }
        atoa(mcand, b->coe, deg+1);
        b->deg = deg;
        justifyp(b);
}

void
polyrem(poly x) {
        int j;
   for(j=0; j <= x->deg; j++) {
      modg(p, x->coe[j]);
   }
   justifyp(x);
}

fullmod(poly x) {
   unsigned int j, m;
   
   polyrem(x);
   if(x->deg < r) return;
   for(j=r; j <= x->deg; j++) { /* Perform reduction mod x^r-1. */
         m = j%r;
         addg(x->coe[j], x->coe[m]);
         if(gcompg(x->coe[m], p) >= 0) subg(p, x->coe[m]);
   }
   x->deg = r-1;
   justifyp(x);
}

void
mulmod(poly x, poly y) {
        mulp(x, y); fullmod(y);
}

void
just(giant g) {
   while((g->n[g->sign-1] == 0) && (g->sign > 0)) --g->sign;
}

void
unstuff_partial(giant g, poly y, int words){
        int j;
        for(j=0; j < y->deg; j++) {
                bcopy((g->n) + j*words, y->coe[j]->n, words*sizeof(short));
      y->coe[j]->sign = words;
      just(y->coe[j]);
   }
}

void
stuff(poly x, giant g, int words) {
        int deg = x->deg, j, coedigs;

   g->sign = words*(1 + deg);
   for(j=0; j <= deg; j++) {
                coedigs = (x->coe[j])->sign;
                bcopy(x->coe[j]->n, (g->n) + j*words, coedigs*sizeof(short));
                bzero((g->n) + (j*words+coedigs), 
                                sizeof(short)*(words-coedigs));
        }
   just(g);
}

void
binarysegmul(poly x, poly y) {
   int bits = bitlen(p), xwords, ywords, words;
   
   xwords = (2*bits + log_2(x->deg+1) + 32 + 15)/16;
   ywords = (2*bits + log_2(y->deg+1) + 32 + 15)/16;
   if(xwords > ywords) words = xwords; else words = ywords;
   stuff(x, globx, words);
   stuff(y, globy, words);
   mulg(globx, globy);
   gtog(y->coe[y->deg], globx);  /* Save high coeff. */
   y->deg += x->deg;
   gtog(globx, y->coe[y->deg]);  /* Move high coeff. */
   unstuff_partial(globy, y, words);
   mulg(x->coe[x->deg], y->coe[y->deg]); /* resolve high coeff. */
   justifyp(y);
}

binarysegsquare(poly y) {
   int bits = bitlen(p), words;
      words = (2*bits + log_2(y->deg+1) + 32 + 15)/16;
      
   stuff(y, globy, words);
   squareg(globy);
   gtog(y->coe[y->deg], globx);  /* Save high coeff. */
   y->deg += y->deg;
   gtog(globx, y->coe[y->deg]);  /* Move high coeff. */
   unstuff_partial(globy, y, words);
   mulg(y->coe[y->deg], y->coe[y->deg]); /* resolve high coeff. */
   justifyp(y);
}

void
mulp(poly x, poly y)
/*  y *= x. */
{
   int n, degx = x->deg, degy = y->deg;

   if((degx < P_BREAK) || (degy < P_BREAK)) {
                grammarmulp(x,y);
                justifyp(y);
                return;
        }
   if(x==y) binarysegsquare(y);
   	else binarysegmul(x, y);
}

void powerpoly(poly x, giant n)
/* Perform windowed ladder. */
{       int pos, code;
        ptop(x, pbuff);  
	ptop(pbuff, pbuff2);
	mulmod(pbuff2, pbuff2); mulmod(pbuff, pbuff2);
        pos = bitlen(n)-2;
        while(pos >= 0) {
		mulmod(x, x);
		if(pos==0) {
			if(bitval(n, pos) != 0) {
				mulmod(pbuff, x);
			}
			break;
		}
		code = (bitval(n, pos) != 0) * 2 + (bitval(n, pos-1) != 0);
		switch(code) {
			case 0: mulmod(x,x); break;
			case 1: mulmod(x,x); 
				mulmod(pbuff, x);
				break;
			case 2: mulmod(pbuff, x); 
				mulmod(x,x); break;
			case 3: mulmod(x,x); mulmod(pbuff2, x); break;
		}
		pos -= 2;
        }
}



giant *
newa(int n) {
        giant *p = (giant *)malloc(n*sizeof(giant));
        int j;
        for(j=0; j<n; j++) {
                p[j] = newgiant(MAX_DIGS);
        }
        return(p);
}

init_polys(unsigned int r) {
    p1 = newpoly(4*r+20);
    pbuff = newpoly(4*r + 20);
    pbuff2 = newpoly(4*r+20);
    globx = newgiant(0);
    globy = newgiant(0);
    mcand = (giant *)newa(MAX_COEFFS);

}

void
polyout(poly x) {
   int j;
   for(j=0; j <= x->deg; j++) {printf("%d: ",j); gout(x->coe[j]);}
}

int
ord(giant p, unsigned int r) {
    unsigned int ct = 0;
    long long x,y=1;

    gtog(p, tmp);    
    x = idivg(r, tmp);
    while(1) {
        y = (x*y) % r;
        ++ct;
        if(y==1) break;
    }
    return ct;
}

main(int argc, char **argv) {
        unsigned int phi, a, j, rem, v;
        giant am = newgiant(0);

       init_tools(1);
       p = newgiant(0); 
       coe = newgiant(0);
       tmp = newgiant(0);

       printf("Enter a candidate p which is not a prime power:\n");
       fflush(stdout);
       gin(p);
       v = bitlen(p); v = v*v;
       for(r=v+1; ; r++){
          if(!primeq(r)) continue;
          if(ord(p, r) <= v) continue;
          break;
       }
       phi = r-1;
       printf("r phi v: %d %d %d\n", r, phi, v); fflush(stdout);
       init_polys(r);
       gtog(p, tmp);
       rem = idivg(r, tmp);
       for(a = 1; a <  phi; a++) {
        	p1->deg = 1;
        	itog(1, p1->coe[1]);
                itog(a, am); negg(am); addg(p, am);
                gtog(am, p1->coe[0]);
		powerpoly(p1, p);
                if(p1->deg > rem) goto COMPOSITE;
                if(!isone(p1->coe[rem])) goto COMPOSITE;
                if(gcompg(p1->coe[0], am) != 0) goto COMPOSITE;
	        for(j=1; j < rem; j++) {
			if(!isZero(p1->coe[j])) goto COMPOSITE;
		}
                printf("%d ", a); fflush(stdout);
                if(a % 8 == 0) {printf("\n"); fflush(stdout);}
	}
        printf("\np is PRIME\n"); gout(p);
        exit(0);
COMPOSITE:
        printf("\np is COMPOSITE: %d \n", a); gout(p);
        exit(0);
}


      		
        