//-------------------------------------------------------------------------
// NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
// Released October 2005
//
// See SOFTWARE DISTRIBUTION POLICY for software distribution policies
//
// The LIS source code and documentation are in the public domain,
// available without fee for educational, research, non-commercial and
// commercial purposes.  Users may distribute the binary or source
// code to third parties provided this statement appears on all copies and
// that no charge is made for such copies.
//
// NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
// SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
// IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
// LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
//
// See COPYRIGHT.TXT for copyright details.
//
//-------------------------------------------------------------------------
/*
 * SUMMARY:      RootBrent.c - Determine surface temperature iteratively
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * LAST-MOD: Mon Jan 24 12:06:38 2000 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Determine surface temperature iteratively using the Brent
 *               method.  
 * DESCRIP-END.
 * FUNCTIONS:    RootBrent()
 * COMMENTS:     
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vicNl.h"


/*****
#include <stdarg.h>
#include "settings.h"
#include "brent.h"
#include "massenergy.h"
#include "DHSVMerror.h"
*****/

#define MAXTRIES 5
#define MAXITER 1000
#define MACHEPS 3e-8
#define TSTEP   10
#define T       1e-7   

/*****************************************************************************
  GENERAL DOCUMENTATION FOR THIS MODULE
  -------------------------------------

  Source: Brent, R. P., 1973, Algorithms for minimization without derivatives,
                        Prentice Hall, Inc., Englewood Cliffs, New Jersey
			Chapter 4
  This source includes an implementation of the algorithm in ALGOL-60, which
  was translated into C for this application.

  The method is also discussed in:
  Press, W. H., S. A. Teukolsky, W. T. Vetterling, B. P. Flannery, 1992,
                Numerical Recipes in FORTRAN, The art of scientific computing,
		Second edition, Cambridge University Press
  (Be aware that this book discusses a Brent method for minimization (brent), 
  and one for root finding (zbrent).  The latter one is similar to the one 
  implemented here and is also copied from Brent [1973].)

  The function returns the surface temperature, TSurf, for which the sum
  of the energy balance terms is zero, with TSurf in the interval 
  [MinTSurf, MaxTSurf].  The surface temperature is calculated to within
  a tolerance (6 * MACHEPS * |TSurf| + 2 * T), where MACHEPS is the relative
  machine precision and T is a positive tolerance, as specified in brent.h.

  The function assures that f(MinTSurf) and f(MaxTSurf) have opposite signs.
  If this is not the case the program will abort.  In addition the program
  will perform not more than a certain number of iterations, as specified
  in brent.h, and will abort if more iterations are needed.
******************************************************************************/
  
/*****************************************************************************
  Function name: RootBrent()

  Purpose      : Calculate the surface temperature in the absence of snow

  Required     :
    int y                 - Row number of current pixel
    int x                 - Column number of current pixel 
    float LowerBound      - Lower bound for root
    float UpperBound      - Upper bound for root
    float (*Function)(float Estimate, va_list ap)
    ...                   - Variable arguments 
                            The number and order of arguments has to be
                            appropriate for the Function pointed to, since
                            the list of arguments after Nargs will be passed
                            on to Function.
                            See the appropriate Function for the correct
                            arguments. 

  Returns      :
    float b               - Effective surface temperature (C)

  Modifies     : none

  Comments     :
*****************************************************************************/
float root_brent(float LowerBound, float UpperBound, 
                float (*Function)(float Estimate, va_list ap), ...)
{
  const char *Routine = "RootBrent";
  char ErrorString[MAXSTRING+1];
  va_list ap;                   /* Used in traversing variable argument list
                                 */ 
  float a;
  float b;
  float c;
  float d;
  float e;
  float fa;
  float fb;
  float fc;
  float m;
  float p;
  float q;
  float r;
  float s;
  float tol;
  int i;
  int j;
  int eval = 0;

  /* sprintf(ErrorString, "%s", Routine); */

  /* initialize variable argujment list */

    
  a = LowerBound;
  b = UpperBound;;
  va_start(ap, Function);
  fa = Function(a, ap);
  eval++;
  va_start(ap, Function);
  fb = Function(b, ap);
  eval++;
  
  /*  if root not bracketed attempt to bracket the root */

  j = 0;
  while ((fa * fb) >= 0  && j < MAXTRIES) {
    a -= TSTEP;
    b += TSTEP;
    va_start(ap, Function);
    fa = Function(a, ap);
    eval++;
    va_start(ap, Function);
    fb = Function(b, ap);
    eval++;
    j++;
  }
  if ((fa * fb) >= 0) {
    sprintf(ErrorString,"ERR: %s -- first error in root_brent -> %f * %f >= 0",
        Routine, fa, fb);
    lis_log_msgC(ErrorString);
    /* fprintf(stderr,"Dumping Input Variables - Check for Valid Values\n"); */
    return(MISSING);
  }
  
  fc = fb;

  for (i = 0; i < MAXITER; i++) {

    if (fb*fc > 0) {
      c = a;
      fc = fa;
      d = b - a;
      e = d;
    }
    
    if (fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    
    tol = 2 * MACHEPS * fabs(b) + T;
    m = 0.5 * (c - b);
    
    if (fabs(m) <= tol || fb == 0) {
      va_end(ap);
      return b;
    }
    
    else {
      if (fabs(e) < tol || fabs(fa) <= fabs(fb)) {
	d = m;
	e = d;
      }
      else {
	s = fb/fa;
	
	if (a == c) {
	  
	  /* linear interpolation */
          
	  p = 2 * m * s;
	  q = 1 - s;
	}
	
	else {
	  
	  /* inverse quadratic interpolation */
	  
	  q = fa/fc;
	  r = fb/fc;
	  p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
	  q = (q - 1) * (r - 1) * (s - 1);
	}
	
	if (p > 0)
	  q = -q;
	else
	  p = -p;
	s = e;
	e = d;
	if ((2 * p) < ( 3 * m * q - fabs(tol * q)) && p < fabs(0.5 * s * q))
	  d = p/q;
	else {
	  d = m;
	  e = d;
	}
      }
      a = b;
      fa = fb;
      b += (fabs(d) > tol) ? d : ((m > 0) ? tol : -tol);
      va_start(ap, Function);
      fb = Function(b, ap);      
      eval++;
    }
  }
  lis_log_msgC("ERR: RootBrent -- second error in root_brent -> Too many iterations");
  /* fprintf(stderr,"Dumping Input Variables - Check for Valid Values\n"); */
  return(MISSING);
}

#undef MAXTRIES
#undef MAXITER
#undef MACHEPS
#undef TSTEP
#undef T
