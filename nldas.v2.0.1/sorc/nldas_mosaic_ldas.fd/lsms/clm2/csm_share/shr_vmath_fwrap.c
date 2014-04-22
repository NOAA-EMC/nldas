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

/* =============================================================================
** Fortran wrappers for shr_vmath calls for systems that only
** provide a C interface to the system vector math routines.
** =============================================================================
*/

#if (defined IRIX64)

void shr_vmath_fwrap_vsqrt_(double *X, double *Y, int *n)
{
   vsqrt(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vexp_(double *X, double *Y, int *n)
{
   vexp(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vlog_(double *X, double *Y, int *n)
{
   vlog(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vsin_(double *X, double *Y, int *n)
{
   vsin(X, Y, *n, 1, 1);
}

void shr_vmath_fwrap_vcos_(double *X, double *Y, int *n)
{
   vcos(X, Y, *n, 1, 1);
}

#endif
