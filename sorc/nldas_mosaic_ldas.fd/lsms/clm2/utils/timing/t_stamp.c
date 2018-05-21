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

#include <sys/time.h>     /* gettimeofday */ 
#include <sys/times.h>    /* times */
#include <unistd.h>       /* gettimeofday */

#include "gpt.h"

/*
** t_stamp: Compute timestamp of usr, sys, and wallclock time (seconds)
**
** Output arguments:
**   wall: wallclock
**   usr:  user time
**   sys:  system time
**
** Return value: 0 (success) or -1 (failure)
*/

int t_stamp (double *wall, double *usr, double *sys)
{
  struct timeval tp;         /* argument to gettimeofday */
  struct tms buf;            /* argument to times */

  *usr = 0;
  *sys = 0;

  if (times (&buf) == -1)
    return t_error ("t_stamp: times() failed. Timing bogus\n");

  *usr = buf.tms_utime / (double) ticks_per_sec;
  *sys = buf.tms_stime / (double) ticks_per_sec;

  gettimeofday (&tp, NULL);
  *wall = tp.tv_sec + 1.e-6*tp.tv_usec;

  return 0;
}
