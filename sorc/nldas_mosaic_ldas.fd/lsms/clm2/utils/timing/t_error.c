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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

/*
** t_error: error return routine to print a message and return a failure
** value.
**
** Input arguments:
**   fmt: format string
**   variable list of additional arguments for vfprintf
**
** Return value: -1 (failure)
*/

int t_error (const char *fmt, ...)
{
  va_list args;

#if ( ! defined DISABLE_TIMERS )
  va_start (args, fmt);

  if (fmt != NULL)
    (void) vfprintf (stderr, fmt, args);

  va_end (args);
#endif
  
  return (-1);
}
