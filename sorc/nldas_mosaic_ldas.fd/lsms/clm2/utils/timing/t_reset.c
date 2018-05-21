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

#include <string.h>  /* memset */
#include <stdio.h>

#include "gpt.h"

/*
** t_reset: reset all known timers to 0
**
** Return value: 0 (success) or -1 (failure)
*/

int t_reset ()
{
  int n;             /* index over threads */
  struct node *ptr;  /* linked list index */

#if ( ! defined DISABLE_TIMERS )
  if ( ! t_initialized)
    return t_error ("t_reset: t_initialize has not been called\n");

  /*
  ** Only allow the master thread to reset timers
  */

  if (get_thread_num () != 0)
    return 0;

  for (n = 0; n < numthreads; n++) {
    for (ptr = timers[n]; ptr != NULL; ptr = ptr->next) {
      memset (timers[n], 0, sizeof (struct node));
      printf ("Reset accumulators for timer %s to zero\n", ptr->name);
    }
  }
#endif
  return 0;
}
