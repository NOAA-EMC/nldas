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

#include <stdio.h>
#include "gpt.h"

/*
** t_setoption: set option value to true or false.
**
** Input arguments:
**   option: option to be set
**   val:    value to which option should be set
**
** Return value: 0 (success) or -1 (failure)
*/

int t_setoption (OptionName option, Boolean val)
{
  int n;

#if ( defined DISABLE_TIMERS )
  return 0;
#endif

  if (t_initialized)
    return (t_error ("t_setoption: Options must be set BEFORE t_initialize\n"));

  for (n = 0; n < npossible; n++) {
    if (possible_event[n].name == option) {
      possible_event[n].enabled = val;

      if (val)
	printf ("t_setoption: option enabled:  %s\n", possible_event[n].string);
      else
	printf ("t_setoption: option disabled: %s\n", possible_event[n].string);

      return 0;
    }
  }

  return (t_error ("t_setoption: Option with enum index %d not available\n",
		     option));
}
