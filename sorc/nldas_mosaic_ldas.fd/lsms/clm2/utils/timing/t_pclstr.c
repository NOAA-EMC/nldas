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
#include "gpt.h"

char *t_pclstr (int code)
{

#if ( defined DISABLE_TIMERS )
  return "";
#endif

#ifdef HAVE_PCL
  switch (code) {

  case PCL_SUCCESS: 
    return "Success";
    
  case PCL_NOT_SUPPORTED:
    return "Event not supported";
    
  case PCL_TOO_MANY_EVENTS:
    return "Too many events";
    
  case PCL_TOO_MANY_NESTINGS:
    return "More nesting levels than allowed";
    
  case PCL_ILL_NESTING:
    return "Bad nesting";
    
  case PCL_ILL_EVENT:
    return "Illegal event identifier";
    
  case PCL_MODE_NOT_SUPPORTED:
    return "Mode not supported";
    
  case PCL_FAILURE:
    return "Failure for unspecified reason";
    
  default:
    return "Unknown error code";
    
  }
#endif
}

      
  
