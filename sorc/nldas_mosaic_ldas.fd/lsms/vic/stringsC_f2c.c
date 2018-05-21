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
#include <string.h>
#include <stdio.h>
#include "ftn.h"
#include "lis_headers.h"
//BOP
//
// !ROUTINE: extract_string_f2c
//
// !DESCRIPTION: 
// This routine extracts a substring from the given packed (pack_string_f2c)
// string.
//
// !REVISION HISTORY: 
// 02 Mar 2004; James Geiger :Initial Specification
//
// !INTERFACE:
char * extract_string_f2c(int n, char * string)
//EOP
{
//BOC
   int len;
   char * substring;
   int i, j, count;

   len = strlen(string);

   substring = (char *) lis_malloc(len,"extract_string_f2c");

   count = 0;
   i = 0;
   while ( count < n-1 )
   {
      if ( string[i++] == ' ' )
      {
         ++count;
      }
   } 
   
   j = 0;
   while ( (substring[j] = string[i]) != ' ' )
   {
      ++i;
      ++j;
   }
   substring[j] = '\0';

   return substring;
//EOC
}
