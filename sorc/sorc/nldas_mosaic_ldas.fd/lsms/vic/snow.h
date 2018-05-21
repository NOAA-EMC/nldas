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
 * SUMMARY:      snow.h - header file for DHSVM snow routines
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              nijssen@u.washington.edu
 * ORIG-DATE:    29-Aug-1996 at 16:03:11
 * LAST-MOD: Tue Jul 14 10:38:32 1998 by Keith Aric Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  header file for DHSVM snow routines
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     
 */

#ifndef SNOW_H
#define SNOW_H

#include <stdarg.h>

/* water holding capacity of snow as a fraction of snow-water-equivalent */ 
#define LIQUID_WATER_CAPACITY	0.035

/* multiplier to calculate the amount of available snow interception as
   a function of LAI (m) */
#define LAI_SNOW_MULTIPLIER	0.0005

/* the amount of snow on the canopy that can only be melted off. (m) */
#define MIN_INTERCEPTION_STORAGE	0.005

/* maximum depth of the surface layer in water equivalent (m)
   [default 0.125] */
#define MAX_SURFACE_SWE	0.125

/* density of new fallen snow [50] */
#define NEW_SNOW_DENSITY	50.
 
/* Minimum SWQ at which the snow pack is assumed to fully cover the grid cell (m) */
#define MAX_FULL_COVERAGE_SWQ   0.076

#endif 
