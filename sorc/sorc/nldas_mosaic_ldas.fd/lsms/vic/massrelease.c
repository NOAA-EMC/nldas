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
 * SUMMARY:      MassRelease.c - Calculates mass release of snow from canopy
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Brian Connelly and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     6-Oct-1996 at 15:42:13
 * LAST-MOD: Mon Sep 28 16:21:39 1998 by VIC Administrator <vicadmin@u.washington.edu>
 * DESCRIPTION:  Calculates mass release of snow from canopy
 * DESCRIP-END.
 * FUNCTIONS:    MassRelease()
 * COMMENTS:     
 */

#include <stdio.h>
#include "vicNl.h"
  

/*****************************************************************************
  Function name: MassRelease()

  Purpose      : Calculates mass release of snow from canopy

  Required     :
    float *InterceptedSnow
    float *TempInterceptionStorage
    float *ReleasedMass
    float *Drip 

  Returns      : none

  Modifies     : see under required (i.e. all variables are modified)

  Comments     :
*****************************************************************************/
void MassRelease(float *InterceptedSnow, float *TempInterceptionStorage, 
                 float *ReleasedMass, float *Drip ) 
{
  float TempDrip;
  float TempReleasedMass;
  float Threshold;
  float MaxRelease;
  
  /* If the amount of snow in the canopy is greater than some minimum
     value, MIN_INTERCEPTION_STORAGE, then calculte mass release and Drip */
  
  if (*InterceptedSnow > MIN_INTERCEPTION_STORAGE) {
    Threshold  = 0.10 * *InterceptedSnow;
    MaxRelease = 0.17 * *InterceptedSnow;
    
    /* If the amount of snow_melt after interception, snow_melt, is >= the
       theshhold then there is mass release.  If snow_melt is < the treshhold
       then there is no mass release but that water remains in
       *TempInterceptionStorage which will be augmented during the next
       compute period */
    
    if ((*TempInterceptionStorage) >= Threshold) {
      
      
      *Drip += Threshold;
      *InterceptedSnow -= Threshold;
      *TempInterceptionStorage -= Threshold;
      if (*InterceptedSnow < MIN_INTERCEPTION_STORAGE)
        TempReleasedMass = 0.0;
      else
        TempReleasedMass = min((*InterceptedSnow - MIN_INTERCEPTION_STORAGE),
                               MaxRelease); 
      *ReleasedMass += TempReleasedMass;
      *InterceptedSnow -= TempReleasedMass;
      MassRelease(InterceptedSnow, TempInterceptionStorage, ReleasedMass,
                  Drip); 
    }

    else {
      TempDrip = min(*TempInterceptionStorage, *InterceptedSnow);
      *Drip += TempDrip;
      *InterceptedSnow -= TempDrip;
    }
  }

  /* (*InterceptedSnow < MIN_INTERCEPTION_STORAGE) If the amount of snow in
     the canopy is less than some minimum value, MIN_INTERCEPTION_STORAGE,
     then only melt can occur and there is no mass release. */

  else {
    TempDrip = min(*TempInterceptionStorage, *InterceptedSnow);
    *Drip += TempDrip;
    *InterceptedSnow -= TempDrip;
    *TempInterceptionStorage = 0.0;
  }
}
