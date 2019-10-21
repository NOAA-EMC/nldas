#!  /bin/csh -fx
##=========================================================================
##
##  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  A U.S. Continental-Scale   
##  D                                      L  Land Modeling and Data 
##  A  --LAND DATA ASSIMILATION SCHEMES--  D  Assimilation Project.
##  S                                      A  This is the GSFC-LDAS Code. 
##  LDASLDASLDASLDASLDASLDASLDASLDASLDASLDAS  http://ldas.gsfc.nasa.gov
##
##   GSFC - NCEP - OH - Princeton - Washington - Rutgers
##
##=========================================================================
## comp.csh: 
##
## DESCRIPTION:
##  UNIX script for compiling GSFC-LDAS code.
##
## REVISION HISTORY:
##  1 Oct 1999: Jared Entin; Initial code
##  15 Oct 1999: Paul Houser; General Revision
##  31 Jan 2003: Jon Gottschalck; Updated for CLM2 checks
##  14 Mar 2019  Youlong Xia, WCOSS phase 3
##=========================================================================
cd MAKDEP
gmake
cd ..
cp ./MAKDEP/makdep makdep
gmake clean
gmake











