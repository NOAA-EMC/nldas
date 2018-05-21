#!/bin/sh
#
# Start of script compallip_64.sh         B. Vuong 04-01-1998 
#
# This script is used to compile all F90 subroutines and place in the 
#    library W3LIB90.
#    Note: 
#          The F90 library is a "portable archive" -vs- the previous
#          "build" library.  Compiler option "-O -64 -i8 -r8" has been included
#          to default to 64-bit integers (CDIR$ INTEGER=64 no longer allowed
#          in F90 codes as well as "-i 64" at compilation time).
#
# Example:  Compile subroutine w3fq02.f and put it into library iplib_64
#           compallip_64.sh w3fq02
#
# rm    $name.o               # delete subroutine old object
# set -xS
#
list=`ls *.f | grep '.f' `
for routine in $list
do
name=`echo $routine | cut -f1 -d. `
echo ' '
echo Compile and place in library iplib_64 subroutine $name.f
#if f90 -c -O -64 -i8 -r8 $name.o $name.f # compile subroutine
if f90 -c $name.o $name.f # compile subroutine

then
ar r ../iplib_64 $name.o
rm    $name.o               # delete subroutine new object
#
echo Subroutine $name.f was compiled and placed in library iplib_64
echo ' '
else
echo ' '
echo Compile error, subroutine $name.f was not placed in library iplib_64
echo ' '
exit
fi
done
#
# End of script compallip_64.sh
