#!/bin/sh
#
# Start of script compallip_32.sh         B. Vuong 04-01-1998 
#
# This script is used to compile all F90 subroutines and place in the 
#    library W3LIB90.
#    Note: 
#          The F90 library is a "portable archive" -vs- the previous
#          "build" library.  Compiler option "-O -n32" has been included
#          to default to 32-bit integers.
#
# Example:  Compile subroutine w3fq02.f and put it into library iplib_32
#           compallip_32.sh w3fq02
#
# rm    $name.o               # delete subroutine old object
# set -xS
#
rm ../iplib_32_alpha
list=`ls *.f | grep '.f' `
for routine in $list
do
name=`echo $routine | cut -f1 -d. `
echo ' '
echo Compile and place in library iplib_32 subroutine $name.f
#if f90 -c -O -n32 $name.o $name.f # compile subroutine
if f90 -assume byterecl -convert big_endian -arch ev67 -tune ev67 -O2 -c $name.o $name.f # compile subroutine

then
ar r ../iplib_32_alpha $name.o
rm    $name.o               # delete subroutine new object
#
echo Subroutine $name.f was compiled and placed in library iplib_32
echo ' '
else
echo ' '
echo Compile error, subroutine $name.f was not placed in library iplib_32
echo ' '
exit
fi
done
#
# End of script compallip_32.sh
