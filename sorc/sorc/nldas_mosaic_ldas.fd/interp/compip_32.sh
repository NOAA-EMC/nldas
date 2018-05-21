#!/bin/sh
# Start of script compip_32.sh          M.Farley  98-03-26
#
echo Compile and place in library iplib_32 subroutine $1.f
#
# Example:  Compile subroutine q9yi32.f and put it into library iplib_32 
#
#           compip_32.sh q9yi32
#
rm    $1.o               # delete subroutine old object
# This Compile options will generate a 32-Bit Object code
if f90 -c -O -n32 $1.f # compile subroutine
then
ar r ../iplib_32 $1.o
rm    $1.o               # delete subroutine new object
# rm    $1.l               # delete listing file
#
echo Subroutine $1.f was compiled and placed in library iplib_32
else
echo Compile error, subroutine $1.f was not placed in library iplib_32
fi
#
# End of script compip_32.sh
