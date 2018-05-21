#!/bin/sh
# Start of script compip_64.sh          M.Farley  98-03-26
#
echo Compile and place in library iplib_64 subroutine $1.f
#
# Example:  Compile subroutine q9yi64.f and put it into library iplib_64
#
#           compip_64.sh q9yi64
#
rm    $1.o               # delete subroutine old object
# This Compile options will generate a 64-Bit Object code
if f90 -c -O -64 -r8 -i8 $1.f # compile subroutine
then
ar r ../iplib_64 $1.o
rm    $1.o               # delete subroutine new object
# rm    $1.l               # delete listing file
#
echo Subroutine $1.f was compiled and placed in library iplib_64
else
echo Compile error, subroutine $1.f was not placed in library iplib_64
fi
#
# End of script compip_64.sh
