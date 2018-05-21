#!/bin/ksh 
########################################
# This script calculates the Julian Day
########################################
set -x

date=$1

if [ $# -eq 0 ]; then
  julian=`date +%j`
else
  let "yy=$date/10000"
  let "mm=($date-$yy*10000)/100"
  let "dd=$date%100"
  case $mm in 
   1) julian=$dd;;
   2) julian=`expr $dd + 31`;;
   3) julian=`expr $dd + 31 + 28`;;
   4) julian=`expr $dd + 31 + 28 + 31`;;
   5) julian=`expr $dd + 31 + 28 + 31 + 30`;;
   6) julian=`expr $dd + 31 + 28 + 31 + 30 + 31`;;
   7) julian=`expr $dd + 31 + 28 + 31 + 30 + 31 + 30`;;
   8) julian=`expr $dd + 31 + 28 + 31 + 30 + 31 + 30 + 31`;;
   9) julian=`expr $dd + 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31`;;
  10) julian=`expr $dd + 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30`;;
  11) julian=`expr $dd + 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31`;;
  12) julian=`expr $dd + 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30`;;
  esac
  
  leap=`expr $yy % 4`
  if [ $leap -eq 0 -a $mm -gt 2 ]; then
     julian=`expr $julian + 1`
  fi
fi
echo $julian
