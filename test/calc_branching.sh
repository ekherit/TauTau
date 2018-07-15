#!/bin/bash
cat tau.dec | grep -v '^#' | grep -v -i decay | awk 'BEGIN {a=0; OFMT ="%10.7f     "}{a=a+$1;n=n+1;print a, $1, $0} END { OFMT ="%10.8f"; print "Sum of ", n, " channels = ", a}'
echo -e "\n\n"
echo "Sorted branching:"
cat tau.dec | grep -v '^#' | grep -v -i decay | sort -nr | 
  awk '
  BEGIN { 
    print "---------------------------------------------------------------------------";
  } 
  {
     a=a+$1; 
     n=n+1; 
     print $0
  } 
  END { 
    print "---------------------------------------------------------------------------"; 
    OFMT ="%10.8f"; 
    print "Sum of ", n, " channels = ", a
  }'

