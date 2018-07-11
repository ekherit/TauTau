#!/bin/bash
cat tau.dec | grep -v '#' | grep -v -i decay | awk 'BEGIN {a=0; OFMT ="%10.7f     "} {a=a+$1; print a, $1, $0}'
