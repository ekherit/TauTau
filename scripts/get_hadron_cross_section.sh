#!/bin/bash
#run it in monte carlo with hadrons
cd /besfs/groups/tauqcd/zhangjy/bos704/subjob/PubSimRec/hadrons
for i in 3.*; do echo W=$i; cat $i/obsxs.dat | awk '{a=a+$2;}; END {print "sigma=",a}'; done;

