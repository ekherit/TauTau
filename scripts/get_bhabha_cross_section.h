#!/bin/bash
#run it in monte carlo with hadrons
cd /besfs/groups/tauqcd/zhangjy/bos704/subjob/PubSimRec/lumi/bhabha
for i in 3.* 
do 
  sigma="`cat $i/CrossSection.txt | tail -2 | head -1`"
  echo $i   sigma
done
