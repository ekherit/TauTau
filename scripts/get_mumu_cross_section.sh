#!/bin/bash
cd /besfs/groups/tauqcd/zhangjy/bos704/subjob/PubSimRec/bkgd/mumu
for i in 3.* 
do 
  sigma="`cat $i/CrossSection.txt | tail -2 | head -1`"
  echo $i "      " $sigma
done
