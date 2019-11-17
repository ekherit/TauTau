#!/bin/bash
#run it in monte carlo with hadrons
cd /besfs/groups/tauqcd/zhangjy/bos704/subjob/PubSimRec/lumi/bhabha
for i in 3.* 
do 
  sigma=`cat $i | head -3  | tail -1`
  echo $i  $sigma
done
