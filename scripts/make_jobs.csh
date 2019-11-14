#!/bin/bash
CDIR=`pwd`
WDIR=$HOME/batch4/TAU2018
cd $WDIR
TDIR=$WDIR/$1
cd $TDIR
mkdir data
$TAUTAUROOT/scripts/make_job_option.py $WDIR/data --config=$TAUTAUROOT/share/all_scan_points_ems3.txt 
cd ../
mkdir mc
cd mc
mkdir bhabha
cd bhabha
$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/lumi/bhabha
