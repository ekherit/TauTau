#!/bin/bash
CDIR=`pwd`
WDIR=$HOME/batch4/TAU2018
cd $WDIR
TDIR=$WDIR/$1
mkdir $TDIR
cd $TDIR

#data
$TAUTAUROOT/scripts/make_job_option.py $WDIR/data data --config=$TAUTAUROOT/share/all_scan_points_ems3.txt 

#monte carlo
mkdir mc

$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/sigmc       mc/signal --W=1 --combine='\d\.\d+'

$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/lumi/bhabha     mc/bb --W=1 --combine='\d\.\d+'
$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/lumi/digam      mc/gg --W=1 --combine='\d\.\d+'
$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/hadrons    mc/hadrons --W=1 --combine='\d\.\d+'
$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/bkgd/mumu     mc/uu   --W=1 --combine='\d\.\d+'
$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/bkgd/pipi     mc/pipi --W=1 --combine='\d\.\d+'
mkdir mc/galuga
#$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/bkgd/galuga/EEee        mc/galuga/ee   --W=1 --combine='\d\.\d+'
$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/bkgd/galuga/EEee     mc/galuga/ee   --W=1 --combine='\d\.\d+'
$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/bkgd/galuga/EEkk     mc/galuga/KK   --W=1 --combine='\d\.\d+'
$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/bkgd/galuga/EEpipi   mc/galuga/pipi --W=1 --combine='\d\.\d+'
$TAUTAUROOT/scripts/make_job_option.py /besfs/groups/tauqcd/zhangjy/bos704/MC/bkgd/galuga/EEuu     mc/galuga/uu   --W=1 --combine='\d\.\d+'

find -name "submit.sh" -exec echo {} \;  > submit_all.sh
find -name "submit.csh" -exec echo {} \; > submit_all.csh
chmod u+x submit_all.sh
chmod u+x submit_all.csh
