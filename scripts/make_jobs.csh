#!/bin/bash
WDIR=$HOME/batch4/TAU2018
cd $WDIR
mkdir $TDIR
$TAUTAUROOT/scripts/make_job_option.py $WDIR/data --config=$TAUTAUROOT/share/all_scan_points_ems3.txt $1/data
