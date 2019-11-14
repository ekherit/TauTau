setenv WDIR  $HOME/batch4/TAU2018
setenv TDIR  $WDIR/$1
mkdir $TDIR
$TAUTAUROOT/scripts/make_job_option.py $WDIR/data --config=$TAUTAUROOT/share/all_scan_points_ems3.txt $TDIR/data
