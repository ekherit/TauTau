#!/bin/tcsh 

######## set cvs
setenv CVS_RSH "/usr/bin/ssh"
setenv CVSROOT ':pserver:bes3@koala.ihep.ac.cn:/bes/bes'

set BOSS_VERSION=6.6.2
#echo '--------for boss 6.5.5_build  environments ---'
source /ihepbatch/bes/$USER/$BOSS_VERSION/cmthome/setupCMT.csh
source /ihepbatch/bes/$USER/$BOSS_VERSION/cmthome/setup.csh
source /ihepbatch/bes/$USER/$BOSS_VERSION/cmthome/setupCVS.csh
source /ihepbatch/bes/$USER/$BOSS_VERSION/TestRelease/*/cmt/setup.csh
source /ihepbatch/bes/$USER/$BOSS_VERSION/JPsi/*/cmt/setup.csh
echo "Welcome to BOSS $BOSS_VERSION"
#echo 'CMTPATH:   '$CMTPATH
#echo 'CMTROOT:   '$CMTROOT $CMTBIN
