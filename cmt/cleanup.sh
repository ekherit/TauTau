if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/.ihep.ac.cn/bes3/offline/ExternalLib/contrib/CMT/v1r20p20081118; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=TauMass -version=0.1 -path=/afs/ihep.ac.cn/users/n/nikolaev/batch/6.6.2 $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

