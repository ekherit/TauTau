# echo "setup TauTau TauTau-00-00-01 in /afs/ihep.ac.cn/users/n/nikolaev/batch/7.0.3"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtTauTautempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtTauTautempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=TauTau -version=TauTau-00-00-01 -path=/afs/ihep.ac.cn/users/n/nikolaev/batch/7.0.3  -no_cleanup $* >${cmtTauTautempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=TauTau -version=TauTau-00-00-01 -path=/afs/ihep.ac.cn/users/n/nikolaev/batch/7.0.3  -no_cleanup $* >${cmtTauTautempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtTauTautempfile}
  unset cmtTauTautempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtTauTautempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtTauTautempfile}
unset cmtTauTautempfile
return $cmtsetupstatus

export LD_LIBRARY_PATH=/ihepbatch/bes/nikolaev/7.0.3/TauTau/TauTau-00-00-01/x86_64-slc6-gcc46-opt:$LD_LIBRARY_PATH
