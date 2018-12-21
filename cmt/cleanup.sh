# echo "cleanup TauTau  in /home/nikolaev"

if test "${CMTROOT}" = ""; then
  CMTROOT=/home/nikolaev/BOSS/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtTauTautempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtTauTautempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=TauTau -version= -path=/home/nikolaev  $* >${cmtTauTautempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=TauTau -version= -path=/home/nikolaev  $* >${cmtTauTautempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtTauTautempfile}
  unset cmtTauTautempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtTauTautempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtTauTautempfile}
unset cmtTauTautempfile
return $cmtcleanupstatus

