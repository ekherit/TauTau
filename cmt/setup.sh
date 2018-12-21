# echo "setup TauTau  in /home/nikolaev"

if test "${CMTROOT}" = ""; then
  CMTROOT=/home/nikolaev/BOSS/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtTauTautempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtTauTautempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=TauTau -version= -path=/home/nikolaev  -no_cleanup $* >${cmtTauTautempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=TauTau -version= -path=/home/nikolaev  -no_cleanup $* >${cmtTauTautempfile}"
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

