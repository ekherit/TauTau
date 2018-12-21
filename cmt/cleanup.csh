# echo "cleanup TauTau  in /home/nikolaev"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /home/nikolaev/BOSS/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtTauTautempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtTauTautempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=TauTau -version= -path=/home/nikolaev  $* >${cmtTauTautempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=TauTau -version= -path=/home/nikolaev  $* >${cmtTauTautempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtTauTautempfile}
  unset cmtTauTautempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtTauTautempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtTauTautempfile}
unset cmtTauTautempfile
exit $cmtcleanupstatus

