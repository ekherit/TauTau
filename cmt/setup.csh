# echo "setup TauTau  in /home/nikolaev"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /home/nikolaev/BOSS/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtTauTautempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtTauTautempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=TauTau -version= -path=/home/nikolaev  -no_cleanup $* >${cmtTauTautempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=TauTau -version= -path=/home/nikolaev  -no_cleanup $* >${cmtTauTautempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtTauTautempfile}
  unset cmtTauTautempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtTauTautempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtTauTautempfile}
unset cmtTauTautempfile
exit $cmtsetupstatus

