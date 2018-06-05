# echo "cleanup TauTau TauTau-00-00-01 in /afs/ihep.ac.cn/users/n/nikolaev/batch/7.0.3"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtTauTautempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtTauTautempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=TauTau -version=TauTau-00-00-01 -path=/afs/ihep.ac.cn/users/n/nikolaev/batch/7.0.3  $* >${cmtTauTautempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=TauTau -version=TauTau-00-00-01 -path=/afs/ihep.ac.cn/users/n/nikolaev/batch/7.0.3  $* >${cmtTauTautempfile}"
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

