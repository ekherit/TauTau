# echo "setup TauTau TauTau-00-00-01 in /afs/ihep.ac.cn/users/n/nikolaev/batch/7.0.3"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtTauTautempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtTauTautempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=TauTau -version=TauTau-00-00-01 -path=/afs/ihep.ac.cn/users/n/nikolaev/batch/7.0.3  -no_cleanup $* >${cmtTauTautempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=TauTau -version=TauTau-00-00-01 -path=/afs/ihep.ac.cn/users/n/nikolaev/batch/7.0.3  -no_cleanup $* >${cmtTauTautempfile}"
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

setenv LD_LIBRARY_PATH /ihepbatch/bes/nikolaev/7.0.3/TauTau/TauTau-00-00-01/x86_64-slc6-gcc46-opt:$LD_LIBRARY_PATH
