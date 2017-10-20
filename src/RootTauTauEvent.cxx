// =====================================================================================
//
//       Filename:  RootEvent.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 16:59:57
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "RootEvent.h"
#include "PhysConst.h"
#include "Utils.h"

void RootEvent::fill(int i,  EvtRecTrack * Trk)
{
  if(Trk->isMucTrackValid())
  {
    RecMucTrack *mucTrk = Trk->mucTrack();
    T.depth[i]= mucTrk->depth();
  }
  else 
  {
    T.depth[i] = - 9999;
  }
}
