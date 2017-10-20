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

#include "RootTauTauEvent.h"
#include "PhysConst.h"
#include "Utils.h"

RootTauTauEvent::~RootTauTauEvent(void)
{
}

void RootTauTauEvent::fill(int i,  EvtRecTrack * track)
{
  if(track->isMdcTrackValid())
  {
    RecMdcKalTrack * mdc = track->mdcKalTrack();
    T.id[i] = mdc->trackId(); //id of the track
    T.q[i] =  mdc->charge(); //charge of the track
    //T.E[i] = mdc->energy();
    T.E[i] = 0;
    T.p[i] = mdc->p();
    std::cout << " i = " << i << " p = " << mdc->p() << std::endl;
    T.px[i] = mdc->px();
    T.py[i] = mdc->py();
    T.pz[i] = mdc->pz();
    T.pt[i] = mdc->pxy();
    T.theta[i]= mdc->theta();
    T.phi[i] = mdc->phi();
    T.x[i] = mdc->x();
    T.y[i] = mdc->y();
    T.z[i] = mdc->z(); 
    T.r[i] = 0; 

    double rho,z,phi;
    calculate_vertex(track->mdcTrack(),rho,z,phi);
    T.vxy[i] = rho;
    T.vz[i] = z; 
    T.vphi[i] = phi; 
  }
  else
  {
    T.id[i] = -999;
    T.q[i]  = -999;
    T.E[i]  = -999; 
    T.p[i]  = -999; 
    T.px[i] = -999; 
    T.py[i] = -999; 
    T.pz[i] = -999; 
    T.pt[i] = -999; 
    T.theta[i]= - 999; 
    T.phi[i] =  -999;
    T.x[i] = -999;
    T.y[i] = -999;
    T.z[i] = -999; 
    T.r[i] = -999; 

    T.vxy[i] = -999;
    T.vz[i] = -999; 
    T.vphi[i] = -999; 
  }
  if(track->isMucTrackValid())
  {
    RecMucTrack *mucTrk = track->mucTrack();
    T.depth[i]= mucTrk->depth();
  }
  else 
  {
    T.depth[i] = - 9999;
  }
}
