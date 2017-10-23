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

const int UNSET_VALUE = -999;

void RootTauTauEvent::fill(int i,  EvtRecTrack * track)
{
  if(track->isMdcTrackValid())
  {
    RecMdcKalTrack * mdc = track->mdcKalTrack();
    T.id[i] = mdc->trackId(); //id of the track
    T.q[i] =  mdc->charge(); //charge of the track
    T.p[i] = mdc->p();
    if(track->isEmcShowerValid())
    {
      T.E[i] = track->emcShower()->energy();
      T.Ep[i] = T.E[i]/T.p[i];
    }
    else
    {
      T.E[i] = UNSET_VALUE;
      T.Ep[i] = UNSET_VALUE;
    }
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
    T.id[i] = UNSET_VALUE;
    T.q[i]  = UNSET_VALUE;
    T.E[i]  = UNSET_VALUE; 
    T.p[i]  = UNSET_VALUE; 
    T.px[i] = UNSET_VALUE; 
    T.py[i] = UNSET_VALUE; 
    T.pz[i] = UNSET_VALUE; 
    T.pt[i] = UNSET_VALUE; 
    T.theta[i]= UNSET_VALUE; 
    T.phi[i] =  UNSET_VALUE;
    T.x[i] = UNSET_VALUE;
    T.y[i] = UNSET_VALUE;
    T.z[i] = UNSET_VALUE; 
    T.r[i] = UNSET_VALUE; 
    T.vxy[i] = UNSET_VALUE;
    T.vz[i] = UNSET_VALUE; 
    T.vphi[i] = UNSET_VALUE; 
  }
  if(track->isMucTrackValid())
  {
    RecMucTrack *mucTrk = track->mucTrack();
    T.depth[i]= mucTrk->depth();
    T.Nmuhit[i] = mucTrk->numHits();
  }
  else 
  {
    T.depth[i] =  UNSET_VALUE;
    T.Nmuhit[i] = UNSET_VALUE;
  }
}
