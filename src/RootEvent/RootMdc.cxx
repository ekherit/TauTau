// =====================================================================================
//
//       Filename:  RootMdc.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 17:02:36
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================
#include "RootMdc.h"

void RootMdc::init_tuple(void)
{
	//add_to_tuple(tuple);
}


void RootMdc::init(void)
{
}


void RootMdc::fill(int i, EvtRecTrack * track)
{
	if(!track->isMdcTrackValid()) return; 
	//RecMdcTrack  *mdcTrk = (*track)->mdcTrack();
  RecMdcKalTrack * mdcTrk = track->mdcKalTrack();
	id[i] = mdcTrk->trackId();
	q[i] = mdcTrk->charge(); 
	p[i] = mdcTrk->p();
	px[i]= mdcTrk->px();
	py[i]= mdcTrk->py();
	pz[i]= mdcTrk->pz();
	theta[i]= mdcTrk->theta();
	phi[i] = mdcTrk->phi();
	x[i]  = mdcTrk->x();
	y[i]  = mdcTrk->y();
	z[i]  = mdcTrk->z();
	x[i]  = mdcTrk->x();
	y[i]  = mdcTrk->y();
	z[i]  = mdcTrk->z();
	double rvxy=0,rvz=0,rvphi=0;
	//calculate_vertex(track->mdcTrack(),rvxy,rvz,rvphi); 
	vxy[i] = rvxy;
	vz[i]  = rvz; 
	vphi[i] = rvphi; 
	if(track->isEmcShowerValid())
	{
		RecEmcShower *emcTrk = track->emcShower();
		E[i] = emcTrk->energy();
	}
	else
	{
		E[i] = 0;
	}
  if(track->isMucTrackValid())
  {
    RecMucTrack *mucTrk = track->mucTrack();
    depth[i]= mucTrk->depth();
  }
  else 
  {
    depth[i] = - 9999;
  }
}
