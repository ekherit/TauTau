// =====================================================================================
//
//       Filename:  Utils.h
//
//    Description:  Some usefull functions
//
//        Version:  1.0
//        Created:  19.10.2015 21:12:17
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================


#pragma once

#include <exception>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

//#include "GaudiKernel/IDataProviderSvc.h"
//#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Bootstrap.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/Helix.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

//#include "utils.h"
#include "TypeDefs.h"
#include "SelectionConfig.h"
#include "PhysConst.h"


inline HepLorentzVector getTotalMomentum(double Wcm = BEAM_CENTER_MASS_ENERGY)
{
	return HepLorentzVector(Wcm*tan(0.5*BEPC_CROSSING_ANGLE),0,0,Wcm/cos(0.5*BEPC_CROSSING_ANGLE));
	//I think formula below is wrong 
	//return HepLorentzVector(Wcm*sin(0.5*BEPC_CROSSING_ANGLE),0,0,Wcm);
}

inline void calculate_vertex(RecMdcTrack *mdcTrk, double & ro, double  & z, double phi)
{
  ro = -9999;
  z = -9999;
  phi = -9999;
  /*  Reconstruct the vertex */
  Hep3Vector xorigin(0,0,0);
  IVertexDbSvc*  vtxsvc;
  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if(vtxsvc->isVertexValid())
  {
    double* dbv = vtxsvc->PrimaryVertex(); 
    double*  vv = vtxsvc->SigmaPrimaryVertex();  
    xorigin.setX(dbv[0]);
    xorigin.setY(dbv[1]);
    xorigin.setZ(dbv[2]);
  }
  /* Vertex game. copy from rhophi analysis */
  double phi0=mdcTrk->helix(1);
  double xv=xorigin.x();
  double yv=xorigin.y();
  //double Rxy=(mdc.x[i]-xv)*cos(phi0)+(mdc.y[i]-yv)*sin(phi0);
  //mdc.r[i]=Rxy;
  HepVector a = mdcTrk->helix();
  HepSymMatrix Ea = mdcTrk->err();
  HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
  HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
  VFHelix helixip(point0,a,Ea); 
  helixip.pivot(IP);
  HepVector vecipa = helixip.a();
  double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
  double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
  double  Rvphi0=vecipa[1];
  ro=Rvxy0;
  z=Rvz0;
  phi=Rvphi0;
}


inline std::list<EvtRecTrack*> createGoodChargedTrackList(
		SelectionConfig & cfg, 
		SmartDataPtr<EvtRecEvent>    & evtRecEvent, 
		SmartDataPtr<EvtRecTrackCol> & evtRecTrkCol
		)
{
  std::list<EvtRecTrack*> good_charged_tracks;
  for(unsigned i = 0; i < evtRecEvent->totalCharged(); i++)
  {
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isMdcTrackValid()) continue;  //use only valid charged tracks
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();  //main drift chambe
    //calculate interaction point distance
    double rvxy,rvz,rvphi;
    calculate_vertex(mdcTrk,rvxy,rvz,rvphi); //find distance to interaction point
    bool IP_track = fabs(rvz)< cfg.IP_MAX_Z && fabs(rvxy)<cfg.IP_MAX_RHO;  //tracks begin near interaction point
    bool good_track = IP_track && fabs(cos(mdcTrk->theta()))<cfg.MAX_COS_THETA; //track is good
    if(good_track) good_charged_tracks.push_back(*itTrk);
  }
	return good_charged_tracks;
}

inline std::list<EvtRecTrack*> createGoodEmcChargedTrackList(
		SelectionConfig & cfg, 
    const std::list<EvtRecTrack*> & good_charged_tracks
		)
{
  std::list<EvtRecTrack*> emc_good_charged_tracks;
  for( std::list<EvtRecTrack*>::iterator it = good_charged_tracks.begin();
      it!=good_charged_tracks.end();
      ++it)
  {
    EvtRecTrack * track = *it;
    if(!track->isMdcTrackValid()) continue;  //use only valid charged tracks
    if(!track->isEmcShowerValid()) continue; //charged track must have energy deposition in EMC
    RecMdcTrack * mdcTrk = track->mdcTrack();  //main drift chambe
    RecEmcShower *emcTrk = track->emcShower(); //Electro Magnet Calorimeer
    emc_good_charged_tracks.push_back(track);
  }
  return emc_good_charged_tracks;
}

inline std::list<EvtRecTrack*> createGoodNeutralTrackList(
		SelectionConfig & cfg, 
		SmartDataPtr<EvtRecEvent>    & evtRecEvent, 
		SmartDataPtr<EvtRecTrackCol> & evtRecTrkCol
		)
{
	std::list<EvtRecTrack*> good_neutral_tracks;
	//collect good neutral track
	for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isEmcShowerValid()) continue; //keep only valid neutral tracks
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		double c =  fabs(cos(emcTrk->theta())); //abs cos theta
		double E  =  emcTrk->energy();
		bool hit_barrel = (c <= cfg.EMC_BARREL_MAX_COS_THETA);
		bool hit_endcup = (cfg.EMC_ENDCUP_MIN_COS_THETA <=c) && (c <= cfg.EMC_ENDCUP_MAX_COS_THETA);
		//barrel and endcup calorimeters have different energy threshold
		bool barrel_good_track = hit_barrel && (E > cfg.EMC_BARREL_MIN_ENERGY);
		bool endcup_good_track = hit_endcup && (E > cfg.EMC_ENDCUP_MIN_ENERGY);
		if(barrel_good_track  || endcup_good_track) 
		{
			//cout << "Energy of good neutral track: " << E << endl;
			good_neutral_tracks.push_back(*itTrk);
		}
	}
	return good_neutral_tracks;
}


inline double angle_to_close_charged(
    RecEmcShower * emcTrk,
		SmartDataPtr<EvtRecEvent>    & evtRecEvent, 
		SmartDataPtr<EvtRecTrackCol> & evtRecTrkCol
    )
{
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
    double dthe = 200.;
    double dphi = 200.;
    double dang = 200.; 
    for( int j = 0; j< evtRecEvent->totalCharged(); j++)
    {
      EvtRecTrackIterator itChargedTrk=evtRecTrkCol->begin() + j;
      if(!(*itChargedTrk)->isExtTrackValid()) continue;
      RecExtTrack *extTrk = (*itChargedTrk)->extTrack();
      if(extTrk->emcVolumeNumber() == -1) continue;
      Hep3Vector extpos = extTrk->emcPosition();
      double angd = extpos.angle(emcpos);
      double thed = extpos.theta() - emcpos.theta();
      double phid = extpos.deltaPhi(emcpos);
      thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      if(angd < dang)
      {
        dang = angd;
        dthe = thed;
        dphi = phid;
      }
      if(dang>=200) continue;
    }
    return dang;
}


inline std::list<EvtRecTrack*> createGoodNeutralTrackList2(
		SelectionConfig & cfg, 
		SmartDataPtr<EvtRecEvent>    & evtRecEvent, 
		SmartDataPtr<EvtRecTrackCol> & evtRecTrkCol
		)
{
	std::list<EvtRecTrack*> good_neutral_tracks;
	//collect good neutral track
	for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isEmcShowerValid()) continue; //keep only valid neutral tracks
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
    double theta = emcTrk->theta();
    double phi = emcTrk->phi();
		double c =  fabs(cos(theta)); //abs cos theta
		double E  =  emcTrk->energy();
		bool hit_barrel = (c <= cfg.EMC_BARREL_MAX_COS_THETA);
		bool hit_endcup = (cfg.EMC_ENDCUP_MIN_COS_THETA <=c) && (c <= cfg.EMC_ENDCUP_MAX_COS_THETA);
		//barrel and endcup calorimeters have different energy threshold
		bool barrel_good_track = hit_barrel && (E > cfg.EMC_BARREL_MIN_ENERGY);
		bool endcup_good_track = hit_endcup && (E > cfg.EMC_ENDCUP_MIN_ENERGY);

    double angle_to_charged  =  180/(CLHEP::pi) * angle_to_close_charged(emcTrk,evtRecEvent, evtRecTrkCol);
    bool close_charged_track = fabs(angle_to_charged) < cfg.NEUTRAL_CLOSE_CHARGED_ANGLE;

    if(close_charged_track) continue;

		if(barrel_good_track  || endcup_good_track) 
		{
			good_neutral_tracks.push_back(*itTrk);
		}
	}
	return good_neutral_tracks;
}


/*  
 *
 *  E1 < E2  true if E1<E2
 *  none < E2 always true
 *  E1 < none always false
 *  none < none always false
 *   */
inline bool EmcEnergyOrder(EvtRecTrack * track1, EvtRecTrack *track2)
{
    if(!track2->isEmcShowerValid())  return false;
    RecEmcShower *emcTrk1 = track2->emcShower();
    if(!track1->isEmcShowerValid())  return true;
    RecEmcShower *emcTrk2 = track1->emcShower();
    return emcTrk1->energy() < emcTrk2->energy();
}

inline bool ChargeOrder(EvtRecTrack * track1, EvtRecTrack *track2)
{
    if(!track2->isMdcTrackValid()) return false; 
    RecMdcTrack * mdc2 = track2->mdcTrack();  
    if(!track1->isMdcTrackValid()) return true; 
    RecMdcTrack * mdc1 = track1->mdcTrack();  
    return mdc1->charge() < mdc2->charge();
}

inline int GetCharge(EvtRecTrack * track)
{
  if(!track->isMdcTrackValid()) return -999;
  RecMdcTrack * mdc = track->mdcTrack();  
  return mdc->charge();
}

inline double GetMomentum(EvtRecTrack * track)
{
  if(!track->isMdcTrackValid()) return -999;
  RecMdcTrack * mdc = track->mdcTrack();  
  return mdc->p();
}

inline Hep3Vector GetHep3Vector(EvtRecTrack * track)
{
  if(!track->isMdcTrackValid()) return Hep3Vector();
  return Hep3Vector(track->mdcKalTrack()->px(), track->mdcKalTrack()->py(),track->mdcKalTrack()->pz());
}


inline double Acoplanarity(double phi1, double phi2)
{
    double dphi = (phi2-phi1);
    if(dphi > +M_PI) dphi = dphi - 2*M_PI;
    if(dphi < -M_PI) dphi = dphi + 2*M_PI;
    return  M_PI - fabs(dphi);
  //return fmod(phi2-phi1 + 4*M_PI, 2*M_PI) - M_PI;
}
