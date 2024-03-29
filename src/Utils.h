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
#include <functional>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include <TMatrixDEigen.h>

#include "TypeDefs.h"
#include "SelectionConfig.h"
#include "PhysConst.h"

#include "Vertex.h"

template < class T> inline T sq(T value) { return value*value; }

inline HepLorentzVector getTotalMomentum(double Wcm = BEAM_CENTER_MASS_ENERGY)
{
  double a_2 = 0.5*BEPC_CROSSING_ANGLE;
	return HepLorentzVector(Wcm*tan(a_2),0,0,Wcm/cos(a_2)); //check this formula 2019-01-25
}

inline HepLorentzVector getInitialFourMomentum(double Wcm = BEAM_CENTER_MASS_ENERGY)
{
  double a_2 = 0.5*BEPC_CROSSING_ANGLE;
	return HepLorentzVector(Wcm*tan(a_2),0,0,Wcm/cos(a_2)); //check this formula 2019-01-25
}


/*
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
    Vertex_t vtx(mdcTrk);
    //calculate interaction point distance
    //double rvxy,rvz,rvphi;
    //calculate_vertex(mdcTrk,rvxy,rvz,rvphi); //find distance to interaction point
    bool IP_track = fabs(vtx.z)< cfg.IP_MAX_Z && fabs(vtx.rho)<cfg.IP_MAX_RHO;  //tracks begin near interaction point
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
  for( std::list<EvtRecTrack*>::const_iterator it = good_charged_tracks.begin();
      it!=good_charged_tracks.end();
      ++it)
  {
    EvtRecTrack * track = *it;
    //if(!track->isMdcTrackValid()) continue;  //use only valid charged tracks
    if(!track->isEmcShowerValid()) continue; //charged track must have energy deposition in EMC
    //RecMdcTrack * mdcTrk = track->mdcTrack();  //main drift chambe
    //RecEmcShower *emcTrk = track->emcShower(); //Electro Magnet Calorimeer
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

*/


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
    RecEmcShower *emc2 = track2->emcShower();
    if(!track1->isEmcShowerValid())  return true;
    RecEmcShower *emc1 = track1->emcShower();
    return emc1->energy() < emc2->energy();
}

inline bool ChargeOrder(EvtRecTrack * track1, EvtRecTrack *track2)
{
    if(!track2->isMdcTrackValid()) return false; 
    RecMdcTrack * mdc2 = track2->mdcTrack();  
    if(!track1->isMdcTrackValid()) return true; 
    RecMdcTrack * mdc1 = track1->mdcTrack();  
    return mdc1->charge() < mdc2->charge();
}

inline bool PtOrder(EvtRecTrack * track1, EvtRecTrack *track2)
{
    if(!track2->isMdcTrackValid()) return false; 
    RecMdcTrack * mdc2 = track2->mdcTrack();  
    if(!track1->isMdcTrackValid()) return true; 
    RecMdcTrack * mdc1 = track1->mdcTrack();  
    return mdc1->pxy() < mdc2->pxy();
}

inline bool MomentumOrder(EvtRecTrack * track1, EvtRecTrack *track2)
{
    if(!track2->isMdcTrackValid()) return false; 
    RecMdcTrack * mdc2 = track2->mdcTrack();  
    if(!track1->isMdcTrackValid()) return true; 
    RecMdcTrack * mdc1 = track1->mdcTrack();  
    return mdc1->p() < mdc2->p();
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

inline double Acoplanarity(EvtRecTrack * t1, EvtRecTrack *t2)
{
  RecMdcKalTrack * mdcKal1 = t1->mdcKalTrack();
  RecMdcKalTrack * mdcKal2 = t2->mdcKalTrack();
  return Acoplanarity( mdcKal1->phi(), mdcKal2->phi());
}

inline double Acolinearity(EvtRecTrack * t1, EvtRecTrack *t2)
{
  RecMdcKalTrack * mdcKal1 = t1->mdcKalTrack();
  RecMdcKalTrack * mdcKal2 = t2->mdcKalTrack();
  return mdcKal1->p3().angle(mdcKal2->p3());
}



inline std::vector<double> getSphericityEigenvalues(const std::vector<EvtRecTrack*> & T)
{
  TMatrixD S(3,3); //sphericity tensor
  //I am not sure that values will be zeroes after creation
  for(int i=0;i<3;++i) for(int j=0;j<3;++j) S[i][j] =  0;
  double p2sum = 0;
  for(int track = 0;track < T.size(); ++track)
  {
    RecMdcKalTrack * mdcKal = T[track]->mdcKalTrack();
    p2sum += sq(mdcKal->p());
    for(int i=0;i<3;++i) for(int j=0;j<3;++j) 
    {
      S[i][j] += mdcKal->p3()[i] * mdcKal->p3()[j];
    }
  }
  //normalize sphericity tenzor
  for(int i=0;i<3;++i) for(int j=0;j<3;++j) S[i][j]/=p2sum;
  //calculate eigen values
  TMatrixDEigen Seigen(S);
  const TVectorD & v = Seigen.GetEigenValuesRe();
  //copy result to my vector
  std::vector<double> result(3);
  for(int i=0;i<3;i++) result[i]=v[i];
  //sort eigenvalues in descent order
  std::sort(result.rbegin(), result.rend());
  return result;
}

inline double Sphericity(std::vector<double> & v)
{
  return 1.5*(v[1]+v[2]);
}

inline double Aplanarity(std::vector<double> & v)
{
  return 1.5*v[2];
}


template < class T>
bool in( const T & value, const T & min, const T & max ) { return  min < value && value < max; }

/* Have to write own accumulator - just copy from cppreference.com */
/*
template<class InputIt, class T, class UnaryOperation, class BinaryOperation>
T accumulate(InputIt first, InputIt last, const UnaryOperation & F, T init, BinaryOperation op=std::plus)
{
    for (; first != last; ++first)
    {
        init = op(init, F(*first)); 
    }
    return init;
}

struct  ShiftedInvariantMassOperator
{
  InvariantMassOperator(double s=0, double e) : shift(s), error(e) {} 
  double operator(const HepLorentzVector & v1, const HepLorentzVector &v2)
  {
    return (v1+v2).mag()-shift;
  }

  double operator(const std::pair<HepLorentzVector, HepLorentzVector> & p)
  {
    return this->(pair.first, pair.second);
  }

  double operator(const std::pair<HepLorentzVector*, HepLorentzVector*> & p)
  {
    return this->(*(pair.first), *(pair.second));
  }
};


template< class PairContainer>
inline double GetPairChi2(const  PairContainers &   pairs, double data_expected=1, double data_error = 1, )
{
  return accumulate(pairs.begin(),pairs.end(), ShiftedInvariantMassOperator(Mexp));
};
*/


//std::vector<std::vector<EvtRecTrack*> > SplitByCharge(const  std::vector<EvtRecTrack*> & Tc) {
//  std::vector<EvtRecTrack*> Tcn, Tcp;// negative and positive charged tracks
//  std::vector<std::vector<EvtRecTrack*> >  T;
//  T.resize(2);
//  T[0].reserve(Tc.size());
//  T[1].reserve(Tc.size());
//  for(auto track : Tc) {
//    int charge = GetCharge(track);
//    if( charge == -1) T[0].push_back(track);
//    if( charge == +1 ) T[1].push_back(track);
//  }
//  //Sort by transverse momentum desccending
//  std::sort(T[0].rbegin(), T[0].rend(),PtOrder);
//  std::sort(T[1].rbegin(), T[1].rend(),PtOrder);
//  return T;
//};

/*

std::vector<EvtRecTrack * > JoinTracks(const std::vector< std::vector<EvtRecTrack*> & T) {
  std::vector<EvtRecTrack * >  R; 
  bool not_empty=true;
  for(size_t idx=0; not_empty; ++idx) {
    not_empty = false;
    for(size_t i =0; i!=T.size(); ++i) {
      if(idx < T[i].size()) {
        not_empty = true;
        R.push_back(T[i]);
      }
    }
  }
  return R;
};
*/
