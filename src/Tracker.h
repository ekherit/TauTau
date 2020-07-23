/*
 * =====================================================================================
 *
 *       Filename:  Tracker.h
 *
 *    Description:  Manipulation with tracks
 *
 *        Version:  1.0
 *        Created:  28.01.2019 18:12:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#pragma once

#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "Vertex.h"

#include <list>
#include <vector>
#include <limits>
#include <cassert>



struct Tracker
{
  typedef std::list<EvtRecTrack*> List;
  typedef std::vector<EvtRecTrack*> Vector;
  SmartDataPtr<EvtRecEvent>    evtRecEvent;
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol;
  Tracker(SmartDataPtr<EvtRecEvent> evt,  SmartDataPtr<EvtRecTrackCol> trk) : 
    evtRecEvent(evt), evtRecTrkCol (trk)
  {
  };

  template<typename Container>
  inline Container  GetCentralTracks(const double zmax, const double rhomax, const bool use_db)
  {
    Container result;
    for(unsigned i = 0; i < evtRecEvent->totalCharged(); i++)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
      if(!(*itTrk)->isMdcTrackValid()) continue;   //use only valid charged tracks
      RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();  //main drift chamber
      Vertex_t vtx(mdcTrk);
      if(use_db) vtx.use_db(mdcTrk);
      bool is_close_z = fabs(vtx.z) < zmax; 
      bool is_close_xy = fabs(vtx.rho) < rhomax;  
      if(is_close_z && is_close_xy) result.push_back(*itTrk);
    }
    return result;
  }



  long GetNtrackCharged(void)  { return  evtRecEvent->totalCharged(); };
  long GetNtrackNeutral(void)  { return  evtRecEvent->totalNeutral(); };
  long GetNtrackTotal(void)    { return  evtRecEvent->totalTracks(); };

  template<typename Container>
    Container GetNeutralTracks(double Emin=0)
    {
      Container result;
      for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++) {
        EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
        if(!(*itTrk)->isEmcShowerValid()) continue; //keep only valid neutral tracks
        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        if(emcTrk->energy() > Emin) result.push_back(*itTrk);
      }
      return result;
    }

  template<typename Container>
    Container GetGoodNeutralTracks(void)
    {
      Container input = GetNeutralTracks<Container>(0.025); //Get preliminary neutral tracks
      Container result;
      for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it) {
        assert((*it)->isEmcShowerValid());
        RecEmcShower *emcTrk = (*it)->emcShower();
        double c =  fabs(cos(emcTrk->theta())); //abs cos theta
        double E  =  emcTrk->energy();
        double t  =  emcTrk->time();
        double angle = 180/(CLHEP::pi) * AngleToCloseCharged(emcTrk);
        bool hit_barrel = (c < 0.8);
        bool hit_endcup = (0.86 <c) && (c < 0.92);
        //barrel and endcup calorimeters have different energy threshold
        bool barrel_good_track = hit_barrel && (E > 0.025);
        bool endcup_good_track = hit_endcup && (E > 0.050);
        //bool good_time = 0 < t && t < 14;
        bool no_close_charged = angle > 20;
        if(//good_time 
           // && 
            (barrel_good_track || endcup_good_track) 
            && 
            no_close_charged
          ) {
          result.push_back(*it);
        }
      }
      return result;
    }

  //minimum and maximum neutral tracks energy of all reconstructed tracks
  double MaxNeutralTracksEnergy(void)
  {
    double Emax = 0;
    for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
      if(!(*itTrk)->isEmcShowerValid()) continue; //keep only valid neutral tracks
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      double E = emcTrk->energy();
      if ( E > Emax )  Emax = E;
    }
    return Emax;
  };

  inline double AngleToCloseCharged(RecEmcShower * emcTrk)
  {
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
    double dthe = 200.;
    double dphi = 200.;
    double dang = 200.; 
    for( int j = 0; j< evtRecEvent->totalCharged(); j++) {
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
      if(angd < dang) {
        dang = angd;
        dthe = thed;
        dphi = phid;
      }
      if(dang>=200) continue;
    }
    return dang;
  }

  double MinNeutralTracksEnergy(void)
  {
    double Emin = std::numeric_limits<double>::max(); //unfound value
    for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      if(!(*itTrk)->isEmcShowerValid()) continue; //keep only valid neutral tracks
      double E = emcTrk->energy();
      if (E < Emin)  Emin = E;
    }
    return Emin;
  };

  double GetTotalNeutralTracksEnergy(void)
  {
    double Etotal=0;
    for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      if(!(*itTrk)->isEmcShowerValid()) continue; //keep only valid neutral tracks
      double E = emcTrk->energy();
      Etotal+=E;
    }
    return Etotal;
  }
};

template <typename Container>
double GetTotalNeutralTracksEnergy(const Container & input, double Emin=0)
{
  double Etotal=0;
  for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    assert((*it)->isEmcShowerValid());
    RecEmcShower *emcTrk = (*it)->emcShower();
    double E = emcTrk->energy();
    if ( E> Emin ) Etotal+=E;
  }
  return Etotal;
}

template <typename Container>
double MinNeutralTracksEnergy(const Container & input)
{
  double Emin = std::numeric_limits<double>::max(); //unfound value
  for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    assert((*it)->isEmcShowerValid());
    RecEmcShower *emcTrk = (*it)->emcShower();
    double E = emcTrk->energy();
    if (E < Emin)  Emin = E;
  }
  return Emin;
}

template <typename Container>
double MaxNeutralTracksEnergy(const Container & input)
{
  double Emax = -1; //unfound value
  for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    assert((*it)->isEmcShowerValid());
    RecEmcShower *emcTrk = (*it)->emcShower();
    double E = emcTrk->energy();
    if (E > Emax)  Emax = E;
  }
  return Emax;
}


template<typename Container>
std::vector<HepLorentzVector>  GetMdcLorentzVector(const Container & input, double mass = 0.5109989461e-3 /* electron mass in GeV*/)
{
  std::vector<HepLorentzVector> result;
  result.reserve(input.size());
  for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    //assert( (*it)->isMdcTrackValid(), "ERROR: GetMdcLorentzVector: No MDC information in a track");
    assert( (*it)->isMdcTrackValid());
    RecMdcKalTrack * t = (*it)->mdcKalTrack();
    result.push_back( t->p4(mass) );
  }
  return result;
};

template<typename Container>
std::vector<HepLorentzVector>  GetEmcLorentzVector(const Container & input)
{
  std::vector<HepLorentzVector> result;
  result.reserve(input.size());
  for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    //assert((*it)->isEmcShowerValid(), "ERROR: GetEmcLorentzVector: No EMC information in a track");
    assert((*it)->isEmcShowerValid());
    RecEmcShower *emcTrk = (*it)->emcShower();
    double E = emcTrk->energy();
    HepLorentzVector p(0,0,E,E);  //px,py,pz, E
    //correct direction
    p.setTheta(emcTrk->theta()); 
    p.setPhi(emcTrk->phi());
    result.push_back(p);
  }
  return result;
};

template<typename Container>
int GetTotalCharge(Container  &  input)
{
  int charge=0;
  for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    //assert( (*it)->isMdcTrackValid(), "ERROR: GetTotalCharge: No MDC information in a track");
    assert( (*it)->isMdcTrackValid());
    RecMdcTrack * mdcTrk = (*it)->mdcTrack();
    charge +=  int(mdcTrk->charge());
  }
  return charge;
}

template <typename Container >
HepLorentzVector GetTotalFourMomentum( const Container  & input )
{
  //assert(std::is_same<typename Container::value_type, HepLorentzVector>, "ERROR: GetTotalFourMomentum: Expect conatainer with HepLorentzVector");
  HepLorentzVector result(0,0,0,0);
  for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    result += *it;
  }
  return result;
};

template <typename Container >
Hep3Vector GetTotalMomentum( const Container  & input )
{
  //assert(std::is_same<typename Container::value_type, HepLorentzVector>, "ERROR: GetTotalMomentum: Expect container with HepLorentzVector");
  Hep3Vector result(0.0,0.0,0.0);
  for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    result += it->vect();
  }
  return result;
};

template <typename Container >
Hep3Vector GetTotalTransverseMomentum( const Container  & input )
{
  //assert(std::is_same<typename Container::value_type, HepLorentzVector>, "ERROR: GetTotalMomentum: Expect container with HepLorentzVector");
  Hep3Vector result = GetTotalMomentum(input);
  result.setZ(0);
  return result;
};

template<typename Container>
inline Container FilterMdcTracksByCosTheta(Container  &  input, const double MAX_COS_THETA)
{
  Container result;
  for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    //assert( (*it)->isMdcTrackValid(), "ERROR: FilterMdcTracksByCosTheta:  No MDC information in track");
    assert( (*it)->isMdcTrackValid());
    RecMdcTrack * mdcTrk = (*it)->mdcTrack();
    if(fabs(cos(mdcTrk->theta())) < MAX_COS_THETA) result.push_back(*it); 
  }
  return result;
}

template<typename Container>
inline Container FilterMdcTracksByEmcEnergy(const Container  &  input, const double MIN_EMC_ENERGY /* minimum allowed emc energy for charged tracks */)
{
  Container result;
  for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    if(!(*it)->isEmcShowerValid()) continue; //charged track must have energy deposition in EMC
    RecEmcShower *emcTrk = (*it)->emcShower();
    if ( emcTrk->energy() > MIN_EMC_ENERGY ) result.push_back(*it);
  }
  return result;
}

template<typename Container>
std::vector< Container >  SplitByCharge(const Container & input) {
  std::vector<Container> result(3);
  /* result[0] - negative charged particles
   * result[1] - positive charged particles
   * result[2] - other (no MDC or neutarls)
   */
  for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it) {
    size_t idx;
    if((*it)->isMdcTrackValid()) {
      double q =(*it)->mdcTrack()->charge();
      idx = q < 0 ? 0 : (q > 0 ? 1 : 2 );
    } else {
      idx=2;
    }
    result[idx].push_back(*it);
  }
  return result;
}

template <typename Container>
Container Zip(const Container & C1, const Container & C2, bool is_add_remains = false) {
  Container result;
  size_t npairs = std::min( C1.size(), C2.size()); //number of pairs
  for(int i=0; i != npairs; ++i) {
    result.push_back(C1[i]);
    result.push_back(C2[i]);
  }
  if(is_add_remains) { //add unpaired tracks back
      const std::vector<EvtRecTrack*> & tmp = C1.size() > C2.size() ?  C1 :  C2;
      for(int i = npairs; i < tmp.size() ; ++i)  result.push_back(tmp[i]);
  }
  return result;
}



