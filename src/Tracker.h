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
//#include <type_traits>
//
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


  long GetNtrackCharged(void)  { return  evtRecEvent->totalCharged(); };
  long GetNtrackNeutral(void)  { return  evtRecEvent->totalNeutral(); };
  long GetNtrackTotal(void)    { return  evtRecEvent->totalTracks(); };

  template<typename Container>
  Container GetNeutralTracks(double Emin)
  {
    Container result;
    for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
      if(!(*itTrk)->isEmcShowerValid()) continue; //keep only valid neutral tracks
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      if(emcTrk->energy() > Emin) result.push_back(*itTrk);
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

  template <typename Container>
    double GetTotalNeutralTracksEnergy(const Container & input, double Emin=0)
    {
      double Etotal=0;
      for(typename Container::const_iterator it = input.begin(); it!=input.end(); ++it)
      {
        //EvtRecTrack * track = *it;
        //assert((*it)->isEmcShowerValid(), "ERROR: GetTotalNeutralTracksEnergy: No EMC information for a track");
        assert((*it)->isEmcShowerValid());
        RecEmcShower *emcTrk = (*it)->emcShower();
        double E = emcTrk->energy();
        if ( E> Emin ) Etotal+=E;
      }
      return Etotal;
    }
};



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
