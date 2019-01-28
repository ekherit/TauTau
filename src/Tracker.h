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
struct Tracker
{
  typedef std::list<EvtRecTrack*> List;
  typedef std::vector<EvtRecTrack*> Vector;
  SmartDataPtr<EvtRecEvent>    evtRecEvent;
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol;
  Tracker(void) : 
    evtRecEvent(SmartDataPtr<EvtRecEvent>(eventSvc(), EventModel::EvtRec::EvtRecEvent)),
    evtRecTrkCol(SmartDataPtr<EvtRecTrackCol>(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol)

  {
  };

  template<typename Container>
  inline Container  GetCentralTracks(const double zmax=10.0, const double rhomax=1.0, const bool use_db=true)
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

  long GetNtrackCharged(void) const { return  evtRecEvent->totalCharged(); };
  long GetNtrackNeutral(void) const { return  evtRecEvent->totalNeutral(); };
  long GetNtrackTotal(void)   const { return  evtRecEvent->evtRecEvent->totalTracks(); };

  template<typename Container>
  Container GetNeutralTracks()
  {
    Container result;
    for(int i = evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
    {
      EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
      if(!(*itTrk)->isEmcShowerValid()) continue; //keep only valid neutral tracks
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      double theta = emcTrk->theta();
      double phi = emcTrk->phi();
      double c =  fabs(cos(theta)); //abs cos theta
      double E  =  emcTrk->energy();
      bool hit_barrel = (c < 0.83);
      bool hit_endcup = (0.83 <=c); //&&(c <= 0.93);
      //barrel and endcup calorimeters have different energy threshold
      bool barrel_good_track = hit_barrel && (E > 0.025);
      bool endcup_good_track = hit_endcup && (E > 0.050);
      if(barrel_good_track  || endcup_good_track) 
      {
        good_neutral_tracks.push_back(*itTrk);
      }
    }
    return good_neutral_tracks;
  }
};
