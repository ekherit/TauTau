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

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "Vertex.h"

#include <list>
#include <vector>
struct Tracker
{
  typedef std::list<EvtRecTrack*> List_t;
  typedef std::vector<EvtRecTrack*> Vector_t;
  SmartDataPtr<EvtRecEvent>    evtRecEvent;
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol;
  Tracker(void) : 
    evtRecEvent(SmartDataPtr<EvtRecEvent>(eventSvc(), EventModel::EvtRec::EvtRecEvent)),
    evtRecTrkCol(SmartDataPtr<EvtRecTrackCol>(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol)

  {
  };

  template<typename Container>
  inline Container  get_centeral_tracks(const double zmax=10.0, const double rhomax=1.0, const bool use_db=true)
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

};
