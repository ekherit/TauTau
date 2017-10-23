/*
 * =====================================================================================
 *
 *       Filename:  TauTau.cxx
 *
 *    Description:  Multihadron event selection for j/psi and psi prime resonance.
 *
 *        Version:  1.0
 *        Created:  04/27/2010 02:50:11 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"

#include "McTruth/McParticle.h"

#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include "TauTau.h"

#include "Utils.h"


#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"

inline double sq(double x) { return x*x; }

TauTau::TauTau(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator)
{
  declareProperty("CENTER_MASS_ENERGY", cfg.CENTER_MASS_ENERGY=1.777*2); //GeV
  declareProperty("MIN_CHARGED_TRACKS", cfg.MIN_CHARGED_TRACKS=2); 
  declareProperty("MAX_CHARGED_TRACKS", cfg.MAX_CHARGED_TRACKS=2); 
  //good charged track configuration
  declareProperty("IP_MAX_Z",      cfg.IP_MAX_Z = 10.0); //cm
  declareProperty("IP_MAX_RHO",    cfg.IP_MAX_RHO = 1.0); //cm
  declareProperty("MAX_COS_THETA", cfg.MAX_COS_THETA = 0.8);
}


StatusCode TauTau::initialize(void)
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in initialize()" << endmsg;
  nproceed_events=0;
  nwrited_events = 0;
  ntautau_events=0;
  nbhabha_events=0;
  ngg_events=0;
  StatusCode status;
  try
  {
    fEvent.make_tuple(this, "FILE1/event","Signal tau tau events");
  }
	catch(std::runtime_error & error)
	{
		log << MSG::ERROR << error.what() << endmsg;
		return StatusCode::FAILURE;
	}
  return StatusCode::SUCCESS;
}

StatusCode TauTau::execute()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  int runNo=eventHeader->runNumber();
  int event=eventHeader->eventNumber();
  time_t t=eventHeader->time();
  bool isprint=false;
  if(nproceed_events<10) isprint=true;
  if(10 <= nproceed_events && nproceed_events < 100 && nproceed_events % 10 ==0) isprint=true;
  if(100 <= nproceed_events && nproceed_events < 1000 && nproceed_events % 100 ==0) isprint = true;
  if(1000 <= nproceed_events && nproceed_events < 10000 && nproceed_events % 1000 ==0) isprint = true;
  if(10000 <= nproceed_events && nproceed_events % 10000 ==0) isprint = true;
  if(isprint)
  {
    std::cout << "proceed event: " << setw(15) << nproceed_events;
    std::cout << "  tau-mu:" << setw(15) << ntautau_events;
    std::cout << "  ee:" << setw(15) << nbhabha_events;
    std::cout << "  gg:" << setw(15) << ngg_events;
    std::cout << std::endl;
  }
  nproceed_events++;

  //  Get information about reconstructed events
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);

	std::list<EvtRecTrack*> good_neutral_tracks;
	std::list<EvtRecTrack*> good_charged_tracks;

  //fill initial value of the selected event
  //come from IP and cos(theta)<max_cos_theta
  good_charged_tracks=createGoodChargedTrackList(cfg, evtRecEvent, evtRecTrkCol);
  //good neutral tracks
  good_neutral_tracks=createGoodNeutralTrackList2(cfg, evtRecEvent, evtRecTrkCol);

  //std::cout << "Number of good charged tracks = " << good_charged_tracks.size() << std::endl;
  //std::cout << "Number of good neutral tracks = " << good_neutral_tracks.size() << std::endl;

  //filter good charged tracks keep only tracks with emcShower
  std::list<EvtRecTrack*> emc_good_charged_tracks;
  for( std::list<EvtRecTrack*>::iterator it=good_charged_tracks.begin();
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
  //std::cout << "Number of good charged tracks with emc = " << good_charged_tracks.size() << std::endl;

  emc_good_charged_tracks.sort(EmcEnergyOrder); //sort over deposited energy in EMC
  emc_good_charged_tracks.reverse(); //begin from hier energie

  /*  Select exactly 2 charged tracks
   *  And no good neutral tracks
   *  */
  if(    good_neutral_tracks.size() == 0 
      && cfg.MIN_CHARGED_TRACKS <= emc_good_charged_tracks.size()
      && emc_good_charged_tracks.size() <= cfg.MAX_CHARGED_TRACKS 
    ) //GOES to SIGNAL TauTau selection
  {
    //sort tracks on charge (lover charge goes first)
    emc_good_charged_tracks.sort(ChargeOrder);
    std::vector<EvtRecTrack*> Tracks
      (
        emc_good_charged_tracks.begin(), 
        emc_good_charged_tracks.end()
      );
    RecMdcTrack * mdc[2] =  { Tracks[0]->mdcTrack(), Tracks[1]->mdcTrack() };

    //opposite charge of this two tracks
    if( mdc[0]->charge()*mdc[1]->charge() >= 0) return StatusCode::SUCCESS;

    //std::cout << "Before fEvent.run " << std::endl;
    fEvent.run = eventHeader->runNumber();
    fEvent.event = eventHeader->eventNumber();
    fEvent.time = eventHeader->time();
    //std::cout << "eventTime = " << eventHeader->time() << std::endl;
    fEvent.channel = 0;
    fEvent.ntrack = 2;
    fEvent.Pid.init();
    for(int i=0;i<Tracks.size();++i)
    {
      fEvent.Pid.fill(i,Tracks[i]);
      fEvent.fill(i,Tracks[i]);
      //std::cout << i << "    " << GetCharge(Tracks[i]) << std::endl;
      //std::cout << i << "    " << GetMomentum(Tracks[i]) << std::endl;
    }

    fEvent.acop = M_PI - (fEvent.T.phi[1]  -  fEvent.T.phi[0] );

    //Hep3Vector p[2] = { GetHep3Vector(Tracks[0]), GetHep3Vector(Tracks[1])};
    Hep3Vector p[2] = { Tracks[0]->mdcKalTrack()->p3(), Tracks[1]->mdcKalTrack()->p3() };
    Hep3Vector psum  = p[0] +  p[1];
    Hep3Vector ptsum = p[0] +  p[1];
    ptsum.setZ(0);
    double Emis = cfg.CENTER_MASS_ENERGY - fEvent.T.E[0]-fEvent.T.E[1];

    fEvent.ptem = ptsum.mag() / Emis;
    fEvent.acol = (p[1].cross(p[0]).mag()/(p[1].mag()*p[0].mag()));
    fEvent.M2 = 0;
    fEvent.write();
    //double E[2] = { sqrt( 
  }
  else //gamma-gamma selection
  {
  }
  return StatusCode::SUCCESS;
}

StatusCode TauTau::finalize()
{
  std::cout << "Event proceed: " << nproceed_events << std::endl;
  std::cout << "Event selected: " << nwrited_events << std::endl;
  std::cout << "Tau candidates: " << ntautau_events << endl;
  std::cout << "Bhabha candidates: " << nbhabha_events << endl;
  std::cout << "Gamma-Gamma candidates: " << ngg_events << endl;
  std::cout << "Selection efficiency: " << nwrited_events/double(nproceed_events) << std::endl;
  return StatusCode::SUCCESS;
}
// for particle id look /ihepbatch/bes/alex/workarea/Analysis/Physics/PsiPrime/G2MuMuAlg-00-00-01/PipiJpsiAlg/src
