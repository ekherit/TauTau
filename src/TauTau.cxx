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
  declareProperty("MAX_COS_THETA", cfg.MAX_COS_THETA = 0.93);

  declareProperty("EMC_ENDCUP_MIN_ENERGY",    cfg.EMC_ENDCUP_MIN_ENERGY = 0.05);
  declareProperty("EMC_BARREL_MIN_ENERGY",    cfg.EMC_BARREL_MIN_ENERGY = 0.025);

  //endcup calorimeter
  declareProperty("EMC_ENDCUP_MIN_COS_THETA", cfg.EMC_ENDCUP_MIN_COS_THETA = 0.86);
  declareProperty("EMC_ENDCUP_MAX_COS_THETA", cfg.EMC_ENDCUP_MAX_COS_THETA = 0.92);
  //barrel calorimeter
  declareProperty("EMC_BARREL_MAX_COS_THETA", cfg.EMC_BARREL_MAX_COS_THETA = 0.8);
  declareProperty("NEUTRAL_CLOSE_CHARGED_ANGLE",    cfg.NEUTRAL_CLOSE_CHARGED_ANGLE = 10);

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
    fEvent.make_tuple(this, "FILE1/tt","Signal tau tau events");
    fGG.make_tuple(this, "FILE1/gg","Two gamma (luminocity) events");
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
  if(nproceed_events==0) std::cout << "Wcm = " << cfg.CENTER_MASS_ENERGY << " GeV" << std::endl;
  if(10 <= nproceed_events && nproceed_events < 100 && nproceed_events % 10 ==0) isprint=true;
  if(100 <= nproceed_events && nproceed_events < 1000 && nproceed_events % 100 ==0) isprint = true;
  if(1000 <= nproceed_events && nproceed_events < 10000 && nproceed_events % 1000 ==0) isprint = true;
  if(10000 <= nproceed_events && nproceed_events % 10000 ==0) isprint = true;
  if(isprint)
  {
    std::cout << "proceed event: " << setw(15) << nproceed_events;
    std::cout << "  ττ:" << setw(15) << ntautau_events;
    std::cout << "  ee:" << setw(15) << nbhabha_events;
    std::cout << "  γγ:" << setw(15) << ngg_events;
    std::cout << std::endl;
  }
  nproceed_events++;

  //  Get information about reconstructed events
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),  EventModel::MC::McParticleCol);

	std::list<EvtRecTrack*> good_neutral_tracks;
	std::list<EvtRecTrack*> good_charged_tracks;

  //fill initial value of the selected event
  //come from IP and cos(theta)<max_cos_theta //define in utils
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

  //SELECTION
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

    //SELECTION charge
    //opposite charge of this two tracks
    if( mdc[0]->charge()*mdc[1]->charge() >= 0) goto SKIP_TAUTAU;

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
      if(eventHeader->runNumber() < 0)
      {
        fEvent.McTruth.fill(i,Tracks[i],mcParticleCol);
      }
      //std::cout << i << "    " << GetCharge(Tracks[i]) << std::endl;
      //std::cout << i << "    " << GetMomentum(Tracks[i]) << std::endl;
    }

    // Calculate acoplanarity
    double phi[2] = {fEvent.T.phi[0],  fEvent.T.phi[1]};
    double dphi = (phi[1]-phi[0]);
    if(dphi>M_PI) dphi = dphi-2*M_PI;
    if(dphi<-M_PI) dphi = dphi+2*M_PI;
    fEvent.acop = M_PI - fabs(dphi);

    //calculae missing energy and ptem
    //Hep3Vector p[2] = { GetHep3Vector(Tracks[0]), GetHep3Vector(Tracks[1])};
    Hep3Vector p[2] = { Tracks[0]->mdcKalTrack()->p3(), Tracks[1]->mdcKalTrack()->p3() };
    Hep3Vector psum  = p[0] +  p[1];
    Hep3Vector ptsum = p[0] +  p[1];
    ptsum.setZ(0);
    double Emis = cfg.CENTER_MASS_ENERGY - fEvent.T.p[0]-fEvent.T.p[1];

    fEvent.ptsum =  ptsum.mag();
    fEvent.ptem = ptsum.mag() / Emis;
    fEvent.acol = (p[1].cross(p[0]).mag()/(p[1].mag()*p[0].mag()));
    fEvent.M2 = 0;
    bool select=true;
    //SELECTION
    for(int i=0;i<2;++i)
    {
      select &=                     fabs(cos(fEvent.T.theta[i])) < 0.8; //goes to barrel
      select &=    0  < fEvent.T.p[i]      && fEvent.T.p[i]      < 2.0;
      select &=  0.05 < fEvent.T.Ep[i]     && fEvent.T.Ep[i]     < 1.1;
      select &=   2.5 < fEvent.Pid.ftof[i] && fEvent.Pid.ftof[i] < 5.5;
      select &=  0.00 < fEvent.ptem        && fEvent.ptem        < 1.1;
    }
    if(select)
    {
      fEvent.write();
      ntautau_events++;
    }
  }
SKIP_TAUTAU:
  //GAMMA GAMMA LUMINOCITY SELECTION
  fGG.N0 = good_neutral_tracks.size();
  fGG.Nq = good_charged_tracks.size();
  if
  ( 
        2 <= fGG.N0 && fGG.N0 <= fGG.NEUTRAL_TRACKS_NUMBER
      && fGG.Nq <= fGG.CHARGED_TRACKS_NUMBER
  )
  {

    /*
    std::cout << fGG.N0 << " " << fGG.Nq << std::endl;
    std::cout << " N0 cut = " << fGG.NEUTRAL_TRACKS_NUMBER;
    std::cout << " Nq cut = " << fGG.CHARGED_TRACKS_NUMBER;
    std::cout << " E/B min = " << fGG.EEB_MIN_CUT;
    std::cout << " E/B max = " << fGG.EEB_MAX_CUT;
    std::cout << " cos = " << fGG.COS_THETA_CUT;
    std::cout << " delta theta = " << fGG.DELTA_THETA_CUT;
    std::cout << " delta phi min = " << fGG.MIN_DELTA_PHI_CUT;
    std::cout << " delta phi max = " << fGG.MAX_DELTA_PHI_CUT;
    std::cout << " EMC_BAR_MIN = " << cfg.EMC_BARREL_MIN_ENERGY;
    std::cout << " EMC_END_MIN = " << cfg.EMC_ENDCUP_MIN_ENERGY;
    std::cout << " END_END_MIN_COS = " << cfg.EMC_ENDCUP_MIN_COS_THETA;
    std::cout << " END_END_MAX_COS = " << cfg.EMC_ENDCUP_MAX_COS_THETA;
    std::cout << " END_BARREL_MAX_COS = " << cfg.EMC_BARREL_MAX_COS_THETA;
    std::cout << std::endl;
    */


    good_neutral_tracks.sort(EmcEnergyOrder);
    good_neutral_tracks.reverse();
    std::vector<EvtRecTrack*> T
      (
        good_neutral_tracks.begin(), 
        good_neutral_tracks.end()
      );
    fGG.Nq = good_charged_tracks.size();
    bool keep=true;
    double E[2];
    double x[2]; //E/Ebeam
    double theta[2];
    double phi[2];
    for(int i=0;i<2;i++)
    {
      RecEmcShower * emc = T[i]->emcShower();
      E[i] = emc->energy();
      x[i] = 2.0*E[i] / cfg.CENTER_MASS_ENERGY;
      phi[i] = emc->phi();
      theta[i] = emc->theta();
      keep = keep && ( fGG.EEB_MIN_CUT < x[i]  && x[i] < fGG.EEB_MAX_CUT );
      keep = keep && fabs(cos(theta[i])) < fGG.COS_THETA_CUT;
      fGG.E[i] = E[i];
      fGG.E_Eb[i] = x[i];
      fGG.phi[i] = phi[i];
      fGG.theta[i] = theta[i];
    }
    fGG.delta_theta =  theta[0] + theta[1] - M_PI;
    fGG.delta_phi = fabs(phi[1]  - phi[0]) - M_PI;
    keep = keep && fabs( fGG.delta_theta) < fGG.DELTA_THETA_CUT;
    keep = keep && fGG.delta_phi > fGG.MIN_DELTA_PHI_CUT;
    keep = keep && fGG.delta_phi < fGG.MAX_DELTA_PHI_CUT;
    //std::cout << fGG.delta_theta << " " << fGG.delta_phi << std::endl;
    if(keep) 
    {
      fGG.write();
      ngg_events++;
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode TauTau::finalize()
{
  std::cout << "Event proceed: " << nproceed_events << std::endl;
  std::cout << "Event selected: " << nwrited_events << std::endl;
  std::cout << "ττ candidates: " << ntautau_events << endl;
  std::cout << "Bhabha candidates: " << nbhabha_events << endl;
  std::cout << "γγ candidates: " << ngg_events << endl;
  std::cout << "Selection efficiency: " << nwrited_events/double(nproceed_events) << std::endl;
  return StatusCode::SUCCESS;
}
// for particle id look /ihepbatch/bes/alex/workarea/Analysis/Physics/PsiPrime/G2MuMuAlg-00-00-01/PipiJpsiAlg/src
