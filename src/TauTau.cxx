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
#include "combinator.h"

#include "Utils.h"

#include "Tracker.h"

#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"

#include "ParticleID/ParticleID.h"


inline double sq(double x) { return x*x; }

TauTau::TauTau(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator)
{
  declareProperty("CENTER_MASS_ENERGY"         , cfg.CENTER_MASS_ENERGY         = 1.777*2); //GeV

  declareProperty("MIN_CHARGED_TRACKS"         , cfg.MIN_CHARGED_TRACKS         = 2);
  declareProperty("MAX_CHARGED_TRACKS"         , cfg.MAX_CHARGED_TRACKS         = 3);

  declareProperty("MIN_NEUTRAL_TRACKS"         , cfg.MIN_NEUTRAL_TRACKS         = 0);
  declareProperty("MAX_NEUTRAL_TRACKS"         , cfg.MAX_NEUTRAL_TRACKS         = 4);

  declareProperty("IP_MAX_Z"                   , cfg.IP_MAX_Z                   = 15.0); //cm
  declareProperty("IP_MAX_RHO"                 , cfg.IP_MAX_RHO                 = 1.5); //cm

  declareProperty("USE_VERTEX_DB"              , cfg.USE_VERTEX_DB              = 1);

  declareProperty("MAX_COS_THETA_FOR_CHARGED"  , cfg.MAX_COS_THETA_FOR_CHARGED  = 0.93);
  declareProperty("MIN_EMC_ENERGY_FOR_CHARGED" , cfg.MIN_EMC_ENERGY_FOR_CHARGED = 0.0); //GeV

  declareProperty("MIN_EMC_ENERGY_FOR_NEUTRAL" , cfg.MIN_EMC_ENERGY_FOR_NEUTRAL = 0.01); //GeV

  declareProperty("MAX_MOMENTUM"               , cfg.MAX_MOMENTUM               = 1.2); //GeV

  declareProperty("MIN_TRANSVERSE_MOMENTUM"    , cfg.MIN_TRANSVERSE_MOMENTUM    = 0.05); //GeV
  //declareProperty("MAX_TRANSVERSE_MOMENTUM"    , cfg.MAX_TRANSVERSE_MOMENTUM    = 1.5); //GeV


  //Will not use
  //declareProperty("MIN_EP_RATIO"               , cfg.MIN_EP_RATIO               = 0.0);
  //declareProperty("MAX_EP_RATIO"               , cfg.MAX_EP_RATIO               = 10);

  declareProperty("MIN_PTEM"                     , cfg.MIN_PTEM                   = 0);
  //declareProperty("MAX_PTEM"                   , cfg.MAX_PTEM                   = 2);

  declareProperty("MIN_TOF"                    , cfg.MIN_TOF                    = 1);
  declareProperty("MAX_TOF"                    , cfg.MAX_TOF                    = 6);

  declareProperty("DELTA_MJPSI"                , cfg.DELTA_MJPSI                   = 0.2);
  declareProperty("TEST_COMBINATIONS", cfg.TEST_COMBINATIONS=0);

}



StatusCode TauTau::initialize(void)
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in initialize()" << endmsg;

  if(cfg.TEST_COMBINATIONS==1)
  {
    // Test for  pair combination
    std::vector<int> A;
    for(int i=1;i<6;++i) A.push_back(i);
    std::cout << "Test for homogenios array: \n";
    print_array(" A = ", A);
    std::cout << "\n";
    std::vector<std::vector<std::pair<int,int> > > R;
    make_unique_pairs(A.begin(),A.end(),R);
    print(R);
    std::cout << std::flush;

    std::vector<char> B;
    for(int i=0;i<3;i++) B.push_back('a'+i);
    std::vector<std::vector<std::pair<int,char> > > Rab;
    std::cout << "Test for homogenios array: \n";
    print_array(" A = ", A);
    std::cout << "\n";
    print_array(" B = ", B);
    std::cout << "\n";
    make_unique_pairs(A.begin(),A.end(),B.begin(),B.end(),Rab);
    print(Rab);
    std::cout << std::flush;
    exit(0);
  }

  //init counters
  nproceed_events=0;
  nwritten_events=0;
  ntautau_events=0;
  nbhabha_events=0;
  ngg_events=0;
  StatusCode status;
  try
  {
    fTT.make_tuple(this,    "FILE1/tt","Signal tau tau events");
    fGG.make_tuple(this,    "FILE1/gg","Two gamma (luminosity) events");
    fBB.make_tuple(this,    "FILE1/bb","Bhabha (luminosity) events");
    fJobInfo.make_tuple(this,  "FILE1/info","Job information");
    fJobInfo.begin_time = time(nullptr);
  }
  catch(std::runtime_error & error)
  {
    log << MSG::ERROR << error.what() << endmsg;
    return StatusCode::FAILURE;
  }
  cfg.print();
  return StatusCode::SUCCESS;
}

StatusCode TauTau::execute()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  int runNo=eventHeader->runNumber();
  if(runNo<0) fJobInfo.type = 1; //monte carlo
  int event=eventHeader->eventNumber();
  time_t t=eventHeader->time();
  bool isprint=false;
  if(nproceed_events<10) isprint=true;
  if(nproceed_events==0) 
  {
    cfg.print_relevant();
  }
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
  SmartDataPtr<EvtRecEvent>            evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol>        evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),  EventModel::MC::McParticleCol);


  Tracker tracker(evtRecEvent, evtRecTrkCol); //helper class for excracting information about tracks
  Tracker::Vector  central_tracks = tracker.GetCentralTracks<Tracker::Vector>(cfg.IP_MAX_Z, cfg.IP_MAX_RHO, cfg.USE_VERTEX_DB);
  Tracker::Vector Tc              = FilterMdcTracksByEmcEnergy(central_tracks, cfg.MIN_EMC_ENERGY_FOR_CHARGED);
  Tracker::Vector Tn              = tracker.GetNeutralTracks<Tracker::Vector>(cfg.MIN_EMC_ENERGY_FOR_NEUTRAL);
  Tracker::Vector Tgn             = tracker.GetGoodNeutralTracks<Tracker::Vector>();

  //std::sort(Tc.begin(),Tc.end(), ChargeOrder);

  /* ****************** TAU PAIR SELECTION **********************************/
  //if ( Tc.size() == central_tracks.size()  && fTT.pass(cfg, eventHeader.ptr(), mcParticleCol.ptr(), Tc,Tn,Tgn))  //all central tracks has energy deposite in EMS
  //if ( fTT.pass(cfg, eventHeader.ptr(), mcParticleCol.ptr(), Tc,Tn,Tgn))
  if ( fTT.pass(cfg, eventHeader.ptr(), mcParticleCol.ptr(), Tc,Tn,Tgn))
  {
    fTT.nctrack   = tracker.GetNtrackCharged(); //save total number of all reconstructed charged tracks
    fTT.nciptrack = central_tracks.size(); //fill the total number of central tracks

    fTT.nntrack   = tracker.GetNtrackNeutral(); //save total number of all reconstructed neutral tracks
    fTT.enmin     = tracker.MinNeutralTracksEnergy();
    fTT.enmax     = tracker.MaxNeutralTracksEnergy();
    fTT.entot     = tracker.GetTotalNeutralTracksEnergy();

    fTT.McTruth.flag1=eventHeader->flag1();
    fTT.McTruth.flag2=eventHeader->flag2();
    fTT.write();
    ntautau_events++;
    nwritten_events++;
  }


  /* *******************  SELECT  DIGAMMA  EVENTS ********************************** */
  //see selection detail in GammaGammaEvent.h
  if(fGG.pass(eventHeader.ptr(), cfg.CENTER_MASS_ENERGY, Tc,Tgn) )
  {
    fGG.write();
    ngg_events++;
    nwritten_events++;
  }

  /* *******************  SELECT  BHABHA EVENTS ********************************** */
  //see selection detail in BhabhaEvent.h
  if(fBB.pass(eventHeader.ptr(), Tc,Tn))
  {
    //std::cout << "Before BB write" << std::endl;
    fBB.nctrack   = tracker.GetNtrackCharged(); //save total number of all reconstructed charged tracks
    fBB.nntrack   = tracker.GetNtrackNeutral(); //save total number of all reconstructed neutral tracks
    fBB.enmin     = tracker.MinNeutralTracksEnergy();
    fBB.enmax     = tracker.MaxNeutralTracksEnergy();
    fBB.entot     = tracker.GetTotalNeutralTracksEnergy();
    fBB.write();
    //std::cout << "After BB write" << std::endl;
    nbhabha_events++;
    nwritten_events++;
  }

  return StatusCode::SUCCESS;
}

StatusCode TauTau::finalize()
{
  cfg.print_relevant();
  std::cout << "Event proceed: "        << nproceed_events                        << std::endl;
  std::cout << "ττ candidates: "        << ntautau_events                         << std::endl;
  std::cout << "Bhabha candidates: "    << nbhabha_events                         << std::endl;
  std::cout << "γγ candidates: "        << ngg_events                             << std::endl;
  std::cout << "Total written: "        << nwritten_events                        << std::endl;
  std::cout << "Selection efficiency: " << ntautau_events/double(nproceed_events) << std::endl;
  fJobInfo.end_time = time(nullptr);
  fJobInfo.N = nproceed_events;
  fJobInfo.n = nwritten_events;
  fJobInfo.write();
  return StatusCode::SUCCESS;
}
// for particle id look /ihepbatch/bes/alex/workarea/Analysis/Physics/PsiPrime/G2MuMuAlg-00-00-01/PipiJpsiAlg/src
