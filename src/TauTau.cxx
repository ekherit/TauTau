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

#include "pair_comb.h"

inline double sq(double x) { return x*x; }

TauTau::TauTau(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator)
{
  declareProperty("CENTER_MASS_ENERGY"         , cfg.CENTER_MASS_ENERGY = 1.777*2); //GeV
  declareProperty("IP_MAX_Z"                   , cfg.IP_MAX_Z           = 10.0); //cm
  declareProperty("IP_MAX_RHO"                 , cfg.IP_MAX_RHO         = 1.0); //cm
  declareProperty("USE_VERTEX_DB"              , cfg.USE_VERTEX_DB      = 1);
  declareProperty("MAX_COS_THETA_FOR_CHARGED"  , cfg.MAX_COS_THETA_FOR_CHARGED = 0.93);
  declareProperty("MIN_EMC_ENERGY_FOR_CHARGED" , cfg.MIN_EMC_ENERGY_FOR_CHARGED=0.025); //GeV

  declareProperty("MIN_EMC_ENERGY_FOR_NEUTRAL" , cfg.MIN_EMC_ENERGY_FOR_NEUTRAL=0.025); //GeV

  declareProperty("MIN_MOMENTUM"               , cfg.MIN_MOMENTUM         = 0.1); //GeV
  declareProperty("MAX_MOMENTUM"               , cfg.MAX_MOMENTUM         = 1.5); //GeV

  declareProperty("MIN_TRANSVERSE_MOMENTUM"    , cfg.MIN_TRANSVERSE_MOMENTUM         = 0.1); //GeV
  declareProperty("MAX_TRANSVERSE_MOMENTUM"    , cfg.MAX_TRANSVERSE_MOMENTUM         = 1.5); //GeV


  declareProperty("MIN_EP_RATIO"               , cfg.MIN_EP_RATIO         = 0.05); 
  declareProperty("MAX_EP_RATIO"               , cfg.MAX_EP_RATIO         = 1.1);

  declareProperty("MIN_PTEM"               , cfg.MIN_PTEM         = 0.0); 
  declareProperty("MAX_PTEM"               , cfg.MAX_PTEM         = 1.5);
  //declareProperty("MIN_CHARGED_TRACKS", cfg.MIN_CHARGED_TRACKS=2); 
  //declareProperty("MAX_CHARGED_TRACKS", cfg.MAX_CHARGED_TRACKS=2); 


  //good netural tracks
  //declareProperty("EMC_ENDCUP_MIN_ENERGY",    cfg.EMC_ENDCUP_MIN_ENERGY = 0.05);
  //declareProperty("EMC_BARREL_MIN_ENERGY",    cfg.EMC_BARREL_MIN_ENERGY = 0.025);

  ////endcup calorimeter
  //declareProperty("EMC_ENDCUP_MIN_COS_THETA", cfg.EMC_ENDCUP_MIN_COS_THETA = 0.83); //was 0.86
  //declareProperty("EMC_ENDCUP_MAX_COS_THETA", cfg.EMC_ENDCUP_MAX_COS_THETA = 0.92);

  ////barrel calorimeter
  //declareProperty("EMC_BARREL_MAX_COS_THETA", cfg.EMC_BARREL_MAX_COS_THETA = 0.83); //was 0.8
  //declareProperty("NEUTRAL_CLOSE_CHARGED_ANGLE",    cfg.NEUTRAL_CLOSE_CHARGED_ANGLE = 10);

  //for tau+ -> pi+pi0
  declareProperty("GAMMA_GAMMA_MIN_INV_MASS", cfg.GAMMA_GAMMA_MIN_INV_MASS = 0.10);
  declareProperty("GAMMA_GAMMA_MAX_INV_MASS", cfg.GAMMA_GAMMA_MIN_INV_MASS = 0.20);
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
  cfg.print();
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
  if(nproceed_events==0) 
  {
    std::cout << "########################################################################################################\n";
    std::cout << "########################################################################################################\n";
    std::cout << "########################################################################################################\n";
    std::cout << "\n";
    cfg.print_relevant();
    std::cout << "########################################################################################################\n";
    std::cout << "########################################################################################################\n";
    std::cout << "########################################################################################################\n";
    //std::cout << "Wcm = " << cfg.CENTER_MASS_ENERGY << " GeV" << std::endl;
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
  std::sort(Tc.begin(),Tc.end(), ChargeOrder);

  fEvent.nciptrack = central_tracks.size(); //fill the total number of central tracks
  fEvent.nctrack   = tracker.GetNtrackCharged(); //save total number of all reconstructed charged tracks
  fEvent.nntrack   = tracker.GetNtrackNeutral(); //save total number of all reconstructed neutral tracks
  fEvent.Enmin     = tracker.MinNeutralTracksEnergy();
  fEvent.Enmax     = tracker.MaxNeutralTracksEnergy();

  //TAU TAU SELECTION
  if(
      Tc.size() == central_tracks.size()  //all central tracks has energy deposite in EMS
    &&
      GetTotalCharge(Tc) == 0  // opposite charge for tracks
    && 
      (    Tc.size() == 2  
        || Tc.size() == 4 
        || Tc.size() == 6 //till 3 charged particle decay for each tau
      ) 
    && 
      (    Tn.size() == 0 
        || Tn.size() == 2 
        || Tn.size() == 4 
        || Tn.size() == 6 
        || Tn.size() == 8 //till 2pi0 for each tau decay
      )
    )
  {
    bool select=true;
    fEvent.ngood_charged_track = Tc.size();
    fEvent.ngood_neutral_track = Tn.size();

    fEvent.Pid.init();

    std::vector<HepLorentzVector> Pc = GetMdcLorentzVector(Tc); //lorentz vector for charged tracks (electron hypoteza)
    std::vector<HepLorentzVector> Pn = GetEmcLorentzVector(Tn);

    HepLorentzVector Psum = GetTotalFourMomentum(Pc);
    Hep3Vector p3sum      = GetTotalMomentum(Pc);
    Hep3Vector ptsum      = GetTotalTransverseMomentum(Pc);
    //double psum = Psum.e(); //total energy of charged tracks ( electron hypoteza)

    for(int i=0;i<Tc.size();++i)
    {
      fEvent.Pid.fill(i, Tc[i]);
      fEvent.fill(i, Tc[i]);
      if(eventHeader->runNumber() < 0)
      {
        fEvent.McTruth.fill(i,Tc[i],mcParticleCol);
      }
    }
    fEvent.M2    = Psum.mag2();
    fEvent.Entot = tracker.GetTotalNeutralTracksEnergy(Tn);
    fEvent.Emis  = cfg.CENTER_MASS_ENERGY - Psum.e() - fEvent.Entot;
    fEvent.ptsum = ptsum.mag();
    fEvent.ptem  = fEvent.ptsum / fEvent.Emis;

    //acoplanarity and acolinearity for momentum with higher transverse momentum
    fEvent.acop = Acoplanarity(Tc[0], Tc[1]);
    fEvent.acol = Acolinearity(Tc[0], Tc[1]);

    std::vector<double> V = getSphericityEigenvalues(Tc);
    fEvent.S = Sphericity(V);
    fEvent.A = Aplanarity(V);
    fEvent.lambda1 = V[0];
    fEvent.lambda2 = V[1];
    fEvent.lambda3 = V[2];

    //find best pi0 combination
    //create combination list
    typedef std::list < std::pair<HepLorentzVector*, HepLorentzVector*> > comb_t;
    typedef std::vector< comb_t > comb_list_t;
    //std::vector<HepLorentzVector*> 
    //comb_list_t pi0_cmb_list;
    //make_unique_pairs(Pn.begin(),Pn.end(),pi0_cmb_list);
    comb_list_t pi0_cmb_list = make_combination_list(Pn); 
    //loop over all combinations
    comb_list_t::iterator best_comb=pi0_cmb_list.begin();
    double chi2_mass=1e100;
    for(comb_list_t::iterator it=pi0_cmb_list.begin(); it!=pi0_cmb_list.end(); ++it)
    {
      double chi2=0;
      for(comb_t::iterator it_pair = it->begin(); it_pair!=it->end(); ++it_pair)
      {
        //calculate invariant mass
        double m = (*(it_pair->first) +  *(it_pair->second)).mag();
        //add to chi square
        chi2+=pow(m-PI0_MASS,2.0);
      }
      if( chi2 < chi2_mass ) 
      {
        chi2_mass = chi2;
        best_comb = it;
      }
    }
    fEvent.npi0 = Pn.size()/2;
    //std::cout << "Before fEvent.Mpi0 filling " << std::endl;
    //std::cout << "neutral size = " << Pn.size() << std::endl;
    fEvent.Nrho = fEvent.npi0*Tc.size();
    if(pi0_cmb_list.size()!=0)
    {
      int idx=0;
      for(comb_t::iterator it_pair = best_comb->begin(); it_pair!=best_comb->end(); ++it_pair)
      {
        double m = (*(it_pair->first) +  *(it_pair->second)).mag(); //again calculate the pi0 mass
        fEvent.Mpi0[idx] = m;
        select &= fabs(m - PI0_MASS) <  0.03; //selection of the pi0
        //now create all combination to tie pi0 with charged tracks
        for(int i = 0; i<Tc.size(); ++i)
        {
          //RecMdcKalTrack * mdcTrk = track->mdcKalTrack();
          HepLorentzVector p = Tc[i]->mdcKalTrack()->p4(PION_MASS); //was error I should use PI+mass
          double Mrho = (p + *(it_pair->first) + *(it_pair->second)).mag();
          fEvent.Mrho[idx*Tc.size()+i] = Mrho;
        }
        idx++;
      }
    }
    select &=( MIN_PTEM < fEvent.ptem  && fEvent.ptem   < MAX_PTEM);
    for(int i=0;i<Tc.size();++i)
    {
      select &= ( cfg.MIN_MOMENTUM             < fEvent.T.p[i]      && fEvent.T.p[i]      < cfg.MAX_MOMENTUM);
      select &= ( cfg.MIN_TRANSVERSE_MOMENTUM  < fEvent.T.pt[i]     && fEvent.T.pt[i]     < cfg.MAX_TRANSVERSE_MOMENTUM);
      select &= ( cfg.MIN_EP_RATIO             < fEvent.T.Ep[i]     && fEvent.T.Ep[i]     < cfg.MAX_EP_RATIO);
      select &= ( cfg.MIN_TOF                  < fEvent.Pid.ftof[i] && fEvent.Pid.ftof[i] < cfg.MAX_TOF);
    }
    if(select)
    {
      fEvent.run   = eventHeader->runNumber();
      fEvent.event = eventHeader->eventNumber();
      fEvent.time  = eventHeader->time();
      fEvent.channel = 0;
      fEvent.write();
      ntautau_events++;
      //std::cout << " EMC_BAR_MIN = " << cfg.EMC_BARREL_MIN_ENERGY;
      //std::cout << " EMC_END_MIN = " << cfg.EMC_ENDCUP_MIN_ENERGY;
      //std::cout << " END_END_MIN_COS = " << cfg.EMC_ENDCUP_MIN_COS_THETA;
      //std::cout << " END_END_MAX_COS = " << cfg.EMC_ENDCUP_MAX_COS_THETA;
      //std::cout << " END_BARREL_MAX_COS = " << cfg.EMC_BARREL_MAX_COS_THETA;
    }
  }
SKIP_TAUTAU:
  //GAMMA GAMMA LUMINOCITY SELECTION
  fGG.N0 = Tn.size();
  fGG.Nq = Tc.size();
  if
  ( 
                2 <= fGG.N0 
        && fGG.N0 <= fGG.NEUTRAL_TRACKS_NUMBER
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


    std::sort(Tn.begin(),Tn.end(), EmcEnergyOrder);
    std::reverse(Tn.begin(),Tn.end());

    //std::vector<EvtRecTrack*> T
    //  (
    //    good_neutral_tracks.begin(), 
    //    good_neutral_tracks.end()
    //  );
    //fGG.Nq = good_charged_tracks.size();
    bool keep=true;
    double E[2];
    double x[2]; //E/Ebeam
    double theta[2];
    double phi[2];
    for(int i=0;i<2;i++)
    {
      RecEmcShower * emc = Tn[i]->emcShower();
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
