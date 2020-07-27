// =====================================================================================
//
//       Filename:  RootEvent.cxx
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 16:59:57
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#include "CLHEP/Vector/LorentzVector.h"

//#include "McTruth/McParticle.h"
//#include "EventModel/EventHeader.h"

#include "TauTauEvent.h"
#include "PhysConst.h"
#include "Utils.h"
#include "Tracker.h"
#include "combinator.h"
#include "pair_comb.h"

TauTauEvent::~TauTauEvent(void)
{
}

const int UNSET_VALUE = -999;

void TauTauEvent::fill(int i,  EvtRecTrack * track)
{
  if(track->isMdcTrackValid())
  {
    //RecMdcKalTrack * mdc = track->mdcKalTrack();
    RecMdcTrack * mdc = track->mdcTrack();
    T.id[i] = mdc->trackId(); //id of the track
    T.q[i] =  mdc->charge(); //charge of the track
    T.p[i] = mdc->p();
    if(track->isEmcShowerValid())
    {
      T.E[i] = track->emcShower()->energy();
      T.Ep[i] = T.E[i]/T.p[i];
      T.temc[i] = track->emcShower()->time();
    }
    else
    {
      T.E[i] = UNSET_VALUE;
      T.Ep[i] = UNSET_VALUE;
      T.temc[i] = UNSET_VALUE;
    }
    T.px[i] = mdc->px();
    T.py[i] = mdc->py();
    T.pz[i] = mdc->pz();
    T.pt[i] = mdc->pxy();
    T.theta[i]= mdc->theta();
    T.phi[i] = mdc->phi();
    T.x[i] = mdc->x();
    T.y[i] = mdc->y();
    T.z[i] = mdc->z(); 
    T.r[i] = sqrt(mdc->x()*mdc->x() + mdc->y()*mdc->y());

    Vertex_t vtx(track->mdcTrack());
    vtx.use_db(track->mdcTrack());
    T.vxy[i] = vtx.rho;
    T.vz[i] = vtx.z;
    T.vphi[i] = vtx.phi;
  }
  else
  {
    T.id[i] = UNSET_VALUE;
    T.q[i]  = UNSET_VALUE;
    T.E[i]  = UNSET_VALUE; 
    T.p[i]  = UNSET_VALUE; 
    T.px[i] = UNSET_VALUE; 
    T.py[i] = UNSET_VALUE; 
    T.pz[i] = UNSET_VALUE; 
    T.pt[i] = UNSET_VALUE; 
    T.theta[i]= UNSET_VALUE; 
    T.phi[i] =  UNSET_VALUE;
    T.x[i] = UNSET_VALUE;
    T.y[i] = UNSET_VALUE;
    T.z[i] = UNSET_VALUE; 
    T.r[i] = UNSET_VALUE; 
    T.vxy[i] = UNSET_VALUE;
    T.vz[i] = UNSET_VALUE; 
    T.vphi[i] = UNSET_VALUE; 
  }
  if(track->isMucTrackValid())
  {
    RecMucTrack *mucTrk = track->mucTrack();
    T.depth[i]= mucTrk->depth();
    T.Nmuhit[i] = mucTrk->numHits();
  }
  else 
  {
    T.depth[i] =  UNSET_VALUE;
    T.Nmuhit[i] = UNSET_VALUE;
  }
}

//fill the neutral tracks
void TauTauEvent::nfill(int i,  EvtRecTrack * track)
{
  RecEmcShower * emc = track->emcShower();
  double E = emc->energy();
  double phi = emc->phi();
  double theta = emc->theta();
  Tn.E[i] = E;
  Tn.id[i] = emc->trackId(); //id of the track
  Tn.q[i] = 0; //charge of the track
  Tn.theta[i]= theta;
  Tn.phi[i] = phi;
  Tn.p[i] = E;
  Tn.Ep[i] = 1;
  double st = sin(theta);
  Tn.px[i] = E*st*cos(phi);
  Tn.py[i] = E*st*sin(phi);
  Tn.pz[i] = E*cos(theta);
  Tn.pt[i] = E*st;
  Tn.x[i] = 0;
  Tn.y[i] = 0;
  Tn.z[i] = 0;
  Tn.r[i] = 0; 
  Tn.vxy[i] = 0;
  Tn.vz[i] = 0;
  Tn.vphi[i] = 0;
  Tn.temc[i] = emc->time();
  if(track->isMucTrackValid())
  {
    RecMucTrack *mucTrk = track->mucTrack();
    Tn.depth[i]= mucTrk->depth();
    Tn.Nmuhit[i] = mucTrk->numHits();
  }
  else 
  {
    Tn.depth[i] =  UNSET_VALUE;
    Tn.Nmuhit[i] = UNSET_VALUE;
  }
}

void FillIndexedArrayWithFourMomentum(NTuple::Array<double>  & d4p, const HepLorentzVector  & dP) {
    d4p[0] = dP.e();
    d4p[1] = dP.vect().mag();
    d4p[2] = dP.perp().mag();
    d4p[3] = dP.vect().theta();
    d4p[4] = dP.vect().phi();
    d4p[5] = dP.m2();
};
typedef std::pair<const HepLorentzVector*, const HepLorentzVector*> LorentzPairPtr;
typedef std::list < LorentzPairPtr > comb_t;
typedef std::vector< comb_t > comb_list_t;

inline bool CloseToMpi0Order(const LorentzPairPtr & p1, const LorentzPairPtr & p2) {
  double m1 = (*(p1.first) +  *(p1.second)).mag();
  double m2 = (*(p2.first) +  *(p2.second)).mag();
  return fabs(m1-PI0_MASS) < fabs(m2-PI0_MASS);
}

inline bool CloseToMrhoOrder(double m1,double m2) {
  return fabs(m1-RHO_MASS) < fabs(m2-RHO_MASS);
}

inline bool CloseToMrhoOrder(const std::pair<double, int> & p1, const std::pair<double, int> & p2) {
  return fabs(p1.first-RHO_MASS) < fabs(p2.first-RHO_MASS);
}

comb_t MakePi0List(const std::vector<HepLorentzVector> & Png) {
    assert( Png.size() >=2 );
    //create combination list
    comb_list_t pi0_cmb_list = make_combination_list(Png); //all pi0 combination 
    //loop over all combinations
    comb_list_t::iterator best_comb=pi0_cmb_list.begin();
    double chi2_mass=std::numeric_limits<double>::max();
    //Find best combinations"
    for(comb_list_t::iterator it=pi0_cmb_list.begin(); it!=pi0_cmb_list.end(); ++it) {
      double chi2=0;
      for(comb_t::iterator it_pair = it->begin(); it_pair!=it->end(); ++it_pair) {
        //calculate invariant mass
        double m = (*(it_pair->first) +  *(it_pair->second)).mag();
        //add to chi square
        chi2+=pow(m-PI0_MASS,2.0);
      }
      if( chi2 < chi2_mass ) {
        chi2_mass = chi2;
        best_comb = it;
      }
    }
    std::sort(best_comb->begin(), best_comb->end(), CloseToMpi0Order);
    return *best_comb;
};

void TauTauEvent::fillPi0Rho(const comb_t & pi0s, std::vector<EvtRecTrack*> & Tq) {
    size_t npi0 = pi0s.size() < MAX_PI0_NUMBER ? pi0s.size() : MAX_PI0_NUMBER;
    Npi0 = npi0;

    std::vector<int> NotEidxs;  //with idxes in the Tq which are not electrons
    for(size_t i=0;i < Tq.size(); ++i) {
      RecEmcShower * emc = Tq[i]->emcShower();
      RecMdcKalTrack * mdc = Tq[i]->mdcKalTrack();
      double E = emc->energy();
      double p = mdc->p();
      if( E/p < 0.8 ) {
    //  T.push_back(Tq[i]);
        NotEidxs.push_back(i);
      }
    };

    std::vector< std::pair<double, int> > Rhos;
    Rhos.reserve(npi0*NotEidxs.size());
    size_t idx=0;
    for(comb_t::iterator it_pair = pi0s.begin(); idx < npi0; ++it_pair, ++idx) {
      double m = (*(it_pair->first) +  *(it_pair->second)).mag(); //again calculate the pi0 mass
      //fill  Mpi0 mass
      Mpi0[idx] = m;

      //Now loop over Not electron and calculate Rho mass
      for(int i = 0; i<NotEidxs.size(); ++i) {
        HepLorentzVector p = Tq[RhoIdxs[i]]->mdcKalTrack()->p4(PION_MASS);
        double mrho = (p + *(it_pair->first) + *(it_pair->second)).mag();
        Rhos.push_back(  std::pair<double, int>(mrho, NotEidxs[i]) );
      }
    }
    std::sort(Rhos.begin(), Rhos.end(), CloseToMrhoOrder);
    size_t nrho = Rhos.size() < MAX_RHO_NUMBER ? Rhos.size() : MAX_RHO_NUMBER;
    Nrho = nrho;
    for(size_t i=0;i <nrho; ++i) {
      Mrho[i] = Rhos[i].first;
      IdxRho[i] = Rhos[i].second;
    }
}
//bool TauTauEvent::pass(const SelectionConfig & cfg, const Event::EventHeader *  eventHeader, const Event::McParticleCol * mcParticleCol,  const  std::vector<EvtRecTrack*>  & Tc, const  std::vector<EvtRecTrack*>  & Tn, const  std::vector<EvtRecTrack*>  & Tgn) 
bool TauTauEvent::pass(const SelectionConfig & cfg, const Event::EventHeader *  eventHeader, const Event::McParticleCol * mcParticleCol,  const  Tracker & tracker) 
{
  //central tracks T_{ineracion point}, T_ip
  Tracker::Vector  Tcc     = tracker.GetCentralTracks<Tracker::Vector>(cfg.IP_MAX_Z, cfg.IP_MAX_RHO, cfg.USE_VERTEX_DB);
  Tracker::Vector  Tce     = FilterMdcTracksByEmcEnergy(Tip, 0);
  Tracker::Vector  Tcg     = FilterMdcTracksByCosTheta(Tce, 0.93); //good charged tracks with emc and from IP
  Tracker::Vector  Tn      = tracker.GetNeutralTracks<Tracker::Vector>(0); //all neutral tracks in the event
  Tracker::Vector  Tng     = tracker.GetGoodNeutralTracks<Tracker::Vector>();
  std::sort(Tn.rbegin(),Tn.rend(), EmcEnergyOrder); //sort by deccending energy
  std::sort(Tng.rbegin(),Tng.rend(), EmcEnergyOrder);

  enmin     = tracker.MinNeutralTracksEnergy();
  enmax     = tracker.MaxNeutralTracksEnergy();
  entot     = tracker.GetTotalNeutralTracksEnergy();


  Nc  = tracker.GetNtrackCharged(); //total number of charged tracks
  Ncc = Tcc.size(); //number of charged tracks from IP
  Nce = Tce.size(); //number of charged tracks with EMC from IP
  Ncg = Tcg.size(); //number of good charged tracks

  Nn     = Tn.size(); //total number of neutral tracks.
  NnE10  = tracker.CountNeutralTracks(0.010); //number of neutrals more then 10 MeV
  NnE25  = tracker.CountNeutralTracks(0.025); //more 25 MeV
  NnE50  = tracker.CountNeutralTracks(0.050); //more 50 MeV
  NnE100 = tracker.CountNeutralTracks(0.100); // more 100 MeV
  Nng    = Tng.size(); //number of good neutrals

  bool select = true;
  select = select && 2 <= Ncg  &&  Ncg <= MAX_CHARGED_TRACKS;
  select = select && 0 <= Nn   &&   Nn <= MAX_NEUTRAL_TRACKS;
  if(!select) return false; 

  //Now devide track into sign of charge Tcs - Tracks charged splited (by charge)
  std::vector< std::vector<EvtRecTrack*>  > Tcs = SplitByCharge(Tcg);
  //Tcs[0] -- negative charged
  //Tcs[1] -- positive charged
  if(Tcs[0].empty() || Tcs[1].empty()) return false; //Must be one opposite charged pair
  //sort it by transverse momentum in dessceding order
  std::sort(Tcs[0].rbegin(), Tcs[0].rend(),MomentumOrder);
  std::sort(Tcs[1].rbegin(), Tcs[1].rend(),MomentumOrder);

  std::vector<EvtRecTrack*> Tq = Zip(Tcs[0],Tcs[1], true);//JoinTracks(T);

  Pid.init();

  std::vector<HepLorentzVector> Pc  = GetMdcLorentzVector(Tq); //lorentz vector for charged tracks (electron hypoteza)
  std::vector<HepLorentzVector> Png = GetEmcLorentzVector(Tng); //loretz vectos of good neutrals
  std::vector<HepLorentzVector> Pn  = GetEmcLorentzVector(Tn); //lorentz vectors of all neutrals

  HepLorentzVector Psum = GetTotalFourMomentum(Pc);   //Total 4 momentum of good emc charged
  HepLorentzVector Pnsum = GetTotalFourMomentum(Pn);  //Total momentum of all neutrals

  Enmin = MinNeutralTracksEnergy(Tn);
  Enmax = MaxNeutralTracksEnergy(Tn);
  Entot = GetTotalNeutralTracksEnergy(Tn,0.0);
  Entot10 = GetTotalNeutralTracksEnergy(Tn,0.010);
  Entot25 = GetTotalNeutralTracksEnergy(Tn,0.025);
  Entot50 = GetTotalNeutralTracksEnergy(Tn,0.050);
  Entot100 = GetTotalNeutralTracksEnergy(Tn,0.100);


  Emis  = cfg.CENTER_MASS_ENERGY - Psum.e() - Entot; 
  Emis0 = cfg.CENTER_MASS_ENERGY - Psum.e();
  Emis1 = cfg.CENTER_MASS_ENERGY - (Pc[0]+Pc[1]).e() - Entot;
  Emis2 = cfg.CENTER_MASS_ENERGY - (Pc[0]+Pc[1]).e();

  ptsum = Psum.perp();
  ptsum2 = (Pc[0]+Pc[1]).perp();

  psum = Psum.vect();
  psum2 = (Pc[0]+Pc[1]).vect();

  ptem  = ptsum / Emis;
  ptem10  = ptsum / (cfg.CENTER_MASS_ENERGY - Psum.e() - Entot10);
  ptem25  = ptsum / (cfg.CENTER_MASS_ENERGY - Psum.e() - Entot25);
  ptem50  = ptsum / (cfg.CENTER_MASS_ENERGY - Psum.e() - Entot50);

  FillIndexedArrayWithFourMomentum(P4, Psum);
  HepLorentzVector P0 = getInitialFourMomentum(cfg.CENTER_MASS_ENERGY); //initial Momentum

  ndp = 6;
  HepLorentzVector dP  = Psum - P0;
  for(size_t i=0; i<d4p.size(); ++i) {
    FillIndexedArrayWithFourMomentum(d4p[i], dP);
    if( i >= Pn.size()) dP += Pn[i];
  }

  HepLorentzVector dPg  = Psum - P0;
  for(size_t i=0; i<d4pg.size(); ++i) {
    FillIndexedArrayWithFourMomentum(d4pg[i], dPg);
    if( i >= Png.size()) dP += Png[i];
  }

  //Save good charged tracks from IP point with EMC signal
  for(size_t i=0;i<Tq.size();++i) {
    Pid.fill(i, Tq[i]);
    fill(i, Tq[i]);
    if(eventHeader->runNumber() < 0) {
      McTruth.fill(i,Tq[i],mcParticleCol);
    }
  }
  for(size_t i=0;i<Tn.size();++i) {
    nfill(i, Tn[i]);
  }

  M2    = Psum.mag2();

  //std::cout << " Before acop calculation" << std::endl;
  //acoplanarity and acolinearity for momentum with higher transverse momentum
  acop = Acoplanarity(Tq[0], Tq[1]);
  acol = Acolinearity(Tq[0], Tq[1]);
  
  cos_theta_mis  = Psum.vect().z()/Psum.vect.mag();
  //std::cout << "cos_theta_mis = " << cos_theta_mis << std::endl;

  //calculate sphericity 
  std::vector<double> V = getSphericityEigenvalues(Tq);
  S = Sphericity(V);
  A = Aplanarity(V);
  lambda1 = V[0];
  lambda2 = V[1];
  lambda3 = V[2];

  select =  select && ( 0 < ptem); 
  //for(int i=0;i<Tq.size();++i) {
  for(int i=0;i<2;++i) { //look only for 2 tracks
    select = select && (T.p[i] < 1.2);
    select = select && (0.05  < T.pt[i] );
    select = select && (1 < Pid.ftof[i] && Pid.ftof[i] < 6);
  }

  //Making pi0
  if (Png.size() >= 2) {
    comb_t pi0s = MakePi0List(Png); //make all gamma gamma combinations orderer by close to pi0
    fillPi0Rho(pi0s); //fill tuple by resuls. In this function we keep only most close to pi0 and rho
  }

  run   = eventHeader->runNumber();
  event = eventHeader->eventNumber();
  time  = eventHeader->time();
  channel = 0;
  return select; 
}
