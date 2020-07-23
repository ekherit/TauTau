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
    }
    else
    {
      T.E[i] = UNSET_VALUE;
      T.Ep[i] = UNSET_VALUE;
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

bool TauTauEvent::pass(const SelectionConfig & cfg, const Event::EventHeader *  eventHeader, const Event::McParticleCol * mcParticleCol,  const  std::vector<EvtRecTrack*>  & Tc, const  std::vector<EvtRecTrack*>  & Tn, const  std::vector<EvtRecTrack*>  & Tgn) 
{
  //Tc - хорошие заряженные центральные треки, имеющией энерговыделение в калориметре
  //Tn - кластеры с энерговыделением больше порогового в любом месте калориметра
  //Tgn - хорошие нейтральные кластеры
  bool select = true;
  select = select && cfg.MIN_CHARGED_TRACKS <= Tc.size()  &&  Tc.size() <= cfg.MAX_CHARGED_TRACKS;
  select = select && cfg.MIN_NEUTRAL_TRACKS <= Tn.size()  &&  Tn.size() <= cfg.MAX_NEUTRAL_TRACKS;
  if(!select) return false;  //number of tracks cut

  //Now devide track into sign of charge
  std::vector< std::vector<EvtRecTrack*>  > Ts = SplitByCharge(Tc);
  //T[0] -- negative charged
  //T[1] -- positive charged
  if(Ts[0].empty() || Ts[1].empty()) return false; //Must be one opposite charged pair
  //sort it by transverse momentum in dessceding order
  std::sort(Ts[0].rbegin(), Ts[0].rend(),PtOrder);
  std::sort(Ts[1].rbegin(), Ts[1].rend(),PtOrder);

  std::vector<EvtRecTrack*> Tq = Zip(Ts[0],Ts[1], true);//JoinTracks(T);

  Nc = Tq.size();
  Nn = Tn.size(); 

  Pid.init();

  std::vector<HepLorentzVector> Pc = GetMdcLorentzVector(Tq); //lorentz vector for charged tracks (electron hypoteza)

  HepLorentzVector Psum = GetTotalFourMomentum(Pc);
  Hep3Vector P3sum      = GetTotalMomentum(Pc);
  Hep3Vector Ptsum      = GetTotalTransverseMomentum(Pc);

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
  Enmin = MinNeutralTracksEnergy(Tn);
  Enmax = MaxNeutralTracksEnergy(Tn);
  Entot = GetTotalNeutralTracksEnergy(Tn);
  Emis  = cfg.CENTER_MASS_ENERGY - Psum.e() - Entot;
  ptsum = Ptsum.mag();
  ptem  = ptsum / Emis;

  //std::cout << " Before acop calculation" << std::endl;
  //acoplanarity and acolinearity for momentum with higher transverse momentum
  acop = Acoplanarity(Tq[0], Tq[1]);
  acol = Acolinearity(Tq[0], Tq[1]);
  
  cos_theta_mis  = Ptsum.z()/P3sum.mag();

  //calculate sphericity 
  std::vector<double> V = getSphericityEigenvalues(Tq);
  S = Sphericity(V);
  A = Aplanarity(V);
  lambda1 = V[0];
  lambda2 = V[1];
  lambda3 = V[2];

  select =  select && (cfg.MIN_PTEM < ptem); 
  for(int i=0;i<2;++i) {
    select = select && (T.p[i] < cfg.MAX_MOMENTUM);
    select = select && (cfg.MIN_TRANSVERSE_MOMENTUM  < T.pt[i] );
    select = select && (cfg.MIN_TOF                  < Pid.ftof[i] && Pid.ftof[i] < cfg.MAX_TOF);
  }

  std::vector<HepLorentzVector> Pn = GetEmcLorentzVector(Tgn);
  if (Pn.size() >= 2) {
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
    //std::cout << "Before making all combination" << std::endl;
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
    npi0 = Pn.size()/2;
    Nrho = npi0*Tq.size();
    //std::cout << "Find good pi0" << std::endl;
    if(pi0_cmb_list.size()!=0) {
      bool has_good_pi0 = false;
      int idx=0;
      for(comb_t::iterator it_pair = best_comb->begin(); it_pair!=best_comb->end(); ++it_pair) {
        double m = (*(it_pair->first) +  *(it_pair->second)).mag(); //again calculate the pi0 mass
        Mpi0[idx] = m;
        has_good_pi0 = has_good_pi0 || ( fabs(m - PI0_MASS) <  0.03); //selection of the pi0
        //now create all combination to tie pi0 with charged tracks
        for(int i = 0; i<Tq.size(); ++i) {
          //RecMdcKalTrack * mdcTrk = track->mdcKalTrack();
          HepLorentzVector p = Tq[i]->mdcKalTrack()->p4(PION_MASS);
          Mrho[idx*Tq.size()+i] = (p + *(it_pair->first) + *(it_pair->second)).mag();
        }
        idx++;
      }
      select &= has_good_pi0; //supress some events completely without pi0
    }
  } else if( Tgn.size() == 1) //chi_c2 -> Jpsi gamma
  {
    //TAU TAU SELECTION and chi_c1 -> Jpsi gamma selection
    double Mjpsi = M2 > 0 ? sqrt(M2) : 0;
    select = select && ( fabs(Mjpsi - JPSI_MASS) < 0.02 );
  } 
  run   = eventHeader->runNumber();
  event = eventHeader->eventNumber();
  time  = eventHeader->time();
  channel = 0;
  return select; 
}
