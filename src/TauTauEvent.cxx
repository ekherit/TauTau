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

#include "RootTauTauEvent.h"
#include "PhysConst.h"
#include "Utils.h"

RootTauTauEvent::~RootTauTauEvent(void)
{
}

const int UNSET_VALUE = -999;

void RootTauTauEvent::fill(int i,  EvtRecTrack * track)
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
    T.r[i] = 0; 

    Vertex_t vtx(track->mdcTrack());
    vtx.use_db(track->mdcTrack());
    //double rho,z,phi;
    //calculate_vertex(track->mdcTrack(),rho,z,phi);
    //T.vxy[i] = rho;
    //T.vz[i] = z; 
    //T.vphi[i] = phi; 
    T.vxy[i] = vtx.rho;
    T.z[i] = vtx.z;
    T.phi[i] = vtx.phi;
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

bool RootTauTauEvent::pass(const SelectionConfig & cfg, const SmartDataPtr<Event::EventHeader> & eventHeader const  std::vector<EvtRecTrack*>  & Tc, const  std::vector<EvtRecTrack*>  & Tn, const  std::vector<EvtRecTrack*>  & Tgn) 
{
  bool result = true;
  result = result && cfg.MIN_CHARGED_TRACKS <= Tc.size()  &&  Tc.size() <= cfg.MAX_CHARGED_TRACKS; //two charged tracks
  result = result && GetTotalCharge(Tc) == 0; //opposite charged tracks
  //result = result && Tn.size() == Tgn.size(); 
  result = result && 
    (
     Tc.size() == 2  
     || Tc.size() == 4 
     || Tc.size() == 6
    ); //till 3 charged particle decay for each tau
  result = result && 
    ( 
     Tn.size() == 0 
     //|| Tn.size() == 1 //this is for eta_c2  -> J/psi gamma
     || (Tgn.size() == 2  && Tn.size() == 2)
     || (Tgn.size() == 4  && Tn.size() == 2)
    );
  if(!result) resurn false;
  bool select=true;
  ngood_charged_track = Tc.size();
  ngood_neutral_track = Tn.size();

  Pid.init();

  std::vector<HepLorentzVector> Pc = GetMdcLorentzVector(Tc); //lorentz vector for charged tracks (electron hypoteza)
  std::vector<HepLorentzVector> Pn = GetEmcLorentzVector(Tgn);

  HepLorentzVector Psum = GetTotalFourMomentum(Pc);
  Hep3Vector P3sum      = GetTotalMomentum(Pc);
  Hep3Vector Ptsum      = GetTotalTransverseMomentum(Pc);

  for(int i=0;i<Tc.size();++i)
  {
    Pid.fill(i, Tc[i]);
    fill(i, Tc[i]);
    if(eventHeader->runNumber() < 0)
    {
      McTruth.fill(i,Tc[i],mcParticleCol);
    }
  }
  M2    = Psum.mag2();
  Entot = GetTotalNeutralTracksEnergy(Tn);
  Emis  = cfg.CENTER_MASS_ENERGY - Psum.e() - Entot;
  ptsum = Ptsum.mag();
  ptem  = ptsum / Emis;

  //acoplanarity and acolinearity for momentum with higher transverse momentum
  acop = Acoplanarity(Tc[0], Tc[1]);
  acol = Acolinearity(Tc[0], Tc[1]);

  //calculate sphericity 
  std::vector<double> V = getSphericityEigenvalues(Tc);
  S = Sphericity(V);
  A = Aplanarity(V);
  lambda1 = V[0];
  lambda2 = V[1];
  lambda3 = V[2];

  if (Tgn.size() % 2 == 0) {
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
    npi0 = Pn.size()/2;
    Nrho = npi0*Tc.size();
    if(pi0_cmb_list.size()!=0)
    {
      bool has_good_pi0 = false;
      int idx=0;
      for(comb_t::iterator it_pair = best_comb->begin(); it_pair!=best_comb->end(); ++it_pair)
      {
        double m = (*(it_pair->first) +  *(it_pair->second)).mag(); //again calculate the pi0 mass
        Mpi0[idx] = m;
        has_good_pi0 = has_good_pi0 || ( fabs(m - PI0_MASS) <  0.03); //selection of the pi0
        //now create all combination to tie pi0 with charged tracks
        for(int i = 0; i<Tc.size(); ++i)
        {
          //RecMdcKalTrack * mdcTrk = track->mdcKalTrack();
          HepLorentzVector p = Tc[i]->mdcKalTrack()->p4(PION_MASS);
          Mrho[idx*Tc.size()+i] = (p + *(it_pair->first) + *(it_pair->second)).mag();
        }
        idx++;
      }
      select &= has_good_pi0; //supress some events completely without pi0
    }
    select =  select && ( cfg.MIN_PTEM < ptem  && ptem   < cfg.MAX_PTEM);
    for(int i=0;i<Tc.size();++i)
    {
      select = select && ( cfg.MIN_MOMENTUM             < T.p[i]      && T.p[i]      < cfg.MAX_MOMENTUM);
      select = select && ( cfg.MIN_TRANSVERSE_MOMENTUM  < T.pt[i]     && T.pt[i]     < cfg.MAX_TRANSVERSE_MOMENTUM);
      select = select && ( cfg.MIN_EP_RATIO             < T.Ep[i]     && T.Ep[i]     < cfg.MAX_EP_RATIO);
      select = select && ( cfg.MIN_TOF                  < Pid.ftof[i] && Pid.ftof[i] < cfg.MAX_TOF);
    }
  } else if( Tgn.size() == 1) //chi_c2 -> Jpsi gamma
  {
    double Mjpsi = M2 > 0 ? sqrt(M2) : 0;
    select = select && ( fabs(Mjpsi - JPSI_MASS) < 0.02 );
  } else  {
    select = false;
  }
  result = select;
  fEvent.run   = eventHeader->runNumber();
  fEvent.event = eventHeader->eventNumber();
  fEvent.time  = eventHeader->time();
  fEvent.channel = 0;
  return result;
}
