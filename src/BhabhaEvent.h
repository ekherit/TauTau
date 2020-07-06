// =====================================================================================
//
//       Filename:  GammaGammaEvent.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 16:52:08
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once

#include <vector>

#include "EventModel/EventHeader.h"

#include "RootEvent/RootTuple.h" 
#include "CLHEP/Vector/LorentzVector.h"

//#include "SelectionConfig.h"

#include "Utils.h"
#include "Tracker.h"
// =====================================================================================
//        Class:  BhabhaEvent
//  Description:  Bhabha event to measure luminosity  
// =====================================================================================

class BhabhaEvent : public RootTuple
{
  public:
    virtual ~BhabhaEvent(void) {};
    BhabhaEvent(void)
    {
      //default preselection selections
      MAX_CHARGED_TRACKS_NUMBER = 3;
      MAX_NEUTRAL_TRACKS_NUMBER = 15;
      COS_THETA_CUT = 0.93; 
      DELTA_THETA_CUT = 0.04;
      MIN_DELTA_PHI_CUT = -0.06;
      MAX_DELTA_PHI_CUT = 0.06;
      MIN_PEB_CUT=0.8;
      MAX_PEB_CUT=1.2;
      MIN_EEB_CUT=0.8;
      MAX_EEB_CUT=1.2;
    }
    long   MAX_CHARGED_TRACKS_NUMBER;
    long   MAX_NEUTRAL_TRACKS_NUMBER;
    double DELTA_THETA_CUT;
    double MIN_DELTA_PHI_CUT;
    double MAX_DELTA_PHI_CUT;
    double COS_THETA_CUT;
    double MIN_PEB_CUT;
    double MAX_PEB_CUT;
    double MIN_EEB_CUT;
    double MAX_EEB_CUT;


    NTuple::Item<long> run;           //run number
    NTuple::Item<long> event;         //event number 
    NTuple::Item<long> time;          //time of the event

    RootTracks T;  //track information (momentum, vertex, muon depth...)

    NTuple::Item<long> nntrack;       //raw number of neutral tracks
    NTuple::Item<long> nctrack;       //raw number of charged tracks

    NTuple::Item<double>  Enmin; //minimum energy of all neutral tracks (not only selected)
    NTuple::Item<double>  Enmax; //maximum energy of all neutral tracks (not only selected)
    NTuple::Item<double>  Entot; //total energy of neutral tracks

    NTuple::Item<long> N0;            //number of neutral tracks with energy more then threshold
    NTuple::Item<long> Nq;            //number of charged tracks in event
    NTuple::Item<double> delta_phi;   //delta phi
    NTuple::Item<double> delta_theta; //delta theta
    NTuple::Item<double> acol;        //Acolinearity
    NTuple::Array<double> E_Eb;       //E/Ebeam
    NTuple::Array<double> p_Eb;       //p/Ebeam

    virtual void bind_tuple(void)
    {
      tuple->addItem ("run", run);
      tuple->addItem ("event", event);
      tuple->addItem ("time", time);
      tuple->addItem ("nntrack", nntrack); //raw number of neutral tracks
      tuple->addItem ("nctrack", nctrack); //raw number of charged tracks
      tuple->addItem ("Enmin",Enmin);
      tuple->addItem ("Enmax",Enmax);
      tuple->addItem ("Entot",Entot);

      tuple->addItem ("Nn", N0, 0,5);  //number of neutral tracks with energy more theshold
      tuple->addItem ("Nc", Nq, 0,3); 
      tuple->addItem ("dphi", delta_phi); 
      tuple->addItem ("dtheta", delta_theta); 
      tuple->addItem("acol",acol);
      tuple->addIndexedItem ("E_Eb", Nq, E_Eb);
      tuple->addIndexedItem ("p_Eb", Nq, p_Eb);
      T.add_to_tuple (tuple,Nq); 
    };

    virtual void init(void) {};

    virtual void fill(int i,  EvtRecTrack * track, double CENTER_MASS_ENERGY=3.6) 
    {
      const int UNSET_VALUE = -999;
      if(track->isMdcTrackValid())
      {
        //RecMdcKalTrack * mdc = track->mdcKalTrack();
        RecMdcTrack * mdc = track->mdcTrack();
        T.id[i] = mdc->trackId(); //id of the track
        T.q[i] =  mdc->charge(); //charge of the track
        T.p[i] = mdc->p();
        p_Eb[i] = 2 * T.p[i] / CENTER_MASS_ENERGY;
        if(track->isEmcShowerValid())
        {
          T.E[i] = track->emcShower()->energy();
          T.Ep[i] = T.E[i]/T.p[i];
          E_Eb[i] = 2 * T.E[i] / CENTER_MASS_ENERGY;
        }
        else
        {
          T.E[i] = UNSET_VALUE;
          T.Ep[i] = UNSET_VALUE;
          E_Eb[i] = UNSET_VALUE;
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


    bool pass(const Event::EventHeader * eventHeader, const  std::vector<EvtRecTrack*>  & Tc, const  std::vector<EvtRecTrack*>  & Tn) 
    {
      bool result = true;
      //fill event number for printing
      run   = eventHeader->runNumber();
      event = eventHeader->eventNumber();
      time  = eventHeader->time();

      result = result && Tc.size() >= 2;
      result = result && Tc.size() <= MAX_CHARGED_TRACKS_NUMBER;
      result = result && Tn.size() <= MAX_NEUTRAL_TRACKS_NUMBER;
      if(!result) return false; //earlier rejection 
      N0 = Tn.size();
      Nq = Tc.size(); 

      std::vector< std::vector<EvtRecTrack*> >  Ts = SplitByCharge(Tc);

      //Alwais there is event one pair
      for(int i=0;i<2;++i) result = result && !(Ts[i].empty()); 
      if(!result) return false;

      //sort by EMC energy order
      for(int i=0;i<3;++i) std::sort(Ts[i].rbegin(), Ts[i].rend(),EmcEnergyOrder);

      std::vector<EvtRecTrack*> Tr = Zip(Ts[0],Ts[1]);

      for(int i=0;i<Tr.size(); ++i) fill(i, Tr[i]); 


      acol = Acolinearity(Tr[0], Tr[1]);
      delta_theta =  T.theta[0] + T.theta[1] - M_PI;
      delta_phi =    fabs(T.phi[1]  - T.phi[0]) - M_PI;
      //print_event(result);

      //Calculate pass result
      result = result && fabs( delta_theta) < DELTA_THETA_CUT;
      //print_event(result);
      result = result && delta_phi > MIN_DELTA_PHI_CUT;
      //print_event(result);
      result = result && delta_phi < MAX_DELTA_PHI_CUT;
      //print_event(result);
      for(int i=0;i<2;++i) {
        result = result && cos(T.theta[i]) < COS_THETA_CUT;
        result = result && (MIN_EEB_CUT < E_Eb[i])  && (E_Eb[i]  < MAX_EEB_CUT);
      }
      //print_event(result);
      return result;
    }

    void print_event(bool pass) {
      std::cout << "event = " << event << " N0 = " << N0 << " Nq = " << Nq << " Δθ = " << delta_theta << "  Δφ=" << delta_phi << " θ[0] = " << T.theta[0] << " θ[1] = " << T.theta[1];
      std::cout << "  E[0]/Eb = " << E_Eb[0] << "   E[1]/Eb = " << E_Eb[1];
      std::cout << " pass = " << pass;
      std::cout << endl;
    }
};
