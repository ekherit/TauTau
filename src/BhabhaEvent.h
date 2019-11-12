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
      //default selections
      MAX_CHARGED_TRACKS_NUMBER = 3;
      MAX_NEUTRAL_TRACKS_NUMBER = 5;
      COS_THETA_CUT = 0.8; //barrel
      DELTA_THETA_CUT = 0.05;
      MIN_DELTA_PHI_CUT = -0.06;
      MAX_DELTA_PHI_CUT = 0.02;
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

    NTuple::Item<long> N0;            //number of neutral tracks
    NTuple::Item<long> Nq;            //number of charged tracks in event
    NTuple::Item<double> delta_phi;   //delta phi
    NTuple::Item<double> delta_theta; //delta theta
    NTuple::Array<double> E_Eb;       //E/Ebeam
    NTuple::Array<double> p_Eb;      // p/Ebeam
    tuple->addItem("acol",acol);

    virtual void bind_tuple(void)
    {
      tuple->addItem ("run", run);
      tuple->addItem ("event", event);
      tuple->addItem ("time", time);
      tuple->addItem ("N0", N0, 0,MAX_NEUTRAL_TRACKS_NUMBER); 
      tuple->addItem ("Nq", Nq, 0,MAX_CHARGED_TRACKS_NUMBER); 
      tuple->addItem ("dphi", delta_phi); 
      tuple->addItem ("dtheta", delta_theta); 
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


    bool pass(const SmartDataPtr<Event::EventHeader> & eventHeader, const  std::vector<EvtRecTrack*>  & Tc, const  std::vector<EvtRecTrack*>  & Tn) 
    {
      bool result = true;
      result = result && Tc.size() < MAX_CHARGED_TRACKS_NUMBER;
      result = result && Tn.size() < MAX_NEUTRAL_TRACKS_NUMBER;
      if(!result) return false; //earlier rejection 
      N0 = Tn.size();
      Nq = Tc.size(); 
      //split track list into positive and negative
      std::vector<EvtRecTrack*> negative_tracks;
      std::vector<EvtRecTrack*> positive_tracks;
      for(int i=0;i<Tc.size();++i) {
        RecMdcTrack * mdc = Tc[i]->mdcTrack();  
        if(mdc->charge() < 0) negative_tracks.push_back(Tc[i]);
        if(mdc->charge() > 0) positive_tracks.push_back(Tc[i]);
      }
      //sort them on Momentum order
      std::sort(negative_tracks.begin(),negative_tracks.end(), MomentumOrder);
      std::sort(positive_tracks.begin(),positive_tracks.end(), MomentumOrder);

      std::vector<EvtRecTrack*> Tr; //result tracks for filling the tuple
      size_t npairs = std::min( negative_tracks.size(), positive_tracks.size()); //number of pairs
      for(int i=0; i < npairs; ++i) {
        Tr.push_back(negative_tracks[i]); //negative charged track goes first
        Tr.push_back(positive_tracks[i]);
      }
      std::vector<EvtRecTrack*> & tmp_tracks = negative_tracks.size() > positive_tracks.size() ?  negative_tracks :  positive_tracks;
      for(int i = npairs; i < tmp_trakcs.size() : ++i)  Tr.push_back(tmp_tracks[i]);


      //now fill the tuple
      for(int i=0;i<Tr.size(); ++i) fill(i, Tr[i]); 

      acol = Acolinearity(Tr[0], Tr[1]);
      delta_theta =  T.theta[0] + T.theta[1] - M_PI;
      delta_phi =    fabs(T.phi[1]  - T.phi[0]) - M_PI;

      //Calculate pass result
      result = result && fabs( delta_theta) < DELTA_THETA_CUT;
      result = result && delta_phi > MIN_DELTA_PHI_CUT;
      result = result && delta_phi < MAX_DELTA_PHI_CUT;
      for(int i=0;i<2;++i) {
        result = result && T[i].theta[i] < COS_THETA_CUT;
        result = result && MIN_EEB_CUT < T[i].E_Eb  && T[i].E_EB  < MAX_EEB_CUT;
      }
      run   = eventHeader->runNumber();
      event = eventHeader->eventNumber();
      time  = eventHeader->time();
      return result;
    }
};
