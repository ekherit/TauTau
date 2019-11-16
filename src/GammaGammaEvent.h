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

#include "Utils.h"
#include "SelectionConfig.h"

// =====================================================================================
//        Class:  GammaGammaEvent
//  Description:  Digamma events to measure luminosity
// =====================================================================================

class GammaGammaEvent : public RootTuple
{
  public:
    virtual ~GammaGammaEvent(void) {};
    GammaGammaEvent(void)
    {
      //default selections
      MAX_CHARGED_TRACKS_NUMBER = 0;
      MAX_NEUTRAL_TRACKS_NUMBER = 5;
      DELTA_THETA_CUT = 0.05;
      MIN_DELTA_PHI_CUT = -0.06;
      MAX_DELTA_PHI_CUT = 0.02;
      COS_THETA_CUT = 0.8; //barrel
      EEB_MIN_CUT = 0.8;
      EEB_MAX_CUT = 1.2;
    }
    long   MAX_CHARGED_TRACKS_NUMBER;
    long   MAX_NEUTRAL_TRACKS_NUMBER;
    double DELTA_THETA_CUT;
    double MIN_DELTA_PHI_CUT;
    double MAX_DELTA_PHI_CUT;
    double COS_THETA_CUT;
    double EEB_MIN_CUT;
    double EEB_MAX_CUT;


    NTuple::Item<long> run;           //run number
    NTuple::Item<long> event;         //event number 
    NTuple::Item<long> time;          //time of the event

    NTuple::Item<long> N0;            //number of neutral tracks
    NTuple::Item<long> Nq;            //number of charged tracks in event
    NTuple::Item<double> delta_phi;   //delta phi
    NTuple::Item<double> delta_theta; //delta theta
    NTuple::Array<double> E;          //energy deposited
    NTuple::Array<double> phi;        //phi
    NTuple::Array<double> theta;      //theta
    NTuple::Array<double> E_Eb;       //E/Ebeam


    //void init_tuple(Algorithm * algo, const char * dir, const char * title)
    //{
    //  RootTuple::init_tuple(algo,dir,title);
    //}
    virtual void bind_tuple(void)
    {
      tuple->addItem ("run", run);
      tuple->addItem ("event", event);
      tuple->addItem ("time", time);
      tuple->addItem ("N0", N0, 0,5); 
      tuple->addItem ("Nq", Nq); 
      tuple->addItem ("dphi", delta_phi); 
      tuple->addItem ("dtheta", delta_theta); 
      tuple->addIndexedItem ("E", N0, E);
      tuple->addIndexedItem ("E_Eb", N0, E_Eb);
      tuple->addIndexedItem ("phi", N0, phi);
      tuple->addIndexedItem ("theta", N0, theta);
    };

    virtual void init(void) {};

    virtual void fill(int i,  EvtRecTrack * track, double CENTER_MASS_ENERGY=3.6)
    {
      if(track->isEmcShowerValid())
      {
        E[i] = track->emcShower()->energy();
        E_Eb[i] = 2 * E[i] / CENTER_MASS_ENERGY;
      }
    }

    bool pass(const Event::EventHeader * eventHeader, double W, const  std::vector<EvtRecTrack*>  & Tc /* charged tracks */, std::vector<EvtRecTrack*>  Tn /* netural tracks */) 
    {
      bool result = true;
      N0 = Tn.size();
      Nq = Tc.size();
      result = result &&  2 <= N0;
      result = result && N0 <= MAX_NEUTRAL_TRACKS_NUMBER;
      result = result && Nq <= MAX_CHARGED_TRACKS_NUMBER;
      if(!result) return false;
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

      std::sort(Tn.rbegin(),Tn.rend(), EmcEnergyOrder);
      //std::reverse(Tn.begin(),Tn.end());
      bool keep=true;
      double E[2];
      double x[2]; //E/Ebeam
      for(int i=0;i<2;i++) {
        RecEmcShower * emc = Tn[i]->emcShower();
        E[i] = emc->energy();
        x[i] = 2.0*E[i] / W;
        phi[i] = emc->phi();
        theta[i] = emc->theta();
        keep = keep && ( EEB_MIN_CUT < x[i]  && x[i] < EEB_MAX_CUT );
        keep = keep && fabs(cos(theta[i])) < COS_THETA_CUT;
        E[i] = E[i];
        E_Eb[i] = x[i];
      }
      delta_theta =  theta[0] + theta[1] - M_PI;
      delta_phi = fabs(phi[1]  - phi[0]) - M_PI;
      keep = keep && fabs( delta_theta) < DELTA_THETA_CUT;
      keep = keep && delta_phi > MIN_DELTA_PHI_CUT;
      keep = keep && delta_phi < MAX_DELTA_PHI_CUT;
      result = keep;

      run   = eventHeader->runNumber();
      event = eventHeader->eventNumber();
      time  = eventHeader->time();

      return result;
    }
};

