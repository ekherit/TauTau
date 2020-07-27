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
      MAX_CHARGED_TRACKS_NUMBER = 0; //number of maximum good charged tracks
      MAX_NEUTRAL_TRACKS_NUMBER = 11;
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
    NTuple::Array<double> t;       //time


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
      tuple->addIndexedItem ("t", N0, t);
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

      std::sort(Tn.rbegin(),Tn.rend(), EmcEnergyOrder);
      bool keep=true;
      for(int i=0;i<2;i++) {
        RecEmcShower * emc = Tn[i]->emcShower();
        E[i] = emc->energy();
        E_Eb[i] = 2.0*E[i] / W;
        phi[i] = emc->phi();
        theta[i] = emc->theta();
        t[i] = emc->time();
        keep = keep && ( (EEB_MIN_CUT < E_Eb[i])  && (E_Eb[i] < EEB_MAX_CUT ) );
        keep = keep && fabs(cos(theta[i])) < COS_THETA_CUT;
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

