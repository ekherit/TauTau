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

#include "RootEvent/RootTuple.h" #include "CLHEP/Vector/LorentzVector.h"

#include "SelectionConfig.h"

// =====================================================================================
//        Class:  GammaGammaEvent
//  Description:   Two gamma   
// =====================================================================================

class GammaGammaEvent : public RootTuple
{
  public:
    virtual ~GammaGammaEvent(void) {};
    GammaGammaEvent(void)
    {
      //default selections
      CHARGED_TRACKS_NUMBER = 0;
      NEUTRAL_TRACKS_NUMBER = 2;
      DELTA_THETA_CUT = 0.05;
      MIN_DELTA_PHI_CUT = -0.06;
      MAX_DELTA_PHI_CUT = 0.02;
      COS_THETA_CUT = 0.8; //barrel
      EEB_MIN_CUT = 0.8;
      EEB_MAX_CUT = 1.2;
    }
    long   CHARGED_TRACKS_NUMBER;
    long   NEUTRAL_TRACKS_NUMBER;
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
      tuple->addItem ("N0", N0, 0,2); 
      tuple->addItem ("Nq", Nq, 0,2); 
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
};
