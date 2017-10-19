// =====================================================================================
//
//       Filename:  RootTauTauEvent.h
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

#include "RootTuple.h"
#include "RootMass.h"
#include "RootTrack.h"

#include "CLHEP/Vector/LorentzVector.h"

// =====================================================================================
//        Class:  RootEvent
//  Description:  Main event information supposed to used for selection
// =====================================================================================
struct RootTauTauEvent : public RootTuple
{
	NTuple::Item<long>    run; //run number
	NTuple::Item<long>    event; //event number 
	NTuple::Item<long>    time; //time of the event
	NTuple::Item<long>    ngood_charged_track;     //number of good charged tracks in event
	NTuple::Item<long>    ngood_neutral_track;     //number of good neutral tracks in event

	NTuple::Item<long>    channel;     //J/psi decay channel 0 -- kaons, 1 -- muons,  10 -- muK,  11 - Kmu
	Track_t T;

	NTuple::Array<double> pid_chi2[5];      //probability from ParticleID
	NTuple::Array<double> pid_prob[5];      //particle id probability
	NTuple::Array<double> prob;      //probability from ParticleID

	virtual void init(void);
	virtual void init_tuple(void);
	virtual void fill(int i,  EvtRecTrack * track);
};
