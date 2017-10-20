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

enum CHANNEL
{
   EE  =  0,  MUE  =  10,  KE  =  20,
   EMU =  1,  MUMU =  11,  KMU =  21,
   EPI =  2,  MUPI =  12,  KPI =  22,
   EK  =  3,  MUK  =  13,  KK  =  23 
};

struct RootTauTauEvent : public RootTuple
{
	NTuple::Item<long>    run; //run number
	NTuple::Item<long>    event; //event number 
	NTuple::Item<long>    time; //time of the event
	NTuple::Item<long>    ngood_charged_track;     //number of good charged tracks in event
	NTuple::Item<long>    ngood_neutral_track;     //number of good neutral tracks in event

  /*
   *  00 - ee  |  10 - μe  |  20 - Ke
   *  01 - eμ  |  11 - μμ  |  21 - Kμ
   *  02 - eπ  |  12 - μπ  |  22 - Kπ
   *  03 - eK  |  13 - μK  |  23 - KK
   */ 
	NTuple::Item<long> channel; //used channel
	NTuple::Item<long> ntrack;  //number of tracks (must be 2)
	RootTracks T; //track information (momentum, vertex, muon depth...)
  RootPid Pid; //particle id for track

	virtual void init_tuple(void)
  {
    tuple->addItem ("channel", channel);
    tuple->addItem ("ntrack", ntrack, 2,2);
    T.add_to_tuple(tuple,ntrack); 
    Pid.add_to_tuple(tuple,ntrack); 
  };
	virtual void init(void)
	virtual void fill(int i,  EvtRecTrack * track);
};
