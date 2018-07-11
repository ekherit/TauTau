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

#include "RootEvent/RootTuple.h"
#include "RootEvent/RootTrack.h"

#include "RootEvent/RootPid.h"
#include "RootEvent/RootMcTruth.h"

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

class RootTauTauEvent : public RootTuple
{
  public:
  virtual ~RootTauTauEvent(void);
	NTuple::Item<long>    run; //run number
	NTuple::Item<long>    event; //event number 
	NTuple::Item<long>    time; //time of the event
	NTuple::Item<long>    ngood_charged_track;     //number of good charged tracks in event
	NTuple::Item<long>    ngood_neutral_track;     //number of good neutral tracks in event
	NTuple::Item<long>    nctrack;     //number of total charged tracks
  NTuple::Item<long>    nntrack;     //number of total neutral tracks

  /*
   *  00 - ee  |  10 - μe  |  20 - Ke
   *  01 - eμ  |  11 - μμ  |  21 - Kμ
   *  02 - eπ  |  12 - μπ  |  22 - Kπ
   *  03 - eK  |  13 - μK  |  23 - KK
   */ 
	NTuple::Item<long> channel; //used channel
	NTuple::Item<long> ntrack;  //number of tracks (must be 2)
	RootTracks T;  //track information (momentum, vertex, muon depth...)
	RootTracks Tn; //neutral track
  RootPid Pid;   //particle id for charged track
  RootMcTruth  McTruth;
  RootMcTruth nMcTruth; // mc truth for nuetral tracks
  NTuple::Item<double>  ptsum;
  NTuple::Item<double>  ptem;
  NTuple::Item<double>  acop;
  NTuple::Item<double>  acol;
  NTuple::Item<double>  M2;

  //void init_tuple(Algorithm * algo, const char * dir, const char * title)
  //{
  //  RootTuple::init_tuple(algo,dir,title);
  //}
	virtual void bind_tuple(void)
  {
    tuple->addItem ("run", run);
    tuple->addItem ("event", event);
    tuple->addItem ("time", time);
    tuple->addItem ("channel", channel);
    tuple->addItem ("ntrack", ntrack, 0,3);
    tuple->addItem ("Nc", ngood_charged_track, 0, 3);
    tuple->addItem ("Nn", ngood_neutral_track, 0, 6);
    T.add_to_tuple (tuple,ngood_charged_track); 
    Tn.add_to_tuple(tuple,ngood_neutral_track,"n");
    Pid.add_to_tuple(tuple,ntrack); 
    McTruth.add_to_tuple(tuple,ntrack);
    nMcTruth.add_to_tuple(tuple,ntrack,"n");
    tuple->addItem("acop",acop);
    tuple->addItem("acol",acol);
    tuple->addItem("ptem",ptem);
    tuple->addItem("ptsum",ptsum);
    tuple->addItem("M2",M2);
  };
	virtual void init(void) {};
	virtual void fill(int i,  EvtRecTrack * track);
};
