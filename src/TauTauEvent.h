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
#include "RootEvent/RootMcTruth.h"
#include "RootEvent/RootPid.h"


#include "SelectionConfig.h"

//class EventModel::MC::McParticleCol;
//class Event::EventHeader;

#include "McTruth/McParticle.h"
#include "EventModel/EventHeader.h"

class EvtRecTrack;

// =====================================================================================
//        Class:  RootEvent
//  Description:  Main event information supposed to used for selection
// =====================================================================================
class TauTauEvent : public RootTuple
{
  public:
    virtual ~TauTauEvent(void);
    NTuple::Item<long>    run; //run number
    NTuple::Item<long>    event; //event number 
    NTuple::Item<long>    time; //time of the event
    NTuple::Item<long>    nctrack; //total number of charged tracks;
    NTuple::Item<long>    nntrack; //total number of neutral tracks;

    NTuple::Item<double>  enmin; //min energy of all neutral tracks
    NTuple::Item<double>  enmax; //max energy of all neutral tracks


    NTuple::Item<long>    nciptrack; //total number of charged tracks came from interaction points

    NTuple::Item<long>    Nc;     //number of good charged tracks in event
    NTuple::Item<long>    Nn;     //number of good neutral tracks in event
    NTuple::Item<long>    npi0; //number of pi0

    NTuple::Item<long> channel; //used channel
    //NTuple::Item<long> ntrack;  //number of tracks (must be 2)
    RootTracks T;  //track information (momentum, vertex, muon depth...)
    RootTracks Tn; //neutral track
    RootPid Pid;   //particle id for charged track
    //NTuple::Array<double> pid; //particle id
    //NTuple::Array<double> mother_pid;
    RootMcTruth  McTruth; // mc truth
    //RootMcTruth nMcTruth; // mc truth for nuetral tracks
    NTuple::Item<double>  ptsum;
    NTuple::Item<double>  ptem;
    NTuple::Item<double>  acop;
    NTuple::Item<double>  acol;
    NTuple::Item<double>  cos_theta_mis;
    NTuple::Item<double>  M2;
    NTuple::Item<double>  S; //sphericity (lambda2+lambda3)*1.5
    NTuple::Item<double>  A; //aplanarity 1.5*lambda3
    NTuple::Item<double>  lambda1; //eigen values of sphericity tenzor
    NTuple::Item<double>  lambda2; //
    NTuple::Item<double>  lambda3; //
    NTuple::Item<double>  Emis; //

    NTuple::Item<double>  Enmin; //minimum energy of accepted neutral tracks
    NTuple::Item<double>  Enmax; //maximum energy of accepted neutral tracks
    NTuple::Item<double>  Entot; //total energy of the acceptd neutral tracks

    NTuple::Array<double> Mpi0;

    NTuple::Item<long>    Nrho; //number of rho candidates
    NTuple::Array<double> Mrho; //this mass of rho candidates


    virtual void bind_tuple(void)
    {
      tuple->addItem ("run", run);
      tuple->addItem ("event", event);
      tuple->addItem ("time", time);
      tuple->addItem ("channel", channel);
      tuple->addItem ("nctrack",   nctrack); 
      tuple->addItem ("nntrack",   nntrack); 

      tuple->addItem ("enmin",enmin);
      tuple->addItem ("enmax",enmax);

      tuple->addItem ("nciptrack", nciptrack, 0,8); //total number of charged tracks come from interaction points
      tuple->addItem ("Nc", Nc, 6);
      tuple->addItem ("Nn", Nn, 0, 8);
      tuple->addItem ("Npi0", npi0, 0, 4);
      tuple->addItem ("Nrho", Nrho, 0, 24);


      tuple->addItem ("Enmin",Enmin);
      tuple->addItem ("Enmax",Enmax);
      tuple->addItem ("Entot",Entot);

      T.add_to_tuple (tuple,Nc); 
      Tn.add_to_tuple(tuple,Nn,"n");
      Pid.add_to_tuple(tuple,Nc); 
      McTruth.add_to_tuple(tuple,Nc);

      //nMcTruth.add_to_tuple(tuple,ntrack,"n");
      tuple->addItem("acop",acop);
      tuple->addItem("acol",acol);
      tuple->addItem("ptem",ptem);
      tuple->addItem("ptsum",ptsum);
      tuple->addItem("Emis",Emis);
      tuple->addItem("cos_theta_mis", cos_theta_mis);
      tuple->addItem("M2",M2);
      tuple->addItem("S",S);
      tuple->addItem("A",A);
      tuple->addItem("l1",lambda1);
      tuple->addItem("l2",lambda2);
      tuple->addItem("l3",lambda3);
      tuple->addIndexedItem("Mpi0", npi0, Mpi0);
      tuple->addIndexedItem("Mrho", Nrho, Mrho);
    };
    virtual void init(void) {};
    virtual void fill(int i,  EvtRecTrack * track);

    bool pass(const  std::vector<EvtRecTrack*>  & Tc, const  std::vector<EvtRecTrack*>  & Tn, const SelectionConfig & cfg);
    bool pass(const SelectionConfig & cfg, const Event::EventHeader * eventHeader, const Event::McParticleCol *,  const  std::vector<EvtRecTrack*>  & Tc, const  std::vector<EvtRecTrack*>  & Tn, const  std::vector<EvtRecTrack*>  & Tgn);
};
