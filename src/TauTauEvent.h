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

#include "Tracker.h"

#include "SelectionConfig.h"

//class EventModel::MC::McParticleCol;
//class Event::EventHeader;

#include "McTruth/McParticle.h"
#include "EventModel/EventHeader.h"

class EvtRecTrack;

typedef std::pair<HepLorentzVector*, HepLorentzVector*> LorentzPairPtr;
typedef std::list < LorentzPairPtr > comb_t;

// =====================================================================================
//        Class:  RootEvent
//  Description:  Main event information supposed to used for selection
// =====================================================================================

class TauTauEvent : public RootTuple
{
  const int MAX_NEUTRAL_TRACKS; 
  const int MAX_CHARGED_TRACKS; 
  const int MAX_PI0_NUMBER;   
  const int MAX_RHO_NUMBER;  
  public:
    TauTauEvent(void);
    virtual ~TauTauEvent(void);
    NTuple::Item<long>    run; //run number
    NTuple::Item<long>    event; //event number 
    NTuple::Item<long>    time; //time of the event
    //NTuple::Item<long>    nctrack; //total number of charged tracks;
    //NTuple::Item<long>    nntrack; //total number of neutral tracks;

    NTuple::Item<long>    Nn;    //total number of neutral tracks
    NTuple::Item<long>    NnE10; //number of neutral tracks with energy more then 10 MeV
    NTuple::Item<long>    NnE25; //number of neutral tracks with energy more then 25 MeV
    NTuple::Item<long>    NnE50; //number of neutral tracks with energy more then 50 MeV
    NTuple::Item<long>    NnE100; //number of neutral tracks with energy more then 100 MeV
    NTuple::Item<long>    Nng; //number of good neutral tracks

    NTuple::Item<long>    Nc;  //Total number of charged tracks
    NTuple::Item<long>    Ncc; //Number of central charged tracks
    NTuple::Item<long>    Nce; //Number of central charged tracks with emc
    NTuple::Item<long>    Ncg; //Number of good central charged tracks with emc


    NTuple::Item<double>  enmin; //min energy of all neutral tracks
    NTuple::Item<double>  enmax; //max energy of all neutral tracks
    NTuple::Item<double>  entot; //total energy of all neutral tracks



    NTuple::Item<long> channel; //used channel
    //NTuple::Item<long> ntrack;  //number of tracks (must be 2)
    RootTracks T;  //track information (momentum, vertex, muon depth...)
    RootTracks Tn; //neutral track
    RootPid Pid;   //particle id for charged track
    //NTuple::Array<double> pid; //particle id
    //NTuple::Array<double> mother_pid;
    RootMcTruth  McTruth; // mc truth
    //RootMcTruth nMcTruth; // mc truth for nuetral tracks
    NTuple::Item<double>  acop;
    NTuple::Item<double>  acol;
    NTuple::Item<double>  cos_theta_mis;
    NTuple::Item<double>  M2;
    NTuple::Item<double>  S; //sphericity (lambda2+lambda3)*1.5
    NTuple::Item<double>  A; //aplanarity 1.5*lambda3
    NTuple::Item<double>  lambda1; //eigen values of sphericity tenzor
    NTuple::Item<double>  lambda2; //
    NTuple::Item<double>  lambda3; //

    NTuple::Item<double>  Emis;  //Wcm  - sum |p_i|  - sum En 
    NTuple::Item<double>  Emis0; //Wcm - sum |p_i| //no neutral tracks in Emis
    NTuple::Item<double>  Emis1; //Wcm - p_0 - p_1 - sum En 
    NTuple::Item<double>  Emis2; //Wcm - p_0 - p_1 


    NTuple::Item<double>  psum; //Total momentum in the event
    NTuple::Item<double>  psum2; //Momentum of two tracks 

    NTuple::Item<double>  ptsum; //sum pt_i  - for all Ncg charged tracks
    NTuple::Item<double>  ptsum2; //|pt_0 + pt_1|

    NTuple::Item<double>  ptem; // ptsum/Emis total
    NTuple::Item<double>  ptem10; // ptsum/Emis10
    NTuple::Item<double>  ptem25; // ptsum/Emis25
    NTuple::Item<double>  ptem50; // ptsum/Emis50


    NTuple::Item<double>  Enmin; //minimum energy of accepted neutral tracks
    NTuple::Item<double>  Enmax; //maximum energy of accepted neutral tracks
    NTuple::Item<double>  Entot; //total energy of the acceptd neutral tracks
    NTuple::Item<double>  Entot10; //total energy of photons with E>10 MeV
    NTuple::Item<double>  Entot25; //total energy of photons with E>25 MeV
    NTuple::Item<double>  Entot50; //total energy of photons with E>50 MeV
    NTuple::Item<double>  Entot100; //total energy of photons with E>100 MeV


    NTuple::Item<long>    Npi0; //number of pi0
    NTuple::Array<double> Mpi0;

    NTuple::Item<long>    Nrho; //number of rho candidates
    NTuple::Array<double> Mrho; //this mass of rho candidates
    NTuple::Array<double> IdxRho; //Indexes in the charged particles T, which used for Rho

    NTuple::Item<long>    ndp; // ndp=2 (Energy and 3d-momentum)

    std::vector< NTuple::Array<double> *> d4p; 
    std::vector< NTuple::Array<double> *> d4pg; 
    NTuple::Array<double> P4; //total four momentum of all charged

    virtual void bind_tuple(void)
    {
      tuple->addItem ("run", run);
      tuple->addItem ("event", event);
      tuple->addItem ("time", time);
      tuple->addItem ("channel", channel);
      //tuple->addItem ("nctrack",   nctrack); 
      //tuple->addItem ("nntrack",   nntrack); 

      tuple->addItem ("enmin",enmin);
      tuple->addItem ("enmax",enmax);
      tuple->addItem ("entot",entot);

      tuple->addItem ("Nc", Nc); //total number of charged tracks
      tuple->addItem ("Ncc", Ncc); //central charged tracks
      tuple->addItem ("Nce", Nce); //central charged tracks with EMC
      tuple->addItem ("Ncg", Ncg, 0, MAX_CHARGED_TRACKS); //Good charged tracks with EMC

      tuple->addItem ("Nn", Nn, 0, MAX_NEUTRAL_TRACKS); //total number of neutral tracks (raw)
      tuple->addItem ("NnE10",  NnE10);
      tuple->addItem ("NnE25",  NnE25);
      tuple->addItem ("NnE50",  NnE50);
      tuple->addItem ("NnE100", NnE100);
      tuple->addItem ("Nng",    Nng);


      tuple->addItem ("Npi0", Npi0, 0, MAX_PI0_NUMBER);
      tuple->addItem ("Nrho", Nrho, 0, MAX_RHO_NUMBER);


      tuple->addItem ("Enmin",Enmin);
      tuple->addItem ("Enmax",Enmax);
      tuple->addItem ("Entot",Entot);
      tuple->addItem ("Entot10",Entot10);
      tuple->addItem ("Entot25",Entot25);
      tuple->addItem ("Entot50",Entot50);
      tuple->addItem ("Entot100",Entot100);

      tuple->addItem("Emis",Emis);
      tuple->addItem("Emis0",Emis0);
      tuple->addItem("Emis1",Emis1);
      tuple->addItem("Emis2",Emis2);

      tuple->addItem("psum",psum);
      tuple->addItem("psum2",psum2);

      tuple->addItem("ptsum",ptsum);
      tuple->addItem("ptsum2",ptsum2);

      tuple->addItem("ptem",ptem);
      tuple->addItem("ptem10",ptem10);
      tuple->addItem("ptem25",ptem25);
      tuple->addItem("ptem50",ptem50);


      T.add_to_tuple (tuple,Ncg); 
      Tn.add_to_tuple(tuple,Nn,"n");
      Pid.add_to_tuple(tuple,Ncg); 
      McTruth.add_to_tuple(tuple,Ncg);

      //nMcTruth.add_to_tuple(tuple,ntrack,"n");
      tuple->addItem("acop",acop);
      tuple->addItem("acol",acol);
      tuple->addItem("cos_theta_mis", cos_theta_mis);
      tuple->addItem("M2",M2);
      tuple->addItem("S",S);
      tuple->addItem("A",A);
      tuple->addItem("l1",lambda1);
      tuple->addItem("l2",lambda2);
      tuple->addItem("l3",lambda3);
      tuple->addIndexedItem("Mpi0", Npi0, Mpi0);
      tuple->addIndexedItem("Mrho", Nrho, Mrho);
      tuple->addIndexedItem("IdxRho", Nrho, IdxRho);

      d4p.resize(4); //for  0, 1,2,3 photons
      d4pg.resize(4); //for 

      for(size_t i = 0;i< d4p.size(); ++i) {
        d4p[i] = new NTuple::Array<double>(); 
      }
      for(size_t i = 0;i< d4pg.size(); ++i) {
        d4pg[i] = new NTuple::Array<double>(); 
      }

      tuple->addItem ("ndp", ndp, 0, 6); 
      // 0 - energy
      // 1 - momentum
      // 2 - transverse momentum
      // 3 - theta
      // 4 - phi
      // 5 - m2
      tuple->addIndexedItem("Psum",ndp, P4);
      char buf[1024];
      for(size_t i=0; i< d4p.size(); ++i) {
        sprintf(buf,"d4p%d",i);
        tuple->addIndexedItem(buf, ndp, *(d4p[i]));
      }
      for(size_t i=0; i< d4pg.size(); ++i) {
        sprintf(buf,"d4pg%d",i);
        tuple->addIndexedItem(buf, ndp, *(d4pg[i]));
      }
    };
    virtual void init(void) {};
    virtual void fill(int i,  EvtRecTrack * track);
    virtual void nfill(int i,  EvtRecTrack * track);

    void fillPi0Rho(const comb_t & pi0s, std::vector<EvtRecTrack*> & );

    //bool pass(const  std::vector<EvtRecTrack*>  & Tc, const  std::vector<EvtRecTrack*>  & Tn, const SelectionConfig & cfg);
    //bool pass(const SelectionConfig & cfg, const Event::EventHeader * eventHeader, const Event::McParticleCol *,  const  std::vector<EvtRecTrack*>  & Tc, const  std::vector<EvtRecTrack*>  & Tn, const  std::vector<EvtRecTrack*>  & Tgn);
      bool pass(const SelectionConfig & cfg, const Event::EventHeader *  eventHeader, const Event::McParticleCol * mcParticleCol,  const  Tracker & tracker);
};

