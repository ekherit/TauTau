// =====================================================================================
//
//       Filename:  RootMcTruth.h
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

#include "RootTuple.h"
#include "McParticle.h"

//#include "CLHEP/Vector/LorentzVector.h"

struct RootMcTruth
{
	NTuple::Item<long> * ntrack; 
	NTuple::Array<long> pid; //particle id
	NTuple::Array<long> mother_pid;
  std::list<int> blacklist;
  std::list<int> whitelist;
	virtual void add_to_tuple(NTuple::Tuple * tuple, NTuple::Item<long> & ntrk, const std::string & prefix="")
  {
    ntrack = & ntrk;
    tuple->addIndexedItem (prefix+"pid",  *ntrack, pid);
    tuple->addIndexedItem (prefix+"mpid",  *ntrack, mother_pid);
  }

	virtual void init(void){};

  void fill(int i, EvtRecTrack * track)
  {
  }

  void set_blacklist(std::list<int> bl) 
  {
    blacklist = bl;
  }
  void set_whitelist(std::list<int> bl) 
  {
    whitelist = bl;
  }

  void add_black_item(int id)
  {
    blacklist.push_back(id);
  }
  void add_white_item(int id)
  {
    whitelist.push_back(id);
  }

	virtual void fill(int i, EvtRecTrack * track,  Event::McParticleCol * mcParticleCol)
  {
    RecMdcKalTrack * mdc = track->mdcKalTrack();
    double x=1e100;
    bool found=false;
    for(Event::McParticleCol::iterator ip=mcParticleCol->begin(); ip!=mcParticleCol->end(); ++ip)
    {
      Event::McParticle * p = *ip;
      int mc_track_id = p->trackIndex();
      int pdg_id  = p->particleProperty();
      if(blacklist.find(pdg_id) != blacklist.end())  continue; //skip some particles
      //or keep specific particles if whitelist is not empty
      if( !whitelist.empty() && whitelist.find(pdg_id)==whitelist.end() ) continue; 
      //HepLorentzVector P4 = p->initialFourMomentum();
      Hep3Vector p_mc = p->initialFourMomentum().vect();
      Hep3Vector delta_p = mdc->p3() - p_mc;
      double y = sqrt( delta_p.mag2() / ( p_mc.mag() * mdc->p3().mag() ) );
      if(y<x && y<0.1)
      {
        found=true;
        x=y;
        pid[i] = p->particleProperty();
        mother_pid[i] = p->mother().particleProperty();
        //std::cout << y << " " << mdc->charge() << " " << pid [i] << " " << mdc->p() << " " << p_mc.mag()  << std::endl;
        //break;
      }
    }
    if (!found) 
    {
      pid[i] = 0;
      mother_pid[i] = 0;
    }
  }
};
