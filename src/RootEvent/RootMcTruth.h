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
#include <set>

#include "RootTuple.h"
#include "McParticle.h"

//#include "CLHEP/Vector/LorentzVector.h"

struct RootMcTruth
{
  NTuple::Item<long> flag1; //tag the number of root particle decay mode in *.dec file
  NTuple::Item<long> flag2; //tag the type of the root particle decay
	NTuple::Array<double> pid; //particle id
	NTuple::Array<double> mother_pid;
  std::set<int> blacklist;
  std::set<int> whitelist;
	virtual void add_to_tuple(NTuple::Tuple * tuple, NTuple::Item<long> & ntrk, const std::string & prefix="")
  {
    //ntrack = & ntrk;
    tuple->addItem ("flag1", flag1);
    tuple->addItem ("flag2", flag2);
    tuple->addIndexedItem (prefix+"pid",  ntrk, pid);
    tuple->addIndexedItem (prefix+"mother_pid", ntrk, mother_pid);
    //add black item for McTruth
    //remove quarks
    add_black_item(1); //d
    add_black_item(-1); //dbar
    add_black_item(2); //u
    add_black_item(-2); //ubar
    add_black_item(3); //s
    add_black_item(-3); //sbar
    add_black_item(4); //c
    add_black_item(-4); //cbar
    add_black_item(5); //b
    add_black_item(-5); //bbar
    add_black_item(6); //t
    add_black_item(-6); //tbar
    add_black_item(7); //bp
    add_black_item(-7); //bpbar
    add_black_item(8); //tp
    add_black_item(-8); //tpbar

    add_black_item(12); //nu_e
    add_black_item(-12); //anti nu_e
    add_black_item(14); //nu_mu
    add_black_item(-14); //anti nu_mu
    add_black_item(16); //nu_tau
    add_black_item(-16); //anti nu_tau
    add_black_item(18); //nu_taup
    add_black_item(-18); //anti nu_taup

    add_black_item(21); //gluon
    add_black_item(-21); //gluon
  }

	virtual void init(void){};

  void fill(int i, EvtRecTrack * track)
  {
  }

  void set_blacklist(std::set<int> bl) 
  {
    blacklist = bl;
  }
  void set_whitelist(std::set<int> bl) 
  {
    whitelist = bl;
  }

  void add_black_item(int id)
  {
    blacklist.insert(id);
  }
  void add_white_item(int id)
  {
    whitelist.insert(id);
  }

  virtual void fill(int i, EvtRecTrack * track,  const Event::McParticleCol * mcParticleCol)
  {
    RecMdcKalTrack * mdc = track->mdcKalTrack();
    double x=1e100;
    bool found=false;
    for(Event::McParticleCol::const_iterator ip=mcParticleCol->begin(); ip!=mcParticleCol->end(); ++ip)
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
      }
    }
    if (!found) 
    {
      pid[i] = 0;
      mother_pid[i] = 0;
    }
  }
};
