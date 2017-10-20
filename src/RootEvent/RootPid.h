/*
 * =====================================================================================
 *
 *       Filename:  RootPid.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  19.10.2017 16:15:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#pragma once

#include "GaudiKernel/NTuple.h"
#include "ParticleID/ParticleID.h"

struct RootPid
{
	NTuple::Array<long> * ntrack; 
	NTuple::Array<double> chi2[5]; 
	NTuple::Array<double> prob[5]; 

  ParticleID * PID;

	virtual void add_to_tuple(NTuple::Tuple * tuple, NTuple::Itel<long> & ntrk)
  {
    ntrack = & ntrk;
		tuple->addIndexedItem ("pid_chi2_e",  *ntrack, chi2[ELECTRON]);
		tuple->addIndexedItem ("pid_prob_e",  *ntrack, prob[ELECTRON]);

		tuple->addIndexedItem ("pid_chi2_mu",  *ntrack, chi2[MUON]);
		tuple->addIndexedItem ("pid_prob_mu",  *ntrack, prob[MUON]);

		tuple->addIndexedItem ("pid_chi2_pi",  *ntrack, chi2[PION]);
		tuple->addIndexedItem ("pid_prob_pi",  *ntrack, prob[PION]);

		tuple->addIndexedItem ("pid_chi2_K",  *ntrack, chi2[KAON]);
		tuple->addIndexedItem ("pid_prob_K",  *ntrack, prob[KAON]);

		tuple->addIndexedItem ("pid_chi2_p",  *ntrack, chi2[KAON]);
		tuple->addIndexedItem ("pid_prob_p",  *ntrack, prob[KAON]);
  }

  void init(void)
  {
    PID = ParticleID::instance();
		PID->init();
  }

  void fill(int i, EvtRecTrack * track)
  {
    PID->setRecTrack(track);
    PID->setMethod(PID->methodProbability());
    PID->setChiMinCut(4);
    PID->usePidSys
          (
             PID->useDedx() | 
             PID->useTof1() | 
             PID->useTof2() | 
             PID->useEmc()  | 
             PID->useMuc()
          );

    PID->identify(PID->all()); 
    PID->calculate();
    for(int pid = 0; pid <5; ++pid)
    {
      chi2[pid][i] = PID->chi(pid);
      prob[pid][i] = PID->prob(pid);
    }
  }
}
