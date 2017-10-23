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
#include "DstEvent/TofHitStatus.h"

enum PID 
{
  ELECTRON=0, 
  MUON=1, 
  PION=2, 
  KAON=3, 
  PROTON=4
};

struct RootPid
{
	NTuple::Item<long> * ntrack; 
	NTuple::Array<double> chi2_pid[5];  //this if from package PID
	NTuple::Array<double> chi2_dedx[5]; 
	NTuple::Array<double> chi2_tof[5]; 

  ParticleID * PID;

	virtual void add_to_tuple(NTuple::Tuple * tuple, NTuple::Item<long> & ntrk)
  {
    ntrack = & ntrk;
    const char * chan[] = {"e","mu","pi","K","p"};
    for( int pid = 0;pid< 5;++pid)
    {
      char item_name[1024];
      sprintf(item_name,"chi2_pid_%s",chan[pid]);
      tuple->addIndexedItem (item_name,  *ntrack, chi2_pid[pid]);

      sprintf(item_name,"chi2_dedx_%s",chan[pid]);
      tuple->addIndexedItem (item_name,  *ntrack, chi2_dedx[pid]);

      sprintf(item_name,"chi2_tof_%s", chan[pid]);
      tuple->addIndexedItem (item_name,  *ntrack, chi2_tof[pid]);
    }
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
      chi2_pid[pid][i]= PID->chi(pid);
    }
    //my pid information
    if(track->isMdcTrackValid() && track->isMdcDedxValid() && track->isTofTrackValid())
    {
      RecMdcDedx  * dedx = track->mdcDedx();
      chi2_dedx[ELECTRON][i] = chi2_dedx[ELECTRON][i] + pow(dedx->chiE() ,2);
      chi2_dedx[MUON][i]     = chi2_dedx[MUON][i]     + pow(dedx->chiMu(),2);
      chi2_dedx[PION][i]     = chi2_dedx[PION][i]     + pow(dedx->chiPi(),2);
      chi2_dedx[KAON][i]     = chi2_dedx[KAON][i]     + pow(dedx->chiK() ,2);
      chi2_dedx[PROTON][i]   = chi2_dedx[PROTON][i]   + pow(dedx->chiP() ,2);

      SmartRefVector<RecTofTrack> tofs = track->tofTrack();
      double  atof[5] = {0,0,0,0,0};
      double dtof2[5] = {0,0,0,0,0};
      unsigned nlayers = 0;
      for(SmartRefVector<RecTofTrack>::iterator iter_tof = tofs.begin();
          iter_tof!=tofs.end(); ++iter_tof)
      {
        TofHitStatus *status = new TofHitStatus;
        status->setStatus((*iter_tof)->status());
        double tof  = (*iter_tof)->tof();
        double texp[5];
        texp[ELECTRON] = (*iter_tof)->texpElectron();
        texp[MUON]     = (*iter_tof)->texpMuon();
        texp[PION]     = (*iter_tof)->texpPion();
        texp[KAON]     = (*iter_tof)->texpKaon();
        texp[PROTON]   = (*iter_tof)->texpProton();
        for(int pid = 0; pid < 5; ++pid)
        {
          atof[pid] += tof - texp[pid];
          dtof2[pid] += (*iter_tof)->errtof() * (*iter_tof)->errtof();
        }
        nlayers++;
      }

      for(int pid = 0; pid < 5; ++pid)
      {
        atof[pid] /= nlayers;
        dtof2[pid] /= nlayers;
        chi2_tof[pid][i] = atof[pid] * atof[pid] / dtof2[pid];
      }
    }
  }
};
