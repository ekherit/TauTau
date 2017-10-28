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
	NTuple::Array<double> fDedx[5]; 

	NTuple::Array<double> ftof; 
	NTuple::Array<double> tof_sigma; 
	NTuple::Array<double> tof_exp[5]; 
	NTuple::Array<double> tof_delta[5]; 

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

      sprintf(item_name,"tof_exp_%s", chan[pid]);
      tuple->addIndexedItem (item_name,  *ntrack, tof_exp[pid]);

      sprintf(item_name,"delta_tof_%s", chan[pid]);
      tuple->addIndexedItem (item_name,  *ntrack, tof_delta[pid]);

      sprintf(item_name,"dedx_%s", chan[pid]);
      tuple->addIndexedItem (item_name,  *ntrack, fDedx[pid]);
    }
    tuple->addIndexedItem ("tof",  *ntrack, ftof);
    tuple->addIndexedItem ("tof_sigma",  *ntrack, tof_sigma);
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
      fDedx[ELECTRON][i] = dedx->chiE() ;
      fDedx[MUON][i]     = dedx->chiMu();
      fDedx[PION][i]     = dedx->chiPi();
      fDedx[KAON][i]     = dedx->chiK() ;
      fDedx[PROTON][i]   = dedx->chiP() ;

      chi2_dedx[ELECTRON][i] = pow(dedx->chiE() ,2.0);
      chi2_dedx[MUON][i]     = pow(dedx->chiMu(),2.0);
      chi2_dedx[PION][i]     = pow(dedx->chiPi(),2.0);
      chi2_dedx[KAON][i]     = pow(dedx->chiK() ,2.0);
      chi2_dedx[PROTON][i]   = pow(dedx->chiP() ,2.0);

      SmartRefVector<RecTofTrack> tofs = track->tofTrack();
      double tof_sum=0;
      double dtof2_sum=0;//sum of the errors in the tof
      double dtof_sum[5] = {0,0,0,0,0}; //sum of the differencies in tof1 and tof2
      double tof_exp_sum[5] = {0,0,0,0,0};
      unsigned nlayers = 0;
      for(SmartRefVector<RecTofTrack>::iterator iter_tof = tofs.begin();
          iter_tof!=tofs.end(); ++iter_tof)
      {
        TofHitStatus *status = new TofHitStatus;
        status->setStatus((*iter_tof)->status());
        if( !(status->is_counter()) ) continue;
        //if( status->layer()!=0 ) continue;

        double t  = (*iter_tof)->tof();
        ftof[i] = t;
        tof_sum += t;
        double texp[5];
        texp[ELECTRON] = (*iter_tof)->texpElectron();
        texp[MUON]     = (*iter_tof)->texpMuon();
        texp[PION]     = (*iter_tof)->texpPion();
        texp[KAON]     = (*iter_tof)->texpKaon();
        texp[PROTON]   = (*iter_tof)->texpProton();

        tof_exp[ELECTRON][i] = (*iter_tof)->texpElectron();
        tof_exp[MUON][i]     = (*iter_tof)->texpMuon();
        tof_exp[PION][i]     = (*iter_tof)->texpPion();
        tof_exp[KAON][i]     = (*iter_tof)->texpKaon();
        tof_exp[PROTON][i]   = (*iter_tof)->texpProton();

        dtof2_sum += (*iter_tof)->errtof() * (*iter_tof)->errtof();
        for(int pid = 0; pid < 5; ++pid) 
        {
          tof_exp_sum[pid] += texp[pid];
          dtof_sum[pid] += (t - texp[pid]);
        }

        nlayers++;
        //std::cout <<  t << "  " << (*iter_tof)->errtof() << "  barrel=" <<  status->is_barrel();
        //std::cout << "  counter=" << status->is_counter();
        //std::cout << " layer = " <<  status->layer();
        //std::cout <<std::endl;
      }
      if(nlayers>0)
      {
        //std::cout <<"new track" << std::endl;
        ftof[i] = tof_sum/nlayers;
        tof_sigma[i] = sqrt(dtof2_sum/nlayers);
        //std::cout << ftof[i] << "  " << tof_sigma[i] << std::endl;
        for(int pid = 0; pid < 5; ++pid)
        {
          tof_exp[pid][i] = tof_exp_sum[pid]/nlayers;
          tof_delta[pid][i] = dtof_sum[pid]/nlayers;
          chi2_tof[pid][i] = pow ( tof_delta[pid][i]/tof_sigma[i], 2); ;
        }
      }
      else
      {
        ftof[i] = 99;
        tof_sigma[i] = 99;
        for(int pid = 0; pid < 5; ++pid)
        {
          tof_exp[pid][i] = 99; 
          tof_delta[pid][i] = 99; 
          chi2_tof[pid][i] = 999;
        }
      }

    }
    //std::cout << "new track" << std::endl;
  }
};
