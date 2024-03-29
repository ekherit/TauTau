/*
 * =====================================================================================
 *
 *       Filename:  JPsi.h
 *
 *    Description:  Multihadron event selection for j/psi and psi prime resonance.
 *
 *        Version:  1.0
 *        Created:  04/27/2010 02:47:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */
#pragma once

#ifndef IBN_TAUTAU_H
#define IBN_TAUTAU_H
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"

#include "SelectionConfig.h"

#include "JobInfo.h"
#include "TauTauEvent.h" 
#include "GammaGammaEvent.h"
#include "BhabhaEvent.h"

class TauTau : public Algorithm 
{
  public:
    TauTau(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();  
  private:
    SelectionConfig cfg;

    //counters
    long int nproceed_events;  //total event proceed
    long int nwritten_events;  //number of written events into all tuples
    long int ntautau_events;   //number of selected (and written) tau tau events
    long int nbhabha_events;   //number of Bhabha events to measure luminosity
    long int ngg_events;       //number of Digamma events to measure luminosity

    //storages with root tuple and selection code
    TauTauEvent     fTT;       //tau tau events
    GammaGammaEvent fGG;       //gamma gamma events for luminosity
    BhabhaEvent     fBB;       //Bhabha events for luminosity
    JobInfo         fJobInfo;     //job information
};

#endif
