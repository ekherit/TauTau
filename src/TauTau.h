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

#ifndef IBN_TAUEMU_H
#define IBN_TAUEMU_H
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"

#include <TMatrixD.h>
#include <vector>
#include <algorithm>

#include "EventModel/EventHeader.h"

#include "ibn/averager.h"

//#include "RootEvent/RootMdc.h"

#include "SelectionConfig.h"

#include "RootTauTauEvent.h"
#include "GammaGammaEvent.h"

class TauTau : public Algorithm 
{
  public:
    TauTau(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();  
  private:
    SelectionConfig cfg;
    long int nproceed_events;
    long int nwrited_events;
    long int ntautau_events; //number of mu ivents for track #1
    long int nbhabha_events; 
    long int ngg_events; 
    RootTauTauEvent fEvent; //tau tau events
    GammaGammaEvent fGG;   //gamma gamma events for luminocity
};

#endif
