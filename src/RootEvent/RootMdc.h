// =====================================================================================
//
//       Filename:  RootMCTopo.h
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
#include "RootTrack.h"
struct RootMdc : public RootTracks, public RootTuple
{
	virtual void init(void);
	virtual void init_tuple(void);
	virtual void fill(int i,  EvtRecTrack * track);
};
