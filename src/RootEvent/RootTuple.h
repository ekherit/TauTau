// =====================================================================================
//
//       Filename:  RootTuple.h
//
//    Description:  Interface class for root tuples
//
//        Version:  1.0
//        Created:  27.10.2015 15:24:06
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once 

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "EvtRecEvent/EvtRecTrack.h"

struct RootTuple
{
	public:
    RootTuple(void) {};
		NTuple::Tuple * tuple; //tuple
		virtual ~RootTuple(void){};
		virtual void init(void)=0;
		virtual void bind_tuple(void)=0;
    virtual void make_tuple(Algorithm * algo, const char * dir, const char * title)
    {
      NTuplePtr nt(algo->ntupleSvc(), dir);
      if(nt) 
      {
        tuple = nt;
      }
      else
      {
        tuple = algo->ntupleSvc()->book(dir, CLID_ColumnWiseTuple, title);
        if(tuple)
        {
          bind_tuple();
        }
        else
        {
          char buf[1024];
          sprintf(buf, "   Cannot book N-tuple: %s %s %x", dir, title, long(tuple));
          throw std::runtime_error(buf);
        }
      }
    }
		virtual void fill(EvtRecTrack * track) {};
    virtual void write(void) { tuple->write(); }
};
