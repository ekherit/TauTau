// =====================================================================================
//
//       Filename:  RootTrack
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

struct RootTracks
{
	NTuple::Array<double>   id; //id of the track
	NTuple::Array<double> q; //charge of the track
	NTuple::Array<double> E;
	NTuple::Array<double> Ep;//E/p ratio
	NTuple::Array<double> p;
	NTuple::Array<double> px;
	NTuple::Array<double> py;
	NTuple::Array<double> pz;
	NTuple::Array<double> pt;
	NTuple::Array<double> theta;
	NTuple::Array<double> phi;
	NTuple::Array<double> x, y, z, r; 
	NTuple::Array<double> vxy, vz, vphi; //poca coordinate of track
  NTuple::Array<double> depth; //depth in muon system
  NTuple::Array<double> Nmuhit; //number of muon hits
  NTuple::Array<double> temc; //time of EMC
	virtual void add_to_tuple(NTuple::Tuple * tuple, NTuple::Item<long> & ntrack, const std::string & prefix="")
	{
		//tuple->addItem ("ntrack", ntrack); 
		tuple->addIndexedItem (prefix+"trackid",   ntrack, id);
		tuple->addIndexedItem (prefix+"q",         ntrack, q);
		tuple->addIndexedItem (prefix+"E",         ntrack, E);
		tuple->addIndexedItem (prefix+"Ep",        ntrack, Ep);
		tuple->addIndexedItem (prefix+"p",         ntrack, p);
		tuple->addIndexedItem (prefix+"px",        ntrack, px);
		tuple->addIndexedItem (prefix+"py",        ntrack, py);
		tuple->addIndexedItem (prefix+"pz",        ntrack, pz);
		tuple->addIndexedItem (prefix+"pt",        ntrack, pt);
		tuple->addIndexedItem (prefix+"theta",     ntrack, theta);
		tuple->addIndexedItem (prefix+"phi",       ntrack, phi);
		tuple->addIndexedItem (prefix+"x",         ntrack, x);
		tuple->addIndexedItem (prefix+"y",         ntrack, y);
		tuple->addIndexedItem (prefix+"z",         ntrack, z);
		tuple->addIndexedItem (prefix+"r",         ntrack, r);
		tuple->addIndexedItem (prefix+"vxy",       ntrack, vxy);
		tuple->addIndexedItem (prefix+"vz",        ntrack, vz);
		tuple->addIndexedItem (prefix+"vphi",      ntrack, vphi);
		tuple->addIndexedItem (prefix+"depth",     ntrack, depth);
		tuple->addIndexedItem (prefix+"Nmuhit",    ntrack, Nmuhit);
		tuple->addIndexedItem (prefix+"temc",      ntrack, temc);
	}
};
