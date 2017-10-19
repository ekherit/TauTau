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

struct RootTrack
{
  static long MAX_TRACK_NUMBER;
  static long MIN_TRACK_NUMBER;
	NTuple::Item<long>    ntrack;  //size of 
	NTuple::Array<long>   trackId; //id of the track
	NTuple::Array<double> q; //charge of the track
	NTuple::Array<double> E;
	NTuple::Array<double> p;
	NTuple::Array<double> px;
	NTuple::Array<double> py;
	NTuple::Array<double> pz;
	NTuple::Array<double> pt; //transvese momentum
	NTuple::Array<double> theta,phi;
	NTuple::Array<double> x, y, z, r; //poca coordinate of track
	NTuple::Array<double> vxy, vz, vphi; //poca coordinate of track
  NTuple::Array<double> depth; //depth in muon system
	virtual void add_to_tuple(NTuple::Tuple * tuple)
	{
		//tuple->addItem ("ntrack", ntrack, MIN_TRACK_NUMBER,MAX_TRACK_NUMBER); 
		tuple->addItem ("ntrack", ntrack); 
		tuple->addIndexedItem ("trackId",   ntrack, trackId);
		tuple->addIndexedItem ("q",     ntrack, q);
		tuple->addIndexedItem ("E",     ntrack, E);
		tuple->addIndexedItem ("p",     ntrack, p);
		tuple->addIndexedItem ("px",    ntrack, px);
		tuple->addIndexedItem ("py",    ntrack, py);
		tuple->addIndexedItem ("pz",    ntrack, pz);
		tuple->addIndexedItem ("pt",    ntrack, pt);
		tuple->addIndexedItem ("theta", ntrack, theta);
		tuple->addIndexedItem ("phi",   ntrack, phi);
		tuple->addIndexedItem ("x",     ntrack, x);
		tuple->addIndexedItem ("y",     ntrack, y);
		tuple->addIndexedItem ("z",     ntrack, z);
		tuple->addIndexedItem ("r",     ntrack, r);
		tuple->addIndexedItem ("vxy",   ntrack, vxy);
		tuple->addIndexedItem ("vz",    ntrack, vz);
		tuple->addIndexedItem ("vphi",  ntrack, vphi);
		tuple->addIndexedItem ("depth",  ntrack, depth);
	}
};
