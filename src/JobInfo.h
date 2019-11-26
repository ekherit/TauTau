/*
 * =====================================================================================
 *
 *       Filename:  Info.h
 *
 *    Description:  Job information
 *
 *        Version:  1.0
 *        Created:  26.11.2019 11:38:19
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#pragma once
#include "RootEvent/RootTuple.h" 

class JobInfo : public RootTuple
{
  public:
    virtual ~JobInfo(void) {};
    JobInfo(void) { }
    NTuple::Item<long> N;             //total number of event proceed
    NTuple::Item<long> n;             //writen number of events
    NTuple::Item<long> begin_time;    //job beging unixtime
    NTuple::Item<long> end_time;      //job end unixtime
    NTuple::Item<long> type;          //job type 0 - data, 1 - MC

    virtual void bind_tuple(void)
    {
      tuple->addItem ("N", N);
      tuple->addItem ("n", n);
      tuple->addItem ("begin_time", begin_time);
      tuple->addItem ("end_time", end_time);
      tuple->addItem ("type",   type);
    };
    virtual void init(void) { };
};
