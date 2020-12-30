/*
 * =====================================================================================
 *
 *       Filename:  ScanPoint.h
 *
 *    Description: Discription of the scanpoint 
 *
 *        Version:  1.0
 *        Created:  30.12.2020 09:48:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#ifndef IBN_TAU_SCAN_POINT_H
#define IBN_TAU_SCAN_POINT_H
#include <string>
#include <memory>
#include <regex>
#include <list>
#include <vector>

#include <TTree.h>

#include "../ibn/valer.h"


struct AcceleratorInfo_t 
{
  ibn::valer<double> energy;        //beam c.m. energy
  ibn::valer<double> energy_spread; //beam energy spread
  ibn::valer<double> luminosity;    //accelerator luminosity
};

struct DataSample_t  : public AcceleratorInfo_t
{
  std::string title;
  std::shared_ptr<TTree> tree;
  ibn::valer<double> cross_section;
  ibn::valer<double> efficiency; //registration efficiency
  ibn::valer<double> effcor; //correction to efficiency
  long N; //number of selected events in data
  long Nmc; //number of selected events in MC
  long N0mc; //initial number of MC events
  long Draw(std::string varexp, std::string selection="", std::string option="", long nentries = std::numeric_limits<long>::max(), long firstentry=0) const {
    return tree->Draw(varexp.c_str(), selection.c_str(), option.c_str(), nentries, firstentry);
  }
};


struct ScanPoint_t : public AcceleratorInfo_t
{
  std::string title;
  std::list<int> run_list; //list of runs

  DataSample_t tt; //tau tau events. This is the signal or background
  DataSample_t bb; //Bhabha  luminosity
  DataSample_t gg; //Digamma luminosity

  std::list<std::pair<int,double> > runs; //for drawing luminosity per runs

  //std::string selection;
  //std::map<std::string, int> NttMap;
  std::list<std::string> file_list;
  std::list<std::string> regexprs; //regular expessions to match files
  std::string scan_title;


  //std::string cut;
  std::string name;
  std::string root_name;
  std::string tex_name;

  long Draw(std::string varexp, std::string selection="", std::string option="", long nentries = std::numeric_limits<long>::max(), long firstentry=0) const {
    return tt.tree->Draw(varexp.c_str(), selection.c_str(), option.c_str(), nentries, firstentry);
  }

  struct by_regexp
  {
    const Scan_t * scan;
    by_regexp(const Scan_t * scan) : scan(scan) {}
    std::string operator()(std::string file)
    {
      for ( auto & s : *scan )
      {
      }
      std::regex re(R"(\D+(\d+)\.root)");
      std::smatch match;
      std::string result;
      int run=0;
      if(std::regex_match(file, match,re))
      {
        if(match.size()>1)
        {
          run = std::stoi(match[1]);
          for( auto & s : *scan)
          {
            if(std::count(std::begin(s.run_list),std::end(s.run_list), run) > 0)
            {
              result=s.title;
            }
          }
        }
      }
      return result;
    }
  };
  std::string type;
};
typedef std::reference_wrapper<ScanPoint_t> ScanPointRef_t;
typedef std::vector<ScanPoint_t>  Scan_t;
typedef std::reference_wrapper<Scan_t> ScanRef_t;
typedef std::shared_ptr<Scan_t> ScanPtr_t;
using PointSelectionResult_t = ScanPoint_t;

#endif
