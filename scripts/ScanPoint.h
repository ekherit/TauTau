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
#include <string>

#include <TTree.h>

#include "../ibn/valer.h"
#include "utils.h"


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

  long GetEntries(std::string sel) const {
    auto ss = split(sel,"||");
    long N{0};
    for(auto s : split(sel, "||") ) {
      N+=tree->GetEntries(std::string(s).c_str());
    }
    return N;
  };

};

struct ScanPoint_t; 
typedef std::reference_wrapper<ScanPoint_t> ScanPointRef_t;
typedef std::vector<ScanPoint_t>  Scan_t;
typedef std::reference_wrapper<Scan_t> ScanRef_t;
typedef std::shared_ptr<Scan_t> ScanPtr_t;

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

  template<typename Projector>
  auto select(Projector  proj, std::string  sel) const ->  ScanPoint_t{
    ScanPoint_t result = *this; //copy data
    DataSample_t & ds = std::invoke(proj,result);
    ds.N = ds.GetEntries(sel);
    return result;
  };

  /*
  auto select(std::string  sel) const ->  ScanPoint_t {
    ScanPoint_t result = *this; //copy data
    result.tt.N = result.tt.GetEntries(sel);
    return result;
  };
  */

  template<typename ... Projs>
    void print_column(void) const {};

  template<typename Proj1, typename ... Projs>
    void print_column(Proj1 proj1,  Projs... ps) const
    {
      constexpr size_t Nbuf{1024};
      char buf[Nbuf];
      constexpr size_t Ns = sizeof...(ps);

      auto data = std::invoke(proj1,*this);
      if constexpr (std::is_same_v<double,decltype(data)> ){
        snprintf(buf,Nbuf, "%15.6f", data);
      }
      if constexpr (std::is_same_v<int,decltype(data)> || std::is_same_v<unsigned,decltype(data)> ){
        snprintf(buf,Nbuf, "%15d", data);
      }
      if constexpr (std::is_same_v<int,decltype(data)> || std::is_same_v<long,decltype(data)> ){
        snprintf(buf,Nbuf, "%15ld", data);
      }
      std::cout << buf;
      if constexpr (Ns>0) print_column(ps...);
    }

};

using PointSelectionResult_t = ScanPoint_t;

#endif
