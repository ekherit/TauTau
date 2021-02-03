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
#include <algorithm>

#include <TTree.h>

#include "../ibn/valer.h"
#include "utils.h"


struct AcceleratorInfo_t 
{
  ibn::valer<double> energy{0,0};        //beam c.m. energy
  ibn::valer<double> energy_spread{0,0}; //beam energy spread
  ibn::valer<double> luminosity{0,0};    //accelerator luminosity
};

struct DataSample_t  : public AcceleratorInfo_t
{
  std::string title;
  std::shared_ptr<TTree> tree;
  ibn::valer<double> cross_section{0,0};
  ibn::valer<double> efficiency{0,0}; //registration efficiency
  ibn::valer<double> effcor{0,0}; //correction to efficiency
  long N{0}; //number of selected events in data
  long Nmc{0}; //number of selected events in MC
  long N0mc{0}; //initial number of MC events

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
  
  DataSample_t operator+=(const DataSample_t & d) {
    energy += d.energy;
    energy_spread += d.energy_spread;
    luminosity += luminosity;
    N += d.N ;
    Nmc += d.Nmc;
    N0mc += d.N0mc;
    //d.cross_section = d1.cross_section
    efficiency += d.efficiency;
    effcor += d.effcor; 
    //if(!tree)  tree.reset( new TChain ( d.tree->GetName(), d.tree->GetTitle()));

    //if(std::string(tree->GetName()) != std::string(tree->GetName())) {
    //  throw std::runtime_error("Unable to add different chains");
    //};
    //dynamic_pointer_cast<TChain>(tree)->Add( dynamic_cast<TChain*>(d.tree.get()) );
    return d;
  };
};

inline DataSample_t operator+(const DataSample_t & d1, const DataSample_t & d2) {
  DataSample_t d;
  d+=d1;
  d+=d2;
  //d.N = d1.N + d2.N;
  //d.Nmc = d1.Nmc + d2.Nmc;
  //d.N0mc = d1.N0mc + d2.N0mc;
  ////d.cross_section = d1.cross_section
  //d.efficiency = d1.efficiency + d2.efficiency;
  //d.effcor = d1.effcor + d2.effcor;
  //if(std::string(d1.tree->GetName()) != std::string(d2.tree->GetName())) {
  //  throw std::runtime_error("Unable to add different chains");
  //};
  //d.tree.reset( new TChain ( d1.tree->GetName(), d1.tree->GetTitle()));
  //dynamic_pointer_cast<TChain>(d.tree)->Add( dynamic_cast<TChain*>(d1.tree.get()) );
  return d;
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
    long N = ds.GetEntries(sel);
    ds.N = N;
    if(result.luminosity.value < 0) { //Monte Carlo case
      ds.Nmc = N;
      auto & eps = ds.efficiency;
      double N0 = ds.N0mc;
      ds.efficiency.value = double(N)/N0;
      ds.efficiency.error = sqrt( eps.value/N0 * ( 1.0 - eps.value));
    } 
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

  ScanPoint_t & operator+=(const ScanPoint_t & sp) {
    std::copy(sp.run_list.begin(), sp.run_list.end(), std::back_inserter(run_list));
    std::copy(sp.file_list.begin(), sp.file_list.end(), std::back_inserter(file_list));
    tt += sp.tt;
    gg += sp.gg;
    bb += sp.bb;
    return *this;
  };

};

inline ScanPoint_t operator+(const ScanPoint_t & sp1, const ScanPoint_t & sp2) {
  ScanPoint_t sp;
  std::merge( sp1.run_list.begin(), sp1.run_list.end(), 
              sp2.run_list.begin(), sp2.run_list.end(),
              std::back_inserter(sp.run_list)
      );

  std::merge( sp1.file_list.begin(), sp1.file_list.end(), 
              sp2.file_list.begin(), sp2.file_list.end(),
              std::back_inserter(sp.file_list)
      );
  sp.tt = sp1.tt + sp2.tt;
  sp.gg = sp1.gg + sp2.gg;
  sp.bb = sp1.bb + sp2.bb;
  return sp;
};

template <typename Container >
inline ScanPoint_t fold(const Container & cont) {
  static_assert(  std::is_same_v< typename Container::value_type,  ScanPoint_t > );
  ScanPoint_t result;
  for ( const auto & sp:  cont ){
    result += sp;
  };
  return result;
};


using PointSelectionResult_t = ScanPoint_t;

#endif
