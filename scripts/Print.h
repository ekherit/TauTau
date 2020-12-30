/*
 * =====================================================================================
 *
 *       Filename:  Print.h
 *
 *    Description: Some printer functions 
 *
 *        Version:  1.0
 *        Created:  30.12.2020 11:43:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#ifndef IBN_TAU_PRINT_H
#define IBN_TAU_PRINT_H
#include <string>
#include <vector>
#include <iostream>
#include <regex>
#include <list>

#include "ScanPoint.h"

class printer
{
  struct print_cfg_t
  {
    std::string head;
    std::string fmt;
  };
  std::vector<print_cfg_t> formats;
  std::string result;
  std::string Head;
  char buf[65535];
  int idx;
  public:
  printer(void) {idx=0;}
  printer & add(std::string h, std::string f) 
  {
    formats.push_back({h,f});
    static const std::regex re(R"(%([^\d]?\d+)(?:\.\d+)?[[:alpha:]])");
    std::smatch sm;
    int size;
    std::string head_format;
    if(std::regex_search(f,sm,re)) head_format = "%"+std::string(sm[1])+"s";
    else std::cerr << "ERROR: bad print format for printer" << std::endl;
    sprintf(buf,head_format.c_str(),h.c_str());
    Head+=buf;
    return *this;
  };
  printer & add(std::string f) 
  {
    return add("",f);
  };
  printer & operator()(std::string h, std::string f) { return add(h,f); }
  template <typename Data>
  printer & operator%(const Data  & d) 
  {
    if(idx==0) result="";
    sprintf(buf,formats[idx].fmt.c_str(),d);
    result+=buf;
    ++idx%=formats.size();
    return *this;
  }
  void print(std::ostream & os)
  {
    os << result;
  }
  const std::string  & head(void)  const { return Head; }
};

std::ostream & operator<<(std::ostream & os, printer & p)
{
  p.print(os);
  return os;
};

std::string get_run_formula(const std::list<int> & runs)
{
  if(runs.empty()) return "";
  std::string formula;
  for ( auto r = std::begin(runs); r!=std::end(runs); ++r)
  {
    if( r == std::begin(runs) ) formula = std::to_string(*r); //start sequence
    else
    {
      if ( *r == *std::prev(r) + 1 )  //different in run number is 1
      {
        if( r == std::prev(std::end(runs)) || *std::next(r) != *r + 1 )  //this is last item or next item is not next run and we need to finish sequency
        {
          formula +="-"+std::to_string(*r); 
        }
        continue;
      }
      else
      {
        formula+=" "+std::to_string(*r);
      }
    }
  }
  return formula;
}

void print(const std::vector<ScanPoint_t> & SPL)
{
  int point_number = 1;
  printer pr;
  //pr.add("#point",     "%-7d");
  pr.add("#name",       "%-10s");
  pr.add("W,GeV",      "%15.6f");
  pr.add("dW,GeV",     "%15.6f");
//  pr.add("E-Mtau,MeV", "%15.3f");
  pr.add("Sw,GeV",     "%15.6f");
  pr.add("dSw,GeV",    "%15.6f");
  pr.add("L,pb^-1",    "%15.6f");
  pr.add("dL,pb^-1",   "%15.6f");
  pr("  run list", "  %-10s");
  std::cout << pr.head() << std::endl;
  for (const auto & p : SPL) 
    std::cout <<  pr 
   //   % point_number++ 
      % p.title.c_str() 
      % p.energy.value % p.energy.error 
//      % ((0.5*p.energy.value-MTAU)/MeV) 
      % p.energy_spread.value
      % p.energy_spread.error
      % p.luminosity.value 
      % p.luminosity.error
      % get_run_formula(p.run_list).c_str() << std::endl;
}

#endif
