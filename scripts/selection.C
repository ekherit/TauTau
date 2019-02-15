/*
 * =====================================================================================
 *
 *       Filename:  selection.C
 *
 *    Description:  Tau tau pair selection
 *
 *        Version:  1.0
 *        Created:  15.06.2018 09:26:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <regex>
#include <string>
#include <list>
#include <algorithm>

#include <variant>

#include "time.h"
#include "stdlib.h"

#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TTree.h>



//gSystem->AddLinkedLibs("-lfmt");

const double GeV=1.0;
const double MeV=1e-3*GeV;

const double MTAU=1776.86*MeV;
const double MPI=0.13957061*GeV; 
const double MPI0=0.134977*GeV;

/* ==================== WORKING WITH UTF-8 ================================================== */

template <class String > 
size_t count_utf8_symbols( String s )
{
  return std::count_if(begin(s),end(s),[](char c) { return (c & 0xc0) != 0x80; } );
}
template <class String > 
size_t count_utf8_extra_byte( String s )
{
  return s.size() - std::count_if(begin(s),end(s),[](char c) { return (c & 0xc0) != 0x80; } );
}
template< typename D>
void print_utf(int width, D d)
{
  std::cout << std::setw(width) << d;
};

template<>
void print_utf(int width, std::string s)
{

  std::cout << std::setw(width+count_utf8_extra_byte(std::string(s))) << s;
};

template<>
void print_utf(int width, const char * s)
{
  print_utf(width, std::string(s));
};

/* ========================================================================================== */

//Get recursive file list
std::list<std::string> get_recursive_file_list(std::string dirname)
{
  std::list<std::string> flst;
  TSystemDirectory dir(dirname.c_str(), dirname.c_str());
  for(auto  file: * dir.GetListOfFiles())
  {
    std::string file_name(file->GetName());
    if(file->IsFolder())
    {
      if(file_name != "." && file_name != "..")
      {
        flst.merge(get_recursive_file_list(dirname+"/"+file_name));
      }
    }
    else
    {
      flst.push_back(dirname+"/"+file_name);
    }
  }
  flst.sort();
  return flst;
}

std::list<std::string> filter_file_list(const std::list<std::string> LST, std::string regexpr=R"(.+\.root)")
{
  std::list<std::string> flst;
  std::regex file_re(regexpr); //regular expression to filter files
  std::smatch file_match;
  for( auto & file : LST)
  {
    if(std::regex_match (file,file_match,file_re))
    {
      flst.push_back(file);
    }
  }
  return flst;
}


template<class Discriminator>
std::map<std::string, std::list<std::string> > combine(const std::list<std::string> LST, Discriminator D)
{
  std::map<std::string, std::list<std::string> > fmap;
  for(const std::string & file : LST)
  {
    fmap[D(file)].push_back(file);
  }
  return fmap;
}

std::map<std::string, std::list<std::string> > combine(const std::list<std::string> LST, std::string regexpr=R"(^\D*(\d+\.?\d+).root$)")
{
  std::regex re(regexpr);
  return combine(LST, [&re](const std::string & file)
                          {  
                            std::string result;
                            std::smatch match;
                            if(std::regex_match(file, match,re))
                            {
                              if(match.size()>1)
                              {
                                result =  match[1];
                              }
                            }
                            return result;
                          }
                 );
};


//template < typename Map>
//void print(const Map & fmap)
void print(const std::map<std::string, std::list<std::string> > & fmap)
{
  for(const auto & point :  fmap)
  {
    for(const auto & file : point.second)
    {
      std::cout << point.first << ": " << file << std::endl;
    }
    std::cout << std::endl;
  }
}

TChain * get_chain(std::string name, std::string newname, std::string title, std::string filename)
{
  TChain * chain = new TChain(name.c_str(), title.c_str());
  chain->AddFile(filename.c_str());
  chain->SetName(newname.c_str());
  return chain;
}

TChain * get_chain(const char * name, const char * newname, std::string title, const char * file_name)
{
  TChain * chain = new TChain(name, title.c_str());
  chain->AddFile(file_name);
  chain->SetName(newname);
  return chain;
}

TChain * get_chain(const char * name, const char * newname, std::string title, int run_begin, int run_end, std::string dir="")
{
  TChain * chain = new TChain(name, title.c_str());
  for(int run = run_begin; run<=run_end; ++run)
  {
    char buf[1024];
    if(dir=="") sprintf(buf, "%d.root",run);
    else sprintf(buf, "%s/%d.root", dir.c_str(), run);
    chain->AddFile(buf);
  }
  chain->SetName(newname);
  return chain;
}


struct Restriction_t
{
  std::string name;
  double min;
  double max;
};

Restriction_t  restrict(std::string value, double min, double max)
{
  return {value, min, max};
};

struct ParticleID_t
{
  std::string name;
  std::vector<std::string> cuts;
  std::vector<Restriction_t> restrictions;
};

struct PointSelectionResult_t
{
  std::string name;
  std::string root_name;
  std::string tex_name;
  long Ntt=0; //number of tau tau events;
  long Ngg=0; //number of gamma gamma events for luminocity
  double W=0;   //point energy
  double dW=0.01;  
  double L=0;
  double dL=0.01;
  double eps=1.0; //registration efficiencey
  double eps_error=0; //registration efficiencey error
  double effcor = 1.0; // efficiency correction
  double effcor_error = 0; //error of the efficiency correction
};

//struct SelectionResult_t
//{
//  std::vector<point_t> data;
//};

struct ScanPoint_t; 

typedef std::vector<ScanPoint_t>  Scan_t;

std::map<std::string, std::string > SelMap;

struct ScanPoint_t
{
  std::string title;
  int begin_run;
  int end_run;
  double W;
  double dW;
  TTree * tt;
  TTree * gg;
  double L=0;
  std::list<std::pair<int,double> > runs;
  int Ntt;
  int Ngg;
  std::string selection;
  std::map<std::string, int> NttMap;
  double ppi;
  double ppi_max;
  double ppi_min;
  double eps=1.0; //registration efficiency
  double eps_error=1.0; //registration efficiency error
  double effcor=1.0; //registration efficiency correction
  double Sw, dSw; //energy spread
  std::list<int> run_list;
  std::list<std::string> file_list;
  std::list<std::string> regexprs; //regular expessions to match files
  std::string scan_title;

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
};


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
  pr.add("#point",     "%-7d");
  pr.add("name",       "%8s");
  pr.add("W,GeV",      "%10.6f");
  pr.add("dW,GeV",     "%10.6f");
  pr.add("E-Mtau,MeV", "%15.3f");
  pr.add("Sw,MeV",     "%10.3f");
  pr.add("dSw,MeV",    "%10.3f");
  pr.add("L,pb^-1",    "%10.3f");
  pr("  run list", "  %-10s");
  std::cout << pr.head() << std::endl;
  for (const auto & p : SPL) 
    std::cout <<  pr 
      % point_number++ 
      % p.title.c_str() 
      % p.W % p.dW 
      % ((0.5*p.W-MTAU)/MeV) 
      % p.Sw 
      % p.dSw 
      % p.L 
      % get_run_formula(p.run_list).c_str() << std::endl;
}


struct PointDiscriminatorByRunList
{
  const Scan_t * scan;
  PointDiscriminatorByRunList(const Scan_t * scan) : scan(scan) {}
  std::string operator()(std::string file)
  {
    std::regex re(R"(.*\D+(\d+)\.root)");
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



//create run list from single line, where runs in format [run1] [run2] [run3-run4] [run5]
std::list<int> get_run_list(std::string line)
{
  std::list<int> run_list;
  if(line.empty()) return run_list;
  std::string str = line;
  bool is_match=true;
  static std::regex run_range_re(R"(^\s*(\d+)\s*-\s*(\d+))"); //look for the run range
  static std::regex run_re(R"(^\s*\d+)"); //look for the single run
  std::smatch match;
  while(is_match)
  {
    is_match = false;
    if(std::regex_search(str,match,run_range_re)) 
    {
      is_match = true;
      //match[0] -- <begin_run>-<end_run>
      //match[1] -- <begin_run>
      //match[2] -- <end_run>
      //std::cout << match[1] << '-'<<match[2] << std::endl;
      int begin_run = stod(match[1]);
      int end_run = stod(match[2]);
      for(int run = begin_run; run<=end_run; ++run) run_list.push_back(run);
      str = match.suffix();
    }
    else if( std::regex_search(str,match,run_re))
    {
      is_match = true;
      run_list.push_back(stod(match[0]));
      str = match.suffix();
    }
  }
  return run_list;
}

Scan_t read_privalov_runtable(std::string filename)
{
  Scan_t theScan;
  std::fstream ifs(filename);
  if(!ifs) { std::cerr << "Unable to open file: " << filename << std::endl; return theScan; }
  std::string point_name;
  std::string point_energy; //beam energy in GeV
  std::string point_runs;
  //std::regex re(R"((\d+-\d+)");
  //std::regex run_range_re(R"(^\s*(\d+)\s*-\s*(\d+))"); //look for the run range
  //std::regex run_re(R"(^\s*\d+)"); //look for the single run
  std::regex energy_re(R"(\d+\.?\d+)"); //for extract energy
  std::smatch match;
  while( ifs >> point_name >> point_energy) 
  { 
    std::getline(ifs,point_runs); 
    ScanPoint_t sp;
    sp.run_list = get_run_list(point_runs);
    if(sp.run_list.empty()) continue;
    sp.title = point_name;
    if(std::regex_search(point_energy, match, energy_re))
    {
      sp.W = stod(match[0])*2.0;
    }
    else
    {
      sp.W = 0;
    }
    theScan.emplace_back(std::move(sp));
  };
  return theScan;
}

std::vector<std::string> parse_line(std::string regex)
{
  std::vector<std::string> result;
  return result;
}

///* eat comment line */
//struct cline
//{
//  std::string comments="#";
//  cline(std::string s) : comments(s) {}
//};
//
///* eat comment line */
//struct eat
//{
//  std::string eat_list=" \t";
//  cline(std::string s) : eat_list(s) {}
//};
//
//std::istream & oprator>>(std::istream & is, eat e)
//{
//  do
//  {
//    bool is_eat = e.eat_list
//
//    for(int i=0;i<e.eat_list.size();++i) es_eat &
//
//
//  } while(is_eat)
//  return is;
//}

Scan_t read_my_runtable(std::string filename)
{
  Scan_t theScan;
  std::fstream ifs(filename);
  if(!ifs) { std::cerr << "Unable to open file: " << filename << std::endl; return theScan; }
  std::string point_name;
  std::string point_runs;
  double point_energy; //beam energy in GeV
  double point_energy_error;
  double point_spread;
  double point_spread_error;
  double point_lum;
  std::regex comment_re(R"(^\s*#.*)");
  std::string line;
  std::smatch sm;
  while( std::getline(ifs,line)) 
  { 
    if ( line.find_first_of('#') != std::string::npos ) continue;
    //std::cout << line << std::endl;
    std::istringstream iss(line);
    iss >> point_name >> point_energy >> point_energy_error >> point_spread >> point_spread_error >> point_lum;
    std::getline(iss,point_runs); 
    ScanPoint_t sp;
    sp.run_list = get_run_list(point_runs);
    if(sp.run_list.empty()) 
    {
      std::istringstream iss2(point_runs);
      std::string run_combine_regexpr;
      while(iss2>>run_combine_regexpr) sp.regexprs.push_back(run_combine_regexpr);
    }
    else for(auto & r : sp.run_list) sp.regexprs.push_back(std::to_string(r));
    sp.title = point_name;
    sp.W = point_energy;
    sp.dW = point_energy_error;
    sp.Sw = point_spread;
    sp.dSw = point_spread_error;
    sp.L = point_lum;
    theScan.emplace_back(std::move(sp));
  };
  return theScan;
}

Scan_t read_simple_runtable(std::string filename)
{
  Scan_t theScan;
  std::fstream ifs(filename);
  if(!ifs) { std::cerr << "Unable to open file: " << filename << std::endl; return theScan; }
  std::string point_name;
  std::string point_runs;
  double point_energy; //beam energy in GeV
  double point_energy_error;
  double point_spread;
  double point_spread_error;
  double point_lum;
  std::smatch sm;
  std::regex sharp_re(R"(^\s*#.*)"); //for extract energy
  std::regex files_re(R"(\w+\.root)");
  std::string line;
  while(std::getline(ifs,line)) 
  { 
    if(regex_match(line,sm,sharp_re)) continue;
    std::stringstream iss(line);
    //std::cout << point_name << std::endl;
    iss >> point_name >> point_energy >> point_energy_error >> point_spread >> point_spread_error >> point_lum;
    std::getline(iss,point_runs); 
    ScanPoint_t sp;
    for (std::sregex_iterator it(point_runs.begin(), point_runs.end(), files_re); it != std::sregex_iterator(); ++it)
    {
      sp.file_list.push_back(it->str());
    }
    sp.begin_run=-1;
    sp.end_run=-1;
    //std::cout << point_name << " " << point_energy << " " << std::endl;
    //sp.run_list = get_run_list(point_runs);
    //if(sp.run_list.empty()) continue;
    sp.title = point_name;
    sp.W = point_energy;
    sp.dW = point_energy_error;
    sp.Sw = point_spread;
    sp.dSw = point_spread_error;
    sp.L = point_lum;
    theScan.emplace_back(std::move(sp));
  };
  return theScan;
}

const char * make_alias(int channel, const char * templ )
{
  return templ;
}

//substitute in string s by substring according to regexpr
std::string sub(std::string s, std::string regexpr, std::string substr)
{
  std::string result;
  std::regex re(regexpr);
  std::regex_replace(std::back_inserter(result), s.begin(), s.end(), re, substr);
  return result;
};


void set_alias(TTree * tt, double W)
{
  tt->SetAlias("MM","(p[0]+p[1])**2 - (px[0]+px[1])**2 - (py[0]+py[1])**2 - (pz[0]+pz[1])**2");
  tt->SetAlias("Epi0","sqrt(p[0]**2 + 0.1396**2)");
  tt->SetAlias("Epi1","sqrt(p[1]**2 + 0.1396**2)");
  tt->SetAlias("Emu0","sqrt(p[0]**2 + 0.1056**2)");
  tt->SetAlias("Emu1","sqrt(p[1]**2 + 0.1056**2)");
  tt->SetAlias("EK0","sqrt(p[0]**2 + 0.4937**2)");
  tt->SetAlias("EK1","sqrt(p[1]**2 + 0.4937**2)");
  tt->SetAlias("psum","sqrt((px[0]+px[1])**2 + (py[0]+py[1])**2 + (pz[0]+pz[1])**2)");
  tt->SetAlias("M2mumu","(Emu0+Emu1)**2-psum*psum");
  tt->SetAlias("M2pipi","(Epi0+Epi1)**2-psum*psum");
  tt->SetAlias("M2emu","(p[0]+Emu1)**2-psum*psum");
  tt->SetAlias("M2mue","(p[1]+Emu0)**2-psum*psum");
  tt->SetAlias("M2KK","(EK0+EK1)**2-psum*psum");

  char Eb[1024];
  sprintf(Eb,"%5.3f*1",W*0.5);
  tt->SetAlias("Eb",Eb);
  tt->SetAlias("Emis","(2*Eb-p[0]-p[1])");
  tt->SetAlias("cos_theta_mis","(pz[0]+pz[1])/Emis");
  tt->SetAlias("cos_theta_mis2","(pz[0]+pz[1])/hypot(hypot(px[0]+px[1], py[0]+py[1]), pz[0]+pz[1])");
  tt->SetAlias("acol2","(px[0]*px[1]+py[0]*py[1]+pz[0]*pz[1])/(p[0]*p[1])");
  tt->SetAlias("MM2","Emis**2 - (px[0]+px[1])**2 - (py[0]+py[1])**2 - (pz[0]+pz[1])**2");

  auto set_alias = [](TTree * tt, std::string track_name_prefix, int track, std::string selection_template) -> void
  {
    std::string track_name = track_name_prefix + std::to_string(track);
    tt->SetAlias(track_name.c_str(), sub(selection_template,R"(#)",std::to_string(track)).c_str());
  };

  for(int i=0;i<4;i++)
  {
    set_alias(tt,"Mrho_cut", i, "0.5 < Mrho[#] && Mrho[#] < 1.0");
  }

  /*
  //Now setup particle id's
  for(int track=0;track<2;++track)
  {
    set_alias(tt, "e",   track, "0.8 < Ep[#] && Ep[#]<1.05 && chi2_dedx_e[#]<2 && abs(delta_tof_e[#])<0.2");
    set_alias(tt, "u",   track, "!e# && 0.1 < E[#] && E[#] < 0.3 && depth[#]>0 && chi2_dedx_mu[#] < 5  && abs(delta_tof_mu[#]) < 0.3 && Ep[#]<0.8");
    set_alias(tt, "pi",  track, "!(e#||u#) && 0.8 < p[#] && p[#] < 1.05 && Ep[#]<0.6 && chi2_dedx_pi[#] < 5 && abs(delta_tof_pi[#])<0.3");
    set_alias(tt, "K",   track, "!(e#||u#||pi#) && 0.8 < p[#] && p[#] < 1.05 && Ep[#]<0.6 && chi2_dedx_K[#] < 3 && abs(delta_tof_K[#])<0.3");
    set_alias(tt, "rho", track, "pi# && Mrho_cut# && Mpi0[#]>0.1128 && Mpi0[#]<0.1464");
  }
  */

  //define channels
  auto make_xy  = [](TTree * tt,  std::string x, std::string y)
  {
    auto s1 = sub("((x0 && y1) || (y0  && x1))", "x", x);
    auto alias = sub(s1, "y", y);
    tt->SetAlias( (x+y).c_str(), alias.c_str());
  };
  for (auto & x : {"e","u", "pi", "K"} )
    for (auto & y : {"e","u", "pi", "K"} )
      make_xy(tt, x,y);

  //auto make_xx  = [](TTree * tt,  std::string x, std::string extra_cut = "")
  //{
  //  auto alias = sub("x0 && x1", "x", x);
  //  tt->SetAlias((x+x).c_str(), (alias + (extra_cut == "" ? "" : "&& (" + extra_cut + ")")).c_str() );
  //};



  tt->SetAlias("pil", "(pi0 && (u1 || e1 ))");
  tt->SetAlias("lpi", "(pi1 && (u0 || e0 ))");

  //tt->SetAlias("Xrho",   "!rho0 && (e0 || u0 || pi0 || K0) && rho1 && !(e1 || u1 || K0)");
  //tt->SetAlias("rhoX",   "!rho1 && (e1 || u1 || pi1 || K1) && rho0 && !(e0 || u0 || K0)");
  tt->SetAlias("Xrho",   "!rho0 && rho1");
  tt->SetAlias("rhoX",   "!rho1 && rho0");
  tt->SetAlias("rhorho", "!(e0 || u0 || K0) &&  !(e1 || u1 || K1) && Mpi0[0]>0.1128 && Mpi0[0]<0.1464 && Mpi0[1]>0.1128 && Mpi0[1]<0.1464  && pi0 && pi1  && ((Mrho_cut1 && Mrho_cut2) || (Mrho_cut1 && Mrho_cut3))");
  tt->SetAlias("lpipi0_cut","(Nc==2 && Npi0==1 && ( (pil && Mrho[0]>0.5 && Mrho[0]<1.0 && ) || ((lpi && (Mrho[1]>0.5 && Mrho[1]<1.0))) ) && acop>1.4)");
  //tt->SetAlias("pipi0pipi0_cut","Nc==2 && Npi0==2 && pipi && ptem>0.6");
  tt->SetAlias("all_cut",   "lpipi0_cut || eu_cut || epi_cut");
}

Scan_t read_data(std::string data_dir, const Scan_t & cfg, std::string filter=R"(\.root$)")
{
  std::cout << "Reading data directory \"" << data_dir << "\":\n";
  auto fl  = filter_file_list(get_recursive_file_list(data_dir)); //get file list *.root
  //for(auto & f : fl) std::cout << f << std::endl;
  auto ml  = combine(fl,PointDiscriminatorByRunList(&cfg)); //combine files into points
  //for(auto & m : ml) 
  //{
  //  std::cout << m.first << ": ";
  //  for(auto & f: m.second) std::cout << f<<" "; 
  //  std::cout << std::endl;
  //}
  Scan_t scan;
  scan.reserve(cfg.size());
  for(auto & point : cfg ) 
  {
    ScanPoint_t sp;
    sp = point;
    sp.file_list = ml[point.title];
    if(sp.file_list.empty()) 
    {
      std::cout << "WARNING: no data for point: " << point.title << std::endl;
      continue;
    }
    scan.push_back(sp);
  }
  //read TChain
  int index = 0;
  for(auto & sp : scan)
  {
    //create chain
    //sp.tt= new TChain((, sp.title.c_str());
    //read tau tau events
    auto tt = new TChain("tt", sp.title.c_str());
    for(auto & file : sp.file_list) tt->AddFile(file.c_str());
    //read gamma gamma events
    auto gg = new TChain("gg", sp.title.c_str());
    for(auto & file : sp.file_list) gg->AddFile(file.c_str());
    //change chain names
    tt->SetName(("tt"+std::to_string(index)).c_str());
    gg->SetName(("gg"+std::to_string(index)).c_str());
    set_alias(tt,sp.W);
    sp.tt = tt;
    sp.gg = gg;
  }
  print(scan);
  return scan;
}

std::vector<ScanPoint_t> read_data2(std::string data_dir, std::string privalov_runtable)
{
  auto cfg = read_privalov_runtable(privalov_runtable); //read scan configuration
  auto fl  = filter_file_list(get_recursive_file_list(data_dir)); //get file list *.root
  auto ml  = combine(fl,PointDiscriminatorByRunList(&cfg)); //combine files into points
  Scan_t scan;
  scan.reserve(cfg.size());
  for(auto & point : cfg ) 
  {
    ScanPoint_t sp;
    sp = point;
    sp.file_list = ml[point.title];
    if(sp.file_list.empty()) 
    {
      std::cout << "WARNING: no data for point: " << point.title << std::endl;
      continue;
    }
    scan.push_back(sp);
  }
  //read TChain
  int index = 0;
  for(auto & sp : scan)
  {
    //create chain
    //sp.tt= new TChain((, sp.title.c_str());
    //read tau tau events
    auto tt = new TChain("tt", sp.title.c_str());
    for(auto & file : sp.file_list) tt->AddFile(file.c_str());
    //read gamma gamma events
    auto gg = new TChain("tt", sp.title.c_str());
    for(auto & file : sp.file_list) gg->AddFile(file.c_str());
    //change chain names
    tt->SetName(("tt"+std::to_string(index)).c_str());
    gg->SetName(("gg"+std::to_string(index)).c_str());
    set_alias(tt,sp.W);
    sp.tt = tt;
    sp.gg = gg;
  }
  return scan;
}

std::vector<ScanPoint_t> read_data3(std::string data_dir, std::string cfg_file)
{
  Scan_t scan = read_simple_runtable(cfg_file); //read scan configuration
  print(scan);
  int index = 0;
  for(auto & sp : scan)
  {
    //create chain
    //sp.tt= new TChain((, sp.title.c_str());
    //read tau tau events
    auto tt = new TChain("tt", sp.title.c_str());
    for(auto & file : sp.file_list) tt->AddFile((data_dir+"/"+file).c_str());
    //read gamma gamma events
    auto gg = new TChain("tt", sp.title.c_str());
    for(auto & file : sp.file_list) gg->AddFile((data_dir+"/"+file).c_str());
    //change chain names
    tt->SetName(("tt"+std::to_string(index)).c_str());
    gg->SetName(("gg"+std::to_string(index)).c_str());
    set_alias(tt,sp.W);
    sp.tt = tt;
    sp.gg = gg;
  }
  return scan;
}

#include <regex>
std::vector<ScanPoint_t> read_mc(std::string  dirname=".", std::string regexpr=R"(.+\.root)")
{
  std::vector<ScanPoint_t> P;
  TSystemDirectory dir(dirname.c_str(), dirname.c_str());
  std::regex file_re(regexpr); //regular expression to filter files
  std::regex energy_re(R"(^\D*(\d+\.?\d+).root$)");
  std::smatch file_match;
  std::smatch energy_match; 
  int point=0;
  for(auto  file: * dir.GetListOfFiles())
  {
    std::string file_name(file->GetName());
    if(std::regex_match (file_name,file_match,file_re))
    {
      //extracting the energy
      if(std::regex_match(file_name,energy_match,energy_re))
      {
        double W = std::stod(energy_match[1]);
        //std::cout << file_name <<  "  W = " << W << "  " << W-MTAU << std::endl;
        if(W*0.5-MTAU > -0.010 &&  W*0.5-MTAU < 0.030) 
        {
          P.push_back({"P"+std::to_string(point), 55116, 55155,  W, 1e-5,0,0});
          P.back().tt = get_chain("tt",("tt"+std::to_string(point)).c_str(), "MC GALUGA", (dirname+"/"+file_name).c_str());
          P.back().gg = get_chain("gg",("gg"+std::to_string(point)).c_str(), "gg lum", (dirname+"/"+file_name).c_str());
          set_alias(P.back().tt, P.back().W);
        }
      }
      else
      {
        std::cout << "Unable to extract energy from filename: " << file_name << std::endl;
      }
    }
  }
  //std::cout << "THE END" << std::endl;
  return P;
}


void draw_lum_per_run2(const std::vector<ScanPoint_t> & SPL)
{
  TMultiGraph * mg = new TMultiGraph;
  int i=0;
  std::vector<TGraphErrors*> g(SPL.size());
  int point=0;
  for (const auto & p : SPL)
  {
    g[point] = new TGraphErrors;
    g[point]->SetMarkerStyle(21);
    g[point]->SetMarkerColor(point+1);
    for(auto & ri : p.runs)
    {
      auto run = ri.first;
      auto lum = ri.second;
      auto chain = get_chain("gg", "gg","g", run, run);
      double Ngg = chain->GetEntries();
      double dNgg = sqrt(Ngg);
      int n = g[point]->GetN();
      g[point]->SetPoint(n,run, Ngg/lum);
      g[point]->SetPointError(n, 0, dNgg/lum);
      i++;
    }
    mg->Add(g[point],"p");
    point++;
  }
  mg->Draw("a");
  mg->GetXaxis()->SetTitle("run");
  mg->GetYaxis()->SetTitle("N_{e^{+}e^{-1} #rightarrow #gamma #gamma} / L_{online}, nb");
}

void  draw_lum_per_run(const char * dir, int begin_run, int end_run)
{
  auto * g = new TGraphErrors;
  for ( int run = begin_run; run<=end_run; ++run)
  {
    auto chain = get_chain("gg", "gg","g", run, run);
    double Ngg = chain->GetEntries();
    double dNgg = sqrt(Ngg);
    int n = g->GetN();
    g->SetPoint(n,run, Ngg);
    g->SetPointError(n, 0, dNgg);
  }
  g->Draw("a*");
}

void draw_tau_per_run(const std::vector<ScanPoint_t> & SPL, const char * selection = "")
{
  TMultiGraph * mg = new TMultiGraph;
  int i=0;
  std::vector<TGraphErrors*> g(SPL.size());
  int point=0;
  for (const auto & p : SPL)
  {
    g[point] = new TGraphErrors;
    g[point]->SetMarkerStyle(21);
    g[point]->SetMarkerColor(point+1);
    for(auto & ri : p.runs)
    {
      auto run = ri.first;
      auto lum = ri.second;
      auto chain = get_chain("gg", "gg","g", run, run);
      auto tt = get_chain("tt","tt","t",run,run);
      double Ngg = chain->GetEntries();
      double dNgg = sqrt(Ngg);
      double L = Ngg/10*pow(p.W/MTAU*0.5,2);
      double dL = dNgg/10*pow(p.W/MTAU*0.5,2);
      double Ntt = tt->GetEntries(selection);
      double dNtt = sqrt(Ntt);
      double sigma = Ntt/L;
      double dsigma = Ntt/L* sqrt(1./Ntt + 1.0/Ngg); 
      int n = g[point]->GetN();
      g[point]->SetPoint(n,run, sigma);
      g[point]->SetPointError(n, 0, dsigma);
      i++;
    }
    mg->Add(g[point],"p");
    point++;
  }
  mg->Draw("a");
  mg->GetXaxis()->SetTitle("run");
  mg->GetYaxis()->SetTitle("N_{#tau #tau} / L_{online}, nb");
}

TGraphErrors * draw_result(const char * selection, const std::vector<ScanPoint_t> & Points)
{
  TGraphErrors * g = new TGraphErrors;
  long totalNtt=0;
  for(int i=0; i<Points.size();++i)
  {
    double Ntt = Points[i].tt->GetEntries(selection);
    double Ngg = Points[i].gg->GetEntries(selection);
    double xs  = 0;
    double dxs = 0;
    totalNtt += Ntt;
    if(Ngg != 0 )
    {
      xs = Ntt/Ngg;
      dxs = sqrt( Ntt/(Ngg*Ngg) + pow(Ntt/(Ngg*Ngg), 2.0)*Ngg );
    }
    g->SetPoint(i, Points[i].W/2.0-MTAU, xs);
    g->SetPointError(i, Points[i].dW/2.0, dxs);
    std::cout << i << " " << Points[i].W/2.0-MTAU << "  " << Ngg << "  " << Ntt << "   " <<  xs << std::endl;
  }
  std::cout << "Total number of tau-tau candidates:" << totalNtt << std::endl;
  g->SetMarkerStyle(21);
  g->Draw("ap");

  TCanvas * cacop = new TCanvas("acop","acop");
  cacop->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cacop->cd(i+1);
    Points[i].tt->Draw("acop",selection);
  }
  TCanvas * cptem = new TCanvas("ptem","ptem");
  cptem->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cptem->cd(i+1);
    Points[i].tt->Draw("ptem",selection);
  }
  TCanvas * cp = new TCanvas("p","p");
  cp->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cp->cd(i+1);
    Points[i].tt->Draw("p",selection);
  }
  return g;
}



void add_last(std::vector<ScanPoint_t> & P)
{
  for ( auto & p : P)
  {
    p.NttMap[p.selection]=p.Ntt;
  }
}
void set_last(std::vector<ScanPoint_t> & P)
{
  for ( auto & p : P)
  {
    p.NttMap.clear();
    p.NttMap[p.selection]=p.Ntt;
  }
}

void fit(std::vector<ScanPoint_t> & P, const char * filename="scan.txt", bool nofit=false, std::string title="")
{
  long totalNtt=0;
  long totalNgg = 0;
  double totalL=0;
  std::string total_title="";
  for(int i=0; i<P.size();++i)
  {
    for(auto & item: P[i].NttMap)
    {
      totalNtt += item.second;
    }
    if(P[i].Ngg==0) P[i].Ngg = P[i].gg->GetEntries();
    totalNgg+=P[i].Ngg;
    totalL += P[i].L;
  }
  //total_title = total_title + " " + P[i].selection;
  double sigma_gg;
  if( totalNgg == 0 || totalL == 0) sigma_gg = 1;
  else sigma_gg = totalNgg == 0 ? 1 : totalNgg/totalL;
  std::ofstream ofs(filename);
  ofs << "#Selection: " << P[0].selection << std::endl;
  ofs << "#Aliases: " << std::endl;
  ofs << "#    e0 = " << P[0].tt->GetAlias("e0") << std::endl;
  ofs << "#    e1 = " << P[0].tt->GetAlias("e1") << std::endl;
  ofs << "#    u0 = " << P[0].tt->GetAlias("u0") << std::endl;
  ofs << "#    u1 = " << P[0].tt->GetAlias("u1") << std::endl;
  ofs << "#    pi0 = " << P[0].tt->GetAlias("pi0") << std::endl;
  ofs << "#    pi1 = " << P[0].tt->GetAlias("pi1") << std::endl;
  ofs << "#" << std::setw(4) << " " << std::setw(15) << "L, nb-1" << std::setw(10) << "dL, nb-1" << std::setw(15) << "W, MeV" << std::setw(15) << "dW, MeV" << std::setw(10) << "SW, MeV" << std::setw(10) << "dSW, MeV";
  ofs << std::setw(10) << "Ntt" << std::setw(10) << "Nee" << std::setw(10) << "Ngg" << std::setw(10) << "effcor" << std::endl;
  for(int i=0; i<P.size();++i)
  {
    int Ntt=0;
    for(auto & item: P[i].NttMap) Ntt+=item.second;
    int Ngg = P[i].Ngg;
    double L = Ngg==0 ? 1 : Ngg/(sigma_gg*pow(P[i].W/(2*MTAU),2.0));
    ofs << std::setw(5) << i <<  std::setw(15) << L*1000 << std::setw(10) << 10 << std::setw(15) << P[i].W/MeV  << std::setw(15) << P[i].dW/MeV;
    ofs << std::setw(10) << 1.256 << " " << std::setw(10) << 0.019;
    ofs << std::setw(10) << Ntt << std::setw(10) <<  1 << std::setw(10) << Ngg <<  std::setw(10) << P[i].effcor << std::endl;
  }
  if(!nofit)
  {
    char command[65536];
    //if(title=="") title=total_title;
    sprintf(command, "taufit --tau-spread=1.256 --title='sigma: %s' '%s' --output '%s.txt' &", title.c_str(), filename,filename);
    //sprintf(command, "taufit --tau-spread=1.256  '%s' --output '%s.txt' &",  filename,filename);
    system(command);
  }
}

void fit_last(std::vector<ScanPoint_t> & P, const char * filename="scan.txt", bool nofit=false, std::string title="")
{
  set_last(P);
  fit(P,filename,nofit,title);
}



void mccmp(const std::vector<ScanPoint_t> & P, const char * selection, const char * cut, const char * his="")
{
  std::string title = std::string("mc_data: ") + selection + " cut= " + cut;
  TCanvas * c = new TCanvas("data_mc_diff", title.c_str());
  TChain * chain = new TChain("data_mc_diff","data_mc_diff");
  set_alias((TTree*)chain, MTAU*2);
  for(int i=1;i<P.size();++i) chain->Add((TChain*)P[i].tt);
  std::string sel(selection);
  sel+=">>h1(";
  sel+=his;
  sel+=")";
  chain->Draw(sel.c_str(),cut,"goff");
  auto hdata  = chain->GetHistogram();
  P[P.size()-1].tt->SetLineColor(kRed);
  sel=selection;
  sel+=">>h2(";
  sel+=his;
  sel+=")";
  P[P.size()-1].tt->Draw(sel.c_str(), cut, "goff");
  auto hmc = P[P.size()-1].tt->GetHistogram();
  hmc->Draw();
  hmc->GetXaxis()->SetTitle(selection);
  hmc->SetLineWidth(2);
  hdata->Scale(double(hmc->GetEntries())/hdata->Integral());
  hdata->SetLineWidth(2);
  hdata->Draw("SAME");
}

TGraphErrors * draw_result2(const std::vector<ScanPoint_t> & Points, const char * selection)
{
  TGraphErrors * g = new TGraphErrors;
  long totalNtt=0;
  std::ofstream ofs("scan.txt");
  for(int i=0; i<Points.size();++i)
  {
    double Ntt = Points[i].tt->GetEntries(selection);
    double Ngg = Points[i].gg->GetEntries();
    double xs  = 0;
    double dxs = 0;
    totalNtt += Ntt;
    if(Ngg != 0 )
    {
      xs = Ntt/Ngg;
      dxs = sqrt( Ntt/(Ngg*Ngg) + pow(Ntt/(Ngg*Ngg), 2.0)*Ngg );
    }
    g->SetPoint(i, Points[i].W/2.0-MTAU, xs);
    g->SetPointError(i, Points[i].dW/2.0, dxs);
    std::cout << i << " " << Points[i].W/2.0-MTAU << "  " << Ngg << "  " << Ntt << "   " <<  xs << std::endl;
    auto & P = Points[i];
    ofs << std::setw(5) << i <<  std::setw(15) << 20000 << "  " << 10 << std::setw(15) << P.W/MeV  << std::setw(15) << P.dW/MeV;
    ofs << std::setw(10) << 1.256 << " " << std::setw(10) << 0.019;
    ofs << std::setw(10) << Ntt << std::setw(10) << " " << 1 << "  " << Ngg <<  " " << 1 << std::endl;
  }
  std::cout << "Total number of tau-tau candidates:" << totalNtt << std::endl;
  TCanvas * c = new TCanvas("cross","cross");
  c->cd();
  g->SetMarkerStyle(21);
  g->Draw("ap");

  std::vector<const char *> var = { "acop", "ptem", "p" ,"acol"};

  TCanvas * can = new TCanvas("param","param");
  can->Divide(var.size(),Points.size());

  for( int i=0;i<Points.size();++i) 
  { 
    for(int v=0;v<var.size();v++)
    {
      can->cd(1+i*var.size()+v);
      Points[i].tt->Draw(var[v],selection);
    }
  }
  return g;
}


struct ChannelSelection_t
{
  std::string title; //channel name  (temrinal)
  std::string root_title; //title for root to print on canvas
  std::string cut;
};

struct ChannelSelectionResult_t : public ChannelSelection_t
{
  long Ntt=0; //total number of events for all points
  long Ngg=0; //total number of gg events
  std::vector<PointSelectionResult_t> Points;
  ChannelSelectionResult_t(void){}
  ChannelSelectionResult_t(const ChannelSelection_t & cs, const std::vector<PointSelectionResult_t> & p) : ChannelSelection_t(cs)
  {
    Points = p;
    for(auto & p : Points)
    {
      Ntt+=p.Ntt;
      Ngg+=p.Ngg;
    } 
  }
  ChannelSelection_t & operator=(const std::vector<PointSelectionResult_t> & P)
  {
    Ntt=0;
    Ngg=0;
    Points = P;
    for(auto & p : Points)
    {
      Ntt+=p.Ntt;
      Ngg+=p.Ngg;
    } 
    return *this;
  }
  ChannelSelection_t & operator=(const ChannelSelection_t & cs )
  {
    title      = cs.title;
    root_title = cs.root_title;
    cut        = cs.cut;
    return *this;
  }
};

struct Selection
{
  std::string name;
  std::string common_cut;
  std::vector<ParticleID_t> pid;
  std::vector<ChannelSelection_t> sel;
  ChannelSelection_t       & operator[](int i)         { return sel[i]; }
  const ChannelSelection_t & operator[](int i) const   { return sel[i]; }
};



std::vector<ChannelSelection_t> SELECTION
{
  {"eμ", "e#mu"     , "Nc==2 && Nn==0 && eu   && ptem>0.25"}   ,
  {"eπ", "e#pi"     , "Nc==2 && Nn==0 && epi  && ptem>0.3"}    ,
  {"μπ", "#mu#pi"   , "Nc==2 && Nn==0 && upi  && ptem>0.62"}   ,
  {"ee", "ee"       , "Nc==2 && Nn==0 && ee   && ptem>0.5  && acop<2.7 && p[0]<0.9 && p[1]<0.9"}     ,
//{"ee", "ee"       , "ee && Nn==0 && Nc==2 && ptem>0.4 &&  S>0.05 && abs(cos_theta_mis2)<0.6"},
  {"μμ", "#mu#mu"   , "Nc==2 && Nn==0 && uu   && ptem>0.5 && acop<2.7 && p[0]<0.9 && p[1]<0.9"}     ,
  //{"μμ", "#mu#mu"   , "Nc==2 && Nn==0 && uu   && ptem>0.4 && abs(cos_theta_mis2) <0.7"}     ,
  {"ππ", "#pi#pi"   , "Nc==2 && Nn==0 && pipi && ptem>0.55"}    ,
  {"eK", "eK"       , "Nc==2 && Nn==0 && eK"}                   ,
  {"πK", "#piK"     , "Nc==2 && Nn==0 && piK  && ptem >0.4"}    ,
  {"μK", "#uK"      , "Nc==2 && Nn==0 && uK   && ptem >0.15"}   ,
  {"KK", "KK"       , "Nc==2 && Nn==0 && KK   && ptem >0.15"}   ,
  {"Xρ", "X#rho"    , "Nc==2 && Npi0 ==1 && (Xrho || rhoX) && acop>1.4"} ,
  {"ρρ", "#rho#rho" , "Nc==2 && Npi0 ==2 && rhorho && ptem>0.5"}
};

std::vector<ChannelSelection_t> SELECTION11
{
  {"eμ", "e#mu"     , "eu   && Nc==2 && Nn==0 && ptem>0.1 && acop>16./180.*3.1415 && acop<(180-16.)/180.*3.1415"}     ,
  {"eπ", "e#pi"     , "epi  && Nc==2 && Nn==0 && ptem>0.22"}    ,
  {"μπ", "#mu#pi"   , "upi  && Nc==2 && Nn==0 && ptem>0.25 && p[0]<0.94 && p[1]<0.94"}    ,
  {"ee", "ee"       , "ee   && Nc==2 && Nn==0 && ptem>0.33"}     ,
  {"μμ", "#mu#mu"   , "uu   && Nc==2 && Nn==0 && ptem>0.35"}     ,
  {"ππ", "#pi#pi"   , "pipi && Nc==2 && Nn==0 && ptem>0.33 && M2 >0.8"}   ,
  {"KK", "KK"       , "KK   && Nc==2 && Nn==0 && ptem>0.33 && M2 >0.8"}   ,
  {"eK", "eK"       , "eK   && Nc==2 && Nn==0 && ptem>0.22"}     ,
  {"πK", "#piK"     , "piK  && Nc==2 && Nn==0 && ptem>0.33 && M2 > 0.8"}    
};

std::vector<ChannelSelection_t> SELECTION2
{
  {"eμ", "e#mu"     , "Nc==2 && Nn==0 && eu   && ptem>0.25" }  ,
  {"eπ", "e#pi"     , "nciptrack ==2 && Nc==2 && Nn==0 && epi  && ptem>0.19"}   ,
  {"μπ", "#mu#pi"   , "nciptrack ==2 && Nc==2 && Nn==0 && upi  && ptem>0.22"}   ,
  {"ee", "ee"       , "nntrack==0 && nciptrack ==2 && Nc==2 && Nn==0 && ee   && ptem>0.3"}  ,
  {"μμ", "#mu#mu"   , "nciptrack ==2 && Nc==2 && Nn==0 && uu   && ptem>0.25"}     ,
  {"ππ", "#pi#pi"   , "nciptrack ==2 && Nc==2 && Nn==0 && pipi && ptem>0.2"}    ,
  {"eK", "eK"       , "nciptrack ==2 && Nc==2 && Nn==0 && eK"}                   ,
  {"πK", "#piK"     , "nntrack==0 && nciptrack ==2 && Nc==2 && Nn==0 && piK  && ptem >1"}    ,
  {"μK", "#uK"      , "nciptrack ==2 && Nc==2 && Nn==0 && uK   && ptem >0.2"}   ,
  {"KK", "KK"       , "nciptrack ==2 && Nc==2 && Nn==0 && KK   && ptem >0.16"}   
};

std::string SEL3_GLOBAL = "Nc==2 && Nn==0 && p[0] < 1.1 && p[1] < 1.1  && abs(cos_theta_mis2)<0.8 && cos(theta[0])<0.8 && cos(theta[1])<0.8";
std::vector<ChannelSelection_t> SELECTION3
{
  {"eμ", "e#mu"     , "eu   && ptem > 0.3"},
  {"eπ", "e#pi"     , "epi  && ptem > 0.3"},
  {"μπ", "#mu#pi"   , "upi  && ptem > 0.3"},
  {"ee", "ee"       , "ee   && ptem > 0.4"},
  {"μμ", "#mu#mu"   , "uu   && ptem > 0.3"},
  {"ππ", "#pi#pi"   , "pipi && ptem > 0.3"},
  {"eK", "eK"       , "eK   && ptem > 0.0"},
  {"πK", "#piK"     , "piK  && ptem > 0.3"},
  {"μK", "#uK"      , "uK   && ptem > 0.3"},
  {"KK", "KK"       , "KK   && ptem > 0.1"}   
};


std::vector<ChannelSelection_t> & DEFAULT_SELECTION=SELECTION;




std::string to_string(time_t t, int TZ)
{
  std::string s;
  return s;
};


void print_event_info(ScanPoint_t  & p, const char * selection, const char * title = "")
{
  long n = p.tt->Draw("run:event:time",selection,"goff");
  setenv("TZ", "Asia/Shanghai",1);
  auto runs   = p.tt->GetV1();
  auto events = p.tt->GetV2();
  auto times  = p.tt->GetV3();
  std::ofstream file("tau_events.txt",std::ios_base::app);
  std::cout << '#' << std::setw(7) << "run" << std::setw(20)<< "eventId" << std::setw(10) << "point" << std::setw(10) << "channel" << std::setw(10) << "pnt_evt" << std::setw(30) << " time" << std::endl;;
  file << '#' << std::setw(7) << "run" << std::setw(20)<< "eventId" << std::setw(10) << "point" << std::setw(10) << "channel" << std::setw(10) << "pnt_evt" << std::setw(30) << " time" << std::endl;;
  for(int i=0;i<n;++i)
  {
    time_t t = time_t(times[i]);
    std::cout << std::setw(8) << runs[i] << std::setw(20)<< long(events[i]) << std::setw(10) << p.title << std::setw(10) << title << std::setw(10) << i<< std::setw(30) << ctime(&t);
    file <<  std::setw(8) << runs[i] << std::setw(20)<< long(events[i]) << std::setw(10) << p.title << std::setw(10) << title << std::setw(10) << i<< std::setw(30) << ctime(&t);
  }
  unsetenv("TZ");
}

void print_event_info(std::vector<ScanPoint_t> & P, const char * selection, const char * title)
{
  for(auto & p: P )
  {
    print_event_info(p,selection,title);
  }
};

void print_event_info(std::vector<ScanPoint_t> & P, std::vector<ChannelSelection_t> & SEL)
{
  for( auto & s : SEL)
  {
    print_event_info(P,s.cut.c_str(), s.title.c_str());
  }
};


void set_pid(std::vector<ScanPoint_t> & DATA, const std::vector<ParticleID_t> & PID)
{
  for( auto & data : DATA)
  {
    for(int track=0;track<2;++track)
    {
      for(auto & pid: PID)
      {
        std::string alias_value;
        std::string name = pid.name + std::to_string(track);
        for (auto  cut: pid.cuts)
        {
          cut = sub(cut,R"(#)",std::to_string(track));
          //std::cout << name << " pid_cut = " <<  cut << std::endl; 
          if(alias_value=="") alias_value = cut; 
          else alias_value += " && " + cut;
        }
        for (auto & r : pid.restrictions)
        {
            std::string value_name = r.name;

            std::string pid_name= pid.name;
            if(pid_name == "u") pid_name = "mu";

            if ( value_name == "chi2_dedx"  || value_name == "delta_tof")  value_name += "_"+pid_name;
            value_name += "["+std::to_string(track)+"]";
            std::string cut = std::to_string(r.min) + " < " + value_name + " && " + value_name  + " < "  +  std::to_string(r.max);
            //std::cout << cut << std::endl;
            if(alias_value=="") alias_value = cut; 
            else alias_value += " && " + cut;
        };
//        std::cout << " alias: " << name << " = " << alias_value << std::endl;
        data.tt->SetAlias(name.c_str(), alias_value.c_str());
      }
    }
  };
};



std::vector<PointSelectionResult_t> new_select(const std::vector<ScanPoint_t> & P, const std::string & sel)
{
  std::vector<PointSelectionResult_t> R(P.size());
  for(int i=0;i<P.size();++i)
  {
    auto & r       = R[i];
    auto & p       = P[i];
    r.Ntt          = p.tt->GetEntries(sel.c_str());
    if(p.gg) r.Ngg = p.gg->GetEntries();
    r.name         = p.title;
    r.root_name    = p.title;
    r.tex_name     = p.title;
    r.W            = p.W;
    r.dW           = p.dW;
    r.L            = p.L;
  }
  return R;
}

std::vector<PointSelectionResult_t> draw(const std::vector<ScanPoint_t> & P, const std::string & sel, const std::string & var, std::string gopt="")
{
  auto c = new TCanvas;
  int canvas_pos_x = 0;
  int canvas_width_x = 1920*2/5.0*P.size();
  int canvas_width_y = 1080*0.5;
  static int canvas_pos_y = 0;
  c->SetWindowPosition(canvas_pos_x,canvas_pos_y);
  canvas_pos_y+=canvas_width_y+60;
  canvas_pos_y%=(1080*2);
  c->SetWindowSize(canvas_width_x,canvas_width_y);
  c->SetTitle(sel.c_str());
  c->Divide(P.size(),1);
  std::vector<PointSelectionResult_t> R(P.size());
  for(int i=0;i<P.size();++i)
  {
    c->cd(i+1);
    auto & r       = R[i];
    auto & p       = P[i];
    r.Ntt          = p.tt->Draw(var.c_str(),sel.c_str(),gopt.c_str());  
    std::string result_title = p.title  + ": " + std::to_string(r.Ntt) + " events";
    p.tt->GetHistogram()->SetTitle(result_title.c_str());
    if(p.gg) r.Ngg = p.gg->GetEntries();
    r.name         = p.title;
    r.root_name    = p.title;
    r.tex_name     = p.title;
    r.W            = p.W;
    r.dW           = p.dW;
    r.L            = p.L;
  }
  return R;
}

std::vector<PointSelectionResult_t> draw(std::vector<ScanPoint_t> & DATA, const Selection & S, int i, const std::string & var,  std::string extracut="", std::string gopt="")
{
  set_pid(DATA, S.pid);
  std::string cut = S[i].cut;
  if(S.common_cut!="") cut += " && " + S.common_cut;
  if(extracut!="")     cut += " && " + extracut;
  return draw(DATA, cut, var, gopt);
};


template<typename Iterator, typename Func>
Iterator find_best(Iterator it_begin, Iterator it_end, Func F)
{
  using ItemType = typename Iterator::value_type;
  using Type = typename std::result_of<Func(ItemType)>::type;
  Type chi2  = std::numeric_limits<Type>::max();
  Iterator result=it_end;
  for(auto it = it_begin; it!=it_end; ++it)
  {
    auto x = F(*it);
    if( x < chi2 ) 
    {
      result = it;
      chi2 = x;
    }
  }
  return result;
};

template<typename Container, typename Func>
typename Container::iterator find_best(Container & c, Func F)
{
  return find_best(std::begin(c), std::end(c), F);
}

ChannelSelectionResult_t  fold(const std::vector<ChannelSelectionResult_t>  & SR, std::string name = "all")
{
  ChannelSelectionResult_t result;
  result.title = name;
  if(SR.empty()) return result;
  result.Points.resize(SR[0].Points.size());
  result.Ntt=0;
  for(int i=0;i<result.Points.size();++i) 
  {
    auto & rp   = result.Points[i];
    rp.eps = 0;
    double eps_error2=0;
    for(const auto & r : SR )  
    {
      auto & p    = r.Points[i];
      result.Ntt += p.Ntt;
      rp.eps     += p.eps;
      eps_error2 += p.eps_error*p.eps_error;
      rp.Ntt     += p.Ntt;
      rp.Ngg      = p.Ngg;
      rp.L        = p.L;
      rp.dL       = p.dL;
      rp.W        = p.W;
      rp.dW       = p.dW;
    }
    rp.eps_error = sqrt(eps_error2);
  }
  auto & reference_point = * find_best(result.Points, [](const auto & p) { return  std::abs(p.W*0.5-MTAU); } );
  //normalize efficiency correction to threshold efficiency
  for(int i=0;i<result.Points.size();++i)
  {
    auto & rp   = result.Points[i];
    rp.effcor = rp.eps / reference_point.eps;
    rp.effcor_error = rp.eps_error / reference_point.eps;
  }
  reference_point.effcor_error = 0;
  return result;
};

void print(const ChannelSelectionResult_t & sr, int opt=1 , int first_column_width=10, int last_column_width=5, int  column_width= 8, int vline_width=6)
{
  int N = sr.Points.size(); //number of points
  int hline_width = first_column_width + last_column_width + vline_width + N*column_width+3;
  auto hline = [&hline_width](std::string  symb="─", int width = 0) 
  { 
    int w = width == 0 ? hline_width : width;
    for(int i=0;i<w;++i)
    {
      std::cout << symb;
    }
    std::cout << std::endl;
  };
  auto vline = [&vline_width](char symb='-', int width = 0) { std::cout << std::setw(width == 0? (vline_width-1) : (width-1)) << symb; };

  if(opt<0) hline(); //preparing for total
  switch (opt)
  {
    case 0:
      hline("━");
      std::cout << std::setw(first_column_width) << "CHNNL/PNT";
      for(int i=0;i<N;++i) std::cout << std::setw(column_width)  << sr.Points[i].name;
      std::cout << std::setw(vline_width) << " │ " << std::setw(last_column_width) << "TOTAL";
      std::cout << std::endl;
      hline();
      std::cout << std::setw(first_column_width+count_utf8_extra_byte(std::string("ΔE,MeV"))) << "ΔE,MeV";
      for(int i=0;i<N;++i) std::cout << std::setw(column_width) << (sr.Points[i].W*0.5-MTAU)*1e3;
      std::cout << std::setw(vline_width) << " │ ";
      std::cout << std::endl;
      hline();
    default:
      std::cout << std::setw(first_column_width+count_utf8_extra_byte(sr.title)) << sr.title;
      for(int i=0; i < N; i++) std::cout << std::setw(column_width) <<  sr.Points[i].Ntt;
      std::cout << std::setw(vline_width) << " │ " << std::setw(last_column_width) <<  sr.Ntt;
      std::cout << std::endl;
  }
  if(opt<0) hline("━");
};


struct PrintConfig_t
{
  int title_width=10;
  int data_width=15;
  int total_width=5;
  int vline_width=6;
  std::string thline=""; //top horizontal line
  std::string bhline=""; //bottom horizontal line
  PrintConfig_t() = default;
  void hline(int N, std::string symb="─")
  {
    int hline_width = title_width + total_width + vline_width + N*data_width+3;
    //int w = width == 0 ? hline_width : width;
    for(int i=0;i<hline_width;++i)
    {
      std::cout << symb;
    }
    std::cout << std::endl;
  }
};

PrintConfig_t PCFG;

template < typename Item, typename Ftitle , typename Fdata, typename Ftotal>
void print_smth(const std::vector<Item> pts, PrintConfig_t cfg, Ftitle ftitle, Fdata fdata, Ftotal ftotal)
{
  int N = pts.size(); //number of points
  int hline_width = cfg.title_width + cfg.total_width + cfg.vline_width + N*cfg.data_width+3;
  auto hline = [&hline_width](std::string  symb="─", int width = 0) 
  { 
    int w = width == 0 ? hline_width : width;
    for(int i=0;i<w;++i)
    {
      std::cout << symb;
    }
    std::cout << std::endl;
  };
  auto vline_width = cfg.vline_width;
  auto vline = [&vline_width](char symb='-', int width = 0) { std::cout << std::setw(width == 0? (vline_width-1) : (width-1)) << symb; };
  print_utf(cfg.title_width,ftitle());
  for(int i=0;i<N;++i) print_utf(cfg.data_width, fdata(pts[i]));
  std::cout << std::setw(vline_width) << " │ ";
  print_utf(cfg.total_width, ftotal(pts));
  std::cout << std::endl;
};

void print_head(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg)
{
  print_smth(pts, cfg, [](){ return "CHNL/PNT"; }, [](auto & p) { return p.name; }, [](auto & s) { return "TOTAL"; } );
};

void print_eps(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg=PrintConfig_t())
{
  print_smth(pts, cfg, [](){ return "ε"; }, [](auto & p) { return p.eps; }, [](auto & s) { return ""; } );
};

void print_Ntt(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg, std::string title="Ntt")
{
  print_smth(pts, cfg, 
      [&title](){ return title; }, 
      [](auto & p) { return p.Ntt; }, 
      [](auto & s) { 
        long N=0; 
        for(auto & x: s)  N+=x.Ntt;  
        return N; 
        } 
      );
};

template<typename Fdata, typename Ftotal>
void print_smth(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg, Fdata fd, Ftotal ft)
{
  if(SR.empty()) return;
  auto hline = [&SR, &cfg](std::string s="─") { cfg.hline(SR[0].Points.size(),s); };
  hline("━");
  print_head(SR[0].Points,cfg);
  hline();
  auto prn = [&](std::string title, auto & sr) {
    print_smth(sr.Points, cfg, 
        [&title](){ return title; }, 
        fd, 
        ft
        );
  };
  for( auto &  sr : SR) 
     prn(sr.title,sr);
  hline();
  auto f = fold(SR);
  prn("all", f);
  hline("━");
};

void print_Ntt(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  print_smth(SR,cfg, 
      [](auto & p) { 
        return p.Ntt; 
      },
      [](auto & s) { 
        double sum=0;
        for(auto & x: s)  sum+=x.Ntt;  
        return sum; 
      } 
  );
};

void print_efficiency(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  if(SR.empty()) return;
  std::string format_str = "%4.2f ± %4.2f";
  cfg.total_width = format_str.length();
  auto hline = [&SR, &cfg](std::string s="─") { cfg.hline(SR[0].Points.size(),s); };
  hline("━");
  print_smth(SR[0].Points, cfg, [](){ return "CHNL/PNT"; }, [](auto & p) { return p.name; }, [](auto & s) { return "AVERAGE"; } );
  hline();
  auto prn = [&](std::string title, auto & sr) {
    print_smth(sr.Points, cfg, 
          [&title](){ return title; }, 
          [&format_str](auto & p) { 
            char buf[1024];
            sprintf(buf,format_str.c_str(), p.eps*100,p.eps_error*100);
            return std::string(buf); 
          },
          [&format_str](auto & P) { 
            double sum=0;
            double error2_sum=0;
            for(auto & p: P)  
            {
              sum+=p.eps;  
              error2_sum+=p.eps_error*p.eps_error;
            }
            double average = sum/P.size(); 
            double error = sqrt(error2_sum)/P.size();
            char buf[1024];
            sprintf(buf,format_str.c_str(), average*100,error*100);
            return std::string(buf); 
          } 
        );
  };
  for( auto &  sr : SR) 
     prn(sr.title,sr);
  hline();
  auto f = fold(SR);
  prn("all", f);
  hline("━");
};

void print_effcor(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  if(SR.empty()) return;
  std::string format_str = "%4.2f ± %4.2f";
  cfg.total_width = format_str.length();
  auto hline = [&SR, &cfg](std::string s="─") { cfg.hline(SR[0].Points.size(),s); };
  hline("━");
  print_smth(SR[0].Points, cfg, [](){ return "CHNL/PNT"; }, [](auto & p) { return p.name; }, [](auto & s) { return "AVERAGE"; } );
  hline();
  auto prn = [&](std::string title, auto & sr) {
    print_smth(sr.Points, cfg, 
          [&title](){ return title; }, 
          [&format_str](auto & p) { 
            char buf[1024];
            sprintf(buf,format_str.c_str(), p.effcor,p.effcor_error);
            return std::string(buf); 
          },
          [&format_str](auto & P) { 
            double sum=0;
            double error2_sum=0;
            for(auto & p: P)  
            {
              sum+=p.effcor;  
              error2_sum+=p.effcor_error*p.effcor_error;
            }
            double average = sum/P.size(); 
            double error = sqrt(error2_sum)/P.size();
            char buf[1024];
            sprintf(buf,format_str.c_str(), average,error);
            return std::string(buf); 
          } 
        );
  };
  for( auto &  sr : SR) 
     prn(sr.title,sr);
  hline();
  auto f = fold(SR);
  prn("all", f);
  hline("━");
};


void print(const std::vector<PointSelectionResult_t> & Points, int opt=1 , int first_column_width=10, int last_column_width=5, int  column_width= 8, int vline_width=6)
{
  ChannelSelectionResult_t sr;
  sr = Points;
  print(sr,0);
}


void print(const  std::vector<ChannelSelectionResult_t> & SR )
{
  PrintConfig_t cfg;
  print_Ntt(SR,cfg);
};

ChannelSelectionResult_t new_select(const std::vector<ScanPoint_t> & P, const ChannelSelection_t & S, std::string extra_cut="")
{
  if(extra_cut!="") extra_cut = " && " + extra_cut;
  ChannelSelectionResult_t r;
  r = S; //copy selection config
  r  = new_select(P, S.cut + extra_cut);
  print(r,0);
  return r;
};

std::vector<ChannelSelectionResult_t> new_select( std::vector<ScanPoint_t> & P, const std::vector<ChannelSelection_t> & S, std::string extra_cut="")
{
  std::vector<ChannelSelectionResult_t> R(S.size());
  if(extra_cut!="") extra_cut = " && " + extra_cut;
  for(int i=0; i<S.size(); ++i)
  {
    auto & s = S[i];
    auto & r = R[i];
    r = s; //save current selection for channel
    r.cut = s.cut + extra_cut;
    r = new_select(P, s.cut + extra_cut);
    print(r,i);
  };
  print(fold(R,"all"),-1);
  return R;
};

std::vector<ChannelSelectionResult_t> new_select(std::vector<ScanPoint_t> & DATA, Selection & SEL, std::string extracut="")
{
  set_pid(DATA, SEL.pid);
  if(extracut!="") extracut = " && " + extracut;
  return new_select(DATA, SEL.sel, SEL.common_cut + extracut);
};

ChannelSelectionResult_t new_select(std::vector<ScanPoint_t> & DATA, const Selection & SEL, int i, std::string extracut="")
{
  set_pid(DATA, SEL.pid);
  if(extracut!="") extracut = " && " + extracut;
  return new_select(DATA, SEL.sel[i], SEL.common_cut + extracut);
};

std::vector<ChannelSelectionResult_t> new_select(std::vector<ScanPoint_t> & DATA, const Selection & SEL, std::vector<int> idx, std::string extracut="")
{
  std::vector<ChannelSelectionResult_t> R(idx.size());
  set_pid(DATA, SEL.pid);
  std::string cut = SEL.common_cut;
  if(extracut != "")
  {
    if(SEL.common_cut !="") cut = SEL.common_cut +"&&"+extracut;
    else cut = extracut;
  }
  if(cut!="") cut = "&&"+cut;

  int j=0;
  for(auto i : idx)
  {
    auto & s = SEL.sel[i];
    auto & r = R[i];
    r = s;
    r.cut = SEL.sel[i].cut +  cut;
    r = new_select(DATA, SEL.sel[i].cut +  cut);
    print(r,j++);
  }
  print(fold(R,"all"),-1);
  return R;
};


std::vector<PointSelectionResult_t> measure_efficiency(const std::vector<ScanPoint_t> & P, std::string cut)
{
  return new_select(P,cut);
}

std::vector<PointSelectionResult_t> & set_efficiency(std::vector<PointSelectionResult_t>  & psr,  const std::vector<PointSelectionResult_t>  & eff, long N0)
{
  int i0=0;
  const double UNSET = std::numeric_limits<double>::max();
  double chi2_thr = UNSET;
  for(int i = 0; i< psr.size(); ++i) {
    auto & r = psr[i];
    int jc; //correspoinding point
    double chi2 = UNSET;
    for( int j = 0;j<eff.size(); ++j) {
      if (  std::abs(eff[j].W-r.W) < chi2 ) {
        jc = j;
        chi2 = std::abs(eff[j].W-r.W); 
      }
    }
    assert ( chi2 == UNSET);
    double eps = double(eff[jc].Ntt)/N0;
    double eps_error =  sqrt( eps/N0 * ( 1.0 - eps));
    //std::cout << i << " " << r.name << " " << "  W=" << r.W << "  eps = " << eps << " +- " << eps_error << std::endl;
    if( std::abs(r.W*0.5 - MTAU) <  chi2_thr )  {
      chi2_thr = std::abs(r.W*0.5 - MTAU);
      i0 = i;
    }
    r.eps = eps;
    r.eps_error = eps_error;

  }
  auto & reference_point = * find_best(psr, [](const auto & p) { return  std::abs(p.W*0.5-MTAU); } );
  //normalize efficiency correction to threshold efficiency
  for(auto & rp : psr) {
    rp.effcor = rp.eps / reference_point.eps;
    rp.effcor_error = rp.eps_error / reference_point.eps;
  }
  reference_point.effcor_error = 0;
  return psr;
}

std::vector<ChannelSelectionResult_t> & set_efficiency( std::vector<ChannelSelectionResult_t> & SR, const Scan_t & MC, long N0)
{
  for(auto & sr : SR)
  {
    auto eff = measure_efficiency(MC,sr.cut);
    set_efficiency(sr.Points, eff, N0);
  };
  print_efficiency(SR);
  return SR;
};


void draw_hist_vs_energies(Scan_t & scan, const char * var, int Nbin, double xmin, double xmax, const char *selection)
{
  std::vector<int> color={kBlack, kRed, kBlue,kGreen+2};
  std::vector<int> line = {1,2,3};
  TLegend * l = new TLegend(0.8,0.8,1.0,1.0);
  for(int i=0;i<scan.size();++i)
  {
    std::string s = var;
    int c = i % color.size();
    int l = i / color.size();
    std::cout << i << " " << l << " " << c << std::endl;
    scan[i].tt->SetLineColor(color[c]);
    scan[i].tt->SetLineStyle(line[l]);
    scan[i].tt->SetLineWidth(3);
    s+=" >> h"+std::to_string(i)+"("+std::to_string(Nbin)+","+std::to_string(xmin)+","+std::to_string(xmax)+")";
    if(i==0) scan[i].tt->Draw(s.c_str(), selection);
    else scan[i].tt->Draw(s.c_str(),selection,"same");
  }
}






void compare(const ScanPoint_t & p1, const ScanPoint_t &p2, const char * varexp, const char * selection, const char * gopt)
{
  auto  c = new TCanvas;
  c->SetWindowSize(1200, 1800);
  c->Divide(1,2);
  c->cd(1);
  p1.tt->Draw(varexp,selection,gopt);
  c->cd(2);
  p2.tt->Draw(varexp,selection,gopt);
}


void correct_efficiency(Scan_t & D /*  data  */, Scan_t & MC /*  MC signal */)
{
  for(auto & P : D)
  {
    //if( (P.W*0.5-MTAU)/MeV < 1 ) P.effcor = 1.0;
    //double effcor = 1.0;
    double energy_distance = 1e10;
    for ( auto  & mc : MC)
    {
      if ( fabs(P.W-mc.W) < energy_distance )
      {
        std::cout << "mc.effcor = " << mc.effcor << " " << P.W << " " << mc.W << std::endl;
        P.effcor = mc.effcor;
        energy_distance = fabs(P.W-mc.W);
      }
    }
  }

  std::cout << "#Correct efficiency: " << std::endl;
  std::cout << "#" << std::setw(4) << " " << std::setw(10) << "W,MeV" << std::setw(10) << "E, MeV" << std::setw(15) << "E-Mtau, MeV" << std::setw(10) << "Ntt" << std::setw(10) << "effcor" << std::endl;
  for( int i = 0;i<D.size(); i++)
  {
    std::cout << std::setw(5) << i+1 << std::setw(10) << D[i].W/MeV << std::setw(10) << D[i].W*0.5/MeV << std::setw(15) << (D[i].W*0.5-MTAU)/MeV <<  std::setw(10) << D[i].Ntt << std::setw(10) << D[i].effcor << std::endl;
  }
}

void set_color_and_line(TAttLine *obj, int index)
{
  std::vector<int> color={kBlack, kRed, kBlue,kGreen+2};
  std::vector<int> line = {1,2,3};
  int c = index % color.size();
  int l = index / color.size();
  obj->SetLineColor(color[c]);
  obj->SetLineStyle(line[l]);
}

void compare(std::vector<ScanPoint_t*> SP, std::vector<std::string> varexp, const char * selection, std::string gopt)
{
  auto  c = new TCanvas;
  int win_size=600;
  c->SetWindowSize(win_size*(SP.size()+1), win_size*varexp.size());
  int index =0;
  c->Divide(SP.size()+1,varexp.size());
  for(auto & ve : varexp)
  {
    TLegend * l = new TLegend(0.8,0.9,1.0,1.0);
    auto hs = new THStack(ve.c_str(),ve.c_str());
    int draw_index=0;
    for(auto sp :  SP)
    {
      std::string  gopt_tmp = gopt;
      sp->tt->SetLineWidth(1);
      set_color_and_line(sp->tt, draw_index);
      c->cd(index+1);
      sp->tt->Draw(ve.c_str(),selection, gopt.c_str());
      auto h = sp->tt->GetHistogram();
      h->SetName(("h"+sp->title+std::to_string(index)).c_str());
      h->SetTitle((sp->title).c_str());
      hs->Add(h);
      index++;
      l->AddEntry(sp->tt,sp->title.c_str());
      draw_index++;
    }
    c->cd(index+1);
    hs->Draw();
    l->Draw();
    index++;
  }
}

void compare(std::vector<ScanPoint_t*> SP, const char * varexp_str, const char * selection, std::string gopt)
{
  compare(SP,std::vector<std::string>{varexp_str}, selection, gopt);
}

 std::string test_replace(std::string particle_name_prefix, std::string selection_template, int track)
{
  std::string track_str = std::to_string(track);
  std::string particle_name = particle_name_prefix + track_str;
  std::regex re(R"(\[\*\])");
  std::string selection;
  std::regex_replace(std::back_inserter(selection), selection_template.begin(), selection_template.end(), re,"["+track_str+"]");
  return selection;
};


void save(const ChannelSelectionResult_t & sr, std::string  filename="scan.txt")
{
  std::cout << "Saving selection: " << sr.title << " to file: " << filename << std::endl;
  std::stringstream os;
  os << "#" << std::setw(4) << " " << std::setw(15) << "L, nb-1" << std::setw(10) << "dL, nb-1" << std::setw(15) << "W, MeV" << std::setw(15) << "dW, MeV" << std::setw(10) << "SW, MeV" << std::setw(10) << "dSW, MeV";
  os << std::setw(10) << "Ntt" << std::setw(10) << "Nee" << std::setw(10) << "Ngg" << std::setw(10) << "effcor" << "\n";
  int idx=0;
  for(const auto & p : sr.Points)
  {
    idx++;
    if(p.name=="") os << std::setw(5) << idx;
    else           os << std::setw(5) << p.name;
    os << std::setw(15) << p.L*1000 << std::setw(10) << 10 << std::setw(15) << p.W/MeV  << std::setw(15) << p.dW/MeV;
    os << std::setw(10) << 1.256 << " " << std::setw(10) << 0.019;
    os << std::setw(10) << p.Ntt << std::setw(10) <<  1 << std::setw(10) << p.Ngg <<  std::setw(10) << p.effcor << "\n";
  }
  std::cout << os.str();
  std::ofstream ofs(filename);
  ofs << os.str();
}

void fit(const ChannelSelectionResult_t & sr, std::string  filename="scan.txt", std::string title="")
{
  save(sr,filename);
  char command[65536];
  sprintf(command, "taufit --tau-spread=1.256 --title='sigma: %s' '%s' --output '%s.txt' &", title.c_str(), filename.c_str(),filename.c_str());
  system(command);
}

void fit(const std::vector<PointSelectionResult_t>  & ps, std::string  filename="scan.txt", std::string title="")
{
  ChannelSelectionResult_t sr;
  sr = ps;
  save(sr,filename);
  char command[65536];
  sprintf(command, "taufit --tau-spread=1.256 --title='sigma: %s' '%s' --output '%s.txt' &", title.c_str(), filename.c_str(),filename.c_str());
  system(command);
}


#include <TH1F.h>
void find_optimal_depth(ScanPoint_t & sp, double kmin = 0, double kmax=100, double dk=0.1)
{
  TGraph * g =new TGraph;
  for(double k = kmin; k<kmax; k+=dk)
  {
    char buf[1024];
    sprintf(buf,"depth[0]-p[0]*%f >> h", k);
    sp.tt->Draw(buf,"depth[0]>0 && depth[0]<100 && abs(pid[0])==13 && p[0]<1.1 && depth[0]-56.76*p[0] >-40 && depth[0]-56.76*p[0]<-10");
    TH1 * h = sp.tt->GetHistogram();
    double rms = h->GetRMS();
    double sigma = h->GetStdDev();
    //std::cout << k << " " << rms << "  " << sigma << std::endl;
    g->SetPoint(g->GetN(), k, sigma);
    delete h;
  }
  g->Draw("a*");
  double minimum = TMath::MinElement(g->GetN(), g->GetY());
  std::cout << "optimal k = " << minimum << std::endl;
  std::cout << "g->GetMinimum() " << g->GetMinimum() << std::endl;
}


void show_depth(ScanPoint_t & sp, double k=58.6, double k0=-40)
{
  auto c = new TCanvas;
  c->Divide(1,2);
  TTree * t = sp.tt;
  c->cd(1);
  t->SetMarkerStyle(21);
  t->SetMarkerColor(kBlue);
  t->Draw("depth[0]:p[0]","depth[0]>0 && depth[0]<100 && abs(pid[0]) == 13");
  t->SetMarkerStyle(22);
  t->SetMarkerColor(kRed);
  t->Draw("depth[0]:p[0]","depth[0]>0 && depth[0]<100 && abs(pid[0]) == 211", "same");
  auto l  = new TLine(0.6, k0 + 0.6*k,  1.05, k0+1.05*k);
  l->SetLineWidth(3);
  l->Draw();

  //std::string mu = "depth[0]>0 && depth[0]<100 && abs(pid[0]) == 13 &&   depth[0]-p[0]*" + std::to_string(k) + " >= " + std::to_string(k0);
  //std::string pi = "depth[0]>0 && depth[0]<100 && abs(pid[0]) == 211 &&  depth[0]-p[0]*" + std::to_string(k) + " < " + std::to_string(k0);
  std::string mu = "depth[0]>0 && depth[0]<60 && abs(pid[0]) == 13";
  std::string pi = "depth[0]>0 && depth[0]<60 && abs(pid[0]) == 211";
  std::string varexp = "depth[0] - p[0]*" + std::to_string(k);
  c->cd(2);
  t->SetMarkerStyle(21);
  t->SetMarkerColor(kBlue);
  t->SetLineColor(kBlue);
  t->Draw(varexp.c_str(),mu.c_str());
  t->SetMarkerStyle(22);
  t->SetMarkerColor(kRed);
  t->SetLineColor(kRed);
  t->Draw(varexp.c_str(),pi.c_str(),"same");

}

void show_ep(ScanPoint_t & sp, std::string var = "Ep[0]", std::string cut = "", std::string opt= "")
{
  TTree * t = sp.tt;
  auto draw = [&t](std::string var, std::string sel, std::string opt, int color, int style, int size)
  {
    t->SetLineColor(color); 
    t->SetLineWidth(size);
    t->SetMarkerStyle(style);
    t->SetMarkerColor(color);
    t->SetMarkerSize(size);
    t->Draw(var.c_str(),sel.c_str(),opt.c_str());
  };
  if( cut!="" ) cut+=" && ";
  draw(var, cut + "abs(pid[0]) ==  11", opt +     "", kBlue,  21, 1);
  draw(var, cut + "abs(pid[0]) ==  13", opt + "SAME", kRed,   21, 1);
  draw(var, cut + "abs(pid[0]) == 211", opt + "SAME", kGreen-3, 21, 1);
};

void find_Ep(ScanPoint_t &sp, int Nbin=0, double parmin=0, double parmax=0)
{
  std::string var="Ep[0]";
  std::string varexpr = var;
  std::string cut = "abs(pid[0])==11";
  if(Nbin>0) varexpr += " >> h("+std::to_string(Nbin)+")";
  if(parmin!=0) cut += " && " + var + " > " + std::to_string(parmin);
  if(parmax!=0) cut += " && " + var + " < " + std::to_string(parmax);
  sp.tt->Draw(var.c_str(),cut.c_str());
  TH1 * h = sp.tt->GetHistogram();
  double ymax = h->GetMaximum();
  int nmax1 = h->FindFirstBinAbove(ymax-1);
  int nmax2 = h->FindLastBinAbove(ymax-1);
  int nmax = (nmax1+nmax2)*0.5;
  double xmin = h->FindFirstBinAbove(ymax/2);
  double xmax = h->FindLastBinAbove(ymax/2);
  std::cout << "maximum at " << nmax1 << " " << nmax2<< "  half width: " << xmin << " " << xmax << std::endl;
  std::cout << "maximum at " << h->GetBinCenter(nmax1) << " " << h->GetBinCenter(nmax2) << "  half width: " << h->GetBinCenter(xmin) << " " << h->GetBinCenter(xmax) << std::endl;
  double I0 = h->Integral();
  std::cout << "Integreal: " << I0 << std::endl;
  TGraph * g  =new TGraph;
  std::map<std::string, std::pair<double, double> > Mrange;
  for(int bin =0; bin<100; bin++)
  {
    double I = h->Integral(nmax-bin, nmax+bin);
    double F = I/I0;
    double xmin = h->GetBinCenter(nmax-bin);
    double xmax = h->GetBinCenter(nmax+bin);
    if( F<0.5) Mrange["0.5"] = {xmin,xmax}; 
    if( F<0.86) Mrange["0.86"] = {xmin,xmax}; 
    if( F<0.9) Mrange["0.9"] = {xmin,xmax}; 
    if( F<0.95) Mrange["0.95"] = {xmin,xmax}; 
    if( F<0.99) Mrange["0.99"] = {xmin,xmax}; 
    if( F<0.999) Mrange["0.999"] = {xmin,xmax}; 
    
    double x = h->GetBinCenter(nmax-bin);
    g->SetPoint(g->GetN(), x, I/I0);
  }
  new TCanvas;
  g->Draw("a*l");
  for( auto & m :  Mrange )
  {
    auto level = m.first;
    auto range = m.second;
    auto xmin = range.first;
    auto xmax = range.second;
    std::cout << "level  " << level << " : " <<  xmin << " " << xmax << std::endl;
  }
};
