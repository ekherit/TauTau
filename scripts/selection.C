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
#include <numeric>
#include <functional>

#include <variant>

#include <typeinfo>

#include "time.h"
#include "stdlib.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPaveStats.h>
#include <TStyle.h>

#include "Selection.h"
#include "draw_helper.h"
#include "fold.h"

#include "../ibn/valer.h"
#include "../ibn/indexer.h"
#include "utils.h"

#include "pdg_table.h"

constexpr double GeV=1.0;
constexpr double MeV=1e-3*GeV;
constexpr double barn = 1.0;
constexpr double pb = 1e-12*barn;
constexpr double nb = 1e-9*barn;

constexpr double MTAU=1776.86*MeV;
constexpr double MPI=0.13957061*GeV; 
constexpr double MPI0=0.134977*GeV;


constexpr double ME_PDG2011=0.510998910*MeV; //+-0.000000013 
constexpr double ME_PDG2019=0.5109989461*MeV; //+-0.000000013
constexpr double ME_PDG=ME_PDG2019;// mass of electron
constexpr double ME=ME_PDG;// mass of electron

constexpr double SIGMA_TOMSON = 0.665245854*barn; 
constexpr double SIGMA_CONST = SIGMA_TOMSON/pb* ME*ME/2.;  // GeV^2/pb


constexpr double ALPHA_PDG2019=1/137.035999139; 
constexpr double ALPHA = ALPHA_PDG2019;
constexpr double ALPHAPI=ALPHA/TMath::Pi(); //alpha / pi
constexpr double PIALPHA=ALPHA*TMath::Pi(); 


//std::string TAUFIT = "taufit --lum=bes --tau-spread=1.258 --energy-correction=-0.078";
static std::string TAUFIT = "taufit --lum=bes --tau-spread=1.258 --energy-correction=+0.011";


void clean_taufit(void) {
  system("killall -9 taufit");
}

static std::string BB_SEL = "(acol-TMath::Pi())>-0.03 && abs(cos(theta[0]) < 0.8 && abs(cos(theta[1])) < 0.8";
static std::string GG_SEL = "";






struct AcceleratorInfo_t {
  ibn::valer<double> energy;     //beam c.m. energy
  ibn::valer<double> energy_spread; //beam energy spread
  ibn::valer<double> luminosity; 
};

struct DataSample_t  : public AcceleratorInfo_t
{
  std::string title;
  std::shared_ptr<TTree> tree;
  ibn::valer<double> cross_section;
  //ibn::valer<double> energy;     //beam c.m. energy
  //ibn::valer<double> luminosity;
  ibn::valer<double> efficiency; //registration efficiency
  ibn::valer<double> effcor; //correction to efficiency
  long N; //number of selected events in data
  long Nmc; //number of selected events in MC
  long N0mc; //initial number of MC events
  long Draw(std::string varexp, std::string selection="", std::string option="", long nentries = std::numeric_limits<long>::max(), long firstentry=0) const {
    return tree->Draw(varexp.c_str(), selection.c_str(), option.c_str(), nentries, firstentry);
  }
};

struct ScanPoint_t; 
typedef std::reference_wrapper<ScanPoint_t> ScanPointRef_t;

typedef std::vector<ScanPoint_t>  Scan_t;
typedef std::reference_wrapper<Scan_t> ScanRef_t;
typedef std::shared_ptr<Scan_t> ScanPtr_t;

std::map<std::string, std::string > SelMap;


struct Simulation_t {
  ScanRef_t signal;
  std::vector<ScanRef_t> bgs;
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

using PointSelectionResult_t = ScanPoint_t;

//struct PointSelectionResult_t
//{
//  std::string name;
//  std::string root_name;
//  std::string tex_name;
//  std::string cut;
//  long Ntt=0; //number of tau tau events;
//  ibn::valer<double> W;
//  ibn::valer<double> L;
//  ibn::valer<double> eps = {1.0,0};
//  ibn::valer<double> effcor = {1.0,0};
//  DataSample_t bb, gg;
//};


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
      % p.energy.value % p.energy.error 
      % ((0.5*p.energy.value-MTAU)/MeV) 
      % p.energy_spread.value
      % p.energy_spread.error
      % p.luminosity.value 
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
      sp.energy = stod(match[0])*2.0;
    }
    else {
      sp.energy = 0;
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
  ibn::valer<double>point_energy; //GeV
  ibn::valer<double>point_lum; 
  //double point_energy; //beam energy in GeV
  //double point_energy_error;
  double point_spread;
  double point_spread_error;
  std::regex comment_re(R"(^\s*#.*)");
  std::string line;
  std::smatch sm;
  while( std::getline(ifs,line)) 
  { 
    if ( line.find_first_of('#') != std::string::npos ) continue;
    //std::cout << line << std::endl;
    std::istringstream iss(line);
    iss >> point_name >> point_energy.value >> point_energy.error >> point_spread >> point_spread_error >> point_lum.value >> point_lum.error;
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
    sp.energy = point_energy;
    sp.energy_spread = point_spread;
    //sp.energy_spread.error = point_spread_error;
    sp.luminosity = point_lum;
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
  ibn::valer<double> point_energy; //GeV
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
    iss >> point_name >> point_energy.value >> point_energy.error >> point_spread >> point_spread_error >> point_lum;
    std::getline(iss,point_runs); 
    ScanPoint_t sp;
    for (std::sregex_iterator it(point_runs.begin(), point_runs.end(), files_re); it != std::sregex_iterator(); ++it)
    {
      sp.file_list.push_back(it->str());
    }
    //sp.begin_run=-1;
    //sp.end_run=-1;
    //std::cout << point_name << " " << point_energy << " " << std::endl;
    //sp.run_list = get_run_list(point_runs);
    //if(sp.run_list.empty()) continue;
    sp.title = point_name;
    sp.energy = point_energy;
    //sp.dW = point_energy_error;
    //sp.Sw.value = point_spread;
    //sp.Sw.error = point_spread_error;
    sp.energy_spread = point_spread;
    sp.luminosity = point_lum;
    theScan.emplace_back(std::move(sp));
  };
  return theScan;
}

const char * make_alias(int channel, const char * templ )
{
  return templ;
}




void set_alias(TTree * tt, double W, double L=1.0)
{
  tt->SetAlias("good_emc_time", "(0<=ntemc[0]&&ntemc[0]<=14&&0<=ntemc[1]&&ntemc[1]<=14)");
  tt->SetAlias("cgood","(abs(cos(theta[0]))<0.93 && abs(cos(theta[1]))<0.93)");
  tt->SetAlias("barrel","abs(cos(theta[0]))<0.8 && abs(cos(theta[1]))<0.8");
  tt->SetAlias("missed_photon_angle2","( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )");
  tt->SetAlias("missed_photon_angle", "( abs(cos_theta_mis) < 0.8 || ( 0.92 > abs(cos_theta_mis) && abs(cos_theta_mis) > 0.86) )");
  tt->SetAlias("missed_photon_angle3", "( abs(cos(d4p3[3])) < 0.8 || ( 0.92 > abs(cos(d4p3[3])) && abs(cos(d4p3[3])) > 0.86))");
  tt->SetAlias("missed_photon_angle4", "( 0.92 > abs(cos_theta_mis) )");
  tt->SetAlias("missed_photon_angle5", "( ( 0.93 > abs(cos_theta_mis) && abs(cos_theta_mis) > 0.84) || (0.83 > abs(cos_theta_mis)) )");
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
  //tt->SetAlias("NnE100","Sum$(nE>0.1)");
  //tt->SetAlias("NnE50","Sum$(nE>0.050)");
  //tt->SetAlias("NnE25","Sum$(nE>0.025)");
  //tt->SetAlias("Ngbarrel","Sum$( abs(cos(ntheta))<0.8 && nE>0.025)");
  //tt->SetAlias("Ngendcup","Sum$( abs(cos(ntheta))<0.92 && abs(cos(ntheta)) >0.86  && nE>0.05)");
  //tt->SetAlias("Ng","Ngbarrel+Ngendcup");



  tt->SetAlias("lum",myfmt("%6.2f*1",L).c_str());
  char Eb[1024];
  sprintf(Eb,"%5.3f*1",W*0.5);
  tt->SetAlias("Eb",Eb);
  tt->SetAlias("MPI",(std::to_string(MPI)+"*1").c_str());
  tt->SetAlias("Emis2","(2*Eb-p[0]-p[1])");
  //tt->SetAlias("cos_theta_mis","(pz[0]+pz[1])/Emis");
  tt->SetAlias("cos_theta_mis2","(pz[0]+pz[1])/hypot(hypot(px[0]+px[1], py[0]+py[1]), pz[0]+pz[1])");
  tt->SetAlias("acol2","(px[0]*px[1]+py[0]*py[1]+pz[0]*pz[1])/(p[0]*p[1])");
  tt->SetAlias("MM2","Emis**2 - (px[0]+px[1])**2 - (py[0]+py[1])**2 - (pz[0]+pz[1])**2");
  tt->SetAlias("MM2pi","(2*Eb-sqrt(p[0]**2+MPI**2)-sqrt(p[1]**2+MPI**2))**2 - (px[0]+px[1])**2 - (py[0]+py[1])**2 - (pz[0]+pz[1])**2");
  tt->SetAlias("M2pi","(sqrt(p[0]**2+MPI**2)+sqrt(p[1]**2+MPI**2))**2 - (px[0]+px[1])**2 - (py[0]+py[1])**2 - (pz[0]+pz[1])**2");
  tt->SetAlias("k_ptem","(1.0*1.0)");

  tt->SetAlias("ptem2","sqrt((px[0]+px[1]+npx[0]+npx[1])^2 + (py[0]+py[1]+npy[0]+npy[1])^2)/Emis");

  auto set_alias = [](TTree * tt, std::string track_name_prefix, int track, std::string selection_template) -> void
  {
    std::string track_name = track_name_prefix + std::to_string(track);
    tt->SetAlias(track_name.c_str(), sub(selection_template,R"(#)",std::to_string(track)).c_str());
  };

  for(int i=0;i<4;i++)
  {
    set_alias(tt,"Mrho_cut", i, "0.5 < Mrho[#] && Mrho[#] < 1.0");
  }

  //define channels
  auto make_xy  = [](TTree * tt,  std::string x, std::string y)
  {
    auto s1 = sub("((x0&&y1)||(y0&&x1))", "x", x);
    auto alias = sub(s1, "y", y);
    tt->SetAlias( (x+y).c_str(), alias.c_str());
  };
  for (auto & x : {"e","u", "pi", "K","rho","X"} )
    for (auto & y : {"e","u", "pi", "K","rho","X"} )
      make_xy(tt, x,y);

  tt->SetAlias("pil", "(pi0 && (u1 || e1 ))");
  tt->SetAlias("lpi", "(pi1 && (u0 || e0 ))");

  //tt->SetAlias("eX",   "((!e0 && e1)||(!e1 && e0))");
  tt->SetAlias("Xrho",   "((!rho0 && rho1)||(!rho1 && rho0))");
  tt->SetAlias("lpipi0_cut","(Nc==2 && Npi0==1 && ( (pil && Mrho[0]>0.5 && Mrho[0]<1.0 && ) || ((lpi && (Mrho[1]>0.5 && Mrho[1]<1.0))) ) && acop>1.4)");
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

    auto bb = new TChain("bb", sp.title.c_str());
    for(auto & file : sp.file_list) bb->AddFile(file.c_str());
    //change chain names
    tt->SetName(("tt"+std::to_string(index)).c_str());
    bb->SetName(("bb"+std::to_string(index)).c_str());
    set_alias(tt,sp.energy,sp.luminosity);
    sp.tt.tree.reset(tt);
    sp.gg.tree.reset(gg);
    sp.bb.tree.reset(bb);
    sp.tt.energy = sp.energy;
    sp.gg.energy = sp.energy;
    sp.bb.energy = sp.energy;
    sp.type = data_dir;
  }
  print(scan);
  return scan;
}

/*
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
    set_alias(tt,sp.W,sp.L);
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
    set_alias(tt,sp.W,sp.L);
    sp.tt = tt;
    sp.gg = gg;
  }
  return scan;
}
*/

std::vector<ScanPoint_t> read_mc(std::string  dirname=".", Scan_t cfg={}, long N0mc=1e6, std::string regexpr=R"(.+\.root)")
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
    if(std::regex_match (file_name,file_match,file_re)) {
      //extracting the energy
      if(std::regex_match(file_name,energy_match,energy_re)) {
        double W = std::stod(energy_match[1]);
        if(cfg.empty()){
        }
        else {
          //find closest energy
          auto & sp  = *std::min_element(cfg.begin(), cfg.end(), [W](auto a, auto b){ return fabs(a.energy-W)<fabs(b.energy-W); } );
          if( abs(sp.energy - W) < 0.001 ) {
            P.push_back(sp);
            auto & p = P.back();
            P.back().type = dirname;
            p.title = dirname+p.title;
            p.title = sub(p.title,R"(mc/)","");
            p.title = sub(p.title,R"(T)","");
            p.energy = W;
            p.luminosity = -p.luminosity;
            p.tt.tree.reset(get_chain("tt",("tt"+p.title).c_str(), "signal", (dirname+"/"+file_name).c_str()));
            p.gg.tree.reset(get_chain("gg",("gg"+p.title).c_str(), "gg lum", (dirname+"/"+file_name).c_str()));
            p.bb.tree.reset(get_chain("bb",("bb"+p.title).c_str(), "bb lum", (dirname+"/"+file_name).c_str()));
            p.tt.energy = W;
            p.gg.energy = W;
            p.bb.energy = W;
            p.tt.N0mc = N0mc;
            p.gg.N0mc = N0mc;
            p.bb.N0mc = N0mc;
            set_alias(p.tt.tree.get(), p.energy,p.luminosity);
          }
        }
      }
      else
      {
        std::cout << "Unable to extract energy from filename: " << file_name << std::endl;
      }
    }
  }
  std::sort(P.begin(),P.end(),
      [](auto & sp1, auto &sp2) {
        return sp1.energy < sp2.energy;
        });

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
      double L = Ngg/10*pow(p.energy/MTAU*0.5,2);
      double dL = dNgg/10*pow(p.energy/MTAU*0.5,2);
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
    double Ntt = Points[i].tt.tree->GetEntries(selection);
    double Ngg = Points[i].gg.tree->GetEntries(selection);
    double xs  = 0;
    double dxs = 0;
    totalNtt += Ntt;
    if(Ngg != 0 )
    {
      xs = Ntt/Ngg;
      dxs = sqrt( Ntt/(Ngg*Ngg) + pow(Ntt/(Ngg*Ngg), 2.0)*Ngg );
    }
    g->SetPoint(i, 0.5*Points[i].energy-MTAU, xs);
    g->SetPointError(i, 0.5*Points[i].energy.error, dxs);
    std::cout << i << " " << 0.5*Points[i].energy-MTAU << "  " << Ngg << "  " << Ntt << "   " <<  xs << std::endl;
  }
  std::cout << "Total number of tau-tau candidates:" << totalNtt << std::endl;
  g->SetMarkerStyle(21);
  g->Draw("ap");

  TCanvas * cacop = new TCanvas("acop","acop");
  cacop->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cacop->cd(i+1);
    Points[i].tt.Draw("acop",selection);
  }
  TCanvas * cptem = new TCanvas("ptem","ptem");
  cptem->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cptem->cd(i+1);
    Points[i].tt.Draw("ptem",selection);
  }
  TCanvas * cp = new TCanvas("p","p");
  cp->Divide(2,3);
  for( int i=0;i<Points.size();++i) 
  { 
    cp->cd(i+1);
    Points[i].tt.tree->Draw("p",selection);
  }
  return g;
}

void mccmp(const std::vector<ScanPoint_t> & P, const char * selection, const char * cut, const char * his="")
{
  std::string title = std::string("mc_data: ") + selection + " cut= " + cut;
  TCanvas * c = new TCanvas("data_mc_diff", title.c_str());
  TChain * chain = new TChain("data_mc_diff","data_mc_diff");
  set_alias((TTree*)chain, MTAU*2,1);
  for(int i=1;i<P.size();++i) chain->Add((TChain*)P[i].tt.tree.get());
  std::string sel(selection);
  sel+=">>h1(";
  sel+=his;
  sel+=")";
  chain->Draw(sel.c_str(),cut,"goff");
  auto hdata  = chain->GetHistogram();
  P[P.size()-1].tt.tree->SetLineColor(kRed);
  sel=selection;
  sel+=">>h2(";
  sel+=his;
  sel+=")";
  P[P.size()-1].tt.tree->Draw(sel.c_str(), cut, "goff");
  auto hmc = P[P.size()-1].tt.tree->GetHistogram();
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
    double Ntt = Points[i].tt.tree->GetEntries(selection);
    double Ngg = Points[i].gg.tree->GetEntries();
    double Nbb = Points[i].bb.tree->GetEntries(BB_SEL.c_str());
    double xs  = 0;
    double dxs = 0;
    totalNtt += Ntt;
    if(Ngg != 0 )
    {
      xs = Ntt/Ngg;
      dxs = sqrt( Ntt/(Ngg*Ngg) + pow(Ntt/(Ngg*Ngg), 2.0)*Ngg );
    }
    g->SetPoint(i, 0.5*Points[i].energy-MTAU, xs);
    g->SetPointError(i, 0.5*Points[i].energy.error, dxs);
    std::cout << i << " " << Points[i].energy/2.0-MTAU << "  " << Ngg << "  " << Ntt << "   " <<  xs << std::endl;
    auto & P = Points[i];
    ofs << std::setw(5) << i <<  std::setw(15) << 20000 << "  " << 10 << std::setw(15) << P.energy/MeV  << std::setw(15) << P.energy.error/MeV;
    ofs << std::setw(10) << 1.256 << " " << std::setw(10) << 0.019;
    ofs << std::setw(10) << Ntt << std::setw(10) << " " << Nbb << "  " << Ngg <<  " " << 1 << std::endl;
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
      Points[i].tt.tree->Draw(var[v],selection);
    }
  }
  return g;
}



struct ChannelSelectionResult_t : public ChannelSelection_t,  public  std::vector<PointSelectionResult_t> 
{
  long Ntt=0; //total number of events for all points
  long Ngg=0; //total number of gg events
  long Nbb=0; //total number of Bhabha events
  double L=0; //total luminosity
  std::string cut;
  ChannelSelectionResult_t(void){}
  ChannelSelectionResult_t(const ChannelSelectionResult_t & ) = default;
  ChannelSelectionResult_t(const ChannelSelection_t & cs, const std::vector<PointSelectionResult_t> & p) : 
    ChannelSelection_t(cs), 
    std::vector<PointSelectionResult_t>(p) {
      accumulate();
  }

  ChannelSelection_t & operator=(const std::vector<PointSelectionResult_t> & P)
  {
    Ntt=0;
    Ngg=0;
    Nbb=0;
    //*this = P;
    clear();
    for(auto & p : P) push_back(p);
    accumulate();
    return *this;
  }
  ChannelSelection_t & operator=(const ChannelSelection_t & cs )
  {
    title      = cs.title;
    cut        = cs.cut;
    //root_title = cs.root_title;
    return *this;
  }
  private:
  void accumulate(void) {
    for(auto & p : *this) {
      Ntt+=p.tt.N;
      Ngg+=p.gg.N;
      Nbb+=p.bb.N;
      L+=p.luminosity;
      //cut = p.cut;
    } 
  }

};

std::string to_string(time_t t, int TZ)
{
  std::string s;
  return s;
};


void print_event_info(ScanPoint_t  & p, const char * selection, const char * title = "")
{
  long n = p.tt.tree->Draw("run:event:time",selection,"goff");
  setenv("TZ", "Asia/Shanghai",1);
  auto runs   = p.tt.tree->GetV1();
  auto events = p.tt.tree->GetV2();
  auto times  = p.tt.tree->GetV3();
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
            if(pid_name == "PI") pid_name = "pi";

            if ( value_name == "chi2_dedx"  || value_name == "delta_tof")  value_name += "_"+pid_name;
            value_name += "["+std::to_string(track)+"]";
            auto rm_zero = [](const std::string  input) -> std::string
            {
              std::string tmp = input;
              std::reverse(tmp.begin(),tmp.end());
              if( tmp.find('.') == std::string::npos) return input;
              std::string result;
              bool stop_skip=false;
              for(int i=0;i<tmp.size();++i)
              {
                if(tmp[i] != '0' || stop_skip)
                {
                  stop_skip = true;
                  result+=tmp[i];
                }
              }
              std::reverse(result.begin(),result.end());
              return result;
            };
            std::string min_str = rm_zero(std::to_string(r.min));
            std::string max_str = rm_zero(std::to_string(r.max));

            //std::string cut = std::to_string(r.min) + " < " + value_name + " && " + value_name  + " < "  +  std::to_string(r.max);
            std::string cut = "("+min_str + "<" + value_name + "&&" + value_name  + "<"  +  max_str +")";
            //std::string cut = buf;
            //std::cout << cut << std::endl;
            if(alias_value=="") alias_value = cut; 
            else alias_value += "&&" + cut;
        };
//        std::cout << " alias: " << name << " = " << alias_value << std::endl;
        data.tt.tree->SetAlias(name.c_str(), alias_value.c_str());
      }
    }
  };
};

void set_kptem(std::vector<ScanPoint_t> & DATA, double kptem)
{
  for( auto & data : DATA)
  {
    data.tt.tree->SetAlias("k_ptem",(std::to_string(kptem)+"*1").c_str());
  };
};

void set_pid_kptem(std::vector<ScanPoint_t> & DATA, const std::vector<ParticleID_t> & PID, double kptem)
{
  set_pid(DATA,PID);
  set_kptem(DATA,kptem);
};


std::vector<PointSelectionResult_t> select(const std::vector<ScanPoint_t> & P, const std::string  sel)
{
  std::vector<PointSelectionResult_t> R(P.size());
  for(int i=0;i<P.size();++i)
  {
    auto & r       = R[i];
    auto & p       = P[i];
    r.tt.N         = p.tt.tree->GetEntries(sel.c_str());
    r.tt.N0mc      = p.tt.N0mc > 0 ? p.tt.N0mc : 0;
    r.bb           = p.bb;
    r.gg           = p.gg;
    r.name         = p.title;
    r.root_name    = p.title;
    r.tex_name     = p.title;
    r.energy       = p.energy;
    r.luminosity   = p.luminosity;
    //r.cut          = sel.c_str();
  }
  return R;
}

ChannelSelectionResult_t select_tmp(const std::vector<ScanPoint_t> & P, const std::string  sel)
{
  //std::vector<PointSelectionResult_t> R(P.size());
  ChannelSelectionResult_t R;
  R.cut = sel;
  long sum=0;
  for(auto & sp : P){
    R.push_back(sp); //copy ScanPoints
    auto & r = R.front();
    //auto & r       = R[i];
    //auto & p       = P[i];
    long N         = sp.tt.tree->GetEntries(sel.c_str());
    r.tt.N         =  N;
    r.tt.N0mc      = sp.tt.N0mc > 0 ? sp.tt.N0mc : 0;
    r.bb           = sp.bb;
    r.gg           = sp.gg;
    r.name         = sp.title;
    r.root_name    = sp.title;
    r.tex_name     = sp.title;
    r.energy       = sp.energy;
    r.luminosity   = sp.luminosity;
    sum+=N;
    std::cout << N <<  " " ;
  }
  std::cout << " sum = " << sum << std::endl;
  return R;
}


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

std::vector<std::vector<PointSelectionResult_t>> draw2(const std::vector<std::vector<ScanPoint_t> *> P, const std::string & var, const std::string & sel, std::string gopt="")
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
  std::vector<std::vector<PointSelectionResult_t>> Result(P.size());
  //find the Scan with bigest number of points
  int nS = std::distance(P.begin(), std::max_element(P.begin(),P.end(),[](auto s1, auto s2){ return s1->size() < s2->size(); }));
  std::vector<ScanPoint_t> & S = *P[nS];
  c->Divide(S.size() ,1);

  auto draw = [&](int s, int color) {
    auto & scan = *P[s];
    Result[s].resize(scan.size());
    for(int i=0;i<scan.size();++i) {
      int nc = i==nS ? 
               1+i :
               1+std::distance(S.begin(), find_best(S, [&scan,&i](const auto & p) { return  std::abs(p.energy-scan[i].energy); } ));
      c->cd(nc);
      auto & r       = Result[s][i]; //point result
      auto & p       = scan[i];      //current point
      p.tt.tree->SetLineColor(s);
      p.tt.tree->SetMarkerColor(s);
      if ( i == nS ) r.tt.N = p.tt.tree->Draw(var.c_str(),sel.c_str(),gopt.c_str());  
      else r.tt.N = p.tt.tree->Draw(var.c_str(),sel.c_str(),(gopt+"SAME").c_str());  
      std::string result_title = p.title  + ": " + std::to_string(r.tt.N) + " events";
      p.tt.tree->GetHistogram()->SetTitle(result_title.c_str());
      //if(p.gg) r.Ngg = p.gg->GetEntries();
      //if(p.bb) r.Nbb = p.gg->GetEntries(BB_SEL.c_str());
      r.name         = p.title;
      r.root_name    = p.title;
      r.tex_name     = p.title;
      r.energy       = p.energy;
      r.luminosity   = p.luminosity;
    }
  };
  for(int s = 0; s < P.size(); ++s)
  {
    draw(s,s+1);
  }
  return Result;
}

std::vector<PointSelectionResult_t> draw(const std::vector<ScanPoint_t> & P, const std::string  var, const std::string  sel, std::string gopt="")
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
    r.tt.N         = p.tt.tree->Draw(var.c_str(),sel.c_str(),gopt.c_str());  
    std::string result_title = p.title  + ": " + std::to_string(r.tt.N) + " events";
    p.tt.tree->GetHistogram()->SetTitle(result_title.c_str());
    //if(p.gg) r.Ngg = p.gg->GetEntries();
    //if(p.bb) r.Nbb = p.gg->GetEntries();
    r.name         = p.title;
    r.root_name    = p.title;
    r.tex_name     = p.title;
    r.energy       = p.energy;
    r.luminosity   = p.luminosity;
  }
  return R;
}

std::vector<PointSelectionResult_t> draw(std::vector<ScanPoint_t> & DATA, const Selection_t & S, int i, const std::string  var,  std::string extracut="", std::string gopt="") 
{
  set_pid(DATA, S.pid);
  std::string cut = S[i].cut;
  if(S.cut!="") cut += " && " + S.cut;
  if(extracut!="")     cut += " && " + extracut;
  return draw(DATA, var, cut, gopt);
};



//fold all points together and draw the result
TH1 *  fold(const std::vector<ScanPoint_t> & P, std::string  var, std::string  Sel, std::string gopt="", int Nbin=100, double Min=0, double Max=0)
{
  //find commond xmin and xmax and number of bins
  std::vector<std::unique_ptr<TH1>> vH(P.size());
  std::vector<double> vweight; 
  //calculate weights
  auto vFc  = [](double v) -> double { 
    if(v==0) return PIALPHA; 
    return PIALPHA/(1. - exp(-PIALPHA/v));
  };

  //tau pair cross section with Coloumb correction
  auto sigma_tree_fc = [&vFc]( double s) -> double
  {
    if( s < 4*MTAU*MTAU) return 0;
    double v= sqrt(1.0 - pow(MTAU*2,2)/s);
    return SIGMA_CONST*(3.-v*v)/2./s*vFc(v);
  };

  for(auto & p : P) 
  {
    double sigma = sigma_tree_fc(p.energy*p.energy);
    if(sigma == 0) sigma  = sigma_tree_fc((2*MTAU)*(2*MTAU));
    //std::cout << p.L << " " <<  sigma << std::endl;
    vweight.push_back(p.luminosity * sigma);
  }
  double wsum = std::accumulate(vweight.begin(), vweight.end(),0.0);
  if(wsum == 0) {
    std::cerr << "ERROR: zero total sum L * sigma" << std::endl;
    std::exit(1);
  }
  for(auto & w : vweight) w*=1e-6;
  std::string sel = Sel;

  auto draw = [&](const ScanPoint_t & p, std::string title, double weight = - 1.0) -> long {
    if(weight>=0) {
      if(Sel=="") sel = std::to_string(weight);
      else sel = std::to_string(weight)+"*("+Sel+")";
    }
    else sel= Sel;
    //std::cout << p.title << " " << sel<< std::endl;
    return p.tt.tree->Draw(title.c_str(),sel.c_str(),"goff");  
  };
  if(Min==Max) { //determine minumum and maximum
    std::vector<double> vmin(P.size()), vmax(P.size());
    std::vector<int> vNbins(P.size());
    for(int i=0;i<P.size();++i)
    {
      auto & p       = P[i];
      std::string htitle = "H2341dodk"+std::to_string(i)+"("+std::to_string(Nbin);
      htitle+=")";
      long N = draw(p, var+">>"+htitle, -vweight[i]);  
      auto h  = (TH1*) p.tt.tree->GetHistogram();
      int Nbins = h->GetNbinsX();
      vNbins[i]  = Nbins;
      double min = h->GetBinCenter(0)-h->GetBinWidth(0)*0.5;
      double max = h->GetBinCenter(Nbins-1)+h->GetBinWidth(Nbins-1)*0.5;
      vmin[i] = min;
      vmax[i] = max;
      //std::cout << htitle << " " << Nbins << " " << min << " " << max << " N = " << N << std::endl;
      delete h;
    }
    Min = *std::min_element(vmin.begin(),vmin.end());
    Max = *std::max_element(vmax.begin(),vmax.end());
   // int Nbins = *std::max_element(vNbins.begin(),vNbins.end());
  }

  //now draw histograms with common bins and range
  for(int i=0;i<P.size();++i) {
    //new TCanvas;
    auto & p       = P[i];
    std::string htitle = "HMoqek231kc"+std::to_string(i);
    std::string hconfig = "("+std::to_string(Nbin)+"," + std::to_string(Min)+"," + std::to_string(Max) + ")";
    //p.tt->Draw((var+">>"+htitle + hconfig).c_str(),sel.c_str());  
    draw(p,(var+">>"+htitle+hconfig).c_str(), -vweight[i]);
    vH[i].reset((TH1*) p.tt.tree->GetHistogram());
  }

  TH1 * H = new TH1F(("his_"+var+std::to_string(HISTO_INDEX)).c_str(), var.c_str(), Nbin, Min,Max);
  for(auto & h : vH) {
    H->Add(h.get());
   // delete h; //remove anused histogram
  }
  HISTO_INDEX++;
  return H;
};


TH1 *  fold_and_draw(const std::vector<ScanPoint_t> & P, std::string  var, std::string  sel, std::string gopt="")
{
  auto c = get_new_tailed_canvas(var + " " + sel);
  auto H = fold(P,var,sel);
  if(gopt == "NORM") H->DrawNormalized();
  else H->Draw();
  return H;
}

void cmp(const Scan_t & S1, const Scan_t & S2, std::string var, std::string sel="", std::string gopt="", int Nbin=100, double Min=0, double Max=0)
{
  auto c  = get_new_tailed_canvas(var + " " + sel);
  auto h1 = fold(S1,var,sel,gopt,Nbin,Min,Max);
  auto h2 = fold(S2,var,sel,gopt,Nbin,Min,Max);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  if(gopt=="NORM") {
    h1->DrawNormalized();
    h2->DrawNormalized("same");
  } else  {
    h1->Draw();
    h2->Draw("same");
  }
};


void cmp(const Scan_t & S1, const Scan_t & S2, const Selection_t & Sel, int i, std::string var, std::string extracut="", std::string gopt = "", int Nbin=100, double Min=0, double Max=0)
{
  std::cout << Sel[i].title << std::endl;
  //std::string sel = Sel.common_cut()  + "&&" + Sel[i].cut + (extracut == "" ? "" : ("&&" + extracut));
  //std::string sel = Sel  && Sel[i].cut + (extracut == "" ? "" : ("&&" + extracut));
  cmp(S1,S2,var,Sel[i] && extracut,gopt,Nbin,Min,Max);
};

TH1 * fold(const Scan_t & D, const Selection_t & Sel,  std::string var, std::string extracut, std::string gopt, int Nbin, double Min, double Max) {
  std::string name = "Hlastfold" + std::to_string(HISTO_INDEX);
  TH1 * H = new TH1F(name.c_str(),name.c_str() ,Nbin,Min,Max);
  for(size_t i =0;i<Sel.size();++i) {
    std::cout << Sel[i].title << ", " << std::flush;
    //std::string sel = Sel.common_cut()  + "&&" + Sel[i].cut + (extracut == "" ? "" : ("&&" + extracut));
    auto h = fold(D, var,Sel[i]&&extracut,gopt,Nbin,Min,Max);
    H->Add(h);
    delete h;
  }
  std::cout << std::endl;
  return H;
}



void cmp(const Scan_t & D1, const Scan_t & D2, const Selection_t & SEL,  std::string var, std::string extracut, std::string gopt, int Nbin, double Min, double Max) {
  std::string title = (var + ":" + extracut);
  auto c  = get_new_tailed_canvas(title.c_str());
  auto h1 = fold(D1,SEL,var,extracut,gopt, Nbin,Min,Max);
  auto h2 = fold(D2,SEL,var,extracut,gopt, Nbin,Min,Max);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h1->SetTitle(title.c_str());
  std::smatch s;
  if( std::regex_match(gopt, s, std::regex("NORM|norm")) ) {
    h1->DrawNormalized("E");
    h2->DrawNormalized("Esame");
  } else  {
    h1->Draw("E");
    h2->Draw("Esame");
  }
}

std::vector<TH1*> cmp(std::vector<std::reference_wrapper<Scan_t>>  SCANS, const Selection_t & SEL,  std::string var, std::string extracut, std::string gopt, int Nbin, double Min, double Max) {
  std::vector<TH1*> H;
  std::string title = (var + ":" + extracut);
  auto c  = get_new_tailed_canvas(title.c_str());
  gStyle->SetOptStat(0);
  std::vector<int> color={kRed, kBlue, kBlack, kGreen+2};
  std::vector<int> line = {1,2,3};
  std::smatch sm;
  bool is_norm = std::regex_match(gopt, sm, std::regex("NORM|norm")); 
  for(unsigned i=0;i<SCANS.size();++i) {
    std::cout << "Proceeding " << title.c_str() << ": ";
    auto h = fold(SCANS[i],SEL,var,extracut,gopt, Nbin,Min,Max);
    h->SetTitle(title.c_str());
    h->SetLineWidth(3);
    h->SetLineColor(color[i % color.size()]);
    h->SetMarkerColor(color[i % color.size()]);
    //h->SetLineStyle(line[(i/color.size())%line.size()]);
    H.push_back(h);
  }
  std::vector<double> ymax;
  std::vector<double> ymin;
  for(unsigned i = 0; i<H.size(); ++i) {
    TH1 * h;
    if(is_norm) {
      h=(TH1*)H[i]->DrawNormalized(gopt.c_str());
    } else {
      H[i]->Draw(gopt.c_str());
      h=H[i];
    }
    if(i==0) gopt+="SAME";
    int bin_max = h->GetMaximumBin();
    int bin_min = h->GetMinimumBin();
    ymax.push_back(h->GetBinContent(bin_max) + 2*h->GetBinError(bin_max));
    ymin.push_back(h->GetBinContent(bin_min) - 2*h->GetBinError(bin_min));
  }
  double max = *std::max_element(ymax.begin(),ymax.end());
  double min = *std::max_element(ymin.begin(),ymin.end());
  std::string title2d=("h2d"+std::to_string(HISTO_INDEX));
  TH2 * h2d = new TH2F(title2d.c_str(),title2d.c_str(), Nbin,Min,Max,1000,min,max);
  h2d->Draw();
  for(auto h : H) {
    if(is_norm) h->DrawNormalized((gopt+ (h == H[0] ? " HIST " : "") ).c_str());
    else  h->Draw(gopt.c_str());
  }
  //for(unsigned i=0;i<H.size(); ++i) {
  //  if(is_norm) H[i]->DrawNormalized((gopt+ (i==0 ? " HIST " : "") ).c_str());
  //  else  H[i]->Draw(gopt.c_str());
  //}
  h2d->SetTitle(var.c_str());
  h2d->GetXaxis()->SetTitle(var.c_str());

  TLegend * l = new TLegend(0.6,0.8,1.0,1.0);
  l->SetHeader(var.c_str());
  for(unsigned i=0;i<SCANS.size();++i) {
    Scan_t  &  s = SCANS[i];
    l->AddEntry(H[i],s[0].type.c_str(), "lp");
  }
  l->Draw();
  return H;
}

TCanvas* cmp(std::vector<std::reference_wrapper<Scan_t>>  SCANS, std::string var, std::string extracut, std::string gopt, int Nbin, double Min, double Max) {
  std::vector<TH1*> H;
  std::string title = (var + ":" + extracut);
  auto c  = get_new_tailed_canvas(title.c_str());
  gStyle->SetOptStat(0);
  std::vector<int> color={kRed, kBlack, kBlue, kGreen+2, kYellow,kCyan};
  std::vector<int> line = {1,2,3};
  std::smatch sm;
  bool is_norm = std::regex_match(gopt, sm, std::regex("NORM|norm")); 
  for(unsigned i=0;i<SCANS.size();++i) {
    auto h = fold(SCANS[i],var,extracut,gopt, Nbin,Min,Max);
    h->SetTitle(title.c_str());
    h->SetLineWidth(3);
    h->SetLineColor(color[i % color.size()]);
    h->SetMarkerColor(color[i % color.size()]);
    h->SetMarkerStyle(i+20);
    h->SetLineStyle(line[(i/color.size())%line.size()]);
    H.push_back(h);
  }
  std::vector<double> ymax;
  std::vector<double> ymin;
  for(unsigned i = 0; i<H.size(); ++i) {
    TH1 * h;
    if(is_norm) {
      h=(TH1*)H[i]->DrawNormalized(gopt.c_str());
    } else {
      H[i]->Draw(gopt.c_str());
      h=H[i];
    }
    if(i==0) gopt+="SAME";
    int bin_max = h->GetMaximumBin();
    int bin_min = h->GetMinimumBin();
    ymax.push_back(h->GetBinContent(bin_max) + 2*h->GetBinError(bin_max));
    ymin.push_back(h->GetBinContent(bin_min) - 2*h->GetBinError(bin_min));
    std::cout << ymax.back() << std::endl;
  }
  double max = *std::max_element(ymax.begin(),ymax.end());
  double min = *std::max_element(ymin.begin(),ymin.end());
  std::string title2d=("h2d"+std::to_string(HISTO_INDEX));
  TH2 * h2d = new TH2F(title2d.c_str(),title2d.c_str(), Nbin,Min,Max,1000,min,max);
  h2d->Draw();
  for(unsigned i=0;i<H.size(); ++i) {
    //if(is_norm) H[i]->DrawNormalized((gopt+ (i==0 ? " HIST " : "") ).c_str());
    H[i]->SetLineWidth(3);
    if(is_norm) H[i]->DrawNormalized((gopt+ " HIST").c_str());
    //if(is_norm) H[i]->DrawNormalized((gopt+ " HIST PLC").c_str());
    else  H[i]->Draw(gopt.c_str());
  }
  h2d->SetTitle(var.c_str());
  h2d->GetXaxis()->SetTitle(var.c_str());

  TLegend * l = new TLegend(0.6,0.8,1.0,1.0);
  for(unsigned i=0;i<SCANS.size();++i) {
    //typeof(SCANS[i]);
    Scan_t  &  s = SCANS[i];
    l->AddEntry(H[i],s[0].type.c_str(), "lp");
  }
  l->Draw();
  c->SaveAs("tmp.root");
  return c;
}

ChannelSelectionResult_t  fold(const std::vector<ChannelSelectionResult_t>  & SR, std::string name = "all")
{
  ChannelSelectionResult_t result;
  result.title = name;
  if(SR.empty()) return result;
  result.resize(SR[0].size());
  result.Ntt=0;
  for(int i=0;i<result.size();++i) 
  {
    auto & rp   = result[i];
    rp.tt.efficiency.value = 0;
    rp.tt.efficiency.error = 0;
    double eps_error2=0;
    for(const auto & r : SR )  
    {
      auto & p    = r[i];
      result.Ntt += p.tt.N;
      rp.tt.efficiency += p.tt.efficiency;
      //rp.eps     += p.eps;
      //eps_error2 += p.eps.error*p.eps.error;
      rp.tt.N      += p.tt.N;
      rp.gg      = p.gg;
      rp.bb      = p.bb;
      rp.energy  = p.energy;
      rp.luminosity = p.luminosity;
    }
    //rp.eps_error = sqrt(eps_error2);
  }
  auto & reference_point = * find_best(result, [](const auto & p) { return  std::abs(p.energy*0.5-MTAU+0.5); } );
  //normalize efficiency correction to threshold efficiency
  for(size_t i=0;i<result.size();++i)
  {
    auto & rp   = result[i];
    rp.tt.effcor = rp.tt.efficiency / reference_point.tt.efficiency.value;
    //rp.effcor.error = rp.eps.error / reference_point.eps;
  }
  reference_point.tt.effcor.error = 0;
  return result;
};

ChannelSelectionResult_t  fold(std::vector<std::vector<ChannelSelectionResult_t>> & SR, std::string name = "all")
{
  std::vector<ChannelSelectionResult_t> VSR;
  for(const auto & vsr : SR) VSR.push_back(fold(vsr));
  return fold(VSR);
};

ChannelSelectionResult_t  fold(const std::vector<ChannelSelectionResult_t>  & s1, const std::vector<ChannelSelectionResult_t>  & s2)
{
  std::vector<ChannelSelectionResult_t> VSR;
  VSR.push_back(fold(s1));
  VSR.push_back(fold(s2));
  return fold(VSR);
};

void print(const ChannelSelectionResult_t & sr, int opt=1 , int first_column_width=10, int last_column_width=5, int  column_width= 8, int vline_width=6)
{
  int N = sr.size(); //number of points
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
      for(int i=0;i<N;++i) std::cout << std::setw(column_width)  << sr[i].name;
      std::cout << std::setw(vline_width) << " │ " << std::setw(last_column_width) << "TOTAL";
      std::cout << std::endl;
      hline();
      std::cout << std::setw(first_column_width+count_utf8_extra_byte(std::string("ΔE,MeV"))) << "ΔE,MeV";
      for(int i=0;i<N;++i) std::cout << std::setw(column_width) << (sr[i].energy*0.5-MTAU)*1e3;
      std::cout << std::setw(vline_width) << " │ ";
      std::cout << std::endl;
      hline();
    default:
      std::cout << std::setw(first_column_width+count_utf8_extra_byte(sr.title)) << sr.title;
      for(int i=0; i < N; i++) std::cout << std::setw(column_width) <<  sr[i].tt.N;
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
  void hline(int N, std::string symb="─", std::string name="")
  {
    int hline_width = title_width + total_width + vline_width + N*data_width+3;
    if(name=="") {
      //int w = width == 0 ? hline_width : width;
      for(int i=0;i<hline_width;++i) std::cout << symb;
    }
    else {
      int f = (hline_width - 2 - name.length())/2;
      for(int  i =0; i< f; ++i) std::cout << symb;
      std::cout << " " << name << " ";
      for(int  i =0; i< f; ++i) std::cout << symb;
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
  print_smth(pts, cfg, [](){ return "ε"; }, [](auto & p) { return p.tt.efficiency.value; }, [](auto & s) { return ""; } );
};

void print_Ntt(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg, std::string title="Ntt")
{
  print_smth(pts, cfg, 
      [&title](){ return title; }, 
      [](auto & p) { return p.tt.N; }, 
      [](auto & s) { 
        long N=0; 
        for(auto & x: s)  N+=x.tt.N;  
        return N; 
        } 
      );
};

void print_Ngg(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg, std::string title="Ngg")
{
  print_smth(pts, cfg, 
      [&title](){ return title; }, 
      [](auto & p) { return p.gg.N; }, 
      [](auto & s) { 
        long N=0; 
        for(auto & x: s)  N+=x.gg.N;  
        return N; 
        } 
      );
};

void print_Nbb(const std::vector<PointSelectionResult_t> pts, PrintConfig_t cfg, std::string title="Nbb")
{
  print_smth(pts, cfg, 
      [&title](){ return title; }, 
      [](auto & p) { return p.bb.N; }, 
      [](auto & s) { 
        long N=0; 
        for(auto & x: s)  N+=x.bb.N;  
        return N; 
        } 
      );
};

template<typename Fdata, typename Ftotal>
void print_smth(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg, Fdata fd, Ftotal ft)
{
  if(SR.empty()) return;
  auto hline = [&SR, &cfg](std::string s="─") { cfg.hline(SR[0].size(),s); };
  hline("━");
  print_head(SR[0],cfg);
  hline();
  auto prn = [&](std::string title, auto & sr) {
    print_smth(sr, cfg, 
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
        return p.tt.N; 
      },
      [](auto & s) { 
        double sum=0;
        for(auto & x: s)  sum+=x.tt.N;  
        return sum; 
      } 
  );
};

void print_Ngg(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  print_smth(SR,cfg, 
      [](auto & p) { 
        return p.gg.N; 
      },
      [](auto & s) { 
        double sum=0;
        for(auto & x: s)  sum+=x.gg.N;  
        return sum; 
      } 
  );
};

void print_Nbb(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  print_smth(SR,cfg, 
      [](auto & p) { 
        return p.bb.N; 
      },
      [](auto & s) { 
        double sum=0;
        for(auto & x: s)  sum+=x.bb.N;  
        return sum; 
      } 
  );
};


void print_efficiency(const std::vector<ChannelSelectionResult_t> SR, PrintConfig_t cfg = PCFG)
{
  if(SR.empty()) return;
  std::string format_str = "%4.3f ± %4.3f";
  cfg.total_width = format_str.length();
  auto hline = [&SR, &cfg](std::string s="─",std::string title="") { cfg.hline(SR[0].size(),s,title); };
  hline("━","REGISTRATION EFFICIENCY");
  print_smth(SR[0], cfg, [](){ return "CHNL/PNT"; }, [](auto & p) { return p.name; }, [](auto & s) { return "AVERAGE"; } );
  hline();
  auto prn = [&](std::string title, auto & sr) {
    print_smth(sr, cfg, 
          [&title](){ return title; }, 
          [&format_str](auto & p) { 
            char buf[1024];
            sprintf(buf,format_str.c_str(), p.tt.efficiency.value*100,p.tt.efficiency.error*100);
            return std::string(buf); 
          },
          [&format_str](auto & P) { 
            double sum=0;
            double error2_sum=0;
            for(auto & p: P)  
            {
              sum+=p.tt.efficiency.value;  
              error2_sum+=p.tt.efficiency.error*p.tt.efficiency.error;
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
  std::string format_str = " %5.4f ± %5.4f ";
  cfg.total_width = format_str.length();
  auto hline = [&SR, &cfg](std::string s="─",std::string title="") { cfg.hline(SR[0].size(),s,title); };
  hline("━","CORRECTION TO EFFICIENCY");
  print_smth(SR[0], cfg, [](){ return "CHNL/PNT"; }, [](auto & p) { return p.name; }, [](auto & s) { return "AVERAGE"; } );
  hline();
  auto prn = [&](std::string title, auto & sr) {
    print_smth(sr, cfg, 
          [&title](){ return title; }, 
          [&format_str](auto & p) { 
            char buf[1024];
            sprintf(buf,format_str.c_str(), p.tt.effcor.value,p.tt.effcor.error);
            if(p.tt.effcor.error == 0 && p.tt.effcor==1) 
            {
              return std::string("      1      ");
            }
            return std::string(buf); 
          },
          [&format_str](auto & P) { 
            double sum=0;
            double error2_sum=0;
            for(auto & p: P)  
            {
              sum+=p.tt.effcor.value;  
              error2_sum+=p.tt.effcor.error*p.tt.effcor.error;
            }
            double average = sum/P.size(); 
            double error = sqrt(error2_sum)/P.size();
            char buf[1024];
            sprintf(buf,format_str.c_str(), average,error);
            return std::string(buf); 
          } 
        );
  };
  auto prn2 = [&](std::string title, auto & sr) {
    print_smth(sr, cfg, 
          [&title](){ return title; }, 
          [&format_str](auto & p) { 
            char buf[1024];
            sprintf(buf,format_str.c_str(), (p.tt.effcor.value-1)*100,p.tt.effcor.error*100);
            if(p.tt.effcor.error == 0 && p.tt.effcor==1) 
            {
              return std::string("      0      ");
            }
            return std::string(buf); 
          },
          [&format_str](auto & P) { 
            double sum=0;
            double error2_sum=0;
            for(auto & p: P)  
            {
              sum+=p.tt.effcor.value;  
              error2_sum+=p.tt.effcor.error*p.tt.effcor.error;
            }
            double average = sum/P.size(); 
            double error = sqrt(error2_sum)/P.size();
            char buf[1024];
            sprintf(buf,format_str.c_str(), (average-1.0)*100,error*100);
            return std::string(buf); 
          } 
        );
  };
  for( auto &  sr : SR) 
     prn(sr.title,sr);
  hline();
  auto f = fold(SR);
  prn("all", f);
  prn2("(cor-1)x100", f);
  hline("━");
};


void print(const std::vector<PointSelectionResult_t> & Points, int opt=1 , int first_column_width=10, int last_column_width=5, int  column_width= 8, int vline_width=6)
{
  ChannelSelectionResult_t sr;
  sr = Points;
  print(sr,0);
}

//void print(const  std::vector<ChannelSelectionResult_t> & SR )
//{
//  PrintConfig_t cfg;
//  print_Ntt(SR,cfg);
//  print_Ngg(SR[0],cfg);
//  print_Nbb(SR[0],cfg);
//};

std::string print_tex(const std::vector<ChannelSelectionResult_t> & SR,std::string ResultTitle="", std::string fit_file="")
{
  std::ostringstream os;
  os << R"(
  \documentclass[a4paper,12pt]{article}
  \usepackage[T2A]{fontenc}
  \usepackage[utf8]{inputenc} 
  \usepackage[english,russian]{babel}
  \usepackage{graphics}
  \usepackage{epsfig}
  \usepackage{amsmath,amssymb}
  \usepackage{ctable}
  \usepackage[
    a4paper, 
    mag=1000,
    nofoot, 
    nohead,
    left=0.5cm, 
    right=0.5cm, 
    top=0.5cm, 
    bottom=0.5cm
  ]{geometry}
  \renewcommand{\arraystretch}{1.2}
  \newcommand{\myC}[1]{\multicolumn{1}{c}{#1}}
  \newcommand{\myCF}[1]{\multicolumn{1}{c}{\hspace{4ex}#1\hspace{4ex}}}
  \newcommand{\myL}[1]{\multicolumn{1}{l}{#1}}
  \newcommand{\myR}[1]{\multicolumn{1}{r}{#1}}
  \begin{document}
  )" << "\n";
  os << R"(\centering \textbf{ \Large {)" << ResultTitle << R"(}}\\)" << "\n";
  os << R"(\today)" << "\n";



  os << R"(\begin{table}[h!])" << "\n";
  os << R"(\caption{Количество отобранных $\tau$ пар})" << "\n";
  os << R"(\centering)" << "\n";
  os << R"(\begin{tabular}[h!]{)";
  os << "c";
  for(int i=0; i< SR[0].size(); ++i) os << "r";
  os << "|"<<"r}" << R"(\\)";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  int col_width=5;
  os << std::setw(20) << "channel";
  for( auto & p : SR[0]) {
    std::string name = R"(\myCF{)" + p.name + R"(})";
    os << " & " << std::setw(col_width) << name;
  }
  os << std::setw(col_width) << " & " << "total"  << R"(\\ \hline)";
  os << "\n";

  os << std::setw(20) << R"($\int L$, $pb^{-1}$)";
  for( auto & p : SR[0]) {
    os << " & ";
    char buf[1024];
    sprintf(buf,"%5.1f", p.luminosity.value);
    os << std::setw(col_width) << buf;
  }
  os << " & ";
  {
    char buf[1024];
    sprintf(buf,"%5.1f", (SR[0].L));
    os << std::setw(col_width) << buf;
  }
  os<< R"(\\ \hline)" <<  "\n";

  os << std::setw(20) << R"($E-M_{\tau}^{\text{PDG}}$, MeV)";
  for( auto & p : SR[0]) {
    os << " & ";
    char buf[1024];
    sprintf(buf,"%3.2f", (p.energy.value*0.5 - MTAU)*1000);
    os << std::setw(col_width) << R"(\myR{)" <<  buf << "}";
  }
  os << " & " << R"(\\ \hline)" <<  "\n";

  

  os << R"(\renewcommand{\arraystretch}{1.1})" << "\n";
  for (auto & sr : SR) {
    std::string title = sr.title; 
    title = sub(title,"μ",R"(\mu)");
    title = sub(title,"π",R"(\pi)");
    title = sub(title,"ρ",R"(\rho)");
    title = "$ " + title + " $";
    os << std::setw(20) << title << " & ";
    for(auto & p : sr) {
      os << std::setw(col_width) << p.tt.N << " & ";
    }
    os << std::setw(col_width) << sr.Ntt << R"(\\)" << "\n";
  }
  os << R"(\hline)" << "\n";
  auto total = fold(SR);
  os << std::setw(20) << "all" << " & "; 
  for( auto & p : total ) {
    os << std::setw(col_width) << p.tt.N << " & ";
  }
  os << std::setw(col_width) << total.Ntt << R"(\\)" << "\n";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  os << R"(\end{tabular})" << "\n";
  os << R"(\end{table})" << "\n";

  os << R"(\begin{table}[h!])" << "\n";
  os << R"(\caption{Эффективность регистрации $\tau$ пар})" << "\n";
  os << R"(\centering)" << "\n";
  os << R"(\resizebox{\textwidth}{!}{)" << "\n";
  os << R"(\begin{tabular}[h!]{)";
  os << "c";
  for(int i=0; i< SR[0].size(); ++i) os << "r";
  //os << "|"<<"r";
  os << "}" << R"(\\)";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  col_width=20;
  os << std::setw(20) << "channel";
  for( auto & p : SR[0]) {
    std::string name = R"(\myC{)" + p.name + "}";
    os << " & " << std::setw(col_width) << name;
  }
  //os << std::setw(col_width) << "average ";
  os << R"(\\ \hline)";
  os << "\n";
  os << R"(\renewcommand{\arraystretch}{1.1})" << "\n";
  for (auto & sr : SR) {
    std::string title = sr.title; 
    title = sub(title,"μ",R"(\mu)");
    title = sub(title,"π",R"(\pi)");
    title = sub(title,"ρ",R"(\rho)");
    title = "$ " + title + " $";
    os << std::setw(20) << title;
    for(auto & p : sr) {
      os << " & ";
      char buf[1024];
      sprintf(buf, "$%4.3f \\pm %4.3f$", p.tt.efficiency.value*100, p.tt.efficiency.error*100);
      os << std::setw(col_width) << buf;
    }
    //os << std::setw(col_width) << sr.Ntt;
    os << R"(\\)" << "\n";
  }
  os << R"(\hline)" << "\n";
  os << std::setw(20) << "all"; 
  for( auto & p : total ) {
      os << " & ";
      char buf[1024];
      sprintf(buf, "$%4.3f \\pm %4.3f$", p.tt.efficiency.value*100, p.tt.efficiency.error*100);
      os << std::setw(col_width) << buf;
  }
  //os << std::setw(col_width) << total.Ntt;
  os << R"(\\)" << "\n";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  os << R"(\end{tabular})" << "\n";
  os << R"(})" << "\n";
  os << R"(\end{table})" << "\n";

  os << R"(\begin{table}[h!])" << "\n";
  os << R"(\caption{Поправка к эффективности регистрации $\tau$ пар})" << "\n";
  os << R"(\centering)" << "\n";
  os << R"(\resizebox{\textwidth}{!}{)" << "\n";
  os << R"(\begin{tabular}[h!]{)";
  os << "c";
  for(int i=0; i< SR[0].size(); ++i) os << "r";
  //os << "|"<<"r";
  os << "}" << R"(\\)";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  col_width=20;
  os << std::setw(20) << "channel";
  for( auto & p : SR[0]){
    std::string name = R"(\myC{)" + p.name + "}";
    os << " & " << std::setw(col_width) << name;
  }
  //os << std::setw(col_width) << "average ";
  os << R"(\\ \hline)";
  os << "\n";
  os << R"(\renewcommand{\arraystretch}{1.1})" << "\n";
  for (auto & sr : SR) {
    std::string title = sr.title; 
    title = sub(title,"μ",R"(\mu)");
    title = sub(title,"π",R"(\pi)");
    title = sub(title,"ρ",R"(\rho)");
    title = "$ " + title + " $";
    os << std::setw(20) << title;
    for(auto & p : sr) {
      os << " & ";
      char buf[1024];
      sprintf(buf, "$%4.3f \\pm %4.3f$", p.tt.effcor.value, p.tt.effcor.error);
      os << std::setw(col_width) << buf;
    }
    //char buf[1024];
    //sprintf(buf, "$%4.3f \\pm %4.3f$", sr.effcor, sr.effcor_error);
    //os << std::setw(col_width) << sr.Ntt;
    os << R"(\\)" << "\n";
  }
  os << R"(\hline)" << "\n";
  os << std::setw(20) << "all"; 
  for( auto & p : total ) {
      os << " & ";
      char buf[1024];
      sprintf(buf, "$%4.3f \\pm %4.3f$", p.tt.effcor.value, p.tt.effcor.error);
      os << std::setw(col_width) << buf;
  }
  //os << std::setw(col_width) << total.Ntt;
  os << R"(\\)" << "\n";
  os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
  os << R"(\end{tabular})" << "\n";
  os << R"(})" << "\n";
  os << R"(\end{table})" << "\n";
  if(fit_file!="")
  {
    auto total = fold(SR);
    os << R"(\begin{table}[h!])" << "\n";
    os << R"(\caption{Количество отобранных $\tau$ пар})" << "\n";
    os << R"(\centering)" << "\n";
    os << R"(\resizebox{\textwidth}{!}{)" << "\n";
    os << R"(\begin{tabular}[h!]{)";
    os << "c";
    for(int i=0; i< SR[0].size(); ++i) os << "r";
    os << "|"<<"r}" << R"(\\)";
    //os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
    int col_width=5;

    auto make_row = [&col_width](
        std::ostream & os, 
        const ChannelSelectionResult_t & sr, 
        std::string title, 
        std::function<std::string(const PointSelectionResult_t & p)> Data,
        std::function<std::string(void)> LastColumn = [](){ return ""; })
    {
      os << std::setw(20) << title;
      for( auto & p : sr) {
        os << " & ";
        os << std::setw(col_width) << Data(p);
      }
      os << " & ";
      os << LastColumn();
      os << R"(\\)" << "\n";
    };


    make_row(os, total, 
        R"(point)", 
        [](const PointSelectionResult_t & p) { 
        //return myfmt(R"(\myCF{%s})", p.name.c_str()); },
        return myfmt(R"(%s)", p.name.c_str()); },
        [&total]() { return "Total"; }
        );
    os<< R"(\hline)" <<  "\n";


    //print luminosity
    make_row(os, total, 
        R"($\int L, pb^{-1}$)", 
        //[](const PointSelectionResult_t & p) { return myfmt(R"(\myR{%4.1f})", p.luminosity.value); },
        [](const PointSelectionResult_t & p) { return myfmt(R"(%4.1f)", p.luminosity.value); },
        [&total]() { return myfmt(R"(%4.1f)", total.L); }
        );

    //print energy for point
    make_row(os, total, 
        R"($E-M_{\tau}$, MeV)", 
        //[](const PointSelectionResult_t & p) { return myfmt(R"(\myR{%3.3f})", (p.energy.value*0.5 - MTAU)*1000.0); }
        [](const PointSelectionResult_t & p) { return myfmt(R"(%3.3f)", (p.energy.value*0.5 - MTAU)*1000.0); }
        );
    os<< R"(\hline)" <<  "\n";

    os << R"(\renewcommand{\arraystretch}{1.1})" << "\n";
    for (auto & sr : SR) {
      std::string title = sr.title; 
      title = sub(title,"μ",R"(\mu)");
      title = sub(title,"π",R"(\pi)");
      title = sub(title,"ρ",R"(\rho)");
      title = "$ " + title + " $";
      make_row(os, sr, title,  
        [](const PointSelectionResult_t & p) { return myfmt(R"(%d)", p.tt.N); },
        [&sr]() { return R"(\textit{)"+std::to_string(sr.Ntt)+"}"; });
    }
    os << R"(\hline)" << "\n";

    //print total event number
    make_row(os, total, 
        "all", 
        [](const PointSelectionResult_t & p) { return myfmt(R"(\textbf{%d})", p.tt.N); },
        [&total]() { return R"(\textbf{)"+std::to_string(total.Ntt)+"}"; }
        );
    os<< R"(\hline)" <<  "\n";

    //print epsilon
    make_row(os, total, 
        R"($\varepsilon$, \%)", 
        [](const PointSelectionResult_t & p) { return myfmt(R"($%5.3f\pm%5.3f$)", p.tt.efficiency.value*100, p.tt.efficiency.error*100); });

    //print epsilon correction
    make_row(os, total, 
        R"($\epsilon^{cor}$)", 
        [](const PointSelectionResult_t & p) { return myfmt(R"($%6.4f\pm%6.4f$)", p.tt.effcor.value, p.tt.effcor.error); });

    make_row(os, total, 
        R"($(\epsilon^{cor}-1)\cdot100\%$)", 
        [](const PointSelectionResult_t & p) { return myfmt(R"($%5.3f\pm%5.3f$)", (p.tt.effcor.value-1)*100.0, p.tt.effcor.error*100.0); });

    //os << R"(\specialrule{.1em}{.05em}{.05em})" << "\n";
    os << R"(\hline)" << "\n";
    os << R"(\end{tabular})" << "\n";
    os << "}\n";
    os << R"(\end{table})" << "\n";


    os << R"(\begin{figure}[h!])" << "\n";
    os << R"(\centering)" << "\n";
    os << R"(\includegraphics[width=0.8\textwidth]{)" << fit_file << "}\n";
    os << R"(\caption{)" <<sub(fit_file,"_",R"(\_)")  << "}\n";
    os << R"(\end{figure})" << "\n";
  }
  os << R"(\end{document})";
  os << "\n";
  return os.str();
};

void make_tex(std::string latex, std::string filename = "scan.tex")
{
  std::ofstream ofs(filename);
  ofs << latex;
  ofs.close();
  char command[65536];
  std::cout << "Compiling latex file..." << filename << std::endl;
  sprintf(command, "pdflatex %s > /dev/null", filename.c_str());
  system(command);
};


ChannelSelectionResult_t select(const Scan_t & P, const ChannelSelection_t & S, std::string extra_cut="")
{
  if(extra_cut!="") extra_cut = " && " + extra_cut;
  ChannelSelectionResult_t r;
  r = S; //copy selection config
  r  = select(P, S.cut + extra_cut);
  print(r,0);
  return r;
};

std::vector<ChannelSelectionResult_t> select(const Scan_t & P, const Selection_t & S, std::string extra_cut="") {
  std::vector<ChannelSelectionResult_t> R(S.size());
  if(extra_cut!="") extra_cut = " && " + extra_cut;
  for(int i=0; i<S.size(); ++i) {
    //std::cout << i << " " << S.name << " " << S.common_cut << std::endl;
    auto & s = S[i];
    auto & r = R[i];
    r = s; //save current selection for channel
    r.cut = s.cut + extra_cut;
    r = select(P, s.cut + extra_cut);
    print(r,i);
  };
  print(fold(R,"all"),-1);
  return R;
};


//ChannelSelectionResult_t select(const std::vector<ScanPoint_t> & DATA, const Selection_t & SEL, int channel, std::string extracut="")
//{
//  return select(DATA, SEL[channel], extracut);
//};

std::vector<ChannelSelectionResult_t> select(const Scan_t & DATA, const Selection_t & SEL, std::vector<int> channels, std::string extracut="") {
  std::vector<ChannelSelectionResult_t> R(channels.size());
  /*
  std::string cut = SEL.common_cut();
  if(extracut != "") {
    if(SEL.common_cut() !="") cut = SEL.common_cut() +"&&"+extracut;
    else cut = extracut;
  }
  if(cut!="") cut = "&&"+cut;
  */

  int j=0;
  for(auto i : channels) {
    auto & s = SEL[i];
    auto & r = R[i];
    r = s;
    r.cut = SEL[i] && extracut;
    r = select(DATA, r.cut);
    print(r,j++);
  }
  print(fold(R,"all"),-1);
  return R;
};


std::vector<PointSelectionResult_t> measure_efficiency(const std::vector<ScanPoint_t> & P, std::string cut)
{
  cut=sub(cut,R"(\s)","");
  std::cout << cut << std::endl;
  auto s = select(P,cut);
  for(auto & sp : s ) {
    auto & eps = sp.tt.efficiency;
    double N0 = sp.tt.N0mc;
    sp.tt.efficiency.value = double(sp.tt.N)/N0;
    sp.tt.efficiency.error = sqrt( eps.value/N0 * ( 1.0 - eps.value));
    //std::cout << sp.tt.N << " " <<  N0 << "  " << sp.tt.efficiency.value << " " << sp.tt.efficiency.error  << std::endl;
    sp.tt.effcor.value = sp.tt.efficiency/s[0].tt.efficiency.value;
    sp.tt.effcor.error = sp.tt.effcor.value*hypot( sp.tt.efficiency.error/sp.tt.efficiency.value, s[0].tt.efficiency.error / s[0].tt.efficiency.value);
  }
  return s;
}

std::vector<PointSelectionResult_t> & set_efficiency(std::vector<PointSelectionResult_t>  & psr,  const std::vector<PointSelectionResult_t>  & eff, long N0)
{
  size_t i0=0;
  const double UNSET = std::numeric_limits<double>::max();
  double chi2_thr = UNSET;
  for(size_t i = 0; i< psr.size(); ++i) {
    auto & r = psr[i];
    int jc; //correspoinding point
    double chi2 = UNSET;
    for( size_t j = 0;j<eff.size(); ++j) {
      if (  std::abs(eff[j].energy-r.energy) < chi2 ) {
        jc = j;
        chi2 = std::abs(eff[j].energy-r.energy); 
      }
    }
    assert ( chi2 == UNSET);
    ibn::valer<double> eps;
    eps.value = (double(eff[jc].tt.N)/N0);
    eps.error = sqrt( eps.value/N0 * ( 1.0 - eps.value));
    //std::cout << i << " " << r.name << " " << "  W=" << r.W << "  eps = " << eps << " +- " << eps_error << std::endl;
    if( std::abs(r.energy*0.5 - MTAU - 0.5) <  chi2_thr )  {
      chi2_thr = std::abs(r.energy*0.5 - MTAU-0.5);
      i0 = i;
    }
    r.tt.efficiency = eps;
  }
  auto & reference_point = * find_best(psr, [](const auto & p) { return  std::abs(p.energy*0.5-MTAU+0.5); } );
  //normalize efficiency correction to threshold efficiency
  for(auto & rp : psr) {
    rp.tt.effcor = rp.tt.efficiency / reference_point.tt.efficiency.value;
    //rp.effcor_error = rp.eps_error / reference_point.eps;
  }
  reference_point.tt.effcor.error = 0;
  return psr;
}

std::vector<ChannelSelectionResult_t> & set_efficiency( std::vector<ChannelSelectionResult_t> & SR, const Scan_t & MC, long N0)
{
  for(auto & sr : SR)
  {
    auto eff = measure_efficiency(MC,sr.cut);
    set_efficiency(sr, eff, N0);
  };
  print_efficiency(SR);
  print_effcor(SR);
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
    scan[i].tt.tree->SetLineColor(color[c]);
    scan[i].tt.tree->SetLineStyle(line[l]);
    scan[i].tt.tree->SetLineWidth(3);
    s+=" >> h"+std::to_string(i)+"("+std::to_string(Nbin)+","+std::to_string(xmin)+","+std::to_string(xmax)+")";
    if(i==0) scan[i].tt.tree->Draw(s.c_str(), selection);
    else scan[i].tt.tree->Draw(s.c_str(),selection,"same");
  }
}






void compare(const ScanPoint_t & p1, const ScanPoint_t &p2, const char * varexp, const char * selection, const char * gopt)
{
  auto  c = new TCanvas;
  c->SetWindowSize(1200, 1800);
  c->Divide(1,2);
  c->cd(1);
  p1.tt.tree->Draw(varexp,selection,gopt);
  c->cd(2);
  p2.tt.tree->Draw(varexp,selection,gopt);
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
      if ( fabs(P.energy-mc.energy) < energy_distance )
      {
        std::cout << "mc.effcor = " << mc.tt.effcor << " " << P.energy << " " << mc.energy << std::endl;
        P.tt.effcor = mc.tt.effcor;
        energy_distance = fabs(P.energy-mc.energy);
      }
    }
  }

  std::cout << "#Correct efficiency: " << std::endl;
  std::cout << "#" << std::setw(4) << " " << std::setw(10) << "W,MeV" << std::setw(10) << "E, MeV" << std::setw(15) << "E-Mtau, MeV" << std::setw(10) << "Ntt" << std::setw(10) << "effcor" << std::endl;
  for( int i = 0;i<D.size(); i++)
  {
    std::cout << std::setw(5) << i+1 << std::setw(10) << D[i].energy/MeV << std::setw(10) << D[i].energy*0.5/MeV << std::setw(15) << (D[i].energy*0.5-MTAU)/MeV <<  std::setw(10) << D[i].tt.N << std::setw(10) << D[i].tt.effcor << std::endl;
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

void compare(std::vector<ScanPoint_t*> SP, std::vector<std::string> varexp, std::string selection, std::string gopt)
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
      sp->tt.tree->SetLineWidth(1);
      set_color_and_line(sp->tt.tree.get(), draw_index);
      c->cd(index+1);
      sp->tt.tree->Draw(ve.c_str(),selection.c_str(), gopt.c_str());
      auto h = sp->tt.tree->GetHistogram();
      h->SetName(("h"+sp->title+std::to_string(index)).c_str());
      h->SetTitle((sp->title).c_str());
      hs->Add(h);
      index++;
      l->AddEntry(sp->tt.tree.get(),sp->title.c_str());
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

TCanvas*  compare2(TTree * t1, TTree * t2, std::string var, std::string sel, std::string gopt="", std::string H="", int bin=0,std::string title="", std::string xaxis_title="")
{
  //std::cout << sel << std::endl;
  auto c = new TCanvas;
  t1->SetLineColor(kRed);
  t1->SetLineWidth(2);
  TH1 * h1;
  TH1 * h2;
  std::string name;
  if(H!="") {
    if(bin==0) name =">>"+H+"1";
    else name = ">>H1("+std::to_string(bin)+")";
  }
  t1->Draw((var+name).c_str(), sel.c_str(),gopt.c_str());
  h1 = (TH1*)t1->GetHistogram()->Clone(sub(name, R"([()\*])","_").c_str());
  bin = h1->GetNbinsX();
  double I1 =h1->Integral(); 
  t2->SetLineColor(kBlue);
  t2->SetLineWidth(2);
  if(H!="") name=">>"+H+"2("+std::to_string(bin)+")";
  t2->Draw((var+name).c_str(), sel.c_str(),(gopt+"same").c_str());  
  h2 = (TH1*)t2->GetHistogram()->Clone(sub(name, R"([()\*])","_").c_str());
  double prob = h1->KolmogorovTest(h2);
  std::cout << var << ": " << prob << std::endl;
  double I2 = h2->Integral();
  h2->Scale(I1/I2);
  h1->Draw();
  h1->SetTitle(var.c_str());
  if(title=="") title = var;
  if(xaxis_title=="") xaxis_title = var;
  h1->SetTitle(title.c_str());
  h1->GetXaxis()->SetTitle(xaxis_title.c_str());
  h2->Draw("same");
  TLatex *  l = new TLatex(0.6,0.95,("prob_{Kolmg} = "+ myfmt("%4.2f",prob)).c_str());
  l->SetNDC(kTRUE);
  l->Draw();
  return c;
};

TCanvas*  compare2(ScanPoint_t & P1, ScanPoint_t & P2, std::string var, std::string sel, std::string gopt="", std::string H="", int bin=0,std::string title="", std::string xtitle="")
{
  return compare2(P1.tt.tree.get(), P2.tt.tree.get(), var, sel,gopt, H, bin,title, xtitle);
}

TCanvas * compare2(ScanPoint_t & P1, ScanPoint_t & P2, std::string var, const Selection_t & SEL, std::string gopt="", std::string H="", int bin=0, std::string title="", std::string xtitle="") {
  //std::string global_cut = SEL.common_cut();
  //std::string or_cut="(";
  //for(auto  it=SEL.begin();it!=SEL.end();++it) {
  //  if(it==SEL.begin()) or_cut+="("+it->cut +")";
  //  else or_cut+= "||("+it->cut+")";
  //}
  //or_cut+=")";
  std::string cut;
  for(auto & sel: SEL) cut=cut||sel;
  std::cout << "total cut = " << cut << std::endl;
  if(H=="") H=var;
  return compare2(P1,P2,var,cut,gopt,H,bin,title, xtitle);
}

TCanvas * compare3(ScanPoint_t & P1, ScanPoint_t & P2, std::string var, std::string subregex, const Selection_t & SEL, std::string gopt="", std::string H="", int bin=0, std::string title="", std::string xtitle="")
{
  std::string global_cut = SEL.common_cut();
  std::string or_cut="(";
  for(auto  it=SEL.begin();it!=SEL.end();++it)
  {
    if(it==SEL.begin()) or_cut+="("+it->cut +")";
    else or_cut+= "||("+it->cut+")";
  }
  or_cut+=")";
  std::string cut = global_cut+"&&"+or_cut;
  std::cout << "compare3 cut = " << cut << std::endl;
  std::cout << "remove var from cut " << std::endl;
  std::string r = var + R"(\s*>\s*\d*)";
  std::cout << "regex = " << r << std::endl;
  std::string result = sub(cut,subregex,"1");
  std::cout << " result = " << result << std::endl;
  return compare2(P1,P2,var,result,gopt,H,bin,title, xtitle);
}

TCanvas*  compare2(Scan_t & S1, Scan_t & S2, std::string var, const Selection_t & SEL, std::string gopt="", std::string H="", int bin=0,std::string title="", std::string xtitle="")
{
  if(S1.empty() || S2.empty()) return nullptr;
  std::string global_cut = SEL.common_cut();
  std::string or_cut="(";
  for(auto  it=SEL.begin();it!=SEL.end();++it)
  {
    if(it==SEL.begin()) or_cut+="("+it->cut +")";
    else or_cut+= "||("+it->cut+")";
  }
  or_cut+=")";
  std::string cut = global_cut+"&&"+or_cut;
  auto make_tree = [](TTree * tt) -> TTree*  {
    TTree * tree = tt->CloneTree(0); //copy structure
    tt->GetListOfClones()->Remove(tree); //detach new tree from initial one
    tree->ResetBranchAddresses(); //reset brunch address
    return tree;
  };
  TTree * tree1 = make_tree(S1[0].tt.tree.get());
  TTree * tree2 = make_tree(S2[0].tt.tree.get());
  auto fill_tree = [](TTree* tree, Scan_t & S, std::string cut)
  {
    for(auto & s : S) {
      auto t = s.tt.tree->CopyTree(cut.c_str());
      t->CopyAddresses(tree);
      for(Long64_t i = 0; i < t->GetEntriesFast(); ++i) {
        if(t->GetEntry(i) < 0) break;
        tree->Fill();
      }
    }
  };
  fill_tree(tree1,S1,cut);
  fill_tree(tree2,S2,cut);
  return compare2(tree1, tree2, var, "", gopt,H,bin,title);
};


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  make_compare
 *  Description:  Compare two data samples according to Selection. Histogram will have
 *  number of bins "bin". If bin==0 then autobinning is used.
 * =====================================================================================
 */
void make_compare(ScanPoint_t & P1, ScanPoint_t & P2, Selection_t & sel, int bin = 0)
{
  auto cmp = [&](std::string var,std::string title="", std::string xtitle="") {
    compare2(P1,P2,var,sel,"",var,bin, title,xtitle)->SaveAs(("data_mc_"+var+".pdf").c_str());
  };
  cmp("p","p","GeV");
  cmp("pt","p_{t}", "GeV");
  cmp("tof","TOF", "ns");
  cmp("M2", "M_{inv}^{2}","GeV^{2}");
  cmp("cos(theta)", "cos#theta", "cos#theta");
  cmp("cos_theta_mis2", "cos#theta_{mis}", "cos#theta_{mis}");
  cmp("ptem","ptem","ptem");
  cmp("Mpi0", "M_{#pi^{0}}", "GeV");
  cmp("Mrho", "M_{#rho}", "GeV");
}

/* This compare of data and MC fold all points together */
void make_compare2(std::vector<std::reference_wrapper<Scan_t>>  SCANS, Selection_t & SEL, int bin)
{
  gStyle->Reset("Pub");
  auto compare = [&](std::string var, int Nb, double min,double max, std::string extracut="") {
    cmp(SCANS,SEL, var,extracut,"NORM",Nb, min,max);
  };
  //compare("ptem",bin,0.1,1.1);
  //compare("p",bin,0.1,1.1);
  //compare("pt",bin,0,1.1);
  //compare("tof",bin,0,6);
  //compare("cos(theta)", bin,-1,1);
  //compare("cos_theta_mis2", bin,-1,1);
  //compare("Mpi0", bin, 0.138,0.141,"Npi0==1");
  //compare("Mrho[0]",bin,0.5,1.1,"PI0 && Npi0==1");
  compare("acop",bin,0,TMath::Pi(),"");
  compare("acol",bin,0,TMath::Pi(),"");
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


void save(const ChannelSelectionResult_t & sr, std::string  filename="scan.txt", std::string default_lum="")
{
  std::cout << "Saving selection: " << sr.title << " to file: " << filename << std::endl;
  std::stringstream os;
  char buf[65535];
  sprintf(buf,"%5s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s",
      "#","L,nb^-1","dL,nb^-1","W,MeV", "dW,MeV","Sw,MeV","dSw,MeV","Ntt","Nbb","Ngg","effcor");
  os << buf << '\n';
  int idx=0;
  for(const auto & p : sr) {
    idx++;
    std::string point_name = p.name == "" ? std::to_string(idx) : p.name;
    double lum = p.luminosity*1000.0;
    double lum_error = 10;
    if(default_lum == "bb" || default_lum=="ee") {
      std::cout << "Use Bhabha + Monte Carlo as default luminosity " << std::endl;
      lum = p.bb.luminosity.value;
      lum_error = p.bb.luminosity.error;
    }
    if(default_lum == "gg") {
      std::cout << "Use gamma gamma + Monte Carlo as default luminosity " << std::endl;
      lum = p.gg.luminosity.value;
      lum_error = p.gg.luminosity.error;
    }
    sprintf(buf,"%5s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10ld %10ld %10ld %10.5f %10.5f",
        point_name.c_str(),
        lum,           lum_error,
        p.energy.value/MeV,        p.energy.error/MeV,
        //1.258,             0.017,
        1.258,             0.060,
        p.tt.N,
        p.bb.N,        p.gg.N,
        p.tt.effcor.value, p.tt.effcor.error);
    os << buf <<'\n';
  }
  std::cout << os.str();
  std::ofstream ofs(filename);
  ofs << os.str();
}


void fit(const ChannelSelectionResult_t & sr, std::string  filename="scan.txt", std::string title="", std::string default_lum = "", bool wait = false)
{
  save(sr,filename,default_lum);
  char command[65536];
  std::string basename = sub(filename,R"(\..+)", "");
  std::string output_file = basename + "_fit";
  sprintf(command, (TAUFIT + " --title='taufit: %s' '%s' --output '%s' & ").c_str(), title.c_str(), filename.c_str(),output_file.c_str());
  system(command);
}

void fit(const std::vector<ChannelSelectionResult_t> & sr, std::string  filename="scan.txt", std::string title="", std::string default_lum = "", bool wait = false)
{
  fit(fold(sr),filename,title, default_lum, wait);
}

#include <TH1F.h>
void find_optimal_depth(ScanPoint_t & sp, double kmin = 0, double kmax=100, double dk=0.1)
{
  TGraph * g =new TGraph;
  for(double k = kmin; k<kmax; k+=dk)
  {
    char buf[1024];
    sprintf(buf,"depth[0]-p[0]*%f >> h", k);
    sp.tt.tree->Draw(buf,"depth[0]>0 && depth[0]<100 && abs(pid[0])==13 && p[0]<1.1 && depth[0]-56.76*p[0] >-40 && depth[0]-56.76*p[0]<-10");
    TH1 * h = sp.tt.tree->GetHistogram();
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
  auto t = sp.tt.tree;
  c->cd(1);
  t->SetMarkerSize(1);
  TLegend * leg = new TLegend(0.8,0.8,1.0,1.0);
  t->SetMarkerStyle(7);
  t->SetMarkerColor(kRed);
  t->Draw("depth[0]:p[0]","depth[0]>0 && depth[0]<100 && abs(pid[0]) == 13");
  leg->AddEntry(t->GetHistogram(),"muon", "lp");
  //t->SetMarkerStyle(21);
  t->SetMarkerColor(kGreen);
  t->Draw("depth[0]:p[0]","depth[0]>0 && depth[0]<100 && abs(pid[0]) == 211", "same");
  leg->AddEntry(t->GetHistogram(),"pion", "lp");
  //t->SetMarkerStyle(21);
  t->SetMarkerColor(kBlue);
  t->Draw("depth[0]:p[0]","depth[0]>0 && depth[0]<100 && abs(pid[0]) == 11","same");
  leg->AddEntry(t->GetHistogram(),"electron", "lp");
  auto l  = new TLine(0.6, k0 + 0.6*k,  1.0, k0+1.0*k);
  l->SetLineWidth(3);
  l->Draw();
  leg->Draw();

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
  auto t = sp.tt.tree;
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
  TLegend * l = new TLegend(0.8,0.8,1.0,1.0);
  draw(var, cut + "abs(pid[0]) ==  11", opt +     "", kBlue,  21, 2);
  l->AddEntry(t->GetHistogram(),"electron", "l");
  draw(var, cut + "abs(pid[0]) ==  13", opt + "SAME", kRed,   21, 2);
  l->AddEntry(t->GetHistogram(),"muon", "l");
  draw(var, cut + "abs(pid[0]) == 211", opt + "SAME", kGreen-3, 21, 2);
  l->AddEntry(t->GetHistogram(),"pion", "l");
  l->Draw();
};

void find_Ep(ScanPoint_t &sp, int Nbin=0, double parmin=0, double parmax=0)
{
  std::string var="Ep[0]";
  std::string varexpr = var;
  std::string cut = "abs(pid[0])==11";
  if(Nbin>0) varexpr += " >> h("+std::to_string(Nbin)+")";
  if(parmin!=0) cut += " && " + var + " > " + std::to_string(parmin);
  if(parmax!=0) cut += " && " + var + " < " + std::to_string(parmax);
  sp.tt.tree->Draw(var.c_str(),cut.c_str());
  TH1 * h = sp.tt.tree->GetHistogram();
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


void draw_momentum(Scan_t & mc) {
  gStyle->Reset("Pub");
  mc[0].tt.tree->SetLineColor(kGreen+3);
  mc[0].tt.tree->SetLineWidth(3);
  mc[0].tt.tree->Draw("p[0]","pid[0]==-211 & Nn==0","NORM HIS");
  auto l = new TLegend(0.9,0.9,1.0,1.0);
  //l->AddEntry((TObject*)nullptr,"E=M_{#tau}",""); 

  mc[0].tt.tree->GetHistogram()->GetXaxis()->SetTitle("p, GeV");
  l->AddEntry(mc[0].tt.tree->GetHistogram(),"#pi");

  mc[0].tt.tree->SetLineColor(kBlue);
  mc[0].tt.tree->Draw("p[0]","pid[0]==11 & Nn==0","NORM SAME HIS");
  l->AddEntry(mc[0].tt.tree->GetHistogram(),"e");
  mc[0].tt.tree->SetLineColor(kRed);
  mc[0].tt.tree->Draw("p[0]","pid[0]==13 & Nn==0","NORM SAME HIS");
  l->AddEntry(mc[0].tt.tree->GetHistogram(),"#mu");

  auto l2 = new TLegend(0.9,0.9,1.0,1.0);
  //l2->AddEntry((TObject*)nullptr,"E=M_{#tau}+23 MeV",""); 
  mc[4].tt.tree->SetLineWidth(3);
  mc[4].tt.tree->SetLineStyle(7);
  mc[4].tt.tree->SetLineColor(kGreen+3);
  mc[4].tt.tree->Draw("p[0]","pid[0]==-211 & Nn==0","NORM SAME HIS");
  l2->AddEntry(mc[4].tt.tree->GetHistogram(),"#pi");
  mc[4].tt.tree->SetLineColor(kBlue);
  mc[4].tt.tree->Draw("p[0]","pid[0]==11 & Nn==0","NORM SAME HIS");
  l2->AddEntry(mc[4].tt.tree->GetHistogram(),"e");
  mc[4].tt.tree->SetLineColor(kRed);
  mc[4].tt.tree->Draw("p[0]","pid[0]==13 & Nn==0","NORM SAME HIS");
  l2->AddEntry(mc[4].tt.tree->GetHistogram(),"#mu");
  l->Draw();
  l2->Draw();
}

void draw_Mpi0(Scan_t & mc, Scan_t & data) {
  gStyle->Reset("Pub");
  mc[0].tt.tree->SetLineColor(kRed);
  mc[0].tt.tree->SetLineWidth(3);
  mc[0].tt.tree->Draw("Mpi0[0]","Npi0==1 && 0.1 < Mpi0[0] && Mpi0[0]<0.2 && (PI0 || PI1)", "NORM HIST");
  auto l = new TLegend(0.9,0.9,1.0,1.0);
  l->AddEntry(mc[0].tt.tree->GetHistogram(),"MC");
  data[0].tt.tree->SetLineColor(kBlue);
  data[0].tt.tree->SetLineWidth(3);
  //data[0].tt->Draw("Mpi0[0]","Npi0==1 && 0.1 < Mpi0[0] && Mpi0[0]<0.2 && (PI0 || PI1)", "SAME NORM HIST");
  data[0].tt.tree->Draw("Mpi0[0]","Npi0==1 && 0.1 < Mpi0[0] && Mpi0[0]<0.2", "SAME NORM HIST");
  l->AddEntry(data[0].tt.tree->GetHistogram(),"DATA");
  l->Draw();
}

void draw_Mrho(Scan_t & mc, Scan_t & data) {
  gStyle->Reset("Pub");
  mc[0].tt.tree->SetLineColor(kRed);
  mc[0].tt.tree->SetLineWidth(3);
  const char * var = "abs(Mrho[0]-0.775) < abs(Mrho[1]-0.775) ? Mrho[0] : Mrho[1]";
  const char * sel = "Npi0==1 && 0.12 < Mpi0[0] && Mpi0[0]<0.14 && Nc==2";
  mc[0].tt.tree->Draw(var,sel, "NORM HIST");
  mc[0].tt.tree->GetHistogram()->GetXaxis()->SetTitle("M_{#pi#gamma#gamma}, GeV");
  auto l = new TLegend(0.9,0.9,1.0,1.0);
  l->AddEntry(mc[0].tt.tree->GetHistogram(),"MC");

  data[0].tt.tree->SetLineColor(kBlue);
  data[0].tt.tree->SetLineWidth(3);
  data[0].tt.tree->Draw(var,sel, "SAME NORM HIST");
  l->AddEntry(data[0].tt.tree->GetHistogram(),"DATA");
  l->Draw();
}





template<typename FCX>
void read_cross_section(std::string filename, std::vector<ScanPoint_t> & P, FCX  fcx) {
  std::ifstream ifs(filename);
  if(!ifs) { 
    std::cerr << "Unable to open file with cross section" << std::endl;
    return;
  }
  std::vector<ScanPoint_t>  tmpSP;
  double W; //in GeV
  double sigma; //in nanobarn
  double sigma_error; // in nanobarn
  std::string tmp; //this is for skip +-
  std::cout << std::setw(10) << "Wmc, GeV" << std::setw(15) << "W in file, GeV" << std::setw(10) << "sigma, nb" << std::setw(10) << "error, nb" << std::endl;
  while(ifs >> W >> sigma >> tmp >> sigma_error >> tmp)  {
    ScanPoint_t sp;
    sp.energy = W;
    fcx(sp).cross_section.value = sigma;
    fcx(sp).cross_section.error = sigma_error;
    tmpSP.push_back(sp);
  }

  for(auto & sp : P) {
    auto & p  = *std::min_element(tmpSP.begin(), tmpSP.end(), [&sp](auto a, auto b){ return fabs(a.energy-sp.energy)<fabs(b.energy-sp.energy); } );
    fcx(sp).cross_section.value = fcx(p).cross_section;
    fcx(sp).cross_section.error = fcx(p).cross_section.error;
    if(fabs(sp.energy-p.energy)>1*MeV) std::cerr << "WARNING: To big difference in energy points" << std::endl;
    std::cout << std::setw(10) << sp.energy << std::setw(15) << p.energy << std::setw(10) << fcx(sp).cross_section.value << std::setw(10) << fcx(sp).cross_section.error << std::endl;
  }
}

void read_bhabha_cross_section(std::string filename, std::vector<ScanPoint_t> & P) {
  std::cout << "Cross section for bhabha process: " << std::endl;
  read_cross_section(filename, P, [](ScanPoint_t & sp) -> DataSample_t & { return sp.bb; } );
}
void read_gg_cross_section(std::string filename, std::vector<ScanPoint_t> & P) {
  std::cout << "Cross section for gamma gamma process: " << std::endl;
  read_cross_section(filename, P, [](ScanPoint_t & sp) -> DataSample_t & { return sp.gg; } );
}

void read_mumu_or_pipi_cross_section(std::string filename, std::vector<ScanPoint_t> & P) {
  std::cout << "Cross section for mumu or pipi: " << std::endl;
  read_cross_section(filename, P, [](ScanPoint_t & sp) -> DataSample_t & { return sp.tt; } );
}

void read_pipi_cross_section(std::string filename, std::vector<ScanPoint_t> & P) {
  std::cout << "Cross section for mumu or pipi: " << std::endl;
  read_cross_section(filename, P, [](ScanPoint_t & sp) -> DataSample_t & { return sp.tt; } );
  for(auto & p : P) {
    //there is no cross section for pi+pi- channel in Babayaga generator
    //so I got it from mu+mu- cross section just multiply it by factor of quarks color
    p.tt.cross_section*=3.0; 
  }
}

void read_hadron_cross_section(std::string filename, std::vector<ScanPoint_t> & P) {
  std::cout << "Cross section for hadrons: " << std::endl;
  read_cross_section(filename, P, [](ScanPoint_t & sp) -> DataSample_t & { return sp.tt; } );
  for(auto & p : P) {
    //there is no cross section for pi+pi- channel in Babayaga generator
    //so I got it from mu+mu- cross section just multiply it by factor of quarks color
    //p.tt.cross_section*=22./31.; 
  }
}
void read_tau_cross_section(std::string filename, std::vector<ScanPoint_t> & P) {
  std::cout << "Cross section for hadrons: " << std::endl;
  read_cross_section(filename, P, [](ScanPoint_t & sp) -> DataSample_t & { return sp.tt; } );
  for(auto & p : P) {
    //there is no cross section for pi+pi- channel in Babayaga generator
    //so I got it from mu+mu- cross section just multiply it by factor of quarks color
    p.tt.cross_section*=1.0; 
  }
}

void read_galuga_cross_section(std::string filename, std::map<std::string, Scan_t> & G) {
  std::ifstream ifs(filename);
  if(!ifs) { 
    std::cerr << "Unable to open file with cross section\n";
    return;
  }
  std::clog << "Reading cross section for GALUGA background.\n";
  auto set_cross_section = [](Scan_t & d, auto W, auto cx ) {
    auto & p  = *std::min_element(d.begin(), d.end(), [&d,W](auto a, auto b){ return fabs(a.energy-W)<fabs(b.energy-W); } );
    p.tt.cross_section = cx;
  };
// W, GeV    EEee                            EEuu                        EEpipi                                    EEkk                                  
  ibn::valer<double> E;
  ibn::valer<double> ee;
  ibn::valer<double> uu;
  ibn::valer<double> pipi;
  ibn::valer<double> KK;
  ifs.ignore(4096,'\n');
  char buf[65353];
  sprintf(buf, "%10s %15s %15s  %15s %15s   %15s %15s   %15s %15s\n", "W,GeV", "eeee, nb", "error", "eeuu, nb", "error", "eepipi, nb", "error", "eeKK, nb", "error");
  std::cout << buf;
  while( ifs >> E.value  >> ee.value >> ee.error >> uu.value >> uu.error >> pipi.value >> pipi.error >> KK.value >> KK.error ) {
    E.error=0;
    sprintf(buf, "%10.6f %15.6f %15.6f  %15.6f %15.6f   %15.6f %15.6f   %15.6f %15.6f\n", E.value, ee.value, ee.error, uu.value, uu.error, pipi.value, pipi.error, KK.value, KK.error);
    std::cout << buf;
    for( auto & [ channel, data ] : G ) {
      set_cross_section(G["ee"], E, ee);
      set_cross_section(G["uu"], E, uu);
      set_cross_section(G["pipi"], E, pipi);
      set_cross_section(G["KK"], E, KK);
    }
  }
}

template<typename Projector>
void me(Scan_t & scan, long N0_MC, Projector  proj, std::string sel="") {
  for(auto & sp : scan) {
    auto & ds = proj(sp);
    ds.N0mc = N0_MC;
    ds.Nmc = ds.tree->GetEntries(sel.c_str());
    ds.efficiency.value = double(ds.Nmc)/ds.N0mc;
    ds.efficiency.error = sqrt( ds.efficiency.value * ( 1.0 - ds.efficiency.value )/ds.N0mc );
    //std::cout <<  sp.W << "  " << ds.Nmc << " " << ds.N0mc << " " << ds.efficiency.value << "  " << ds.efficiency.error  <<   "   "  << ds.cross_section.value  <<  std::endl;
  }
}

template<typename Projector>
void measure_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel, Projector proj) {
  me(mc, N0_MC, proj, sel);
  //find close point and select
  for(auto & sp :data) {
    auto & p  = *std::min_element(mc.begin(), mc.end(), [&sp](auto a, auto b){ return fabs(a.energy-sp.energy)<fabs(b.energy-sp.energy); } );
    if(fabs(sp.energy-p.energy)>1*MeV) std::cerr << "WARNING: To big difference in energy points" << std::endl;
    auto & l = proj(sp);
    auto & m = proj(p);
    l.cross_section.value  = m.cross_section.value*pow(m.energy.value/l.energy.value, 2.0);
    l.cross_section.error  = l.cross_section.value*std::hypot( m.cross_section.error/m.cross_section.value, 2.0*l.energy.error/l.energy.value);
    l.efficiency = m.efficiency;
    l.Nmc = m.Nmc;
    l.N0mc = m.N0mc;
    l.effcor = m.effcor;
    l.N = l.tree->GetEntries(sel.c_str());
    double vis_cx = (l.efficiency*l.cross_section); //visible_cross_section;
    double vis_cx_error = std::hypot(l.cross_section.error*l.efficiency, l.cross_section*l.efficiency.error);
    //std::cout << vis_cx << "  " << vis_cx_error << std::endl;
    l.luminosity.value  = l.N / vis_cx;
    l.luminosity.error = l.luminosity.value*sqrt(  1./l.N  +  pow(vis_cx_error/vis_cx,2.0) );
    //std::cout << sp.W << "  " << l.Nmc << " " << l.N0mc << " " << l.efficiency.value << "  " << l.efficiency.error << "  " << l.luminosity.value << " " << l.luminosity.error <<  std::endl;
  }
}

void measure_bhabha_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel) {
  measure_luminosity(data,mc,N0_MC,sel, [](ScanPoint_t &sp) -> DataSample_t & { return sp.bb; } );
}

void measure_gg_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel) {
  measure_luminosity(data,mc,N0_MC,sel, [](ScanPoint_t &sp) -> DataSample_t & { return sp.gg; } );
}

void measure_luminosity(Scan_t & data, Scan_t & bb, Scan_t & gg, long N0_MC) {
  measure_gg_luminosity(data,gg,N0_MC, GG_SEL);
  measure_bhabha_luminosity(data,bb,N0_MC, BB_SEL);
}

//measure luminosity just for monte carlo only
template<typename Projector>
void measure_efficiency(Scan_t & mc, Projector proj, std::string sel) {
  for(auto & sp : mc) {
    //auto & ds = proj(sp);
    auto & ds = std::invoke(proj,sp);
    ds.Nmc = ds.tree->GetEntries(sel.c_str());
    ds.efficiency.value = double(ds.Nmc)/ds.N0mc;
    ds.efficiency.error = sqrt( ds.efficiency.value * ( 1.0 - ds.efficiency.value )/ds.N0mc );
  }
}

void measure_ee_luminosity(Scan_t & mc, std::string sel) {
  //measure_luminosity(mc,sel, [](ScanPoint_t &sp) -> DataSample_t & { return sp.bb; } );
  measure_efficiency(mc,&ScanPoint_t::bb, sel);
}

void measure_gg_luminosity(Scan_t & mc, std::string sel) {
  //measure_luminosity(mc,sel, [](ScanPoint_t &sp) -> DataSample_t & { return sp.gg; } );
  //measure_luminosity(mc,sel, &ScanPoint_t::gg);
  measure_efficiency(mc,&ScanPoint_t::gg, sel);
}

void measure_tt_efficiency(Scan_t & mc, std::string sel) {
  measure_efficiency(mc,&ScanPoint_t::tt, sel);
}

void set_luminosity(Scan_t & DATA, Scan_t & MC ) {
  for(auto & mc : MC ) {
    auto & p  = *std::min_element(DATA.begin(), DATA.end(), [&mc](auto a, auto b){ return fabs(a.energy-mc.energy)<fabs(b.energy-mc.energy); } );
    mc.bb = p.bb;
    mc.gg = p.gg;
  }
}

void set_luminosity(Scan_t & DATA, std::map<std::string, Scan_t>  & G ) {
  for(auto & [name, s] : G ) {
    set_luminosity(DATA,s);
  }
}

void print_luminosity(Scan_t & data) {
  //char buf[1024*16];
  printf("%5s %10.6s %10.6s %10.6s %10.6s %10s %10.6s %10.6s %10s %10s\n",
      "point", "E,GeV", "Lonline,pb^-1","Lgg, pb^-1", "dLgg,pb^-1", "Ngg","Lbb,pb^-1", "dLbb,pb^-1", "Nbb", "Lbb/Lgg"
      );
  ibn::valer<double> Lgg={0,0};
  ibn::valer<double> Lee={0,0};
  double Lonline=0;
  long Ngg=0;
  long Nee=0;
  for(auto & sp : data) {
    Lgg+=sp.gg.luminosity;
    Lee+=sp.bb.luminosity;
    Ngg+=sp.gg.N;
    Nee+=sp.bb.N;
    Lonline+sp.luminosity.value;
    printf("%5s %10.6f %10.6f %10.6f %10.6f %10ld %10.6f %10.6f %10ld %10.6f \n",
        sp.title.c_str(), sp.energy.value, sp.luminosity.value, 
        sp.gg.luminosity.value*1e-3, sp.gg.luminosity.error*1e-3, sp.gg.N,
        sp.bb.luminosity.value*1e-3, sp.bb.luminosity.error*1e-3, sp.bb.N,
        sp.bb.luminosity.value/sp.gg.luminosity.value
        );
  }
    printf("%5s %10s %10.6f %10.6f %10.6f %10ld %10.6f %10.6f %10ld %10.6f %10.6f\n",
        "total", "", Lonline, Lgg.value*1e-3, Lgg.error*1e-3, Ngg,
                     Lee.value*1e-3, Lee.error*1e-3, Nee,
                     (Lee/Lgg).value, (Lee/Lgg).error
        );
}

void draw_luminosity(Scan_t & data) {
  auto make_graph = [&data] ( int color, int marker, int line, auto F ) {
    auto g = new TGraphErrors;
    for(auto & sp : data) {
      int n = g->GetN();
      ibn::valer<double> x = F(sp);
      g->SetPoint(n, sp.energy.value, (x.value-1)*100);
      g->SetPointError(n, sp.energy.error, x.error*100);
    }
    return g;
  };
  auto g = make_graph(kBlack,1,1, [](auto & sp) { return sp.gg.luminosity/sp.bb.luminosity ; } );
  g->SetMarkerSize(1);
  g->SetMarkerStyle(21);
  g->Draw("ap");
  g->GetXaxis()->SetTitle("W_{cm}, GeV");
  g->GetYaxis()->SetTitle("(L_{#gamma#gamma}/L_{ee} - 1)*100 %");
  /*
  auto g = new TGraphErrors;
  for(auto & sp : data) {
    int n = g->GetN();
    ibn::valer<double> x = sp.gg.luminosity/sp.bb.luminosity;
    g->SetPoint(n, sp.W.value, (x.value-1)*100);
    g->SetPointError(n, sp.W.error, x.error*100);
  }
  g->Draw("al*");
  auto gog = new TGraphErrors;
  for(auto & sp : data) {
    int n = gog->GetN();
    ibn::valer<double> x = sp.gg.luminosity/(sp.L*1e3);
    gog->SetPoint(n, sp.W.value, (x.value-1)*100);
    gog->SetPointError(n, sp.W.error, x.error*100);
  }
  gog->SetLineColor(kRed);
  gog->SetMarkerColor(kRed);
  gog->Draw("l*");
  auto gob = new TGraphErrors;
  for(auto & sp : data) {
    int n = gob->GetN();
    ibn::valer<double> x = sp.bb.luminosity/(sp.L*1e3);
    gob->SetPoint(n, sp.W.value, (x.value-1)*100);
    gob->SetPointError(n, sp.W.error, x.error*100);
  }
  gob->SetLineColor(kBlue);
  gob->SetMarkerColor(kBlue);
  gob->Draw("l*");
  */
}

void data_vs_mc_bhabha(const ScanPoint_t & D, const ScanPoint_t & B) {
  gStyle->Reset("Pub");
  int idx=0;
  auto draw = [&](const char * var, int Nbin, double min, double max, double ymax) {
    auto c = new TCanvas;
    c->SetTitle(var);
    auto h = new TH2F(var,var,Nbin, min,max,Nbin,0, ymax);
    h->Draw();
    auto l = new TLegend(0.8,0.8,1.0,1.0);
    h->GetXaxis()->SetTitle(var);
    D.bb.tree->SetLineColor(kBlue); 
    D.bb.tree->SetLineWidth(3);
    D.bb.tree->Draw(var,BB_SEL.c_str(),"SAME NORM HIS");
    l->AddEntry(D.bb.tree->GetHistogram(),"DATA", "lp");
    B.bb.tree->SetLineWidth(2);
    B.bb.tree->SetLineColor(kRed); 
    B.bb.tree->Draw(var,BB_SEL.c_str(),"SAME NORM HIS");
    l->AddEntry(B.bb.tree->GetHistogram(),"MC", "lp");
    l->Draw();
    c->SaveAs(("data_vs_mc_bhabha"+std::to_string(++idx)+".svg").c_str());
  };
  int N=100;
  draw("E_Eb[0]",N,0.75,1.1,0.07);
  draw("Ep[0]",N, 0.8, 1.1,0.06);
  draw("cos(theta)",N,-0.9,0.9,0.06);
  draw("dtheta*180/TMath::Pi()",N,-3,3,0.035);
  draw("dphi*180/TMath::Pi()",N,-3,3,0.025);
  draw("(acol-TMath::Pi())*180/TMath::Pi()",N,-2.5,0.1,0.06);

  gStyle->Reset("Default");
  auto c1 = new TCanvas;
  c1->SetTitle("Bhabha mc");
  B.bb.tree->Draw("dtheta:dphi", BB_SEL.c_str(),"col");

  auto c2 = new TCanvas;
  c2->SetTitle("Bhabha data");
  D.bb.tree->Draw("dtheta:dphi", BB_SEL.c_str(),"col");
}


#include <chrono>
#include <TMinuitMinimizer.h>
#include <Math/WrappedFunction.h>

template<typename Func>
std::tuple<double,double, double> find_maximum(Func func, double a, double b, double prec=-1) {
  TMinuitMinimizer m(ROOT::Minuit::kSimplex,1);
  m.SetFunction(ROOT::Math::WrappedMultiFunction(
        [&](const double *x) { 
          return -func(x[0]);
        }));
  m.SetStrategy(2);
  m.SetPrecision(prec);
  m.SetLimitedVariable(0,"test",0.5*(a+b),0.2*(b-a),a,b);
  m.Minimize();
  m.PrintResults();
  return {m.X()[0],-m.MinValue(), m.Edm()};
}


void find_optimal_cut(ScanPoint_t & SIG, ScanPoint_t & BG, std::string sel, long N0, std::string param, double param_min, double param_max, double prec=1e-6) {
  auto c = new TCanvas;
  double sigxs = 0.1*1e3; //pb
  double bgxs  = 22*1e3; //pb
  auto g = new TGraphErrors;
  int idx=0;
  double error;
  auto fun = [&](double par) {
    char cut[1024];
    sprintf(cut,"%s && (%s%f)",sel.c_str(),  param.c_str(), par);
    long Nsig = SIG.tt.tree->GetEntries(cut);
    long Nbg  = BG.tt.tree ->GetEntries(cut);

    struct mc_t {
      long  N; //number of selected events
      long  N0; //initial MC number
      double sigma; //the cross section
      double L; //the integrated luminosity
      double eps; //registration efficiency
      double eps_error;
      double Nexp; //expected number of events
      double Nvis; //visible number of events
      double Nvis_error; //visible number of events
      mc_t(long n, long n0, double sig, double lum) {
        N = n;
        N0 = n0;
        sigma = sig;
        L = lum;
        eps = double(N)/N0;
        eps_error = sqrt(double(N))/N0;
        Nexp = fabs(sigma*L);
        Nvis =  eps * Nexp;
        Nvis_error = eps_error * Nexp;
      }
    };

    mc_t sig(Nsig,N0,0.1*1e3, SIG.luminosity);
    mc_t    bg(Nbg,N0,24*1e3, SIG.luminosity);

    double ntot = sig.Nvis + bg.Nvis;
    double x = ntot <=0 ? 0 : sig.Nvis/sqrt(ntot);

    error = std::hypot( sig.Nvis_error*(1./sqrt(ntot) - sig.Nvis*0.5*pow(ntot,-1.5)), bg.Nvis_error*sig.Nvis*0.25*pow(ntot,-1.5));
    return x;
  };
  using clock = std::chrono::system_clock;
  auto t1 = clock::now();
  double xm=std::numeric_limits<double>::lowest();
  double fmax=std::numeric_limits<double>::lowest();
  for(double ptem = 0; ptem < 1.0; ptem+=0.05) {
    double f = fun(ptem);
    if(f>fmax) { xm = ptem; fmax = f;}
    g->SetPoint(idx, ptem, f);
    g->SetPointError(idx++, 0, error);
  }
  std::cout << "time1: " << std::chrono::duration_cast<std::chrono::seconds>(clock::now()-t1).count() << " s" << std::endl;
  g->Draw("ac*");
  t1 = clock::now();
  auto [xmax, ymax, edm] = find_maximum(fun, param_min, param_max);
  std::cout << "time2: " << std::chrono::duration_cast<std::chrono::seconds>(clock::now()-t1).count() << " s" << std::endl;
}

void find_optimal_ptem_cut(ScanPoint_t & SIG, ScanPoint_t & BG, std::string sel, long N0, double prec=1e-6) {
  find_optimal_cut(SIG,BG,sel,N0,"ptem>",0,1);
}

struct PhysicsChannel_t {
  double sigma; //cross section
  ScanPoint_t & point;
};
/*

void find_optimal_cut(ScanPoint_t & SIG, ScanPoint_t & BG, std::string sel, long N0, std::string param, double param_min, double param_max, double prec=1e-6) {
  auto c = new TCanvas;
  double sigxs = 0.1*1e3; //pb
  double bgxs  = 22*1e3; //pb
  auto g = new TGraphErrors;
  int idx=0;
  double error;
  auto fun = [&](double par) {
    char cut[1024];
    sprintf(cut,"%s && (%s%f)",sel.c_str(),  param.c_str(), par);
    long Nsig = SIG.tt->GetEntries(cut);
    long Nbg  = BG.tt ->GetEntries(cut);

    struct mc_t {
      long  N; //number of selected events
      long  N0; //initial MC number
      double sigma; //the cross section
      double L; //the integrated luminosity
      double eps; //registration efficiency
      double eps_error;
      double Nexp; //expected number of events
      double Nvis; //visible number of events
      double Nvis_error; //visible number of events
      mc_t(long n, long n0, double sig, double lum) {
        N = n;
        N0 = n0;
        sigma = sig;
        L = lum;
        eps = double(N)/N0;
        eps_error = sqrt(double(N))/N0;
        Nexp = fabs(sigma*L);
        Nvis =  eps * Nexp;
        Nvis_error = eps_error * Nexp;
      }
    };

    mc_t sig(Nsig,N0,0.1*1e3, SIG.L);
    mc_t    bg(Nbg,N0,24*1e3, SIG.L);

    double ntot = sig.Nvis + bg.Nvis;
    double x = ntot <=0 ? 0 : sig.Nvis/sqrt(ntot);

    error = std::hypot( sig.Nvis_error*(1./sqrt(ntot) - sig.Nvis*0.5*pow(ntot,-1.5)), bg.Nvis_error*sig.Nvis*0.25*pow(ntot,-1.5));
    return x;
  };
  using clock = std::chrono::system_clock;
  auto t1 = clock::now();
  double xm=std::numeric_limits<double>::lowest();
  double fmax=std::numeric_limits<double>::lowest();
  for(double ptem = 0; ptem < 1.0; ptem+=0.05) {
    double f = fun(ptem);
    if(f>fmax) { xm = ptem; fmax = f;}
    g->SetPoint(idx, ptem, f);
    g->SetPointError(idx++, 0, error);
  }
  std::cout << "time1: " << std::chrono::duration_cast<std::chrono::seconds>(clock::now()-t1).count() << " s" << std::endl;
  g->Draw("ac*");
  t1 = clock::now();
  auto [xmax, edm] = find_maximum(fun, param_min, param_max);
  std::cout << "time2: " << std::chrono::duration_cast<std::chrono::seconds>(clock::now()-t1).count() << " s" << std::endl;
}

*/


void draw(ScanPoint_t & D, ScanPoint_t & M, std::string var, std::string cut, int Nbin, double xmin, double xmax) {
  auto c = new TCanvas;
  auto hs = new THStack("hs","");

  c->Divide(2,2);

  c->cd(1);
  D.tt.Draw(var, cut, "");
  auto h1 = D.tt.tree->GetHistogram();
  h1->SetFillColor(kRed);
  hs->Add(h1);

  c->cd(2);
  M.tt.Draw(var, cut,"");
  auto h2 = M.tt.tree->GetHistogram();
  h2->SetFillColor(kBlue);
  hs->Add(h2);

  c->cd(3);
  hs->Draw();
}



TH1F *  fold_histogram(std::vector<ScanPointRef_t> DATA, std::string var, std::string cut, std::string opt, int Nbin, double xmin, double xmax) {
  auto c = new TCanvas;
  c->Divide(2,2);
  auto hf = new TH1F("fold","fold", Nbin, xmin,xmax);
  char buf[4096];
  for(auto [idx, ref] : ibn::indexer(DATA) ) {
    c->cd(idx+1);
    auto & s= ref.get();
    sprintf(buf,"%s>>hstack_%ld(%d,%f,%f)", var.c_str(), idx,Nbin,xmin,xmax);
    s.tt.Draw(buf, cut, "");
    auto h = s.tt.tree->GetHistogram();
    int color = idx+1;
    h->SetFillColor(color);
    h->SetLineColor(color);
    std::cout <<  s.tt.cross_section.value  << " " << s.gg.luminosity.value  << " " << s.tt.N0mc << std::endl;
    //hf->Add(h,s.tt.cross_section.value * s.gg.luminosity.value / s.tt.N0mc);
    hf->Add(h);
  }
  return hf;
};

void draw_stack(std::vector<ScanPointRef_t> DATA, std::string var, std::string cut, int Nbin, double xmin, double xmax) {
  auto c = new TCanvas;
  c->Divide(2,2);
  auto hs = new THStack("hs","");

  char buf[4096];

  int idx=0;
  for(auto  ref : DATA) {
    c->cd(idx+1);
    auto & s= ref.get();
    sprintf(buf,"%s>>hstack_%d(%d,%f,%f)", var.c_str(), idx,Nbin,xmin,xmax);
    s.tt.Draw(buf, cut, "");
    auto h = s.tt.tree->GetHistogram();
    int color = idx+1;
    h->SetFillColor(color);
    h->SetLineColor(color);
    hs->Add(h);
    idx++;
  }
  c->cd(idx+1);
  hs->Draw();
}


void draw(ScanPoint_t & DATA, ScanPoint_t & MC, std::vector<ScanPointRef_t> BGs, Selection_t & SEL, int i, std::string var, int Nbin, double xmin, double xmax) {
  std::string cut = SEL.common_cut() + " && " +  SEL[i].cut;
  //std::string cut = "";
  auto hs = new THStack("hs","");
  auto h = fold_histogram(BGs, var, cut, "HIS", Nbin,xmin,xmax);
  h->SetFillColor(kRed);
  hs->Add(h);
  auto c = new TCanvas;
  c->Divide(2,2);
  c->cd(1);
  h->Draw();
  h->SetTitle("Background");
  char buf[4096];
  int idx=0;

  c->cd(2);
  sprintf(buf,"%s>>hsig_%d(%d,%f,%f)", var.c_str(), idx,Nbin,xmin,xmax);
  MC.tt.tree->Draw(buf, cut.c_str(),"HIS");
  auto hmc = MC.tt.tree->GetHistogram();
  hmc->Scale(MC.tt.cross_section.value*MC.gg.luminosity.value/MC.tt.N0mc);
  hmc->SetFillColor(kBlue);

  c->cd(3);
  hs->Add(hmc);

  sprintf(buf,"%s>>data_%d(%d,%f,%f)", var.c_str(), idx,Nbin,xmin,xmax);
  DATA.tt.tree->SetLineColor(kBlack);
  DATA.tt.tree->SetMarkerColor(kBlack);
  DATA.tt.tree->SetLineWidth(3);
  DATA.tt.tree->Draw(buf, cut.c_str(), "E");
  auto dh=DATA.tt.tree->GetHistogram();
  //dh->Draw("same");

  hs->Draw("same");
}


TH1 * select_histogram(const ScanPoint_t & sp, std::string var, std::string cut, int Nbin, double xmin, double xmax, bool is_scale=false) {
  char buf[4096];
  sprintf(buf,"%s>>h%d(%d,%f,%f)", var.c_str(), ++HISTO_INDEX, Nbin,xmin,xmax);
  //std::cout << buf << std::endl;
  sp.tt.tree->Draw(buf, cut.c_str(),"goff");
  auto h = sp.tt.tree->GetHistogram(); if(is_scale) h->Scale(sp.tt.cross_section.value*sp.gg.luminosity.value/sp.tt.N0mc);
  return h;
}

TH1 * select_histogram(const Scan_t & SCAN, std::string var, std::string cut, int Nbin, double xmin, double xmax, bool is_scale=false) {
  std::vector<TH1*> Hs;
  for(auto & sp : SCAN) {
    Hs.push_back( select_histogram(sp,var,cut,Nbin,xmin,xmax, is_scale));
  }
  return fold(Hs, Nbin,xmin,xmax);
};

TH1 * select_histogram(const std::vector<ScanRef_t> Ss, std::string var, std::string cut, int Nbin, double xmin, double xmax, bool is_scale=false) {
  std::vector<TH1*> Hs;
  for(auto & scan : Ss) {
    Hs.push_back( select_histogram(scan,var,cut,Nbin,xmin,xmax, is_scale));
  }
  return fold(Hs,Nbin,xmin,xmax);
};

std::tuple<double,double> my_chi2_histogram_test ( TH1 * h1, TH1 *h2) {
  int Nbin1{h1->GetNbinsX()};
  int Nbin2{h2->GetNbinsX()};
  if(Nbin1!=Nbin2)  {
    std::cerr << "ERROR: Different number of bins: " << Nbin1 << ", " << Nbin2 << " for histogram " << h1->GetName() << " and " << h2->GetName() << std::endl;
    return {-1,-1};
  }
  double chi2{0};
  int Ndf{0};
  for(int i=0;i<Nbin1; ++i) {
    ibn::valer<double> y1{h1->GetBinContent(i), h1->GetBinError(i)};
    ibn::valer<double> y2{h2->GetBinContent(i), h2->GetBinError(i)};
    //if(std::max(y1.error,y2.error) < 2*std::numeric_limits<double>::min()) {
    //  y1.error = 1;
    //};
    static constexpr double epsilon = std::numeric_limits<double>::min();
    if( y1.error< 2*epsilon ) y1.error = 1.0;
    if( y2.error< 2*epsilon ) y2.error = 1.0;
    if( y1.value > 2*epsilon && y2.value > 2*epsilon) {
      double chi2_tmp = pow(y1.value - y2.value,2.0)/(y1.error*y1.error + y2.error*y2.error);
      chi2 += chi2_tmp;
      //std::cout << i << std::setw(10) << y1.value << std::setw(10) << y2.value << std::setw(10) << y1.value - y2.value <<std::setw(10)  << y1.error << std::setw(10) << y2.error << std::setw(10)  << std::hypot(y1.error, y2.error) << " chi2_tmp = " << std::setw(10) << chi2_tmp << " chi2 = " << std::setw(10) << chi2 << std::endl;
      ++Ndf;
    }
  }
  //std::cout << chi2/(2*Ndf) << std::endl;
  return {chi2, TMath::Prob(chi2,2*Ndf)};
}

TCanvas * draw(const Selection_t & SEL, const Scan_t & DATA, const Scan_t & SIGNAL, const std::vector<ScanRef_t> BGs, std::string var, int Nbin, double xmin, double xmax, std::string extracut="") {
  auto  c = new TCanvas;
  gStyle->SetPalette(kOcean);
  auto hs = new THStack(("hstack"+std::to_string(++HISTO_INDEX)).c_str(),var.c_str());
  std::vector<TH1*> Hsig; //signal histogram for different selection channels
  std::vector<TH1*> Hdata;
  std::vector< std::vector<TH1*> > Hbgs(BGs.size()); 

  for(auto & sel : SEL ) { //loop over selection channels
    std::string cut = SEL.common_cut() + " && " +  sel.cut + (extracut == "" ? "" :  (" && " + extracut) );
    //std::string cut = remove_some_cuts(var,SEL.common_cut + " && " +  sel.cut + (extracut == "" ? "" :  (" && " + extracut) ));
    std::cout << cut << std::endl;
    const bool SCALE = true;
    const bool NOSCALE = false;

    auto select = [&var,&cut,&Nbin,&xmin, &xmax] (auto &  sample , bool is_scale=false) {
      return select_histogram(sample,var,cut, Nbin,xmin,xmax, is_scale);
    };

    //select data
    auto hdata = select(DATA);
    Hdata.push_back(hdata);

    //select signal monte carlo
    auto hsig = select(SIGNAL,SCALE);
    Hsig.push_back(hsig);

    //select background monte carlo
    for(size_t i=0;i!=BGs.size();++i) {
      auto hbg = select(BGs[i],SCALE);
      Hbgs[i].push_back(hbg);
    }
  }

  auto lfold = [Nbin,xmin,xmax] (auto & hs) {
    return fold(hs,Nbin,xmin,xmax);
  };

  TH1 * hsig = lfold(Hsig);
  TH1 * hdata = lfold(Hdata);

  std::vector<TH1*> hbgs(BGs.size());
  for(size_t i=0;i!=BGs.size();++i) {
    hbgs[i] = lfold(Hbgs[i]);
    std::cout << "bg " << i << "  " << hbgs[i]->GetEntries() << "  " << hbgs[i]->Integral() << std::endl;
  };
  TH1 * hbg = lfold(hbgs);
  double Nbg_mc = hbg->Integral();
  static constexpr double Nbg_exp = 45.9828;
  double scale_bg = Nbg_exp/Nbg_mc;
  hbg->Scale(scale_bg);
  for(auto & h : hbgs) h->Scale(scale_bg);

  auto set_color = [](auto & o , int color) {
    o->SetFillColor(color);
    o->SetLineColor(color);
    o->SetMarkerColor(color);
  };

  set_color(hsig,kRed);
  set_color(hdata,kBlue);
  set_color(hbg, kBlack);

  //int col=kBlack;
  int col=0;
  std::vector palette{1,8,6,11,40,42,46,29};
  for(auto & h :hbgs) {
    //while(col==2 || col==4 || col == 3 || col == 10) ++col;
    std::cout << "h " << h->GetName() << " color = " << palette[col] << " kRed = " << kRed << "  kBlue = " << kBlue << " kBlack= " << kBlack <<  std::endl;
    set_color(h,palette[col]);
    ++col;
  }
  for(auto & h :hbgs) {
    hs->Add(h,"HIST");
  }
  //hs->Add(hbg,"HIST");

  hsig->SetFillColor(kWhite);
  hs->Add(hsig,"HIST");

  hdata->Draw("E");
  hdata->SetTitle(var.c_str());
  hdata->GetXaxis()->SetTitle(var.c_str());
  hs->Draw("same");
  hdata->SetLineWidth(3);
  hdata->Draw("same E");

  auto lbl = new TLegend(0.6,0.8,1.0,1.0);
  //std::vector<ScanRef_t> BGs ={HADR, BB, UU, GG, GALUGA["ee"], GALUGA["uu"],GALUGA["pipi"], GALUGA["KK"]};
  {
    lbl->AddEntry(hdata,"data","pl");
    lbl->AddEntry(hsig,"MC signal+bg", "plf");
    std::vector bg_leg_name = {"hadrons","Bhabha", "#mu#mu","#gamma#gamma","EEee","EE#mu#mu","EE#pi#pi","EEKK"};
    for(size_t i = 0; i!=hbgs.size(); ++i) {
      lbl->AddEntry(hbgs[i],bg_leg_name[i],"plf");
    }
  }

  lbl->Draw();

  //auto Ndata = hdata->GetEntries();
  auto get_integral = [Nbin](TH1 * h) -> ibn::valer<double> {
    ibn::valer<double> N;
    N.value = h->IntegralAndError(0,Nbin,N.error);
    return N;
  };
  ibn::valer<double> Ndata  = get_integral(hdata);
  ibn::valer<double> Nsig   = get_integral(hsig);
  ibn::valer<double> Nbg    = get_integral(hbg);
  ibn::valer<double> Nsum   = Nsig + Nbg;

  //std::cout << "bg: " << Nbg << std::endl;
  auto hsig_bg = fold({hsig, hbg},Nbin,xmin,xmax);
  double prob_kolm_sig_bg = hdata->KolmogorovTest(hsig_bg);
  //std::cout << "Kolmogorov test (sig+bg): " << prob_kolm_sig_bg << std::endl;
  double prob_kolm_sig = hdata->KolmogorovTest(hsig);
  //std::cout << "Kolmogorov test (sig): " << prob_kolm_sig << std::endl;
  double prob_chi2_sig_bg = hdata->Chi2Test(hsig_bg,"UW");
  //std::cout << "Chi2 test (sig+bg): " << prob_chi2_sig_bg << std::endl;
  double prob_chi2_sig = hdata->Chi2Test(hsig,"UW");
  //std::cout << "Chi2 test (sig): " << prob_chi2_sig << std::endl;
  //
  std::cout << "My chi2 test" << std::endl;
  auto [chi2_my_sig_bg, prob_my_sig_bg] = my_chi2_histogram_test(hsig_bg,hdata);
  auto [chi2_my_sig, prob_my_sig] = my_chi2_histogram_test(hsig,hdata);
  std::cout << "chi2_my_sig_bg = " << chi2_my_sig_bg << std::endl;
  std::cout << "chi2_my_sig = " << chi2_my_sig << std::endl;

  std::cout << "====================  Number of events ===================  " << std::endl;
  std::cout << std::setw(20) <<  "sample"  << std::setw(15) << "N" << std::endl;
  std::cout << "----------------------------------------------------------  " << std::endl;
  std::cout << std::setw(20) <<  "signal"      << std::setw(15) << Nsig     <<  " ± " << Nsig.error  << std::endl;
  std::cout << std::setw(20) <<  "background"  << std::setw(15) << Nbg      <<  " ± " << Nbg.error   << std::endl;
  std::cout << std::setw(20) <<  "sig + bg"    << std::setw(15) << Nsum     <<  " ± " << Nsum.error  << std::endl;
  std::cout << std::setw(20) <<  "data"        << std::setw(15) << Ndata    <<  " ± " << Ndata.error << std::endl;
  std::cout << "====================  Probability Test ===================  " << std::endl;
  std::cout << std::setw(20) << " sample "             << std::setw(15) << " Kolmogorov "     << std::setw(15)   <<   " Chi2 "           <<  std::setw(15) << "My chi2" << std::endl;
  std::cout << "----------------------------------------------------------  " << std::endl;
  std::cout << std::setw(20) << "signal + background"  << std::setw(15) << prob_kolm_sig_bg   << std::setw(15)   <<   prob_chi2_sig_bg   <<  std::setw(15) << prob_my_sig_bg << std::endl;
  std::cout << std::setw(20) << "signal"               << std::setw(15) << prob_kolm_sig      << std::setw(15)   <<   prob_chi2_sig      <<  std::setw(15) << prob_my_sig    << std::endl;
  std::cout << "----------------------------------------------------------  " << std::endl;
  return c;
}

TCanvas * draw(const Selection_t & SEL, const Scan_t & DATA, const Simulation_t & MC, std::string var, int Nbin, double xmin, double xmax, std::string extracut="") {
  auto  c = new TCanvas;
  gStyle->SetPalette(kOcean);
  auto hs = new THStack(("hstack"+std::to_string(++HISTO_INDEX)).c_str(),var.c_str());
  std::vector<TH1*> Hsig; //signal histogram for different selection channels
  std::vector<TH1*> Hdata;
  auto & SIGNAL = MC.signal;
  auto & BGs = MC.bgs;
  std::vector< std::vector<TH1*> > Hbgs(BGs.size()); 

  auto sel = [&]( const Selection_t & sel ) {
      std::string cut = sel.cut + (extracut == "" ? "" :  (" && " + extracut) );
      std::cout << cut << std::endl;
      const bool SCALE = true;
      const bool NOSCALE = false;

      auto select = [&var,&cut,&Nbin,&xmin, &xmax] (auto &  sample , bool is_scale=false) {
        return select_histogram(sample,var,cut, Nbin,xmin,xmax, is_scale);
      };

      //select data
      auto hdata = select(DATA);
      Hdata.push_back(hdata);

      //select signal monte carlo
      auto hsig = select(SIGNAL,SCALE);
      Hsig.push_back(hsig);

      //select background monte carlo
      for(size_t i=0;i!=BGs.size();++i) {
        auto hbg = select(BGs[i],SCALE);
        Hbgs[i].push_back(hbg);
      }
  };
  if(SEL.empty()) sel(SEL);
  else for(auto & s : SEL ) sel(s);

  auto lfold = [Nbin,xmin,xmax] (auto & hs) {
    return fold(hs,Nbin,xmin,xmax);
  };

  TH1 * hsig = lfold(Hsig);
  TH1 * hdata = lfold(Hdata);

  std::vector<TH1*> hbgs(BGs.size());
  for(size_t i=0;i!=BGs.size();++i) {
    hbgs[i] = lfold(Hbgs[i]);
    std::cout << "bg " << i << "  " << hbgs[i]->GetEntries() << "  " << hbgs[i]->Integral() << std::endl;
  };
  TH1 * hbg = lfold(hbgs);
  double Nbg_mc = hbg->Integral();
  static constexpr double Nbg_exp = 45.9828;
  double scale_bg = Nbg_exp/Nbg_mc;
  if(!SEL.empty()){
    hbg->Scale(scale_bg);
    for(auto & h : hbgs) h->Scale(scale_bg);
  }

  auto set_color = [](auto & o , int color) {
    o->SetFillColor(color);
    o->SetLineColor(color);
    o->SetMarkerColor(color);
  };

  set_color(hsig,kRed);
  set_color(hdata,kBlue);
  set_color(hbg, kBlack);

  //int col=kBlack;
  int col=0;
  std::vector palette{1,8,6,11,40,42,46,29};
  for(auto & h :hbgs) {
    //while(col==2 || col==4 || col == 3 || col == 10) ++col;
    std::cout << "h " << h->GetName() << " color = " << palette[col] << " kRed = " << kRed << "  kBlue = " << kBlue << " kBlack= " << kBlack <<  std::endl;
    set_color(h,palette[col]);
    ++col;
  }
  for(auto & h :hbgs) {
    hs->Add(h,"HIST");
  }
  //hs->Add(hbg,"HIST");

  hsig->SetFillColor(kWhite);
  hs->Add(hsig,"HIST");

  hdata->Draw("E");
  hdata->SetTitle(var.c_str());
  hdata->GetXaxis()->SetTitle(var.c_str());
  hs->Draw("same");
  hdata->SetLineWidth(3);
  hdata->Draw("same E");

  auto lbl = new TLegend(0.6,0.8,1.0,1.0);
  //std::vector<ScanRef_t> BGs ={HADR, BB, UU, GG, GALUGA["ee"], GALUGA["uu"],GALUGA["pipi"], GALUGA["KK"]};
  {
    lbl->AddEntry(hdata,"data","pl");
    lbl->AddEntry(hsig,"MC signal+bg", "plf");
    std::vector bg_leg_name = {"hadrons","Bhabha", "#mu#mu","#gamma#gamma","EEee","EE#mu#mu","EE#pi#pi","EEKK"};
    for(size_t i = 0; i!=hbgs.size(); ++i) {
      lbl->AddEntry(hbgs[i],bg_leg_name[i],"plf");
    }
  }

  lbl->Draw();

  //auto Ndata = hdata->GetEntries();
  auto get_integral = [Nbin](TH1 * h) -> ibn::valer<double> {
    ibn::valer<double> N;
    N.value = h->IntegralAndError(0,Nbin,N.error);
    return N;
  };
  ibn::valer<double> Ndata  = get_integral(hdata);
  ibn::valer<double> Nsig   = get_integral(hsig);
  ibn::valer<double> Nbg    = get_integral(hbg);
  ibn::valer<double> Nsum   = Nsig + Nbg;

  //std::cout << "bg: " << Nbg << std::endl;
  auto hsig_bg = fold({hsig, hbg},Nbin,xmin,xmax);
  double prob_kolm_sig_bg = hdata->KolmogorovTest(hsig_bg);
  //std::cout << "Kolmogorov test (sig+bg): " << prob_kolm_sig_bg << std::endl;
  double prob_kolm_sig = hdata->KolmogorovTest(hsig);
  //std::cout << "Kolmogorov test (sig): " << prob_kolm_sig << std::endl;
  double prob_chi2_sig_bg = hdata->Chi2Test(hsig_bg,"UW");
  //std::cout << "Chi2 test (sig+bg): " << prob_chi2_sig_bg << std::endl;
  double prob_chi2_sig = hdata->Chi2Test(hsig,"UW");
  //std::cout << "Chi2 test (sig): " << prob_chi2_sig << std::endl;
  //
  std::cout << "My chi2 test" << std::endl;
  auto [chi2_my_sig_bg, prob_my_sig_bg] = my_chi2_histogram_test(hsig_bg,hdata);
  auto [chi2_my_sig, prob_my_sig] = my_chi2_histogram_test(hsig,hdata);
  std::cout << "chi2_my_sig_bg = " << chi2_my_sig_bg << std::endl;
  std::cout << "chi2_my_sig = " << chi2_my_sig << std::endl;

  std::cout << "====================  Number of events ===================  " << std::endl;
  std::cout << std::setw(20) <<  "sample"  << std::setw(15) << "N" << std::endl;
  std::cout << "----------------------------------------------------------  " << std::endl;
  std::cout << std::setw(20) <<  "signal"      << std::setw(15) << Nsig     <<  " ± " << Nsig.error  << std::endl;
  std::cout << std::setw(20) <<  "background"  << std::setw(15) << Nbg      <<  " ± " << Nbg.error   << std::endl;
  std::cout << std::setw(20) <<  "sig + bg"    << std::setw(15) << Nsum     <<  " ± " << Nsum.error  << std::endl;
  std::cout << std::setw(20) <<  "data"        << std::setw(15) << Ndata    <<  " ± " << Ndata.error << std::endl;
  std::cout << "====================  Probability Test ===================  " << std::endl;
  std::cout << std::setw(20) << " sample "             << std::setw(15) << " Kolmogorov "     << std::setw(15)   <<   " Chi2 "           <<  std::setw(15) << "My chi2" << std::endl;
  std::cout << "----------------------------------------------------------  " << std::endl;
  std::cout << std::setw(20) << "signal + background"  << std::setw(15) << prob_kolm_sig_bg   << std::setw(15)   <<   prob_chi2_sig_bg   <<  std::setw(15) << prob_my_sig_bg << std::endl;
  std::cout << std::setw(20) << "signal"               << std::setw(15) << prob_kolm_sig      << std::setw(15)   <<   prob_chi2_sig      <<  std::setw(15) << prob_my_sig    << std::endl;
  std::cout << "----------------------------------------------------------  " << std::endl;
  return c;
}

//TCanvas * unconstrain_par_and_draw(const ChannelSelection_t & sel, const  Scan_t & DATA, const Simulation_t & SIM, std::string vars, int Nbin, double xmin, double xmax, std::string extracut="", double bg_scale=-1) {
//  //auto  c = new TCanvas(sel.root_title().c_str(), sel.root_title().c_str());
//  auto pars = split(vars,"&&");
//  std::string var(pars[0]);
//  //cout << var << endl;
//  auto c = get_new_tailed_canvas(sel.root_title2() + ":" + var);
//  gStyle->SetPalette(kOcean);
//  gStyle->SetOptStat(kFALSE);
//  auto hs = new THStack(("hstack"+std::to_string(++HISTO_INDEX)).c_str(),var.c_str());
//  std::vector<TH1*> Hsig; //signal histogram for different selection channels
//  std::vector<TH1*> Hdata;
//  std::vector< std::vector<TH1*> > Hbgs(SIM.bgs.size()); 
//
//  {
//    //std::string cut = SEL.common_cut + " && " +  sel.cut + (extracut == "" ? "" :  (" && " + extracut) );
//    if(extracut!="") extracut = "&&" + extracut;
//    std::string cut = remove_some_cuts(vars,sel.cut+extracut);
//    std::cout << cut << std::endl;
//    const bool SCALE = true;
//    const bool NOSCALE = false;
//
//    auto select = [&var,&cut,&Nbin,&xmin, &xmax] (auto &  sample , bool is_scale=false) {
//      return select_histogram(sample,var,cut, Nbin,xmin,xmax, is_scale);
//    };
//
//    //select data
//    auto hdata = select(DATA);
//    Hdata.push_back(hdata);
//
//    //select signal monte carlo
//    auto hsig = select(SIM.signal,SCALE);
//    Hsig.push_back(hsig);
//
//    //select background monte carlo
//    for(size_t i=0;i!=SIM.bgs.size();++i) {
//      auto hbg = select(SIM.bgs[i],SCALE);
//      Hbgs[i].push_back(hbg);
//    }
//  }
//
//  auto lfold = [Nbin,xmin,xmax] (auto & hs) {
//    return fold(hs,Nbin,xmin,xmax);
//  };
//
//  TH1 * hsig = lfold(Hsig);
//  TH1 * hdata = lfold(Hdata);
//
//  std::vector<TH1*> hbgs(SIM.bgs.size());
//  for(size_t i=0;i!=SIM.bgs.size();++i) {
//    hbgs[i] = lfold(Hbgs[i]);
//    std::cout << "bg " << i << "  " << hbgs[i]->GetEntries() << "  " << hbgs[i]->Integral() << std::endl;
//  };
//  TH1 * hbg = lfold(hbgs);
//  /*
//     double Nbg_mc = hbg->Integral();
//     static constexpr double Nbg_exp = 45.9828;
//     double scale_bg = Nbg_exp/Nbg_mc;
//     */
//  double scale_bg = bg_scale >= 0 ? bg_scale  :  1.0;
//  //if(bg_scale>0) scale_bg=bg_scale;
//
//  hbg->Scale(scale_bg);
//  for(auto & h : hbgs) h->Scale(scale_bg);
//
//  auto set_color = [](auto & o , int color) {
//    o->SetFillColor(color);
//    o->SetLineColor(color);
//    o->SetMarkerColor(color);
//  };
//
//  set_color(hsig,kRed);
//  set_color(hdata,kBlue);
//  set_color(hbg, kBlack);
//
//  //int col=kBlack;
//  int col=0;
//  std::vector palette{1,8,6,11,40,42,46,29};
//  for(auto & h :hbgs) {
//    //while(col==2 || col==4 || col == 3 || col == 10) ++col;
//    std::cout << "h " << h->GetName() << " color = " << palette[col] << " kRed = " << kRed << "  kBlue = " << kBlue << " kBlack= " << kBlack <<  std::endl;
//    set_color(h,palette[col]);
//    ++col;
//  }
//  for(auto & h :hbgs) {
//    hs->Add(h,"HIST");
//  }
//  //hs->Add(hbg,"HIST");
//
//  hsig->SetFillColor(kWhite);
//  hs->Add(hsig,"HIST");
//
//  hdata->Draw("E");
//  hdata->SetTitle((sel.root_title() + " : " + var).c_str());
//  hdata->GetXaxis()->SetTitle(var.c_str());
//  hs->Draw("same");
//  hdata->SetLineWidth(3);
//  hdata->Draw("same E");
//
//  auto lbl = new TLegend(0.8,0.6,1.0,1.0);
//  //std::vector<ScanRef_t> BGs ={HADR, BB, UU, GG, GALUGA["ee"], GALUGA["uu"],GALUGA["pipi"], GALUGA["KK"]};
//  {
//    lbl->AddEntry(hdata,"data","pl");
//    lbl->AddEntry(hsig,"MC signal+bg", "plf");
//    std::vector bg_leg_name = {"hadrons","Bhabha", "#mu#mu","#gamma#gamma","EEee","EE#mu#mu","EE#pi#pi","EEKK"};
//    for(size_t i = 0; i!=hbgs.size(); ++i) {
//      lbl->AddEntry(hbgs[i],bg_leg_name[i],"plf");
//    }
//  }
//
//  lbl->Draw();
//
//  //auto Ndata = hdata->GetEntries();
//  auto get_integral = [Nbin](TH1 * h) -> ibn::valer<double> {
//    ibn::valer<double> N;
//    N.value = h->IntegralAndError(0,Nbin,N.error);
//    return N;
//  };
//  ibn::valer<double> Ndata  = get_integral(hdata);
//  ibn::valer<double> Nsig   = get_integral(hsig);
//  ibn::valer<double> Nbg    = get_integral(hbg);
//  ibn::valer<double> Nsum   = Nsig + Nbg;
//
//  //std::cout << "bg: " << Nbg << std::endl;
//  auto hsig_bg = fold({hsig, hbg},Nbin,xmin,xmax);
//  double prob_kolm_sig_bg = hdata->KolmogorovTest(hsig_bg);
//  //std::cout << "Kolmogorov test (sig+bg): " << prob_kolm_sig_bg << std::endl;
//  double prob_kolm_sig = hdata->KolmogorovTest(hsig);
//  //std::cout << "Kolmogorov test (sig): " << prob_kolm_sig << std::endl;
//  double prob_chi2_sig_bg = hdata->Chi2Test(hsig_bg,"UW");
//  //std::cout << "Chi2 test (sig+bg): " << prob_chi2_sig_bg << std::endl;
//  double prob_chi2_sig = hdata->Chi2Test(hsig,"UW");
//  //std::cout << "Chi2 test (sig): " << prob_chi2_sig << std::endl;
//  //
//  std::cout << "My chi2 test" << std::endl;
//  auto [chi2_my_sig_bg, prob_my_sig_bg] = my_chi2_histogram_test(hsig_bg,hdata);
//  auto [chi2_my_sig, prob_my_sig] = my_chi2_histogram_test(hsig,hdata);
//  std::cout << "chi2_my_sig_bg = " << chi2_my_sig_bg << std::endl;
//  std::cout << "chi2_my_sig = " << chi2_my_sig << std::endl;
//
//  std::cout << "====================  Number of events ===================  " << std::endl;
//  std::cout << std::setw(20) <<  "sample"  << std::setw(15) << "N" << std::endl;
//  std::cout << "----------------------------------------------------------  " << std::endl;
//  std::cout << std::setw(20) <<  "signal"      << std::setw(15) << Nsig     <<  " ± " << Nsig.error  << std::endl;
//  std::cout << std::setw(20) <<  "background"  << std::setw(15) << Nbg      <<  " ± " << Nbg.error   << std::endl;
//  std::cout << std::setw(20) <<  "sig + bg"    << std::setw(15) << Nsum     <<  " ± " << Nsum.error  << std::endl;
//  std::cout << std::setw(20) <<  "data"        << std::setw(15) << Ndata    <<  " ± " << Ndata.error << std::endl;
//  std::cout << "====================  Probability Test ===================  " << std::endl;
//  std::cout << std::setw(20) << " sample "             << std::setw(15) << " Kolmogorov "     << std::setw(15)   <<   " Chi2 "           <<  std::setw(15) << "My chi2" << std::endl;
//  std::cout << "----------------------------------------------------------  " << std::endl;
//  std::cout << std::setw(20) << "signal + background"  << std::setw(15) << prob_kolm_sig_bg   << std::setw(15)   <<   prob_chi2_sig_bg   <<  std::setw(15) << prob_my_sig_bg << std::endl;
//  std::cout << std::setw(20) << "signal"               << std::setw(15) << prob_kolm_sig      << std::setw(15)   <<   prob_chi2_sig      <<  std::setw(15) << prob_my_sig    << std::endl;
//  std::cout << "----------------------------------------------------------  " << std::endl;
//  gSystem->ProcessEvents();
//  return c;
//}
//
//void unconstrain_par_and_draw_per_channel(const Selection_t & SEL, const  Scan_t & DATA, const Simulation_t & SIM, std::string var, int Nbin, double xmin, double xmax, std::string extracut="", double bg_scale=-1) {
//  for(const auto & sel : SEL) {
//    unconstrain_par_and_draw(sel,DATA,SIM,var,Nbin,xmin,xmax,extracut, bg_scale);
//  }
//}

TCanvas * unconstrain_par_and_draw(const Selection_t & SEL, const  Scan_t & DATA, const Simulation_t & SIM, const std::string vars, int Nbin, double xmin, double xmax, std::string extracut="", double bg_scale=-1) {
  gStyle->Reset("Pub");
  //auto  c = new TCanvas;
  const std::string var= std::string(split(vars,"&&")[0]);
  auto  c = get_new_tailed_canvas(var);
  gStyle->SetPalette(kOcean);
  auto hs = new THStack(("hstack"+std::to_string(++HISTO_INDEX)).c_str(),var.c_str());
  std::vector<TH1*> Hsig; //signal histogram for different selection channels
  std::vector<TH1*> Hdata;
  std::vector< std::vector<TH1*> > Hbgs(SIM.bgs.size()); 

  auto sel = [&](const Selection_t & sel) {
    //std::string cut = SEL.common_cut + " && " +  sel.cut + (extracut == "" ? "" :  (" && " + extracut) );
    //std::string cut = remove_some_cuts(var,SEL.common_cut() + " && " +  sel.cut + (extracut == "" ? "" :  (" && " + extracut) ));
    std::string cut = remove_some_cuts(vars,sel.cut) + (extracut == "" ? "" :  (" && " + extracut));
    std::cout << cut << std::endl;
    const bool SCALE = true;
    const bool NOSCALE = false;

    auto select = [&var,&cut,&Nbin,&xmin, &xmax] (auto &  sample , bool is_scale=false) {
      return select_histogram(sample,var,cut, Nbin,xmin,xmax, is_scale);
    };

    //select data
    auto hdata = select(DATA);
    Hdata.push_back(hdata);

    //select signal monte carlo
    auto hsig = select(SIM.signal,SCALE);
    Hsig.push_back(hsig);

    //select background monte carlo
    for(size_t i=0;i!=SIM.bgs.size();++i) {
      auto hbg = select(SIM.bgs[i],SCALE);
      Hbgs[i].push_back(hbg);
    }
  };

  if(SEL.empty()) sel(SEL);
  else {
    for(auto & s : SEL ) { sel(s); }
  }

  auto lfold = [Nbin,xmin,xmax] (auto & hs) {
    return fold(hs,Nbin,xmin,xmax);
  };

  TH1 * hsig = lfold(Hsig);
  TH1 * hdata = lfold(Hdata);

  std::vector<TH1*> hbgs(SIM.bgs.size());
  for(size_t i=0;i!=SIM.bgs.size();++i) {
    hbgs[i] = lfold(Hbgs[i]);
    std::cout << "bg " << i << "  " << hbgs[i]->GetEntries() << "  " << hbgs[i]->Integral() << std::endl;
  };
  TH1 * hbg = lfold(hbgs);
  /*
  double Nbg_mc = hbg->Integral();
  static constexpr double Nbg_exp = 45.9828;
  double scale_bg = Nbg_exp/Nbg_mc;
  */
  double scale_bg = bg_scale >= 0 ? bg_scale  :  1.0;
  //if(bg_scale>0) scale_bg=bg_scale;

  hbg->Scale(scale_bg);
  for(auto & h : hbgs) h->Scale(scale_bg);

  auto set_color = [](auto & o , int color) {
    o->SetFillColor(color);
    o->SetLineColor(color);
    o->SetMarkerColor(color);
  };

  set_color(hsig,kRed);
  set_color(hdata,kBlue);
  set_color(hbg, kBlack);

  //int col=kBlack;
  int col=0;
  std::vector palette{1,8,6,11,40,42,46,29};
  for(auto & h :hbgs) {
    //while(col==2 || col==4 || col == 3 || col == 10) ++col;
    std::cout << "h " << h->GetName() << " color = " << palette[col] << " kRed = " << kRed << "  kBlue = " << kBlue << " kBlack= " << kBlack <<  std::endl;
    set_color(h,palette[col]);
    ++col;
  }
  for(auto & h :hbgs) {
    hs->Add(h,"HIST");
  }
  //hs->Add(hbg,"HIST");

  hsig->SetFillColor(kWhite);
  hs->Add(hsig,"HIST");

  hdata->Draw("E");
  hdata->SetTitle(var.c_str());
  hdata->GetXaxis()->SetTitle(var.c_str());
  hs->Draw("same");
  hdata->SetLineWidth(3);
  hdata->Draw("same E");

  auto lbl = new TLegend(0.8,0.5,1.0,1.0);
  //std::vector<ScanRef_t> BGs ={HADR, BB, UU, GG, GALUGA["ee"], GALUGA["uu"],GALUGA["pipi"], GALUGA["KK"]};
  {
    lbl->AddEntry(hdata,"data","pl");
    lbl->AddEntry(hsig,"MC signal+bg", "plf");
    std::vector bg_leg_name = {"hadrons","Bhabha", "#mu#mu","#gamma#gamma","EEee","EE#mu#mu","EE#pi#pi","EEKK"};
    for(size_t i = 0; i!=hbgs.size(); ++i) {
      lbl->AddEntry(hbgs[i],bg_leg_name[i],"plf");
    }
  }

  lbl->Draw();

  //auto Ndata = hdata->GetEntries();
  auto get_integral = [Nbin](TH1 * h) -> ibn::valer<double> {
    ibn::valer<double> N;
    N.value = h->IntegralAndError(0,Nbin,N.error);
    return N;
  };
  ibn::valer<double> Ndata  = get_integral(hdata);
  ibn::valer<double> Nsig   = get_integral(hsig);
  ibn::valer<double> Nbg    = get_integral(hbg);
  ibn::valer<double> Nsum   = Nsig + Nbg;

  //std::cout << "bg: " << Nbg << std::endl;
  auto hsig_bg = fold({hsig, hbg},Nbin,xmin,xmax);
  double prob_kolm_sig_bg = hdata->KolmogorovTest(hsig_bg);
  //std::cout << "Kolmogorov test (sig+bg): " << prob_kolm_sig_bg << std::endl;
  double prob_kolm_sig = hdata->KolmogorovTest(hsig);
  //std::cout << "Kolmogorov test (sig): " << prob_kolm_sig << std::endl;
  double prob_chi2_sig_bg = hdata->Chi2Test(hsig_bg,"UW");
  //std::cout << "Chi2 test (sig+bg): " << prob_chi2_sig_bg << std::endl;
  double prob_chi2_sig = hdata->Chi2Test(hsig,"UW");
  //std::cout << "Chi2 test (sig): " << prob_chi2_sig << std::endl;
  //
  std::cout << "My chi2 test" << std::endl;
  auto [chi2_my_sig_bg, prob_my_sig_bg] = my_chi2_histogram_test(hsig_bg,hdata);
  auto [chi2_my_sig, prob_my_sig] = my_chi2_histogram_test(hsig,hdata);
  std::cout << "chi2_my_sig_bg = " << chi2_my_sig_bg << std::endl;
  std::cout << "chi2_my_sig = " << chi2_my_sig << std::endl;

  std::cout << "====================  Number of events ===================  " << std::endl;
  std::cout << std::setw(20) <<  "sample"  << std::setw(15) << "N" << std::endl;
  std::cout << "----------------------------------------------------------  " << std::endl;
  std::cout << std::setw(20) <<  "signal"      << std::setw(15) << Nsig     <<  " ± " << Nsig.error  << std::endl;
  std::cout << std::setw(20) <<  "background"  << std::setw(15) << Nbg      <<  " ± " << Nbg.error   << std::endl;
  std::cout << std::setw(20) <<  "sig + bg"    << std::setw(15) << Nsum     <<  " ± " << Nsum.error  << std::endl;
  std::cout << std::setw(20) <<  "data"        << std::setw(15) << Ndata    <<  " ± " << Ndata.error << std::endl;
  std::cout << "====================  Probability Test ===================  " << std::endl;
  std::cout << std::setw(20) << " sample "             << std::setw(15) << " Kolmogorov "     << std::setw(15)   <<   " Chi2 "           <<  std::setw(15) << "My chi2" << std::endl;
  std::cout << "----------------------------------------------------------  " << std::endl;
  std::cout << std::setw(20) << "signal + background"  << std::setw(15) << prob_kolm_sig_bg   << std::setw(15)   <<   prob_chi2_sig_bg   <<  std::setw(15) << prob_my_sig_bg << std::endl;
  std::cout << std::setw(20) << "signal"               << std::setw(15) << prob_kolm_sig      << std::setw(15)   <<   prob_chi2_sig      <<  std::setw(15) << prob_my_sig    << std::endl;
  std::cout << "----------------------------------------------------------  " << std::endl;
  return c;
}

void draw_all(Selection_t & SEL, Scan_t & DATA, Scan_t & MC, std::vector<ScanRef_t> BGs, int Nbin=80) {
  gStyle->Reset("Pub");
  auto compare = [&](std::string var, int Nb, double xmin,double xmax, std::string extracut="") {
    std::cout << "Draw " << var << std::endl;
    auto c = draw(SEL,DATA,MC,BGs,var, Nbin, xmin, xmax, extracut);
    gSystem->ProcessEvents();
    c->SetTitle(var.c_str());
  };
  compare("ptem"           , Nbin , 0  , 1.2);
  compare("p"              , Nbin , 0.1   , 1.1);
  compare("pt"             , Nbin , 0.15     , 1.1);
  compare("tof"            , Nbin , 2.2     , 5.5);
  compare("cos(theta)"     , Nbin , -0.9    , 0.9);
  compare("cos_theta_mis2" , Nbin , -1    , 1);
  compare("Mpi0"           , Nbin , 0.12 , 0.141       , "Npi0==1");
  compare("Mrho[0]"        , Nbin , 0.5   , 1.1         , "PI0 && Npi0==1");
  compare("acop"           , Nbin , 0     , 3.1415 , "");
  compare("acol"           , Nbin , 0     , 3.1415 , "");
};

void unconstrain_draw_all(const Selection_t & SEL, const Scan_t & DATA, const Simulation_t & MC, int Nbin=80) {
  gStyle->Reset("Pub");
  auto compare = [&](std::string var, int Nb, double xmin,double xmax, std::string extracut="") {
    std::cout << "Draw " << var << std::endl;
    auto c = unconstrain_par_and_draw(SEL,DATA,MC,var, Nbin, xmin, xmax, extracut);
    gSystem->ProcessEvents();
    c->SetTitle(var.c_str());
  };
  compare("ptem"           , Nbin , 0  , 1.2);
  compare("p"              , Nbin , 0.1   , 1.1);
  compare("pt"             , Nbin , 0.15     , 1.1);
  compare("tof"            , Nbin , 2.2     , 5.5);
  compare("cos(theta)"     , Nbin , -0.9    , 0.9);
  compare("cos_theta_mis2" , Nbin , -1    , 1);
  compare("Mpi0"           , Nbin , 0.12 , 0.141       , "Npi0==1");
  compare("Mrho[0]"        , Nbin , 0.5   , 1.1         , "PI0 && Npi0==1");
  compare("acop"           , Nbin , 0     , 3.1415 , "");
  compare("acol"           , Nbin , 0     , 3.1415 , "");
};

struct HistConfig_t {
  int Nbin = 100;
  double xmin = 0;
  double xmax = 1;
  std::string title;
  std::string name;
  void print(void) const {
    std::cout << Nbin << std::endl;
    std::cout << xmin << std::endl;
    std::cout << xmax << std::endl;
    std::cout << name << std::endl;
    std::cout << title << std::endl;
  }
};


TCanvas * draw_background_vs_energy(const Selection_t & SEL, const std::vector<ScanRef_t> BGs, std::string extracut="" ) {
  //hcfg.print();
  //check number of points in the background
  if(BGs.empty()) return nullptr;
  size_t Npnt = BGs[0].get().size();
  if(! std::all_of(BGs.begin(),BGs.end(), [&Npnt] (const auto & bg) {  return Npnt==bg.get().size(); } )) return nullptr;
  std::vector<ibn::valer<double> > CS(Npnt); //cross section
  std::vector<ibn::valer<double> > W(Npnt); //point energy

  std::vector< std::vector< ibn::valer<double> > > CSs;

  for(size_t pnt = 0; pnt!=Npnt; ++pnt) {
    CS[pnt]=ibn::valer<double>{0,0};
     W[pnt]=ibn::valer<double>{0,0};
  }

  for(auto & bg: BGs) {
    auto SR = fold(select(bg, SEL)); 
    std::vector<ibn::valer<double>> cs(Npnt);
    for(size_t pnt=0; pnt!=Npnt; ++pnt) {
      auto & sp = bg.get()[pnt];
      long N = SR[pnt].tt.N;
      auto & tt = bg.get()[pnt].tt;
      ibn::valer<double> eps{double(N)/tt.N0mc, 0};
      eps.error = sqrt( eps.value * (1.0 - eps.value) /double(tt.N0mc) );
      cs[pnt] =  tt.cross_section  * eps;
      W[pnt]=SR[pnt].energy;
    }
    CSs.push_back(cs);
    for(size_t pnt=0; pnt!=Npnt; ++pnt) {
      CS[pnt]+=cs[pnt];
    }
    std::cout << std::endl;
  }
  //drawing different backgrounds
  auto make_graph = [&W](const auto & cs , int color, int marker=21)  {
    TGraphErrors  * g = new TGraphErrors(cs.size());
    for(size_t pnt = 0; pnt!=cs.size(); ++pnt) {
      g->SetPoint(pnt,(W[pnt].value*0.5 - MTAU)/MeV, cs[pnt].value*nb/pb); //convert to MeV and pb
      g->SetPointError(pnt,W[pnt].error*0.5/MeV, cs[pnt].error*nb/pb);
    }
    g->SetMarkerStyle(marker);
    g->SetMarkerColor(color);
    g->SetLineWidth(3);
    g->SetMarkerSize(2);
    g->SetLineColor(color);
    return g;
  };


  //drawing part
  auto c = new TCanvas;
  gStyle->Reset("Pub");
  gStyle->SetOptFit();

  /*
  TGraphErrors  * g = new TGraphErrors(CS.size());
  for(size_t pnt = 0; pnt!=Npnt; ++pnt) {
    g->SetPoint(pnt,(W[pnt].value*0.5 - MTAU)/MeV, CS[pnt].value*nb/pb); //convert to MeV and pb
    g->SetPointError(pnt,W[pnt].error*0.5/MeV, CS[pnt].error*nb/pb);
  }
  g->SetMarkerStyle(21);
  g->SetMarkerColor(kBlue);
  g->SetLineWidth(3);
  g->SetLineColor(kBlue);
  */
  auto g = make_graph(CS,kBlue,21);
  g->Draw("ap");
  g->SetTitle("Visible background cross section (MC)");
  g->GetYaxis()->SetTitle("pb");
  g->GetXaxis()->SetTitle("E-M_{#tau}^{PDG}, MeV");
  g->Fit("pol1");
  g->GetFunction("pol1")->SetLineColor(kRed);


  auto c2 = new TCanvas;
  //c2->SetLogy();
  TMultiGraph * mg = new TMultiGraph("visible_backgrounds", "Visible backgrounds (MC)");
  std::vector palette{1,8,6,11,40,42,46,29};
  auto lbl = new TLegend(0.6,0.8,1.0,1.0);
    std::vector bg_leg_name = {"hadrons","Bhabha", "#mu#mu","#gamma#gamma","EEee","EE#mu#mu","EE#pi#pi","EEKK"};
  for(size_t i = 0;i<BGs.size(); ++i) {
    auto g = make_graph(CSs[i],palette[i],20+i);
    g->SetMarkerSize(2);
    g->SetLineWidth(3);
    mg->Add(g,"lp");
    lbl->AddEntry(g,bg_leg_name[i],"plf");
  }
  mg->Draw("a");
  lbl->Draw();

  return c;
};

std::tuple<ibn::valer<double>, ibn::valer<double> > select(Scan_t & SIG, std::vector<ScanRef_t> BGs, std::string_view sel="") {

  auto get_entries = [](const ScanPoint_t & sp, std::string_view sel="") -> ibn::valer<double> {
    long N = sp.tt.tree->GetEntries(sel.data());
    ibn::valer<double> eps{double(N)/sp.tt.N0mc};
    eps.error = sqrt( eps.value * (1.0 - eps.value) /double(sp.tt.N0mc) );
    auto cs =  sp.tt.cross_section  * eps;
    auto Nexp = cs*sp.gg.luminosity; 
    return Nexp;
  };

  ibn::valer<double> Nsig(0);
  for (const auto & sp : SIG)  Nsig+= get_entries(sp, sel);

  ibn::valer<double> Nbg(0);
  for (const auto & bg : BGs)  
    for(const auto &sp : bg.get()) Nbg+= get_entries(sp, sel);

  return {Nsig, Nbg};
}

std::tuple<ibn::valer<double>, ibn::valer<double> > select_data_sig(const Scan_t & DATA, const Scan_t & SIG, std::string_view sel="") {

  auto get_entries = [](const ScanPoint_t & sp, std::string_view sel="") -> ibn::valer<double> {
    long N = sp.tt.tree->GetEntries(sel.data());
    ibn::valer<double> eps{double(N)/sp.tt.N0mc};
    eps.error = sqrt( eps.value * (1.0 - eps.value) /double(sp.tt.N0mc) );
    auto cs =  sp.tt.cross_section  * eps;
    auto Nexp = cs*sp.gg.luminosity; 
    return Nexp;
  };

  ibn::valer<double> Nsig(0);
  for (const auto & sp : SIG)  Nsig+= get_entries(sp, sel);

  ibn::valer<double> Ndata(0);
  for (const auto & sp : DATA)  {
    ibn::valer<double> n(sp.tt.tree->GetEntries(sel.data()));
    n.error = sqrt(n.value);
    Ndata += n;
  }
  Ndata.error=0;
  return {Nsig,Ndata};
}

std::tuple<ibn::valer<double>, ibn::valer<double> > select(Scan_t & SIG, Scan_t & BG, std::string sel) {
  return select(SIG, {BG} );
}

ibn::valer<double> get_data_sig_value(const Scan_t & DATA, const Scan_t & SIG, const std::string_view cut){
    auto  [Nsig, Ndata]  =  select_data_sig(DATA, SIG, cut);
    ibn::valer<double> x{Nsig.value/sqrt(Ndata)};
    x.error = x.value*std::hypot(Nsig.error/sqrt(Nsig.value), 0.5*Ndata.error/sqrt(Ndata.value));
    return x;
};

ibn::valer<double> get_sig_bg_value(Simulation_t & sim, const std::string_view cut){
    auto  [Nsig, Nbg]  =  select(sim.signal, sim.bgs, cut);
    ibn::valer<double> x{Nsig.value/sqrt(Nsig.value + Nbg.value)};
    x.error = x.value*0.5/(Nsig.value+Nbg.value)*std::hypot( Nbg.error,  Nsig.error*(2*Nbg.value/Nsig.value+1) );
    return x;
};

TGraphErrors * find_optimal_cut_new(Simulation_t & sim, std::string_view sel, std::string_view parfmt, double min, double max, double step) {
  auto g = new TGraphErrors;
  long idx{0};
  for( double par = min; par<max; par+=step ) {
    char par_cut[1024];
    sprintf(par_cut, parfmt.data(), par);  
    std::string cut = std::string(sel) + "&&" + par_cut;
    auto x = get_sig_bg_value(sim, cut);
    g -> SetPoint(idx, par, x.value);
    g -> SetPointError(idx, 0, x.error);
    ++idx;
  }
  //g->Draw("a*");
  g->SetLineWidth(2);
  g->SetMarkerSize(2);
  g->SetMarkerStyle(21);
  return g;
}


std::tuple<TGraph*, double, double, double>  scan_param(const Scan_t &DATA, const Scan_t & SIG, std::string_view sel, std::string_view parfmt, double min, double max, double step) {
  auto g = new TGraphErrors;
  long idx{0};
  for( double par = min; par<max; par+=step ) {
    char par_cut[1024];
    sprintf(par_cut, parfmt.data(), par);  
    std::string cut = std::string(sel) + "&&" + par_cut;
    auto x = get_data_sig_value(DATA,SIG, cut);
    if(std::isfinite(x.value)) { 
      g -> SetPoint(idx, par, x.value);
      g -> SetPointError(idx, 0, x.error);
      ++idx;
    }
  }
  auto  [xmax1, ymax1] = get_maximum(*g);
  std::cout << "xmax = " << xmax1 << "  ymax " << ymax1 << std::endl;

  auto [xmax, ymax, edm] = find_maximum(
      [&DATA,&SIG,&parfmt, &sel] (double par) 
      {  
        char par_cut[1024];
        sprintf(par_cut, parfmt.data(), par);  
        std::string cut = std::string(sel) + "&&" + par_cut;
        auto x = get_data_sig_value(DATA,SIG, cut);
        return x.value;
      },
      std::max(xmax1-2*step,min),
      std::min(max, xmax1+2*step),
      1e-3
      );
  //std::cout << "xmax = " << maximum << "  edm = " << edm << std::endl;
  //g->Draw("a*");
  g->SetLineWidth(2);
  g->SetMarkerSize(2);
  g->SetMarkerStyle(21);
  //g->Print();
  return {g, xmax, ymax, edm};
}


void find_optimal_cuts(Simulation_t & sim, Selection_t & SEL, int i, std::string par, std::string_view parfmt, int Nbin, double min, double max) {
  gStyle->Reset("Pub");
  auto c = new TCanvas(SEL[i].root_title().c_str(), SEL[i].root_title().c_str());
  std::string initial_cut =  SEL.common_cut();
  if(SEL[i].cut!="") initial_cut + "&&" + SEL[i].cut;
  std::string cut =  remove_some_cuts(par, initial_cut);
  //}
}

void find_optimal_cuts(Simulation_t & sim, Selection_t & SEL, int i, std::string par, std::string_view parfmt, double min, double max, double step) {
  gStyle->Reset("Pub");
  //auto cptem = new TCanvas;
  //auto gptem = find_optimal_cut(SIG,BGs, SEL.common_cut + " && Nn==0 && eu && acop>0.1","ptem>%.3f", 0, 0.3,0.01); 
  //for(size_t channel = 0; channel<SEL.size(); ++i) {
  auto c = new TCanvas(SEL[i].root_title().c_str(), SEL[i].root_title().c_str());
  std::string initial_cut =  SEL.common_cut();
  if(SEL[i].cut!="") initial_cut += "&&" + SEL[i].cut;
  std::string cut =  remove_some_cuts(par, initial_cut);
  std::cout << SEL[i].title << "  : " << par <<  "  cut = " <<  cut << std::endl;
  auto gptem = find_optimal_cut_new(sim, cut, parfmt, min, max, step);
  gptem->Draw("ap");
  //}
}

void find_optimal_cuts(Simulation_t & sim, Selection_t & SEL, std::string par, std::string_view parfmt, double min, double max, double step) {
  for(size_t ch = 0; ch!=SEL.size(); ++ch) {
    find_optimal_cuts(sim, SEL, ch, par, parfmt, min,max,step);
  }
}

void find_optimal_cuts(Simulation_t & sim, Selection_t & SEL) {
  find_optimal_cuts(sim, SEL, "ptem", "ptem>%.3f", 0, 1.0, 0.05);
}

void scan_param(const Scan_t & DATA, const Scan_t & SIG, Selection_t & SEL, int i, std::string par, std::string_view parfmt, double min, double max, double step) {
  gStyle->Reset("Pub");
  //auto c = new TCanvas(SEL[i].root_title.c_str(), SEL[i].root_title.c_str());
  new TCanvas;
  std::string initial_cut =  SEL.common_cut();
  std::cout << "SEL[i].cut = " << SEL[i].cut << std::endl;
  std::regex re(R"(\s*)");
  std::smatch sm;
  if(!std::regex_match(SEL[i].cut, sm, re)) {
    initial_cut += ("&&" + SEL[i].cut);
  }
  //std::cout << "initial_cut = " << initial_cut << std::endl;
  std::string cut =  remove_some_cuts(par, initial_cut);
  //std::cout << SEL[i].title << "  : " << par <<  "  cut = " <<  cut << std::endl;
  auto [g,xmax,ymax,edm] = scan_param(DATA,SIG, cut, parfmt, min, max, step);
  g->Draw("ap");
  //}
}

TCanvas * scan_param(const Scan_t & DATA, const Scan_t & SIG, Selection_t & SEL, std::string par, std::string_view parfmt, double min, double max, double step) {
  gStyle->Reset("Pub");
  auto c = new TCanvas(("scan_param_"+par).c_str(), ("Scan parameter + " + par ).c_str());
  TMultiGraph * mg = new TMultiGraph;
  std::cout << "=========== Scaning " << par << "  ======================== " << std::endl;
  std::cout << "common_cut = " << SEL.common_cut() << std::endl;
  TLegend * l = new TLegend(0.6,0.8,1.0,1.0);
  l->SetHeader(par.c_str());
  for(size_t i = 0;i<SEL.size(); ++i) {
    std::cout << SEL[i].title << " cut = " <<  SEL[i].cut << std::endl;
    std::string cut =  remove_some_cuts(par, SEL[i].cut);
    auto [g,xmax,ymax, edm] = scan_param(DATA,SIG, cut, parfmt, min, max, step);
    g->SetLineColor(i+1);
    g->SetMarkerColor(i+1);
    mg->Add(g);
    char buf[1024];
    sprintf(buf, "%s : %s = %.3f, max = %.3f", SEL[i].root_title().c_str(), par.c_str(), xmax, ymax);
    l->AddEntry(g,buf,"lp");
    if(i==0)  {
      mg->Draw("apl");
      l->Draw();
    }
    c->Modified();
    c->Update();
  }
  return c;
}


TCanvas * unconstrain_par_and_draw2d(const Selection_t & SEL, const std::vector<ScanRef_t> bgs, std::string var, int Nbinx, double xmin, double xmax, int Nbiny, double ymin, double ymax, std::string extracut="") {
  auto  c = new TCanvas;
  gStyle->SetPalette(kOcean);

  auto local_select = [Nbinx, xmin, xmax, Nbiny, ymin, ymax, &var]( const ScanPoint_t & sp , std::string cut)  {
    std::string name = "h"+std::to_string(++HISTO_INDEX);
    char expr[1024];
    sprintf(expr, "%s>>%s(%d,%f,%f,%d,%f,%f)", var.c_str(), name.c_str(), Nbinx, xmin,xmax, Nbiny, ymin, ymax);
    sp.tt.tree->Draw(expr,cut.c_str(),"goff");
    return std::shared_ptr<TH2>((TH2*) sp.tt.tree->GetHistogram());
    //TH2* h2  = new TH2F(name.c_str(), var.c_str(), Nbinx, xmin, xmax, Nbiny, ymin, ymax); 
    //h2->Draw("goff");
  };

  using his2ds_t  = std::vector<std::shared_ptr<TH2>>;
  //std::vector< his2ds_t > hsels;
  his2ds_t hsels;
  for(auto & sel : SEL ) { //loop over selection channels
    //std::string cut = SEL.common_cut + " && " +  sel.cut + (extracut == "" ? "" :  (" && " + extracut) );
    std::string cut = remove_some_cuts(var,SEL.common_cut() + " && " +  sel.cut + (extracut == "" ? "" :  (" && " + extracut) ));
    std::cout << cut << std::endl;
    his2ds_t hbgs;
    for(auto & scan : bgs ) {
      his2ds_t hpts;
      for(auto & sp : scan.get() ) {
        hpts.push_back(local_select(sp,cut));
      }
      hbgs.push_back(fold(hpts));
    }
    hsels.push_back(fold(hbgs));
  }
  auto h = fold(hsels);
  h->Draw("colz");
  return c;
}


TCanvas * draw_efficiency(const std::vector<PointSelectionResult_t> & SR) {
  auto * c = get_new_tailed_canvas("efficiency");
  TGraphErrors * g = new TGraphErrors;
  int i=0;
  for(auto & sp : SR ) {
    std::cout << sp.energy.value << " " << sp.tt.efficiency.value << std::endl;
    g->SetPoint(i, sp.energy.value, sp.tt.efficiency.value);
    g->SetPointError(i, sp.energy.error, sp.tt.efficiency.error);
    ++i;
  }
  g->Draw("a*");
  return c;
}


void check_common_cuts(const Scan_t & SCAN, std::string initial_cut) {
  std::map<std::string, Scan_t> R;
  R["base"] = measure_efficiency(SCAN,initial_cut);
  auto proceed = [&](std::string s) {
    R[s] = measure_efficiency(SCAN, remove_some_cuts(s, initial_cut));
  };
  proceed("Nc");
  proceed("Nn");
  proceed("p");
  proceed("pt");
  proceed("theta");
  proceed("cos_theta_mis2");
  std::map<std::string, TGraphErrors *> mgraph;
  TMultiGraph * mg0 = new TMultiGraph;
  TMultiGraph * mg = new TMultiGraph;
  int color=0;
  int marker=20;
  for(auto & [par, sr] : R) {
    auto g = new TGraphErrors;
    auto g0 = new TGraphErrors;
    ++color;
    ++marker;
    g->SetLineWidth(3);
    g->SetMarkerSize(2);
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(marker);

    g0->SetLineWidth(3);
    g0->SetMarkerSize(2);
    g0->SetLineColor(color);
    g0->SetMarkerColor(color);
    g0->SetMarkerStyle(marker);
    int i=0;
    for(auto & sp : sr) {
      g->SetPoint(i, sp.energy.value, sp.tt.effcor.value);
      g->SetPointError(i, sp.energy.error, sp.tt.effcor.error);
      g0->SetPoint(i, sp.energy.value, sp.tt.efficiency.value);
      g0->SetPointError(i, sp.energy.error, sp.tt.efficiency.error);
      ++i;
    }
    mgraph[par]=g;
    mg->Add(g,"lp");
    mg0->Add(g0,"lp");
  }
  auto c_eff = get_new_tailed_canvas("efficiency");
  mg0->Draw("a");
  auto c_cor = get_new_tailed_canvas("effcor");
  mg->Draw("a");
  auto c_dif = get_new_tailed_canvas("effcor variation");
  TMultiGraph * mg2 = new TMultiGraph;
  auto & base = R["base"];
  color=0;
  marker=20;
  TLegend * l = new TLegend(0.8,0.8,1.0,1.0);
  std::ofstream ofs("effcor_variation_cut.txt");
  char buf[128];
  for(auto & [par, sr] : R) {
    auto g = new TGraphErrors;
    ++color;
    ++marker;
    g->SetLineWidth(2);
    g->SetMarkerSize(3);
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(marker);
    int i=0;
    sprintf(buf,"%20s",par.c_str());
    ofs << buf;
    for(auto & sp : sr) {
      double cor =  sp.tt.effcor.value/base[i].tt.effcor.value;
      double error = cor*hypot(base[i].tt.effcor.error/ base[i].tt.effcor.value, sp.tt.effcor.error/sp.tt.effcor.value);
      g->SetPoint(i, sp.energy.value, cor);
      g->SetPointError(i, sp.energy.error,error);
      ++i;
      sprintf(buf,"%15.8f  %15.8f", cor, error);
      ofs << buf;
    }
    ofs << endl;
    mgraph[par]=g;
    l->AddEntry(g,par.c_str(), "lp");
    mg2->Add(g,"lp");
  }
  mg2->Draw("a");
  l->Draw();
}

TCanvas * draw_ee_efficiency(const std::vector<ScanPoint_t> & SR) {
  auto * c = get_new_tailed_canvas("efficiency");
  TGraphErrors * g = new TGraphErrors;
  int i=0;
  for(auto & sp : SR ) {
    std::cout << sp.energy.value << " " << sp.bb.efficiency.value << std::endl;
    g->SetPoint(i, sp.energy.value, sp.bb.efficiency.value);
    g->SetPointError(i, sp.energy.error, sp.bb.efficiency.error);
    ++i;
  }
  g->Draw("a*");
  g->GetXaxis()->SetTitle("W, GeV");
  g->GetYaxis()->SetTitle("#varepsilon^{ee}");
  return c;
}

TCanvas * draw_gg_efficiency(const std::vector<ScanPoint_t> & SR, std::string gopt="") {
  TCanvas * c = nullptr;
  if(gopt!="same") {
    c = get_new_tailed_canvas("#gamma#gamma efficiency");
  }
  TGraphErrors * g = new TGraphErrors;
  int i=0;
  for(auto & sp : SR ) {
    std::cout << sp.energy.value << " " << sp.gg.efficiency.value << std::endl;
    g->SetPoint(i, sp.energy.value, sp.gg.efficiency.value);
    g->SetPointError(i, sp.energy.error, sp.gg.efficiency.error);
    ++i;
  }
  if(gopt!="same") {
    g->Draw("a*");
    g->GetXaxis()->SetTitle("W, GeV");
    g->GetYaxis()->SetTitle("#varepsilon^{#gamma#gamma}");
  }
  else {
    g->Draw("*");
  }
  return c;
}

template<typename Proj1, typename Proj2>
TGraphErrors * efficiency_compare(std::vector<ScanPoint_t> & MC1, std::vector<ScanPoint_t> & MC2, Proj1 proj1, Proj2 proj2, std::string cut1, std::string cut2) {
  measure_efficiency(MC1, proj1, cut1);
  measure_efficiency(MC2, proj2, cut2);
  TGraphErrors * g = new TGraphErrors;
  //int i=0;
  for(auto & sp1 :  MC1 ) {
    for(auto & sp2 : MC2) {
      if( fabs(sp1.energy.value -  sp2.energy.value) < 0.0001 ) {
        auto & ds1 = std::invoke(proj1,sp1);
        auto & ds2 = std::invoke(proj2,sp2);
        ibn::valer<double> ratio = ds1.efficiency/ds2.efficiency;
        ratio.error =  ratio.value*hypot(ds1.efficiency.error/ds1.efficiency.value, ds2.efficiency.error/ds2.efficiency.value );
        AddPoint(g,sp1.energy, ratio);
        //g->SetPoint(i, sp1.energy.value, ratio.value);
        //g->SetPointError(i, sp1.energy.error, ratio.error);
        //++i;
      }
    }
  }
  return g;
}

TCanvas * draw_tt_ee_compare(std::vector<ScanPoint_t> & SIG, std::vector<ScanPoint_t> & BB, std::string sig_cut, std::string bb_cut) {
  auto g = efficiency_compare(SIG,BB,&ScanPoint_t::tt, &ScanPoint_t::bb, sig_cut, bb_cut);
  auto c = get_new_tailed_canvas("Efficiency compare: #varepsilon_{#tau#tau} / #varepsilon_{ee}");
  g->Draw("a*");
  g->GetXaxis()->SetTitle("W, GeV");
  g->GetYaxis()->SetTitle("#varepsilon_{#tau#tau} / #varepsilon_{ee}");
  gStyle->SetOptFit(kTRUE);
  g->Fit("pol1");
  return c;
}

TCanvas * draw_tt_gg_compare(std::vector<ScanPoint_t> & SIG, std::vector<ScanPoint_t> & GG, std::string sig_cut, std::string bb_cut) {
  auto g = efficiency_compare(SIG,GG,&ScanPoint_t::tt, &ScanPoint_t::gg, sig_cut, bb_cut);
  auto c = get_new_tailed_canvas("Efficiency compare: #varepsilon_{#tau#tau} / #varepsilon_{#gamma #gamma}");
  g->Draw("a*");
  g->GetXaxis()->SetTitle("W, GeV");
  g->GetYaxis()->SetTitle("#varepsilon_{#tau#tau} / #varepsilon_{#gamma#gamma}");
  gStyle->SetOptFit(kTRUE);
  g->Fit("pol1");
  return c;
}

TCanvas * draw_ee_gg_compare(std::vector<ScanPoint_t> & BB, std::vector<ScanPoint_t> & GG, std::string sig_cut, std::string bb_cut) {
  auto g = efficiency_compare(BB,GG,&ScanPoint_t::bb, &ScanPoint_t::gg, sig_cut, bb_cut);
  auto c = get_new_tailed_canvas("Efficiency compare: #varepsilon_{ee} / #varepsilon_{#gamma #gamma}");
  g->Draw("a*");
  g->GetXaxis()->SetTitle("W, GeV");
  g->GetYaxis()->SetTitle("#varepsilon_{ee} / #varepsilon_{#gamma#gamma}");
  gStyle->SetOptFit(kTRUE);
  g->Fit("pol1");
  return c;
}

TCanvas * draw_gg_ee_compare(std::vector<ScanPoint_t> & GG, std::vector<ScanPoint_t> & BB, std::string gg_cut, std::string ee_cut) {
  auto g = efficiency_compare(GG,BB,&ScanPoint_t::gg,&ScanPoint_t::bb, gg_cut,  ee_cut);
  auto c = get_new_tailed_canvas("Efficiency compare: #varepsilon_{gg} / #varepsilon_{ee}");
  g->Draw("a*");
  g->GetXaxis()->SetTitle("W, GeV");
  g->GetYaxis()->SetTitle("#varepsilon_{gg} / #varepsilon_{ee}");
  gStyle->SetOptFit(kTRUE);
  g->Fit("pol1");
  return c;
}

//TCanvas * draw_tt_ee_compare(const std::vector<ScanPoint_t> & SIG, std::vector<ScanPoint_t> & BB, std::string sig_cut, std::string bb_cut) {
//  auto S = measure_efficiency(SIG,sig_cut);
//  measure_ee_luminosity(BB,bb_cut);
//  TCanvas * c = nullptr;
//  c = get_new_tailed_canvas("#gamma#gamma efficiency");
//  TGraphErrors * g = new TGraphErrors;
//  int i=0;
//  for(auto & sig : S ) {
//    for(auto & bb : BB) {
//      if( fabs(sig.energy.value -  bb.energy.value) < 0.0001 ) {
//        ibn::valer<double> ratio = sig.tt.efficiency/bb.bb.efficiency;
//        ratio.error =  ratio.value*hypot(sig.tt.efficiency.error/sig.tt.efficiency.value, bb.bb.efficiency.error/bb.bb.efficiency.value );
//        //std::cout << bb.bb.efficiency.error << std::endl;
//        //std::cout << sig.tt.efficiency.error << std::endl;
//        //std::cout << ratio.error << std::endl;
//        //std::cout << sig.energy.value << " " << bb.energy.value << " " << ratio.value << " " << ratio.error  << std::endl;
//        g->SetPoint(i, sig.energy.value, ratio.value);
//        g->SetPointError(i, sig.energy.error, ratio.error);
//        ++i;
//      }
//    }
//  }
//  g->Draw("a*");
//  g->GetXaxis()->SetTitle("W, GeV");
//  g->GetYaxis()->SetTitle("#varepsilon^{#tau#tau}/#varepsilon^{ee}");
//  gStyle->SetOptFit(kTRUE);
//  g->Fit("pol1");
//  return c;
//}

void compare2D_base ( Scan_t & MC, Scan_t & D, std::string var, std::string cut="" ){
  MC[0].tt.tree->SetMarkerColor(kBlue); 
  MC[0].tt.tree->SetLineColor(kBlue); 
  MC[0].tt.tree->SetLineWidth(3); 
  MC[0].tt.tree->SetMarkerStyle(21); 
  long Nsig = MC[0].tt.tree->Draw(var.c_str(),  cut.c_str(),"");

  D[0].tt.tree->SetMarkerColor(kRed); 
  D[0].tt.tree->SetLineColor(kRed); 
  D[0].tt.tree->SetLineWidth(3); 
  D[0].tt.tree->SetMarkerStyle(21); 
  long n = D[0].tt.tree->Draw(var.c_str(), cut.c_str(),"same"); 
  std::cout << "Nsig = " << Nsig << "  Nbg  = " << n << std::endl;
}
TCanvas * compare2D ( Scan_t & MC, Scan_t & D, std::string var, std::string cut="" ){
  auto c = new TCanvas; 
  compare2D_base(MC, D, var,cut );
  return c;
}

//draw ptem:cos_theta_mis2 Emis:cos_theta_mis2  acop:cos_theta_mis2
TCanvas * compare2D_2 ( Scan_t & MC, Scan_t & D, std::string cut="" ){
  auto c = new TCanvas("compare2D_2c","compare2D_2",  3600,2000); 
  c->Divide(4,2);
  int idx=0;
  c->cd(++idx);
  compare2D_base(MC, D, "ptem:cos_theta_mis2",cut);
  c->cd(++idx);
  compare2D_base(MC, D, "Emis:cos_theta_mis2",cut);
  c->cd(++idx);
  compare2D_base(MC, D, "ptem:acop",cut);
  c->cd(++idx);
  compare2D_base(MC, D, "Emis:acop",cut);
  c->cd(++idx);
  compare2D_base(MC, D, "acop:cos_theta_mis2",cut);
  c->cd(++idx);
  compare2D_base(MC, D, "acol:cos_theta_mis2",cut);
  c->cd(++idx);
  compare2D_base(MC, D, "M2:ptem",cut);
  c->cd(++idx);
  compare2D_base(MC, D, "ptem:Enmax",cut);
  return c;
}


// Рисует график зависимости эффективности регистрации от энергии с разными катами
void check_common_cuts2(const Scan_t & SCAN, const std::vector<std::string> & cuts, std::string file_prefix="", const std::vector<std::string>  titles ={}) {
  //std::multimap<std::string, Scan_t> R;
  std::vector<Scan_t> R;
  for(auto & cut : cuts) {
    auto c = cut;
    //c.erase(std::remove(c.begin(), c.end(), ' '), c.end());
    //std::cout << c <<std::endl;
    //R.insert({cut, measure_efficiency(SCAN,cut)});
    R.push_back(measure_efficiency(SCAN,c));
  }

  std::map<std::string, TGraphErrors *> mgraph;
  TMultiGraph * mg0 = new TMultiGraph;
  TMultiGraph * mg = new TMultiGraph;
  //int color=0;
  int marker=20;
  std::vector<int> color={kBlack, kRed, kBlue, kGreen+2, kOrange+7, kPink+7, kMagenta, kBlue-6, kCyan};
  int idx=0;
  auto escale  = [](double W) {
    return 0.5*(W-2.0*MTAU)*1000.0;
  };
  auto descale  = [](double dW) {
    return 0.5*dW*1000.0;
  };
  for(auto & sr : R) {
    auto g = new TGraphErrors;
    auto g0 = new TGraphErrors;
    ++marker;
    g->SetLineWidth(3);
    g->SetMarkerSize(2);
    g->SetLineColor(color[idx%color.size()]);
    g->SetMarkerColor(color[idx%color.size()]);
    g->SetMarkerStyle(marker);

    g0->SetLineWidth(3);
    g0->SetMarkerSize(2);
    g0->SetLineColor(color[idx%color.size()]);
    g0->SetMarkerColor(color[idx%color.size()]);
    g0->SetMarkerStyle(marker);
    int i=0;
    for(auto & sp : sr) {
      g->SetPoint(i, escale(sp.energy.value), (sp.tt.effcor.value-1.0)*100.0);
      g->SetPointError(i, descale(sp.energy.error), sp.tt.effcor.error*100.0);
      g0->SetPoint(i, escale(sp.energy.value), sp.tt.efficiency.value*100.0);
      g0->SetPointError(i, descale(sp.energy.error), sp.tt.efficiency.error*100.0);
      ++i;
    }
    mgraph[cuts[idx]]=g;
    mg->Add(g,"lp");
    mg0->Add(g0,"lp");
    ++idx;
  }
  TMultiGraph * mg2 = new TMultiGraph;
  auto & base = *R.begin();
  marker=20;
  //TLegend * l = new TLegend(0.8,0.8,1.0,1.0);
  TLegend * l = new TLegend(0.33,0.15, 0.66,0.85);
  std::ofstream ofs("effcor_variation_cut.txt");
  char buf[65535];
  idx=0;
  for(auto & sr : R) {
    //std::cout << idx << std::endl;
    std::string par = cuts[idx];
    //std::cout << par << std::endl;
    auto g = new TGraphErrors;
    ++marker;
    g->SetLineWidth(3);
    g->SetMarkerSize(3);
    g->SetLineColor(color[idx%color.size()]);
    g->SetMarkerColor(color[idx%color.size()]);
    g->SetMarkerStyle(marker);
    int i=0;
    sprintf(buf,"%-20s",par.c_str());
    ofs << buf;
    std::cout << buf;
    //std::cout << " base.size = " << base.size() << std::endl;
    for(auto & sp : sr) {
 //     std::cout << sp.energy.value << std::endl;
      double cor =  sp.tt.effcor.value/base[i].tt.effcor.value;
      double error = cor*hypot(base[i].tt.effcor.error/ base[i].tt.effcor.value, sp.tt.effcor.error/sp.tt.effcor.value);
      g->SetPoint(i, escale(sp.energy.value), (cor-1.0)*100.0);
      g->SetPointError(i, descale(sp.energy.error),error*100.0);
      ++i;
      sprintf(buf,"%15.8f  %15.8f     ", (cor-1.0)*100.0, error*100);
      ofs << buf;
      std::cout << buf;
    }
    ofs << std::endl;
    std::cout << std::endl;
    mgraph[par]=g;
    if(idx<titles.size()) {
        l->AddEntry(g,titles[idx].c_str(), "lp");
    } else {
      if(par.size() < 30) {
        l->AddEntry(g,par.c_str(), "lp");
      } else {
        std::string shrinked_name = par.substr(0,27)+"...";
        l->AddEntry(g, shrinked_name.c_str(), "lp");
      }
    }
    mg2->Add(g,"lp");
    ++idx;
  }
  auto c_eff = get_new_tailed_canvas("efficiency");
  mg0->Draw("a");
  mg0->GetXaxis()->SetTitle("E-M_{#tau}, MeV");
  mg0->GetYaxis()->SetTitle("#varepsilon^{cut}_{i}, %");
  mg0->GetYaxis()->SetDecimals();
  l->Draw();

  auto c_cor = get_new_tailed_canvas("effcor");
  mg->Draw("a");
  mg->GetXaxis()->SetTitle("E-M_{#tau}, MeV");
  mg->GetYaxis()->SetTitle("#varepsilon^{cut}_{i}/#varepsilon^{cut}_{M#tau}-1, %");
  mg->GetYaxis()->SetDecimals();
  l->Draw();

  auto c_dif = get_new_tailed_canvas("effcor variation");
  mg2->Draw("a");
  mg2->GetYaxis()->SetTitle("#varepsilon^{cut}_{i}/#varepsilon^{cut}_{M#tau} / #varepsilon^{0}_{i}/#varepsilon^{0}_{M#tau}  - 1, %");
  mg2->GetYaxis()->SetDecimals();
  mg2->GetXaxis()->SetTitle("E-M_{#tau}, MeV");
  l->Draw();
  if(file_prefix!="") {
    c_eff->SaveAs((file_prefix+"eff.pdf").c_str());
    c_cor->SaveAs((file_prefix+"effcor.pdf").c_str());
    c_dif->SaveAs((file_prefix+"effcor-diff.pdf").c_str());
  }
}

void registration_efficiency_sys(const Scan_t & SCAN, const Selection_t & SEL,  std::string file_prefix="sys-") {
  std::vector< std::pair<std::string, std::string > > pars
  {
      //{""                         , "nocut"}             ,
      //{"NnE50==0"                 , "N_{n}=0"}           ,
      //{"Ncg==2"                   , "N_{cg}=2"}           ,
      //{"Ncc==2"                   , "N_{cc}=2"}           ,
      //{"barrel"                   , "|cos(#theta)|<0.8"}            ,
      ////{"q[0]==-1 && q[1]==1"      , "opposite charge"}   ,
      //{"pt[0]>0.2 && pt[1]>0.2"   , "p_{t}>0.2"}            ,
      //{"p[0]<1.1 && p[1]<1.1"     , "p<1.1"}             ,
      //{"missed_photon_angle"      , "cos(#theta_{mis})"} ,
      //{"2.5 < tof && tof< 5.5"    , "2.5<tof<5.5 ns"}       ,
      //{"ptem50>0.25"              , "ptem>0.25"}         ,
      //{"ptem50<1.1"               , "ptem<1.1"}         ,
      //{"E[0]>0.025 && E[1]>0.025" , "E>0.1"}             ,
      //{"z[0]<5 && z[1]<5"         , "|z|<5"}             ,
      //{"vxy[0]<0.5 && vxy[1]<0.5" , "#rho<0.5"}          ,
      //{"Nng==2"                   , "N^{good}_{#gamma} = 2" },
      //{"eX"                       , "PID"}                ,
      //{SEL[0].cut                 , "eX"}               ,
      //{SEL[1].cut                 , "e#rho"}               ,
      {SEL[0].cut || SEL[1].cut,  "eX + e#rho"}               ,
  };
  std::vector<std::string> cuts, titles;
  for(auto & [cut, title] : pars) {
    cuts.push_back(cut);
    titles.push_back(title);
  }
  check_common_cuts2(SCAN, cuts,file_prefix, titles); 
}



struct Config {
  std::string name = "";
  std::string lum = "gg"; //gg or ee
  std::string energy_version ="ems3";
  std::string gglum_extra_cut = "";
  std::string eelum_extra_cut = "";
  std::string taufit_extra_opt = "";
  //double ecor=+0.011;
  //double spread=1.258;
  //bool free_energy=true;
  //bool free_luminosity=true;
  //bool free_effcor=true;
  //bool minos=true;
  //bool pdgshift=true;
} CFG;

void clear_canvas(std::string regex="") {
  //auto lst = ;
  std::vector<TCanvas*> to_delete;
  std::regex re(regex);
  std::smatch sm;
  for(auto l : *(gROOT->GetListOfCanvases())) {
    std::string name = l->GetName();
    if(std::regex_search(name, sm, re)) {
      to_delete.push_back((TCanvas*)l);
    }
  }
  for(auto l : to_delete) delete l;
}

//void parameter_example(Selection_t & S)
//{
//  fold_and_draw(DATA,"ptem","Nc==2 && Nn==0 &&" + S.common_cut(),"NORM");
//  fold_and_draw(SIGNAL,"ptem","Nc==2 && Nn==0 &&" + S.common_cut(),"NORM SAME");
//}



//TCanvas * draw_tt_ee_compare(std::string sig_cut, std::string bb_cut) {
//  return draw_tt_ee_compare(SIGNAL, BB, sig_cut, bb_cut);
//}
//
//TCanvas * draw_tt_gg_compare(std::string cut1, std::string cut2) {
//  return draw_tt_gg_compare(SIGNAL, GG, cut1, cut2);
//}
//
//TCanvas * draw_gg_ee_compare(std::string cut1, std::string cut2) {
//  return draw_gg_ee_compare(GG, BB, cut1, cut2);
//}


void test_for_mixing(Scan_t & D, std::string extracut = "") {
  TCanvas * c = new TCanvas;
  c->Divide(2,2);
  int i=0;
  auto get_graph  = [&](std::string cut1, std::string cut2) {
    TGraphErrors * g1  = new TGraphErrors;
    int i=0;
    for(auto & d : D) {
      std::string cut = cut1;
      std::string truth = cut2;
      double N0  = d.tt.tree->GetEntries(cut.c_str());
      double N1  =  d.tt.tree->GetEntries((cut+truth).c_str());
      g1->SetPoint(i, d.energy.value, 100*N1/N0);
      g1->SetPointError(i, d.energy.error, 100*N1/N0*sqrt(1.0/N1 + 1.0/N0));
      ++i;
    }
    return g1;
  };
  gStyle->SetOptFit();

  auto make_and_draw_graph = [&](std::string cut1, std::string cut2, int col=1, std::string gopt="a*") {
    TGraphErrors * g1 = get_graph(cut1,cut2);
    g1->SetMarkerColor(col);
    g1->SetLineColor(col);
    g1->Draw(gopt.c_str());
    g1->GetXaxis()->SetTitle("W_{cm}, MeV");
    g1->GetYaxis()->SetTitle("missid, %");
    g1->Fit("pol1");
  };


  c->cd(1);
  make_and_draw_graph("e0 && X1", " && pid[1] == -11", kBlack,"a*");
  c->cd(2);
  make_and_draw_graph("e1 && X0", " && pid[0] == 11", kRed, "*a");

  c->cd(3);
  make_and_draw_graph("e0 && X1", " && pid[0] != 11", kBlack,"a*");
  c->cd(4);
  make_and_draw_graph("e1 && X0", " && pid[1] != -11", kRed, "a*");

}

void count(Scan_t & D, std::string cut) {
  long N = D[0].tt.tree->Draw("pid[0]:pid[1]:mother_pid[0]:mother_pid[1]", cut.c_str(),"goff");
  double * pid0 = D[0].tt.tree->GetV1();
  double * pid1 = D[0].tt.tree->GetV2();
  double * mpid0 = D[0].tt.tree->GetV3();
  double * mpid1 = D[0].tt.tree->GetV4();
  std::map<std::string, long> m; //basic count
  std::map<std::pair<int, int> , long> om;
  std::set<std::pair<long, std::string>, std::greater<std::pair<long, std::string>> > sm;
  std::set<std::pair<long, std::pair<int, int>>,  std::greater<std::pair<long, std::pair<int, int>>> > som;
  for(int i=0;i<N;++i) {
    if(mpid0[i]==15 && mpid1[i]==-15) {
      if(pid0[i] == 11 && pid1[i]==-11) ++m["ee"];

      if(pid0[i] == 11 && pid1[i]==-13) ++m["eu"];
      if(pid0[i] == 13 && pid1[i]==-11) ++m["eu"];

      if(pid0[i] == 11 && pid1[i]==211) ++m["epi"];
      if(pid0[i] == -211 && pid1[i]==-11) ++m["epi"];

      if(pid0[i] == 11 && pid1[i]==321) ++m["eK"];
      if(pid0[i] == -321 && pid1[i]==-11) ++m["eK"];

      if(pid0[i] == 13 && pid1[i]==-13) ++m["uu"];

      if(pid0[i] == 13 && pid1[i]==211) ++m["upi"];
      if(pid0[i] == -211 && pid1[i]==-13) ++m["upi"];

      if(pid0[i] == 13 && pid1[i]==321) ++m["uK"];
      if(pid0[i] == -321 && pid1[i]==-13) ++m["uK"];


      if(pid0[i] == -211 && pid1[i]==211) ++m["pipi"];

      if(pid0[i] == -211 && pid1[i]==321) ++m["piK"];
      if(pid0[i] == -321 && pid1[i]==211) ++m["piK"];

      if(pid0[i] == -321 && pid1[i]==321) ++m["KK"];

      switch(int(pid0[i])) {
        case 11:
        case 13:
        case -211:
        case -321:
          break;
        default:
          ++om[{int(pid0[i]),int(pid1[i])}];
          ++m["other"];
      }
      switch(int(pid1[i])) {
        case -11:
        case -13:
        case 211:
        case 321:
          break;
        default:
          ++om[{int(pid0[i]),int(pid1[i])}];
          ++m["other"];
      }
    }
    else {
      std::cout << mpid0[i] << " " << mpid1[i] <<  "  " << pid0[i] << "  " << pid1[i] << std::endl;
    }
  }
  long n=0;
  long nother=0;
  for(auto  & [channel, count] : m) {
    sm.insert( {count, channel} );
    n+=count;
  }

  for(auto & [p, count] : om) {
    som.insert( {count, p} );
    nother+=count;
  }


  int width=70;
  auto print_double_line = [&width] () {
    for(int i=0;i<70; ++i) std::cout << "=";
    std::cout << "\n";
  };
  /*
  print_double_line();
  std::cout << myfmt("%16s %20d %20.3f%%", "Total number", N, 100.0) << std::endl;
  print_double_line();
  for(auto  & [channel, count] : m) {
    std::cout << myfmt("%16s %20d %20.3f%%", channel.c_str(), count, 100*double(count)/N) << std::endl;
  }
  print_double_line();
  std::cout << myfmt("%16s %20d %20.3f%%", "Channels above", n, 100*double(n)/N) << std::endl;
  std::cout << myfmt("%16s %20d %20.3f%%", "Other channels", nother, 100*double(nother)/N) << std::endl;
  //std::cout << "Other channels:  " << nother << std::endl;
  print_double_line();
  for(auto & [p, count] : om) {
    std::cout << myfmt("%7d,%7d %20d %20.3f%%", p.first, p.second, count, 100*double(count)/N) << std::endl;
  }
  */


  print_double_line();
  std::cout << myfmt("%16s %20d %20.3f%%", "Total number", N, 100.0) << std::endl;
  print_double_line();
  for(auto  & [count,channel] : sm) {
    std::cout << myfmt("%16s %20d %20.3f%%", channel.c_str(), count, 100*double(count)/N) << std::endl;
  }
  print_double_line();
  std::cout << myfmt("%16s %20d %20.3f%%", "Channels above", n, 100*double(n)/N) << std::endl;
  std::cout << myfmt("%16s %20d %20.3f%%", "Other channels", nother, 100*double(nother)/N) << std::endl;
  //std::cout << "Other channels:  " << nother << std::endl;
  print_double_line();
  for(auto & [count, p] : som) {
    std::cout << myfmt("%7d,%7d %20d %20.3f%%", p.first, p.second, count, 100*double(count)/N) << std::endl;
  }
  print_double_line();

}

struct event_t{
  int mpid0;
  int mpid1;
  int pid0;
  int pid1;
  int q0;
  int q1;
  double c0;
  double c1;
};

inline  bool operator<(const event_t & e1, const event_t & e2) {
  if      ( e1.mpid0 <  e2.mpid0 ) return true;
  else if ( e2.mpid0 <  e1.mpid0 ) return false;
    else if ( e1.mpid0 == e2.mpid0 ) {
      if      ( e1.mpid1 <  e2.mpid1 ) return true;
      else if ( e2.mpid1 <  e1.mpid1 ) return false;
        else if ( e1.mpid1 == e2.mpid1 ) {
          if      ( e1.pid0 < e2.pid0 )  return true;
          else if ( e2.pid0 < e1.pid0 )  return false;
          else if ( e1.pid0 == e2.pid0 ) return e1.pid1 < e2.pid1;
        }
    }
    std::cout << "ERROR: Unable to compare " << std::endl;
    return false;
}

size_t count_multibyte_symbols(const std::string & s) {
  size_t count=0;
  for(size_t i=0;i<s.size();++i) {
    if( ((s[i] >> 6) & 0b11) == 0b10  ) {
      count++;
    }
  }
  return count;
}

size_t utf_size(const std::string  & s) {
  return s.size()-count_multibyte_symbols(s);
}

std::string adjast_utf_string(std::string  s, int width=0) {
  int count=0;
  for(size_t i=0;i<s.size();++i) {
    if( ((s[i] >> 6) & 0b11) == 0b10  ) {
      count++;
    }
  }
  width-=s.size();
  if(width<0) width=0;
  return s+std::string(width+count,' ');
}

std::map<std::string, std::vector<double> >  make_vectors(Scan_t & D, std::string cut) {
  std::vector<double> pid0;
  std::vector<double> pid1;
  std::vector<double> mpid0;
  std::vector<double> mpid1;
  std::vector<double> q0;
  std::vector<double> q1;
  std::map<std::string, std::vector<double> > v;
  for(const auto &  d : D) {
    auto & t  = d.tt.tree;
    long N = t->Draw("pid[0]:pid[1]:mother_pid[0]:mother_pid[1]", cut.c_str(),"goff");
    for(long i =0;i<N;++i) {
      v["pid0"].push_back(t->GetV1()[i]);
      v["pid1"].push_back(t->GetV2()[i]);
      v["mpid0"].push_back(t->GetV3()[i]);
      v["mpid1"].push_back(t->GetV4()[i]);
    }
    N = t->Draw("q[0]:q[1]:cos(theta[0]):cos(theta[1])", cut.c_str(),"goff");
    for(long i =0;i<N;++i) {
      v["q0"].push_back(t->GetV1()[i]);
      v["q1"].push_back(t->GetV2()[i]);
      v["c0"].push_back(t->GetV3()[i]);
      v["c1"].push_back(t->GetV4()[i]);
    }
  }
  return v;
}

void count2(Scan_t & D, std::string cut="", std::string filter="") {
  /*
  auto & t  = D.tt.tree;

  long N = t->Draw("pid[0]:pid[1]:mother_pid[0]:mother_pid[1]", cut.c_str(),"goff");
  auto make_vector = [](double * ptr, long N) {
    return std::vector<double>(ptr, ptr+N);
  };
  auto pid0 = make_vector(t->GetV1(), N);
  auto pid1 = make_vector(t->GetV2(), N);
  auto mpid0 = make_vector(t->GetV3(), N);
  auto mpid1 = make_vector(t->GetV4(), N);

  long N1 = t->Draw("pid[0]:pid[1]:cos(theta[0]):cos(theta[1])", cut.c_str(),"goff");
  auto p0 = make_vector(t->GetV1(), N1);
  auto p1 = make_vector(t->GetV2(), N1);
  auto c0 = make_vector(t->GetV3(), N1);
  auto c1 = make_vector(t->GetV4(), N1);

  long N2 = t->Draw("pid[0]:pid[1]:q[0]:q[1]", cut.c_str(),"goff");
  auto pp0 = make_vector(t->GetV1(), N2);
  auto pp1 = make_vector(t->GetV2(), N2);
  auto q0 = make_vector(t->GetV3(), N2);
  auto q1 = make_vector(t->GetV4(), N2);

  //check pid
  if(N != N1 || N1 != N2 ) {
    std::cout << "ERROR: different length of the arrays" << std::endl;
  }
  */
  auto V = make_vectors(D,cut);
  auto pid0 = V["pid0"]; 
  auto pid1 = V["pid1"]; 
  auto mpid0 = V["mpid0"]; 
  auto mpid1 = V["mpid1"]; 
  auto q0 = V["q0"]; 
  auto q1 = V["q1"]; 
  auto c0 = V["c0"]; 
  auto c1 = V["c1"]; 

  std::map<event_t, long> m;
  std::multimap<long, event_t, std::greater<long> > rm;
  std::multimap<long, event_t, std::greater<long> > rm_unknown;
  std::multimap<long, event_t, std::greater<long> > rm_tt,rm_not_tt, rm_tt_eX, rm_tt_not_eX, head;

  long N  = pid0.size();
  for(long i =0; i<N; ++i) {
   // std::cout << pid0[i] << " " << pid1[i] << " " << mpid0[i] << " " << mpid1[i] << " " << c0[i] << " " << c1[i] << " "<< q0[i] << " " <<  q1[i] << std::endl;
    //if( pid0[i] != p0[i] || pid1[i]!=p1[i] || pp0[i]!=p0[i] || pp1[i]!=p1[i] ) {
    //  std::cout << "ERROR: mismatch of pid" << std::endl;
    //}
    m[{int(mpid0[i]), int(mpid1[i]), int(pid0[i]), int(pid1[i]),int(q0[i]), int(q1[i]), c0[i], c1[i]}]++;
  }

  for(auto & [e,n] : m) {
    //std::cout << e.mpid0 << " " << e.mpid1 << " " << e.pid0 << " " << e.pid1 << " " << e.q0 << " " << e.q1 << "   :   " << n << endl;
    if( e.mpid0 == 0 || e.mpid1 == 0 || e.pid0 ==0  || e.pid1 == 0) {
      rm_unknown.insert({n,e});
    } else {
      rm.insert({n,e});
      if(filter.find("tt")!= std::string::npos && e.mpid0==15 && e.mpid1==-15) {
        rm_tt.insert({n,e});
        if(filter.find("eX") != std::string::npos && ( (e.pid0 == 11 && e.pid1!=-11 ) || (e.pid0 != 11 && e.pid1==-11 ) )) {
          rm_tt_eX.insert({n,e});
        } else {
          rm_tt_not_eX.insert({n,e});
        }
      }
      else {
        rm_not_tt.insert({n,e});
      }
    }

  }

  auto make_names = [](const event_t & e ) -> std::tuple<std::string, std::string, std::string, std::string> {
    auto im0 = PdgTable.find(e.mpid0);
    auto im1 = PdgTable.find(e.mpid1);
    auto ip0 = PdgTable.find( e.pid0);
    auto ip1 = PdgTable.find( e.pid1);
    std::string m0;
    std::string m1;
    std::string p0;
    std::string p1;

    if(im0==PdgTable.end()) m0=std::to_string(e.mpid0);
    else m0=im0->second;

    if(im1==PdgTable.end()) m1=std::to_string(e.mpid1);
    else m1=im1->second;

    if(ip0==PdgTable.end()) p0=std::to_string(e.pid0);
    else p0=ip0->second;

    if(ip1==PdgTable.end()) p1=std::to_string(e.pid1);
    else p1=ip1->second;
    return {m0,m1,p0,p1};
  };

  auto make_process = [&make_names](const event_t & e) {
    auto [m0,m1,p0,p1] = make_names(e);
    return m0+m1+" → "+p0+p1;
  };

  double total_identified  = 0;
  double total_unknown  = 0;
  double total;
  for(auto & [n,e] :rm) { total_identified+=n; }
  for(auto & [n,e] :rm_unknown) { total_unknown+=n; }
  total = total_unknown + total_identified;
  int width=80;
  auto print_double_line = [&width] () {
    for(int i=0;i<width; ++i) std::cout << "=";
    std::cout << "\n";
  };
  auto print_line = [&width] () {
    for(int i=0;i<width; ++i) std::cout << "-";
    std::cout << "\n";
  };

  print_double_line();
  auto print = [&](auto & rm, double ntotal=1, double ntotal_identified=1,   std::string title = "") {
    char buf[1024];
    if(rm.empty()) {
      print_double_line();
      sprintf(buf, "%-30s %6s  %10s  %10s \n", "channel", "N", "% of total", "% of idtfy");
      std::cout << buf;
      width = strlen(buf);
      print_line();
      return;
    }

    //calculate and print digest for topolist 
    long sum = 0;
    for(auto & [n,e] : rm) sum+=n; 
    std::string fmt = "%-30s %6ld  %10.4f%% %10.4f%%\n";
    title = adjast_utf_string("Total " + title,30);
    sprintf(buf,fmt.c_str(), title.c_str(), sum, sum*100.0/ntotal, sum*100.0/ntotal_identified);
    std::cout << buf;
    print_line();

    //print stats for each channel
    for(auto & [n,e] : rm) {
      auto process = make_process(e);
      process = adjast_utf_string(process,30);
      double p = double(n)/ntotal*100.0;
      double pi = double(n)/ntotal_identified*100.0;
      sprintf(buf, fmt.c_str(), process.c_str(), n, p, pi);
      std::cout << buf;
//      std::cout << "         "  <<  e.mpid0 << " " << e.mpid1 << " " << e.pid0 << " " << e.pid1 << " " << e.q0 << " " << e.q1 << "   :   " << n;
      //std::cout << endl;
    }
    print_double_line();
  };
  char buf[1024];
  //print(rm, total, total_identified);
  sprintf(buf, "%-20s %20.0f\n","Total events:",  total);
  std::cout << buf;
  print(head);
  //print(rm_tt,total,total_identified,"τ⁻τ⁺");
  print(rm_tt_eX,total,total_identified,"τ⁻τ⁺ → eX");
  print(rm_tt_not_eX,total,total_identified,"τ⁻τ⁺ → xx");
  print(rm_not_tt, total,total_identified, "not τ⁻τ⁺");
  print(rm_unknown, total, total_identified, "unknown");
}

void cmpall(Scan_t DATA, Scan_t SIGNAL, Selection_t &S) {
  //cmp({SIGNAL , DATA} , S , "ptem"           , "" , "NORM" , 40 , 0   , 1.1);
  //cmp({SIGNAL , DATA} , S , "cos_theta_mis2" , "" , "NORM" , 40 , -1  , +1);
  //cmp({SIGNAL , DATA} , S , "acop"           , "" , "NORM" , 40 , 0   , TMath::Pi());
  //cmp({SIGNAL , DATA} , S , "acol"           , "" , "NORM" , 40 , 0   , 1);
  //cmp({SIGNAL , DATA} , S , "p"              , "" , "NORM" , 40 , 0   , 1.1);
  //cmp({SIGNAL , DATA} , S , "Emis"           , "" , "NORM" , 40 , 1.0 , 3.5);
  //cmp({SIGNAL , DATA} , S , "cos(theta)"     , "" , "NORM" , 40 , -1  , -1);
  cmp({SIGNAL , DATA} , S , "Mpi0"     , "" , "NORM" , 40 , 0.11  , 0.15);
  cmp({SIGNAL , DATA} , S , "Mrho"     , "" , "NORM" , 40 , 0.5  , 1.1);
}
