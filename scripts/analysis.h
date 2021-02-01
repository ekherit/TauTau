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
#include <TFitResult.h>

#include "../ibn/valer.h"
#include "../ibn/indexer.h"

#include "PhysConst.h"

#include "Selection.h"
#include "draw_helper.h"
#include "fold.h"
#include "ScanPoint.h"
#include "RunTable.h"
#include "Read.h"
#include "Print.h"

#include "utils.h"

#include "pdg_table.h"



//std::string TAUFIT = "taufit --lum=bes --tau-spread=1.258 --energy-correction=-0.078";
static std::string TAUFIT = "taufit --lum=bes --tau-spread=1.258 --energy-correction=+0.011";
//static std::string PSIFIT = "psifit --free-lum --free-effcor --pdg-shift-file=shift.txt ";
static std::string PSIFIT = "psifit --free-lum --free-effcor --free-mpdg --minos ";


void clean_taufit(void) {
  system("killall -9 taufit");
}

static std::string BB_SEL = "(acol-TMath::Pi())>-0.03 && abs(cos(theta[0]) < 0.8 && abs(cos(theta[1])) < 0.8";
static std::string GG_SEL = "";


std::map<std::string, std::string > SelMap;


struct Analysis {
  const std::string WORKDIR="./";
  const std::string SHAREDIR="/home/nikolaev/tauscan/TauTau/share/";
  std::string TAUFIT = "taufit --minos --pdgshift --lum=default --tau-spread=1.2299 --ems-cmenergy-shift=-0.0201182 --free-energy --free-luminosity --free-effcor --draw-tau-mass-precision=4  --draw-diff";

  std::string runtable_name = SHAREDIR+"/all_scan_points_ems32.txt";

  Scan_t ALLRUNTABLE      = read_my_runtable(runtable_name);
  Scan_t RUNTABLE         = read_my_runtable(runtable_name, DATA_TAU);
  Scan_t JPSI_RUNTABLE    = read_my_runtable(runtable_name, DATA_JPSI);
  Scan_t PSIP_RUNTABLE    = read_my_runtable(runtable_name, DATA_PSIP);
  Scan_t DATA             = read_data(WORKDIR+"data", RUNTABLE);
  Scan_t JPSI             = read_mh(WORKDIR+"mhdata",JPSI_RUNTABLE);
  Scan_t PSIP             = read_mh(WORKDIR+"mhdata",PSIP_RUNTABLE);
  Scan_t JPSIMC           = read_mc_mh(WORKDIR+"mcmh", JPSI_RUNTABLE);
  Scan_t PSIPMC           = read_mc_mh(WORKDIR+"mcmh", PSIP_RUNTABLE);
  Scan_t JPSIGG           = read_data(WORKDIR+"data", JPSI_RUNTABLE);
  Scan_t JPSIGGMC         = read_mc(WORKDIR+"mc/gg", JPSI_RUNTABLE,1e5);
  Scan_t PSIGG            = read_data(WORKDIR+"data", PSIP_RUNTABLE);
  Scan_t PSIGGMC          = read_mc(WORKDIR+"mc/gg", PSIP_RUNTABLE,1e5);
  Scan_t SIGNAL           = read_mc(WORKDIR+"mc/signal", RUNTABLE, 1e6);

  //two-gamma background
  Scan_t EEEE       = read_mc(WORKDIR+"mc/galuga/ee"   , RUNTABLE, 1e6);
  Scan_t EEUU       = read_mc(WORKDIR+"mc/galuga/uu"   , RUNTABLE, 1e6);
  Scan_t EEPIPI     = read_mc(WORKDIR+"mc/galuga/pipi"   , RUNTABLE, 1e6);
  Scan_t EEKK       = read_mc(WORKDIR+"mc/galuga/KK"   , RUNTABLE, 1e6);

  std::map<std::string, ScanRef_t> GALUGA;

  //hadronic background
  Scan_t HADR       = read_mc(WORKDIR+"mc/hadrons", RUNTABLE, 1e6);
  //bhabha background
  Scan_t BB         = read_mc(WORKDIR+"mc/bb", RUNTABLE, 1e6);
  //uu background
  Scan_t UU         = read_mc(WORKDIR+"mc/uu", RUNTABLE, 1e6);
  //pipi background (unused)
  Scan_t PIPI       = read_mc(WORKDIR+"mc/pipi", RUNTABLE, 1e6);
  // for luminocity measurement
  Scan_t GG         = read_mc(WORKDIR+"mc/gg", RUNTABLE, 1e6);

  //list of all backgrounds
  std::vector<ScanRef_t> BGs ={HADR, BB, UU, GG, EEEE, EEUU, EEPIPI, EEKK};

  struct Simulation_t {
    ScanRef_t signal;
    std::vector<ScanRef_t> bgs;
  };

  Simulation_t MC = { SIGNAL, BGs }; 

  std::string BB_SEL = "(acol-TMath::Pi())>-0.04 && abs(cos(theta[0])) < 0.8 && abs(cos(theta[1])) < 0.8 && Ep[0]>0.8 && Ep[1]>0.8 && abs(z[0])<10 && abs(z[1])<10 && vxy[0]<1.0 && vxy[1]<1.0";
  std::string GG_SEL = "";
  std::string MH_SEL = "ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2";
  std::vector<ParticleID_t> PID;
  /*
  = 
  {
    {
      {"e", {"Ep[#]>0.80"}, 
        { 
          //restrict("Ep",         0.8,  1.1), 
          restrict("chi2_dedx",   0 ,  3.0),  
          restrict("delta_tof", -0.3,  0.3) 
        }}, 
      {"u", { "!e#" 
              ,"depth[#]-p[#]*58.61 > -35"
                ,"Nmuhit[#]>=2"
            },
      { 
        {"depth",        0, 100},
        restrict("E",          0.1, 0.3), 
        restrict("Ep",           0, 0.7), 
        restrict("chi2_dedx",    0, 3.0),  
        restrict("delta_tof", -0.3, 0.3) 
      }},
      {"pi", { "!(e#||u#)" },
        {
          restrict("Ep",           0, 0.7), 
          restrict("chi2_dedx",    0, 3.0),  
          restrict("delta_tof", -0.3, 0.3)
            ,restrict("p", 0.73, 1.1) 
        }}, 
      {"PI", { "!(e#||u#)" }, 
        {
          restrict("Ep",           0, 0.7), 
          restrict("chi2_dedx",    0, 3.0),  
          restrict("delta_tof", -0.3, 0.3)
        }}, 
      {"K", { "!(e#||u#||pi#)" },
        { 
          restrict("Ep",           0, 0.6), 
          restrict("chi2_dedx",    0, 3.0),  
          restrict("delta_tof", -0.3, 0.3)
        }},
      {"rho", {" PI# && 0.12 < Mpi0[0] && Mpi0[0] < 0.14"},
        { 
          restrict("Mrho", 0.6,1.0),
        }}, 

      {"X", { "chi2_dedx_e[#]>3" 
            },
      { 
        restrict("Ep",           0, 0.8), 
      }},
    }
  };
  */

  Selection_t SEL;
  /*
  =
  {
    "sel12", //selection name
    "Ncg==2"
      "&& Ncg==Ncc"
      "&& abs(vz[0])<10 && abs(vz[1])<10"
      "&& vxy[0]<1.0 && vxy[1]<1.0"
      "&& E[0]>0.025 && E[1]>0.025"
      "&& q[0]==-1 && q[1]==1"
      "&& p[0] < 1.1 && p[1] < 1.1"
      "&& pt[0] > 0.2 && pt[1] > 0.2"
      "&& cgood"
      "&& 2.5 < tof[0] && tof[0] < 5.5 && 2.5 < tof[1] && tof[1] < 5.5"
      "&& ptem50>0.25 && ptem50<1.1"
      , 
    PID,
    { 
      {"eX",      "NnE50==0 && eX && missed_photon_angle"},
      {"eρ",     "Nng==2   && eX && Mpi0[0] < 0.14 && Mpi0[0]>0.12 && good_emc_time" },
    },
  };
  */

  Analysis(void) {
    std::cout << "In default constructor\n";
    Init();
  };

  Analysis(std::string workdir, std::string rtname) : 
    WORKDIR{workdir}, 
    runtable_name(SHAREDIR+rtname) 
  {
    std::cout << "In rtname\n";
    Init();
  };

  Analysis(std::string workdir, std::string rtname, const  Selection_t & sel) : 
    WORKDIR{workdir},
    runtable_name{SHAREDIR+rtname}, 
    SEL(sel) 
  {
    std::cout << "In rtname\n";
    Init();
  };


  /*
  Analysis(const Analysis & an) {
    std::cout << "In constructor\n";
    *this = an;
    Init();
  };
  */

  void Init(void) {
    std::cout << "Apply common cuts to SEL\n";
    SEL.apply_common_cut();
    std::vector<ScanRef_t> LUM_MCs = {BB,GG};
    //std::vector<ScanRef_t> BG_MCs =  {HADR, UU, PIPI};
    //std::vector<ScanRef_t> BGall_MCs =  BG_MCs;

    std::map<std::string, ScanRef_t> GALUGA = {
      {"ee", EEEE},
      {"uu",  EEUU},
      {"pipi", EEPIPI},
      {"KK", EEKK},
    };

    //for( auto & p: GALUGA) BGall_MCs.push_back(p.second);

    read_tau_cross_section(SHAREDIR+"/tau_cross_section.txt", SIGNAL);
    read_bhabha_cross_section(SHAREDIR+"/bhabha_cross_section.txt", BB);
    read_gg_cross_section(SHAREDIR+"/gg_cross_section.txt", GG);
    read_galuga_cross_section(SHAREDIR+"/galuga_cross_section.txt", GALUGA);
    read_mumu_or_pipi_cross_section(SHAREDIR+"/mumu_cross_section.txt", UU);
    read_pipi_cross_section(SHAREDIR+"/pipi_cross_section.txt", PIPI);
    read_hadron_cross_section(SHAREDIR+"/hadron_cross_section.txt", HADR);

    std::cout << "Setting  PID" << std::endl;
    set_pid(DATA       , PID);
    set_pid(SIGNAL         , PID);
    for(auto d : BGs)       set_pid(d,PID);
    for(auto d : LUM_MCs)   set_pid(d,PID);

    std::cout << "Measuring luminosity" << std::endl;
    measure_luminosity(DATA,BB,GG,1e6);
    set_luminosity(DATA,SIGNAL);
    for(auto d : BGs) set_luminosity(DATA,d);
    set_gg_luminosity(DATA,BB);
  };

  void measure_luminosity(void) {

  };

  void set_pid(std::vector<ScanPoint_t> & D, const std::vector<ParticleID_t> & Pid);

  template<typename Projector> 
  void measure_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel, Projector proj);
  void measure_bhabha_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel); 
  void measure_gg_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel); 

  void measure_luminosity(Scan_t & data, Scan_t & bb, Scan_t & gg, long N0_MC) {
    measure_gg_luminosity(data,gg,N0_MC, GG_SEL);
    measure_bhabha_luminosity(data,bb,N0_MC, BB_SEL);
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

  void set_gg_luminosity(Scan_t & DATA, Scan_t & MC ) {
    for(auto & mc : MC ) {
      auto & p  = *std::min_element(DATA.begin(), DATA.end(), [&mc](auto a, auto b){ return fabs(a.energy-mc.energy)<fabs(b.energy-mc.energy); } );
      mc.gg = p.gg;
    }
  }

  //Apply selection cut "sel" to ProjSampe "sample" (tt,bb,gg) and return new Scan_t
  template<typename ProjSample> 
  auto select(const Scan_t & P, const std::string & sel,  ProjSample sample) const ->  Scan_t;

  template<typename ProjSample> 
  auto select(const Scan_t & P, const ChannelSelection_t & S, ProjSample sample, std::string extra_cut="") const -> ChannelSelectionResult_t;

  template<typename ProjSample> 
  auto select(const Scan_t & P, const Selection_t & S, ProjSample sample, std::string extra_cut="") const -> std::vector<ChannelSelectionResult_t>;

  auto select_tt(const Scan_t & P, const std::string & sel) const ->  Scan_t;
  auto select_bb(const Scan_t & P, const std::string & sel) const ->  Scan_t;
  auto select_gg(const Scan_t & P, const std::string & sel) const ->  Scan_t;


  auto select_tt(const Scan_t & P, const ChannelSelection_t & S, std::string extra_cut="")const  -> ChannelSelectionResult_t;
  auto select_bb(const Scan_t & P, const ChannelSelection_t & S, std::string extra_cut="")const  -> ChannelSelectionResult_t;
  auto select_gg(const Scan_t & P, const ChannelSelection_t & S, std::string extra_cut="")const  -> ChannelSelectionResult_t;

  auto select_tt(const Scan_t & P,        const Selection_t & S, std::string extra_cut="")const  -> std::vector<ChannelSelectionResult_t>;
  auto select_bb(const Scan_t & P,        const Selection_t & S, std::string extra_cut="")const  -> std::vector<ChannelSelectionResult_t>;
  auto select_gg(const Scan_t & P,        const Selection_t & S, std::string extra_cut="")const  -> std::vector<ChannelSelectionResult_t>;

  auto select_tt(std::string extra_cut="") const -> std::vector<ChannelSelectionResult_t> {
    return select_tt(DATA, SEL, extra_cut);
  }

};


void Analysis::set_pid(std::vector<ScanPoint_t> & D, const std::vector<ParticleID_t> & Pid) 
{
  for( auto & data : D) {
    for(int track=0;track<2;++track) {
      for(auto & pid: Pid) {
        std::string alias_value;
        std::string name = pid.name + std::to_string(track);
        for (auto  cut: pid.cuts) {
          cut = sub(cut,R"(#)",std::to_string(track));
          //std::cout << name << " pid_cut = " <<  cut << std::endl; 
          if(alias_value=="") alias_value = cut; 
          else alias_value += " && " + cut;
        }
        for (auto & r : pid.restrictions) {
          std::string value_name = r.name;

          std::string pid_name= pid.name;
          if(pid_name == "u") pid_name = "mu";
          if(pid_name == "PI") pid_name = "pi";

          if ( value_name == "chi2_dedx"  || value_name == "delta_tof")  value_name += "_"+pid_name;
          value_name += "["+std::to_string(track)+"]";
          auto rm_zero = [](const std::string  input) -> std::string {
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

          std::string cut = "("+min_str + "<" + value_name + "&&" + value_name  + "<"  +  max_str +")";
          if(alias_value=="") alias_value = cut; 
          else alias_value += "&&" + cut;
        };
        data.tt.tree->SetAlias(name.c_str(), alias_value.c_str());
      }
    }
  };
};

template<typename Projector>
void Analysis::measure_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel, Projector proj) {
  me(mc, N0_MC, proj, sel);
  //find close point and select
  for(auto & sp :data) {
    auto & p  = *std::min_element(mc.begin(), mc.end(), [&sp](auto a, auto b){ return fabs(a.energy-sp.energy)<fabs(b.energy-sp.energy); } );
    if(fabs(sp.energy-p.energy)>1*MeV) {
      std::cerr << "WARNING: Luminosity measurement:  To big difference in energy points: " << sp.energy.value/MeV << "(data)  and " << p.energy/MeV << " MeV (MC)" << std::endl;
    }
    auto & l = proj(sp);
    auto & m = proj(p);
    l.cross_section.value  = m.cross_section.value*pow(m.energy.value/l.energy.value, 2.0);
    l.cross_section.error  = l.cross_section.value*std::hypot( m.cross_section.error/m.cross_section.value, 2.0*l.energy.error/l.energy.value);
    //std::cout << "l.sigma = " << l.cross_section.value << " +- " << l.cross_section.error << std::endl;
    l.efficiency = m.efficiency;
    l.Nmc = m.Nmc;
    l.N0mc = m.N0mc;
    l.effcor = m.effcor;
    l.N = l.tree->GetEntries(sel.c_str());
    double vis_cx = (l.efficiency*l.cross_section); //visible_cross_section;
    double vis_cx_error = std::hypot(l.cross_section.error*l.efficiency, l.cross_section*l.efficiency.error);
    //std::cout << "Visible cross section : " << vis_cx << "  " << vis_cx_error << std::endl;
    m.luminosity.value  = l.N / vis_cx;
    m.luminosity.error = m.luminosity.value*sqrt(  1./l.N  +  pow(vis_cx_error/vis_cx,2.0) );
    l.luminosity.value = m.luminosity.value;
    l.luminosity.error = m.luminosity.error;
    //std::cout << sp.energy << "  " << l.Nmc << " " << l.N0mc << " " << l.efficiency.value << "  " << l.efficiency.error << "  " << l.luminosity.value << " " << l.luminosity.error <<  "mc" << m.luminosity.error << std::endl;
  }
}

void Analysis::measure_bhabha_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel) {
  measure_luminosity(data,mc,N0_MC,sel, [](ScanPoint_t &sp) -> DataSample_t & { return sp.bb; } );
}

void Analysis::measure_gg_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel) {
  measure_luminosity(data,mc,N0_MC,sel, [](ScanPoint_t &sp) -> DataSample_t & { return sp.gg; } );
}


template<typename Sample>
auto Analysis::select(const Scan_t & P, const std::string & sel, Sample sample) const ->  Scan_t {
  std::vector<ScanPoint_t> R;
  for(const auto & p : P) {
    R.push_back(p.select(sample, sel));
  };
  return R;
}

auto Analysis::select_tt(const Scan_t & P, const std::string & sel) const ->  Scan_t {
  return select(P, sel, &ScanPoint_t::tt);
}

auto Analysis::select_bb(const Scan_t & P, const std::string & sel) const ->  Scan_t {
  return select(P, sel, &ScanPoint_t::tt);
}

auto Analysis::select_gg(const Scan_t & P, const std::string & sel) const ->  Scan_t {
  return select(P, sel, &ScanPoint_t::tt);
}

template<typename Proj>
auto Analysis::select(const Scan_t & P, const ChannelSelection_t & S, Proj sample, std::string extra_cut) const -> ChannelSelectionResult_t 
{
  ChannelSelectionResult_t r;
  r = S; //copy selection config
  r  = select(P, S.cut + extra_cut, sample);
  return r;
};


auto Analysis::select_tt(const Scan_t & P, const ChannelSelection_t & S, std::string extra_cut) const -> ChannelSelectionResult_t 
{
  return select(P,S,&ScanPoint_t::tt);
};

template<typename Proj>
std::vector<ChannelSelectionResult_t> Analysis::select(const Scan_t & P, const Selection_t & S, Proj sample, std::string extra_cut) const {
  std::vector<ChannelSelectionResult_t> R;
  int i=0;
  for(const ChannelSelection_t & s : S) {
    ChannelSelectionResult_t r = select(P, s, sample, extra_cut);
    print(r,i);
    R.push_back(r);
    ++i;
  }
  print(fold(R,"all"),-1);
  return R;
};

auto Analysis::select_tt(const Scan_t & P,        const Selection_t & S, std::string extra_cut)const  -> std::vector<ChannelSelectionResult_t> {
  return select(P,S,&ScanPoint_t::tt,extra_cut);
}


