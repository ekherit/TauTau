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

struct Analysis 
{
  const std::string WORKDIR="./";
  const std::string SHAREDIR="/home/nikolaev/tauscan/TauTau/share/";
  double EMS_ENERGY_SHIFT = -0.0201182;
  ibn::valer<double> TAU_SPREAD = {1.2299, 0.034};

  std::string TAUFIT = "taufit --minos --pdgshift --lum=default  --free-energy --free-luminosity --free-effcor --free-spread --draw-tau-mass-precision=4  --draw-diff ";

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
  std::vector<Selection_t> SEL;

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

  Analysis(
      std::string workdir, 
      std::string rtname,  
      const std::vector<ParticleID_t> & pid, 
      std::string common_cut,
      const  std::vector< std::tuple<std::string, std::string > > & sel
      ) : 
    WORKDIR{workdir},
    runtable_name{SHAREDIR+rtname},
    PID(pid),
    SEL(make_selections(pid,common_cut, sel) ) 
    {
      Init();
    };


  void Init(void);  

  void measure_tau_luminosity(void) {
    measure_luminosity(DATA,BB,GG,1e6);
    set_luminosity(DATA,SIGNAL);
    for(auto d : BGs) set_luminosity(DATA,d);
    set_gg_luminosity(DATA,BB);
  };

  void set_pid(std::vector<ScanPoint_t> & D, const std::vector<ParticleID_t> & Pid);

  template<typename Projector> void measure_efficiency(Scan_t & scan, long N0_MC, Projector  proj, std::string sel="");

  template<typename Projector> void measure_luminosity        ( Scan_t & data , Scan_t & mc          , long N0_MC      , std::string sel   , Projector proj);

  void measure_bhabha_luminosity ( Scan_t & data , Scan_t & mc          , long N0_MC      , std::string sel);
  void measure_gg_luminosity     ( Scan_t & data , Scan_t & mc          , long N0_MC      , std::string sel);
  void measure_luminosity        ( Scan_t & data , Scan_t & bb          , Scan_t & gg     , long N0_MC);
  void set_luminosity            ( Scan_t & DATA , Scan_t & MC );
  void set_luminosity            ( Scan_t & DATA , std::map<std::string , Scan_t>  & G );
  void set_gg_luminosity         ( Scan_t & DATA , Scan_t & MC );

  //Apply selection cut "sel" to ProjSampe "sample" (tt,bb,gg) and return new Scan_t
  template<typename ProjSample> auto select(const Scan_t & P, const std::string & sel,  ProjSample sample) const ->  Scan_t;
  template<typename ProjSample> auto select(const Scan_t & P, const Selection_t & S, ProjSample sample, std::string extra_cut="") const -> SelectionResult_t;
  template<typename ProjSample> auto select(const Scan_t & P, const std::vector<Selection_t> & S, ProjSample sample, std::string extra_cut="") const -> std::vector<SelectionResult_t>;

  auto select_tt(const Scan_t & P, const std::string & sel) const ->  Scan_t;
  auto select_bb(const Scan_t & P, const std::string & sel) const ->  Scan_t;
  auto select_gg(const Scan_t & P, const std::string & sel) const ->  Scan_t;


  auto select_tt(const Scan_t & P, const Selection_t & S, std::string extra_cut="")const  -> SelectionResult_t;
  auto select_bb(const Scan_t & P, const Selection_t & S, std::string extra_cut="")const  -> SelectionResult_t;
  auto select_gg(const Scan_t & P, const Selection_t & S, std::string extra_cut="")const  -> SelectionResult_t;

  auto select_tt(const Scan_t & P, const std::vector<Selection_t> & S, std::string extra_cut="")const  -> std::vector<SelectionResult_t>;
  auto select_bb(const Scan_t & P, const std::vector<Selection_t> & S, std::string extra_cut="")const  -> std::vector<SelectionResult_t>;
  auto select_gg(const Scan_t & P, const std::vector<Selection_t> & S, std::string extra_cut="")const  -> std::vector<SelectionResult_t>;

  auto select_tt(std::string extra_cut="") const -> std::vector<SelectionResult_t> { return select_tt(DATA, SEL, extra_cut); }

  template<typename Proj> void set_effcor    (Scan_t  & P, Proj proj, double Wref);
  template<typename Proj> void set_efficiency(Scan_t  & psr,  const Scan_t & mc, Proj proj);
  template<typename Proj> void set_effcor    (std::vector<SelectionResult_t>  & P, Proj proj, double Wref);
  template<typename Proj> void set_efficiency(std::vector<SelectionResult_t>  & D /*data*/,  const std::vector<SelectionResult_t> & M /*Monte Carlo*/, Proj proj);

  inline void do_tau(std::string extra_cut = "") {
    if(!is_tau_luminosity_measured)  {
      std::clog << "Measure tau luminosity...\n";
      measure_tau_luminosity();
      is_tau_luminosity_measured = true;
    } else {
      std::clog << "Tau luminosity already measured before\n";
    };
    print_luminosity(DATA);
    std::clog << "Selecting tau tau events for data...\n";
    auto result = select_tt(DATA,SEL,extra_cut);
    std::clog << "Selecting tau tau events for MC...\n";
    auto mc   = select_tt(SIGNAL,SEL,extra_cut);
    set_effcor(mc, &ScanPoint_t::tt, 2*MTAU-1.0*MeV);
    set_efficiency(result,mc, &ScanPoint_t::tt);
    print_efficiency(result);
    print_effcor(result);
    save(fold(result), "test.txt");
    fit("test.txt","test_title");
  };

  void save(const SelectionResult_t & sr, std::string  filename="scan.txt", std::string default_lum="gg")  const;
  void fit(std::string  filename, std::string title="", std::string default_lum = "", bool wait = false) const;

  private:
  bool is_tau_luminosity_measured = false;
};

inline void Analysis::Init(void) {
  std::cout << "Apply common cuts to SEL\n";
  //SEL.apply_common_cut();
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
};

inline void Analysis::set_pid(std::vector<ScanPoint_t> & D, const std::vector<ParticleID_t> & Pid) 
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
inline void Analysis::measure_efficiency(Scan_t & scan, long N0_MC, Projector  proj, std::string sel) {
  for(auto & sp : scan) {
    auto & ds = proj(sp);
    ds.N0mc = N0_MC;
    ds.Nmc = ds.tree->GetEntries(sel.c_str());
    ds.efficiency.value = double(ds.Nmc)/ds.N0mc;
    ds.efficiency.error = sqrt( ds.efficiency.value * ( 1.0 - ds.efficiency.value )/ds.N0mc );
  }
}

template<typename Projector>
inline void Analysis::measure_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel, Projector proj) {
  measure_efficiency(mc, N0_MC, proj, sel);
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

inline void Analysis::measure_bhabha_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel) {
  measure_luminosity(data,mc,N0_MC,sel, [](ScanPoint_t &sp) -> DataSample_t & { return sp.bb; } );
}

inline void Analysis::measure_gg_luminosity(Scan_t & data, Scan_t & mc, long N0_MC, std::string sel) {
  measure_luminosity(data,mc,N0_MC,sel, [](ScanPoint_t &sp) -> DataSample_t & { return sp.gg; } );
}


inline  void Analysis::measure_luminosity(Scan_t & data, Scan_t & bb, Scan_t & gg, long N0_MC) {
    measure_gg_luminosity(data,gg,N0_MC, GG_SEL);
    measure_bhabha_luminosity(data,bb,N0_MC, BB_SEL);
  }

inline  void Analysis::set_luminosity(Scan_t & DATA, Scan_t & MC ) 
{
  for(auto & mc : MC ) {
    auto & p  = *std::min_element(DATA.begin(), DATA.end(), [&mc](auto a, auto b){ return fabs(a.energy-mc.energy)<fabs(b.energy-mc.energy); } );
    mc.bb = p.bb;
    mc.gg = p.gg;
  }
}

inline  void Analysis::set_luminosity(Scan_t & DATA, std::map<std::string, Scan_t>  & G )
{
  for(auto & [name, s] : G ) {
    set_luminosity(DATA,s);
  }
}

inline  void Analysis::set_gg_luminosity(Scan_t & DATA, Scan_t & MC )
{
  for(auto & mc : MC ) {
    auto & p  = *std::min_element(DATA.begin(), DATA.end(), [&mc](auto a, auto b){ return fabs(a.energy-mc.energy)<fabs(b.energy-mc.energy); } );
    mc.gg = p.gg;
  }
}


template<typename Sample>
inline auto Analysis::select(const Scan_t & P, const std::string & sel, Sample sample) const ->  Scan_t
{
  std::vector<ScanPoint_t> R;
  for(const auto & p : P) {
    R.push_back(p.select(sample, sel));
  };
  return R;
}

inline auto Analysis::select_tt(const Scan_t & P, const std::string & sel) const ->  Scan_t {
  return select(P, sel, &ScanPoint_t::tt);
}

inline auto Analysis::select_bb(const Scan_t & P, const std::string & sel) const ->  Scan_t {
  return select(P, sel, &ScanPoint_t::tt);
}

inline auto Analysis::select_gg(const Scan_t & P, const std::string & sel) const ->  Scan_t {
  return select(P, sel, &ScanPoint_t::tt);
}

template<typename Proj>
auto Analysis::select(const Scan_t & P, const Selection_t & S, Proj sample, std::string extra_cut) const -> SelectionResult_t 
{
  //std::cout << "Analysis::select_tt(Scan_t, Selection_t, Proj, std::string)" << std::endl;
  SelectionResult_t r;
  r = S; //copy selection config
  //Scan_t sps   //std::cout << "after sps" << std::endl;
  r = select(P, S.cut + extra_cut, sample);
  //std::cout << "after r = sps" << std::endl;
  return r;
};


inline auto Analysis::select_tt(const Scan_t & P, const Selection_t & S, std::string extra_cut) const -> SelectionResult_t 
{
  //std::cout << "Analysis::select_tt(Scan_t, Selection_t, std::string)" << std::endl;
  return select(P,S,&ScanPoint_t::tt);
};

template<typename Proj>
inline std::vector<SelectionResult_t> Analysis::select(const Scan_t & P, const std::vector<Selection_t> & S, Proj sample, std::string extra_cut) const {
  std::vector<SelectionResult_t> R;
  int i=0;
  for(const Selection_t & s : S) {
    SelectionResult_t r = select(P, s, sample, extra_cut);
    print(r,i);
    R.push_back(r);
    ++i;
  }
  print(fold(R,"all"),-1);
  return R;
};

inline auto Analysis::select_tt(const Scan_t & P, const std::vector<Selection_t> & S, std::string extra_cut)const  -> std::vector<SelectionResult_t> {
  return select(P,S,&ScanPoint_t::tt,extra_cut);
}


template<typename Proj>
void  Analysis::set_effcor(Scan_t  & P, Proj proj, double Wref)  {
  //determine reference point
  std::cout << "Find reference point\n";
  DataSample_t  & ds0 = std::invoke(proj, *find_minimum(P, [&Wref](const auto & p) { std::cout << " p.energy = " << p.energy << "   Wref = " << Wref << std::endl; return  std::abs(p.energy - Wref); } ) );
  std::cout << "End of finding reference point \n";
  for(auto & sp : P) {
    DataSample_t & ds = std::invoke(proj, sp);
    ds.effcor.value = ds.efficiency/ds0.efficiency.value;
    ds.effcor.error = ds.efficiency.error/ds0.efficiency.value;  //we should not take into account the efficiency error of the reference point
    //, ds0.efficiency.error / ds0.efficiency.value);
  };
  ds0.effcor.error = 0;
};

  /* Теперь надо присвоить эффективность регистрации для данных, а также
   * поправки к эффективности регистрации 
   * Для этого надо для точки по энергии найти ближайшую точку в моделировании
   * */
template< typename Proj> 
void Analysis::set_efficiency(Scan_t  & psr,  const Scan_t & mc, Proj proj)
{
  for(size_t i = 0; i< psr.size(); ++i) {
    ScanPoint_t & dsp = psr[i]; //data scan point
    const ScanPoint_t & msp = *find_minimum(mc, [&dsp](const ScanPoint_t & p)->double { return std::abs(p.energy.value-dsp.energy.value);});
    DataSample_t & d = std::invoke(proj, dsp);
    const DataSample_t & m = std::invoke(proj, msp);
    d.efficiency = m.efficiency;
    d.effcor     = m.effcor;
  }
}


template<typename Proj>
inline void  Analysis::set_effcor(std::vector<SelectionResult_t>  & P, Proj proj, double Wref)  {
  for(size_t i=0;i!=P.size();++i) {
    set_effcor(P[i], proj, Wref);
  }
};

template< typename Proj> 
void Analysis::set_efficiency(std::vector<SelectionResult_t> & D /*data*/,  const std::vector<SelectionResult_t> & M /*Monte Carlo*/, Proj proj)
{
  assert(D.size() == M.size());
  for(size_t i=0;i!=D.size();++i) {
    set_efficiency(D[i], M[i], proj);
  }
}



inline void Analysis::save(const SelectionResult_t & sr, std::string  filename, std::string default_lum) const
{
  std::cout << "Saving selection: " << sr.title << " to file: " << filename << std::endl;
  std::stringstream os;
  os << ibn::mformat("15.6", "#", "L,nb^-1","dL,nb^-1","W,MeV", "dW,MeV","Sw,MeV","dSw,MeV","Ntt","Nbb","Ngg","effcor");
  os << "\n";
  auto lum = [&default_lum](const ScanPoint_t & sp) -> ibn::valer<double> {
    if(default_lum == "bb" || default_lum=="ee") { return sp.bb.luminosity; }
    if(default_lum == "gg") { return sp.gg.luminosity; }
    return sp.luminosity;//*(nb/pb);
  };

  int idx=0;
  for(const auto & p : sr) {
    std::string name = (p.name == "" ? std::to_string(++idx) : p.name);
    os << ibn::mformat("15.6", 
        name.c_str(),
        lum(p).value, lum(p).error,
        p.energy.value/MeV,  p.energy.error/MeV,
        1.258,             0.060,
        p.tt.N,
        p.bb.N,        
        p.gg.N,
        p.tt.effcor.value, p.tt.effcor.error);
    os << "\n";
  }
  std::cout << os.str();
  std::ofstream ofs(filename);
  ofs << os.str();
}


inline void Analysis::fit(std::string  filename, std::string title, std::string default_lum, bool wait) const
{

  std::string command= ibn::format("%s --tau-spread=%.6f --tau-spread-error=%.6f --ems-cmenergy-shift=%.6f  %s &", TAUFIT.c_str(), TAU_SPREAD.value, TAU_SPREAD.error, EMS_ENERGY_SHIFT, filename.c_str()); 
  std::cout << "Doing fit:\n";
  std::cout << command << "\n";
  system(command.c_str());
}

