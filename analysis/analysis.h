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

#include <fmt/core.h>

#include "time.h"
#include "stdlib.h"

#include <TStyle.h>

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
#include "Draw.h"

#include "utils.h"

#include "pdg_table.h"

template<typename Proj> 
void copy_luminosity(Proj proj, const Scan_t & IN, Scan_t & OUT) 
{
  for(auto & out : OUT) {
    const auto & in  = *std::min_element(IN.begin(), IN.end(), [&out](const auto & in1, const auto & in2){ return fabs(in1.energy-out.energy)<fabs(in2.energy-out.energy); } );
    std::invoke(proj,out) = std::invoke(proj,in);
  };
};

struct Sample_t {
  Scan_t data;
  Scan_t mc;

  template<typename Projector>
    inline void measure_efficiency(Projector  proj, std::string sel) {
      for(auto & sp : mc) {
        auto & ds = std::invoke(proj,sp);
        ds.Nmc = ds.tree->GetEntries(sel.c_str());
        ds.efficiency.value = double(ds.Nmc)/ds.N0mc;
        ds.efficiency.error = sqrt( ds.efficiency.value * ( 1.0 - ds.efficiency.value )/ds.N0mc );
      }
    }

  template<typename Projector>
  inline void measure_luminosity(Projector proj, std::string sel) {
    measure_efficiency(proj, sel);
    //find close point and select
    for(auto & sp :data) {
      auto & p  = *std::min_element(mc.begin(), mc.end(), [&sp](auto a, auto b){ return fabs(a.energy-sp.energy)<fabs(b.energy-sp.energy); } );
      if(fabs(sp.energy-p.energy)>1*MeV) {
        std::cerr << "WARNING: Luminosity measurement:  To big difference in energy points: " << sp.energy.value/MeV << "(data)  and " << p.energy/MeV << " MeV (MC)" << std::endl;
      }
      auto & l = std::invoke(proj,sp);
      auto & m = std::invoke(proj, p);
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
    copy_luminosity(proj, data,mc);
  }
};

struct Process_t {
  Sample_t signal; //signal 
  Sample_t gg; //gamma-gamma luminosity
  Sample_t bb; //bhabha luminosity
  std::vector<ScanRef_t> bgs; //Monte Carlo for backgrounds
};

struct Analysis 
{
  const std::string WORKDIR="./";
  const std::string SHAREDIR="/home/nikolaev/tauscan/TauTau/share/";
  // double EMS_ENERGY_SHIFT = -0.0201182;
  // ibn::valer<double> TAU_SPREAD = {1.2299, 0.034};
  double EMS_ENERGY_SHIFT = -0.02138719;
  ibn::valer<double> TAU_SPREAD = {1.22760117, 0.034};

  std::string TAUFIT = "taufit --minos --pdgshift --lum=default  --free-energy --free-luminosity --free-effcor --free-spread --draw-tau-mass-precision=4  --draw-diff ";

  std::string runtable_name = SHAREDIR+"/all_scan_points_ems32.txt";

  Scan_t  TAU_RUNTABLE     = read_my_runtable(runtable_name, DATA_TAU);
  Scan_t JPSI_RUNTABLE    = read_my_runtable(runtable_name,  DATA_JPSI);
  Scan_t PSIP_RUNTABLE    = read_my_runtable(runtable_name,  DATA_PSIP);

  //two-gamma background
  Scan_t EEEE       = read_mc(WORKDIR+"mc/galuga/ee"   , TAU_RUNTABLE, 1e6);
  Scan_t EEUU       = read_mc(WORKDIR+"mc/galuga/uu"   , TAU_RUNTABLE, 1e6);
  Scan_t EEPIPI     = read_mc(WORKDIR+"mc/galuga/pipi"   , TAU_RUNTABLE, 1e6);
  Scan_t EEKK       = read_mc(WORKDIR+"mc/galuga/KK"   , TAU_RUNTABLE, 1e6);

  //std::map<std::string, ScanRef_t> GALUGA;

  //hadronic background
  Scan_t HADR       = read_mc(WORKDIR+"mc/hadrons", TAU_RUNTABLE, 1e6);
  //bhabha background
  Scan_t BB         = read_mc(WORKDIR+"mc/bb", TAU_RUNTABLE, 1e6);
  //uu background
  Scan_t UU         = read_mc(WORKDIR+"mc/uu", TAU_RUNTABLE, 1e6);
  //pipi background (unused)
  Scan_t PIPI       = read_mc(WORKDIR+"mc/pipi", TAU_RUNTABLE, 1e6);
  // for luminocity measurement
  Scan_t GG         = read_mc(WORKDIR+"mc/gg", TAU_RUNTABLE, 1e6);

  //list of all backgrounds
  std::vector<ScanRef_t> BGs ={HADR, BB, UU, GG, EEEE, EEUU, EEPIPI, EEKK};


  Process_t TAU = {
    .signal = {
      .data =  read_data(WORKDIR+"data/tt", TAU_RUNTABLE),
      .mc   =  read_mc(WORKDIR+"mc/tt", TAU_RUNTABLE, 1e6),
    },
    .gg =  {
      .data =  read_data(WORKDIR+"data/tt", TAU_RUNTABLE),
      .mc   =  read_mc(WORKDIR+"mc/gg", TAU_RUNTABLE, 1e6)
    },
    .bb =  {
      .data =  read_data(WORKDIR+"data/tt", TAU_RUNTABLE),
      .mc   =  read_mc(WORKDIR+"mc/bb", TAU_RUNTABLE, 1e6)
    },
    .bgs = {HADR, BB, UU, GG, EEEE, EEUU, EEPIPI, EEKK}
  };

  Process_t JPSI = {
    .signal = {
      .data =  read_mh(WORKDIR+"data/mh",JPSI_RUNTABLE),
      .mc   =  read_mc_mh(WORKDIR+"mc/mh", JPSI_RUNTABLE),
    },
    .gg =  {
      .data =  read_data(WORKDIR+"data/tt", JPSI_RUNTABLE),
      .mc   =  read_mc(WORKDIR+"mc/gg", JPSI_RUNTABLE,1e5)
    },
  };

  Process_t PSI = {
    .signal = {
      .data =  read_mh(WORKDIR+"data/mh",PSIP_RUNTABLE),
      .mc   =  read_mc_mh(WORKDIR+"mc/mh", PSIP_RUNTABLE),
    },
    .gg =  {
      .data =  read_data(WORKDIR+"data/tt", PSIP_RUNTABLE),
      .mc   =  read_mc(WORKDIR+"mc/gg", PSIP_RUNTABLE,1e5)
    },
  };

  std::string GG_SEL = "";
  std::string BB_SEL = "(acol-TMath::Pi())>-0.04 && abs(cos(theta[0])) < 0.8 && abs(cos(theta[1])) < 0.8 && Ep[0]>0.8 && Ep[1]>0.8 && abs(z[0])<10 && abs(z[1])<10 && vxy[0]<1.0 && vxy[1]<1.0";
  std::string MH_SEL = "ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2 && maxNmuhit==0";
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
      std::string gg_sel,
      std::string bb_sel,
      std::string mh_sel,
      const std::vector<ParticleID_t> & pid, 
      std::vector<std::string> common_cut,
      const  std::vector< std::tuple<std::string, std::string > > & sel
      ) : 
    WORKDIR{workdir},
    runtable_name{SHAREDIR+rtname},
    GG_SEL{gg_sel},
    BB_SEL{bb_sel},
    MH_SEL{mh_sel},
    PID(pid),
    SEL(make_selections(pid,common_cut, sel) ) 
    {
      Init();
      std::cout << "Runtable: " << std::endl;
      print(JPSI_RUNTABLE);
      print(TAU_RUNTABLE);
      print(PSIP_RUNTABLE);
      std::cout << "GG_SEL = " << GG_SEL << "\n";
      std::cout << "BB_SEL = " << BB_SEL << "\n";
      std::cout << "MH_SEL = " << MH_SEL << "\n";
      std::cout << "TAU_SEL:  \n"; 
      //this lambda remove surrounding ( ) 
      auto nake = [](std::string_view a) -> std::string_view {
          size_t b=0;
          size_t e=a.size();
          for(auto it = a.begin(); it!=a.end();++it) {
              if(*it==' ' || *it=='(') ++b;
              else break;
          }
          for(auto it = a.rbegin(); it!=a.rend();++it) {
              if(*it==' ' || *it==')') --e;
              else  break;
          }
          return a.substr(b,e-b);
      };
      auto shift = [](int N=1) {
          std::string s;
          for(int i=0;i!=N;++i) {
              s+="     ";
          };
          return s;
      };
      std::cout << shift(1) << "Common tau-tau cuts:\n";
      for(auto & s : common_cut) {
          std::cout << shift(2) << '\"'<<s << "\"\n";
      }
      std::cout << shift(1) << "Channel specific cuts:\n";
      for(auto & s: sel) {
          std::cout << shift(3) << std::get<0>(s) <<  "\n";
          for(auto item : split(std::get<1>(s))) {
          std::cout << shift(4) << "\"" << nake(item) << "\"\n";
          }
      };
    };


  void Init(void);  


  void measure_tau_luminosity(void) {
    TAU.gg.measure_luminosity(&ScanPoint_t::gg, GG_SEL);
    TAU.bb.measure_luminosity(&ScanPoint_t::bb, GG_SEL);
    ::copy_luminosity(&ScanPoint_t::gg, TAU.gg.data, TAU.signal.mc);
    ::copy_luminosity(&ScanPoint_t::gg, TAU.bb.data, TAU.signal.mc);
    for(auto d : BGs) {
      ::copy_luminosity(&ScanPoint_t::gg, TAU.gg.data, d);
      ::copy_luminosity(&ScanPoint_t::gg, TAU.bb.data, d);
    };
    ::copy_luminosity(&ScanPoint_t::gg, TAU.gg.data, TAU.signal.data);
    ::copy_luminosity(&ScanPoint_t::bb, TAU.bb.data, TAU.signal.data);
  };

  void set_pid(std::vector<ScanPoint_t> & D, const std::vector<ParticleID_t> & Pid);


  //template<typename Projector> void measure_luminosity( Scan_t & data, Scan_t & mc, std::string sel, Projector proj);

  //copy luminosity porj from  D1 to D2 
  template<typename Proj> void copy_luminosity (const Scan_t & D1, Scan_t & D2, Proj proj) const;

  template<typename Proj> void set_default_luminosity (Scan_t & D, Proj proj) const;


  //Apply selection cut "sel" to ProjSampe "sample" (tt,bb,gg) and return new Scan_t
  template<typename Proj> auto select(Proj proj, const Scan_t & P, const std::string & sel) const ->  Scan_t;
  template<typename Proj> auto select(Proj proj, const Scan_t & P, const Selection_t & S, std::string extra_cut="") const -> SelectionResult_t;
  template<typename Proj> auto select(Proj proj, const Scan_t & P, const std::vector<Selection_t> & S, std::string extra_cut="") const -> std::vector<SelectionResult_t>;

  template<typename Proj> void set_effcor    (Proj proj, Scan_t  & P, double Wref);

  template<typename Proj> void copy_efficiency(Proj proj, const Scan_t & mc, Scan_t  & psr) const;

  template<typename Proj> void set_effcor    (Proj proj, std::vector<SelectionResult_t>  & P, double Wref);
  template<typename Proj> void copy_efficiency(Proj proj, const std::vector<SelectionResult_t> & M /*Monte Carlo*/, std::vector<SelectionResult_t>  & D /*data*/) const;

  inline void do_tau(std::string extra_cut = "") 
  {
    if(!is_tau_luminosity_measured)  {
      std::clog << "Measure tau luminosity...\n";
      measure_tau_luminosity();
      is_tau_luminosity_measured = true;
    } else {
      std::clog << "Tau luminosity already measured before\n";
    };
    print_luminosity(TAU.signal.data);
    std::clog << "Selecting tau tau events for data...\n";
    auto result = select(&ScanPoint_t::tt, TAU.signal.data,SEL,extra_cut);
    std::clog << "Selecting tau tau events for MC...\n";
    auto mc   = select(&ScanPoint_t::tt, TAU.signal.mc, SEL,extra_cut);
    set_effcor( &ScanPoint_t::tt, mc, 2*MTAU-1.0*MeV);
    copy_efficiency(&ScanPoint_t::tt, mc,result );
    print_efficiency(result);
    print_effcor(result);
    save(fold(result), "test.txt");
    fit("test.txt","test_title");
  };

  void measure_luminosity(Process_t & P) {
    read_gg_cross_section("../TauTau/share/gg_cross_section.txt", P.gg.mc);
    //measure_luminosity(P.gg.data, P.gg.mc, GG_SEL, &ScanPoint_t::gg);
    P.gg.measure_luminosity(&ScanPoint_t::gg, GG_SEL);
    copy_luminosity(P.gg.data, P.signal.data, &ScanPoint_t::gg);
    copy_luminosity(P.gg.data, P.signal.mc,   &ScanPoint_t::gg);
    set_default_luminosity(P.signal.data,     &ScanPoint_t::gg);
    print(P.signal.data);
  }

  Scan_t  res(Process_t & P, std::string extra_cut = "") {
    measure_luminosity(P);
    std::cout << "Selecting multihadronic events" << std::endl;
    auto res =   select(&ScanPoint_t::tt,  P.signal.data,  MH_SEL && extra_cut);
	std::cout << "In res function : type res: " << typeid(res).name() << std::endl;
    auto mcres = select(&ScanPoint_t::tt,  P.signal.mc,    MH_SEL && extra_cut);
    print(mcres);
	//fmt::print("Efficiency:\n");
	//print_efficiency(mcres);
    set_effcor( &ScanPoint_t::tt, mcres,0.0);  //set correction to efficiency
    copy_efficiency( &ScanPoint_t::tt, mcres, res);
    print(res);
    return res;
  };

  void res(std::string suffix="test") {
    std::string jpsiname = "jpsi_"+suffix+".txt";
    std::string psiname = "psip_"+suffix+".txt";
    std::string resname = "res_"+suffix+".txt";
    auto jpsi = res(JPSI);
    auto psi  = res(PSI);
    save(jpsi,jpsiname);
    save(psi,psiname);
    system(("cat "+jpsiname + " " + psiname + " > " +  resname).c_str());
  };

  void save(const Scan_t & sr, std::string  filename="scan.txt", std::string default_lum="gg")  const;

  void fit(std::string  filename, std::string title="", std::string default_lum = "", bool wait = false) const;

  private:
  bool is_tau_luminosity_measured = false;
};

inline void Analysis::Init(void) {
  //std::cout << "Apply common cuts to SEL\n";
  std::vector<ScanRef_t> LUM_MCs = {BB,GG};

  std::map<std::string, ScanRef_t> GALUGA = {
    {"ee", EEEE},
    {"uu",  EEUU},
    {"pipi", EEPIPI},
    {"KK", EEKK},
  };

  read_tau_cross_section(SHAREDIR+"/tau_cross_section.txt", TAU.signal.mc);
  read_bhabha_cross_section(SHAREDIR+"/bhabha_cross_section.txt", TAU.bb.mc);
  read_gg_cross_section(SHAREDIR+"/gg_cross_section.txt", TAU.gg.mc);

  read_bhabha_cross_section(SHAREDIR+"/bhabha_cross_section.txt", BB);
  read_gg_cross_section(SHAREDIR+"/gg_cross_section.txt", GG);
  read_galuga_cross_section(SHAREDIR+"/galuga_cross_section.txt", GALUGA);
  read_mumu_or_pipi_cross_section(SHAREDIR+"/mumu_cross_section.txt", UU);
  read_pipi_cross_section(SHAREDIR+"/pipi_cross_section.txt", PIPI);
  read_hadron_cross_section(SHAREDIR+"/hadron_cross_section.txt", HADR);

  std::cout << "Setting particle identification aliases (PID)" << std::endl;
  set_pid(TAU.signal.data, PID);
  set_pid(TAU.signal.mc, PID);
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
 
template<typename Proj> void Analysis::copy_luminosity(const Scan_t & IN, Scan_t & OUT, Proj proj) const {
  for(auto & out : OUT) {
    const auto & in  = *std::min_element(IN.begin(), IN.end(), [&out](const auto & in1, const auto & in2){ return fabs(in1.energy-out.energy)<fabs(in2.energy-out.energy); } );
    std::invoke(proj,out) = std::invoke(proj,in);
  };
};

template<typename Proj> void Analysis::set_default_luminosity(Scan_t & D, Proj proj) const {
  for(auto & d : D) {
    d.luminosity = std::invoke(proj,d).luminosity;
  };
};



template<typename Proj>
inline auto Analysis::select(Proj proj, const Scan_t & P, const std::string & sel) const ->  Scan_t
{
  std::vector<ScanPoint_t> R;
  for(const auto & p : P) {
    R.push_back(p.select(proj, sel));
  };
  return R;
}


template<typename Proj>
auto Analysis::select(Proj proj, const Scan_t & P, const Selection_t & S, std::string extra_cut) const -> SelectionResult_t 
{
  //std::cout << "Analysis::select_tt(Scan_t, Selection_t, Proj, std::string)" << std::endl;
  SelectionResult_t r;
  r = S; //copy selection config
  //Scan_t sps   //std::cout << "after sps" << std::endl;
  r = select(proj, P, S.cut + extra_cut);
  //std::cout << "after r = sps" << std::endl;
  return r;
};

template<typename Proj>
inline std::vector<SelectionResult_t> Analysis::select(Proj proj, const Scan_t & P, const std::vector<Selection_t> & S, std::string extra_cut) const {
  std::vector<SelectionResult_t> R;
  int i=0;
  for(const Selection_t & s : S) {
    SelectionResult_t r = select(proj, P, s, extra_cut);
    print(r,i);
    R.push_back(r);
    ++i;
  }
  print(fold(R,"all"),-1);
  return R;
};


template<typename Proj>
void  Analysis::set_effcor( Proj proj, Scan_t  & P,double Wref)  {
  //determine reference point
  DataSample_t  & ds0 = std::invoke(proj, *find_minimum(P, [&Wref](const auto & p) { return  std::abs(p.energy - Wref); } ) );
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
inline void Analysis::copy_efficiency(Proj proj, const Scan_t & mc, Scan_t  & psr) const
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
inline void  Analysis::set_effcor( Proj proj, std::vector<SelectionResult_t>  & P,double Wref)  {
  for(size_t i=0;i!=P.size();++i) {
    set_effcor( proj, P[i],Wref);
  }
};

template< typename Proj> 
void Analysis::copy_efficiency(Proj proj, const std::vector<SelectionResult_t> & M /*Monte Carlo*/, std::vector<SelectionResult_t> & D /*data*/) const
{
  assert(D.size() == M.size());
  for(size_t i=0;i!=D.size();++i) {
    copy_efficiency(proj, M[i], D[i]);
  }
}



inline void Analysis::save(const Scan_t & sr, std::string  filename, std::string default_lum) const
{
  std::cout << "Saving  to file: " << filename << std::endl;
  std::stringstream os;
  os << ibn::mformat("15.6", "#", "L,nb^-1","dL,nb^-1","W,MeV", "dW,MeV","Sw,MeV","dSw,MeV","Nsig","Nbb","Ngg","cor" "error", "eps", "error");
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
        TAU_SPREAD.value,    TAU_SPREAD.error,
        p.tt.N,
        p.bb.N,        
        p.gg.N,
        p.tt.effcor.value, p.tt.effcor.error, 
		p.tt.efficiency.value, p.tt.efficiency.error);
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

#include <TF1.h>
inline void draw_mhadr_efficiency( Analysis & A, std::string cut="") {
    auto jpsi = A.res(A.JPSI, cut);
    double eps_jpsi{0};
    for(const auto & p : jpsi) {
        eps_jpsi+=p.tt.efficiency.value;
    };
    eps_jpsi/=jpsi.size();
    auto gj = draw_efficiency(jpsi);


    TF1 * f_eps_jpsi = new TF1("f_eps_jpsi","[0]*(1 + (x-3.096)*0.01*1000*[1])");
    f_eps_jpsi->SetParameter(0,eps_jpsi);
    f_eps_jpsi->SetParameter(1,0);

    gStyle->SetOptFit();
    gj->Fit("f_eps_jpsi");

    auto psi  =  A.res(A.PSI, cut);
    double eps_psi{0};
    for(const auto & p : psi) {
        eps_psi+=p.tt.efficiency.value;
    };
    eps_psi/=psi.size();

    auto gp = draw_efficiency(psi);

    TF1 * f_eps_psi = new TF1("f_eps_psi","[0]*(1 + (x-3.686)*0.01*1000*[1])");
    f_eps_psi->SetParameter(0,eps_psi);
    f_eps_psi->SetParameter(1,0);
    gp->Fit("f_eps_psi");
};

inline void draw_mhadr_efficiency_vs_pt( Analysis & A, std::string cut="") {
    auto average_epsilon = [](auto & jpsi) {
        return std::accumulate(jpsi.begin(), jpsi.end(), 0.0, [](const auto sum, const auto & p2) { return sum + p2.tt.efficiency.value; } ) / jpsi.size();
    };

    TGraphErrors * g = new TGraphErrors;
    int i{0};
    for(double ptmin=0.08; ptmin<0.5; ptmin+=0.01) {
        std::string ptcut = cut && ("ptmin>"+std::to_string(ptmin));
        auto res = A.res(A.JPSI, ptcut);

        auto eps_min = std::min_element(res.begin(), res.end(), [](const auto & p1, const auto & p2) { return p1.tt.efficiency.value < p2.tt.efficiency.value;})->tt.efficiency;
        auto eps_max = std::max_element(res.begin(), res.end(), [](const auto & p1, const auto & p2) { return p1.tt.efficiency.value < p2.tt.efficiency.value;})->tt.efficiency;
        double eps_average = average_epsilon(res);
        std::cout << ptcut << "   " << eps_min.value << "  " << eps_max.value << "   " << eps_average << std::endl;

        ibn::valer<double> v = (eps_max-eps_min)/eps_average;
        g->SetPoint(i, ptmin, v.value);
        g->SetPointError(i, 0, v.error);
        ++i;
    }
    g->Draw("a*l");
};


template<typename Proj>
inline std::tuple<TGraphErrors*, TGraphErrors*> mhadr_eps_variate_var( Analysis & A, Proj proj, std::string var, double xmin, double xmax, double dx, std::string cut="" ) {
    auto average_epsilon = [](auto & jpsi) -> double{
        return std::accumulate(jpsi.begin(), jpsi.end(), 0.0, [](const auto sum, const auto & p2) { return sum + p2.tt.efficiency.value; } ) / jpsi.size();
    };

    auto average_energy = [](auto & jpsi) -> double {
        return std::accumulate(jpsi.begin(), jpsi.end(), 0.0, [](const auto sum, const auto & p2) { return sum + p2.tt.energy.value; } ) / jpsi.size();
    };

    TGraphErrors * g = new TGraphErrors;
    TGraphErrors * gminmax = new TGraphErrors;
    int i{0};
    auto & RES = std::invoke(proj, A);
    TF1 f_eps("f_eps","[0]*(1 + (x-[2])*0.01*1000*[1])") ;
    for(double x=xmin; x<xmax; x+=dx) {
        std::string lcut = cut && ibn::format("%s%f", var.c_str(), x);
        auto res = A.res(RES, lcut);
        //auto graph = draw_efficiency(res);
        //auto graph = make_efficiency_graph(res);
        auto graph = make_efficiency_graph(res);
        auto eps_min = std::min_element(res.begin(), res.end(), [](const auto & p1, const auto & p2) { return p1.tt.efficiency.value < p2.tt.efficiency.value;})->tt.efficiency;
        auto eps_max = std::max_element(res.begin(), res.end(), [](const auto & p1, const auto & p2) { return p1.tt.efficiency.value < p2.tt.efficiency.value;})->tt.efficiency;
        double eps_average = average_epsilon(res);
        std::cout << x << "   " << eps_min.value << "  " << eps_max.value << "   " << eps_average << std::endl;
        f_eps.SetParameter(0,eps_average);
        f_eps.SetParameter(1,0);
        f_eps.SetParameter(2, average_energy(res));
        f_eps.FixParameter(2, average_energy(res));
        graph->Fit("f_eps","0");
        g->SetPoint(i, x, f_eps.GetParameter(1));
        g->SetPointError(i, 0, f_eps.GetParError(1));

        ibn::valer<double> v = (eps_max-eps_min)/eps_average;
        gminmax->SetPoint(i, x, v.value*100);
        gminmax->SetPointError(i, 0, v.error*100);
        ++i;
    }
    return {g, gminmax};
}


template<typename Proj>
inline void draw_mhadr_efficiency_slope_vs_pt( Analysis & A, Proj proj, std::string cut="", double step=0.01) {

    auto [g_slope, g_minmax] = mhadr_eps_variate_var(A, proj, "ptmin>",0.08,0.5,step,cut);

    TCanvas * c = new TCanvas("draw_mhadr_efficiency_slope_vs_pt", "Efficiency variation vs pt", 1600, 1600);
    c->Divide(1,2);
    c->cd(1);
    g_slope->Draw("a*l");
    g_slope->GetXaxis()->SetTitle("pt, GeV");
    g_slope->GetYaxis()->SetTitle("slope, %/MeV");

    c->cd(2);
    g_minmax->Draw("*la");
    g_minmax->GetXaxis()->SetTitle("pt, GeV");
    g_minmax->GetYaxis()->SetTitle("(#varepsilon_{max}-#varepsilon_{min})/#varepsilon_{average}*100%");

};

template<typename Proj>
inline void draw_mhadr_efficiency_slope_vs_cos( Analysis & A, Proj proj, std::string cut="", double step=0.1) {

    auto [g_slope, g_minmax] = mhadr_eps_variate_var(A, proj, "maxctheta<",0.5,0.93,step,cut);

    TCanvas * c = new TCanvas("draw_mhadr_efficiency_slope_vs_maxcostheta", "Efficiency variation vs maxctheta", 1600, 1600);
    c->Divide(1,2);
    c->cd(1);
    g_slope->Draw("a*l");
    g_slope->GetXaxis()->SetTitle("|cos(#theta)|, GeV");
    g_slope->GetYaxis()->SetTitle("slope, %/MeV");

    c->cd(2);
    g_minmax->Draw("*la");
    g_minmax->GetXaxis()->SetTitle("|cos(#theta)|, GeV");
    g_minmax->GetYaxis()->SetTitle("(#varepsilon_{max}-#varepsilon_{min})/#varepsilon_{average}*100%");

    auto [g_slope2, g_minmax2] = mhadr_eps_variate_var(A, proj, "minctheta>",-0.93,-0.5, step,cut);

    TCanvas * c2 = new TCanvas("draw_mhadr_efficiency_slope_vs_mincostheta", "Efficiency variation vs minctheta", 1600, 1600);
    c2->Divide(1,2);
    c2->cd(1);
    g_slope2->Draw("a*l");
    g_slope2->GetXaxis()->SetTitle("|cos(#theta)|, GeV");
    g_slope2->GetYaxis()->SetTitle("slope, %/MeV");

    c2->cd(2);
    g_minmax2->Draw("*la");
    g_minmax2->GetXaxis()->SetTitle("|cos(#theta)|, GeV");
    g_minmax2->GetYaxis()->SetTitle("(#varepsilon_{max}-#varepsilon_{min})/#varepsilon_{average}*100%");

};

template<typename Proj>
inline void draw_mhadr_efficiency_slope_vs_S( Analysis & A, Proj proj, std::string cut="") {

    auto [g_slope, g_minmax] = mhadr_eps_variate_var(A, proj, "S>",0.0,0.3,0.005,cut);

    TCanvas * c = new TCanvas("draw_mhadr_efficiency_slope_vs_S", "Efficiency variation vs S", 1600, 1600);
    c->Divide(1,2);
    c->cd(1);
    g_slope->Draw("a*l");
    g_slope->GetXaxis()->SetTitle("S, GeV");
    g_slope->GetYaxis()->SetTitle("slope, %/MeV");

    c->cd(2);
    g_minmax->Draw("*la");
    g_minmax->GetXaxis()->SetTitle("S, GeV");
    g_minmax->GetYaxis()->SetTitle("(#varepsilon_{max}-#varepsilon_{min})/#varepsilon_{average}*100%");

};

inline TCanvas * draw_res(const Scan_t & scan, std::string var, std::string cut="") {
    auto c = make_canvas("var");
    //c->Divide(4,2);
    int i=1;
    int color=1;
    TLegend * l = new TLegend(0.8,0.8,1.0,1.0);
    for(auto & sp :scan) {
        std::string gopt="PLC PMC NORM HIS";
        if(i!=1) gopt = "PLC PMC NORM same HIS";
        if(color==10) ++color;
        //sp.tt.tree->SetMarkerColor(color);
        //sp.tt.tree->SetLineColor(color);
        sp.tt.tree->SetLineWidth(3);
        //c->cd(i+1);
        sp.tt.tree->Draw(var.c_str(), cut.c_str(), gopt.c_str());
        auto h = sp.tt.tree->GetHistogram();
        //h->SetMarkerColor(i+1);
        //h->SetLineColor(i+1);
        l->AddEntry(h,("point " + std::to_string(i)).c_str());
        ++i;
        ++color;
        c->Modified();
        c->Update();
    };
    l->Draw();
    c->Modified();
    c->Update();
    return c;
};

#include <TH2F.h>
//рисует разницу в гистограммах для разных точек сканирования
inline TCanvas * draw_res_diff(const Scan_t & scan, std::string var, int Nbin, double xmin, double xmax, std::string cut="") {
    auto c = make_canvas("var");
    //c->Divide(4,2);
    int i=1;
    TLegend * l = new TLegend(0.8,0.8,1.0,1.0);
    std::vector<TH1F*> hs;
    std::vector<int> color;
    for(auto & sp :scan) {
        std::string gopt="PLC PMC NORM HIS";
        if(i!=1) gopt = "PLC PMC NORM same HIS";
        sp.tt.tree->SetLineWidth(3);
        //c->cd(i+1);
        std::string v = ibn::format("%s>>his%d(%d,%f,%f)", var.c_str(), HISTO_INDEX++, Nbin, xmin, xmax);
        sp.tt.tree->Draw(v.c_str(), cut.c_str(), gopt.c_str());
        auto h = (TH1F*)sp.tt.tree->GetHistogram();
        color.push_back(h->GetLineColor());
        //h->SetMarkerColor(i+1);
        //h->SetLineColor(i+1);
        l->AddEntry(h,("point " + std::to_string(i)).c_str());
        ++i;
        c->Modified();
        c->Update();
        hs.push_back(h);
    };
    l->Draw();
    c->Modified();
    c->Update();

    auto c2 = make_canvas("var");
    TH2F * his = new TH2F("haxis","haxis", Nbin, xmin, xmax, Nbin, -0.05, 0.05);
    his->Draw();
    int ref_point = hs.size()/2;
    auto his_ref = (TH1F*) hs[ref_point]->Clone("ref_his");
    TLegend * l2 = new TLegend(0.8,0.8,1.0,1.0);
    for(int i=0;i<hs.size();++i) {
        //if(i==ref_point) continue;
        auto hh = (TH1F*) hs[i]->Clone(ibn::format("hh%d",HISTO_INDEX).c_str());
        hh->Add(his_ref,-1);
        //hh->Divide(his_ref);
        hh->Draw("same");
        l2->AddEntry(hh,("point " + std::to_string(i+1)).c_str());
    }
    l2->Draw();

    return c;
};

inline void draw_res(Analysis & A, std::string var, std::string cut) {
    TCanvas *  c{nullptr};
    gStyle->SetPalette(kRainBow);
    c = draw_res(A.JPSI.signal.data,var, cut);
    c->SetTitle((var + " for data jpsi").c_str());
    c->Write();
    c->SaveAs(("mhadr/mh_"+var+"jpsi_data.pdf").c_str());
    c=draw_res(A.JPSI.signal.mc,  var, cut);
    c->SetTitle((var + " for mc jpsi").c_str());
    c->Write();
    c->SaveAs(("mhadr/mh_"+var+"jpsi_mc.pdf").c_str());
    c= draw_res(A.PSI.signal.data, var, cut);
    c->SetTitle((var + " for data psi(2S)").c_str());
    c->Write();
    c->SaveAs(("mhadr/mh_"+var+"psi_data.pdf").c_str());
    c= draw_res(A.PSI.signal.mc,   var, cut);
    c->SetTitle((var + " for mc   psi(2S)").c_str());
    c->Write();
    c->SaveAs(("mhadr/mh_"+var+"psi_mc.pdf").c_str());
};

#include <TFile.h>
inline void draw_res( Analysis & A) {
    TFile f("tmp.root", "RECREATE");

    //draw_res(A, "ptmin","S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2");
    //draw_res(A, "maxctheta","ptmin>0.2 && S>0.06 && minctheta>-0.8 && Nchc>2");
    //draw_res(A, "minctheta","ptmin>0.2 && S>0.06 && maxctheta<0.8 && Nchc>2");
    //draw_res(A, "S","ptmin>0.2 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2");
    //draw_res(A, "Nchc","ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8");
    draw_res(A, "Nchgcemc","");
    //f.Write();
    f.Close();
}
