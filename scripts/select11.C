/*
 * =====================================================================================
 *
 *       Filename:  select.C
 *
 *    Description:  Selection configuration of tau events
 *
 *        Version:  11
 *        Created:  10.02.2019 14:27:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#include "selection.C"
#include "sys.C"
#include <map>
#include "Config.h"

std::string LOCAL_TAUFIT = "taufit --minos --pdgshift --lum=default --tau-spread=1.258 --energy-correction=+0.011 --free-energy --free-luminosity --free-effcor --draw-tau-mass-precision=4";

//const char * runtable_name = "../TauTau/share/scan_points_ems3_privalov_lum.txt";
//const char * runtable_name = "scan_points_ems3_bhabha_lum.txt";
const char * runtable_name = "../TauTau/share/scan_points_ems3_online_lum.txt";
//const char * runtable_name = "../TauTau/share/all_scan_points_ems3.txt";

auto RUNTABLE  = read_my_runtable(runtable_name);


constexpr long N0MC  = 1e6;

//read data from directory "data". Runs are combined into points according to runtable 
auto DATA        = read_data("data", RUNTABLE);

//Monte Carlo simulation of the signal
auto SIGNAL          = read_mc("mc/signal", RUNTABLE, N0MC);
//auto SIGNAL          = read_mc("mc/signal_sys", RUNTABLE, N0MC);

//background GALUGA
std::map<std::string, Scan_t> GALUGA =
{
  {"ee"   , read_mc("mc/galuga/ee"   , RUNTABLE, N0MC)} ,
  {"uu"   , read_mc("mc/galuga/uu"   , RUNTABLE, N0MC)} ,
  {"pipi" , read_mc("mc/galuga/pipi" , RUNTABLE, N0MC)} ,
  {"KK"   , read_mc("mc/galuga/KK"   , RUNTABLE, N0MC)}
};

//hadronic background
auto HADR       = read_mc("mc/hadrons", RUNTABLE, N0MC);
//bhabha background
auto BB         = read_mc("mc/bb", RUNTABLE, N0MC);
auto UU         = read_mc("mc/uu", RUNTABLE, N0MC);
auto PIPI         = read_mc("mc/pipi", RUNTABLE, N0MC);

// for luminocity measurement
auto GG         = read_mc("mc/gg", RUNTABLE, N0MC);

std::vector<ScanRef_t> BGs ={HADR, BB, UU, GG, GALUGA["ee"], GALUGA["uu"],GALUGA["pipi"], GALUGA["KK"]};

Simulation_t MC = { SIGNAL, BGs }; 

//std::string LOCAL_BB_SEL = "(acol-TMath::Pi())>-0.03";
std::string LOCAL_BB_SEL = "(acol-TMath::Pi())>-0.04 && abs(cos(theta[0])) < 0.8 && abs(cos(theta[1])) < 0.8 && Ep[0]>0.8 && Ep[1]>0.8";
std::string LOCAL_GG_SEL = "";

//particla identification configuration
std::vector<ParticleID_t> PID = 
{
  {
    {"e", {"Ep[#]>0.8"}, 
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
    {"PI", { "!(e#||u#)" }, //this pi meson is like previose pi but without moment cut. It will be used in rho id
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
//      ,restrict("p", 0.8, 1.05) 
      }},
      {"rho", {" PI# && 0.12 < Mpi0[0] && Mpi0[0] < 0.14"},
      { 
        restrict("Mrho", 0.6,1.0),
      }}, 

      {"X", { "chi2_dedx_e[#]>3" 
      //{"X", { "!e#"
          },
      { 
        restrict("Ep",           0, 0.8), 
      }},
  }
};


Selection_t SEL11 =
{
  "sel11", //selection name
  //common cuts
  "Nc==2"
  "&& q[0]==-1 && q[1]==1"
  "&& p[0] < 1.1 && p[1] < 1.1"
  "&& pt[0] > 0.2 && pt[1] > 0.2"
  "&& barrel"
  "&& ( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )"
  , 
  PID,
  { 
    {"eX",  "Nn<=2 &&  eX   && ptem>0.25"},
  },
};


auto & SEL = SEL11;
std::vector<SelRef_t> SELS = {SEL11};



void doall(Scan_t & D = DATA/* data */ , Selection_t & S=SEL, Scan_t & M = SIGNAL, Config cfg = CFG)
{
  std::cout << "Configurateion: " << cfg.name << "  luminosity: " << cfg.lum << std::endl;
  S.apply_common_cut();
  set_kptem(D,1.0); 
  //GG_SEL = LOCAL_GG_SEL + (cfg.gglum_extra_cut != "" ? "&&"+cfg.gglum_extra_cut : "");
  //BB_SEL = LOCAL_BB_SEL + (cfg.eelum_extra_cut != "" ? "&&"+cfg.eelum_extra_cut : "");
  GG_SEL = LOCAL_GG_SEL && cfg.gglum_extra_cut;
  BB_SEL = LOCAL_BB_SEL && cfg.eelum_extra_cut;
  TAUFIT = LOCAL_TAUFIT + " " + cfg.taufit_extra_opt;
  measure_luminosity(D, BB, GG, N0MC);
  auto result = select(D,S); 
  print_luminosity(D);
  set_kptem(M,1.0);
  set_efficiency(result,M,N0MC);
  std::string file_to_fit  = S.title;
  file_to_fit += "_" + cfg.energy_version;
  file_to_fit += "_"+cfg.lum;
  if(cfg.lum == "ee" && cfg.eelum_extra_cut != "") {
    std::string c = cfg.eelum_extra_cut;
    c.erase(std::remove(c.begin(), c.end(), ' '), c.end());
    c.erase(std::remove(c.begin(), c.end(), '&'), c.end());
    file_to_fit += "_"+c;
  }
  if(cfg.lum == "gg" && cfg.gglum_extra_cut != "") {
    std::string c = cfg.gglum_extra_cut;
    c.erase(std::remove(c.begin(), c.end(), ' '), c.end());
    c.erase(std::remove(c.begin(), c.end(), '&'), c.end());
    file_to_fit += "_"+c;
  }
  if(cfg.name != "") file_to_fit += "_"+cfg.name;
  std::string prefix=file_to_fit;
  file_to_fit += ".txt";
  std::cout << "Ouput file: " << file_to_fit << "\n";
  fit(result, file_to_fit,  S.title, cfg.lum);
  //sleep(10);
  make_tex(print_tex(result, file_to_fit, prefix + "_fit.pdf"),prefix + ".tex");
};




void select() 
{
  BB_SEL = LOCAL_BB_SEL;
  GG_SEL = LOCAL_GG_SEL;
  TAUFIT = LOCAL_TAUFIT;
  double Kptem = 1.0;
  for(auto & S : SELS) {
    S.get().apply_common_cut();
  }

  std::vector<ScanRef_t> LUM_MCs = {BB,GG};
  std::vector<ScanRef_t> BG_MCs =  {HADR, UU, PIPI};
  std::vector<ScanRef_t> BGall_MCs =  BG_MCs;
  for( auto & p: GALUGA) BGall_MCs.push_back(p.second);

  read_tau_cross_section("../TauTau/share/tau_cross_section.txt", SIGNAL);
  read_bhabha_cross_section("../TauTau/share/bhabha_cross_section.txt", BB);
  read_gg_cross_section("../TauTau/share/gg_cross_section.txt", GG);
  read_galuga_cross_section("../TauTau/share/galuga_cross_section.txt", GALUGA);
  read_mumu_or_pipi_cross_section("../TauTau/share/mumu_cross_section.txt", UU);
  read_pipi_cross_section("../TauTau/share/pipi_cross_section.txt", PIPI);
  read_hadron_cross_section("../TauTau/share/hadron_cross_section.txt", HADR);

  set_pid_kptem(DATA       , PID , Kptem);
  set_pid_kptem(SIGNAL         , PID , Kptem);
  for(auto d : BGall_MCs) set_pid_kptem(d,PID,Kptem);
  for(auto d : LUM_MCs)   set_pid_kptem(d,PID,Kptem);

  measure_luminosity(DATA,BB,GG,1e6);
  set_luminosity(DATA,SIGNAL);
  for(auto d : BGall_MCs) set_luminosity(DATA,d);

  /*
  set_pid_kptem(HADR       , PID , Kptem);
  set_pid_kptem(BB         , PID , Kptem);
  set_pid_kptem(GG         , PID , Kptem);
  set_pid_kptem(UU         , PID , Kptem);
  set_pid_kptem(PIPI         , PID , Kptem);
  */
 //doall(SEL,1,DATA,SIGNAL,"sel","gg"); 
 //sys("cos_theta_mis3.sys","cos_theta_mis2, GeV");
}



