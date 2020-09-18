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

//ems3
std::string LOCAL_TAUFIT = "taufit --minos --pdgshift --lum=default --tau-spread=1.244 --ems-cmenergy-shift=+0.072 --free-energy --free-luminosity --free-effcor --draw-tau-mass-precision=4 --base_shift=1.23091e-01 --draw-diff";
//ems2
//std::string LOCAL_TAUFIT = "taufit --minos --pdgshift --lum=default --tau-spread=1.239      --energy-correction=+0.0391396 --free-energy --free-luminosity --free-effcor --draw-tau-mass-precision=4 --base_shift=1.23091e-01 --draw-diff";


//std::string runtable_name = "../TauTau/share/all_scan_points_ems3.txt";
//std::string runtable_name = "ems31.txt";
//std::string runtable_name = "ems32.txt";
//std::string runtable_name = "ems3_with_exp.txt";
//std::string runtable_name = "../TauTau/share/all_scan_points_ems2.txt";
std::string runtable_name = "ems20.txt";
auto ALLRUNTABLE      = read_my_runtable(runtable_name);

auto RUNTABLE      = read_my_runtable(runtable_name, DATA_TAU);
//auto RUNTABLE      = read_my_runtable("../TauTau/share/all_scan_points_ems3.txt");
auto JPSI_RUNTABLE =  read_my_runtable(runtable_name, DATA_JPSI);
auto PSIP_RUNTABLE =  read_my_runtable(runtable_name, DATA_PSIP);

//#ifdef EMS2
////auto RUNTABLE      = read_my_runtable("../TauTau/share/tau_scan_points_ems2.txt");
//auto JPSI_RUNTABLE =  read_my_runtable(rutab, DATA_JPSI);
//auto PSIP_RUNTABLE =  read_my_runtable("../TauTau/share/psip_scan_points_ems2.txt", DATA_PSIP);
//#endif


constexpr long N0MC  = 1e6;

//read data from directory "data". Runs are combined into points according to runtable 
auto DATA      = read_data("data", RUNTABLE);
auto JPSI      = read_mh("mhdata",JPSI_RUNTABLE);
auto PSIP      = read_mh("mhdata",PSIP_RUNTABLE);

auto JPSIMC    = read_mc_mh("mcmh", JPSI_RUNTABLE);
auto PSIPMC    = read_mc_mh("mcmh", PSIP_RUNTABLE);

auto JPSILUM   = read_privalov_lum("mhlum",JPSI_RUNTABLE);
auto PSIPLUM   = read_privalov_lum("mhlum",PSIP_RUNTABLE);
auto JPSILUMMC   = read_mc_mhlum("mhlum",JPSI_RUNTABLE,1e5);
auto PSIPLUMMC  = read_mc_mhlum("mhlum",PSIP_RUNTABLE,1e5);

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
std::string LOCAL_BB_SEL = "(acol-TMath::Pi())>-0.04 && abs(cos(theta[0])) < 0.8 && abs(cos(theta[1])) < 0.8 && Ep[0]>0.8 && Ep[1]>0.8 && abs(z[0])<10 && abs(z[1])<10 && vxy[0]<1.0 && vxy[1]<1.0";
std::string LOCAL_GG_SEL = "";
std::string LOCAL_MH_SEL = "Echmin>0.05 && Nchc==Nchgcemc && ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2";
std::string PRIVALOV_GG_SEL="fabs(cos(theta[0]))<0.8 && fabs(cos(theta[1]))<0.8 && theta[0]+theta[1]-TMath::Pi()<0.055 && theta[0]+theta[1]-TMath::Pi()>-0.06 && fabs(phi[0]-phi[1])-TMath::Pi()<0.014 && fabs(phi[0]-phi[1])-TMath::Pi()>-0.054";

//particla identification configuration
std::vector<ParticleID_t> PID = 
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


Selection_t SEL12 =
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

Selection_t SELEE =
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
    //{"ee",      "NnE50==0 && ee && missed_photon_angle"},
    {"uu",      "NnE50==0 && uu && missed_photon_angle"},
  },
};

Selection_t SELN =
{
  "seln", //selection name
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
    {"ee",      "NnE25==0 && ee && Emis>1.8 && ptem25 > 0.34*(1+cos_theta_mis2^2) && acop>0.1"},
    {"μμ",      "NnE25==0 && uu && ptem25>0.2"},
    {"ππ",      "NnE25==0 && pipi && ptem25>0.15 &&  1.7 < Emis && Emis < 1.95"},
    {"μπ",      "NnE25==0 && upi && barrel"},
    {"πρ",      "Nng==2 && Npi0 == 1 && pirho && acop > 1.2 && Emis>1.5 && Emis <1.8"},
    {"μρ",      "Nng==2 && Npi0 == 1 && urho && acop>0.1"},
    //{"Xρ",      "Nng==2   && XX && Mpi0[0] < 0.14 && Mpi0[0]>0.12" },
    //{"XX",      "NnE50==0 && XX" },
  },
};

Selection_t SEL_HIGH_BG =
{
  "sel_high_bg", //selection name
  "Ncg<=3"
  "&& Ncg==Ncc"
  "&& E[0]>0.025 && E[1]>0.025"
  "&& q[0]==-1 && q[1]==1"
  "&& cgood"
  "&& abs(vz[0])<12 && abs(vz[1])<12"  // cosm
  "&& vxy[0]<1.2 && vxy[1]<1.2"        //cosm
  "&& 1 < tof[0] && tof[0] < 6 && 1 < tof[1] && tof[1] < 6" //cosm
  "&& p[0] < 1.2 && p[1] < 1.2" //ISR
  "&& pt[0] > 0.05 && pt[1] > 0.05" //beam
  "&& ptem50>0.20 && ptem50<1.1" //twogam
  , 

  PID,
  { 
    {"eX",      "NnE50==0 && eX && abs(cos_theta_mis)<0.95"},
    {"eρ",     "(Nng==2 || Nng==3) && eX && good_emc_time" },
  },
};


auto & SEL = SEL12;
std::vector<SelRef_t> SELS = {SEL12, SEL_HIGH_BG, SELN,SELEE};



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

#include <sys/stat.h>
#include <sys/types.h>
#include <sys/unistd.h>
void dores(std::string suffix) {
  //system("rm shift.txt");
  std::cout << "Proceeding J/psi" << std::endl;
  std::cout << "Selecting multihadronic events" << std::endl;
  auto jpsi = select(JPSI, LOCAL_MH_SEL);
  std::cout << "measuring efficiency" << std::endl;
  auto mcjpsi = measure_efficiency(JPSIMC, LOCAL_MH_SEL);
  print(jpsi);
  set_efficiency(jpsi,mcjpsi,1e6);
  std::string jpsi_filename = "jpsi"+suffix+".txt";
  save_and_fitmh(jpsi,jpsi_filename,"default", "");

  std::cout << "Proceeding psi(2S). Selecting multihadronics ... ";
  auto psip = select(PSIP, LOCAL_MH_SEL);
  std::cout << "measuring efficiency ... ";
  auto mcpsip = measure_efficiency(PSIPMC, LOCAL_MH_SEL);
  print(psip);
  set_efficiency(psip,mcpsip,1e6);
  std::cout << std::endl;
  std::string psip_filename = "psip"+suffix+".txt";
  save_and_fitmh(psip,psip_filename,"default", "");
  //fitmh(psip,psip_filename);
  //struct stat st;
  //stat("shift.txt", &st);
  //time_t t0;
  //t0 = st.st_mtim.tv_sec;
  //do {
  //  sleep(1);
  //  stat("shift.txt", &st);
  //} while( st.st_mtim.tv_sec == t0 );
  //t0 = st.st_mtim.tv_sec;
  //do {
  //  sleep(1);
  //  stat("shift.txt", &st);
  //} while( st.st_mtim.tv_sec == t0 );

  std::string res_filename = "res"+suffix+".txt";
  system(("cat " + jpsi_filename + "  " + psip_filename + " > " + res_filename).c_str());
  system((PSIFIT + " --both 1 " + res_filename).c_str());
}


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


  std::cout << "Setting  PID" << std::endl;
  set_pid_kptem(DATA       , PID , Kptem);
  set_pid_kptem(SIGNAL         , PID , Kptem);
  for(auto d : BGall_MCs) set_pid_kptem(d,PID,Kptem);
  for(auto d : LUM_MCs)   set_pid_kptem(d,PID,Kptem);

  std::cout << "Measuring luminosity" << std::endl;
  measure_luminosity(DATA,BB,GG,1e6);
  set_luminosity(DATA,SIGNAL);
  for(auto d : BGall_MCs) set_luminosity(DATA,d);
  set_gg_luminosity(DATA,BB);

  std::cout << "Add cross section for gg luminosity (Privalov) " << std::endl;
  read_privalov_gg_cross_section("../TauTau/share/privalov_gg_cross_section.txt", JPSILUMMC);
  read_privalov_gg_cross_section("../TauTau/share/privalov_gg_cross_section.txt", PSIPLUMMC);
  std::cout << "Measuring privalov luminosity" << std::endl;
  measure_gg_luminosity(JPSILUM, JPSILUMMC, 1e5, PRIVALOV_GG_SEL);
  measure_gg_luminosity(PSIPLUM, PSIPLUMMC, 1e5, PRIVALOV_GG_SEL);
  set_res_gg_luminosity(JPSI, JPSILUM);
  set_res_gg_luminosity(PSIP, PSIPLUM);

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



