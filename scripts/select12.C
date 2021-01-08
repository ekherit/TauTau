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

//ems3 mh noemc mylum2
std::string LOCAL_TAUFIT = "taufit --minos --pdgshift --lum=default --tau-spread=1.2299 --ems-cmenergy-shift=-0.0201182 --free-energy --free-luminosity --free-effcor --draw-tau-mass-precision=4  --draw-diff";
//ems2
//std::string LOCAL_TAUFIT = "taufit --minos --pdgshift --lum=default --tau-spread=1.2266 --ems-cmenergy-shift=0.0613 --free-energy --free-luminosity --free-effcor --draw-tau-mass-precision=4  --draw-diff --base-shift=1.19568e-01 ";
//ems32 base
//std::string LOCAL_TAUFIT = "taufit --minos --pdgshift --lum=default --tau-spread=1.238 --ems-cmenergy-shift=-0.0128 --free-energy --free-luminosity --free-effcor --draw-tau-mass-precision=4  --draw-diff";

//ems2
//std::string LOCAL_TAUFIT = "taufit --minos --pdgshift --lum=default --tau-spread=1.235 --ems-cmenergy-shift=0.0687 --free-energy --free-luminosity --free-effcor --draw-tau-mass-precision=4  --draw-diff --base-shift=1.19056e-01 ";

//std::string LOCAL_TAUFIT = "taufit --minos --pdgshift --lum=default --tau-spread=1.248 --ems-cmenergy-shift=-0.0107 --free-energy --free-luminosity --free-effcor --draw-tau-mass-precision=4 --base-shift=1.21898e-01 --draw-diff";
//ems3
//std::string LOCAL_TAUFIT = "taufit --minos --pdgshift --lum=default --tau-spread=1.248 --ems-cmenergy-shift=-0.0107 --free-energy --free-luminosity --free-effcor --draw-tau-mass-precision=4 --base-shift=1.23091e-01 --draw-diff";
//ems2
//std::string LOCAL_TAUFIT = "taufit --minos --pdgshift --lum=default --tau-spread=1.239      --energy-correction=+0.0391396 --free-energy --free-luminosity --free-effcor --draw-tau-mass-precision=4 --base-shift=1.23091e-01 --draw-diff";


//std::string runtable_name = "../TauTau/share/all_scan_points_ems2.txt";
std::string runtable_name = "../TauTau/share/all_scan_points_ems32.txt";
auto ALLRUNTABLE      = read_my_runtable(runtable_name);

auto RUNTABLE      = read_my_runtable(runtable_name, DATA_TAU);
auto JPSI_RUNTABLE =  read_my_runtable(runtable_name, DATA_JPSI);
auto PSIP_RUNTABLE =  read_my_runtable(runtable_name, DATA_PSIP);
//auto RES_RUNTABLE  =  read_my_runtable(runtable_name, DATA_RES);



constexpr long N0MC  = 1e6;

//read data from directory "data". Runs are combined into points according to runtable 
auto DATA      = read_data("data", RUNTABLE);


//auto RES      = read_mh("mhdata", RES_RUNTABLE);
//auto RESMC    = read_mc_mh("mcmh", RES_RUNTABLE, N0MC);
//auto RESGG    = read_data("data", RES_RUNTABLE);
//auto RESGGMC  = read_mc("mc/gg", RES_RUNTABLE,1e5);


auto JPSI      = read_mh("mhdata",JPSI_RUNTABLE);
auto PSIP      = read_mh("mhdata",PSIP_RUNTABLE);
auto JPSIMC    = read_mc_mh("mcmh", JPSI_RUNTABLE);
auto PSIPMC    = read_mc_mh("mcmh", PSIP_RUNTABLE);

auto JPSIGG  = read_data("data", JPSI_RUNTABLE);
auto JPSIGGMC  = read_mc("mc/gg", JPSI_RUNTABLE,1e5);

auto PSIGG  = read_data("data", PSIP_RUNTABLE);
auto PSIGGMC  = read_mc("mc/gg", PSIP_RUNTABLE,1e5);

//
//auto JPSILUM   = read_privalov_lum("mhlum",JPSI_RUNTABLE);
//auto PSIPLUM   = read_privalov_lum("mhlum",PSIP_RUNTABLE);
//auto JPSILUMMC   = read_mc_mhlum("mhlum",JPSI_RUNTABLE,1e5);
//auto PSIPLUMMC  = read_mc_mhlum("mhlum",PSIP_RUNTABLE,1e5);

//Monte Carlo simulation of the signal
auto SIGNAL          = read_mc("mc/signal", RUNTABLE, N0MC);
//auto SIGNAL          = read_mc("mc/signal_sys", RUNTABLE, N0MC);

//two-gamma background
Scan_t EEEE       = read_mc("mc/galuga/ee"   , RUNTABLE, 1e6);
Scan_t EEUU       = read_mc("mc/galuga/uu"   , RUNTABLE, 1e6);
Scan_t EEPIPI     = read_mc("mc/galuga/pipi"   , RUNTABLE, 1e6);
Scan_t EEKK       = read_mc("mc/galuga/KK"   , RUNTABLE, 1e6);

std::map<std::string, ScanRef_t> GALUGA = {
  {"ee", EEEE},
  {"uu",  EEUU},
  {"pipi", EEPIPI},
  {"KK", EEKK},
};
/*
std::map<std::string, Scan_t> GALUGA =
{
  {"ee"   , read_mc("mc/galuga/ee"   , RUNTABLE, N0MC)} ,
  {"uu"   , read_mc("mc/galuga/uu"   , RUNTABLE, N0MC)} ,
  {"pipi" , read_mc("mc/galuga/pipi" , RUNTABLE, N0MC)} ,
  {"KK"   , read_mc("mc/galuga/KK"   , RUNTABLE, N0MC)}
}; 
*/

//hadronic background
auto HADR       = read_mc("mc/hadrons", RUNTABLE, N0MC);
//bhabha background
auto BB         = read_mc("mc/bb", RUNTABLE, N0MC);
auto UU         = read_mc("mc/uu", RUNTABLE, N0MC);
auto PIPI         = read_mc("mc/pipi", RUNTABLE, N0MC);

// for luminocity measurement
auto GG         = read_mc("mc/gg", RUNTABLE, N0MC);

std::vector<ScanRef_t> BGs ={HADR, BB, UU, GG, EEEE, EEUU, EEPIPI, EEKK};

Simulation_t MC = { SIGNAL, BGs }; 

//std::string LOCAL_BB_SEL = "(acol-TMath::Pi())>-0.03";
std::string LOCAL_BB_SEL = "(acol-TMath::Pi())>-0.04 && abs(cos(theta[0])) < 0.8 && abs(cos(theta[1])) < 0.8 && Ep[0]>0.8 && Ep[1]>0.8 && abs(z[0])<10 && abs(z[1])<10 && vxy[0]<1.0 && vxy[1]<1.0";
std::string LOCAL_GG_SEL = "";
//studi systematic in GG luminosity
//std::string LOCAL_GG_SEL = "fabs(cos(theta[0]))<0.7 && fabs(cos(theta[1]))<0.7 && abs(dtheta)<0.03 && E_Eb<1.025 && 0.85 < E_Eb && -0.04 < dphi && dphi< 0.0 && N0==2";
//std::string LOCAL_GG_SEL = "N0==2";

//std::string LOCAL_MH_SEL = "Echmin>0.05 && Nchc==Nchgcemc && ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2";
//std::string LOCAL_MH_SEL = "Echmin>0.1 && Nchc==Nchgcemc && ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2";
//std::string LOCAL_MH_SEL = "Echmin>0.05 && Nchc==Nchgcemc && ptmin>0.25 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2";
//std::string local_mh_sel = "echmin>0.05 && nchc==nchgcemc && ptmin>0.2 && s>0.1 && maxctheta<0.8 && minctheta>-0.8 && nchc>2";
//std::string LOCAL_MH_SEL = "Echmin>0.05 && Nchc==Nchgcemc && ptmin>0.2 && S>0.06 && maxctheta<0.7 && minctheta>-0.7 && Nchc>2";
//std::string LOCAL_MH_SEL = "Echmin>0.05 && Nchc==Nchgcemc && ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>3";
//теперь это новый отбор многоадронных
std::string LOCAL_MH_SEL = "ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2";
//std::string LOCAL_MH_SEL = "S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2";
//std::string LOCAL_MH_SEL = "ptmin>0.19 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2";
//std::string LOCAL_MH_SEL = "ptmin>0.2 && S>0.15 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2";
//std::string LOCAL_MH_SEL = "ptmin>0.2 && S>0.0 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2";
//std::string LOCAL_MH_SEL = "ptmin>0.2 && S>0.06 && maxctheta<0.75 && minctheta>-0.75 && Nchc>2";
//std::string LOCAL_MH_SEL = "ptmin>0.2 && S>0.06 && maxctheta<0.85 && minctheta>-0.85 && Nchc>2";
//std::string LOCAL_MH_SEL = "ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>3";
//std::string LOCAL_MH_SEL = "ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2 && zmax < 5 && zmin>-5 && rhomax<0.5 && pmax < 1.9";

//std::string PRIVALOV_GG_SEL="fabs(cos(theta[0]))<0.8 && fabs(cos(theta[1]))<0.8 && theta[0]+theta[1]-TMath::Pi()<0.055 && theta[0]+theta[1]-TMath::Pi()>-0.06 && fabs(phi[0]-phi[1])-TMath::Pi()<0.014 && fabs(phi[0]-phi[1])-TMath::Pi()>-0.054";
//std::string PRIVALOV_GG_SEL="fabs(cos(theta[0]))<0.7 && fabs(cos(theta[1]))<0.7 && theta[0]+theta[1]-TMath::Pi()<0.04 && theta[0]+theta[1]-TMath::Pi()>-0.04 && fabs(phi[0]-phi[1])-TMath::Pi()<0.01 && fabs(phi[0]-phi[1])-TMath::Pi()>-0.05";
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
//    {"μμ",      "NnE50==0 && uu "},
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

void reslum(Scan_t & data, Scan_t & lum, Scan_t mclum) {
  read_gg_cross_section("../TauTau/share/gg_cross_section.txt", mclum);
  cout << "Measuring res luminosity " << std::endl;
  measure_gg_luminosity(lum, mclum,1e5,LOCAL_GG_SEL);
  set_res_gg_luminosity(data, lum);
  print(data);
}

Scan_t  res(std::string fname, Scan_t & data, Scan_t &mcdata, Scan_t & lum, Scan_t mclum) {
  reslum(data,lum, mclum);
  std::cout << "Selecting multihadronic events" << std::endl;
  auto res = select(data, LOCAL_MH_SEL);
  auto mcres = measure_efficiency(mcdata, LOCAL_MH_SEL);
  set_efficiency(res,mcres,1e6);
  print(res);
  savemh(res,fname);
  return res;
};

void res(std::string suffix="_test") {
  std::string jpsiname = "jpsi"+suffix+".txt";
  std::string psiname = "psip"+suffix+".txt";
  std::string resname = "res"+suffix+".txt";
  res(jpsiname, JPSI, JPSIMC, JPSIGG, JPSIGGMC);
  res(psiname, PSIP, PSIPMC, PSIGG, PSIGGMC);
  system(("cat "+jpsiname + " " + psiname + " > " +  resname).c_str());
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
  //read_galuga_cross_section("../TauTau/share/galuga_cross_section.txt", GALUGA);
  read_galuga_cross_section("../TauTau/share/galuga_cross_section.txt", EEEE, EEUU,EEPIPI,EEKK);
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

  //std::cout << "Add cross section for gg luminosity (Privalov) " << std::endl;
  //read_privalov_gg_cross_section("../TauTau/share/privalov_gg_cross_section.txt", JPSILUMMC);
  //read_privalov_gg_cross_section("../TauTau/share/privalov_gg_cross_section.txt", PSIPLUMMC);
  //std::cout << "Measuring privalov luminosity" << std::endl;
  //measure_gg_luminosity(JPSILUM, JPSILUMMC, 1e6, PRIVALOV_GG_SEL);
  //measure_gg_luminosity(PSIPLUM, PSIPLUMMC, 1e6, PRIVALOV_GG_SEL);
  //set_res_gg_luminosity(JPSI, JPSILUM);
  //set_res_gg_luminosity(PSIP, PSIPLUM);

  /*
  set_pid_kptem(HADR       , PID , Kptem);
  set_pid_kptem(BB         , PID , Kptem);
  set_pid_kptem(GG         , PID , Kptem);
  set_pid_kptem(UU         , PID , Kptem);
  set_pid_kptem(PIPI         , PID , Kptem);
  */
 //doall(SEL,1,DATA,SIGNAL,"sel","gg"); 
 //sys("cos_theta_mis3.sys","cos_theta_mis2, GeV");
  //set_luminosity(JPSILUM
}

void calc_online_lum(const Scan_t & SP, std::string run_info_filename) {
  std::ifstream ifs(run_info_filename);
  if(!ifs) {
    std::cerr << "Unalbe to open file: " << run_info_filename << "\n" << std::endl;
  }
  std::regex comment_re(R"(^\s*#.*)"); //for extract energy
  std::smatch match;
  std::string line;
  std::map<int, double> mL;
  while( std::getline(ifs, line) ) {
    if(std::regex_match(line, match, comment_re) ) continue;
    std::istringstream iss(line);
    int run;
    double L;
    iss >> run >> L;
    mL[run]=L;
    //std::cout << run << " " << L << std::endl;
  };

  double Ltot{0};
  for(auto & sp : SP) {
    double Lsum{0};
    for(auto run : sp.run_list) {
      auto it = mL.find(run);
      if(it!=mL.end()) {
        //std::cout << sp.name << "   " <<  run << "  " <<  it->first << "  " << it->second << std::endl;
        Lsum+=it->second;
      }
    }
    printf("%20s %20.3f MeV %20.5f pb^-1\n", sp.title.c_str(), sp.energy/MeV, Lsum*1e-3);
    //std::cout << sp.title << "  " << sp.energy/MeV << "          " <<  Lsum*1e-3 << " pb^-1"  << std::endl;
    Ltot+=Lsum;
  }
    printf("%20s %20s     %20.5f pb^-1\n", "Total", "", Ltot*1e-3);

};


void lum_mc_data_compare(Scan_t & data = DATA, Scan_t & G = GG, Scan_t &  B = BB, int Nbin=100, int Nmax=0) {
  int point = 5;
  auto & b = B[point];
  auto & g = G[point];
  auto & d = data[point];
  compare(b , d , &ScanPoint_t::bb , "E_Eb"       , LOCAL_BB_SEL && "Nc==2" , {.file_name="lum_bb_data_mc_cmp_Eeb.pdf", .title="E/E_{beam}"    , .xaxis_title="E/E_{beam}, rad"    ,        .Nbin=Nbin, .Nmax=Nmax});
  compare(b , d , &ScanPoint_t::bb , "dtheta"     , LOCAL_BB_SEL            , {.file_name="lum_bb_data_mc_cmp_dtheta.pdf", .title="#Delta #theta" , .xaxis_title="#Delta #theta, rad" ,     .Nbin=Nbin, .Nmax=Nmax});
  compare(b , d , &ScanPoint_t::bb , "dphi"       , LOCAL_BB_SEL            , {.file_name="lum_bb_data_mc_cmp_dphi.pdf", .title="#Delta #phi"   , .xaxis_title="#Delta #phi, rad"   ,       .Nbin=Nbin, .Nmax=Nmax});
  compare(b , d , &ScanPoint_t::bb , "cos(theta)" , LOCAL_BB_SEL && "Nc==2" , {.file_name="lum_bb_data_mc_cmp_cos_theta.pdf", .title="cos(#theta)"   , .xaxis_title="cos #theta"    ,       .Nbin=Nbin, .Nmax=Nmax});

  compare(g , d , &ScanPoint_t::gg , "E_Eb"       , LOCAL_GG_SEL && "N0==2" , {.file_name="lum_gg_data_mc_cmp_Eeb.pdf"       , .title="E/E_{beam}"    , .xaxis_title="E/E_{beam}"    ,      .Nbin=Nbin, .Nmax=Nmax});
  compare(g , d , &ScanPoint_t::gg , "dtheta"     , LOCAL_GG_SEL && "N0==2" , {.file_name="lum_gg_data_mc_cmp_dtheta.pdf"    , .title="#Delta #theta" , .xaxis_title="#Delta #theta, rad" , .Nbin=Nbin, .Nmax=Nmax});
  compare(g , d , &ScanPoint_t::gg , "dphi"       , LOCAL_GG_SEL && "N0==2" , {.file_name="lum_gg_data_mc_cmp_dphi.pdf"      , .title="#Delta #phi"   , .xaxis_title="#Delta #phi, rad"   , .Nbin=Nbin, .Nmax=Nmax});
  compare(g , d , &ScanPoint_t::gg , "cos(theta)" , LOCAL_GG_SEL && "N0==2" , {.file_name="lum_gg_data_mc_cmp_cos_theta.pdf" , .title="cos(#theta)"   , .xaxis_title="cos #theta"    ,      .Nbin=Nbin, .Nmax=Nmax});
}
void mh_mc_data_compare(int Nbin=100) {
  auto & mc = JPSIMC[0];
  auto & data = JPSI[0];
  compare(mc, data, &ScanPoint_t::tt , "S", LOCAL_MH_SEL , {.file_name="mh_S_mc_data_cmp.pdf", .title="Sphericity"    , .xaxis_title="S", .Nbin=Nbin});
  compare(mc, data, &ScanPoint_t::tt , "Nchc", LOCAL_MH_SEL , {.file_name="mh_Nq_mc_data_cmp.pdf", .title="Multiplicity"    , .xaxis_title="N_{q}", .Nbin=Nbin});
  compare(mc, data, &ScanPoint_t::tt , "cos(rtheta)", LOCAL_MH_SEL , {.file_name="mh_cos_theta_mc_data_cmp.pdf", .title="cos(theta)"    , .xaxis_title="cos #theta", .Nbin=Nbin});
  compare(mc, data, &ScanPoint_t::tt , "rp", LOCAL_MH_SEL , {.file_name="mh_p_mc_data_cmp.pdf", .title="p"    , .xaxis_title="p", .Nbin=Nbin});
}

#include "/home/nikolaev/work/ibn/averager.h"
void calc_cbs_spread(const Scan_t & S) {
  new TCanvas;
  ibn::phys_averager A;
  double sum{0};
  TGraphErrors * g = new TGraphErrors;
  int i=0;
  for(const auto & sp : S) {
    std::cout << sp.energy << " " << sp.energy_spread*1e6 << " " << sp.energy_spread.error*1e6 << std::endl;
    sum+=sp.energy_spread*1e6;
    A.add(sp.energy_spread*1e6, sp.energy_spread.error*1e6);
    g->SetPoint(i, sp.energy, sp.energy_spread*1e3);
    g->SetPointError(i, 0, sp.energy_spread.error*1e3);
    ++i;
  }
  std::cout << sum/S.size() << std::endl;
  std::cout << "chi2=" << A.chi2() << std::endl;
  std::cout << "Average spread: " << A.average() << " +-  " << A.sigma_average() << std::endl;
  g->Draw("a*");
  g->Fit("pol0");
}
