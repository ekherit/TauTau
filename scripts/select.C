/*
 * =====================================================================================
 *
 *       Filename:  select.C
 *
 *    Description:  Selection configuration of tau events
 *
 *        Version:  1.0
 *        Created:  10.02.2019 14:27:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#include "../selection.C"
#include <map>

auto DATA        = read_data("data", read_my_runtable("../scan_points.txt"));
auto DATAp       = read_data("data-old", read_my_runtable("../scan_points.txt"));
auto DATA5       = read_data("../sel5/data", read_my_runtable("../scan_points.txt"));
auto DATA703     = read_data("data703", read_my_runtable("../scan_points.txt"));
auto DATA11      = read_data3("../tau2011","../tau2011/runtable.txt");
auto MC          = read_mc("mc/signal");

std::map<std::string, Scan_t> GALUGA =
{
  {"ee",   read_mc("mc/EEee")},
  {"uu",   read_mc("mc/EEuu")},
  {"pipi", read_mc("mc/EEpipi")},
  {"KK",   read_mc("mc/EEkk")}
};

/*  bhabha */
auto BB          = read_mc("mc/BB");
auto UU          = read_mc("mc/uu");
auto HADR        = read_mc("mc/hadrons");
auto HADR704     = read_mc("mc/hadrons704");

std::vector<ParticleID_t> PID = 
{
  {
    {"e", {""}, 
      { 
        restrict("Ep",         0.8,  1.05), 
        restrict("chi2_dedx",   0 ,  2.0),  
        restrict("delta_tof", -0.2,  0.2) 
      }}, 
    {"u", { "!e#" },
      { 
        {"depth",        0, 100},
        restrict("E",          0.1, 0.3), 
        restrict("Ep",           0, 0.8), 
        restrict("chi2_dedx",    0, 5.0),  
        restrict("delta_tof", -0.3, 0.3) 
      }},
    {"pi", { "!(e#||u#)" },
      {
        restrict("Ep",           0, 0.6), 
        restrict("chi2_dedx",    0, 5.0),  
        restrict("delta_tof", -0.3, 0.3),
        restrict("p", 0.8, 1.05) 
      }}, 
    {"K", { "!(e#||u#||pi#)" },
      { 
        restrict("Ep",           0, 0.6), 
        restrict("chi2_dedx",    0, 3.0),  
        restrict("delta_tof", -0.3, 0.3),
        restrict("p", 0.8, 1.05) 
      }},
    {"rho", { "pi# && Mrho_cut#"},
      { 
        restrict("Mpi0",    0.1128, 0.1464), 
      }} 
  }
};

std::vector<ParticleID_t> PID2 = 
{
  {
    {"e", {""}, 
      { 
        restrict("Ep",         0.85,  1.05), 
        restrict("chi2_dedx",   0 ,  2.0),  
        restrict("delta_tof", -0.2,  0.2) 
      }}, 
    {"u", { "!e#" 
           ,"depth[#]-p[#]*58.61 > -35"
          },
      { 
        {"depth",        0, 100},
        restrict("E",          0.15, 0.25), 
        restrict("Ep",           0, 0.6), 
        restrict("chi2_dedx",    0, 5.0),  
        restrict("delta_tof", -0.3, 0.3) 
      }},
    {"pi", { "!(e#||u#)" },
      {
        restrict("Ep",           0, 0.6), 
        restrict("chi2_dedx",    0, 5.0),  
        restrict("delta_tof", -0.3, 0.3)
 //     ,restrict("p", 0.8, 1.05) 
      }}, 
    {"K", { "!(e#||u#||pi#)" },
      { 
        restrict("Ep",           0, 0.6), 
        restrict("chi2_dedx",    0, 3.0),  
        restrict("delta_tof", -0.3, 0.3)
//      ,restrict("p", 0.8, 1.05) 
      }},
    {"rho", { "pi# && 0.1128 < Mpi0[0] && Mpi0[0] < 0.1464 && 0.6 < Mrho[#] && Mrho[#]<1.0"},
      { 
      }}, 
    {"L", { "(e# || u#)"},
      { 
      }} 
  }
};

Selection SELECTION6 = {
  "SEL6", //selection name
  //extra cut
  "Nc==2 && Nn==0"
  "&& p[0] < 1.05 && p[1] < 1.05"
  "&& MM2>1.0"
  "&& abs(cos_theta_mis2)<0.8"
  "&& E[0]>0.05 && E[1]>0.05"
 // " && cos(theta[0])<0.8 && cos(theta[1])<0.8"
  , 
  //PID
  PID2,
  //kinematic selection for different channel
  {  {"eμ1", "e#mu"     , "e0 && u1   && ptem > 0.1"}
    ,{"eμ2", "e#mu"     , "e1 && u0   && ptem > 0.1"}
    ,{"eμ", "e#mu"      , "eu   && ptem > 0.1"}
    ,{"eπ", "e#pi"      , "epi  && ptem > 0.1"}
    ,{"e-π+", "e#pi"    , "e0 && pi1   && ptem > 0.1"}
    ,{"π-e+", "e#pi"    , "e1 && pi0   && ptem > 0.1"}
    ,{"ee", "ee"        , "ee   && ptem > 0.5 && (abs(cos_theta_mis2)<0.8 ||  abs(cos_theta_mis2) >0.86 && abs(cos_theta_mis2)<0.9) && MM2>1.0 && abs(cos(theta[0]))<0.8  && abs(cos(theta[1])) < 0.8"}
    ,{"μπ", "#mu#pi"    , "upi  && ptem > 0.5 && (abs(cos_theta_mis2)<0.8 ||  abs(cos_theta_mis2) >0.86 && abs(cos_theta_mis2)<0.9)"}
    ,{"μμ", "#mu#mu"    , "uu   && ptem > 0.3 && abs(cos_theta_mis2)<0.8 && MM2>1.0"}
    ,{"ππ", "#pi#pi"    , "pipi && ptem > 0.5 && abs(cos_theta_mis2)<0.8 && abs(cos(theta[0]))<0.8  && abs(cos(theta[1])) < 0.8 && MM2>1.0"}
    //,{"eK", "eK"       , "eK   && ptem > 0.0"},
    //,{"πK", "#piK"     , "piK  && ptem > 0.4 && abs(cos_theta_mis2)<0.8&& abs(cos(theta[0]))<0.8  && abs(cos(theta[1])) < 0.8"},
    //,{"μK", "#uK"      , "uK   && ptem > 0.4 && abs(cos_theta_mis2)<0.8&& abs(cos(theta[0]))<0.8  && abs(cos(theta[1])) < 0.8"},
    //,{"KK", "KK"       , "KK   && ptem > 0.1"}   
    //,{"K(π,μ,K)", "K(#pi,#mu,K)", "(KK || piK || uK) && ptem > 0.4 && abs(cos_theta_mis2)<0.8"}   
    ,{"Xρ", "X#rho"     , "Nc==2 && Npi0 ==1 && (Xrho || rhoX) && acop>1.4"} 
    ,{"ρρ", "#rho#rho"  , "Nc==2 && Npi0 ==2 && rhorho && ptem>0.5"}
  }
};

Selection SEL =
{
  "SEL", //selection name
  //common cuts
  "Nc==2"
  "&& p[0] < 1.07 && p[1] < 1.07"
  "&& E[0]>0.05 && E[1]>0.05"
  "&& cos(theta[0])<0.8 && cos(theta[1])<0.8"
  "&& abs(cos_theta_mis2)<0.8"
  "&& MM<3.5"
  , 
  //PID
  PID2,
  //kinematic selection for different channel
  {  {"e-μ+", "e^{-}#mu^{+}",    "Nn==0 && e0 && u1   && ptem > 0.2"}
    ,{"e+μ-", "e^{+}#mu^{-}",    "Nn==0 && e1 && u0   && ptem > 0.2"}
    ,{"e-π+", "e^{-}#pi^{+}",    "Nn==0 && e0 && pi1  && ptem > 0.2 && MM2>1"}
    ,{"e+π-", "e^{+}#pi^{-}",    "Nn==0 && e1 && pi0  && ptem > 0.3 && MM2>1"}
    ,{"ρ-e+","#rho^{-}e^{+}",    "Nn==2 && Npi0 == 1 &&  rho0 && e1  && ptem>0.2"}
    ,{"ρ+e-","#rho^{+}e^{-}",    "Nn==2 && Npi0 == 1 &&  rho1 && e0  && ptem>0.05"}
    ,{"ρ-μ+","#rho^{-}#mu^{+}",  "Nn==2 && Npi0 == 1 &&  rho0 && u1 && !rho1 && ptem>0.3 && M2>1.0"}
    ,{"ρ+μ-","#rho^{+}#mu^{-}",  "Nn==2 && Npi0 == 1 &&  rho1 && u0 && !rho0 && ptem>0.05"}
    ,{"μπ",     "#mu#pi"    ,        "Nn==0 && upi  && ptem > 0.3 && M2>1.0"}
//    ,{"ρπ","#rho^{-}#pi^{+}",    "Nn==2 && Npi0 == 1 &&  ((rho0 && pi1 && !rho1) || (rho1 && pi0 && !rho0)) && ptem>0.4 && MM<3.0"}
    ,{"e+e-", "e^{+}e^{-}"       ,"Nn==0 && ee  &&  ptem>0.5 && MM2>2"}
    ,{"μ+μ-", "#mu^{-}#mu^{+}",  "Nn==0 && uu   && ptem > 0.3 && abs(cos_theta_mis2)<0.8 && MM2>1.0"}
    ,{"π+π-", "#pi^{-}#pi^{+}"    , "pipi && ptem > 0.5 &&  MM2>2.0 && M2>2.0"}
    ,{"ρ+ρ-", "#rho^{-}#rho^{+}",  "Nn==4 && Npi0==2  && rho0 && rho1 && ptem > 0.2 && MM2>1.0 && M2>1.0"}
  }
};

Selection TEST =
{
  "TEST", //selection name
  //common cuts
  "Nc==2 && Nn==0"
  "&& p[0] < 1.05 && p[1] < 1.05"
  "&& cos(theta[0])<0.8 && cos(theta[1])<0.8"
  "&& abs(cos_theta_mis2)<0.8"
  , 
  //PID
  PID2,
  //kinematic selection for different channel
  {
    {"e0!e1", "e0!e1"     , "e0 && !e1  && Ep[1]<0.6 && ptem>0.1"},
    {"e1!e0", "e1!e0"     , "e1 && !e0  && Ep[0]<0.6 && ptem>0.1"},
  }
};


void select() 
{
  set_pid(DATA,PID2);
  set_pid(MC,PID2);
  for( auto & p: GALUGA)  set_pid(p.second,PID2);
}
