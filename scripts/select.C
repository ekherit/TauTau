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

#include "selection.C"
#include <map>

auto DATA        = read_data("data", read_my_runtable("../scan_points.txt"));
//auto DATAp       = read_data("data-old", read_my_runtable("../scan_points.txt"));
//auto DATA5       = read_data("../sel5/data", read_my_runtable("../scan_points.txt"));
//auto DATA703     = read_data("data703", read_my_runtable("../scan_points.txt"));
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
        restrict("Ep",         0.8,  1.1), 
        restrict("chi2_dedx",   0 ,  3.0),  
        restrict("delta_tof", -0.3,  0.3) 
      }}, 
    {"u", { "!e#" 
           ,"depth[#]-p[#]*58.61 > -35"
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
 //     ,restrict("p", 0.8, 1.05) 
      }}, 
    {"K", { "!(e#||u#||pi#)" },
      { 
        restrict("Ep",           0, 0.6), 
        restrict("chi2_dedx",    0, 3.0),  
        restrict("delta_tof", -0.3, 0.3)
//      ,restrict("p", 0.8, 1.05) 
      }},
 //   {"rho", { "pi# && 0.1128 < Mpi0[0] && Mpi0[0] < 0.1464 && 0.6 < Mrho[#] && Mrho[#]<1.0"},
      {"rho", { "pi# && 0.12 < Mpi0[0] && Mpi0[0] < 0.14 && 0.6 < Mrho[#] && Mrho[#]<1.0"},
      { 
      }}, 
    {"L", { "(e# || u#)"},
      { 
      }} 
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
  PID,
  //kinematic selection for different channel
  {  {"e-μ+",   "Nn==0 && e0 && u1   && ptem > 0.2"}
    ,{"e+μ-",   "Nn==0 && e1 && u0   && ptem > 0.2"}
    ,{"e-π+",   "Nn==0 && e0 && pi1  && ptem > 0.2 && MM2>1"}
    ,{"e+π-",   "Nn==0 && e1 && pi0  && ptem > 0.3 && MM2>1"}
    ,{"ρ-e+",   "Nn==2 && Npi0 == 1 &&  rho0 && e1  && ptem>0.2"}
    ,{"ρ+e-",   "Nn==2 && Npi0 == 1 &&  rho1 && e0  && ptem>0.05"}
    ,{"ρ-μ+",   "Nn==2 && Npi0 == 1 &&  rho0 && u1 && !rho1 && ptem>0.3 && M2>1.0"}
    ,{"ρ+μ-",   "Nn==2 && Npi0 == 1 &&  rho1 && u0 && !rho0 && ptem>0.05"}
    ,{"μπ",     "Nn==0 && upi  && ptem > 0.3 && M2>1.0"}
//  ,{"ρπ",    "Nn==2 && Npi0 == 1 &&  ((rho0 && pi1 && !rho1) || (rho1 && pi0 && !rho0)) && ptem>0.4 && MM<3.0"}
    ,{"e+e-",   "Nn==0 && ee  &&  ptem>0.5 && MM2>2"}
    ,{"μ+μ-",   "Nn==0 && uu   && ptem > 0.3 && abs(cos_theta_mis2)<0.8 && MM2>1.0"}
    ,{"π+π-",   "pipi && ptem > 0.5 &&  MM2>2.0 && M2>2.0"}
    ,{"ρ+ρ-",   "Nn==4 && Npi0==2  && rho0 && rho1 && ptem > 0.2 && MM2>1.0 && M2>1.0"}
  }
};

Selection SEL2 =
{
  "SEL2", //selection name
  //common cuts
  "Nc==2"
  "&& p[0] < 1.00 && p[1] < 1.00"
  "&& E[0]>0.05 && E[1]>0.05"
  "&& pt[0]>0.2 && pt[1]>0.2"
  "&& M2<3.0"
  //"&& cos(theta[0])<0.8 && cos(theta[1])<0.8"
  //"&& abs(cos_theta_mis2)<0.8"
  //"&& MM<3.5"
  //"&& ptem>0.3"
  //"&& acop>0.1"
//  "&& MM2>1.0"
//  "&& M2>1.0"
  , 
  PID,
  //kinematic selection for different channel
  { 
     {"eμ",   "Nn==0 && eu  && ptem>0.15"}
    ,{"eπ",   "Nn==0 && epi && ptem>0.25"}
    //{"e!e",  "Nn==0 && ((e0 && (!e1 && Ep[1]<0.7)) || (e1 && (!e0 && Ep[0]<0.7)))"}
    ,{"eρ",    "Nn==2 && Npi0 == 1 && erho && ptem>0.15"}
    ,{"μπ",    "Nn==0 && upi && ptem>0.3"}
////    ,{"μ-π+",    "Nn==0 && (u0 && !e1) && (E[1]>0.3 || E[1]<0.1)"}
////    ,{"π-μ+",    "Nn==0 && (u1 && !e0) && (E[0]>0.3 || E[0]<0.1)"}
//    ,{"μρ",    "Nn==2 && Npi0 ==1 && urho"}
////    ,{"Xrho",  "Nn==2 && Npi0 == 1 && Xrho"}
////  ,{"ρπ",    "Nn==2 && Npi0 == 1 &&  ((rho0 && pi1 && !rho1) || (rho1 && pi0 && !rho0)) && ptem>0.4 && MM<3.0"}
    ,{"ee", "Nn==0 && ee       && ptem > 0.3 && abs(cos_theta_mis2)<0.8 && acop<2.8"}
    ,{"μμ", "Nn==0 && uu       && ptem > 0.2 && abs(cos_theta_mis2)<0.8 && acop<2.8"}
    ,{"ππ", "Nn==0 && pipi     && ptem > 0.3 &&  M2>1.0 && abs(cos_theta_mis2) < 0.8 && acop<2.8"}
    ,{"ρρ", "Nn==4 && Npi0==2  && rhorho && ptem>0.3 && M2>1.0 && abs(cos_theta_mis2)<0.8 && acop<2.8"}
//    //,{"!e!e",   "Nn==0 && (!e1 && Ep[1]<0.7) && (!e0 && Ep[0]<0.7)"}
  }
};

Selection SEL3 =
{
  "SEL3", //selection name
  //common cuts
  "Nc==2"
  "&& p[0] < 1.00 && p[1] < 1.00"
  "&& E[0]>0.025 && E[1]>0.025"
  "&& pt[0]>0.2 && pt[1]>0.2"
  "&& M2<3.5"
  "&& abs(cos(theta[0]))<0.93"
  "&& abs(cos(theta[1]))<0.93"
  "&& ( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )"
  "&& tof<4.5"
  , 
  PID,
  //kinematic selection for different channel
  { 
     {"eμ",   "Nn==0 && eu  && ptem>0.15*k_ptem"}
    ,{"eπ",   "Nn==0 && epi && ptem>0.20*k_ptem"}
   ,{"eρ",   "Nn==2 && Npi0 == 1 && erho && ptem>0.1*k_ptem"}
   ,{"μπ",   "Nn==0 && upi && ptem>0.4*k_ptem"}
   ,{"μρ",   "Nn==2 && Npi0 == 1 && urho && ptem>0.3*k_ptem"}
   ,{"ee",   "Nn==0 && ee && ptem > 0.4*k_ptem && M2>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
   ,{"μμ",   "Nn==0 && uu && ptem > 0.15*k_ptem &&  barrel && abs(cos_theta_mis2)<0.8"}
   ,{"ππ",   "Nn==0 && pipi&& ptem > 0.25*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
   ,{"ρρ",   "Nn==4 && Npi0==2  && rhorho && ptem>0.4*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
  }
};

Selection SEL_XX =
{
  "SEL3", //selection name
  //common cuts
  "Nc==2"
  "&& p[0] < 1.00 && p[1] < 1.00"
  "&& E[0]>0.05 && E[1]>0.05"
  "&& pt[0]>0.2 && pt[1]>0.2"
  "&& abs(cos(theta[0]))<0.93"
  "&& abs(cos(theta[1]))<0.93"
  "&& M2<3.5"
  "&& M2>0.5"
  "&& abs(cos_theta_mis2) < 0.8"
  "&& acop<2.8"
  , 
  PID,
  //kinematic selection for different channel
  { 
    // {"ee", "Nn==0 && ee       && ptem > 0.3"}
    {"μμ",  "Nn==0 && uu       && ptem > 0.2"}
    ,{"ππ", "Nn==0 && pipi     && ptem > 0.3"}
    ,{"ρρ", "Nn==4 && Npi0==2  && rhorho && ptem>0.3"}
  }
};

Selection SEL4 =
{
  "SEL3", //selection name
  //common cuts
  "Nc==2"
  "&& p[0] < 1.00 && p[1] < 1.00"
  "&& E[0]>0.025 && E[1]>0.025"
  "&& pt[0]>0.2 && pt[1]>0.2"
  "&& M2<3.5"
  "&& abs(cos(theta[0]))<0.93"
  "&& abs(cos(theta[1]))<0.93"
  //"&& ( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )"
  , 
  PID,
  //kinematic selection for different channel
  { 
     {"eμ",   "Nn==0 && eu  && ptem>0.3*k_ptem"}
    ,{"eπ",   "Nn==0 && epi && ptem>0.3*k_ptem"}
    ,{"eρ",   "Nn==2 && Npi0 == 1 && erho"}
    ,{"μπ",    "Nn==0 && upi && ptem>0.3*k_ptem"}
    ,{"μρ", "Nn==2 && Npi0 == 1 && urho && ptem>0.3*k_ptem"}
    ,{"ee", "Nn==0 && ee && ptem > 0.3*k_ptem && M2>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
    ,{"μμ", "Nn==0 && uu && ptem > 0.3*k_ptem &&  barrel && abs(cos_theta_mis2)<0.8"}
    ,{"ππ", "Nn==0 && pipi&& ptem > 0.3*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
    ,{"ρρ", "Nn==4 && Npi0==2  && rhorho && ptem>0.3*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
  }
};

Selection SEL5 =
{
  "SEL3", //selection name
  //common cuts
  "Nc==2"
  "&& p[0] < 1.00 && p[1] < 1.00"
  "&& E[0]>0.025 && E[1]>0.025"
  "&& pt[0]>0.2 && pt[1]>0.2"
//  "&& M2<3.5"
  "&& abs(cos(theta[0]))<0.93"
  "&& abs(cos(theta[1]))<0.93"
 // "&& ( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )"
  , 
  PID,
  //kinematic selection for different channel
  { 
     {"eμ",   "Nn==0 && eu  && ptem>0.2*k_ptem"}
    ,{"eπ",   "Nn==0 && epi && ptem>0.20*k_ptem"}
    ,{"eρ",   "Nn==2 && Npi0 == 1 && erho && ptem>0.2*k_ptem" }
    ,{"μπ",   "Nn==0 && upi && ptem>0.2*k_ptem"}
    ,{"μρ",   "Nn==2 && Npi0 == 1 && urho && ptem>0.2*k_ptem"}
    //,{"ee",   "Nn==0 && ee && ptem > 0.4*k_ptem && M2>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
    //,{"μμ",   "Nn==0 && uu && ptem > 0.15*k_ptem &&  barrel && abs(cos_theta_mis2)<0.8"}
    //,{"ππ",   "Nn==0 && pipi&& ptem > 0.25*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
    //,{"ρρ",   "Nn==4 && Npi0==2  && rhorho && ptem>0.4*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
  }
};


void select() 
{
  set_pid(DATA,PID);
  set_pid(MC,PID);
  set_pid(HADR704,PID);
  for( auto & p: GALUGA)  set_pid(p.second,PID);
}

