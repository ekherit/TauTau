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

auto DATA        = read_data("data", read_my_runtable("scan_points.txt"));
auto MC          = read_mc("mc/signal");

std::map<std::string, Scan_t> GALUGA =
{
  {"ee",   read_mc("mc/galuga/EEee")},
  {"uu",   read_mc("mc/galuga/EEuu")},
  {"pipi", read_mc("mc/galuga/EEpipi")},
  {"KK",   read_mc("mc/galuga/EEkk")}
};

auto HADR       = read_mc("mc/hadrons");

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


Selection SEL7 =
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
      {"eμ",   "Nn==0 && eu  && ptem>0.1*k_ptem"}
     ,{"eπ",   "Nn==0 && epi && ptem>0.20*k_ptem"}
     ,{"eρ",   "Nn==2 && Npi0 == 1 && erho && ptem>0.1*k_ptem"}
     ,{"μπ",   "Nn==0 && upi && ptem>0.4*k_ptem"}
     ,{"μρ",   "Nn==2 && Npi0 == 1 && urho && ptem>0.3*k_ptem"}
     ,{"ee",   "Nn==0 && ee && ptem > 0.4*k_ptem && M2>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
     ,{"μμ",   "Nn==0 && uu && ptem > 0.15*k_ptem &&  barrel && abs(cos_theta_mis2)<0.8"}
     ,{"ππ",   "Nn==0 && pipi&& ptem > 0.25*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
     ,{"πρ",   "Nn==2 && Npi0 == 1 && pirho && ptem>0.4*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
     ,{"ρρ",   "Nn==4 && Npi0==2  && rhorho && ptem>0.4*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8 && 0.12 < Mpi0[1] && Mpi0[1] < 0.14"}
     ,{"XX",   "Nn==4 && Npi0==2 && pipi && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14  && ((0.6 <Mrho[0] && Mrho[0]<1.0 && 0.6<Mrho[3] && Mrho[3]<1.0)  || (0.6 <Mrho[1] && Mrho[1]<1.0 && 0.6<Mrho[2] && Mrho[2]<1.0))  && ptem>0.4*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8 && 0.12 < Mpi0[1] && Mpi0[1] < 0.14"}
  }
};

void doall(Selection & S, double kptem=1.0, Scan_t & D = DATA/* data */ , Scan_t & M = MC/* signal Monte Carlo */, std::string name="sel7" /* the name of do all */)
{
  set_kptem(D,kptem); 
  auto result = new_select(D,S); 
  set_kptem(M,kptem);
  set_efficiency(result,M,1000000);
  fit(result, name + ".txt",  S.name);
  sleep(10);
  make_tex(print_tex(result,S.name, name + "_fit.pdf"),name + ".tex");
};

void select7() 
{
  set_pid(DATA,PID);
  set_kptem(DATA,1);
  set_pid(MC,PID);
  set_kptem(MC,1);
  set_pid(HADR,PID);
  set_kptem(HADR,1);
  for( auto & p: GALUGA) 
  { 
    set_kptem(p.second,1);
    set_pid(p.second,PID);
  }
}

