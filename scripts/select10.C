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
#include "sys.C"
#include <map>

std::string TAUFIT = "taufit --lum=default --tau-spread=1.258 --energy-correction=+0.011 --free-energy --free-luminosity --free-effcor";

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
std::string LOCAL_BB_SEL = "(acol-TMath::Pi())>-0.04 && abs(cos(theta[0])) < 0.8 && abs(cos(theta[1])) < 0.8 && Ep[0]>0.7 && Ep[1]>0.7";
std::string LOCAL_GG_SEL = "";

//particla identification configuration
std::vector<ParticleID_t> PID = 
{
  {
    {"e", {}, 
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
        restrict("Ep",           0, 0.7), 
      }},
  }
};

////selection configuration
//Selection_t SELCFG =
//{
//  "SEL8", //selection name
//  //common cuts
//  "Nc==2"
//  "&& p[0] < 1.1 && p[1] < 1.1"
//  "&& E[0]>0.1 && E[1]>0.1"
//  "&& pt[0]>0.2 && pt[1]>0.2"
//  "&& abs(cos(theta[0]))<0.93"
//  "&& abs(cos(theta[1]))<0.93"
//  "&& ( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )"
////  "&&  abs(cos_theta_mis2) < 0.6"
////  "&& Emis > 0.5 && Emis < 3.0"
////  "&& M2<3.5"
//  , 
//  PID,
//  //kinematic selection for different channel
//  { 
//   //     {"eμ",   "Nn==0 && eu  && ptem>0.15*k_ptem"},
//      {"eμ1",   "Nn==0 && eu  && ptem>0.15*k_ptem && acop>0.1"},
//  //     {"eπ",   "Nn==0 && epi && ptem>0.2*k_ptem"},
//      {"eπ1",   "Nn==0 && epi && ptem>0.15*k_ptem && acop>0.1"},
//  //      ,{"eπ2",   "Nn==0 && epi && ptem>k_ptem*(0.15 +0.05*cos_theta_mis2^2) && Emis>1.5 && Emis<2.6"}
//  //      ,{"eπ3",   "Nn==0 && epi"}
//  //     ,{"eπ4",   "Nn==0 && epi && ptem>k_ptem*(0.15 +0.05*cos_theta_mis2^2)"}
//  //      ,{"eπ5",   "Nn==0 && epi && ptem>0.2 && acop>0.1"}
//  //      {"eρ",   "Nn==2 && Npi0 == 1 && erho"},
//  //      {"eρ1",   "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*(0.15 +0.05*cos_theta_mis2^2)"},
//  //      {"eρ2",   "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.15"},
//  //      {"eρ3",   "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.15 && acop>0.1"},
//  //      {"eρ4",   "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.17"},
//          {"eρ5",   "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.15 &&  abs(cos_theta_mis2)<0.8 && acop>0.1"},
//  //        {"μπ0",   "Nn==0 && upi"},
//  //        {"μπ1",   "Nn==0 && upi && ptem>0.3*k_ptem"},
//  //        {"μπ2",   "Nn==0 && upi && ptem>0.15*k_ptem"},
//  //        {"μπ3",   "Nn==0 && upi && ptem>0.15*k_ptem &&  abs(cos_theta_mis2)<0.8 && abs(cos(theta[0]))<0.8 && abs(cos(theta[1]))<0.8"},
//  //        {"μπ4",   "Nn==0 && upi && ptem>0.15*k_ptem*(1+cos_theta_mis2*cos_theta_mis2)"},
//  //        {"μπ5",   "Nn==0 && upi && ptem>0.21*k_ptem*(1+0.5*cos_theta_mis2*cos_theta_mis2)"},
//          {"μπ6",   "Nn==0 && upi && ptem>0.25*k_ptem && acop>0.1"},
//  //  {"μρ",   "Nn==2 && Npi0 == 1 && urho"},
//  //  {"μρ1",   "Nn==2 && Npi0 == 1 && urho && ptem>0.1*k_ptem"},
//  //  {"μρ2",   "Nn==2 && Npi0 == 1 && urho && ptem>0.15*k_ptem"},
//  //  {"μρ3",   "Nn==2 && Npi0 == 1 && urho && ptem>0.20*k_ptem"},
//  //  {"μρ4",   "Nn==2 && Npi0 == 1 && urho && ptem>0.25*k_ptem"},
//    {"μρ5",   "Nn==2 && Npi0 == 1 && urho && ptem>0.25*k_ptem && acop>0.1"},
// //    {"ee",   "Nn==0 && ee && ptem > 0.3*k_ptem && M2>1.0 && barrel && abs(cos_theta_mis2)<0.8"},
////     {"ee",   "Nn==0 && ee"},
//   //  {"ee0",   "Nn==0 && ee && ptem > 0.1"},
// //    {"ee1",   "Nn==0 && ee && ptem > 0.15"},
//  //   {"ee2",   "Nn==0 && ee && ptem > 0.2"},
//   //  {"ee3",   "Nn==0 && ee && ptem > 0.25"},
////     {"ee4",   "Nn==0 && ee && ptem > 0.3 +1.1*cos_theta_mis2^2"},
//    // {"ee5",   "Nn==0 && ee && ptem > 0.25*(1.0 + 0.5*cos_theta_mis2*cos_theta_mis2)"},
//     //{"ee6",   "Nn==0 && ee && ptem > 0.3 +0.3*cos_theta_mis2^2 && acop>0.1"},
//     //{"ee7",   "Nn==0 && ee && ptem > 0.3 +0.44*cos_theta_mis2^2 && acop>0.1"},
//     //{"ee8",   "Nn==0 && ee && Emis>1.8 && ptem > 0.32 +0.3*cos_theta_mis2^2 && acop>0.1"},
//     {"ee9",   "Nn==0 && ee && Emis>1.8 && ptem > k_ptem*0.34*(1+cos_theta_mis2^2) && acop>0.1"},
// //   {"μμ0",   "Nn==0 && uu"},
// //   {"μμ1",   "Nn==0 && uu && ptem > 0.1*k_ptem"},
// //   {"μμ2",   "Nn==0 && uu && ptem > 0.15*k_ptem"},
// //   {"μμ3",   "Nn==0 && uu && ptem > 0.2*k_ptem && acop>0.1"},
//  //  {"μμ4",   "Nn==0 && uu && ptem > 0.25*k_ptem"},
//  //  {"μμ5",   "Nn==0 && uu && ptem > 0.15*k_ptem &&  barrel && abs(cos_theta_mis2)<0.8"},
//  //  {"μμ6",   "Nn==0 && uu && ptem > 0.15*k_ptem*(1+cos_theta_mis2^2)"},
//  //  {"μμ7",   "Nn==0 && uu && ptem > 0.18*k_ptem*(1+cos_theta_mis2^2)"},
//  //  {"μμ8",   "Nn==0 && uu && ptem > 0.18*k_ptem*(1+cos_theta_mis2^2)"},
////    {"μμ9",   "Nn==0 && uu && ptem > k_ptem*(0.18 + 0.1*cos_theta_mis2^2) && acop>0.1 "},
////     {"μμ",   "Nn==0 && uu && ptem > 0.2*k_ptem && Emis>1.1 && abs(cos_theta_mis2)<0.8"},
////    {"μμ10",   "Nn==0 && uu && ptem>0.1 && acop>0.1"},
////    {"μμ11",   "Nn==0 && uu && E[0]>0.15 && E[0]<0.25 && E[1]>0.15 && E[1]<0.25 && ptem>0.2"},
//      {"μμ12",   "Nn==0 && uu && ptem>0.25*k_ptem"},
// //   // ,{"ππo",   "Nn==0 && pipi&& ptem > 0.25*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
// //    {"ππ",   "Nn==0 && pipi "},
//     {"ππ",   "Nn==0 && pipi && ptem>0.15*k_ptem &&  1.7 < Emis && Emis < 1.95"},
//     //{"πρo",   "Nn==2 && Npi0 == 1 && pirho && ptem>0.4*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"},
//     //{"πρ",   "Nn==2 && Npi0 == 1 && pirho"},
//     //{"πρ2",   "Nn==2 && Npi0 == 1 && pirho && ptem>0.4*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8 && Emis >1.4 && Emis<2.0"},
//     {"πρ",   "Nn==2 && Npi0 == 1 && pirho && acop > 1.2 && Emis>1.5 && Emis <1.8"},
//
////    {"ρρ",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14  && ((0.6 <Mrho[0] && Mrho[0]<1.0 && 0.6<Mrho[3] && Mrho[3]<1.0)  || (0.6 <Mrho[1] && Mrho[1]<1.0 && 0.6<Mrho[2] && Mrho[2]<1.0))  && ptem>0.4*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8 && 0.12 < Mpi0[1] && Mpi0[1] < 0.14 && Emis>1.35"}
//  //  {"ρρ",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14  && ((0.6 <Mrho[0] && Mrho[0]<1.0 && 0.6<Mrho[3] && Mrho[3]<1.0)  || (0.6 <Mrho[1] && Mrho[1]<1.0 && 0.6<Mrho[2] && Mrho[2]<1.0))"},
//  //  {"ρρ0",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14"},
//  //  {"ρρ1",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14  && ((0.6 <Mrho[0] && Mrho[0]<1.0 && 0.6<Mrho[3] && Mrho[3]<1.0))"},
//  //  {"ρρ2",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14  && ((0.6 <Mrho[1] && Mrho[1]<1.0 && 0.6<Mrho[2] && Mrho[2]<1.0)) && 0.12 < Mpi0[1] && Mpi0[1] < 0.14"},
//    //{"ρρ3",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14 &&  abs((Mrho[0]-0.767)/0.15)<1 && abs((Mrho[3]-0.767)/0.15)<1 && ptem>k_ptem*0.5*(1.0+cos_theta_mis2^2)"},
//    //{"ρρ4",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14 &&  abs((Mrho[1]-0.767)/0.15)<1 && abs((Mrho[2]-0.767)/0.15)<1 && ptem>k_ptem*0.5*(1.0+cos_theta_mis2^2)"},
//    //{"ρρ5",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14 &&  abs((Mrho[0]-0.767)/0.15)<1 && abs((Mrho[3]-0.767)/0.15)<1"},
//    //{"ρρ6",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14 &&  abs((Mrho[1]-0.767)/0.15)<1 && abs((Mrho[2]-0.767)/0.15)<1"},
//  },
//};

Selection_t SEL8 =
{
  "SEL8", //selection name
  //common cuts
  "Nc==2"
  "&& p[0] < 1.1 && p[1] < 1.1"
//  "&& E[0]>0.025 && E[1]>0.025"
  "&& pt[0]>0.2 && pt[1]>0.2"
  "&& abs(cos(theta[0]))<0.8"
  "&& abs(cos(theta[1]))<0.8"
  "&& ( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )"
  , 
  PID,
  { //kinematic selection for different channel
/*  0 */ {"eμ",  "Nn==0 && eu  && ptem>0.15*k_ptem"},
/*  1 */ {"eπ",  "Nn==0 && epi && ptem>0.15*k_ptem"},
/*  2 */ {"eρ",  "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.15 &&  abs(cos_theta_mis2)<0.8 && acop>0.1"},
/*  3 */ {"μπ",  "Nn==0 && upi && ptem>0.25*k_ptem && acop>0.1"},
/*  4 */ {"μρ",  "Nn==2 && Npi0 == 1 && urho && ptem>0.25*k_ptem && acop>0.1"},
/*  5 */ {"ee",  "Nn==0 && ee && Emis>1.8 && ptem > k_ptem*0.34*(1+cos_theta_mis2^2) && acop>0.1"},
/*  6 */ {"μμ",  "Nn==0 && uu && ptem>0.2*k_ptem"},
/*  7 */ {"ππ",  "Nn==0 && pipi && ptem>0.15*k_ptem &&  1.7 < Emis && Emis < 1.95"},
/*  8 */ {"πρ",  "Nn==2 && Npi0 == 1 && pirho && acop > 1.2 && Emis>1.5 && Emis <1.8 && ptem>0.15*k_ptem"},
  },
};


Selection_t SELTEST =
{
  "SELTEST", //selection name
  //common cuts
  "Nc==2"
  "&& p[0] < 1.1 && p[1] < 1.1"
  "&& pt[0]>0.2 && pt[1]>0.2"
  "&& barrel"
  "&& ( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )"
  , 
  PID,
  { //kinematic selection for different channel
// {"eX",  "Nn==0 &&  eX   && ptem>0.21*(1+0.5*cos_theta_mis2^2)"},
 {"eX",  "Nn==0 &&  eX   && ptem>0.23"},
// {"eX",  "Nn==0 &&  eX   && ptem>0.15 && Enmax<0.001"},
/*  0 */ //{"uX",  "Nn==0 &&  uX   && ptem>0.3 && Emis<2.2 && abs(cos_theta_mis2) < 0.6"},
    //{"eu",  "Nn==0 &&  eu   && ptem>0.1*(1+0.5*cos_theta_mis2^2) "},
//    {"eu",  "Nn==0 &&  eu   && ptem>0.15 "},
//    {"eπ",  "Nn==0 &&  epi   && ptem>0.15"},
//    {"eX*",  "Nn==0 && eX && !eu && !epi && ptem>0.25 && acol>0.05"},
//     {"eρ",  "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.15 &&  abs(cos_theta_mis2)<0.8"},
//         {"μρ",  "Nn==2 && Npi0 == 1 && urho && 1.4<Emis && Emis < 2.2"},
//          {"μπ",  "Nn==0 && upi && ptem>0.25*k_ptem && acop>0.1"},
//         {"ee",  "Nn==0 && ee && Emis>1.8 && ptem > k_ptem*0.34*(1+cos_theta_mis2^2) && acop>0.1"},
//         {"uu",  "Nn==0 &&  uu  && ptem>0.18*(1+0.2*cos_theta_mis2^2)"}, 
//         {"ππ","Nn==0 && pipi && ptem>0.18 && Emis>1.6 && Emis<2.0"},
//         {"πρ",  "Nn==2 && Npi0 == 1 && pirho && acop > 1.2 && Emis>1.5 && Emis <1.8"},
  },
};

//Selection_t SELRHO =
//{
//  "SEL8", //selection name
//  //common cuts
//  "Nc==2"
// // "&& p[0] < 1.1 && p[1] < 1.1"
//  "&& p[0] < 1.05 && p[1] < 1.05"
//  "&& E[0]>0.025 && E[1]>0.025"
//  "&& pt[0]>0.2 && pt[1]>0.2"
//  //"&& abs(cos(theta[0]))<0.93"
//  //"&& abs(cos(theta[1]))<0.93"
//  "&& abs(cos(theta[0]))<0.8"
//  "&& abs(cos(theta[1]))<0.8"
//  "&& ( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )"
//  , 
//  PID,
//  //kinematic selection for different channel
//  { 
//    {"eρ",  "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.15 &&  abs(cos_theta_mis2)<0.8 && acop>0.1"},
//    {"μρ",  "Nn==2 && Npi0 == 1 && urho && ptem>0.25*k_ptem && acop>0.1"},
//    {"πρ",  "Nn==2 && Npi0 == 1 && pirho && acop > 1.2 && Emis>1.5 && Emis <1.8 && ptem>0.15*k_ptem"},
//  },
//};

auto & SEL = SEL8;
std::vector<SelRef_t> SELS = {SEL8,SELTEST};

//void doall(Selection_t & S=SEL, double kptem=1.0, Scan_t & D = DATA/* data */ , Scan_t & M = SIGNAL/* signal Monte Carlo */, std::string name="sel" /* the name of do all */, std::string default_lum="gg")
void doall(Scan_t & D = DATA/* data */ , Selection_t & S=SEL, double kptem=1.0,  Scan_t & M = SIGNAL/* signal Monte Carlo */, std::string name="sel" /* the name of do all */, std::string default_lum="gg")
{
  S.apply_common_cut();
  set_kptem(D,kptem); 
  measure_luminosity(D, BB, GG,1e6);
  auto result = select(D,S); 
  print_luminosity(D);
  set_kptem(M,kptem);
  set_efficiency(result,M,1000000);
  fit(result, name + ".txt",  S.title, default_lum);
  sleep(10);
  make_tex(print_tex(result,S.title, name + "_fit.pdf"),name + ".tex");
};

void cmpall(Selection_t &S=SEL) {
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



void select() 
{
  BB_SEL = LOCAL_BB_SEL;
  GG_SEL = LOCAL_GG_SEL;
  TAUFIT_STR = TAUFIT;
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

void parameter_example(Selection_t & S=SEL)
{
  fold_and_draw(DATA,"ptem","Nc==2 && Nn==0 &&" + S.common_cut(),"NORM");
  fold_and_draw(SIGNAL,"ptem","Nc==2 && Nn==0 &&" + S.common_cut(),"NORM SAME");
}



TCanvas * draw_tt_ee_compare(std::string sig_cut, std::string bb_cut) {
  return draw_tt_ee_compare(SIGNAL, BB, sig_cut, bb_cut);
}

TCanvas * draw_tt_gg_compare(std::string cut1, std::string cut2) {
  return draw_tt_gg_compare(SIGNAL, GG, cut1, cut2);
}

TCanvas * draw_gg_ee_compare(std::string cut1, std::string cut2) {
  return draw_gg_ee_compare(GG, BB, cut1, cut2);
}


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
  long N = D[0].tt.tree->Draw("pid[0]:pid[1]", cut.c_str(),"goff");
  double * pid0 = D[0].tt.tree->GetV1();
  double * pid1 = D[0].tt.tree->GetV2();
  std::map<std::string, long> m;
  std::map<std::pair<int, int> , long> om;
  std::set<std::pair<long, std::string>, std::greater<std::pair<long, std::string>> > sm;
  std::set<std::pair<long, std::pair<int, int>>,  std::greater<std::pair<long, std::pair<int, int>>> > som;
  for(int i=0;i<N;++i) {
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

//TCanvas * draw_tt_gg_compare(std::string cut1, std::string cut2) {
//  return draw_tt_gg_compare(SIGNAL, GG, cut1, cut2);
//}
