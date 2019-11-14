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

const char * runtable_name = "scan_points_ems3_privalov_lum.txt";
//const char * runtable_name = "scan_points_ems3.txt";
//const char * runtable_name = "scan_points_ems2_privalov_lum.txt";

auto RUNTABLE  = read_my_runtable(runtable_name);

//read data from directory "data". Runs are combined into points according to runtable 
auto DATA        = read_data("data", RUNTABLE);

//Monte Carlo simulation of the signal
auto MC          = read_mc("mc/signal", RUNTABLE);

//background GALUGA
std::map<std::string, Scan_t> GALUGA =
{
  {"ee"   , read_mc("mc/galuga/EEee"   , RUNTABLE)} ,
  {"uu"   , read_mc("mc/galuga/EEuu"   , RUNTABLE)} ,
  {"pipi" , read_mc("mc/galuga/EEpipi" , RUNTABLE)} ,
  {"KK"   , read_mc("mc/galuga/EEkk"   , RUNTABLE)}
};

//hadronic background
auto HADR       = read_mc("mc/hadrons", RUNTABLE);
//bhabha background
auto BB         = read_mc("mc/bhabha", RUNTABLE);

// for luminocity measurement
auto GG         = read_mc("mc/gg", RUNTABLE);


//std::string LOCAL_BB_SEL = "(acol-TMath::Pi())>-0.03";
std::string LOCAL_BB_SEL = "(acol-TMath::Pi())>-0.04 && abs(cos(theta[0])) < 0.8 && abs(cos(theta[1])) < 0.8";

//particla identification configuration
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
  }
};

//selection configuration
Selection SELCFG =
{
  "SEL8", //selection name
  //common cuts
  "Nc==2"
  "&& p[0] < 1.1 && p[1] < 1.1"
  "&& E[0]>0.1 && E[1]>0.1"
  "&& pt[0]>0.2 && pt[1]>0.2"
  "&& abs(cos(theta[0]))<0.93"
  "&& abs(cos(theta[1]))<0.93"
  "&& ( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )"
//  "&&  abs(cos_theta_mis2) < 0.6"
//  "&& Emis > 0.5 && Emis < 3.0"
//  "&& M2<3.5"
  , 
  PID,
  //kinematic selection for different channel
  { 
   //     {"eμ",   "Nn==0 && eu  && ptem>0.15*k_ptem"},
      {"eμ1",   "Nn==0 && eu  && ptem>0.15*k_ptem && acop>0.1"},
  //     {"eπ",   "Nn==0 && epi && ptem>0.2*k_ptem"},
      {"eπ1",   "Nn==0 && epi && ptem>0.15*k_ptem && acop>0.1"},
  //      ,{"eπ2",   "Nn==0 && epi && ptem>k_ptem*(0.15 +0.05*cos_theta_mis2^2) && Emis>1.5 && Emis<2.6"}
  //      ,{"eπ3",   "Nn==0 && epi"}
  //     ,{"eπ4",   "Nn==0 && epi && ptem>k_ptem*(0.15 +0.05*cos_theta_mis2^2)"}
  //      ,{"eπ5",   "Nn==0 && epi && ptem>0.2 && acop>0.1"}
  //      {"eρ",   "Nn==2 && Npi0 == 1 && erho"},
  //      {"eρ1",   "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*(0.15 +0.05*cos_theta_mis2^2)"},
  //      {"eρ2",   "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.15"},
  //      {"eρ3",   "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.15 && acop>0.1"},
  //      {"eρ4",   "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.17"},
          {"eρ5",   "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.15 &&  abs(cos_theta_mis2)<0.8 && acop>0.1"},
  //        {"μπ0",   "Nn==0 && upi"},
  //        {"μπ1",   "Nn==0 && upi && ptem>0.3*k_ptem"},
  //        {"μπ2",   "Nn==0 && upi && ptem>0.15*k_ptem"},
  //        {"μπ3",   "Nn==0 && upi && ptem>0.15*k_ptem &&  abs(cos_theta_mis2)<0.8 && abs(cos(theta[0]))<0.8 && abs(cos(theta[1]))<0.8"},
  //        {"μπ4",   "Nn==0 && upi && ptem>0.15*k_ptem*(1+cos_theta_mis2*cos_theta_mis2)"},
  //        {"μπ5",   "Nn==0 && upi && ptem>0.21*k_ptem*(1+0.5*cos_theta_mis2*cos_theta_mis2)"},
          {"μπ6",   "Nn==0 && upi && ptem>0.25*k_ptem && acop>0.1"},
  //  {"μρ",   "Nn==2 && Npi0 == 1 && urho"},
  //  {"μρ1",   "Nn==2 && Npi0 == 1 && urho && ptem>0.1*k_ptem"},
  //  {"μρ2",   "Nn==2 && Npi0 == 1 && urho && ptem>0.15*k_ptem"},
  //  {"μρ3",   "Nn==2 && Npi0 == 1 && urho && ptem>0.20*k_ptem"},
  //  {"μρ4",   "Nn==2 && Npi0 == 1 && urho && ptem>0.25*k_ptem"},
    {"μρ5",   "Nn==2 && Npi0 == 1 && urho && ptem>0.25*k_ptem && acop>0.1"},
 //    {"ee",   "Nn==0 && ee && ptem > 0.3*k_ptem && M2>1.0 && barrel && abs(cos_theta_mis2)<0.8"},
//     {"ee",   "Nn==0 && ee"},
   //  {"ee0",   "Nn==0 && ee && ptem > 0.1"},
 //    {"ee1",   "Nn==0 && ee && ptem > 0.15"},
  //   {"ee2",   "Nn==0 && ee && ptem > 0.2"},
   //  {"ee3",   "Nn==0 && ee && ptem > 0.25"},
//     {"ee4",   "Nn==0 && ee && ptem > 0.3 +1.1*cos_theta_mis2^2"},
    // {"ee5",   "Nn==0 && ee && ptem > 0.25*(1.0 + 0.5*cos_theta_mis2*cos_theta_mis2)"},
     //{"ee6",   "Nn==0 && ee && ptem > 0.3 +0.3*cos_theta_mis2^2 && acop>0.1"},
     //{"ee7",   "Nn==0 && ee && ptem > 0.3 +0.44*cos_theta_mis2^2 && acop>0.1"},
     //{"ee8",   "Nn==0 && ee && Emis>1.8 && ptem > 0.32 +0.3*cos_theta_mis2^2 && acop>0.1"},
     {"ee9",   "Nn==0 && ee && Emis>1.8 && ptem > k_ptem*0.34*(1+cos_theta_mis2^2) && acop>0.1"},
 //   {"μμ0",   "Nn==0 && uu"},
 //   {"μμ1",   "Nn==0 && uu && ptem > 0.1*k_ptem"},
 //   {"μμ2",   "Nn==0 && uu && ptem > 0.15*k_ptem"},
 //   {"μμ3",   "Nn==0 && uu && ptem > 0.2*k_ptem && acop>0.1"},
  //  {"μμ4",   "Nn==0 && uu && ptem > 0.25*k_ptem"},
  //  {"μμ5",   "Nn==0 && uu && ptem > 0.15*k_ptem &&  barrel && abs(cos_theta_mis2)<0.8"},
  //  {"μμ6",   "Nn==0 && uu && ptem > 0.15*k_ptem*(1+cos_theta_mis2^2)"},
  //  {"μμ7",   "Nn==0 && uu && ptem > 0.18*k_ptem*(1+cos_theta_mis2^2)"},
  //  {"μμ8",   "Nn==0 && uu && ptem > 0.18*k_ptem*(1+cos_theta_mis2^2)"},
//    {"μμ9",   "Nn==0 && uu && ptem > k_ptem*(0.18 + 0.1*cos_theta_mis2^2) && acop>0.1 "},
//     {"μμ",   "Nn==0 && uu && ptem > 0.2*k_ptem && Emis>1.1 && abs(cos_theta_mis2)<0.8"},
//    {"μμ10",   "Nn==0 && uu && ptem>0.1 && acop>0.1"},
//    {"μμ11",   "Nn==0 && uu && E[0]>0.15 && E[0]<0.25 && E[1]>0.15 && E[1]<0.25 && ptem>0.2"},
      {"μμ12",   "Nn==0 && uu && ptem>0.25*k_ptem"},
 //   // ,{"ππo",   "Nn==0 && pipi&& ptem > 0.25*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"}
 //    {"ππ",   "Nn==0 && pipi "},
     {"ππ",   "Nn==0 && pipi && ptem>0.15*k_ptem &&  1.7 < Emis && Emis < 1.95"},
     //{"πρo",   "Nn==2 && Npi0 == 1 && pirho && ptem>0.4*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8"},
     //{"πρ",   "Nn==2 && Npi0 == 1 && pirho"},
     //{"πρ2",   "Nn==2 && Npi0 == 1 && pirho && ptem>0.4*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8 && Emis >1.4 && Emis<2.0"},
     {"πρ",   "Nn==2 && Npi0 == 1 && pirho && acop > 1.2 && Emis>1.5 && Emis <1.8"},

//    {"ρρ",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14  && ((0.6 <Mrho[0] && Mrho[0]<1.0 && 0.6<Mrho[3] && Mrho[3]<1.0)  || (0.6 <Mrho[1] && Mrho[1]<1.0 && 0.6<Mrho[2] && Mrho[2]<1.0))  && ptem>0.4*k_ptem && M2pi>1.0 && barrel && abs(cos_theta_mis2)<0.8 && 0.12 < Mpi0[1] && Mpi0[1] < 0.14 && Emis>1.35"}
  //  {"ρρ",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14  && ((0.6 <Mrho[0] && Mrho[0]<1.0 && 0.6<Mrho[3] && Mrho[3]<1.0)  || (0.6 <Mrho[1] && Mrho[1]<1.0 && 0.6<Mrho[2] && Mrho[2]<1.0))"},
  //  {"ρρ0",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14"},
  //  {"ρρ1",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14  && ((0.6 <Mrho[0] && Mrho[0]<1.0 && 0.6<Mrho[3] && Mrho[3]<1.0))"},
  //  {"ρρ2",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14  && ((0.6 <Mrho[1] && Mrho[1]<1.0 && 0.6<Mrho[2] && Mrho[2]<1.0)) && 0.12 < Mpi0[1] && Mpi0[1] < 0.14"},
    //{"ρρ3",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14 &&  abs((Mrho[0]-0.767)/0.15)<1 && abs((Mrho[3]-0.767)/0.15)<1 && ptem>k_ptem*0.5*(1.0+cos_theta_mis2^2)"},
    //{"ρρ4",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14 &&  abs((Mrho[1]-0.767)/0.15)<1 && abs((Mrho[2]-0.767)/0.15)<1 && ptem>k_ptem*0.5*(1.0+cos_theta_mis2^2)"},
    //{"ρρ5",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14 &&  abs((Mrho[0]-0.767)/0.15)<1 && abs((Mrho[3]-0.767)/0.15)<1"},
    //{"ρρ6",   "Nn==4 && Npi0==2 && PI0 && PI1 && 0.12 < Mpi0[0] && Mpi0[0] < 0.14  && 0.12 < Mpi0[1] && Mpi0[1] < 0.14 &&  abs((Mrho[1]-0.767)/0.15)<1 && abs((Mrho[2]-0.767)/0.15)<1"},
  },
};

Selection SEL8 =
{
  "SEL8", //selection name
  //common cuts
  "Nc==2"
 // "&& p[0] < 1.1 && p[1] < 1.1"
  "&& p[0] < 1.05 && p[1] < 1.05"
  "&& E[0]>0.025 && E[1]>0.025"
  "&& pt[0]>0.2 && pt[1]>0.2"
  //"&& abs(cos(theta[0]))<0.93"
  //"&& abs(cos(theta[1]))<0.93"
  "&& abs(cos(theta[0]))<0.8"
  "&& abs(cos(theta[1]))<0.8"
  "&& ( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )"
  , 
  PID,
  //kinematic selection for different channel
  { 
    {"eμ",  "Nn==0 && eu  && ptem>0.15*k_ptem && acop>0.1"},
    {"eπ",  "Nn==0 && epi && ptem>0.15*k_ptem && acop>0.1"},
    {"eρ",  "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.15 &&  abs(cos_theta_mis2)<0.8 && acop>0.1"},
    {"μπ",  "Nn==0 && upi && ptem>0.25*k_ptem && acop>0.1"},
    {"μρ",  "Nn==2 && Npi0 == 1 && urho && ptem>0.25*k_ptem && acop>0.1"},
    {"ee",  "Nn==0 && ee && Emis>1.8 && ptem > k_ptem*0.34*(1+cos_theta_mis2^2) && acop>0.1"},
    {"μμ",  "Nn==0 && uu && ptem>0.25*k_ptem"},
    {"ππ",  "Nn==0 && pipi && ptem>0.15*k_ptem &&  1.7 < Emis && Emis < 1.95"},
    {"πρ",  "Nn==2 && Npi0 == 1 && pirho && acop > 1.2 && Emis>1.5 && Emis <1.8 && ptem>0.15*k_ptem"},
  },
};

Selection SELRHO =
{
  "SEL8", //selection name
  //common cuts
  "Nc==2"
 // "&& p[0] < 1.1 && p[1] < 1.1"
  "&& p[0] < 1.05 && p[1] < 1.05"
  "&& E[0]>0.025 && E[1]>0.025"
  "&& pt[0]>0.2 && pt[1]>0.2"
  //"&& abs(cos(theta[0]))<0.93"
  //"&& abs(cos(theta[1]))<0.93"
  "&& abs(cos(theta[0]))<0.8"
  "&& abs(cos(theta[1]))<0.8"
  "&& ( abs(cos_theta_mis2) < 0.8 || ( 0.92 > abs(cos_theta_mis2) && abs(cos_theta_mis2) > 0.86) )"
  , 
  PID,
  //kinematic selection for different channel
  { 
    {"eρ",  "Nn==2 && Npi0 == 1 && erho &&  ptem>k_ptem*0.15 &&  abs(cos_theta_mis2)<0.8 && acop>0.1"},
    {"μρ",  "Nn==2 && Npi0 == 1 && urho && ptem>0.25*k_ptem && acop>0.1"},
    {"πρ",  "Nn==2 && Npi0 == 1 && pirho && acop > 1.2 && Emis>1.5 && Emis <1.8 && ptem>0.15*k_ptem"},
  },
};


auto & SEL = SEL8;

void doall(Selection & S=SEL, double kptem=1.0, Scan_t & D = DATA/* data */ , Scan_t & M = MC/* signal Monte Carlo */, std::string name="sel" /* the name of do all */, std::vector<int> skip_list = {})
{
  set_kptem(D,kptem); 
  auto result = new_select(D,S); 

  set_kptem(M,kptem);
  set_efficiency(result,M,1000000);
  fit(result, name + ".txt",  S.name);
  sleep(10);
  make_tex(print_tex(result,S.name, name + "_fit.pdf"),name + ".tex");
};

void cmpall(Selection &S=SEL) {
  //cmp({MC , DATA} , S , "ptem"           , "" , "NORM" , 40 , 0   , 1.1);
  //cmp({MC , DATA} , S , "cos_theta_mis2" , "" , "NORM" , 40 , -1  , +1);
  //cmp({MC , DATA} , S , "acop"           , "" , "NORM" , 40 , 0   , TMath::Pi());
  //cmp({MC , DATA} , S , "acol"           , "" , "NORM" , 40 , 0   , 1);
  //cmp({MC , DATA} , S , "p"              , "" , "NORM" , 40 , 0   , 1.1);
  //cmp({MC , DATA} , S , "Emis"           , "" , "NORM" , 40 , 1.0 , 3.5);
  //cmp({MC , DATA} , S , "cos(theta)"     , "" , "NORM" , 40 , -1  , -1);
  cmp({MC , DATA} , S , "Mpi0"     , "" , "NORM" , 40 , 0.11  , 0.15);
  cmp({MC , DATA} , S , "Mrho"     , "" , "NORM" , 40 , 0.5  , 1.1);
}


void select() 
{
  BB_SEL = LOCAL_BB_SEL;
  double Kptem = 1.0;
  set_pid_kptem(DATA       , PID , Kptem);
  set_pid_kptem(MC         , PID , Kptem);
  set_pid_kptem(HADR       , PID , Kptem);
  set_pid_kptem(BB         , PID , Kptem);
  set_pid_kptem(GG         , PID , Kptem);
  for( auto & p: GALUGA) {
    set_pid_kptem(p.second , PID , Kptem);
  }
}

void parameter_example(Selection & S=SEL)
{
  fold_and_draw(DATA,"ptem","Nc==2 && Nn==0 &&" + S.common_cut,"NORM");
  fold_and_draw(MC,"ptem","Nc==2 && Nn==0 &&" + S.common_cut,"NORM SAME");
}

