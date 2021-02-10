/*
 * =====================================================================================
 *
 *       Filename:  analysis.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08.01.2021 20:02:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */
#include "analysis.h"

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
      }},
      {"rho", {" PI# && 0.12 < Mpi0[0] && Mpi0[0] < 0.14"},
      { 
        restrict("Mrho", 0.6,1.0),
      }}, 

      {"X", { "chi2_dedx_e[#]>3" 
          },
      { 
        restrict("Ep",           0, 0.8), 
      }},
  }
};

Analysis base(
    //WORKDIR
    "./",
    //runtable
    "all_scan_points_ems3.txt",
    //GG selection
    "",
    //BB selection
    "(acol-TMath::Pi())>-0.04 && abs(cos(theta[0])) < 0.8 && abs(cos(theta[1])) < 0.8 && Ep[0]>0.8 && Ep[1]>0.8 && abs(z[0])<10 && abs(z[1])<10 && vxy[0]<1.0 && vxy[1]<1.0", //bb_sel
    //MH selection
    "ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2",
    //PID
    PID,
    //TauTau selection common cuts
    { 
        "Ncg==2 && Ncg==Ncc",
        "abs(vz[0])<10 && abs(vz[1])<10",
        "vxy[0]<1.0 && vxy[1]<1.0",
        "E[0]>0.025 && E[1]>0.025",
        "q[0]==-1 && q[1]==1",
        "p[0] < 1.1 && p[1] < 1.1",
        "pt[0] > 0.2 && pt[1] > 0.2",
        "cgood",
        "2.5 < tof[0] && tof[0] < 5.5 && 2.5 < tof[1] && tof[1] < 5.5",
        "ptem50>0.25 && ptem50<1.1",
    }
    , 
    //channel selections
      { 
        {"eX",      "NnE50==0 && eX && missed_photon_angle"},
        {"eρ",     "Nng==2   && eX && (Mpi0[0] < 0.14 && Mpi0[0]>0.12) && good_emc_time" },
      }
    );


Analysis test(
    //WORKDIR
    "./",
    //runtable
    "all_scan_points_ems3.txt",
    //GG selection
    "",
    //BB selection
    "(acol-TMath::Pi())>-0.04 && abs(cos(theta[0])) < 0.8 && abs(cos(theta[1])) < 0.8 && Ep[0]>0.8 && Ep[1]>0.8 && abs(z[0])<10 && abs(z[1])<10 && vxy[0]<1.0 && vxy[1]<1.0", //bb_sel
    //MH selection
    //"ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2",
    "ptmin>0.2 && S>0.06 && maxctheta<0.8 && minctheta>-0.8 && Nchc>2 && maxNmuhit==0",
    //PID
    PID,
    //TauTau selection common cuts
    { 
        "Ncg==2 && Ncg==Ncc",
        "abs(vz[0])<10 && abs(vz[1])<10",
        "vxy[0]<1.0 && vxy[1]<1.0",
        "E[0]>0.025 && E[1]>0.025",
        "q[0]==-1 && q[1]==1",
        "p[0] < 1.1 && p[1] < 1.1",
        "pt[0] > 0.2 && pt[1] > 0.2",
        "cgood",
        "2.5 < tof[0] && tof[0] < 5.5 && 2.5 < tof[1] && tof[1] < 5.5",
        "ptem50>0.25 && ptem50<1.1",
    }
    , 
    //channel selections
      { 
        {"eX",      "NnE50==0 && eX && missed_photon_angle"},
        {"eρ",     "Nng==2   && eX && (Mpi0[0] < 0.14 && Mpi0[0]>0.12) && good_emc_time" },
      }
    );

void analysis(void) {
};
