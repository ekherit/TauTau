// =====================================================================================
//
//       Filename:  SelectionConfig.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  20.10.2015 12:26:35
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once
#include <iostream>

struct SelectionConfig
{
  double CENTER_MASS_ENERGY;      //center mass energy

  int MIN_CHARGED_TRACKS; //minimum good charged tracks in selection
  int MAX_CHARGED_TRACKS; //maximum good charged tracks in selection

  int MAX_NEUTRAL_TRACKS; //maximum good neutral tracks in selection

  double IP_MAX_RHO, IP_MAX_Z; //interection point cut
  int  USE_VERTEX_DB; //user vertex db information

  double MAX_COS_THETA_FOR_CHARGED; //maximum  cos(theta) for good charged track
  double MIN_EMC_ENERGY_FOR_CHARGED;


  double MIN_EMC_ENERGY_FOR_NEUTRAL;

  double MIN_MOMENTUM;
  double MAX_MOMENTUM;

  double MIN_TRANSVERSE_MOMENTUM;
  double MAX_TRANSVERSE_MOMENTUM;

  double MIN_EP_RATIO;
  double MAX_EP_RATIO;

  double MIN_PTEM;
  double MAX_PTEM;

  double MIN_TOF;
  double MAX_TOF;

  double DELTA_MJPSI;
  //double EMC_ENDCUP_MIN_COS_THETA;
  //double EMC_ENDCUP_MAX_COS_THETA;
  //double EMC_ENDCUP_MIN_ENERGY;
  //double EMC_BARREL_MAX_COS_THETA;
  //double EMC_BARREL_MIN_ENERGY;

  //double NEUTRAL_CLOSE_CHARGED_ANGLE;

  //double MAX_MUON_EP_RATIO;
  //double MAX_KAON_EP_RATIO;

  //double MIN_MUON_EP_RATIO;
  //double MIN_KAON_EP_RATIO;

  //double MAX_PION_MOMENTUM; //maximum pion momentum
  //double MIN_PION_MOMENTUM; //maximum pion momentum

  //double MIN_RECOIL_MASS; //minimum recoil mass cut
  //double MAX_RECOIL_MASS; //minimum recoil mass cut

  //double MIN_KAON_MOMENTUM; //minimum kaon momentum
  //double MAX_KAON_MOMENTUM; //maximum pion momentum

  //double MIN_MUON_MOMENTUM; //minimum kaon momentum
  //double MAX_MUON_MOMENTUM; //maximum pion momentum

  //double MIN_INVARIANT_MASS; //minimum invariant  mass cut
  //double MAX_INVARIANT_MASS; //manimum invariant  mass cut

  //double MIN_KAON_MISSING_MASS;   //minimum kaon missing mass
  //double MAX_KAON_MISSING_MASS;   //minimum kaon missing mass
  //double MIN_MUON_MISSING_MASS;   //minimum muon missing mass
  //double MAX_MUON_MISSING_MASS;   //minimum muon missing mass

  //double MIN_MISSING_MASS; 
  //double MAX_MISSING_MASS; 


  //double GAMMA_GAMMA_MIN_INV_MASS;
  //bool FILL_MDC;
  //bool FILL_EMC;
  //bool FILL_DEDX;
  //bool FILL_TOF;
  //bool FILL_MUC;
  //bool FILL_NEUTRAL;


  int TEST_COMBINATIONS;

  void print(void)
  {
    //std::cout  << "CENTER_MASS_ENERGY = " << CENTER_MASS_ENERGY << std::endl;      //center mass energy
    //std::cout  << "MIN_CHARGED_TRACKS = " << MIN_CHARGED_TRACKS << std::endl; //minimum good charged tracks in selection
    //std::cout  << "MAX_CHARGED_TRACKS = " << MAX_CHARGED_TRACKS << std::endl; //maximum good charged tracks in selection
    //std::cout  << "MAX_NEUTRAL_TRACKS = " << MAX_NEUTRAL_TRACKS << std::endl; //maximum good neutral tracks in selection
    //std::cout  << "IP_MAX_RHO = " << IP_MAX_RHO << std::endl; //std::cout << erection point cut
    //std::cout  << "MAX_COS_THETA = " << MAX_COS_THETA << std::endl; //maximum  cos(theta) for good charged track

    //std::cout  << "EMC_ENDCUP_MIN_COS_THETA = " << EMC_ENDCUP_MIN_COS_THETA << std::endl;
    //std::cout  << "EMC_ENDCUP_MAX_COS_THETA = " << EMC_ENDCUP_MAX_COS_THETA << std::endl;
    //std::cout  << "EMC_ENDCUP_MIN_ENERGY = " << EMC_ENDCUP_MIN_ENERGY << std::endl;
    //std::cout  << "EMC_BARREL_MAX_COS_THETA = " << EMC_BARREL_MAX_COS_THETA << std::endl;
    //std::cout  << "EMC_BARREL_MIN_ENERGY = " << EMC_BARREL_MIN_ENERGY << std::endl;

    //std::cout  << "NEUTRAL_CLOSE_CHARGED_ANGLE = " << NEUTRAL_CLOSE_CHARGED_ANGLE << std::endl;

    //std::cout  << "MAX_MUON_EP_RATIO = " << MAX_MUON_EP_RATIO << std::endl;
    //std::cout  << "MAX_KAON_EP_RATIO = " << MAX_KAON_EP_RATIO << std::endl;

    //std::cout  << "MIN_MUON_EP_RATIO = " << MIN_MUON_EP_RATIO << std::endl;
    //std::cout  << "MIN_KAON_EP_RATIO = " << MIN_KAON_EP_RATIO << std::endl;

    //std::cout  << "MAX_PION_MOMENTUM = " << MAX_PION_MOMENTUM << std::endl; //maximum pion momentum
    //std::cout  << "MIN_PION_MOMENTUM = " << MIN_PION_MOMENTUM << std::endl; //maximum pion momentum

    //std::cout  << "MIN_RECOIL_MASS = " << MIN_RECOIL_MASS << std::endl; //minimum recoil mass cut
    //std::cout  << "MAX_RECOIL_MASS = " << MAX_RECOIL_MASS << std::endl; //minimum recoil mass cut

    //std::cout  << "MIN_KAON_MOMENTUM = " << MIN_KAON_MOMENTUM << std::endl; //minimum kaon momentum
    //std::cout  << "MAX_KAON_MOMENTUM = " << MAX_KAON_MOMENTUM << std::endl; //maximum pion momentum

    //std::cout  << "MIN_MUON_MOMENTUM = " << MIN_MUON_MOMENTUM << std::endl; //minimum kaon momentum
    //std::cout  << "MAX_MUON_MOMENTUM = " << MAX_MUON_MOMENTUM << std::endl; //maximum pion momentum

    //std::cout  << "MIN_INVARIANT_MASS = " << MIN_INVARIANT_MASS << std::endl; //minimum invariant  mass cut
    //std::cout  << "MAX_INVARIANT_MASS = " << MAX_INVARIANT_MASS << std::endl; //manimum invariant  mass cut

    //std::cout  << "MIN_KAON_MISSING_MASS = " << MIN_KAON_MISSING_MASS << std::endl;   //minimum kaon missing mass
    //std::cout  << "MAX_KAON_MISSING_MASS = " << MAX_KAON_MISSING_MASS << std::endl;   //minimum kaon missing mass
    //std::cout  << "MIN_MUON_MISSING_MASS = " << MIN_MUON_MISSING_MASS << std::endl;   //minimum muon missing mass
    //std::cout  << "MAX_MUON_MISSING_MASS = " << MAX_MUON_MISSING_MASS << std::endl;   //minimum muon missing mass

    //std::cout  << "MIN_MISSING_MASS = " << MIN_MISSING_MASS << std::endl; 
    //std::cout  << "MAX_MISSING_MASS = " << MAX_MISSING_MASS << std::endl; 

    //std::cout  << "MAX_KIN_CHI2 = " << MAX_KIN_CHI2 << std::endl; //maximum chi2 for kinematic fit
    //std::cout  << "MAX_PID_CHI2 = " << MAX_PID_CHI2 << std::endl; //maximum chi2 for my PID

    //std::cout  << "FILL_MDC = " << FILL_MDC << std::endl;
    //std::cout  << "FILL_EMC = " << FILL_EMC << std::endl;
    //std::cout  << "FILL_DEDX = " << FILL_DEDX << std::endl;
    //std::cout  << "FILL_TOF = " << FILL_TOF << std::endl;
    //std::cout  << "FILL_MUC = " << FILL_MUC << std::endl;
    //std::cout  << "FILL_NEUTRAL = " << FILL_NEUTRAL << std::endl;

  }

  void print_relevant(void)
  {
    std::cout << "########################################################################################################\n";
    std::cout << "########################################################################################################\n";
    std::cout << "                 SELECTION CONFIGURATION                \n";
    std::cout << "########################################################################################################\n";
    std::cout << "\n";
    //std::cout << "TEST = " << TEST << std::endl;      //center mass energy
    std::cout << "CENTER_MASS_ENERGY = "         << CENTER_MASS_ENERGY         << "\n"; 

    std::cout << "MIN_CHARGED_TRACKS = "         << MIN_CHARGED_TRACKS         << "\n";
    std::cout << "MAX_CHARGED_TRACKS = "         << MAX_CHARGED_TRACKS         << "\n";

    std::cout << "IP_MAX_Z = "                   << IP_MAX_Z                   << "\n"; 
    std::cout << "IP_MAX_RHO = "                 << IP_MAX_RHO                 << "\n"; 
    std::cout << "USE_VERTEX_DB = "              << USE_VERTEX_DB              << "\n"; 
    std::cout << "MAX_COS_THETA = "              << MAX_COS_THETA_FOR_CHARGED  << "\n"; 
    std::cout << "MIN_EMC_ENERGY_FOR_CHARGED = " << MIN_EMC_ENERGY_FOR_CHARGED << "\n"; 

    std::cout << "MIN_EMC_ENERGY_FOR_NEUTRAL = " << MIN_EMC_ENERGY_FOR_NEUTRAL << "\n"; 

    std::cout << "MIN_MOMENTUM = "               << MIN_MOMENTUM               << "\n"; 
    std::cout << "MAX_MOMENTUM = "               << MAX_MOMENTUM               << "\n"; 
    std::cout << "MIN_TRANSVERSE_MOMENTUM = "    << MIN_TRANSVERSE_MOMENTUM    << "\n"; 
    std::cout << "MAX_TRANSVERSE_MOMENTUM = "    << MAX_TRANSVERSE_MOMENTUM    << "\n"; 
    std::cout << "MIN_EP_RATIO = "               << MIN_EP_RATIO               << "\n"; 
    std::cout << "MAX_EP_RATIO = "               << MAX_EP_RATIO               << "\n"; 
    std::cout << "MIN_PTEM = "                   << MIN_PTEM                   << "\n"; 
    std::cout << "MAX_PTEM = "                   << MAX_PTEM                   << "\n"; 
    std::cout << "MIN_TOF = "                    << MIN_TOF                    << "\n"; 
    std::cout << "MAX_TOF = "                    << MAX_TOF                    << "\n"; 
    std::cout << "########################################################################################################\n";
  };

};
