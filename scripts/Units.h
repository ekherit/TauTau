/*
 * =====================================================================================
 *
 *       Filename:  Units.h
 *
 *    Description:  Define some unints and physical constants
 *
 *        Version:  1.0
 *        Created:  30.12.2020 10:47:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
 *   Organization:  Budker Insitute of Nuclear Physics
 *
 * =====================================================================================
 */

#ifndef IBN_TAU_UNITS_H
#define IBN_TAU_UNITS_H

#include <TMath.h>

constexpr double GeV=1.0;
constexpr double MeV=1e-3*GeV;
constexpr double keV=1e-6*GeV;
constexpr double barn = 1.0;
constexpr double pb = 1e-12*barn;
constexpr double nb = 1e-9*barn;

constexpr double MTAU=1776.86*MeV;
constexpr double MPION=0.13957061*GeV; 
constexpr double MPI0=0.134977*GeV;

constexpr const double MJPSI = 3096.900*MeV;
constexpr const double MJPSI_ERROR = 6.0*keV;
constexpr const double MPSIP = 3686.097*MeV; 
constexpr const double MPSIP_ERROR =10.0*keV;


constexpr double ME_PDG2011=0.510998910*MeV; //+-0.000000013 
constexpr double ME_PDG2019=0.5109989461*MeV; //+-0.000000013
constexpr double ME_PDG=ME_PDG2019;// mass of electron
constexpr double ME=ME_PDG;// mass of electron

constexpr double SIGMA_TOMSON = 0.665245854*barn; 
constexpr double SIGMA_CONST = SIGMA_TOMSON/pb* ME*ME/2.;  // GeV^2/pb


constexpr double ALPHA_PDG2019=1/137.035999139; 
constexpr double ALPHA = ALPHA_PDG2019;
constexpr double ALPHAPI=ALPHA/TMath::Pi(); //alpha / pi
constexpr double PIALPHA=ALPHA*TMath::Pi(); 

#endif
