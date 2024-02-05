#ifndef NuMIFlux_h
#define NuMIFlux_h

#pragma once

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#include "TChain.h"
#include "TH1.h"
#include "TFile.h"
#include "TSystem.h"
#include "TVector3.h"

#include "FluggNtuple/FluxNtuple.h"
//#include "calcLocationWeights.h"

class NuMIFlux {
public :

  int Nfiles = 0;

  static const int numu  = 56;
  static const int anumu = 55;
  static const int nue   = 53;
  static const int anue  = 52;

  int highest_evtno = 0;
  double NominalPOT = 6e20;
  bool debug = false;
  double fDefaultWeightCorrection = 1./(10000. * TMath::Pi());
  double Ntarget = 4.76e31/56.41e6*256.35*233*1036.8; //TPC active!!!
  double AccumulatedPOT=0.;
  int treeNumber = -1;

  double histMin = 0;
  double histMax = 20;
  int histNbins = 400;

  TChain *cflux;

  FluxNtuple *fluxNtuple;

  TH1D* numuFluxHisto;
  TH1D* anumuFluxHisto;
  TH1D* nueFluxHisto;
  TH1D* anueFluxHisto;
  TH1D* numuCCHisto;
  TH1D* hPOT;
  TTree* outTree;
  // TGraph *genieXsecNumuCC;

  NuMIFlux(string pattern="/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files/rhc/*_7000.root");
  virtual ~NuMIFlux();

  void CalculateFlux(string outfile="NuMIFlux.root");
  TVector3 RandomInTPC();
  TVector3 FromDetToBeam(const TVector3& det);
  double estimate_pots(int highest_potnum);
  int calcEnuWgt( FluxNtuple* decay, const TVector3& xyz, double& enu, double& wgt_xy);

private:

  Float_t nuE;
  Float_t wgt;
  Int_t ptype;
  Int_t ncascade;
  Int_t ntype;
  Int_t pmedium;
  Int_t decaytype;

};

#endif
