#ifndef Dk2NuFlux_h
#define Dk2NuFlux_h

#pragma once

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#include "TChain.h"
#include "TH1.h"
#include "TH2D.h"
#include "TFile.h"
#include "TSystem.h"
#include "TVector3.h"

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/calcLocationWeights.h"

// Krishan: Updated this script to take in the dk2u file format instead of flugg

class Dk2NuFlux {
public :

  int Nfiles = 0;

  static const int numu  =  14;
  static const int anumu = -14;
  static const int nue   =  12;
  static const int anue  = -12;

  int highest_evtno = 0;
  double NominalPOT = 6e20;
  bool debug = false;
  double fDefaultWeightCorrection = 1./(10000. * TMath::Pi());
  double Ntarget = 4.76e31/56.41e6* 256.35*233*1036.8; //TPC active!!!
  double AccumulatedPOT=0.;
  int treeNumber = -1;

  double histMin = 0;
  double histMax = 20;
  int histNbins = 4000;

  TChain *cflux;
  TChain *cflux_meta;

  TH1D* numuFluxHisto;
  TH1D* anumuFluxHisto;
  TH1D* nueFluxHisto;
  TH1D* anueFluxHisto;
  TH1D* numuCCHisto;
  TH1D* anumuCCHisto;
  TH1D* nueCCHisto;
  TH1D* anueCCHisto;
  TH1D* hPOT;
  TTree* outTree;
  //
  // // MIPP
  // TH2D* pionplus_MIPP;
  // TH2D* pionminus_MIPP;
  // TH2D* Kplus_MIPP;
  // TH2D* Kminus_MIPP;
  // // NA49
  // TH2D* pionplus_NA49;
  // TH2D* pionminus_NA49;
  // TH2D* Kplus_NA49;
  // TH2D* Kminus_NA49;

  // TGraph *genieXsecNumuCC;
  // TGraph *genieXsecNumubarCC;
  // TGraph *genieXsecNueCC;
  // TGraph *genieXsecNuebarCC;


  Dk2NuFlux(string pattern="/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_dk2nu_zero_threshold/FHC/g4numiv6_minervame_me000z200i_0_0000.root", string outfile="NuMIFlux.root");

  virtual ~Dk2NuFlux();

  void CalculateFlux();
  TVector3 RandomInTPC();
  TVector3 FromDetToBeam(const TVector3& det);
  double estimate_pots(int highest_potnum);
  int calcEnuWgt(bsim::Dk2Nu* decay, const TVector3& xyz, double& enu, double& wgt_xy);

private:

  Float_t nuE;
  Float_t wgt;
  Float_t wgt_ppfx;
  Int_t ptype;
  Int_t ncascade;
  Int_t ntype;
  Int_t pmedium;
  Int_t decaytype;

  TFile* fout;
};

#endif
