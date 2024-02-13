#ifndef Outputs_h
#define Outputs_h

#pragma once

#include <string>

#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"

struct RootOutput
{
  RootOutput(std::string outname="flux.root");
  ~RootOutput(){
    fout->Close();
  };

  void Write();

  Float_t nuE       = -5.;
  Float_t nudirX    = -5.;
  Float_t nudirY    = -5.;
  Float_t nudirZ    = -5.;

  Float_t wgt       = 1.;
  Float_t wgt_ppfx  = 1.;
  Int_t ptype       = -1000;
  Int_t ncascade    = -1;
  Int_t ntype       = -1000;
  Int_t pmedium     = -1000;
  Int_t decaytype   = -1000;

  Float_t pE        = -5.;
  Float_t pPt       = -5.;
  Float_t pPz       = -5.;
  Float_t pTheta    = -5.;
  Float_t pxF       = -5.;

  TFile* fout;

  TH1D* numuFluxHisto;
  TH1D* anumuFluxHisto;
  TH1D* nueFluxHisto;
  TH1D* anueFluxHisto;
  TH1D* hPOT;
  TTree* outTree;

  double histMin    = 0.;
  double histMax    = 20.;
  int histNbins     = 4000;
};
#endif
