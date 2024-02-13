#ifndef FluggFlux_h
#define FluggFlux_h

#pragma once

#include <iostream>
#include <iomanip>
#include <string>

#include "TChain.h"
#include "TH1.h"
#include "TFile.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "FluggTree.h"
#include "IFlux.h"
#include "Outputs.h"
#include "Constants.h"

class FluggFlux: virtual private NuMI::IFlux
{
public :
  FluggFlux(std::string pattern="", std::string outfile="FluggFlux.root");
  ~FluggFlux() override {};

  void CalculateFlux() override;
  int CalculateWeight(FluggTree* decay,
                      const TVector3& xyz, double& enu, double& wgt_xy);

private:
  TChain* cflux;
  double AccumulatedPOT=0.;

  FluggTree* fluxNtuple;
  RootOutput* fOutput;
};
#endif
