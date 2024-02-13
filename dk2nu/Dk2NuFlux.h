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
#include "IFlux.h"
#include "Outputs.h"
#include "Constants.h"

// Krishan: Updated this script to take in the dk2u file format instead of flugg

class Dk2NuFlux: virtual private NuMI::IFlux
{
public :

  Dk2NuFlux(std::string pattern="", std::string outfile="Dk2NuFlux.root");
  ~Dk2NuFlux() override {};

  void CalculateFlux() override;
  int CalculateWeight(bsim::Dk2Nu* decay,
                      const TVector3& xyz, double& enu, double& wgt_xy);

private:
  TChain* cflux;
  TChain *cflux_meta;
  double AccumulatedPOT=0.;

  RootOutput* fOutput;
};

#endif
