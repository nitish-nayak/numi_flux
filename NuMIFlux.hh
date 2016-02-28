#ifndef NuMIFlux_h
#define NuMIFlux_h

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

  int numu = 56;
  int highest_evtno = 0;
  double NominalPOT = 6e20;
  bool debug = false;
  double fDefaultWeightCorrection = 1./(10000. * TMath::Pi());
  double Ntarget = 4.76e31/56.41e6*256.35*233*1036.8; //TPC active!!!

  double histMin = 0;
  double histMax = 6;
  int histNbins = 120;

  TChain *cflux;

  TH1D* nuFluxHisto;
  TH1D* nuCCHisto;
  TFile* f = new TFile("NuMIFlux.root", "RECREATE");


  NuMIFlux(string pattern="/uboone/data/flux/numi/current/flugg_mn000z200i_20101117.gpcfgrid_lowth/flugg_mn000z200i_20101117.gpcfgrid_lowth_00*.root");
  virtual ~NuMIFlux();

  void CalculateFlux();
  TVector3 RandomInTPC();
  TVector3 FromDetToBeam(const TVector3& det);
  double estimate_pots(int highest_potnum);
  int calcEnuWgt( FluxNtuple& decay, const TVector3& xyz, double& enu, double& wgt_xy);

};

#endif

#ifdef NuMIFlux_cxx

NuMIFlux::NuMIFlux(string pattern) {

  const char* path = "/uboone/app/users/mdeltutt/NuMIFlux";
  if ( path ) {
    TString libs = gSystem->GetDynamicPath();
    libs += ":";
    libs += path;
    //libs += "/lib";
    gSystem->SetDynamicPath(libs.Data());       
    gSystem->Load("FluxNtuple_C.so");
  }

  cflux = new TChain("h10");
  cflux->Add(pattern.c_str());

  Nfiles = cflux->GetNtrees();
  cout << "Number of files: " << Nfiles << endl;

  //Inizialise histos
  nuFluxHisto = new TH1D("nuFluxHisto", "Neutrino Flux; #nu_{#mu} Energy [GeV];#nu_{#mu} / cm^{2} / 6e20 POT",histNbins,histMin,histMax);
  nuCCHisto = new TH1D("nuCCHisto", "numu CC; #nu_{#mu} Energy [GeV]; #nu_{#mu} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
}

NuMIFlux::~NuMIFlux() {

}

#endif // #ifdef NuMIFlux_cxx
