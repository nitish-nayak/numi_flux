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
  TH1D* nuCCHisto;
  TGraph *genieXsecNumuCC;
  TFile* f = new TFile("NuMIFlux.root", "RECREATE");


  NuMIFlux(string pattern="/uboone/data/flux/numi/v2/flugg_mn000z200i_rp11_lowth_pnut_f112c0f093bbird/flugg_mn000z200i_rp11_bs1.1_pnut_lowth_f112c0f093bbird_0*.root");
//"/uboone/data/flux/numi/current/flugg_mn000z200i_20101117.gpcfgrid_lowth/flugg_mn000z200i_20101117.gpcfgrid_lowth_00*.root"
//
  virtual ~NuMIFlux();

  void CalculateFlux();
  TVector3 RandomInTPC();
  TVector3 FromDetToBeam(const TVector3& det);
  double estimate_pots(int highest_potnum);
  int calcEnuWgt( FluxNtuple* decay, const TVector3& xyz, double& enu, double& wgt_xy);

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
  TString titleBase1 = "Neutrino Flux;";
  TString titleBase2 = " Energy [GeV];";
  TString titleBase3 = " / cm^{2} / 6e20 POT";
  // numu
  numuFluxHisto = new TH1D("numuFluxHisto", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
  // anumu
  anumuFluxHisto = new TH1D("anumuFluxHisto", (titleBase1 + "#bar{#nu}_{#mu}" + titleBase2 +"#bar{#nu}_{#mu}" + titleBase3),histNbins,histMin,histMax);
  // nue
  nueFluxHisto = new TH1D("nueFluxHisto", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
  // anue
  anueFluxHisto = new TH1D("anueFluxHisto", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 + "#bar{#nu}_{e}" + titleBase3),histNbins,histMin,histMax);
  nuCCHisto = new TH1D("nuCCHisto", "numu CC; #nu_{#mu} Energy [GeV]; #nu_{#mu} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
}

NuMIFlux::~NuMIFlux() {

}

#endif // #ifdef NuMIFlux_cxx
