#define Outputs_cxx

#include "Outputs.h"


//___________________________________________________________________________
RootOutput::RootOutput(std::string outname)
{
  fout = new TFile(outname.c_str(), "RECREATE");

  //Inizialise histos
  TString titleBase1 = "Neutrino Flux;";
  TString titleBase2 = " Energy [GeV];";
  TString titleBase3 = " / cm^{2} / 6e20 POT";
  // numu
  numuFluxHisto = new TH1D("numuFluxHisto", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
  numuCCHisto = new TH1D("numuCCHisto", "numu CC; #nu_{#mu} Energy [GeV]; #nu_{#mu} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  // anumu
  anumuFluxHisto = new TH1D("anumuFluxHisto", (titleBase1 + "#bar{#nu}_{#mu}" + titleBase2 +"#bar{#nu}_{#mu}" + titleBase3),histNbins,histMin,histMax);
  anumuCCHisto = new TH1D("anumuCCHisto", "numu bar CC; #bar{#nu}_{#mu} Energy [GeV]; #bar{#nu}_{#mu} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  // nue
  nueFluxHisto = new TH1D("nueFluxHisto", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
  nueCCHisto = new TH1D("nueCCHisto", "nue CC; #nu_{e} Energy [GeV]; #nu_{e} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  // anue
  anueFluxHisto = new TH1D("anueFluxHisto", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 + "#bar{#nu}_{e}" + titleBase3),histNbins,histMin,histMax);
  anueCCHisto = new TH1D("anueCCHisto", "nue bar CC; #bar{#nu}_{e} Energy [GeV]; #bar{#nu}_{#mu} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  hPOT = new TH1D("POT", "Total POT", 1, 0, 1);

  outTree = new TTree("outTree", "outTree");
  outTree->Branch("nuE", &nuE, "nuE/F");
  outTree->Branch("wgt", &wgt, "wgt/F");
  outTree->Branch("wgt_ppfx", &wgt_ppfx, "wgt_ppfx/F");
  outTree->Branch("ptype", &ptype, "ptype/I");
  outTree->Branch("ntype", &ntype, "ntype/I");
  outTree->Branch("ncascade", &ncascade, "ncascade/I");
  outTree->Branch("pmedium", &pmedium, "pmedium/I");
  outTree->Branch("decaytype", &decaytype, "decaytype/I");

}

//___________________________________________________________________________
RootOutput::Write()
{
  fout->cd();
  numuFluxHisto  -> Write();
  anumuFluxHisto -> Write();
  nueFluxHisto   -> Write();
  anueFluxHisto  -> Write();

  hPOT->Write();
  outTree->Write();

}

