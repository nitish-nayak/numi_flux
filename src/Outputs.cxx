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

  numuFluxHisto = new TH1D("numuFluxHisto", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
  anumuFluxHisto = new TH1D("anumuFluxHisto", (titleBase1 + "#bar{#nu}_{#mu}" + titleBase2 +"#bar{#nu}_{#mu}" + titleBase3),histNbins,histMin,histMax);
  nueFluxHisto = new TH1D("nueFluxHisto", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
  anueFluxHisto = new TH1D("anueFluxHisto", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 + "#bar{#nu}_{e}" + titleBase3),histNbins,histMin,histMax);

  hPOT = new TH1D("POT", "Total POT", 1, 0, 1);

  outTree = new TTree("outTree", "outTree");
  outTree->Branch("nuE", &nuE, "nuE/F");
  outTree->Branch("nudirX", &nudirX, "nudirX/F");
  outTree->Branch("nudirY", &nudirY, "nudirY/F");
  outTree->Branch("nudirZ", &nudirZ, "nudirZ/F");
  outTree->Branch("wgt", &wgt, "wgt/F");
  outTree->Branch("wgt_ppfx", &wgt_ppfx, "wgt_ppfx/F");
  outTree->Branch("ptype", &ptype, "ptype/I");
  outTree->Branch("ntype", &ntype, "ntype/I");
  outTree->Branch("ncascade", &ncascade, "ncascade/I");
  outTree->Branch("pmedium", &pmedium, "pmedium/I");
  outTree->Branch("decaytype", &decaytype, "decaytype/I");
  outTree->Branch("pE", &pE, "pE/F");
  outTree->Branch("pPt", &pPt, "pPt/F");
  outTree->Branch("pPz", &pPz, "pPz/F");
  outTree->Branch("pTheta", &pTheta, "pTheta/F");
  outTree->Branch("pxF", &pxF, "pxF/F");
}

//___________________________________________________________________________
void RootOutput::Write()
{
  fout->cd();
  numuFluxHisto  -> Write();
  anumuFluxHisto -> Write();
  nueFluxHisto   -> Write();
  anueFluxHisto  -> Write();

  hPOT->Write();
  outTree->Write();

}

