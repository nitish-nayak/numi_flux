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
  outTree->Branch("nuL", &nuL, "nuL/F");
  outTree->Branch("wgt", &wgt, "wgt/F");
  outTree->Branch("wgt_ppfx", &wgt_ppfx, "wgt_ppfx/F");
  outTree->Branch("wgt_tgtatt", &wgt_tgtatt, "wgt_tgtatt/F");
  outTree->Branch("wgt_absorp", &wgt_absorp, "wgt_absorp/F");
  outTree->Branch("wgt_ttpcpion", &wgt_ttpcpion, "wgt_ttpcpion/F");
  outTree->Branch("wgt_ttpckaon", &wgt_ttpckaon, "wgt_ttpckaon/F");
  outTree->Branch("wgt_ttpcnucleon", &wgt_ttpcnucleon, "wgt_ttpcnucleon/F");
  outTree->Branch("wgt_ttncpion", &wgt_ttncpion, "wgt_ttncpion/F");
  outTree->Branch("wgt_ttnucleona", &wgt_ttnucleona, "wgt_ttnucleona/F");
  outTree->Branch("wgt_ttmesoninc", &wgt_ttmesoninc, "wgt_ttmesoninc/F");
  outTree->Branch("wgt_others", &wgt_others, "wgt_others/F");
  outTree->Branch("wgt_ppfxunivs", &wgt_ppfxunivs);

  outTree->Branch("ptype", &ptype, "ptype/I");
  outTree->Branch("gptype", &gptype, "gptype/I");
  outTree->Branch("ntype", &ntype, "ntype/I");
  outTree->Branch("ncascade", &ncascade, "ncascade/I");
  outTree->Branch("pmedium", &pmedium, "pmedium/I");
  outTree->Branch("decaytype", &decaytype, "decaytype/I");
  outTree->Branch("pE", &pE, "pE/F");
  outTree->Branch("pPt", &pPt, "pPt/F");
  outTree->Branch("pPz", &pPz, "pPz/F");
  outTree->Branch("pTheta", &pTheta, "pTheta/F");
  outTree->Branch("pScTheta", &pScTheta, "pScTheta/F");
  outTree->Branch("pxF", &pxF, "pxF/F");
  outTree->Branch("pxF_inc", &pxF_inc, "pxF_inc/F");
  outTree->Branch("pvx", &pvx, "pvx/F");
  outTree->Branch("pvy", &pvy, "pvy/F");
  outTree->Branch("pvz", &pvz, "pvz/F");
  outTree->Branch("gpE", &gpE, "gpE/F");
  outTree->Branch("gpPt", &gpPt, "gpPt/F");
  outTree->Branch("gpPz", &gpPz, "gpPz/F");
  outTree->Branch("gpTheta", &gpTheta, "gpTheta/F");
  outTree->Branch("gpxF", &gpxF, "gpxF/F");
  outTree->Branch("pProc", &pProc);
  outTree->Branch("nvx", &nvx, "nvx/F");
  outTree->Branch("nvy", &nvy, "nvy/F");
  outTree->Branch("nvz", &nvz, "nvz/F");

  outTree->Branch("E_novand", &E_novand, "E_novand/F");
  outTree->Branch("E_minosnd", &E_minosnd, "E_minosnd/F");
  outTree->Branch("E_minerva", &E_minerva, "E_minerva/F");
  outTree->Branch("wgt_novand", &wgt_novand, "wgt_novand/F");
  outTree->Branch("wgt_minosnd", &wgt_minosnd, "wgt_minosnd/F");
  outTree->Branch("wgt_minerva", &wgt_minerva, "wgt_minerva/F");
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

