#include "../include/Analyzer1D.h"
#include "glob.h"

using namespace ROOT;
using namespace ROOT::RDF;
using namespace NuMI;

void analyzer_example(std::string pattern, unsigned int nthreads = 1)
{
  glob_t glob_result;
  int result = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
  if(result != 0){
    std::cerr << "Couldn't find files! Exiting!" << std::endl;
    exit(1);
  }

  nthreads > 1 ? ROOT::EnableImplicitMT(nthreads) : void();

  TChain chain("outTree");
  chain.Add(pattern.c_str());
  RDataFrame df(chain);

  TH1D* dummy_flux = new TH1D("flux", "#nu_{#mu} Flux", 100, 0, 5);

  // define my histogram helper class
  Ana1DHelper helper_cv = Ana1DHelper(*dummy_flux, nthreads);
  Ana1DHelper helper = Ana1DHelper(*dummy_flux, nthreads);

  // fancy indexing is possible but I haven't implemented handling for it yet
  const int Nunivs = 10;
  // auto univWgts = [&Nunivs](float wgt, const RVecF &wgt_ppfxunivs){ return wgt*VecOps::Take(wgt_ppfxunivs, Nunivs); }
  std::vector<Ana1DHelper> helper_univs;
  for(int iUniv = 0; iUniv < Nunivs; iUniv++)
      helper_univs.emplace_back(*dummy_flux, nthreads);

  // define my selections and apply it to my data
  auto selNumu = [](int ntype){ return ntype == 14; };
  std::cout << "Filtering.." << std::endl;
  auto f_df = df.Filter(selNumu, {"ntype"});

  // define my weights
  auto cvWgt = [](float wgt, float wgt_ppfx){ return wgt*wgt_ppfx; };
  f_df = f_df.Define("flux_cv_wgt", cvWgt, {"wgt", "wgt_ppfx"});

  // the move is important because the helper class should not survive the event loop
  // book the action
  std::cout << "Booking univs.." << std::endl;
  std::vector<RResultPtr<TH1D>> r_numu_univs;
  for(int iUniv = 0; iUniv < Nunivs; iUniv++){
    auto univWgt = [iUniv](float wgt, const RVecF &wgt_ppfxunivs){ return wgt*wgt_ppfxunivs[iUniv]; };
    f_df = f_df.Define("flux_univ"+std::to_string(iUniv)+"_wgt",
                       univWgt, {"wgt", "wgt_ppfxunivs"});

    auto r_numu_univ = f_df.Book<float, float>(std::move(helper_univs[iUniv]),
                                               {"nuE", "flux_univ"+std::to_string(iUniv)+"_wgt"});
    r_numu_univs.emplace_back(r_numu_univ);
  }
  RResultPtr<TH1D> r_numu_cv = f_df.Book<float, float>(std::move(helper_cv), {"nuE", "flux_cv_wgt"});
  RResultPtr<TH1D> r_numu = f_df.Book<float, float>(std::move(helper), {"nuE", "wgt"});

  // register some callbacks
  auto c = f_df.Count();
  c.OnPartialResult(100000, [](unsigned int n){ std::cout << "===== " << n << " =====" << std::endl; });

  // booking is all lazy, the event loop on the dataframe is triggered by trying to access the pointer (for everything)
  // trigger it here when I dereference c
  std::cout << "Total Events : " << *c << std::endl;
  std::cout << "CV Flux : " << r_numu_cv->Integral() << std::endl;
  std::cout << "Raw Flux : " << r_numu->Integral() << std::endl;

  // now do pot
  TH1D* hPOT = new TH1D("hPOT", "", 1, 0, 1);
  for(int i = 0; i < (int)glob_result.gl_pathc; i++){
    TFile* fin = TFile::Open(glob_result.gl_pathv[i], "read");
    TH1D* h = (TH1D*)fin->Get("POT");
    hPOT->Add(h);
    fin->Close();
  }
  std::cout << "Total POT : " << hPOT->Integral() << std::endl;

  // the histograms now exist so we can save them
  TFile* fout = new TFile("numu_fluxes.root", "recreate");
  hPOT->Write("POT");
  r_numu_cv->Write("numu_cv");
  r_numu->Write("numu");
  fout->mkdir("univs"); fout->cd("univs");
  for(int i = 0; i < Nunivs; i++)
    r_numu_univs[i]->Write(TString::Format("numu_univ%d", i));

  fout->Close();
}
