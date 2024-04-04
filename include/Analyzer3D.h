#ifndef Analyzer3D_h
#define Analyzer3D_h

#pragma once

#include <iostream>
#include <iomanip>
#include <string>

#include <ROOT/RDataFrame.hxx>
#include <memory>
#include "IAnalyzer.h"

namespace NuMI {

  static int nAna3D;
  class Ana3DHelper : public IAnalyzer<TH3D, 1>, public ROOT::Detail::RDF::RActionImpl<Ana3DHelper>
  {
  public:
    // we need these
    using Result_t = TH3D;
  protected:
    std::vector<std::shared_ptr<TH3D>> fHistos; // one per data processing slot per NDIM
    std::shared_ptr<TH3D> fHistoResult; // one per NDIM
    unsigned int fThreads;
    unsigned int fDim;

  public:
    Ana3DHelper(const TH3D& dummy, unsigned int nThreads=1)
      : IAnalyzer<TH3D, 1>(dummy, nThreads), fThreads(nThreads)
    {
      nAna3D++;
      unsigned int nSlots = ROOT::IsImplicitMTEnabled() ? fThreads : 1;

      // I hate that ROOT is forcing me to do this
      // do x
      int nbinsx = dummy.GetNbinsX();
      TArrayD xarr(nbinsx+1);
      if(!dummy.GetXaxis()->IsVariableBinSize()){
        for(int i = 1; i <= nbinsx+1; i++)
          xarr[i-1] = dummy.GetXaxis()->GetBinLowEdge(i);
      }
      else
        xarr = *(dummy.GetXaxis()->GetXbins());
      // now do y
      int nbinsy = dummy.GetNbinsY();
      TArrayD yarr(nbinsy+1);
      if(!dummy.GetYaxis()->IsVariableBinSize()){
        for(int i = 1; i <= nbinsy+1; i++)
          yarr[i-1] = dummy.GetYaxis()->GetBinLowEdge(i);
      }
      else
        yarr = *(dummy.GetYaxis()->GetXbins());
      // now do z
      int nbinsz = dummy.GetNbinsZ();
      TArrayD zarr(nbinsz+1);
      if(!dummy.GetZaxis()->IsVariableBinSize()){
        for(int i = 1; i <= nbinsz+1; i++)
          zarr[i-1] = dummy.GetZaxis()->GetBinLowEdge(i);
      }
      else
        zarr = *(dummy.GetZaxis()->GetXbins());

      // now initialize my histograms
      fHistoResult = std::make_shared<TH3D>(TString::Format("%s_id%d_result", dummy.GetName(), nAna3D), dummy.GetTitle(),
                                            dummy.GetNbinsX(), xarr.GetArray(), dummy.GetNbinsY(), yarr.GetArray(), dummy.GetNbinsZ(), zarr.GetArray());

      for(unsigned int i = 0; i < nSlots; i++) {
        fHistos.emplace_back(std::make_shared<TH3D>(
                             TString::Format("%s_id%d_slot%d", dummy.GetName(), nAna3D, i), dummy.GetTitle(),
                             dummy.GetNbinsX(), xarr.GetArray(), dummy.GetNbinsY(), yarr.GetArray(), dummy.GetNbinsZ(), zarr.GetArray())
        );
      }
    }
    // don't allow copies
    Ana3DHelper(Ana3DHelper &&) = default;
    Ana3DHelper(const Ana3DHelper &) = delete;

    std::shared_ptr<TH3D> GetResultPtr() const { return fHistoResult; }

    // we need these for df.Book()
    void Initialize() {}
    void InitTask(TTreeReader *, unsigned int) {}

    template <typename... DataTypes>
    void Exec(unsigned int slot, DataTypes... values)
    {
      // first one is val, second one is weight
      // cast everything into double
      std::array<double, sizeof...(DataTypes)> valuesArr{static_cast<double>(values)...};

      if(valuesArr.size() == 3){
        fHistos[slot]->Fill(valuesArr[0], valuesArr[1], valuesArr[2]);
      }
      else if(valuesArr.size() == 4){
        fHistos[slot]->Fill(valuesArr[0], valuesArr[1], valuesArr[2], valuesArr[3]);
      }
      else{
       std::cerr << "Expect atleast three columns (varx, vary, varz) as input!! Exiting!" << std::endl;
       std::exit(1);
      }

    }

    void Finalize() const
    {
      for(unsigned int i = 0; i < fHistos.size(); i++)
        fHistoResult->Add(fHistos[i].get());
    }

    std::string GetActionName(){
      return "Ana3DHelper";
    }
  };
  // end helper
}

#endif
