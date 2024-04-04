#ifndef Analyzer1D_h
#define Analyzer1D_h

#pragma once

#include <iostream>
#include <iomanip>
#include <string>

#include <ROOT/RDataFrame.hxx>
#include <memory>
#include "IAnalyzer.h"

namespace NuMI {

  static int nAna1D;
  class Ana1DHelper : public IAnalyzer<TH1D, 1>, public ROOT::Detail::RDF::RActionImpl<Ana1DHelper>
  {
  public:
    // we need these
    using Result_t = TH1D;
  protected:
    std::vector<std::shared_ptr<TH1D>> fHistos; // one per data processing slot per NDIM
    std::shared_ptr<TH1D> fHistoResult; // one per NDIM
    unsigned int fThreads;
    unsigned int fDim;

  public:
    Ana1DHelper(const TH1D& dummy, unsigned int nThreads=1)
      : IAnalyzer<TH1D, 1>(dummy, nThreads), fThreads(nThreads)
    {
      nAna1D++;
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
      fHistoResult = std::make_shared<TH1D>(TString::Format("%s_id%d_result", dummy.GetName(), nAna1D), dummy.GetTitle(),
                                            dummy.GetNbinsX(), xarr.GetArray());

      for(unsigned int i = 0; i < nSlots; i++) {
        fHistos.emplace_back(std::make_shared<TH1D>(
                             TString::Format("%s_id%d_slot%d", dummy.GetName(), nAna1D, i), dummy.GetTitle(),
                             dummy.GetNbinsX(), xarr.GetArray())
        );
      }
    }
    // don't allow copies
    Ana1DHelper(Ana1DHelper &&) = default;
    Ana1DHelper(const Ana1DHelper &) = delete;

    std::shared_ptr<TH1D> GetResultPtr() const { return fHistoResult; }

    // we need these for df.Book()
    void Initialize() {}
    void InitTask(TTreeReader *, unsigned int) {}

    template <typename... DataTypes>
    void Exec(unsigned int slot, DataTypes... values)
    {
      // first one is val, second one is weight
      // cast everything into double
      std::array<double, sizeof...(DataTypes)> valuesArr{static_cast<double>(values)...};

      if(valuesArr.size() == 1){
        fHistos[slot]->Fill(valuesArr[0]);
      }
      else if(valuesArr.size() == 2){
        fHistos[slot]->Fill(valuesArr[0], valuesArr[1]);
      }
      else{
       std::cerr << "Expect max two columns (var, wgt) as input!! Exiting!" << std::endl;
       std::exit(1);
      }

    }

    void Finalize() const
    {
      for(unsigned int i = 0; i < fHistos.size(); i++)
        fHistoResult->Add(fHistos[i].get());
    }

    std::string GetActionName(){
      return "Ana1DHelper";
    }
  };
  // end helper
}

#endif
