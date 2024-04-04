#ifndef Analyzer2D_h
#define Analyzer2D_h

#pragma once

#include <iostream>
#include <iomanip>
#include <string>

#include <ROOT/RDataFrame.hxx>
#include <memory>
#include "IAnalyzer.h"

namespace NuMI {

  static int nAna2D;
  class Ana2DHelper : public IAnalyzer<TH2D, 1>, public ROOT::Detail::RDF::RActionImpl<Ana2DHelper>
  {
  public:
    // we need these
    using Result_t = TH2D;
  protected:
    std::vector<std::shared_ptr<TH2D>> fHistos; // one per data processing slot per NDIM
    std::shared_ptr<TH2D> fHistoResult; // one per NDIM
    unsigned int fThreads;
    unsigned int fDim;

  public:
    Ana2DHelper(const TH2D& dummy, unsigned int nThreads=1)
      : IAnalyzer<TH2D, 1>(dummy, nThreads), fThreads(nThreads)
    {
      nAna2D++;
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

      // now initialize my histograms
      fHistoResult = std::make_shared<TH2D>(TString::Format("%s_id%d_result", dummy.GetName(), nAna2D), dummy.GetTitle(),
                                            dummy.GetNbinsX(), xarr.GetArray(), dummy.GetNbinsY(), yarr.GetArray());

      for(unsigned int i = 0; i < nSlots; i++) {
        fHistos.emplace_back(std::make_shared<TH2D>(
                             TString::Format("%s_id%d_slot%d", dummy.GetName(), nAna2D, i), dummy.GetTitle(),
                             dummy.GetNbinsX(), xarr.GetArray(), dummy.GetNbinsY(), yarr.GetArray())
        );
      }
    }
    // don't allow copies
    Ana2DHelper(Ana2DHelper &&) = default;
    Ana2DHelper(const Ana2DHelper &) = delete;

    std::shared_ptr<TH2D> GetResultPtr() const { return fHistoResult; }

    // we need these for df.Book()
    void Initialize() {}
    void InitTask(TTreeReader *, unsigned int) {}

    template <typename... DataTypes>
    void Exec(unsigned int slot, DataTypes... values)
    {
      // first one is val, second one is weight
      // cast everything into double
      std::array<double, sizeof...(DataTypes)> valuesArr{static_cast<double>(values)...};

      if(valuesArr.size() == 2){
        fHistos[slot]->Fill(valuesArr[0], valuesArr[1]);
      }
      else if(valuesArr.size() == 3){
        fHistos[slot]->Fill(valuesArr[0], valuesArr[1], valuesArr[2]);
      }
      else{
       std::cerr << "Expect atleast two columns (varx, vary) as input!! Exiting!" << std::endl;
       std::exit(1);
      }

    }

    void Finalize() const
    {
      for(unsigned int i = 0; i < fHistos.size(); i++)
        fHistoResult->Add(fHistos[i].get());
    }

    std::string GetActionName(){
      return "Ana2DHelper";
    }
  };
  // end helper
}

#endif
