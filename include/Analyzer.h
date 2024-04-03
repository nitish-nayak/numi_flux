#ifndef Analyzer_h
#define Analyzer_h

#pragma once

#include "Outputs.h"

#include <iostream>
#include <iomanip>
#include <string>

#include <ROOT/RDataFrame.hxx>
#include <memory>

namespace NuMI {

  static int AnaID;
  class Ana1DHelper : public ROOT::Detail::RDF::RActionImpl<Ana1DHelper>
  {
  public:
    // we need these
    using Result_t = TH1D;
  private:
    std::vector<std::shared_ptr<TH1D>> fHistos; // one per data processing slot
    std::shared_ptr<TH1D> fHistoResult;
    const unsigned int fThreads;

  public:
    Ana1DHelper(const TH1D& dummy, const unsigned int nThreads=1)
      : fThreads(nThreads)
    {
      AnaID++;
      const unsigned int nSlots = ROOT::IsImplicitMTEnabled() ? fThreads : 1;

      const double nbins = dummy.GetNbinsX();
      const double xlow  = dummy.GetBinLowEdge(1);
      const double xhigh = dummy.GetBinLowEdge(nbins+1);
      fHistoResult = std::make_shared<TH1D>(TString::Format("%s_id%d_result", dummy.GetName(), AnaID), dummy.GetTitle(),
                                            nbins, xlow, xhigh);

      for(unsigned int i = 0; i < nSlots; i++) {
        fHistos.emplace_back(std::make_shared<TH1D>(
                             TString::Format("%s_id%d_slot%d", dummy.GetName(), AnaID, i), dummy.GetTitle(),
                             nbins, xlow, xhigh)
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

    void Finalize()
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
