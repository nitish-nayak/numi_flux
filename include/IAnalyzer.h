#ifndef IAnalyzer_h
#define IAnalyzer_h

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <memory>

namespace NuMI {

  template <typename T, unsigned int NDIM>
  class IAnalyzer
  {
  private:
    std::vector<std::vector<std::shared_ptr<T>>> fHistos; // one per data processing slot per NDIM
    std::vector<std::shared_ptr<T>> fHistoResult; // one per NDIM
  protected:
    const unsigned int fDim;

  public:
    IAnalyzer(const T& dummy, const unsigned int nThreads=1)
      : fDim(NDIM)
    {
      if(fDim < 1){
        std::cerr << "Class template requires atleast >= 1 NDIM" << std::endl;
        std::exit(1);
      }
      fHistos = std::vector<std::vector<std::shared_ptr<T>>>(fDim);
      fHistoResult = std::vector<std::shared_ptr<T>>(fDim);
    }
    // don't allow copies
    IAnalyzer(IAnalyzer &&) = default;
    IAnalyzer(const IAnalyzer &) = delete;
  };
  // end helper
}

#endif
