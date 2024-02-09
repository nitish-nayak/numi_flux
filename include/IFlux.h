#ifndef IFlux_h
#define IFlux_h

#pragma once

#include <iostream>
#include <iomanip>
#include <string>

#include "TVector3.h"

class IFlux
{
public:
  IFlux(){};
  virtual ~IFlux() {};

  TVector3 RandomInTPC() const;
  TVector3 FromDetToBeam(const TVector3& det) const;
  double EstimatePOT(int highest_potnum) const;

  virtual void CalculateFlux() = 0;
  int CalculateWeight() = delete;

  void SetDebug(){ fDebug = true; }

protected:
  bool fDebug = false;
};
#endif
