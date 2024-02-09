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
  IFlux() = delete;
  virtual ~IFlux() = delete;

  TVector3 RandomInTPC() const;
  TVector3 FromDetToBeam(const TVector3& det) const;
  double EstimatePOT(int highest_potnum) const;

  virtual void CalculateFlux() const = 0;
  virtual int CalculateWeight() = 0;

  void SetDebug(){ fDebug = true; }

private:
  bool fDebug = false;
}
