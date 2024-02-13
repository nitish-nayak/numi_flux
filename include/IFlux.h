#ifndef IFlux_h
#define IFlux_h

#pragma once

#include "Constants.h"

#include <iostream>
#include <iomanip>
#include <string>

#include "TVector3.h"
#include "TLorentzVector.h"

namespace NuMI {

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

  // define some useful inlines
  //___________________________________________________________________________
  inline double CosTheta(const TLorentzVector& vec)
  {
    return vec.CosTheta();
  }

  //___________________________________________________________________________
  inline double SinTheta(const TLorentzVector& vec)
  {
    return std::sqrt(1. - pow(CosTheta(vec), 2.0));
  }

  //___________________________________________________________________________
  inline double Pt(const TLorentzVector& vec)
  {
    return vec.Pt();
  }

  //___________________________________________________________________________
  inline double Pz(const TLorentzVector& vec)
  {
    return vec.Pz();
  }

  //___________________________________________________________________________
  inline double E(const TLorentzVector& vec)
  {
    return vec.E();
  }

  //___________________________________________________________________________
  inline double xF(const TLorentzVector& inc_vec, const TLorentzVector& prod_vec)
  {
    double Inc_P = inc_vec.P();
    double Prod_P = prod_vec.P();

    double Inc_E = inc_vec.E();
    double Prod_E = prod_vec.E();

    TVector3 Inc_Dir = inc_vec.BoostVector().Unit();
    TVector3 Prod_Dir = prod_vec.BoostVector().Unit();
    double cos_theta = Inc_Dir.Dot(Prod_Dir);

    double Ecm = std::sqrt(inc_vec.M2() + pow(kNUCLEONMASS, 2.) +
                           2.*Inc_E*kNUCLEONMASS);
    double Betacm = std::sqrt(pow(Inc_E, 2.) - inc_vec.M2()) /
                    (Inc_E + kNUCLEONMASS);

    double Gammacm = 1. / std::sqrt(1. - pow(Betacm, 2.));
    double PLcm = Gammacm * (Prod_P*cos_theta - Betacm*Prod_E);

    return PLcm * 2. / Ecm;
  }

  //___________________________________________________________________________
  inline double MassFromGeantCode(const int geant_code)
  {
    double mass = kPIMASS;
    switch ( geant_code ) {
        case kgeant_pionplus:
        case kgeant_pionminus:
            mass = kPIMASS;
            break;
        case kgeant_kaonplus:
        case kgeant_kaonminus:
            mass = kKMASS;
            break;
        case kgeant_k0long:
        case kgeant_k0short:
        case kgeant_k0mix:
            mass = kK0MASS;
            break;
        case kgeant_muplus:
        case kgeant_muminus:
            mass = kMUMASS;
            break;
        case kgeant_proton:
            mass = kPROTONMASS;
            break;
        case kgeant_neutron:
            mass = kNEUTRONMASS;
            break;
        default:
            mass = kOMEGAMASS;
            break;
    }
    return mass;
  }

  //___________________________________________________________________________
  inline double MassFromPdgCode(const int pdg_code)
  {
    double mass = kPIMASS;
    switch ( pdg_code ) {
        case kpdg_pionplus:
        case kpdg_pionminus:
            mass = kPIMASS;
            break;
        case kpdg_kaonplus:
        case kpdg_kaonminus:
            mass = kKMASS;
            break;
        case kpdg_k0long:
        case kpdg_k0short:
        case kpdg_k0mix:
            mass = kK0MASS;
            break;
        case kpdg_muplus:
        case kpdg_muminus:
            mass = kMUMASS;
            break;
        case kpdg_proton:
            mass = kPROTONMASS;
            break;
        case kpdg_neutron:
            mass = kNEUTRONMASS;
            break;
        default:
            mass = kOMEGAMASS;
            break;
    }
    return mass;
  }
}

#endif
