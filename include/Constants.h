#ifndef Constants_h
#define Constants_h

#include "TMath.h"

const double kNominalPOT = 6.e20;
const double kDefaultWeightCorrection = 1./(10000. * TMath::Pi());
const double kNTarget = 4.76e31/56.41e6 * 256.35 * 233 * 1036.8; //TPC active!!!

#ifdef HISTORIC_MASS
  const double kPIMASS      = 0.13957;
  const double kKMASS       = 0.49368;
  const double kK0MASS      = 0.49767;
  const double kMUMASS      = 0.105658389;
  const double kOMEGAMASS   = 1.67245;
#else
  const double kPIMASS      = 0.1395701;     // 0.13957;
  const double kKMASS       = 0.493677;      // 0.49368;
  const double kK0MASS      = 0.497614;      // 0.49767;
  const double kMUMASS      = 0.1056583715;  // 0.105658389;
  const double kOMEGAMASS   = 1.67245;       // 1.67245;
#endif

// from CLHEP/Units/PhysicalConstants.h
// used by Geant as CLHEP::neutron_mass_c2
const double kNEUTRONMASS = 0.93956536;
const int kpdg_nue       =   12;  // extended Geant 53
const int kpdg_nuebar    =  -12;  // extended Geant 52
const int kpdg_numu      =   14;  // extended Geant 56
const int kpdg_numubar   =  -14;  // extended Geant 55
const int kpdg_muplus      =   -13;  // Geant  5
const int kpdg_muminus     =    13;  // Geant  6
const int kpdg_pionplus    =   211;  // Geant  8
const int kpdg_pionminus   =  -211;  // Geant  9
const int kpdg_k0long      =   130;  // Geant 10  ( K0=311, K0S=310 )
const int kpdg_k0short     =   310;  // Geant 16
const int kpdg_k0mix       =   311;
const int kpdg_kaonplus    =   321;  // Geant 11
const int kpdg_kaonminus   =  -321;  // Geant 12
const int kpdg_omegaminus  =  3334;  // Geant 24
const int kpdg_omegaplus   = -3334;  // Geant 32
const int kpdg_neutron     =  2112;
const int kpdg_antineutron = -2112;


const int kgeant_nue       =   53;  // extended Geant 53
const int kgeant_nuebar    =  52;  // extended Geant 52
const int kgeant_numu      =   56;  // extended Geant 56
const int kgeant_numubar   =  55;  // extended Geant 55
const int kgeant_muplus     =   5;  // Geant  5
const int kgeant_muminus    =    6;  // Geant  6
const int kgeant_pionplus   =   8;  // Geant  8
const int kgeant_pionminus  =  9;  // Geant  9
const int kgeant_k0long     =   10;  // Geant 10  ( K0=311, K0S=310 )
const int kgeant_k0short    =   16;  // Geant 16
const int kgeant_k0mix      =   311;
const int kgeant_kaonplus   =   11;  // Geant 11
const int kgeant_kaonminus  =  12;  // Geant 12
const int kgeant_omegaminus =  24;  // Geant 24
const int kgeant_omegaplus  = 32;  // Geant 32

const double kRDET = 100.0;   // set to flux per 100 cm radius
