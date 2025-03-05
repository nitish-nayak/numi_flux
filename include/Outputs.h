#ifndef Outputs_h
#define Outputs_h

#pragma once

#include <string>
#include <vector>

#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"

struct RootOutput
{
  RootOutput(std::string outname="flux.root");
  ~RootOutput(){
    fout->Close();
  };

  void Write();

  Float_t nuE       = -5.;          // neutrino energy
  Float_t nudirX    = -5.;          // neutrino X direction
  Float_t nudirY    = -5.;          // neutrino Y direction
  Float_t nudirZ    = -5.;          // neutrino Z direction
  Float_t nuL       = -5.;          // neutrino propagation length
  // neutrino energies at various locations (should be same as nuE but re-stored for redundancies)
  Float_t E_novand  = -5.;          // neutrino energy
  Float_t E_minosnd = -5.;          // neutrino energy
  Float_t E_minerva = -5.;          // neutrino energy

  Float_t wgt       = 1.;           // weight at location (smeared uboone), includes correction for muon polarization
  Float_t wgt_ppfx  = 1.;           // ppfx CV weight when available (only for dk2nu)
  // breakdowns of ppfx weights
  Float_t wgt_tgtatt = 1.;
  Float_t wgt_absorp = 1.;
  Float_t wgt_ttpcpion = 1.;
  Float_t wgt_ttpckaon = 1.;
  Float_t wgt_ttpcnucleon = 1.;
  Float_t wgt_ttncpion = 1.;
  Float_t wgt_ttnucleona = 1.;
  Float_t wgt_ttmesoninc = 1.;
  Float_t wgt_others = 1.;
  Float_t wgt_novand  = 1.;           // weight at location, includes correction for muon polarization
  Float_t wgt_minosnd = 1.;           // weight at location, includes correction for muon polarization
  Float_t wgt_minerva = 1.;           // weight at location, includes correction for muon polarization
  Float_t wgt_lp1 = 1.;
  Float_t wgt_lp2 = 1.;
  Float_t wgt_th = 1.;
  Float_t wgt_lp3 = 1.;
  Float_t wgt_shs = 1.;

  // get ppfx multiverse weights as well
  std::vector<Float_t> wgt_ppfxunivs;

  Int_t ptype       = -1000;        // parent particle code. uses GEANT codes for flugg, PDG for dk2nu
  Int_t gptype      = -1000;        // grandparent particle code. uses GEANT codes for flugg, PDG for dk2nu
  Int_t ncascade    = -1;           // number of interactions before neutrino, 2 = primary, etc
  Int_t ntype       = -1000;        // neutrino pdg code
  Int_t pmedium     = -1000;        // material codes for neutrino parent projectiles (See MINOS docdb: 9070)
  Int_t decaytype   = -1000;        // type of decay (See MINOS docdb: 9070)

  Float_t pE        = -5.;          // parent energy in lab frame
  Float_t pPt       = -5.;          // parent pT in lab frame
  Float_t pPz       = -5.;          // parent pZ in lab frame
  Float_t pTheta    = -5.;          // parent angle (in radians) wrt beamline (z-axis)
  Float_t pScTheta  = -5.;          // parent scattering angle (in radians) wrt its ancestor direction
  Float_t pxF       = -5.;          // parent Feynman-x, only available for primary in flugg, dk2nu has for any parent
  Float_t pxF_inc   = -5.;          // parent Feynman-x, assuming incident particle is 120 GeV proton beam
  Float_t pvx       = -5.;          // parent production point x (cm)
  Float_t pvy       = -5.;          // parent production point y (cm)
  Float_t pvz       = -5.;          // parent production point z (cm)
  Float_t gpE       = -5.;          // grandparent energy in lab frame
  Float_t gpPt      = -5.;          // grandparent pT in lab frame
  Float_t gpPz      = -5.;          // grandparent pZ in lab frame
  Float_t gpTheta   = -5.;          // grandparent angle (in radians) wrt beamline (z-axis)
  Float_t gpxF      = -5.;          // grandparent Feynman-x, only available for primary in flugg, dk2nu has for any parent
  Float_t nvx       = -5.;          // neutrino production point x (cm)
  Float_t nvy       = -5.;          // neutrino production point y (cm)
  Float_t nvz       = -5.;          // neutrino production point z (cm)

  std::string pProc = "NotFilled";            // parent process name

  TFile* fout;

  TH1D* numuFluxHisto;
  TH1D* anumuFluxHisto;
  TH1D* nueFluxHisto;
  TH1D* anueFluxHisto;
  TH1D* hPOT;
  TTree* outTree;

  double histMin    = 0.;
  double histMax    = 20.;
  int histNbins     = 4000;
};
#endif
