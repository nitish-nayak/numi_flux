#define Dk2NuFlux_cxx

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/calcLocationWeights.h"

// Nitish : Add ability to read in PPFX libraries
#include "MakeReweight.h"

using namespace NeutrinoFluxReweight;

#include <string>
#include <cstddef>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <iomanip>

using namespace std;

#include "TChain.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TRotation.h"
#include "TMath.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include <TParticlePDG.h>
#include <TDatabasePDG.h>

#include "Dk2NuFlux.h"

//___________________________________________________________________________
Dk2NuFlux::Dk2NuFlux(std::string pattern, std::string outfile)
{
  cflux = new TChain("dk2nuTree");
  cflux->Add(pattern.c_str());

  cflux_meta = new TChain("dkmetaTree");
  cflux_meta->Add(pattern.c_str());

  Int_t nfiles = cflux->GetNtrees();
  std::cout << "Number of files: " << nfiles << std::endl;

  fOutput = new RootOutput(outfile);
}

//___________________________________________________________________________
void Dk2NuFlux::CalculateFlux()
{
  // initialize ppfx
  MakeReweight* fPPFXrw = MakeReweight::getInstance();
  std::string fPPFXMode = "ubnumi_cvonly";
  gSystem->Setenv("MODE", fPPFXMode.c_str());
  Int_t fSeed = 84;
  fPPFXrw->setBaseSeed(fSeed); // used in EventWeight ubsim, shouldn't matter for CV
  std::string inputOptions = std::string(std::getenv("NUMIANA_DIR"))+"/dk2nu/ppfx/inputs_"+fPPFXMode+".xml";
  if(!(fPPFXrw->AlreadyInitialized())){
    fPPFXrw->SetOptions(inputOptions);
  }
  std::cout << "Setup PPFX with seed " << fSeed << " from : " << inputOptions << std::endl;


  bsim::Dk2Nu* fDk2Nu = new bsim::Dk2Nu;
  bsim::DkMeta* fDkMeta= new bsim::DkMeta;
  cflux->SetBranchAddress("dk2nu", &fDk2Nu);
  cflux_meta->SetBranchAddress("dkmeta",&fDkMeta);

  //***************************************
  //
  //  Loop over the entries.
  //
  //***************************************

  int treeNumber = -1;
  Long64_t nflux = cflux->GetEntries();
  std::cout << "Total number of entries: " << nflux << std::endl;

  for (Long64_t i=0; i < nflux; ++i ) {

    // Get entry i. fluxNtuple is now filled with entry i info.
    cflux->GetEntry(i);
    cflux_meta->GetEntry(i);

    // Alert the user
    if (i % 100000 == 0) std::cout << "On entry " << i/1.0e6 << " M"<< std::endl;

    if(treeNumber != cflux->GetTreeNumber()) {
      treeNumber = cflux->GetTreeNumber();
      std::cout << "Moving to tree number " << treeNumber << "." << std::endl;
      AccumulatedPOT += fDkMeta->pots;
      std::cout << "AccumulatedPOT: " << AccumulatedPOT << std::endl;
    }

    double wgt_xy = 0.;  // neutrino weight
    double enu    = 0.;  // neutrino energy in lab frame

    // Pick a random point in the TPC (in detector coordinates)
    TVector3 xyz_det = RandomInTPC();
    if (fDebug) std::cout << "xyz_det = [" << xyz_det.X() << ", " <<
                                              xyz_det.Y() << ", " <<
                                              xyz_det.Z() << "]" << std::endl;

    // From detector to beam coordinates
    TVector3 xyz_beam = FromDetToBeam(xyz_det);
    if (fDebug) std::cout << "xyz_beam = [" << xyz_beam.X() << ", " <<
                                               xyz_beam.Y() << ", " <<
                                               xyz_beam.Z() << "]" << std::endl;

    TVector3 nudir_unit = xyz_beam; // get the neutrino direction in beam coordinates
    nudir_unit -= TVector3(fDk2Nu->decay.vx, fDk2Nu->decay.vy, fDk2Nu->decay.vz);
    fOutput->nuL = nudir_unit.Mag();

    nudir_unit = nudir_unit.Unit();

    // Neutrons are not yet implemented so skip for now
    // if (fDk2Nu->decay.ptype == 2112) continue;

    // Calculate the weight
    int ret = CalculateWeight(fDk2Nu, xyz_beam, enu, wgt_xy);
    if (ret != 0) std::cout << "Error with CalculateWeight. Return " << ret << std::endl;
    if (fDebug) std::cout << "wgt_xy " << wgt_xy << std::endl;


    // Calculate the total weight
    double weight = wgt_xy * fDk2Nu->decay.nimpwt * kDefaultWeightCorrection;

    if (std::isnan(weight) == 1) { // catch NaN values
      std::cout << "got a nan: wgt\t" << weight << std::endl;
      weight = 0;
    }

    if (std::isnan(enu) == 1) { // catch NaN values
      std::cout << "got a nan enu:\t" << enu << std::endl;
      enu = 0;
    }

    // Fill the histograms
    switch (fDk2Nu->decay.ntype) {
      case kpdg_numu:
        fOutput->numuFluxHisto->Fill(enu, weight);
        break;
      case kpdg_numubar:
        fOutput->anumuFluxHisto->Fill(enu, weight);
        break;
      case kpdg_nue:
        fOutput->nueFluxHisto->Fill(enu, weight);
        break;
      case kpdg_nuebar:
        fOutput->anueFluxHisto->Fill(enu, weight);
        break;
    }

    // get ppfx weights
    try {
      fPPFXrw->calculateWeights(fDk2Nu, fDkMeta);
      fOutput->wgt_ppfx = fPPFXrw->GetCVWeight();
    } catch (...) {
      std::cout<<"Failed to calculate wgt"<<std::endl;
      fOutput->wgt_ppfx = 1.;
    }

    // parent information
    TLorentzVector pvec = TLorentzVector((fDk2Nu->decay.ppdxdz)*(fDk2Nu->decay.pppz), (fDk2Nu->decay.ppdydz)*(fDk2Nu->decay.pppz),
                                          fDk2Nu->decay.pppz, fDk2Nu->decay.ppenergy);
    double inc_pbeam_mom = std::sqrt(pow(120., 2.) - pow(kPROTONMASS, 2.));
    TLorentzVector gpvec = TLorentzVector(0., 0., inc_pbeam_mom, 120.);
    std::vector<bsim::Ancestor> ancestors = fDk2Nu->ancestor;

    fOutput->pE = NuMI::E(pvec);
    fOutput->pPt = NuMI::Pt(pvec);
    fOutput->pPz = NuMI::Pz(pvec);
    fOutput->pTheta = TMath::ACos(NuMI::CosTheta(pvec));
    fOutput->pxF_inc = NuMI::xF(gpvec, pvec);

    // outsource interaction chain to ppfx
    InteractionChainData icd(fDk2Nu, fDkMeta);
    int nicd = icd.interaction_chain.size();

    if(nicd >= 2) {

      const InteractionData id = icd.interaction_chain[nicd-2];
      double ScPt = id.Pt;
      double ScPz = id.Pz;
      fOutput->pScTheta = TMath::ACos(ScPz/std::sqrt(pow(ScPt, 2.) + pow(ScPz, 2.)));
      fOutput->pxF = id.xF;

      gpvec = TLorentzVector(id.Inc_P4[0], id.Inc_P4[1], id.Inc_P4[2], id.Inc_P4[3]);
      fOutput->gpE = NuMI::E(gpvec);
      fOutput->gpPt = NuMI::Pt(gpvec);
      fOutput->gpPz = NuMI::Pz(gpvec);
      fOutput->gpTheta = TMath::ACos(NuMI::CosTheta(gpvec));

      fOutput->pProc = id.Proc;
      fOutput->pvx   = id.Vtx[0];
      fOutput->pvy   = id.Vtx[1];
      fOutput->pvz   = id.Vtx[2];
    }

    // // for primary particles, keep flugg calculation (even if slightly in-accurate)
    // if(fDk2Nu->tgtexit.tgen == 2)
    //   fOutput->pxF = NuMI::xF(gpvec, pvec);
    // if(fDk2Nu->tgtexit.tgen > 2) {
    //   unsigned int pidx = ancestors.size() - 2;
    //   unsigned int gpidx = ancestors.size() - 3;
    //
    //   double pmass = NuMI::MassFromPdgCode(ancestors[pidx].pdg);
    //   double gpmass = NuMI::MassFromPdgCode(ancestors[gpidx].pdg);
    //
    //   double pmom = std::sqrt(pow(ancestors[pidx].startpx, 2.) + pow(ancestors[pidx].startpy, 2.) +
    //                           pow(ancestors[pidx].startpz, 2.));
    //   double gpmom = std::sqrt(pow(ancestors[gpidx].startpx, 2.) + pow(ancestors[gpidx].startpy, 2.) +
    //                            pow(ancestors[gpidx].startpz, 2.));
    //
    //   pvec = TLorentzVector(ancestors[pidx].startpx, ancestors[pidx].startpy,
    //                         ancestors[pidx].startpz, std::sqrt(pow(pmom, 2.) + pow(pmass, 2.)));
    //   gpvec = TLorentzVector(ancestors[gpidx].startpx, ancestors[gpidx].startpy,
    //                          ancestors[gpidx].startpz, std::sqrt(pow(gpmom, 2.) + pow(gpmass, 2.)));
    //   fOutput->pxF = NuMI::xF(gpvec, pvec);
    // }

    fOutput->nuE = enu;
    fOutput->nudirX = nudir_unit.X();
    fOutput->nudirY = nudir_unit.Y();
    fOutput->nudirZ = nudir_unit.Z();

    fOutput->nvx = fDk2Nu->decay.vx;
    fOutput->nvy = fDk2Nu->decay.vy;
    fOutput->nvz = fDk2Nu->decay.vz;
    fOutput->wgt = weight;
    fOutput->ptype = fDk2Nu->decay.ptype;
    fOutput->ntype = fDk2Nu->decay.ntype;
    fOutput->ncascade = fDk2Nu->tgtexit.tgen;
    fOutput->pmedium = fDk2Nu->decay.ppmedium;
    fOutput->decaytype = fDk2Nu->decay.ndecay;

    (fOutput->outTree)->Fill();

  } // end of loop over the entries


  //***************************************
  //
  // POT scaling
  //
  //***************************************

  double scale =    kNominalPOT/AccumulatedPOT;
  (fOutput->hPOT)->SetBinContent(1, AccumulatedPOT);
  (fOutput->numuFluxHisto)  -> Scale(scale);
  (fOutput->anumuFluxHisto) -> Scale(scale);
  (fOutput->nueFluxHisto)   -> Scale(scale);
  (fOutput->anueFluxHisto)  -> Scale(scale);

  std::cout << std::endl << ">>> TOTAL POT: " << AccumulatedPOT << std::endl << std::endl;

  //***************************************
  //
  // Writing on file
  //
  //***************************************
  fOutput->Write();

  // reset ppfx state
  fPPFXrw->resetInstance();

}

//___________________________________________________________________________
int Dk2NuFlux::CalculateWeight(bsim::Dk2Nu* decay, const TVector3& xyz,double& enu, double& wgt_xy)
{
  // Neutrino Energy and Weight at arbitrary point
  // Based on:
  //   NuMI-NOTE-BEAM-0109 (MINOS DocDB # 109)
  //   Title:   Neutrino Beam Simulation using PAW with Weighted Monte Carlos
  //   Author:  Rick Milburn
  //   Date:    1995-10-01
  // History:
  // jzh  3/21/96 grab R.H.Milburn's weighing routine
  // jzh  5/ 9/96 substantially modify the weighting function use dot product
  //              instead of rotation vecs to get theta get all info except
  //              det from ADAMO banks neutrino parent is in Particle.inc
  //              Add weighting factor for polarized muon decay
  // jzh  4/17/97 convert more code to double precision because of problems
  //              with Enu>30 GeV
  // rwh 10/ 9/08 transliterate function from f77 to C++
  // Original function description:
  //   Real function for use with PAW Ntuple To transform from destination
  //   detector geometry to the unit sphere moving with decaying hadron with
  //   velocity v, BETA=v/c, etc..  For (pseudo)scalar hadrons the decays will
  //   be isotropic in this  sphere so the fractional area (out of 4-pi) is the
  //   fraction of decays that hit the target.  For a given target point and
  //   area, and given x-y components of decay transverse location and slope,
  //   and given decay distance from target ans given decay GAMMA and
  //   rest-frame neutrino energy, the lab energy at the target and the
  //   fractional solid angle in the rest-frame are determined.
  //   For muon decays, correction for non-isotropic nature of decay is done.
  // Arguments:
  //    dk2nu    :: contains current decay information
  //    xyz      :: 3-vector of position to evaluate
  //                in *beam* frame coordinates  (cm units)
  //    enu      :: resulting energy
  //    wgt_xy   :: resulting weight
  // Return:
  //    (int)    :: error code
  // Assumptions:
  //    Energies given in GeV
  //    Particle codes have been translated from GEANT into PDG codes
  // for now ... these masses _should_ come from TDatabasePDG
  // but use these hard-coded values to "exactly" reproduce old code
  //
  // old mass values are v01_07_* and before
  // new mass values (v01_08_* and after) are Geant4 v4.10.3 values
  //
  double xpos = xyz.X();
  double ypos = xyz.Y();
  double zpos = xyz.Z();
  enu    = 0.0;  // don't know what the final value is
  wgt_xy = 0.0;  // but set these in case we return early due to error
  // in principle we should get these from the particle DB
  // but for consistency testing use the hardcoded values
  double parent_mass = NuMI::MassFromPdgCode(decay->decay.ptype);

  double parentp2 = ( decay->decay.pdpx*decay->decay.pdpx +
                      decay->decay.pdpy*decay->decay.pdpy +
                      decay->decay.pdpz*decay->decay.pdpz );

  double parent_energy = TMath::Sqrt( parentp2 + parent_mass*parent_mass);

  double parentp = TMath::Sqrt( parentp2 );

  double gamma     = parent_energy / parent_mass;
  double gamma_sqr = gamma * gamma;

  double beta_mag  = TMath::Sqrt( ( gamma_sqr - 1.0 )/gamma_sqr );

  // Get the neutrino energy in the parent decay CM
  double enuzr = decay->decay.necm;

  // Get angle from parent line of flight to chosen point in beam frame
  double rad = TMath::Sqrt( (xpos-decay->decay.vx)*(xpos-decay->decay.vx) +
                          (ypos-decay->decay.vy)*(ypos-decay->decay.vy) +
                          (zpos-decay->decay.vz)*(zpos-decay->decay.vz) );

  double emrat = 1.0;
  double costh_pardet = -999., theta_pardet = -999.;

  // Boost correction, but only if parent hasn't stopped
  if ( parentp > 0. ) {
    costh_pardet = ( decay->decay.pdpx*(xpos-decay->decay.vx) +
                     decay->decay.pdpy*(ypos-decay->decay.vy) +
                     decay->decay.pdpz*(zpos-decay->decay.vz) ) / ( parentp * rad);

    if ( costh_pardet >  1.0 ) costh_pardet =  1.0;

    if ( costh_pardet < -1.0 ) costh_pardet = -1.0;

    theta_pardet = TMath::ACos(costh_pardet);

    // Weighted neutrino energy in beam, approx, good for small theta
    emrat = 1.0 / ( gamma * ( 1.0 - beta_mag * costh_pardet ));
  }

  enu = emrat * enuzr;  // the energy ... normally

  // Get solid angle/4pi for detector element
  // small angle approximation, fixed by Alex Radovic
  //SAA//  double sangdet = ( kRDET*kRDET /
  //SAA//                   ( (zpos-decay.vz)*(zpos-decay.vz) ) ) / 4.0;
  double sanddetcomp = TMath::Sqrt( ( (xpos-decay->decay.vx)*(xpos-decay->decay.vx) ) +
                                  ( (ypos-decay->decay.vy)*(ypos-decay->decay.vy) ) +
                                  ( (zpos-decay->decay.vz)*(zpos-decay->decay.vz) )   );

  double sangdet = (1.0-TMath::Cos(TMath::ATan( kRDET / sanddetcomp )))/2.0;

  // Weight for solid angle and lorentz boost
  wgt_xy = sangdet * ( emrat * emrat );  // ! the weight ... normally

  // Done for all except polarized muon decay
  // in which case need to modify weight
  // (must be done in double precision)
  // BUT do this only for case of muon decay, not muon capture
  //   until beamline simulation code gets updated these generally show up as
  //   decay.ndecay == 0, but certainly not dkp_mup_nusep or dkp_mum_nusep
  // so was:
  // if ( decay.ptype  == kpdg_muplus      ||
  //      decay.ptype  == kpdg_muminus        ) {
  // now:
  if ( decay->decay.ndecay == bsim::dkp_mup_nusep || decay->decay.ndecay  == bsim::dkp_mum_nusep    ) {
    double beta[3], p_dcm_nu[4], p_nu[3], p_pcm_mp[3], partial;

    // Boost neu neutrino to mu decay CM
    beta[0] = decay->decay.pdpx / parent_energy;
    beta[1] = decay->decay.pdpy / parent_energy;
    beta[2] = decay->decay.pdpz / parent_energy;
    p_nu[0] = (xpos-decay->decay.vx)*enu/rad;
    p_nu[1] = (ypos-decay->decay.vy)*enu/rad;
    p_nu[2] = (zpos-decay->decay.vz)*enu/rad;

    partial = gamma * (beta[0]*p_nu[0] + beta[1]*p_nu[1] + beta[2]*p_nu[2] );
    partial = enu - partial/(gamma+1.0);

    // the following calculation is numerically imprecise
    // especially p_dcm_nu[2] leads to taking the difference of numbers
    //  of order ~10's and getting results of order ~0.02's
    // for g3numi we're starting with floats (ie. good to ~1 part in 10^7)
    p_dcm_nu[0] = p_nu[0] - beta[0]*gamma*partial;
    p_dcm_nu[1] = p_nu[1] - beta[1]*gamma*partial;
    p_dcm_nu[2] = p_nu[2] - beta[2]*gamma*partial;
    p_dcm_nu[3] = TMath::Sqrt( p_dcm_nu[0]*p_dcm_nu[0] +
                               p_dcm_nu[1]*p_dcm_nu[1] +
                               p_dcm_nu[2]*p_dcm_nu[2] );

    // Boost parent of mu to mu production CM
    double particle_energy = decay->decay.ppenergy;
    gamma = particle_energy/parent_mass;

    beta[0] = decay->decay.ppdxdz * decay->decay.pppz / particle_energy;
    beta[1] = decay->decay.ppdydz * decay->decay.pppz / particle_energy;
    beta[2] = decay->decay.pppz / particle_energy;

    partial = gamma * ( beta[0]*decay->decay.muparpx +
                        beta[1]*decay->decay.muparpy +
                        beta[2]*decay->decay.muparpz );
    partial = decay->decay.mupare - partial/(gamma+1.0);

    p_pcm_mp[0] = decay->decay.muparpx - beta[0]*gamma*partial;
    p_pcm_mp[1] = decay->decay.muparpy - beta[1]*gamma*partial;
    p_pcm_mp[2] = decay->decay.muparpz - beta[2]*gamma*partial;

    double p_pcm = TMath::Sqrt ( p_pcm_mp[0]*p_pcm_mp[0] +
                                 p_pcm_mp[1]*p_pcm_mp[1] +
                                 p_pcm_mp[2]*p_pcm_mp[2] );

    const double eps = 1.0e-30;  // ? what value to use

    if ( p_pcm < eps || p_dcm_nu[3] < eps ) {
      return 3; // mu missing parent info?
    }

    // Calc new decay angle w.r.t. (anti)spin direction
    double costh = ( p_dcm_nu[0]*p_pcm_mp[0] +
                        p_dcm_nu[1]*p_pcm_mp[1] +
                        p_dcm_nu[2]*p_pcm_mp[2] ) /
                    ( p_dcm_nu[3]*p_pcm );

    // protect against small excursions
    if ( costh >  1.0 ) costh =  1.0;
    if ( costh < -1.0 ) costh = -1.0;

    // Calc relative weight due to angle difference
    double wgt_ratio = 0.0;
    switch ( decay->decay.ntype ) {
      case kpdg_nue:
      case kpdg_nuebar:
          wgt_ratio = 1.0 - costh;
          break;

      case kpdg_numu:
      case kpdg_numubar:
      {
        double xnu = 2.0 * enuzr / kMUMASS;
        wgt_ratio = ( (3.0-2.0*xnu )  - (1.0-2.0*xnu)*costh ) / (3.0-2.0*xnu);

        if ( wgt_ratio < 0.0 ) {
          std::cerr << "bsim::CalculateWeight encountered serious problem: "
                    << " wgt_ratio " << wgt_ratio
                    << " enu " << enu << " costh " << costh << " xnu " << xnu
                    << " enuzr=decay->decay.necm " << enuzr << " kMUMASS " << kMUMASS
                    << " norig " << decay->decay.norig
                    << " ndecay->decay " << decay->decay.ndecay
                    << " ntype " << decay->decay.ntype
                    << " ptype " << decay->decay.ptype
                    << std::endl;
          enu    = 0;
          wgt_xy = 0;
          return 4; // bad, bad, bad calculation
        }

        break;
      }

      default:
        enu    = 0.0;
        wgt_xy = 0.0;
        return 2; // bad neutrino type for muon decay
      }

      wgt_xy = wgt_xy * wgt_ratio;

  } // ptype is muon

  return 0;

}
