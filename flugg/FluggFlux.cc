#define FluggFlux_cxx

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>

#include "TChain.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TRotation.h"
#include "TMath.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"

#include "FluggFlux.h"

//___________________________________________________________________________
FluggFlux::FluggFlux(std::string pattern, std::string outfile)
{
  TString libs = gSystem->GetDynamicPath();
  libs += ":";
  libs += std::string(std::getenv("NUMIANA_DIR"))+"/flugg";
  gSystem->SetDynamicPath(libs.Data());
  gSystem->Load("FluggTree_C.so");

  cflux = new TChain("h10");
  cflux->Add(pattern.c_str());

  Int_t nfiles = cflux->GetNtrees();
  std::cout << "Number of files: " << nfiles << std::endl;

  fOutput = new RootOutput(outfile);
}

//___________________________________________________________________________
FluggFlux::FluggFlux(bool isfilelist, std::string filelist, std::string outfile)
{
  if (!isfilelist){
    std::cerr << "Require input isfilelist argument to be true" << std::endl;
    std::exit(1);
  }
  TString libs = gSystem->GetDynamicPath();
  libs += ":";
  libs += std::string(std::getenv("NUMIANA_DIR"))+"/flugg";
  gSystem->SetDynamicPath(libs.Data());
  gSystem->Load("FluggTree_C.so");

  cflux = new TChain("h10");
  std::ifstream f_stream;
  f_stream.open(filelist.c_str());
  std::string f_line;
  while (f_stream.good()) {
    std::getline(f_stream, f_line);
    if(f_line.find(".root") > 100000) continue;
    cflux->Add(f_line.c_str());
  }

  Int_t nfiles = cflux->GetNtrees();
  std::cout << "Number of files: " << nfiles << std::endl;

  fOutput = new RootOutput(outfile);
}

//___________________________________________________________________________
void FluggFlux::CalculateFlux()
{

  fluxNtuple = new FluggTree(cflux);

  //***************************************
  //
  //  Loop over the entries.
  //
  //***************************************
  int treeNumber = -1;
  int highest_evtno = 0;
  Long64_t nflux = cflux->GetEntries();
  std::cout << "Total number of entries: " << nflux << std::endl;

  for (Long64_t i=0; i < nflux; ++i ) {

    // Get entry i. fluxNtuple is now filled with entry i info.
    cflux->GetEntry(i);

    // Alert the user
    if (i % 100000 == 0) std::cout << "On entry " << i << std::endl;
    if(treeNumber != cflux->GetTreeNumber()) {
      treeNumber = cflux->GetTreeNumber();
      std::cout << "Moving to tree number " << treeNumber << "." << std::endl;
      AccumulatedPOT += EstimatePOT(highest_evtno);
      std::cout << "AccumulatedPOT: " << AccumulatedPOT << std::endl;
    }

    double wgt_xy = 0.;  // neutrino weight
    double enu    = 0.;  // neutrino energy in lab frame
    // I don't know why the neutrino energies would be different but the fact that these are also computed in the weight function worries me a bit
    // it should be the same as enu above based on the physical quantities used in the calculation (correction from parent COM frame + corrections for muon polarization etc)
    // but there are various failure modes where this gets reset, which I don't fully understand, so I'm just going to recompute things.. should be free of cost anyway
    double enu_novand    = 0.;  // neutrino energy in lab frame
    double enu_minosnd    = 0.;  // neutrino energy in lab frame
    double enu_minerva    = 0.;  // neutrino energy in lab frame
    double wgt_novand    = 0.;  // neutrino energy in lab frame
    double wgt_minosnd    = 0.;  // neutrino energy in lab frame
    double wgt_minerva    = 0.;  // neutrino energy in lab frame

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
    nudir_unit -= TVector3(fluxNtuple->Vx, fluxNtuple->Vy, fluxNtuple->Vz);
    fOutput->nuL = nudir_unit.Mag();

    nudir_unit = nudir_unit.Unit();

    //  the weight
    int ret = CalculateWeight(fluxNtuple, xyz_beam, enu, wgt_xy);
    if (ret != 0) std::cout << "Error with CalculateWeight. Return " << ret << std::endl;
    if (fDebug) std::cout << "wgt_xy " << wgt_xy << std::endl;

    int ret_novand = CalculateWeight(fluxNtuple, kNOvA_ND, enu_novand, wgt_novand);
    if (ret_novand != 0) std::cout << "Error with CalculateWeight. Return " << ret_novand << std::endl;
    if (fDebug) std::cout << "wgt_novand " << wgt_novand << std::endl;

    int ret_minosnd = CalculateWeight(fluxNtuple, kMINOS_ND, enu_minosnd, wgt_minosnd);
    if (ret_minosnd != 0) std::cout << "Error with CalculateWeight. Return " << ret_minosnd << std::endl;
    if (fDebug) std::cout << "wgt_minosnd " << wgt_minosnd << std::endl;

    int ret_minerva = CalculateWeight(fluxNtuple, kMINERvA, enu_minerva, wgt_minerva);
    if (ret_minerva != 0) std::cout << "Error with CalculateWeight. Return " << ret_minerva << std::endl;
    if (fDebug) std::cout << "wgt_minerva " << wgt_minerva << std::endl;

    //  the total weight
    double weight = wgt_xy * fluxNtuple->Nimpwt * kDefaultWeightCorrection;
    double weight_novand = wgt_novand * fluxNtuple->Nimpwt * kDefaultWeightCorrection;
    double weight_minosnd = wgt_minosnd * fluxNtuple->Nimpwt * kDefaultWeightCorrection;
    double weight_minerva = wgt_minerva * fluxNtuple->Nimpwt * kDefaultWeightCorrection;

    // Fill the histograms
    switch (fluxNtuple->Ntype) {
      case kgeant_numu:
        fOutput->numuFluxHisto->Fill(enu, weight);
        fOutput->ntype = kpdg_numu;
        break;
      case kgeant_numubar:
        fOutput->anumuFluxHisto->Fill(enu, weight);
        fOutput->ntype = kpdg_numubar;
        break;
      case kgeant_nue:
        fOutput->nueFluxHisto->Fill(enu, weight);
        fOutput->ntype = kpdg_nue;
        break;
      case kgeant_nuebar:
        fOutput->anueFluxHisto->Fill(enu, weight);
        fOutput->ntype = kpdg_nuebar;
        break;
    }

    // parent information
    TLorentzVector pvec = TLorentzVector((fluxNtuple->ppdxdz)*(fluxNtuple->pppz), (fluxNtuple->ppdydz)*(fluxNtuple->pppz),
                                         fluxNtuple->pppz, fluxNtuple->ppenergy);
    double inc_pbeam_mom = std::sqrt(pow(120., 2.) - pow(kPROTONMASS, 2.));
    TLorentzVector gpvec = TLorentzVector(0., 0., inc_pbeam_mom, 120.);

    fOutput->pE = NuMI::E(pvec);
    fOutput->pPt = NuMI::Pt(pvec);
    fOutput->pPz = NuMI::Pz(pvec);
    fOutput->pTheta = TMath::ACos(NuMI::CosTheta(pvec));
    fOutput->pxF_inc = NuMI::xF(gpvec, pvec);
    fOutput->pvx = fluxNtuple->ppvx;
    fOutput->pvy = fluxNtuple->ppvy;
    fOutput->pvz = fluxNtuple->ppvz;
    // trust xF only for primary particles
    if(fluxNtuple->tgen == 2)
      fOutput->pxF = NuMI::xF(gpvec, pvec);

    // POT stuff
    if ( fluxNtuple->evtno > highest_evtno )
      highest_evtno = fluxNtuple->evtno;

    fOutput->nuE = enu;
    fOutput->nudirX = nudir_unit.X();
    fOutput->nudirY = nudir_unit.Y();
    fOutput->nudirZ = nudir_unit.Z();

    fOutput->nvx = fluxNtuple->Vx;
    fOutput->nvy = fluxNtuple->Vy;
    fOutput->nvz = fluxNtuple->Vz;
    fOutput->wgt = weight;
    fOutput->ptype = fluxNtuple->ptype;
    fOutput->ncascade = fluxNtuple->tgen;
    fOutput->pmedium = fluxNtuple->ppmedium;
    fOutput->decaytype = fluxNtuple->Ndecay;

    fOutput->E_novand = enu_novand;
    fOutput->E_minosnd = enu_minosnd;
    fOutput->E_minerva = enu_minerva;
    fOutput->wgt_novand = weight_novand;
    fOutput->wgt_minosnd = weight_minosnd;
    fOutput->wgt_minerva = weight_minerva;

    (fOutput->outTree)->Fill();

  } // end of loop over the entries


  //***************************************
  //
  // POT scaling
  //
  //***************************************

  AccumulatedPOT += EstimatePOT(highest_evtno); // To account for last tree
  (fOutput->hPOT)->SetBinContent(1, AccumulatedPOT);

  double scale = kNominalPOT/AccumulatedPOT;
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
}

//___________________________________________________________________________
int FluggFlux::CalculateWeight(FluggTree* decay, const TVector3& xyz,
                               double& enu, double& wgt_xy)
{

  // Stolen: https://cdcvs.fnal.gov/redmine/projects/dk2nu/repository/show/trunk/dk2nu
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
  //   detector geometry to the unit sphere moving with decayng hadron with
  //   velocity v, BETA=v/c, etc..  For (pseudo)scalar hadrons the decay will
  //   be isotropic in this  sphere so the fractional area (out of 4-pi) is the
  //   fraction of decay that hit the target.  For a given target point and
  //   area, and given x-y components of decay transverse location and slope,
  //   and given decay distance from target ans given decay GAMMA and
  //   rest-frame neutrino energy, the lab energy at the target and the
  //   fractional solid angle in the rest-frame are determined.
  //   For muon decay, correction for non-isotropic nature of decay is done.

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
  double xpos = xyz.X();
  double ypos = xyz.Y();
  double zpos = xyz.Z();

  enu    = 0.0;  // don't know what the final value is
  wgt_xy = 0.0;  // but set these in case we return early due to error


  // in principle we should get these from the particle DB
  // but for consistency testing use the hardcoded values
  double parent_mass = NuMI::MassFromGeantCode(decay->ptype);

  double parentp2 = ( decay->pdPx*decay->pdPx +
                      decay->pdPy*decay->pdPy +
                      decay->pdPz*decay->pdPz );
  double parent_energy = TMath::Sqrt( parentp2 +
                                     parent_mass*parent_mass);
  double parentp = TMath::Sqrt( parentp2 );

  double gamma     = parent_energy / parent_mass;
  double gamma_sqr = gamma * gamma;
  double beta_mag  = TMath::Sqrt( ( gamma_sqr - 1.0 )/gamma_sqr );

  // Get the neutrino energy in the parent decay CM
  double enuzr = decay->Necm;
  // Get angle from parent line of flight to chosen point in beam frame
  double rad = TMath::Sqrt( (xpos-decay->Vx)*(xpos-decay->Vx) +
                           (ypos-decay->Vy)*(ypos-decay->Vy) +
                           (zpos-decay->Vz)*(zpos-decay->Vz) );

  double emrat = 1.0;
  double costh_pardet = -999., theta_pardet = -999.;

  // boost correction, but only if parent hasn't stopped
  if ( parentp > 0. ) {
      costh_pardet = ( decay->pdPx*(xpos-decay->Vx) +
                      decay->pdPy*(ypos-decay->Vy) +
                      decay->pdPz*(zpos-decay->Vz) )
      / ( parentp * rad);
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
  //SAA//                   ( (zpos-decay->Vz)*(zpos-decay->Vz) ) ) / 4.0;
  double sanddetcomp = TMath::Sqrt( ( (xpos-decay->Vx)*(xpos-decay->Vx) ) +
                                   ( (ypos-decay->Vy)*(ypos-decay->Vy) ) +
                                   ( (zpos-decay->Vz)*(zpos-decay->Vz) )   );
  double sangdet = (1.0-TMath::Cos(TMath::ATan( kRDET / sanddetcomp )))/2.0;

  // Weight for solid angle and lorentz boost
  wgt_xy = sangdet * ( emrat * emrat );  // ! the weight ... normally

  // Done for all except polarized muon decay
  // in which case need to modify weight
  // (must be done in double precision)
  if ( decay->ptype == kgeant_muplus || decay->ptype == kgeant_muminus) {
      double beta[3], p_dcm_nu[4], p_nu[3], p_pcm_mp[3], partial;

      // Boost neu neutrino to mu decay CM
      beta[0] = decay->pdPx / parent_energy;
      beta[1] = decay->pdPy / parent_energy;
      beta[2] = decay->pdPz / parent_energy;
      p_nu[0] = (xpos-decay->Vx)*enu/rad;
      p_nu[1] = (ypos-decay->Vy)*enu/rad;
      p_nu[2] = (zpos-decay->Vz)*enu/rad;
      partial = gamma *
      (beta[0]*p_nu[0] + beta[1]*p_nu[1] + beta[2]*p_nu[2] );
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
      double particle_energy = decay->ppenergy;
      gamma = particle_energy/parent_mass;
      beta[0] = decay->ppdxdz * decay->pppz / particle_energy;
      beta[1] = decay->ppdydz * decay->pppz / particle_energy;
      beta[2] =                    decay->pppz / particle_energy;
      partial = gamma * ( beta[0]*decay->muparpx +
                         beta[1]*decay->muparpy +
                         beta[2]*decay->muparpz );
      partial = decay->mupare - partial/(gamma+1.0);
      p_pcm_mp[0] = decay->muparpx - beta[0]*gamma*partial;
      p_pcm_mp[1] = decay->muparpy - beta[1]*gamma*partial;
      p_pcm_mp[2] = decay->muparpz - beta[2]*gamma*partial;
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
      if ( costh >  1.0 ) costh =  1.0;
      if ( costh < -1.0 ) costh = -1.0;
      // Calc relative weight due to angle difference
      double wgt_ratio = 0.0;
      /*
       if (decay->Ntype == kpdg_nuebar) wgt_ratio = 1.0 - costh;
       if (decay->Ntype == kpdg_numubar) {
       double xnu = 2.0 * enuzr / kMUMASS;
       wgt_ratio = ( (3.0-2.0*xnu )  - (1.0-2.0*xnu)*costh ) / (3.0-2.0*xnu);
       }
       */
      switch ( decay->Ntype ) {
          case kgeant_nue:
          case kgeant_nuebar:
              wgt_ratio = 1.0 - costh;
              break;
          case kgeant_numu:
          case kgeant_numubar:
          {
              double xnu = 2.0 * enuzr / kMUMASS;
              wgt_ratio = ( (3.0-2.0*xnu )  - (1.0-2.0*xnu)*costh ) / (3.0-2.0*xnu);
              break;
          }
          default:
              return 2; // bad neutrino type
      }

      wgt_xy = wgt_xy * wgt_ratio;

  } // ptype is muon

  return 0;
}
