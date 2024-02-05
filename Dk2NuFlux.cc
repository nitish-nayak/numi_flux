#define Dk2NuFlux_cxx
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#include "TChain.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TRotation.h"
#include "TMath.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include <TParticlePDG.h>
#include <TDatabasePDG.h>

#include "Dk2NuFlux.hh"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/calcLocationWeights.h"

Dk2NuFlux::Dk2NuFlux(string pattern) {

  const char* path = "/uboone/app/users/bnayak/flugg_reweight/flugg_pointing/NuMIFlux/FluggNtuple";
  if ( path ) {
    TString libs = gSystem->GetDynamicPath();
    libs += ":";
    libs += path;
    gSystem->SetDynamicPath(libs.Data());
    // gSystem->Load("FluxNtuple_C.so");
  }

  cflux = new TChain("dk2nuTree");
  cflux->Add(pattern.c_str());

  cflux_meta = new TChain("dkmetaTree");
  cflux_meta->Add(pattern.c_str());

  Nfiles = cflux->GetNtrees();
  cout << "Number of files: " << Nfiles << endl;

  //Inizialise histos
  TString titleBase1 = "Neutrino Flux;";
  TString titleBase2 = " Energy [GeV];";
  TString titleBase3 = " / cm^{2} / 6e20 POT";
  // numu
  numuFluxHisto = new TH1D("numuFluxHisto", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
  numuCCHisto = new TH1D("numuCCHisto", "numu CC; #nu_{#mu} Energy [GeV]; #nu_{#mu} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  // anumu
  anumuFluxHisto = new TH1D("anumuFluxHisto", (titleBase1 + "#bar{#nu}_{#mu}" + titleBase2 +"#bar{#nu}_{#mu}" + titleBase3),histNbins,histMin,histMax);
  anumuCCHisto = new TH1D("anumuCCHisto", "numu bar CC; #bar{#nu}_{#mu} Energy [GeV]; #bar{#nu}_{#mu} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  // nue
  nueFluxHisto = new TH1D("nueFluxHisto", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
  nueCCHisto = new TH1D("nueCCHisto", "nue CC; #nu_{e} Energy [GeV]; #nu_{e} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  // anue
  anueFluxHisto = new TH1D("anueFluxHisto", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 + "#bar{#nu}_{e}" + titleBase3),histNbins,histMin,histMax);
  anueCCHisto = new TH1D("anueCCHisto", "nue bar CC; #bar{#nu}_{e} Energy [GeV]; #bar{#nu}_{#mu} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  hPOT = new TH1D("POT", "Total POT", 1, 0, 1);

  outTree = new TTree("outTree", "outTree");
  outTree->Branch("nuE", &nuE, "nuE/F");
  outTree->Branch("wgt", &wgt, "wgt/F");
  outTree->Branch("ptype", &ptype, "ptype/I");
  outTree->Branch("ntype", &ntype, "ntype/I");
  outTree->Branch("ncascade", &ncascade, "ncascade/I");
  outTree->Branch("pmedium", &pmedium, "pmedium/I");
  outTree->Branch("decaytype", &decaytype, "decaytype/I");

}

Dk2NuFlux::~Dk2NuFlux() {

}

void Dk2NuFlux::CalculateFlux(string outfile) {

  bsim::Dk2Nu* fDk2Nu = new bsim::Dk2Nu;
  bsim::DkMeta* fDkMeta= new bsim::DkMeta;
  cflux->SetBranchAddress("dk2nu", &fDk2Nu);
  cflux_meta->SetBranchAddress("dkmeta",&fDkMeta);

  //***************************************
  //
  //  Loop over the entries.
  //
  //***************************************

  Long64_t nflux = cflux->GetEntries();
  std::cout << "Total number of entries: " << nflux << std::endl;
  for (Long64_t i=0; i < nflux; ++i ) {

    // Get entry i. fluxNtuple is now filled with entry i info.
    cflux->GetEntry(i);
    cflux_meta->GetEntry(i);

    // Alert the user
    if (i % 100000 == 0) cout << "On entry " << i/1.0e6 << " M"<< endl;

    if(treeNumber != cflux->GetTreeNumber()) {
      treeNumber = cflux->GetTreeNumber();
      std::cout << "Moving to tree number " << treeNumber << "." << std::endl;
      AccumulatedPOT += fDkMeta->pots;
      cout << "AccumulatedPOT: " << AccumulatedPOT << endl;
    }

    double wgt_xy = 0.;  // neutrino weight
    double enu    = 0.;  // neutrino energy in lab frame

    // Pick a random point in the TPC (in detector coordinates)
    TVector3 xyz_det = RandomInTPC();
    if (debug) cout << "xyz_det = [" << xyz_det.X() << ", " << xyz_det.Y() << ", " << xyz_det.Z() << "]" << endl;

    // From detector to beam coordinates
    TVector3 xyz_beam = FromDetToBeam(xyz_det);
    if (debug) cout << "xyz_beam = [" << xyz_beam.X() << ", " << xyz_beam.Y() << ", " << xyz_beam.Z() << "]" << endl;

    // Neutrons are not yet implemented so skip for now
    // if (fDk2Nu->decay.ptype == 2112) continue;

    // Calculate the weight
    int ret = calcEnuWgt(fDk2Nu, xyz_beam, enu, wgt_xy);
    if (ret != 0) cout << "Error with calcEnuWgt. Return " << ret << endl;
    if (debug) cout << "wgt_xy " << wgt_xy << endl;


    // Calculate the total weight
    double weight = wgt_xy * fDk2Nu->decay.nimpwt * fDefaultWeightCorrection;

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
      case numu:
        numuFluxHisto->Fill(enu, weight);
        break;
      case anumu:
        anumuFluxHisto->Fill(enu, weight);
        break;
      case nue:
        nueFluxHisto->Fill(enu, weight);
        break;
      case anue:
        anueFluxHisto->Fill(enu, weight);
        break;
    }

    // // Now fill the contraint histograms for each file
    // GetConstraints_ThinTarg( fDk2Nu );
    // GetConstraints_ThickTarg(fDk2Nu);
    nuE = enu;
    wgt = weight;
    ptype = fDk2Nu->decay.ptype;
    ntype = fDk2Nu->decay.ntype;
    ncascade = fDk2Nu->tgtexit.tgen;
    pmedium = fDk2Nu->decay.ppmedium;
    decaytype = fDk2Nu->decay.ndecay;

    outTree->Fill();


  } // end of loop over the entries


  //***************************************
  //
  // POT scaling
  //
  //***************************************

  double scale =    NominalPOT/AccumulatedPOT;
  hPOT->SetBinContent(1, AccumulatedPOT);
  numuFluxHisto  -> Scale(scale);
  anumuFluxHisto -> Scale(scale);
  nueFluxHisto   -> Scale(scale);
  anueFluxHisto  -> Scale(scale);

  // // Constrained Histograms
  // pionplus_NA49  -> Scale(scale);
  // pionplus_MIPP  -> Scale(scale);
  // pionminus_NA49 -> Scale(scale);
  // pionminus_MIPP -> Scale(scale);
  // Kplus_NA49     -> Scale(scale);
  // Kplus_MIPP     -> Scale(scale);
  // Kminus_NA49    -> Scale(scale);
  // Kminus_MIPP    -> Scale(scale);
  cout << endl << ">>> TOTAL POT: " << AccumulatedPOT << endl << endl;


  //***************************************
  //
  // Apply now GENIE xsec
  //
  // source /nusoft/app/externals/setup
  // setup genie_xsec R-2_8_0   -q default
  // root -l  $GENIEXSECPATH/xsec_graphs.root
  // >  _file0->cd("nu_mu_Ar40")
  // >  tot_cc->Draw()
  //
  //***************************************

  // const char* genieXsecPath = gSystem->ExpandPathName("$(GENIEXSECPATH)");
  // if ( !genieXsecPath ) {
  //     std::cout << "$(GENIEXSECPATH) not defined." << std::endl;
  //     std::cout << "Please setup *genie_xsec*. (setup genie_xsec R-2_8_0   -q default)." << std::endl;
  // }
  //
  // if ( genieXsecPath ) {
  //     TString genieXsecFileName = genieXsecPath;
  //     genieXsecFileName += "/xsec_graphs.root";
  //     TFile *genieXsecFile = new TFile(genieXsecFileName,"READ");
  //     genieXsecNumuCC    = (TGraph *) genieXsecFile->Get("nu_mu_Ar40/tot_cc");
  //     genieXsecNumubarCC = (TGraph *) genieXsecFile->Get("nu_mu_bar_Ar40/tot_cc");
  //     genieXsecNueCC     = (TGraph *) genieXsecFile->Get("nu_e_Ar40/tot_cc");
  //     genieXsecNuebarCC  = (TGraph *) genieXsecFile->Get("nu_e_bar_Ar40/tot_cc");
  //     genieXsecFile->Close();
  //
  //     // TSpline3* genieXsecSplineNumuCC = new TSpline3("genieXsecSplineNumuCC", genieXsecNumuCC, "", 0,6);
  //
  //     double value;
  //     for(int i=1; i<histNbins+1; i++) {
  //         value = numuFluxHisto->GetBinContent(i);
  //         value *= genieXsecNumuCC->Eval(numuFluxHisto->GetBinCenter(i)); // Eval implies linear interpolation
  //         value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
  //         numuCCHisto->SetBinContent(i, value);
  //
  //         value = anumuFluxHisto->GetBinContent(i);
  //         value *= genieXsecNumubarCC->Eval(anumuFluxHisto->GetBinCenter(i)); // Eval implies linear interpolation
  //         value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
  //         anumuCCHisto->SetBinContent(i, value);
  //
  //         value = nueFluxHisto->GetBinContent(i);
  //         value *= genieXsecNueCC->Eval(nueFluxHisto->GetBinCenter(i)); // Eval implies linear interpolation
  //         value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
  //         nueCCHisto->SetBinContent(i, value);
  //
  //         value = anueFluxHisto->GetBinContent(i);
  //         value *= genieXsecNuebarCC->Eval(anueFluxHisto->GetBinCenter(i)); // Eval implies linear interpolation
  //         value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
  //         anueCCHisto->SetBinContent(i, value);
  //
  //     }
  // } // end if ( genieXsecPath )




  //***************************************
  //
  // Writing on file
  //
  //***************************************

  TFile* f = new TFile(outfile.c_str(), "RECREATE");
  f->cd();
  numuFluxHisto  -> Write();
  anumuFluxHisto -> Write();
  nueFluxHisto   -> Write();
  anueFluxHisto  -> Write();

  hPOT->Write();
  outTree->Write();
  // if ( genieXsecPath ) {
  //     numuCCHisto     -> Write();
  //     genieXsecNumuCC ->SetName("numu_tot_cc");
  //     genieXsecNumuCC -> Write();
  //
  //     anumuCCHisto       -> Write();
  //     genieXsecNumubarCC ->SetName("anumu_tot_cc");
  //     genieXsecNumubarCC -> Write();
  //
  //     nueCCHisto     -> Write();
  //     genieXsecNueCC ->SetName("nue_tot_cc");
  //     genieXsecNueCC -> Write();
  //
  //     anueCCHisto       -> Write();
  //     genieXsecNuebarCC ->SetName("anue_tot_cc");
  //     genieXsecNuebarCC -> Write();
  // }
  f->Close();

}


//___________________________________________________________________________
TVector3 Dk2NuFlux::RandomInTPC() {

  TDatime *d = new TDatime;
  TRandom *r = new TRandom(d->GetTime());

  double xTPC = 256.35;  // cm
  double yTPC = 233.;  // cm
  double zTPC = 1036.8; // cm

  double x = r->Uniform(0., xTPC);
  double y = r->Uniform(-yTPC/2., yTPC/2.);
  double z = r->Uniform(0., zTPC);

  TVector3 det;
  det.SetXYZ(x,y,z);

  delete d;
  delete r;

  return det;
}
//___________________________________________________________________________
TVector3 Dk2NuFlux::FromDetToBeam( const TVector3& det ) {

  TVector3 beam;
  TRotation R;

  //corrected rotation matrix using the 0,0,0 position for MicroBooNE
  //Previous matrix is calculated relative to MiniBooNE, which is not in the centre of the BNB!

  TVector3 newX(0.92103853804025682, 0.0000462540012621546684, -0.38947144863934974);
  TVector3 newY(0.0227135048039241207, 0.99829162468141475, 0.0538324139386641073);
  TVector3 newZ(0.38880857519374290, -0.0584279894529063024, 0.91946400794392302);
  //old matrix
  /*
  TVector3 newX(0.921228671,   0.00136256111, -0.389019125);
  TVector3 newY(0.0226872648,  0.998103714,    0.0572211871);
  TVector3 newZ(0.388359401,  -0.061539578,    0.919450845);
  */

  R.RotateAxes(newX,newY,newZ);
  if (debug) {
    cout << "R_{beam to det} = " << endl;
    cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << endl;
    cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << endl;
    cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << endl;
    cout << endl;
  }
  R.Invert(); // R is now the inverse
  if (debug) {
    cout << "R_{det to beam} = " << endl;
    cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << endl;
    cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << endl;
    cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << endl;
    cout << endl;
  }
  // Now R allows to go from detector to beam coordinates.
  // NuMIDet is vector from NuMI target to uB detector (in beam coordinates)
  // Updated position - leaving old positions here (July 2018)
  //TVector3 NuMIDet (54.499, 74.461,  677.611); // m
  TVector3 NuMIDet (55.02, 72.59,  672.70); //m
  NuMIDet *= 100.; // To have NuMIDet in cm

  beam = R * det + NuMIDet;

  return beam;
}
//___________________________________________________________________________
// DEPRECIATED IN THIS CODE
double Dk2NuFlux::estimate_pots(int highest_potnum) {

  // Stolen: https://cdcvs.fnal.gov/redmine/projects/dk2nu/repository/show/trunk/dk2nu
  // looks like low counts are due to "evtno" not including
  // protons that miss the actual target (hit baffle, etc)
  // this will vary with beam conditions parameters
  // so we should round way up, those generating flugg files
  // aren't going to quantize less than 1000
  // though 500 would probably be okay, 100 would not.
  // also can't use _last_ potnum because muons decay->> don't
  // have theirs set

  // Marco: Trying with 10000
  const Int_t    nquant = 10000; //1000; // 500;  // 100
  const Double_t rquant = nquant;

  Int_t estimate = (TMath::FloorNint((highest_potnum-1)/rquant)+1)*nquant;
  return estimate;
}
//___________________________________________________________________________
int Dk2NuFlux::calcEnuWgt(bsim::Dk2Nu* decay, const TVector3& xyz,double& enu, double& wgt_xy){

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
  const double kRDET = 100.0;   // set to flux per 100 cm radius
  double xpos = xyz.X();
  double ypos = xyz.Y();
  double zpos = xyz.Z();
  enu    = 0.0;  // don't know what the final value is
  wgt_xy = 0.0;  // but set these in case we return early due to error
  // in principle we should get these from the particle DB
  // but for consistency testing use the hardcoded values
  double parent_mass = kPIMASS;
  switch ( decay->decay.ptype ) {

    case kpdg_pionplus:
    case kpdg_pionminus:
        parent_mass = kPIMASS;
        break;

    case kpdg_kaonplus:
    case kpdg_kaonminus:
        parent_mass = kKMASS;
        break;

    case kpdg_k0long:
    case kpdg_k0short:
    case kpdg_k0mix:
            parent_mass = kK0MASS;
            break;

    case kpdg_muplus:
    case kpdg_muminus:
        parent_mass = kMUMASS;
        break;

    case kpdg_omegaminus:
    case kpdg_omegaplus:
        parent_mass = kOMEGAMASS;
        break;

    case kpdg_neutron:
    case kpdg_antineutron:
        parent_mass = kNEUTRONMASS;
        break;

    default:
    std::cerr << "bsim::calcEnuWgt unknown particle type " << decay->decay.ptype
                << std::endl << std::flush;
    enu    = 0.0;
    wgt_xy = 0.0;
    return 1;
  }
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
          std::cerr << "bsim::calcEnuWgt encountered serious problem: "
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
