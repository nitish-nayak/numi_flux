#define NuMIFlux_cxx
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

#include "NuMIFlux.hh"
#include "FluggNtuple/FluxNtuple.h"

//#include "dk2nu.h"
//#include "dkmeta.h"
//#include "calcLocationWeights.cxx"
//#include "FluxNtuple_C.so"
	


void NuMIFlux::CalculateFlux() {
//void NuMIFlux(string pattern="/uboone/data/flux/numi/current/flugg_mn000z200i_20101117.gpcfgrid_lowth/flugg_mn000z200i_20101117.gpcfgrid_lowth_001.root"){

/*
  //const char* path = gSystem->ExpandPathName("$(DK2NU)");
  const char* path = "/uboone/app/users/mdeltutt/NuMIFlux";
  if ( path ) {
    TString libs = gSystem->GetDynamicPath();
    libs += ":";
    libs += path;
    //libs += "/lib";
    gSystem->SetDynamicPath(libs.Data());	
    gSystem->Load("FluxNtuple_C.so");
  }
*/

  //TTree * tree = new TTree();

  FluxNtuple fluxNtuple(cflux);//  = new FluxNtuple(cflux);

  //***************************************
  //
  //  Loop over the entries.
  //
  //***************************************

  Long64_t nflux = cflux->GetEntries();
  std::cout << "Total number of entries: " << nflux << std::endl;
  for (Long64_t i=0; i < nflux; ++i ) {
    if (i % 100000 == 0) cout << "On entry " << i << endl;
    cflux->GetEntry(i);

    if (fluxNtuple.Ntype == numu) {
      double wgt_xy = 0.;  // neutrino weight
      double enu    = 0.;  // neutrino energy in lab frame

      // Pick a random point in the TPC (in detector coordinates)
      TVector3 xyz_det = RandomInTPC();
      if (debug) cout << "xyz_det = [" << xyz_det.X() << ", " << xyz_det.Y() << ", " << xyz_det.Z() << "]" << endl;

      // From detector to beam coordinates
      TVector3 xyz_beam = FromDetToBeam(xyz_det);
      if (debug) cout << "xyz_beam = [" << xyz_beam.X() << ", " << xyz_beam.Y() << ", " << xyz_beam.Z() << "]" << endl;   
 
      // Calculate the weight
      int ret = calcEnuWgt(fluxNtuple, xyz_beam, enu, wgt_xy);     
      if (ret != 0) cout << "Error with calcEnuWgt. Return " << ret << endl;
      if (debug) cout << "wgt_xy " << wgt_xy << endl;

      // Fill the histogram
      double weight = wgt_xy * fluxNtuple.Nimpwt * fDefaultWeightCorrection;
      nuFluxHisto->Fill(enu, weight);

      // POT stuff
      if ( fluxNtuple.evtno > highest_evtno ) 
        highest_evtno = fluxNtuple.evtno;

    }
  } // end of loop over the entries

  double AccumulatedPOT = estimate_pots(highest_evtno) * Nfiles;
  double scale = NominalPOT/AccumulatedPOT;
  nuFluxHisto->Scale(scale);

  //***************************************
  //
  // Apply now GENIE xsec
  //
  // root -l  $GENIEXSECPATH/xsec_graphs.root
  // >  _file0->cd("nu_mu_Ar40")
  // >  tot_cc->Draw()
  //
  //***************************************

  const char* genieXsecPath = gSystem->ExpandPathName("$(GENIEXSECPATH)");
  if ( !genieXsecPath ) {
    std::cout << "$(GENIEXSECPATH) not defined." << std::endl;
    std::cout << "Please setup *genie_xsec*. (setup genie_xsec R-2_8_0   -q default)." << std::endl; 
  }
  TString genieXsecFileName = genieXsecPath;
  genieXsecFileName += "/xsec_graphs.root";
  TFile *genieXsecFile = new TFile(genieXsecFileName,"READ");
  genieXsecFile->cd("nu_mu_Ar40");
  TGraph *genieXsecNumuCC = (TGraph *) gDirectory->Get("tot_cc");
  genieXsecFile->Close();

  // TSpline3* genieXsecSplineNumuCC = new TSpline3("genieXsecSplineNumuCC", genieXsecNumuCC, "", 0,6);


  double value;
  for(int i=1; i<histNbins+1; i++) {
    value = nuFluxHisto->GetBinContent(i);
    value *= genieXsecNumuCC->Eval(nuFluxHisto->GetBinCenter(i)); // Eval implies linear interpolation
    value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
    nuCCHisto->SetBinContent(i,value);

  }




  f->cd();
  nuCCHisto->Write();
  genieXsecNumuCC->Write();
  nuFluxHisto->Write();
  f->Close();

  cout << "POT: " << highest_evtno << endl;

}


//___________________________________________________________________________
TVector3 NuMIFlux::RandomInTPC() {
    
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
    
    return det;
}


//___________________________________________________________________________
TVector3 NuMIFlux::FromDetToBeam( const TVector3& det ) {
    
    TVector3 beam;
    TRotation R;
    TVector3 newX(0.921228671,   0.00136256111, -0.389019125);
    TVector3 newY(0.0226872648,  0.998103714,    0.0572211871);
    TVector3 newZ(0.388359401,  -0.061539578,    0.919450845);
    
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
    TVector3 NuMIDet (54.499, 74.461,  677.611); // m
    NuMIDet *= 100.; // To have NuMIDet in cm
    
    beam = R * det + NuMIDet;
    
    return beam;
}



//___________________________________________________________________________
double NuMIFlux::estimate_pots(int highest_potnum) {

  // looks like low counts are due to "evtno" not including                     
  // protons that miss the actual target (hit baffle, etc)
  // this will vary with beam conditions parameters
  // so we should round way up, those generating flugg files
  // aren't going to quantize less than 1000
  // though 500 would probably be okay, 100 would not.
  // also can't use _last_ potnum because muons decays don't
  // have theirs set

  const Int_t    nquant = 1000; // 500;  // 100                                 
  const Double_t rquant = nquant;

  Int_t estimate = (TMath::FloorNint((highest_potnum-1)/rquant)+1)*nquant;
  return estimate;
}



//___________________________________________________________________________
int NuMIFlux::calcEnuWgt( FluxNtuple& decay, const TVector3& xyz,
                         double& enu, double& wgt_xy)
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
    const double kPIMASS = 0.13957;
    const double kKMASS  = 0.49368;
    const double kK0MASS = 0.49767;
    const double kMUMASS = 0.105658389;
    const double kOMEGAMASS = 1.67245;
    
    /*
     const int kpdg_nue       =   12;  // extended Geant 53
     const int kpdg_nuebar    =  -12;  // extended Geant 52
     const int kpdg_numu      =   14;  // extended Geant 56
     const int kpdg_numubar   =  -14;  // extended Geant 55
     
     const int kpdg_muplus     =   -13;  // Geant  5
     const int kpdg_muminus    =    13;  // Geant  6
     const int kpdg_pionplus   =   211;  // Geant  8
     const int kpdg_pionminus  =  -211;  // Geant  9
     const int kpdg_k0long     =   130;  // Geant 10  ( K0=311, K0S=310 )
     const int kpdg_k0short    =   310;  // Geant 16
     const int kpdg_k0mix      =   311;
     const int kpdg_kaonplus   =   321;  // Geant 11
     const int kpdg_kaonminus  =  -321;  // Geant 12
     const int kpdg_omegaminus =  3334;  // Geant 24
     const int kpdg_omegaplus  = -3334;  // Geant 32
     */
    
    // Marco: redefine (hopefully just for now)
    
    const int kpdg_nue       =   53;  // extended Geant 53
    const int kpdg_nuebar    =  52;  // extended Geant 52
    const int kpdg_numu      =   56;  // extended Geant 56
    const int kpdg_numubar   =  55;  // extended Geant 55
    
    const int kpdg_muplus     =   5;  // Geant  5
    const int kpdg_muminus    =    6;  // Geant  6
    const int kpdg_pionplus   =   8;  // Geant  8
    const int kpdg_pionminus  =  9;  // Geant  9
    const int kpdg_k0long     =   10;  // Geant 10  ( K0=311, K0S=310 )
    const int kpdg_k0short    =   16;  // Geant 16
    const int kpdg_k0mix      =   311;
    const int kpdg_kaonplus   =   11;  // Geant 11
    const int kpdg_kaonminus  =  12;  // Geant 12
    const int kpdg_omegaminus =  24;  // Geant 24
    const int kpdg_omegaplus  = 32;  // Geant 32
    
    // Marco: end
    
    
    
    const double kRDET = 100.0;   // set to flux per 100 cm radius
    
    double xpos = xyz.X();
    double ypos = xyz.Y();
    double zpos = xyz.Z();
    
    enu    = 0.0;  // don't know what the final value is
    wgt_xy = 0.0;  // but set these in case we return early due to error
    
    
    // in principle we should get these from the particle DB
    // but for consistency testing use the hardcoded values
    double parent_mass = kPIMASS;
    
    /*
     if ( decay.ptype == kpdg_pionminus)  parent_mass = kPIMASS;
     if ( decay.ptype == kpdg_kaonminus)  parent_mass = kKMASS;
     if ( decay.ptype == kpdg_k0mix)      parent_mass = kK0MASS;
     if ( decay.ptype == kpdg_muminus)    parent_mass = kMUMASS;
     if ( decay.ptype == kpdg_omegaplus)  parent_mass = kOMEGAMASS;
     */
    switch ( decay.ptype ) {
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
        default:
            std::cerr << "bsim::calcEnuWgt unknown particle type " << decay.ptype
            << std::endl << std::flush;
            assert(0);
            return 1;
    }
    
    
    
    
    double parentp2 = ( decay.pdPx*decay.pdPx +
                       decay.pdPy*decay.pdPy +
                       decay.pdPz*decay.pdPz );
    double parent_energy = TMath::Sqrt( parentp2 +
                                       parent_mass*parent_mass);
    double parentp = TMath::Sqrt( parentp2 );
    
    double gamma     = parent_energy / parent_mass;
    double gamma_sqr = gamma * gamma;
    double beta_mag  = TMath::Sqrt( ( gamma_sqr - 1.0 )/gamma_sqr );
    
    // Get the neutrino energy in the parent decay CM
    double enuzr = decay.Necm;
    // Get angle from parent line of flight to chosen point in beam frame
    double rad = TMath::Sqrt( (xpos-decay.Vx)*(xpos-decay.Vx) +
                             (ypos-decay.Vy)*(ypos-decay.Vy) +
                             (zpos-decay.Vz)*(zpos-decay.Vz) );
    
    double emrat = 1.0;
    double costh_pardet = -999., theta_pardet = -999.;
    
    // boost correction, but only if parent hasn't stopped
    if ( parentp > 0. ) {
        costh_pardet = ( decay.pdPx*(xpos-decay.Vx) +
                        decay.pdPy*(ypos-decay.Vy) +
                        decay.pdPz*(zpos-decay.Vz) )
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
    //SAA//                   ( (zpos-decay.Vz)*(zpos-decay.Vz) ) ) / 4.0;
    double sanddetcomp = TMath::Sqrt( ( (xpos-decay.Vx)*(xpos-decay.Vx) ) +
                                     ( (ypos-decay.Vy)*(ypos-decay.Vy) ) +
                                     ( (zpos-decay.Vz)*(zpos-decay.Vz) )   );
    double sangdet = (1.0-TMath::Cos(TMath::ATan( kRDET / sanddetcomp )))/2.0;
    
    // Weight for solid angle and lorentz boost
    wgt_xy = sangdet * ( emrat * emrat );  // ! the weight ... normally
    
    // Done for all except polarized muon decay
    // in which case need to modify weight
    // (must be done in double precision)
    if ( decay.ptype == kpdg_muplus || decay.ptype == kpdg_muminus) {
        double beta[3], p_dcm_nu[4], p_nu[3], p_pcm_mp[3], partial;
        
        // Boost neu neutrino to mu decay CM
        beta[0] = decay.pdPx / parent_energy;
        beta[1] = decay.pdPy / parent_energy;
        beta[2] = decay.pdPz / parent_energy;
        p_nu[0] = (xpos-decay.Vx)*enu/rad;
        p_nu[1] = (ypos-decay.Vy)*enu/rad;
        p_nu[2] = (zpos-decay.Vz)*enu/rad;
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
        double particle_energy = decay.ppenergy;
        gamma = particle_energy/parent_mass;
        beta[0] = decay.ppdxdz * decay.pppz / particle_energy;
        beta[1] = decay.ppdydz * decay.pppz / particle_energy;
        beta[2] =                    decay.pppz / particle_energy;
        partial = gamma * ( beta[0]*decay.muparpx +
                           beta[1]*decay.muparpy +
                           beta[2]*decay.muparpz );
        partial = decay.mupare - partial/(gamma+1.0);
        p_pcm_mp[0] = decay.muparpx - beta[0]*gamma*partial;
        p_pcm_mp[1] = decay.muparpy - beta[1]*gamma*partial;
        p_pcm_mp[2] = decay.muparpz - beta[2]*gamma*partial;
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
         if (decay.Ntype == kpdg_nuebar) wgt_ratio = 1.0 - costh;
         if (decay.Ntype == kpdg_numubar) {
         double xnu = 2.0 * enuzr / kMUMASS;
         wgt_ratio = ( (3.0-2.0*xnu )  - (1.0-2.0*xnu)*costh ) / (3.0-2.0*xnu);
         }
         */  
        switch ( decay.Ntype ) {
            case kpdg_nue:
            case kpdg_nuebar:
                wgt_ratio = 1.0 - costh;
                break;
            case kpdg_numu:
            case kpdg_numubar:
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
