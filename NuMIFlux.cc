#include <iostream>	
#include <iomanip>
#include <string>
	
using namespace std;
	
#include "TChain.h"
#include "TSystem.h"

//#include "FluxNtuple.h"

//#include "dk2nu.h"
//#include "dkmeta.h"
#include "calcLocationWeights.cxx"

	
#ifndef __CINT__
#endif  // ifndef __CINT__


void NuMIFlux(string pattern="/uboone/data/flux/numi/current/flugg_mn000z200i_20101117.gpcfgrid_lowth/flugg_mn000z200i_20101117.gpcfgrid_lowth_00*.root"){

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


  TChain* cflux = new TChain("h10");
  cflux->Add(pattern.c_str());

  FluxNtuple fluxNtuple(cflux);//  = new FluxNtuple(cflux);
  //cflux->SetBranchAddress("Ntype",&fluxNtuple);


  Long64_t nflux = cflux->GetEntries();
  std::cout << "Total number of entries: " << nflux << std::endl;
  for (Long64_t i=0; i < nflux; ++i ) {
    if (i % 1000 == 0) cout << "On entry " << i << endl;
    cflux->GetEntry(i);
    if ( i < 5 ) cout << "Ntype" << fluxNtuple.Ntype << endl;

    TVector3 xyz;// = new TVector3();
    double wgt_xy = 0;
    double enu = 0;

    // Pick a random point in the TPC
    //xyz = RandomInTPC();

    // From detector to beam coordinates
    TVector3 xyz_beam = FromDetToBeam(xyz);
    cout << xyz_beam.X() << endl;
    
    // Calculate the weight
    int ret = calcEnuWgt(fluxNtuple, xyz_beam, enu, wgt_xy);     
    if (ret != 0) cout << "Error with bsim::calcEnuWgt. Return " << ret << endl;
    cout << "wgt_xy " << wgt_xy << endl;
    break; 
  }

}
