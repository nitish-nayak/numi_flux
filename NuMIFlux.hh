#ifndef NuMIFlux_h
#define NuMIFlux_h

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#include "TChain.h"
#include "TSystem.h"
#include "TVector3.h"

#include "FluxNtuple.h"
//#include "calcLocationWeights.h"

class NuMIFlux {
public :

  int numu = 56;
  int highest_evtno = 0;
  bool debug = false;

  TChain *cflux;

  NuMIFlux(string pattern="/uboone/data/flux/numi/current/flugg_mn000z200i_20101117.gpcfgrid_lowth/flugg_mn000z200i_20101117.gpcfgrid_lowth_001.root");
  virtual ~NuMIFlux();

  void CalculateFlux();
  TVector3 RandomInTPC();
  TVector3 FromDetToBeam(const TVector3& det);
  int calcEnuWgt( FluxNtuple& decay, const TVector3& xyz, double& enu, double& wgt_xy);

};

#endif

#ifdef NuMIFlux_cxx

NuMIFlux::NuMIFlux(string pattern) {

  const char* path = "/uboone/app/users/mdeltutt/NuMIFlux";
  if ( path ) {
    TString libs = gSystem->GetDynamicPath();
    libs += ":";
    libs += path;
    //libs += "/lib";
    gSystem->SetDynamicPath(libs.Data());       
    gSystem->Load("FluxNtuple_C.so");
  }

  cflux = new TChain("h10");
  cflux->Add(pattern.c_str());

  cout << "Number of files: " << cflux->GetNtrees() << endl;


}

NuMIFlux::~NuMIFlux() {

}

#endif // #ifdef NuMIFlux_cxx
