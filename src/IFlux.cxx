#define IFlux_cxx

#include <cstddef>
#include <vector>
#include <stdlib.h>

#include "TRandom.h"
#include "TRotation.h"
#include "TMath.h"
#include "TDatime.h"

#include "IFlux.h"

//___________________________________________________________________________
TVector3 IFlux::RandomInTPC() const
{
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
TVector3 IFlux::FromDetToBeam( const TVector3& det ) const
{

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
  if (fDebug) {
    std::cout << "R_{beam to det} = " << std::endl;
    std::cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << std::endl;
    std::cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << std::endl;
    std::cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << std::endl;
    std::cout << std::endl;
  }
  R.Invert(); // R is now the inverse
  if (fDebug) {
    std::cout << "R_{det to beam} = " << std::endl;
    std::cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << std::endl;
    std::cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << std::endl;
    std::cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << std::endl;
    std::cout << std::endl;
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
double IFlux::EstimatePOT(int highest_potnum) const
{
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
