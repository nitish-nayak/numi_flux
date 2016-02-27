#include <strstream>

void setupNuMIFlux() {

  TString path = "";
  path += " -I";
  path += gSystem->ExpandPathName("$(DK2NU)");
  path += "/tree";
  gSystem->AddIncludePath(path);
  cout << gSystem->GetIncludePath() << endl;
  //gROOT->SetMacroPath("/nusoft/app/externals/dk2nu/v01_01_03/Linux64bit+2.6-2.5/debug-e5/dk2nu/tree:.");


  TString t = ".include ";
  t += "/nusoft/app/externals/dk2nu/v01_01_03/Linux64bit+2.6-2.5/debug-e5/dk2nu/tree";
  gROOT->ProcessLine(t);

/*
  ostrstream s;
  s << ".include " << "/nusoft/app/externals/dk2nu/v01_01_03/Linux64bit+2.6-2.5/debug-e5/dk2nu/tree" << ends;
  //s << ".include " << gSystem->Getenv("MYPATH") << ends;
  gROOT->ProcessLine(s.str());
*/
}
