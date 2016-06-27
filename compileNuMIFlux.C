void compileNuMIFlux(){

  gSystem->Load("FluggNtuple/FluxNtuple_C.so");
  TString t = ".L NuMIFlux.cc+;";
  gROOT->ProcessLine(t);

}
