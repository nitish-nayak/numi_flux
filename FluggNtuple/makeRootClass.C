void makeRootClass() {
  TFile *f = new TFile("/uboone/data/flux/numi/current/flugg_mn000z200i_20101117.gpcfgrid_lowth/flugg_mn000z200i_20101117.gpcfgrid_lowth_001.root");
  TTree *v = (TTree*)f->Get("h10");
  h10->MakeClass("FluxNtuple");
}
