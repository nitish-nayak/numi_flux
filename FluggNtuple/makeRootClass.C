void makeRootClass() {
  TFile *f = new TFile("/uboone/data/flux/numi/v2/flugg_mn000z200i_rp11_lowth_pnut_f112c0f093bbird/flugg_mn000z200i_rp11_bs1.1_pnut_lowth_f112c0f093bbird_0000.root");
  TTree *v = (TTree*)f->Get("h10");
  h10->MakeClass("FluxNtuple");
}
