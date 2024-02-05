void makeRootClass() {
  TFile *f = new TFile("/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files/rhc/flugg_mn000z-200i_rp11_bs1.1_pnut_lowth_f11f093bbird_target_7000.root");
  TTree *v = (TTree*)f->Get("h10");
  v->MakeClass("FluxNtuple");
}
