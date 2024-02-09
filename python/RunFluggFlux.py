import os
import sys
import re
import ROOT

ROOT.gInterpreter.AddIncludePath('../include')
ROOT.gInterpreter.AddIncludePath('../flugg')
ROOT.gSystem.Load("../flugg/FluggFlux_cc.so")
ROOT.gROOT.SetBatch(True)

#  outfile='/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/flugg_rhc.root'
#  f = FluggFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files/rhc/flugg*_70*.root", outfile)
#  f.CalculateFlux()

outfile='/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/flugg_fhc_test.root'
f = ROOT.FluggFlux('/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files/fhc/flugg*_0000.root', outfile)
f.CalculateFlux()
