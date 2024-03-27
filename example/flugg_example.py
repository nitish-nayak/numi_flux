import os
import sys
import re
import ROOT

ROOT.gSystem.Load("../lib/FluggFlux_cc.so")
ROOT.gROOT.SetBatch(True)

#  outfile='/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/flugg_rhc_extrainfo.root'
#  f = ROOT.FluggFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files/rhc/flugg*_70*.root", outfile)
#  f.CalculateFlux()

outfile='/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/flugg_fhc_extrainfo.root'
f = ROOT.FluggFlux('/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files/fhc/flugg*_00*.root', outfile)
f.CalculateFlux()
