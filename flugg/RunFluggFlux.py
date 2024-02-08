import os
import sys
import re
import ROOT
import subprocess

recompile = False
if len(sys.argv) > 1:
    recompile = eval(sys.argv[1])

if recompile:
    subprocess.call(["make", "flugg"])

ROOT.gROOT.ProcessLine('#include "FluggFlux.hh"')
ROOT.gSystem.Load("FluggFlux_cc.so")
ROOT.gROOT.SetBatch(True)

#  f = FluggFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files/rhc/flugg*_70*.root")
#  f.CalculateFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/flugg_rhc.root")

f = FluggFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files/fhc/flugg*_00*.root")
f.CalculateFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/flugg_fhc.root")
