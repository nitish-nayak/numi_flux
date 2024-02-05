import os,sys,string,time
import re
import ROOT

recompile = False
if len(sys.argv) > 1:
    recompile = eval(sys.argv[1])

if recompile:
    ROOT.gROOT.ProcessLine(".L NuMIFlux.cc+")
ROOT.gSystem.Load("NuMIFlux_cc.so")
ROOT.gROOT.SetBatch(True)

from ROOT import NuMIFlux

#  f = NuMIFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files/rhc/flugg*_70*.root")
#  f.CalculateFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/flugg_rhc.root")
f = NuMIFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files/fhc/flugg*_00*.root")
f.CalculateFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/flugg_fhc.root")
