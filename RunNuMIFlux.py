import os,sys,string,time
import re
import ROOT

recompile = eval(sys.argv[1])

if recompile:
    ROOT.gROOT.ProcessLine(".L NuMIFlux.cc+")
ROOT.gSystem.Load("NuMIFlux_cc.so")
ROOT.gROOT.SetBatch(True)

from ROOT import NuMIFlux

f = NuMIFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files/rhc/flugg*_700*.root")
f.CalculateFlux()
