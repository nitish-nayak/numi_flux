import os,sys,string,time
import re
import ROOT

recompile = False
if len(sys.argv) > 1:
    recompile = eval(sys.argv[1])

if recompile:
    ROOT.gROOT.ProcessLine(".L Dk2NuFlux.cc+")
ROOT.gSystem.Load("Dk2NuFlux_cc.so")
ROOT.gROOT.SetBatch(True)

from ROOT import Dk2NuFlux

#  f = Dk2NuFlux("/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_dk2nu_zero_threshold/RHC/*_5*_0000.root")
#  f.CalculateFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_rhc.root")
f = Dk2NuFlux("/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_dk2nu_zero_threshold/FHC/g4numiv6_minervame_me000z200i_2*root")
f.CalculateFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_fhc.root")
