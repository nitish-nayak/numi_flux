import os
import sys
import re
import ROOT
import subprocess

recompile = False
if len(sys.argv) > 1:
    recompile = eval(sys.argv[1])

if recompile:
    subprocess.call(["make", "dk2nu"])

ROOT.gROOT.ProcessLine('#include "Dk2NuFlux.hh"')
ROOT.gSystem.Load("Dk2NuFlux_cc.so")
ROOT.gROOT.SetBatch(True)

#  outfile="/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_fhc_ppfx.root"
#  f = ROOT.Dk2NuFlux("/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_dk2nu_zero_threshold/FHC/g4numiv6_minervame_me000z200i_2*.root", outfile)

outfile="/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_rhc_ppfx.root"
f = ROOT.Dk2NuFlux("/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_dk2nu_zero_threshold/RHC/*_5*_0000.root", outfile)
f.CalculateFlux()
