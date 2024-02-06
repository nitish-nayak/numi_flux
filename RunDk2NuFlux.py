import os,sys,string,time
import re
import ROOT
import subprocess

recompile = False
if len(sys.argv) > 1:
    recompile = eval(sys.argv[1])

if recompile:
    #  ROOT.gSystem.AddIncludePath(" -I${PPFX_DIR}/include");
    #  ROOT.gSystem.AddIncludePath(" -I${BOOSTROOT}");
    #  ROOT.gSystem.AddIncludePath(" -I${DK2NU}/include");
    #
    #  print(ROOT.gSystem.GetMakeSharedLib())
    #  ROOT.gSystem.Load("${PPFX_DIR}/lib/libppfx.so")
    #
    #  ROOT.gROOT.ProcessLine(".L Dk2NuFlux.cc+")
    subprocess.call(["make", "dk2nu"])

ROOT.gROOT.ProcessLine('#include "Dk2NuFlux.hh"')
ROOT.gSystem.Load("Dk2NuFlux_cc.so")
ROOT.gROOT.SetBatch(True)

#  f = Dk2NuFlux("/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_dk2nu_zero_threshold/RHC/*_5*_0000.root")
#  f.CalculateFlux("/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_rhc.root")
outfile="/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_fhc_ppfx.root"
f = ROOT.Dk2NuFlux("/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_dk2nu_zero_threshold/FHC/g4numiv6_minervame_me000z200i_2*.root", outfile)
f.CalculateFlux()
