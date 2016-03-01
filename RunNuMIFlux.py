#!/usr/bin/env python
#

import os,sys,string, time
import ROOT
ROOT.gSystem.Load("FluggNtuple/FluxNtuple_C.so")
ROOT.gSystem.Load("NuMIFlux_cc.so")
from ROOT import NuMIFlux
from glob import glob

fname = glob("/uboone/data/flux/numi/current/flugg_mn000z200i_20101117.gpcfgrid_lowth/flugg_mn000z200i_20101117.gpcfgrid_lowth_00*.root")

f = NuMIFlux()
f.CalculateFlux()


raw_input("Please press enter to exit.")
