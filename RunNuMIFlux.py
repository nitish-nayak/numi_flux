#!/usr/bin/env python

import os,sys,string, time
import ROOT

ROOT.gROOT.ProcessLine(".L NuMIFlux.cc+")
ROOT.gSystem.Load("NuMIFlux_cc.so")
ROOT.gROOT.SetBatch(True)

from ROOT import NuMIFlux
from glob import glob

#  fname = glob("/uboone/data/flux/numi/current/flugg_mn000z200i_20101117.gpcfgrid_lowth/flugg_mn000z200i_20101117.gpcfgrid_lowth_00*.root")

f = NuMIFlux()
f.CalculateFlux()
