import os
import sys
import ROOT

ROOT.gSystem.Load("$NUMIANA_DIR/lib/Dk2NuFlux_cc.so")
ROOT.gROOT.SetBatch(True)

filelist = sys.argv[1]
seed = int(sys.argv[2])

outfile="numi_flux_output.root"
f = ROOT.Dk2NuFlux(True, filelist, outfile)
f.SetSeedPPFX(seed)
f.SetModePPFX("ubnumi_multisim")
f.CalculateFlux()
