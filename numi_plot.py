#!/usr/bin/env python
#
# read a root anatree using Pyroot
#
# run with python pyrootMacroTPCActive.py 14 -1
# 14 is the nu flavour (nu-mu in this case)
# -1 means loop over all the entries

import os,sys,string, time
import ROOT
from math import *
from ROOT import TTree, TObject, TFile, gDirectory, TH1D, TH2D, TH3D, TCanvas, gROOT, TGaxis, gStyle, TColor, TLegend, THStack, TChain, TLatex, TText
#from ROOT import *
from array import array
from glob import glob

# Importing rootlogon.C
ROOT.gROOT.SetMacroPath('~/');
ROOT.gROOT.Macro( os.path.expanduser( 'rootlogon.C' ) )

# Opening root file
#f = TFile("files/NuMIFlux.root")
f = TFile("files/NuMIFlux_anti.root")
f.ls()

h_numu = f.Get("numuFluxHisto")
h_anumu = f.Get("anumuFluxHisto")
h_nue = f.Get("nueFluxHisto")
h_anue = f.Get("anueFluxHisto")

h_numu.SetLineColor(ROOT.kRed+1)
h_numu.GetXaxis().SetRangeUser(0,15)
#h_numu.GetXaxis().SetRangeUser(0,6)
h_numu.SetTitle("")
h_numu.GetXaxis().SetTitle("Neutrino Energy [GeV]")
h_numu.GetYaxis().SetTitle("#Phi(#nu) / 50 MeV / cm^{2} / 6x10^{20} POT")
h_anumu.SetLineColor(ROOT.kBlue+1)
h_nue.SetLineColor(ROOT.kRed+1)
h_nue.SetLineStyle(2)
h_anue.SetLineColor(ROOT.kBlue+1)
h_anue.SetLineStyle(2)




c = TCanvas()
c.SetLogy()


h_numu.Draw()
h_anumu.Draw("same")
h_nue.Draw("same")
h_anue.Draw("same")


leg = TLegend(.5, .55, .9, .85)
leg.SetFillStyle(0);
leg.AddEntry(h_numu,   "#nu_{#mu}",        "l");
leg.AddEntry(h_anumu,  "#bar{#nu}_{#mu}",       "l");
leg.AddEntry(h_nue,   "#nu_{e}",         "l");
leg.AddEntry(h_anue,  "#bar{#nu}_{e}",   "l");
leg.Draw();


#t = TLatex(.9, .95, "MicroBooNE Simulation");
#t.SetTextColor(ROOT.kGray+1);
#t.SetNDC();
#t.SetTextSize(2/30.);
#t.SetTextAlign(32);
#t.Draw();

t2 = TLatex(.51, .48, "#splitline{Off-axis NuMI Flux}{at MicroBooNE}");
t2.SetTextColor(ROOT.kRed+2);
t2.SetNDC();
t2.SetTextSize(1.4/30.);
t2.SetTextAlign(11);
t2.Draw();

t3 = TLatex(.51, .40, "Anti-Neutrino Mode");
t3.SetTextColor(ROOT.kBlack);
t3.SetNDC();
t3.SetTextSize(1.4/30.);
t3.SetTextAlign(11);
t3.Draw();






raw_input("Please press enter to exit.")
