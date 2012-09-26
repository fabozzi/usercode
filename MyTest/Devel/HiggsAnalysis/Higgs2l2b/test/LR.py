import sys
sys.argv.append('-b')
import ROOT
import os, commands
from settings_2012 import *
from SamplesColors import *
from HiggsAnalysis.Higgs2l2b.scaleFactors import scale


ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()        # don't pop up canvases
c1 = ROOT.TCanvas()

masses=["400","600"]

for mh in masses:
    f = ROOT.TFile("tree_NoNorm_Mu/GluGluToHToZZTo2L2Q_M-"+mh+"_8TeV_tree.root")
    t = f.Get("lljjmassTree")
    ROOT.SetOwnership(t, False)
    
    t.Draw("mZZ")
    hmass = t.GetHistogram()
    hmass.SetLineColor(ROOT.kBlue +2)
    hmass.SetLineWidth(2)
    hmass.SetTitle("m_{H} = "+mh+" GeV")
    
    t.Draw("mZZ","LRweight", "SAME")
    hmassRew = t.GetHistogram()
    hmassRew.SetLineColor(ROOT.kOrange +10)
    hmassRew.SetLineWidth(2)
    hmassRew.GetXaxis().SetRangeUser(0,1000)
    
    leg = ROOT.TLegend(.6, .6, .85, .75)
    leg.SetFillColor(0)
    leg.SetTextSize(0.03)
    leg.AddEntry(hmass, "WO LR", "l")
    leg.AddEntry(hmassRew, "WITH LR", "l")
    leg.Draw("SAME")

    c1.Print("LR_"+mh+".png")
    c1.Print("LR_"+mh+".eps")
    
    
