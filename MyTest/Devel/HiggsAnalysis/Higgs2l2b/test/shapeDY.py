
#!/bin/env python

import sys

sys.argv.append('-b')

import ROOT
from settings import *
from SamplesColors import *
#ROOT.gSystem.Load("libRooFitCore")

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gStyle.SetOptStat(0)

c1 = ROOT.TCanvas("name","title", 5)
c1.SetTicks()


    
var = "zllmass"
xmin = 70
xmax = 110

name2012 = "DoubleMuRun2012A"
name2011 = "DoubleMuRun2011AB"


path11 = "/data3/home/decosa/h2l2j_16Apr12/CMSSW_4_2_8/src/HiggsAnalysis/Higgs2l2b/test/NoNorm_Mu/MuRun2011AllEdmNtp.root"
path12 = "/data3/home/decosa/h2l2j_CMSSW_5_2_3_patch2/CMSSW_5_2_3_patch2/src/HiggsAnalysis/Higgs2l2b/test/NoNorm_Mu/MuRun2012AEdmNtp.root"

pad1 = ROOT.TPad()
ROOT.SetOwnership(pad1, False)
pad1.SetPad('pad1','pad1',0.0,0.4,1.0,1.0, 10)
pad1.Draw('same')
pad1.cd()
pad1.SetTopMargin(0.11)
pad1.SetLeftMargin(0.1)
pad1.SetRightMargin(0.05)
pad1.SetBottomMargin(0.0)
pad1.SetTicks()
c1.cd()
    
pad2 = ROOT.TPad()
ROOT.SetOwnership(pad2, False)
pad2.SetPad('pad2','pad2',0.0,0.,1.0,0.4, 10)
pad2.Draw()
pad2.cd()
pad2.SetTopMargin(0.)
pad2.SetLeftMargin(0.1)
pad2.SetRightMargin(0.05)
pad2.SetBottomMargin(0.20)
pad2.SetTicks()

c1.cd()
pad1.cd()

data2011 = ROOT.TFile.Open(path11).Get(var)
data2012 = ROOT.TFile.Open(path12).Get(var)

a = data2011.Integral()
b = data2012.Integral()

data2011.Scale(b/a)

data2011.SetFillStyle(3003)
data2011.SetFillColor(46)
data2011.SetLineColor(46)
data2011.GetXaxis().SetRangeUser(xmin, xmax)
data2011.GetXaxis().SetTitle("m_{#mu#mu}")
data2011.Draw("HISTE")
data2012.SetFillStyle(3004)
data2012.SetFillColor(38)
data2012.SetLineColor(38)
data2012.GetXaxis().SetRangeUser(xmin, xmax)
data2012.GetXaxis().SetTitle("m_{#mu#mu}")
data2012.Draw("HISTEsame")
leg = ROOT.TLegend(.60, .75, .98, .98)
leg.SetFillColor(0)
leg.SetTextSize(0.04)
leg.AddEntry(data2012, name2012, "f")
leg.AddEntry(data2011, name2011, "f")
leg.Draw("same")

pad2.cd()

diff = data2012.Clone("diff")
diff.SetLineColor(30)
#diff.SetMarkerStyle(20)
diff.Add(data2011, -1)
diff.SetTitle("2012-2011")
diff.GetXaxis().SetLabelSize(0.07)
diff.GetYaxis().SetLabelSize(0.07)
diff.GetXaxis().SetTitleSize(0.08)
diff.Draw()
xax=diff.GetXaxis()
line = ROOT.TLine()
line.DrawLine(xmin, 0., xmax, 0.)

name = "zmm_DATA_2011_2012"
c1.Print(name + ".eps")
c1.Print(name + ".png")

