import sys
sys.argv.append('-b')
import ROOT
import os, commands
from settings_2012 import *
from SamplesColors import *
from HiggsAnalysis.Higgs2l2b.scaleFactors import scale

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')

ROOT.gROOT.SetBatch()        # don't pop up canvases

ROOT.gStyle.SetPalette(1);
ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas()


chain = ROOT.TChain("Events")
chain.Add("/data3/scratch/users/decosa/Higgs/Summer12/GluGluToHToZZTo2L2Q_M-275_8TeV/edmntp12Jul12Cleaned/h2l2q_ntuple_clean_0.root")
chain.Add("/data3/scratch/users/decosa/Higgs/Summer12/GluGluToHToZZTo2L2Q_M-275_8TeV/edmntp12Jul12Cleaned/h2l2q_ntuple_clean_1.root")

def get2DPlot(var):
    chain.Draw("muHiggszjjMass:"+var+">>histo(80, 225., 425., 60, 70., 110.)","muHiggsMatch == 1", "COL")
    histo = chain.GetHistogram()
    histo.GetYaxis().SetTitle("m_{jj} [GeV/c^{2}]")
    histo.GetXaxis().SetTitle("m_{ZZ} [GeV/c^{2}]")
    histo.SetTitle("")
    histo.Draw("COL")
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextAlign(11)
    
    latex.DrawLatex(0.1, 0.93, "CMS Simulation 2012");
    latex.DrawLatex(0.79,0.93, "#sqrt{s} = 8 TeV");
    nameFile = "with_KinFit"
    if (var=="muHiggsMass"): nameFile = "wo_KinFit"
    c1.Print(nameFile+".png")
    c1.Print(nameFile+".eps")



get2DPlot("muHiggsMass")
get2DPlot("muHiggsRefitMass")
