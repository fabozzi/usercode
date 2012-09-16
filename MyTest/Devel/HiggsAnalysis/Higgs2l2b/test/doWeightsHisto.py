import sys
sys.argv.append('-b')
import ROOT
import array

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()  # don't pop up canvases


### INSTRUCTIONS
### Use this scripts to create root files containing histograms with weghts
### depending on eta and pt. To use it, please, check the Eta and Pt ranges definition.
### Set the name of the histo.
### Set the values in effArray, it is an array of arrays.
name = "ElIDEff"
title = name

wfile = ROOT.TFile.Open("Weights.root", "RECREATE")

### Define eta and Pt ranges

arrayEta  = array.array("f", [0.00, 0.80, 1.44, 1.56, 2.00, 2.50])
arrayPt = array.array("f", [10, 15, 20, 30, 40, 50, 200])

#Range di Pt [10-15, 15-20, 20-30, 30-40, 40-50, 50-200]         
effArray = [
             [0.989, 1.063, 1.034, 1.019, 1.016, 1.007],   # 0.00 < eta < 0.80  
             [1.125, 1.046, 1.010, 1.008, 1.003, 0.996],   # 0.80 < eta < 1.44 
             [1.313, 1.016, 1.051, 0.992, 0.980, 1.018],   # 1.44 < eta < 1.56
             [0.823, 1.100, 1.007, 0.998, 1.006, 1.002],   # 1.56 < eta < 2.00
             [1.450, 1.173, 1.096, 1.053, 1.028, 1.011],   # 2.00 < eta < 2.50
            ]

histoEff = ROOT.TH2F(name, title, len(arrayEta)-1, arrayEta, len(arrayPt)-1, arrayPt)


for iEta, values in enumerate(effArray):
    [histoEff.SetBinContent(iEta+1, iPt+1, weight) for iPt, weight in enumerate(values)]

wfile.Write()
wfile.Close()
