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

### ******* El ID eff *******
### *************************


def fillFile(name):

    if(name=="ElIDEff"):
        arrayEta  = array.array("f", [0.00, 0.80, 1.44, 1.56, 2.00, 2.50]) # Define eta and Pt ranges
        arrayPt = array.array("f", [10, 15, 20, 30, 40, 50, 200])

        ## #Range di Pt [10-15, 15-20, 20-30, 30-40, 40-50, 50-200]         
        effArray = [
                       [0.989, 1.063, 1.034, 1.019, 1.016, 1.007],   # 0.00 < eta < 0.80  
                       [1.125, 1.046, 1.010, 1.008, 1.003, 0.996],   # 0.80 < eta < 1.44 
                       [1.313, 1.016, 1.051, 0.992, 0.980, 1.018],   # 1.44 < eta < 1.56
                       [0.823, 1.100, 1.007, 0.998, 1.006, 1.002],   # 1.56 < eta < 2.00
                       [1.450, 1.173, 1.096, 1.053, 1.028, 1.011]   # 2.00 < eta < 2.50
                   ]

### ******* Ele8 trigger eff *******
### *************************
    elif(name == "Ele8TriggerEff"):
        arrayEta  = array.array("f", [0.00, 0.80, 1.44, 1.55, 2.00, 2.50])
        arrayPt = array.array("f", [10, 20, 40, 200])
        #Range di Pt [10-20, 20-40, 40-200]         
        effArray = [
             [1.062, 1.005, 0.999],   # 0.00 < eta < 0.80  
             [1.038, 0.995, 0.999],   # 0.80 < eta < 1.44 
             [1.,1.,1.],              # 1.44 < eta < 1.55
             [1.010, 0.997, 0.998],   # 1.55 < eta < 2.00
             [1.013, 0.999, 0.998]   # 2.00 < eta < 2.50
            ]

### ******* Ele17 trigger eff *******
### *************************
    elif(name == "Ele17TriggerEff"):
        arrayEta  = array.array("f", [0.00, 0.80, 1.44, 1.55, 2.00, 2.50]) # Define eta and Pt ranges
        arrayPt = array.array("f", [10, 20, 40, 200])
        ## #Range di Pt [10-20, 20-40, 40-200]         
        effArray = [
            [0.868, 1.003, 0.998],   # 0.00 < eta < 0.80  
            [0.751, 0.994, 0.999],   # 0.80 < eta < 1.44 
            [1.,1.,1.],              # 1.44 < eta < 1.55
            [0.923, 0.994, 0.997],   # 1.55 < eta < 2.00
            [0.886, 1.001, 1.000]   # 2.00 < eta < 2.50
            ]
        
        
    title = name
    
    histoEff = ROOT.TH2F(name, title, len(arrayEta)-1, arrayEta, len(arrayPt)-1, arrayPt)


    for iEta, values in enumerate(effArray):
        [histoEff.SetBinContent(iEta+1, iPt+1, weight) for iPt, weight in enumerate(values)]

    ROOT.SetOwnership(histoEff, True)    



wfile = ROOT.TFile.Open("Weights.root", "UPDATE")

fillFile("Ele17TriggerEff")
fillFile("Ele8TriggerEff")
fillFile("ElIDEff")

wfile.Write()
wfile.Close()
