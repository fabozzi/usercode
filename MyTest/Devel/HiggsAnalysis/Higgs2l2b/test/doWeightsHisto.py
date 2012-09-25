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


names = [  "MuIDEffTight", "EleIDEffLoose", "Ele8TriggerEff", "Ele17TriggerEff" ]

### ******* El ID eff: data/MC SF *******
### *************************************
name_EleWPLoose = "EleIDEffLoose"

arrayEta_EleWPLoose  = array.array("f", [0.00, 0.80, 1.442, 1.556, 2.00, 2.50]) # Define eta and Pt ranges
arrayPt_EleWPLoose = array.array("f", [10, 15, 20, 30, 40, 50, 200])

#Range di Pt            [10-15, 15-20, 20-30, 30-40, 40-50, 50-200]         
effArray_EleWPLoose = [
                        [0.861, 0.950, 1.017, 1.019, 1.016, 1.005],   # 0.00 < eta < 0.80  
                        [0.864, 0.965, 0.998, 1.008, 1.002, 0.993],   # 0.80 < eta < 1.44 
                        [1.019, 1.068, 1.154, 1.028, 0.989, 1.005],   # 1.44 < eta < 1.56
                        [0.823, 0.925, 1.009, 1.010, 1.009, 1.003],   # 1.56 < eta < 2.00
                        [1.109, 1.152, 1.102, 1.066, 1.040, 1.019]    # 2.00 < eta < 2.50
                      ]

errArray_EleWPLoose = [
                        [0.0645, 0.0249, 0.0053, 0.0022, 0.0005, 0.0030],   # 0.00 < eta < 0.80  
                        [0.0657, 0.0237, 0.0062, 0.0026, 0.0019, 0.0037],   # 0.80 < eta < 1.44 
                        [0.2760, 0.0812, 0.0232, 0.0111, 0.0101, 0.0175],   # 1.44 < eta < 1.56
                        [0.0784, 0.0344, 0.0098, 0.0047, 0.0036, 0.0062],   # 1.56 < eta < 2.00
                        [0.1184, 0.0378, 0.0119, 0.0021, 0.0052, 0.0082]    # 2.00 < eta < 2.50
                      ]

### ******* Mu ID eff: data/MC SF *******
### *************************************
name_MuWPTight = "MuIDEffTight"

arrayEta_MuWPTight  = array.array("f", [0.00, 0.80, 2.10, 2.40]) # Define eta and Pt ranges
arrayPt_MuWPTight = array.array("f", [20, 40, 100])

#Range di Pt            [20-40, 40-100]         
effArray_MuWPTight = [
                        [1.00425, 1.00119],   # 0.00 < eta < 0.80  
                        [1.00740, 1.00425],   # 0.80 < eta < 2.10 
                        [1.02159, 1.01404]   # 2.10 < eta < 2.40
                      ]

errArray_MuWPTight = [
                        [0.00042, 0.00039],   # 0.00 < eta < 0.80  
                        [0.00046, 0.00039],   # 0.80 < eta < 2.10 
                        [0.00142, 0.00141]   # 2.10 < eta < 2.40
                      ]


### ******* Ele8 trigger eff *******
### *******************************
name_EleTrig8 = "Ele8TriggerEff"

### Define eta and Pt ranges

arrayEta_EleTrig8  = array.array("f", [0.00, 0.80, 1.4442, 1.566, 2.00, 2.50])
arrayPt_EleTrig8 = array.array("f", [10, 20, 40, 200])

#Range di Pt       [10-20, 20-40, 40-200]         
effArray_EleTrig8 = [
                      [0.9545, 0.9830, 0.9889],   # 0.00 < eta < 0.80  
                      [0.8521, 0.9316, 0.9715],   # 0.80 < eta < 1.44 
                      [1.,1.,1.],                 # 1.44 < eta < 1.55
                      [0.8387, 0.8948, 0.9380],   # 1.56 < eta < 2.00
                      [0.8677, 0.9331, 0.9508]    # 2.00 < eta < 2.50
                    ]

errArray_EleTrig8 = [
                      [0.0038, 0.0005, 0.0004],   # 0.00 < eta < 0.80  
                      [0.0062, 0.0011, 0.0007],   # 0.80 < eta < 1.44 
                      [0.,0.,0.],                 # 1.44 < eta < 1.55
                      [0.0106, 0.0021, 0.0014],   # 1.55 < eta < 2.00
                      [0.0102, 0.0019, 0.0016]    # 2.00 < eta < 2.50
                    ]

### ******* Ele17 trigger eff *******
### *********************************
name_EleTrig17 = "Ele17TriggerEff"

### Define eta and Pt ranges

arrayEta_EleTrig17  = array.array("f", [0.00, 0.80, 1.4442, 1.566, 2.00, 2.50])
arrayPt_EleTrig17 = array.array("f", [10, 20, 40, 200])

#Range di Pt       [10-20, 20-40, 40-200]         
effArray_EleTrig17 = [
                       [0.4735, 0.9856, 0.9913],   # 0.00 < eta < 0.80  
                       [0.3426, 0.9360, 0.9763],   # 0.80 < eta < 1.44 
                       [1.,1.,1.],                 # 1.44 < eta < 1.56
                       [0.4439, 0.9006, 0.9447],   # 1.56 < eta < 2.00
                       [0.4519, 0.9444, 0.9624]    # 2.00 < eta < 2.50
                     ]

errArray_EleTrig17 = [
                       [0.0089, 0.0004, 0.0003],   # 0.00 < eta < 0.80  
                       [0.0082, 0.0011, 0.0006],   # 0.80 < eta < 1.44 
                       [0.,0.,0.],                 # 1.44 < eta < 1.55
                       [0.0142, 0.0020, 0.0014],   # 1.55 < eta < 2.00
                       [0.0148, 0.0018, 0.0014]    # 2.00 < eta < 2.50
                     ]



arraysEta = {
    "MuIDEffTight" : arrayEta_MuWPTight,
    "EleIDEffLoose" : arrayEta_EleWPLoose,
    "Ele8TriggerEff" : arrayEta_EleTrig8,
    "Ele17TriggerEff" : arrayEta_EleTrig17,
   }

arraysPt = {
    "MuIDEffTight" : arrayPt_MuWPTight,
    "EleIDEffLoose" : arrayPt_EleWPLoose,
    "Ele8TriggerEff" : arrayPt_EleTrig8,
    "Ele17TriggerEff" : arrayPt_EleTrig17,
   }

effArrays = {
    "MuIDEffTight" : effArray_MuWPTight,
    "EleIDEffLoose" : effArray_EleWPLoose,
    "Ele8TriggerEff" : effArray_EleTrig8,
    "Ele17TriggerEff" : effArray_EleTrig17,
   }

errArrays = {
    "MuIDEffTight" : errArray_MuWPTight,
    "EleIDEffLoose" : errArray_EleWPLoose,
    "Ele8TriggerEff" : errArray_EleTrig8,
    "Ele17TriggerEff" : errArray_EleTrig17,
   }



#wfile = ROOT.TFile.Open("Weights.root", "UPDATE")
wfile = ROOT.TFile.Open("Weights.root", "RECREATE")

for name in names:

    print name
    print len(arraysEta[name])-1
    print arraysEta[name]
    print len(arraysPt[name])-1
    print arraysPt[name]
    
    title = name
    histoEff = ROOT.TH2F(name, title, len(arraysEta[name])-1, arraysEta[name], len(arraysPt[name])-1, arraysPt[name])

    for iEta, values in enumerate(effArrays[name]):
        [histoEff.SetBinContent(iEta+1, iPt+1, weight) for iPt, weight in enumerate(values)]
        print iEta, values
        print iPt, weight
    for iEta, values in enumerate(errArrays[name]):
        [histoEff.SetBinError(iEta+1, iPt+1, weight) for iPt, weight in enumerate(values)]
    histoEff.Write()
    
#    histoEff = ROOT.TH2F(name, title, len(arrayEta)-1, arrayEta, len(arrayPt)-1, arrayPt)
#
#    for iEta, values in enumerate(effArray):
#        [histoEff.SetBinContent(iEta+1, iPt+1, weight) for iPt, weight in enumerate(values)]

#wfile.Write()
wfile.Close()
