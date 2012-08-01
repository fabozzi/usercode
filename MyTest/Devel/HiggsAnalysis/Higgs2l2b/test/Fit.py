
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

c1 = ROOT.TCanvas()
#c2 = ROOT.TCanvas()
#c3 = ROOT.TCanvas()
if(run == 'A'):intLumi = lumiA
elif(run == 'B'):intLumi = lumiB
else: intLumi = lumiA + lumiB

intLumi = 20000
lumiFact = intLumi/1000
for k, v in scale.iteritems():
        scale[k] = v * lumiFact
	

### *************
### Configuration
### *************

outfilename = "Fit_jpt40_30_lpt_40_20_TCHE_minMET55_3rdJveto"

info = "\n\nLUMINOSITY: "+ str(lumiFact)
#info += "\nLeptons: pt>40/20, Jets pt>40/20\nCut on zzdpt (-40, 50), CSV, no cut on zzdphi and trkMet < 60\nSignal and background shapes taken by 0btag for drellyan, fitting just normalization\nBest candidate choise: candidate with greater leading jet pt\nData sample to fit generated according to MC\n"

info += "\nzzdpt (-40, 50), TCHE, no cut on zzdphi, trkMet < 60\nSig: gauss. Bkg from MC. Fitting just normalization\nBest candidate choise: candidate with candidate with ll mass closer to Z mass\nData sample to fit generated according to the model\n"

xmin = 0
xmax = 400
samples = [ "WW", "WZ",  'T-tW', 'Tbar-tW',"ZZ", "TT", "DYJets" ]
var = "zjjmass_0btag"
var2 = "zjjmass_2btag"
var3 = "zjjmass_sbTT_2btag"

outfile = open(outfilename+".txt", "w")

all_files = dict( [ (s, ROOT.TFile.Open("NoNorm_Mu/"+s+"EdmNtp.root") )for s in samples])
all_files_El = dict( [ (s, ROOT.TFile.Open("NoNorm_El/"+s+"EdmNtp.root") )for s in samples])
all_hist = dict( [ (s, f.Get(var2) ) for s,f in all_files.iteritems()])
all_hist_El = dict( [ (s, f.Get(var2) ) for s,f in all_files_El.iteritems()])

x = ROOT.RooRealVar("x","x",xmin,xmax)
nbins = 40
x.setBins(nbins)

rebin = 10

[h.Add(h2) for h,h2 in zip(all_hist.itervalues(), all_hist_El.itervalues())]

print "ZZ events (all range) prescaling: ", all_hist["ZZ"].Integral()

zz_Hist = all_hist["ZZ"].Clone()
zz_Hist.Rebin(rebin)

#ttSB["TT"].Rebin(5)


[h.Scale(scale[k]) for k, h in all_hist.iteritems()]

### ****************************************************************************
### Check on numer of Signal events expected for this lumi
### ****************************************************************************

print "ZZ events (all range) afterscaling: ", all_hist["ZZ"].Integral()
print "-"*20


print "Events between 60 and 120 GeV"
print "ZZ events: ", all_hist["ZZ"].Integral(60,120)
print "top events: ", (all_hist["TT"].Integral(60, 120) + all_hist["T-tW"].Integral(60,120) + all_hist["Tbar-tW"].Integral(60,120) )
print "ZJets events: ", all_hist["DYJets"].Integral(60,120)
print "-"*20


print "-"*20

### ********************
### Write info in a file
### ********************

outfile.write(info)
outfile.write( "-"*20 )
outfile.write("\nEvents expected from MC")
outfile.write("\nZZ events: "+ str(all_hist["ZZ"].Integral(xmin, xmax)))
outfile.write("\nZJets events: "+ str(all_hist["DYJets"].Integral(xmin, xmax)))
outfile.write("\nTT events: "+  str( (all_hist["TT"].Integral(xmin, xmax) + all_hist["T-tW"].Integral(xmin, xmax) + all_hist["Tbar-tW"].Integral(xmin, xmax) ) ))
outfile.write( "\n"+"-"*20)


outfile.write(info)
outfile.write( "-"*20 )
outfile.write("\nEvents expected from MC between 60 and 120 GeV")
outfile.write("\nZZ events: "+ str(all_hist["ZZ"].Integral(60, 120)))
outfile.write("\nZJets events: "+ str(all_hist["DYJets"].Integral(60, 120)))
outfile.write("\nTT events: "+  str( (all_hist["TT"].Integral(60, 120) + all_hist["T-tW"].Integral(60,120) + all_hist["Tbar-tW"].Integral(60, 120) ) ))

sig =  all_hist["ZZ"].Integral(60, 120)
bkg = all_hist["TT"].Integral(60, 120) + all_hist["T-tW"].Integral(60,120) + all_hist["Tbar-tW"].Integral(60, 120) + all_hist["DYJets"].Integral(60, 120)
significance = sig / ROOT.TMath.Sqrt(sig + bkg)
outfile.write("\n\nSignificance: "+ str(significance))
outfile.write( "\n\n"+"-"*20)

print "ZZ events (all range): ", all_hist["ZZ"].Integral()
print "ZZ events: ", all_hist["ZZ"].Integral(xmin, xmax)
print "top events: ", (all_hist["TT"].Integral(xmin, xmax) + all_hist["T-tW"].Integral(xmin, xmax) + all_hist["Tbar-tW"].Integral(xmin, xmax) )
print "ZJets events: ", all_hist["DYJets"].Integral(xmin, xmax)
print "-"*20


print "Check over ZZ pre TT",  all_hist["ZZ"].Integral(xmin, xmax)
print "Check over DY pre TT", all_hist["DYJets"].Integral(xmin, xmax)

### ****************************************************************************
### TT background from MC, including single top contribution --> JZB sidebands
### ****************************************************************************


ttSB = dict( [ (s, f.Get(var2).Clone() ) for s,f in all_files.iteritems()])
ttSB_El = dict( [ (s, f.Get(var2).Clone() ) for s,f in all_files_El.iteritems()])

print "Check over ZZ after TT pre DY",  all_hist["ZZ"].Integral(xmin, xmax)

#ttSB = dict( [ (s, f.Get(var3) ) for s,f in all_files.iteritems()])
#ttSB_El = dict( [ (s, f.Get(var3) ) for s,f in all_files_El.iteritems()])
[h.Add(h2) for h,h2 in zip(ttSB.itervalues(), ttSB_El.itervalues())]
### Normalize to the same lumi
ttL = {"TT":157.5, "Tbar-tW":7.46, "T-tW":7.46}
L = min(ttL.itervalues())
ttScaling  = dict( [ (s,ttSB[s].Integral()*L/ttL[s]) for s in ttL.keys()] )
[h.Scale(scale[k]) for k, h in zip(ttL.keys(), ttSB.itervalues())]
### Merge TT, Tbar-tW and T-tW channels
ttSB["TT"].Add(ttSB["Tbar-tW"])
ttSB["TT"].Add(ttSB["T-tW"])

#ttSB_ = ttSB["TT"].Clone("tt_SB")

#ttSB_ = ttSB["TT"].Clone("tt_SB")

ttSB_ = all_hist["TT"].Clone("TT")
ttSB_.Add(all_hist["Tbar-tW"])
ttSB_.Add(all_hist["T-tW"])
ttSB_.Rebin(rebin)


#ttSB["TT"].Rebin(rebin)
#ttSB_.Rebin(rebin)
print "TT scaled: ", ttSB["TT"].Integral() 
print "TT scaled (clone): ", ttSB_.Integral()

### ****************************************************************************
### Z+Jets shape from MC --> 0 btag category
### ****************************************************************************

print "Check over ZZ pre DY", all_hist["ZZ"].Integral(xmin, xmax)
print "Check over DY pre DY", all_hist["DYJets"].Integral(xmin, xmax)

#DY = all_files["DYJets"].Get(var2).Clone()
#DY_El = all_files_El["DYJets"].Get(var2).Clone()
#DY.Add(DY_El)
DY = all_hist["DYJets"].Clone("DY")
DY.Rebin(rebin)



print "Check over ZZ after DY", all_hist["ZZ"].Integral(xmin, xmax)
print "Check over DY after DY", all_hist["DYJets"].Integral(xmin, xmax)
print "DY scaled: ", DY.Integral(xmin, xmax)

### ****************************************************************************
### Data Sample
### ****************************************************************************

#all_hist["ZZ"].GetXaxis().SetRangeUser(60.,120.)

#fdata =  ROOT.TFile.Open("NoNorm_Mu/MuRun2011AllEdmNtp.root")
#fdataEl =  ROOT.TFile.Open("NoNorm_El/ElRun2011AllEdmNtp.root")
#data = fdata.Get(var2).Clone()
#dataEl = fdataEl.Get(var2).Clone()
#data.Add(dataEl)
#data.Rebin(rebin)


### **********
### Rebin MC
### **********

print "Check over ZZ preMerge and rebin", all_hist["ZZ"].Integral()
print "Check over DY preMerge and rebin", all_hist["DYJets"].Integral()

[h.Rebin(rebin) for h in all_hist.itervalues()]

print "Check over ZZ preMerge and after rebin", all_hist["ZZ"].Integral()
print "Check over DY preMerge and after rebin", all_hist["DYJets"].Integral()

### ****************************************************************************
### Preparing histos for shape
### ****************************************************************************

#all_DataHist = dict([(s, ROOT.RooDataHist(var, "", ROOT.RooArgList(x), h) ) for s,h in all_hist.iteritems() ] )
#all_DataHist["DYJets"] = ROOT.RooDataHist(var, "", ROOT.RooArgList(x), all_hist["DYJets"])
#ttSBHist = ROOT.RooDataHist(var, "", ROOT.RooArgList(x), ttSB["TT"])


ttSB_.Smooth(50)
ttSBHist = ROOT.RooDataHist(var, "", ROOT.RooArgList(x), ttSB_)
DY.Smooth(50)
DYHist = ROOT.RooDataHist(var, "", ROOT.RooArgList(x), DY)
zz_Hist.Smooth(50)
zz = ROOT.RooDataHist(var, "", ROOT.RooArgList(x), zz_Hist)

print "Check over ZZ preMerge", all_hist["ZZ"].Integral()
print "Check over DY preMerge", all_hist["DYJets"].Integral()

sample = all_hist["TT"].Clone("sample")
sample.Add(all_hist["Tbar-tW"])
sample.Add(all_hist["T-tW"])
sample.Add(all_hist["DYJets"])
#bkgttHist = ROOT.RooDataHist(var, "", ROOT.RooArgList(x), bkgtt)
sample.Add(all_hist["ZZ"])

print "Check over ZZ", all_hist["ZZ"].Integral()
print "Check over TT", all_hist["TT"].Integral()
print "Check over DY", all_hist["DYJets"].Integral()

### ****************************************************************************
### sample to fit (MC)
### ****************************************************************************


#sample.Rebin(rebin)
sample.GetXaxis().SetRangeUser(xmin,xmax)
nEvents = int(sample.Integral())
#sample.Scale(1/nEvents)
sampleHist = ROOT.RooDataHist(var2, "", ROOT.RooArgList(x), sample)
samplePDF =  ROOT.RooHistPdf("samplePdf", "samplePdf", ROOT.RooArgSet(x), sampleHist)

nsamp = ROOT.RooRealVar("nsamp", "sample", nEvents)
esample = ROOT.RooExtendPdf("esample", "esample", samplePDF, nsamp)

nEvents = int(sample.Integral())
print "number of events in MC sample", sample.Integral()
print "number of events to Generate", nEvents
#sampleToFitHist = esample.generateBinned(ROOT.RooArgSet(x), nEvents)

#sampleToFitTH1 = sampleToFit.createHistogram("sampleToFit")
#sampleToFitTH1.Rebin(rebin)
#sampleToFitHist = ROOT.RooDataHist(var, "", ROOT.RooArgList(x), sampleToFitTH1 )
#sampleToFitHist = sampleToFit.binnedClone()

#histdata = ROOT.RooDataHist(var, "", ROOT.RooArgList(x), data)

mean_gen = ROOT.RooRealVar("mean","mean",89.)
sigma_gen = ROOT.RooRealVar("sigma","sigma",11.)

mean = ROOT.RooRealVar("mean","mean",89.)
sigma = ROOT.RooRealVar("sigma","sigma",11., 0., 20.)

gauss_gen= ROOT.RooGaussian("gauss","gaussian PDF",x,mean_gen,sigma_gen) ;

gauss= ROOT.RooGaussian("gauss","gaussian PDF",x,mean,sigma) ;



## mean = ROOT.RooRealVar("mean","mean",89.)
## sigma = ROOT.RooRealVar("sigma","sigma",11., 0., 50.) 
## gauss= ROOT.RooGaussian("gauss","gaussian PDF",x,mean,sigma) ;

## mean = ROOT.RooRealVar("mean","mean",300.)
## sigma = ROOT.RooRealVar("sigma","sigma",11., 0., 50.) 
## gauss= ROOT.RooGaussian("gauss","gaussian PDF",x,mean,sigma) ;

### ****************************************************************************
### Preparing shapes
### ****************************************************************************

zjetsPDF = ROOT.RooHistPdf("zjetsPdf", "zjetsPdf", ROOT.RooArgSet(x), DYHist)
#zjGenPDF = ROOT.RooHistPdf("zjGenPdf", "zjGenPdf", ROOT.RooArgSet(x), DYHist)
TTPDF = ROOT.RooHistPdf("TTPdf", "TTPdf", ROOT.RooArgSet(x), ttSBHist)
zzPDF = ROOT.RooHistPdf("zzPdf", "zzPdf", ROOT.RooArgSet(x), zz)

zjetsPDF.setInterpolationOrder(2)
TTPDF.setInterpolationOrder(2)


bzjets_2 = ROOT.RooRealVar("bzjets", "background ZJets fraction", 1625)
btt_2 = ROOT.RooRealVar("btt", "background tt fraction", 552)
s_2 = ROOT.RooRealVar("s", "signal fraction", 50)

bzjets = ROOT.RooRealVar("bzjets", "background ZJets fraction", 1625, 0, 5000)
btt = ROOT.RooRealVar("btt", "background tt fraction", 522, 0, 2000)
s = ROOT.RooRealVar("s", "signal fraction", 50, 0., 3000.)

## bzjets = ROOT.RooRealVar("bzjets", "background ZJets fraction", 712.392543296, 0, 2000)
## btt = ROOT.RooRealVar("btt", "background tt fraction", 286.366024015, 0, 500)
## s = ROOT.RooRealVar("s", "signal fraction", 19.9752820263, 0., 100.)

## bzjets = ROOT.RooRealVar("bzjets", "background ZJets fraction", 633, 0, 2000)
## btt = ROOT.RooRealVar("btt", "background tt fraction", 383, 0, 500)
## s = ROOT.RooRealVar("s", "signal fraction", 18, 0., 100.)

esig = ROOT.RooExtendPdf("esig", "esig", gauss, s)
esig_2 = ROOT.RooExtendPdf("esig", "esig", gauss_gen, s_2)
#esig = ROOT.RooExtendPdf("esig", "esig", zzPDF, s)

ebkg_zjets = ROOT.RooExtendPdf("ebkg_zjets", "ebkg_zjets", zjetsPDF, bzjets)
ebkg_tt = ROOT.RooExtendPdf("ebkg_tt", "ebkg_tt", TTPDF, btt)

ebkg_zjets_2 = ROOT.RooExtendPdf("ebkg_zjets", "ebkg_zjets", zjetsPDF, bzjets_2)
ebkg_tt_2 = ROOT.RooExtendPdf("ebkg_tt", "ebkg_tt", TTPDF, btt_2)

bkgPDF = ROOT.RooAddPdf("bkgPDF", "bkgPDF", ROOT.RooArgList( ebkg_tt, ebkg_zjets ))
#model = ROOT.RooAddPdf("model", "esig + ebkg_tt",  ROOT.RooArgList(esig, ebkg_tt) )

mean2 = ROOT.RooRealVar("mean2","mean",300)
sigma2 = ROOT.RooRealVar("sigma2","sigma",11., 0., 50.) 
gauss2= ROOT.RooGaussian("gauss2","gaussian PDF",x,mean,sigma) ;
s2 = ROOT.RooRealVar("s", "signal fraction", 800)
esigTest = ROOT.RooExtendPdf("esigTest", "esig", gauss2, s2)

model2 = ROOT.RooAddPdf("model2", "esigTest + ebkg_tt +ebkg_zjets",  ROOT.RooArgList(esig_2, ebkg_tt_2, ebkg_zjets_2) )

model= ROOT.RooAddPdf("model", "esig + ebkg_tt +ebkg_zjets",  ROOT.RooArgList(esig, ebkg_tt, ebkg_zjets) )


rnd = ROOT.TRandom3()

nEvts = rnd.Poisson(nEvents)
print "generated number of events: ", nEvts

sDist = ROOT.TH1F("","",300, 0, 300)
pullDist = ROOT.TH1F("","",100, -10, 10)
#sDist.Rebin(rebin)

#sDist.GetXaxis().SetRangeUser(xmin,xmax)

for i in range(1000):
	ROOT.RooRandom.randomGenerator().SetSeed(i) 
	sampleToFitHist = model2.generateBinned(ROOT.RooArgSet(x), nEvts)
	fr = model.fitTo(sampleToFitHist, ROOT.RooFit.SumW2Error(ROOT.kFALSE))
	sDist.Fill(s.getVal())
	#	pullDist.Fill((s.getVal() - 800)/s.getError())
	pullDist.Fill((s.getVal() - all_hist["ZZ"].Integral(xmin, xmax))/s.getError())
	
sDist.Draw()		   
c1.SaveAs("nsigDistr.eps")
c1.SaveAs("nsigDistr.png")
pullDist.Draw()
c1.SaveAs("pullNsigDistr.eps")
c1.SaveAs("pullNsigDistr.png")
print "-------------- stoys integral "
print str(sDist.Integral())
	
xframe = x.frame()
model.Print("t")
print "-"*100
#model.plotOn(xframe)
#sampleHist.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kGray + 1),  ROOT.RooFit.MarkerColor(ROOT.kGray + 1), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2) )

sampleHist.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kGray + 1),  ROOT.RooFit.MarkerColor(ROOT.kGray + 1))
sampleToFitHist.plotOn(xframe )
#sample.plotOn(xframe)


### ***********************************
### Do a fit and plot dataset the model
### ***********************************
fr = model.fitTo(sampleToFitHist, ROOT.RooFit.SumW2Error(ROOT.kFALSE))


model.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kRed))
model.plotOn(xframe, ROOT.RooFit.Components("esig"),ROOT.RooFit.LineColor(ROOT.kOrange +8) )
model.plotOn(xframe, ROOT.RooFit.Components("ebkg_tt"),ROOT.RooFit.LineColor(ROOT.kMagenta +3))
#model.plotOn(xframe, ROOT.RooFit.Components("ebkg_zjets"),ROOT.RooFit.LineColor(ROOT.kGreen -3))
bkgPDF.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kGreen -3))
nDOF = nbins - 2
print "-"*20
print "Chi2: ", xframe.chiSquare()
print "-"*20
probchi2 = ROOT.TMath.Prob(xframe.chiSquare()*nDOF, nDOF)
print "ProbChi2: ", probchi2
print "-"*20

### ********************
### Write info in a file
### ********************
outfile.write("\n\nModel after fit\n")
outfile.write("\n  s: "+str(s.getVal()) + " +/- " + str(s.getError()))
outfile.write("\n  btt: "+str(btt.getVal()) + " +/- " + str(btt.getError()))
outfile.write("\n  bzjets: "+str(bzjets.getVal()) + " +/- " + str(bzjets.getError()))
sign =  s.getVal() / s.getError()

print "significance: ", sign
relErr =  s.getError()/ s.getVal()
outfile.write("\n\nRelative error: "+ str(relErr))
outfile.write("\n\nSignificance: "+str(sign))
print "s: ",s.Print("t")
print "btt: ",btt.Print()
print "bzjets: ",bzjets.Print()

xframe.Draw()
c1.Print(outfilename+".eps")
c1.Print(outfilename+".png")

del xframe
[f.Close() for f in all_files_El.itervalues()]
f.Close()
#fdataEl.Close()
