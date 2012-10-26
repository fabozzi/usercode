
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
ROOT.TH1.SetDefaultSumw2()


if(run == 'A'):intLumi = lumiA
elif(run == 'B'):intLumi = lumiB
else: intLumi = lumiA + lumiB

lumiFact = intLumi/1000.

#lumiFact = 20

for k, v in scale.iteritems():
    scale[k] = v #* lumiFact

cmd = "mkdir "+folder
os.system(cmd)

sigSamples = ["GluGluToHToZZTo2L2Q_M-350_8TeV"]

drawSignal= False
def makePlot(varKey):

    norm = 0
    print "Variable: ", varKey
    var = varKey
    rebin = Vars[varKey][0]
    xmin = Vars[varKey][1]
    xmax = Vars[varKey][2]
    logScale = Vars[varKey][3]
    label = Vars[varKey][4]
    
    c1 = ROOT.TCanvas()
    m = 0.
    leg = ROOT.TLegend(.6, .7, .9, .9)
    leg.SetFillColor(0)
    leg.SetTextSize(0.03)

    # create a dictionary of histograms for variable var for each sample

    hs = ROOT.THStack("hs",var)
 
    all_hist = dict([ (s, ROOT.TFile.Open("NoNorm_"+ch+"/" + s+"EdmNtp.root").Get(var)) for s in samples])
    [h.Rebin(rebin) for h in all_hist.values()]
    [h.Scale(scale[k]) for k, h in all_hist.iteritems()]
    ## for s in samples:
    ##     print s, scale[s]
        

    sig_samples = dict([ (s, ROOT.TFile.Open("NoNorm_"+ch+"/" + s+"EdmNtp.root").Get(var)) for s in sigSamples])
    [h.Rebin(rebin) for h in sig_samples.values()]
    [h.Scale(scale[k]) for k, h in sig_samples.iteritems()]
    
    if (ch == 'All'):
        all_hist_2 = dict([ (s, ROOT.TFile.Open("NoNorm_Mu/" + s+"EdmNtp.root").Get(var)) for s in samples ])
        [h.Rebin(rebin) for h in all_hist_2.values()]
        [h.Scale(scale[k]) for k, h in all_hist_2.iteritems()]
        [h.Add(h2) for h,h2 in zip(all_hist.itervalues(), all_hist_2.itervalues()) ]


    ## if("DY0" in samples and "DY1" in samples and "DY2" in samples and "DY3" in samples and "DY4" in samples):
    ##     all_hist["DYJetsToLL"] = all_hist["DY0"]
    ##     all_hist["DYJetsToLL"].Add(all_hist["DY1"])
    ##     all_hist["DYJetsToLL"].Add(all_hist["DY2"])
    ##     all_hist["DYJetsToLL"].Add(all_hist["DY3"])
    ##     all_hist["DYJetsToLL"].Add(all_hist["DY4"])
   
    
    #sampleNew = list(samples)
    ##     sampleNew.remove("DY0")
    ##     sampleNew.remove("DY1")
    ##     sampleNew.remove("DY2")
    ##     sampleNew.remove("DY3")
    ##     sampleNew.remove("DY4")
    ##     sampleNew.append("DYJetsToLL")
        
    ### Evaluate number of events after the selection
    

    if data:
        if( run=='A'): dataRun2012 = ROOT.TFile.Open("NoNorm_"+ch+"/"+ch+"Run2012AEdmNtp.root"if ch !="All" else "NoNorm_El/ElRun2012AEdmNtp.root").Get(var)
        elif( run=='B'): dataRun2012 = ROOT.TFile.Open("NoNorm_"+ch+"/"+ch+"Run2012BEdmNtp.root"if ch !="All" else "NoNorm_El/ElRun2012BEdmNtp.root").Get(var)

        else:
            dataRun2012 = ROOT.TFile.Open("NoNorm_"+ch+"/"+ch+"Run2012AllEdmNtp.root"if ch !="All" else "NoNorm_El/ElRun2012AllEdmNtp.root").Get(var)

        if (ch == 'All'):
            if( run=='A'): dataRun2012_2 = ROOT.TFile.Open("NoNorm_Mu/MuRun2012AEdmNtp.root").Get(var)
            elif( run=='B'): dataRun2012_2 = ROOT.TFile.Open("NoNorm_Mu/MuRun2012BEdmNtp.root").Get(var)
            else:
                dataRun2012_2 = ROOT.TFile.Open("NoNorm_Mu/MuRun2012AllEdmNtp.root").Get(var)
            dataRun2012.Add(dataRun2012_2)

        dataRun2012.Rebin(rebin)
        ROOT.SetOwnership(dataRun2012,True)
                  
        #    [h.Scale(1./h.Integral()) if h.Integral()!=0 else h.Scale(1)  for  h in all_hist.itervalues()]
    



    ### Normalization to data 

    if(normalizeToData):
        norm = dataRun2012.Integral()
        
        integ = {}
        integ = dict([ (s, h.Integral()) for s,h in all_hist.iteritems()])    
        normLumi = 0;
        for integral in integ.itervalues():
            normLumi +=integral
        print "normLumi"
        print normLumi
        factNorm = float(norm/normLumi)
    
        #if (norm!=0 and run=='AB' and ch == 'All'):
        if (norm!=0):
            [h.Scale(factNorm) for  h in all_hist.itervalues()]
            #if(ch == 'All'):[h.Scale(factNorm) for  h in all_hist_2.itervalues()]

    if (not normalizeToData):
              [h.Scale(lumiFact) for  h in all_hist.itervalues()]
             # if(ch == 'All'):[h.Scale(lumiFact) for  h in all_hist_2.itervalues()]
   ### END normalization to data

    maximum = max([ h.GetMaximum() for h in all_hist.itervalues()])
    minimum = min([ h.GetMinimum() for h in all_hist.itervalues()])
    
    if(maximum < dataRun2012.GetMaximum()):maximum = dataRun2012.GetMaximum()
    if(minimum > dataRun2012.GetMinimum()):minimum = dataRun2012.GetMinimum()

    lim_max = maximum*1.03

    if(var.startswith("cos") or var.startswith("phi")  or  var.startswith("lepteta") or var.startswith("npv")):lim_max = maximum*1.4
    #if(logScale): lim_max = maximum*100
    if(logScale): lim_max = maximum*10
    
    if(logScale and (var.startswith("npv") or var.startswith("lepteta"))):lim_max = maximum*10
    if(logScale and (var.startswith("leptpt") or var.startswith("jetpt")) ):lim_max = maximum*2
   # if( var.startswith("lepteta")): lim_max = maximum*100
    
    # histogram style setting     
    for i, s in enumerate(samples):
        hist = all_hist[s]
        hist.SetName(hist.GetName()+'_'+ch+'_'+s)
        if(xmin!=xmax): hist.GetXaxis().SetRangeUser(xmin, xmax)
        if(logScale): hist.SetMinimum(0.001)
        if(logScale and minimum!=0): hist.SetMinimum(minimum*0.5)
        hist.SetStats(ROOT.kFALSE)
        hist.SetFillColor(colors[i])
        hist.SetMaximum(lim_max)
        ROOT.SetOwnership(hist,False)
        hist.GetXaxis().SetTitle(label)
        hist.SetTitle(var+"_"+ch)
        hs.Add(hist,'hist')
        leg.AddEntry(hist, s, "f")

    
#    backgrounds = samples[:]
#    backgrounds.remove("ZZ")
#    print "Backgrounds: ", backgrounds
#    print integ
#    #print "All samples: ", samples
#    bkg = 0;
#    for s in backgrounds:
#        bkg +=integ[s]

#    print "Significance: "
#    print  integ["ZZ"]/ROOT.TMath.Sqrt(bkg)
    
    ROOT.SetOwnership(hs,False)
    hs.Draw()
    hs.GetHistogram().GetXaxis().SetTitle(label);
    hs.SetMaximum(lim_max)
   # hs.GetHistogram().SetTitle(var+"_"+ch);
   #    hs.SetTitle(var+"_"+ch);
    hs.SetTitle("");
    hs.Draw()    
    if(xmin!=xmax): hs.GetXaxis().SetRangeUser(xmin, xmax)
    if data :
        dataRun2012.SetMarkerStyle(22)
        dataRun2012.Draw("SAME")
        print "Events in data sample", dataRun2012.Integral()
        if(run == 'A'): leg.AddEntry(dataRun2012, (ch+"Run2012A"), "l")
        if(run == 'B'): leg.AddEntry(dataRun2012, (ch+"Run2012B"), "l")
        if (run == 'AB'): leg.AddEntry(dataRun2012, (ch+"Run2012AB"), "l")


    if(drawSignal):
        hsig =  sig_samples["GluGluToHToZZTo2L2Q_M-350_8TeV"]
        hsig.SetLineColor(ROOT.kPink+3)
        hsig.SetLineWidth(2)
        scaleFact = 1000
        hsig.Scale(scaleFact)
        hsig.Draw("HISTSAME")
        leg.AddEntry(hsig, "H350 (x"+str(scaleFact)+")", "l")
        

    leg.Draw("SAME")
    
    if(logScale): c1.SetLogy()
    name = folder + var+"_"+ch
    #if(logScale): name = folder + var+"_"+ch +"_LOG_"
    if (run == 'A'): name = folder + var+"_"+ch + "Run2012A"
    if (run == 'B'): name = folder + var+"_"+ch + "Run2012B"
    if (run == 'AB'): name = folder + var+"_"+ch + "Run2012"


    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextAlign(11)
    #latex.SetTextAlign(31)
    
    latex.DrawLatex(0.1, 0.93, "CMS preliminary 2012");
    latex.DrawLatex(0.63,0.93, str(intLumi) + " pb^{-1} at #sqrt{s} = 8 TeV");

  
    if(logScale):
        c1.Print(name + "_LOG.eps")
        c1.Print(name + "_LOG.png")
        c1.Print(name + "_LOG.pdf")
        c1.Print(name + "_LOG.C")
    else:
        c1.Print(name + ".eps")
        c1.Print(name + ".png")
        c1.Print(name + ".pdf")
        c1.Print(name + ".C")
        

# make plots for all variables in vars
map(makePlot, Vars.keys())
