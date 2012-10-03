
import sys
sys.argv.append('-b')
import ROOT

from HiggsAnalysis.Higgs2l2b.scaleFactors import scale
from matrix2latex import matrix2latex

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')

ROOT.gROOT.SetBatch()

masses = [200, 210, 220, 230,250, 275, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600]
bkgSamples = [ "WW", "WZ", "ZZ", "TT", "DYJetsToLL_M-50"]
sigSamples = ["GluGluToHToZZTo2L2Q_M-"+str(m)+"_8TeV" for m in masses]
    
### Evaluate the integral between two points
def integ(histo, lower, upper):
    bmin = histo.GetXaxis().FindBin(lower)
    bmax = histo.GetXaxis().FindBin(upper)
    integ = histo.Integral(bmin,bmax)
    integ -= (histo.GetBinContent(bmin)*(lower - histo.GetXaxis().GetBinLowEdge(bmin)))/ histo.GetXaxis().GetBinWidth(bmin)
    integ -= (histo.GetBinContent(bmax)*(histo.GetXaxis().GetBinUpEdge(bmax) -upper))/ histo.GetXaxis().GetBinWidth(bmax)
    return integ


### Evaluate the number of events in a mass window
def evalNumEvts(histosBkg, histosSig, mass):
    
    low = 0.94*mass
    up = 1.10*mass
    signal = "GluGluToHToZZTo2L2Q_M-"+str(mass)+"_8TeV"
    nEvts = dict([(s,integ(h, low, up) ) for s, h in histosBkg.iteritems()])
    nEvts_err = dict([(s,integ(h, low, up) ) for s, h in histosBkg.iteritems()])
    nEvts[signal] =integ(histosSig[signal], low, up)
    nEvts_err[signal] = ROOT.TMath.Sqrt(integ(histosSig[signal], low, up))
    nEvts["Background"] = nEvts["DYJetsToLL_M-50"] + nEvts["TT"]+ nEvts["ZZ"] + nEvts["WZ"] + nEvts["WW"]
    nEvts_err["Background"] = ROOT.TMath.Sqrt(nEvts["DYJetsToLL_M-50"] + nEvts["TT"]+ nEvts["ZZ"] + nEvts["WZ"] + nEvts["WW"])
    nEvts_plus_err =  dict([(s, str.format("{0:.2f}", nEvts[s]) +"$\pm$"+ str.format("{0:.3f}", nEvts_err[s])  ) for s in nEvts.iterkeys() ])
    return nEvts_plus_err
    #print nEvts
    #print "MASS: ", str(mass)
    #print "Event numbers in the mass window: [",str(low),", ",str(up), " ]"

    #for (s, num) in nEvts.iteritems():
        #print s, ": ", str(num)
    #return nEvts


### Evaluate the number of events in all the mass range
def evalNumEvtsAllRange(var, ch):
    ### Take histograms for bkg and signal
    bkg_hist = dict([ (s, ROOT.TFile.Open("NoNorm_"+ch+"/" + s+"EdmNtp.root").Get(var)) for s in bkgSamples])
    [h.Scale(scale[k]) for k, h in bkg_hist.iteritems()]
    sig_hist = dict([ (s, ROOT.TFile.Open("NoNorm_"+ch+"/" + s+"EdmNtp.root").Get(var)) for s in sigSamples])
    print sig_hist.values()
    [h.Scale(scale[k]) for k, h in sig_hist.iteritems()]
    
    ### Evaluate the number of events in the mass windows for signal and bakground
    ### for variable "var" (correspondig to a btag category) and for channel "ch"

    print "CATEGORY: ", var
    print "Channel: ", ch
    print " "
    nEvts = dict([(s,h.Integral() ) for s, h in bkg_hist.iteritems() ])
    nEvts_err = dict([(s,ROOT.TMath.Sqrt(h.Integral()) ) for s, h in bkg_hist.iteritems() ])
    
    nEvts["Background"] = nEvts["DYJetsToLL_M-50"] + nEvts["TT"]+ nEvts["ZZ"] + nEvts["WZ"] + nEvts["WW"]
    nEvts_err["Background"] = ROOT.TMath.Sqrt(nEvts["DYJetsToLL_M-50"] + nEvts["TT"]+ nEvts["ZZ"] + nEvts["WZ"] + nEvts["WW"])

    for s,h in sig_hist.iteritems():
        nEvts[s]=h.Integral()
        nEvts_err[s]=ROOT.TMath.Sqrt(h.Integral())

    nEvts_plus_err =  dict([(s, str.format("{0:.2f}", nEvts[s]) +"$\pm$"+ str.format("{0:.3f}", nEvts_err[s])  ) for s in nEvts.iterkeys() ])
    ## print "Event numbers"
    ## for (s, num) in nEvts.iteritems():
    ##     print s, ": ", str(num)
    ## print "*"*40

    #return nEvts, nEvts_err
    return nEvts_plus_err



### Evaluate the number of events in a mass window for all mass points
def tableEvtsMassWind(var, ch):

    if ch == "Mu": channel = "$2\mu2j$"
    else: channel = "$2e2j$"
    if var == "lljjmass": cat = "2 b-tag"
    elif var == "lljjmass_1btag": cat = "1 b-tag"
    else: cat = "0 b-tag"
    ### Take histograms for bkg and signal
    bkg_hist = dict([ (s, ROOT.TFile.Open("NoNorm_"+ch+"/" + s+"EdmNtp.root").Get(var)) for s in bkgSamples])
    [h.Scale(scale[k]) for k, h in bkg_hist.iteritems()]
    sig_hist = dict([ (s, ROOT.TFile.Open("NoNorm_"+ch+"/" + s+"EdmNtp.root").Get(var)) for s in sigSamples])
    [h.Scale(scale[k]) for k, h in sig_hist.iteritems()]
    
    ### Evaluate the number of events in the mass windows for signal and bakground
    ### for variable "var" (correspondig to a btag category) and for channel "ch"

    print "CATEGORY: ", var
    print "Channel: ", ch
    print " "
    [evalNumEvts(bkg_hist, sig_hist, m) for m in masses]

    print "*"*40

    hr = ['Mass [GeV]','Signal', '$Z$+Jets', '$tt$', 'ZZ/WZ/WW', 'TotalBkg']
    hc = [str(m) for m in masses]
    numbers = []
    for m in masses:
        signal = "GluGluToHToZZTo2L2Q_M-"+str(m)+"_8TeV"
        print signal
        nEvts = evalNumEvts(bkg_hist, sig_hist, m)
        print nEvts
        numbers.append([m,
                        nEvts[signal],
                        nEvts["DYJetsToLL_M-50"],
                        nEvts["TT"],
                        nEvts["ZZ"] + nEvts["WZ"] + nEvts["WW"],
                        nEvts["Background"]
                        ])


        t = matrix2latex(numbers, var+"_"+ch+"_evts", headerColumn=hc, headerRow=hr, format='%.2f', alignment = '|c|c|c|',caption='Expected yields of signal and background with  1~fb$^{-1}$ based on simulation in the '+cat+' category. The numbers show  '+channel+' expectations. Tighter $mZZ$ mass requirements are applied as  (-6\%, +\0\%) of the mass hypothesis.')


def tableEvtsAllRange(ch):
    if ch == "Mu": channel = "muon"
    else: channel = "electron"
    n_2btag = evalNumEvtsAllRange("lljjmass", ch)
    print n_2btag 
    n_1btag = evalNumEvtsAllRange("lljjmass_1btag", ch)
    n_0btag = evalNumEvtsAllRange("lljjmass_0btag", ch)
    #masses = [200, 210, 220, 230,250, 275, 275, 300, 325, 350, 375, 400, 450, 475, 500, 550, 575, 600]
    hr = [' ','0 b-tag yields', '1 b-tag yields','2 b-tag yields']
    hc = ["Background"]
    for m in masses:
        hc.append(str(m)+" GeV")
    #print hc
    numbers = []
    #str.format("{0:.3f}", pi) 
    numbers.append([ n_0btag["Background"],  n_1btag["Background"],  n_2btag["Background"] ])
    for m in masses:
        signal = "GluGluToHToZZTo2L2Q_M-"+str(m)+"_8TeV"
        numbers.append([ n_0btag[signal], n_1btag[signal], n_2btag[signal]])
        #print numbers
    #t = matrix2latex(numbers, "evtsAllRange_"+ch+".py",headerRow=hr)
    
    t = matrix2latex(numbers, "evtsAllRange_"+ch, headerColumn=hc, headerRow=hr, format='%.2f', alignment = '|c|c|c|',caption='List of expected background and signal yields in the '+channel+' channel with 1~fb$^{-1}$of data after all selection and within the $ZZ$ invariant mass range [0,1000]')    
    



### Call the function evalNumEvtsAllRange to get the number of events in the full range for signal and bkg
### evalNumEvtsAllRange(category, channel)
### category = "lljjmass", "lljjmass_1btag" or "lljjmass_0btag"
### channel = "El" or "Mu"


### Call the function evalNumEvtsMassWind to get the number of events in the mass range [-6%, +10%]mH for signal and bkg
### It evaluates numbers for all the mass points
### evalNumEvtsMassWind(category, channel)
### category = "lljjmass", "lljjmass_1btag" or "lljjmass_0btag"
### channel = "El" or "Mu"

channels = ["El", "Mu"]    
for ch in channels:
    tableEvtsAllRange(ch)

    tableEvtsMassWind("lljjmass", ch)
    tableEvtsMassWind("lljjmass_1btag", ch)
    tableEvtsMassWind("lljjmass_0btag", ch)


