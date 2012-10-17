
import sys
sys.argv.append('-b')
import ROOT

from HiggsAnalysis.Higgs2l2b.scaleFactors import scale
from matrix2latex import matrix2latex

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')

ROOT.gROOT.SetBatch()

masses = [200, 210, 220, 230,250, 275, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600]
#masses = [450]
bkgSamples = [ "WW", "WZ", "ZZ", "TT", "DYJetsToLL_M-50"]
sigSamples = ["GluGluToHToZZTo2L2Q_M-"+str(m)+"_8TeV" for m in masses]
lumi = 5.1


ROOT.TH1.SetDefaultSumw2()

class num:
    pass



def integ(histo, mass):
    lower = 0.94*mass
    upper = 1.10*mass
    bmin = histo.GetXaxis().FindBin(lower)
    bmax = histo.GetXaxis().FindBin(upper)
    num.error = ROOT.Double(0.)
    if (mass == -1): num.integral = histo.IntegralAndError(0, 1000, num.error )
    else:
        num.integral = histo.IntegralAndError(bmin,bmax, num.error)
        num.integral -= (histo.GetBinContent(bmin)*(lower - histo.GetXaxis().GetBinLowEdge(bmin)))/ histo.GetXaxis().GetBinWidth(bmin)
        num.integral -= (histo.GetBinContent(bmax)*(histo.GetXaxis().GetBinUpEdge(bmax) -upper))/ histo.GetXaxis().GetBinWidth(bmax)
    return num

def fillEvtsList(ch, btag,  s, mH = -1):
    channel = 0.
    if (ch == "Mu"): channel = 1.
    f = ROOT.TFile("tree_NoNorm_"+ch+"/" + s+"_tree.root")
    tree = f.Get("lljjmassTree")
    histo = ROOT.TH1F("mH", "mH",1000, 0., 1000.)
    cut = "(mJJ>75 && mJJ<105 && lep==" + str(channel) + " && nBTags==" + str(btag) + ")"
    print cut
    tree.Draw("mZZ>> mH", "weight"+"*"+ cut)
    histo.Scale(lumi)
    num = integ(histo, mH)
    print num.integral
    print num.error
    return num 


def evalNumEvts(var, ch, mH = -1):
    btag = 0.
    if(var=="lljjmass_1btag"): btag = 1.
    elif(var=="lljjmass"): btag = 2.
    
    nEvts_err = dict([(s, fillEvtsList(ch, btag, s, mH).error ) for s in bkgSamples ])
    nEvts = dict([(s, fillEvtsList(ch, btag, s, mH).integral ) for s in bkgSamples ])
    nEvts["Background"] = nEvts["DYJetsToLL_M-50"] + nEvts["TT"]+ nEvts["ZZ"] + nEvts["WZ"] + nEvts["WW"]
    nEvts_err["Background"] = ROOT.TMath.Sqrt(
        ROOT.TMath.Power(nEvts_err["DYJetsToLL_M-50"], 2) +
        ROOT.TMath.Power(nEvts_err["TT"], 2)+
        ROOT.TMath.Power(nEvts_err["ZZ"], 2) +
        ROOT.TMath.Power(nEvts_err["WZ"], 2) +
        ROOT.TMath.Power(nEvts_err["WW"], 2)
        )
    nEvts["ZZ/WZ/WW"] = nEvts["ZZ"] + nEvts["WZ"] + nEvts["WW"]
    nEvts_err["ZZ/WZ/WW"] = ROOT.TMath.Sqrt(
        ROOT.TMath.Power(nEvts_err["ZZ"], 2) +
        ROOT.TMath.Power(nEvts_err["WZ"], 2) +
        ROOT.TMath.Power(nEvts_err["WW"], 2)
        )
        
    for s in sigSamples:
        nEvts_err[s] =  fillEvtsList(ch, btag, s, mH).error
        nEvts[s]=  fillEvtsList(ch, btag, s, mH).integral
        
    nEvts_plus_err =  dict([(s, str.format("{0:.2f}", nEvts[s]) +"$\pm$"+ str.format("{0:.2f}", nEvts_err[s])  ) for s in nEvts.iterkeys() ])
    print nEvts
    return nEvts_plus_err


### Evaluate the number of events in the full mass range
def tableEvtsAllRange(ch):
    if ch == "Mu": channel = "muon"
    else: channel = "electron"
    n_2btag = evalNumEvts("lljjmass", ch)
    print n_2btag 
    n_1btag = evalNumEvts("lljjmass_1btag", ch)
    n_0btag = evalNumEvts("lljjmass_0btag", ch)
    hr = [' ','0 b-tag yields', '1 b-tag yields','2 b-tag yields']
    hc = ["Background"]
    for m in masses:
        hc.append(str(m)+" GeV")
    numbers = []
    numbers.append([ n_0btag["Background"],  n_1btag["Background"],  n_2btag["Background"] ])
    for m in masses:
        signal = "GluGluToHToZZTo2L2Q_M-"+str(m)+"_8TeV"
        numbers.append([ n_0btag[signal], n_1btag[signal], n_2btag[signal]])
    
    t = matrix2latex(numbers, "evtsAllRange_"+ch, headerColumn=hc, headerRow=hr, format='%.2f', alignment = '|c|c|c|',caption='List of expected background and signal yields in the '+channel+' channel with 1~fb$^{-1}$of data after all selection and within the $ZZ$ invariant mass range [0,1000]')    
    



### Evaluate the number of events in a mass window for all mass points
def tableEvtsMassWind(var, ch):

    if ch == "Mu": channel = "$2\mu2j$"
    else: channel = "$2e2j$"
    if var == "lljjmass": cat = "2 b-tag"
    elif var == "lljjmass_1btag": cat = "1 b-tag"
    else: cat = "0 b-tag"
  
    hr = ['Mass [GeV]','Signal', '$Z$+Jets', '$tt$', 'ZZ/WZ/WW', 'TotalBkg']
    hc = [str(m) for m in masses]
    numbers = []
    for m in masses:
        nEvts = evalNumEvts(var, ch, m)
        signal = "GluGluToHToZZTo2L2Q_M-"+str(m)+"_8TeV"
        numbers.append([
                        nEvts[signal],
                        nEvts["DYJetsToLL_M-50"],
                        nEvts["TT"],
                        nEvts["ZZ/WZ/WW"],
                        nEvts["Background"]
                        ])


        t = matrix2latex(numbers, var+"_"+ch+"_evts", headerColumn=hc, headerRow=hr, format='%.2f', alignment = '|c|c|c|c|c|',caption='Expected yields of signal and background with  1~fb$^{-1}$ based on simulation in the '+cat+' category. The numbers show  '+channel+' expectations. Tighter $mZZ$ mass requirements are applied as  (-6\%, +\0\%) of the mass hypothesis.')


channels = ["El", "Mu"]    
for ch in channels:
    print "getting number of evts for all the cats"
    tableEvtsAllRange(ch)

    tableEvtsMassWind("lljjmass", ch)
    tableEvtsMassWind("lljjmass_1btag", ch)
    tableEvtsMassWind("lljjmass_0btag", ch)


