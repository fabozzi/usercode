import sys
sys.argv.append('-b')
import ROOT
import os, commands
from variables import Vars

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')

ROOT.gROOT.SetBatch()        # don't pop up canvases

#samples = ["ZZ", "300"]
#samples = [ "ZZ", "DYJets_Summ11", "300", "400", "500"]
samples = [ "TT", "ZZ", "DYJets_Summ11"]

#l = [ ROOT.kTeal - 5, ROOT.kOrange -2, ROOT.kBlue +2, ROOT.kBlue +1, ROOT.kBlue -4, ROOT.kGreen ]
l = [ ROOT.kBlue +2,  ROOT.kRed , ROOT.kGreen +1, ROOT.kBlue +2, ROOT.kBlue +1, ROOT.kBlue -4, ROOT.kGreen ]
scale = { '300':0.00073510 , 'ZZ':0.002798, '400':0.0005031 , '500':0.00021778, 'DYJets_Summ11':0.0840179, 'TT':0.0165 }

intLumi = 4640.6
#intLumi = 20000.

lumiFact = intLumi/1000
for k, v in scale.iteritems():
    scale[k] = v * lumiFact

folder = "zzanalysis/"
cmd = "mkdir "+folder
os.system(cmd)

ch = "El"
# ch = 'Mu'
#ch = 'All'

def makePlot(varKey):

    var = varKey
    rebin = Vars[varKey][0]
    xmin = Vars[varKey][1]
    xmax = Vars[varKey][2]
    logScale = Vars[varKey][3]
    c1 = ROOT.TCanvas()
    m = 0.
    leg = ROOT.TLegend(.72, .72, .98, .98)
    leg.SetFillColor(0)
    leg.SetTextSize(0.04)

    # create a dictionary of histograms for variable var for each sample


    all_hist = dict([ (s, ROOT.TFile.Open("NoNorm_"+ch+"/" + s+"EdmNtp.root" if ch !="All" else "NoNorm_El/" + s+"EdmNtp.root").Get(var)) for s in scale.keys() ])
    dataRun2011B = ROOT.TFile.Open("NoNorm_"+ch+"/"+ch+"Run2011BEdmNtp.root"if ch !="All" else "NoNorm_El/ElRun2011BEdmNtp.root").Get(var)
    dataRun2011A = ROOT.TFile.Open("NoNorm_"+ch+"/"+ch+"Run2011AEdmNtp.root"if ch !="All" else "NoNorm_El/ElRun2011AEdmNtp.root").Get(var)
    #    dataRun2011A = ROOT.TFile.Open("NoNorm_"+ch+"/"+ch+"Run2011A_Oct03EdmNtp.root"if ch !="All" else "NoNorm_El/ElRun2011A_Oct03EdmNtp.root").Get(var)
    dataRun2011A.Add(dataRun2011B)
    dataRun2011A.Rebin(rebin)
    ROOT.SetOwnership(dataRun2011A,True)
    # rebin and scale histograms
    #    [h.Sumw2() for h in all_hist.itervalues()]    
    [h.Rebin(rebin) for h in all_hist.values()]
    ### uncomment the following line to normalize at lumi
    [h.Scale(scale[k]) for k, h in all_hist.iteritems()]
    if (ch == 'All'):
        all_hist_2 = dict([ (s, ROOT.TFile.Open("NoNorm_Mu/" + s+"EdmNtp.root").Get(var)) for s in scale.keys() ])
        [h.Rebin(rebin) for h in all_hist_2.values()]
        [h.Scale(scale[k]) for k, h in all_hist_2.iteritems()]
        [h.Add(h2) for h,h2 in zip(all_hist.itervalues(), all_hist_2.itervalues()) ]
     
        #    [h.Scale(1./h.Integral()) if h.Integral()!=0 else h.Scale(1)  for  h in all_hist.itervalues()]
    maximum = max([ h.GetMaximum() for h in all_hist.itervalues()])
    minimum = min([ h.GetMinimum() for h in all_hist.itervalues()])
    
    # histogram style setting     
    for i, s in enumerate(samples):
        hist = all_hist[s]
        hist.SetName(hist.GetName()+'_'+ch+'_'+s)
        if(xmin!=xmax): hist.GetXaxis().SetRangeUser(xmin, xmax)
        if(logScale): hist.SetMinimum(0.001)
        if(logScale and minimum!=0): hist.SetMinimum(minimum)
        hist.SetStats(ROOT.kFALSE)
        #        hist.SetFillColor(l[i])
        hist.SetLineColor(l[i])
        hist.SetLineWidth(4)
        hist.SetMaximum(maximum*1.5)
        hist.GetXaxis().SetTitle(var)
        hist.Draw("HIST" if i==0 else "HISTSAME") 
        leg.AddEntry(hist, s, "f")

    dataRun2011A.Draw("SAME")
    leg.AddEntry(dataRun2011A, "DATA", "*")
    leg.Draw("SAME")
    if(logScale): c1.SetLogy()
    c1.Print(folder + var+"_"+ch+"_data.eps")
    c1.Print(folder + var+"_"+ch+"_data.png")



# make plots for all variables in vars
map(makePlot, Vars.keys())
