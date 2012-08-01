import sys
sys.argv.append('-b')

import ROOT
import os, commands

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')

ROOT.gROOT.SetBatch()        # don't pop up canvases

Vars = {}
#### Vars[variableName] = (rebin, xmin, xmax)

Vars["phi_eta"]=(  -3., 3.,-3., 3.,)
Vars["pt_phi"]=(  -3., 3.,0., 0.,)
#Vars["pt_eta"]=(  -3., 3.,30., 150.,)
#Vars["zllmass_eta"]=( -3., 3., 70., 110.)

#Vars["mjjmll"]=( 0., 300., 0., 300.)
#Vars["drJ1Lmjj"]= (0., 0., 0., 0.)
#Vars["drJ2Lmjj"]= (0., 0., 0., 0.)
#Vars["jzbmjj"]= (0.,350., -100., 150.)
#Vars["jzbmjjcorr"]= (0.,0., -200., 50.)
#Vars["zzdptmjj"]= (0.,350., -150., 150.)
#Vars["jptmjj"]= (0.,350., 0., 250.)
#Vars["trkMET_MET"]= (0.,0., 0., 0.)
#Vars["csvmjj"]=( 0., 0., 0., 0.)
#Vars["drzzmjj"]=( 0., 0., 0., 0.)
#Vars["ptjjmjj"]=( 0., 0., 0., 0.)
#Vars["drllmjj"]=( 0., 0., 0., 0.)
#Vars["dphijjmjj"]=( 0., 0., 0., 0.)
#Vars["detajjmjj"]=( 0., 0., 0., 0.)
#Vars["ptllptjj"]=( 0., 0., 0., 0.)

print Vars
#samples = [ "ZZ", "DYJets_Summ11"]
#samples = [ "MuRun2012A"]
samples = [ "ElRun2012A"]

#samples = [ "DYJetsToLL_M-50"]

#samples = [ "ZZ", "TT"]
#l = [ ROOT.kTeal - 5, ROOT.kOrange -2, ROOT.kBlue +2, ROOT.kBlue +1, ROOT.kBlue -4, ROOT.kGreen ]
l = [  ROOT.kGreen-3, ROOT.kRed ]
#l = [ ROOT.kRed , ROOT.kMagenta + 4]



rebinX = 1
rebinY = 2
rebinY2 = 2
ch = "El"
#ch = 'Mu'
#ch = 'All'


folder = "2Dplots_"+ch+"/"
if(samples[0]=="DYJetsToLL_M-50"): folder = "2Dplots_MC_"+ch+"/"
cmd = "rm "+folder+"*"
os.system(cmd)
cmd = "mkdir "+folder
os.system(cmd)
name = "phi_"

if(samples[0]=="DYJetsToLL_M-50"): l = [  ROOT.kBlue + 1 ]
    
c1 = ROOT.TCanvas("name","title", 5)
c1.SetTicks()

pad1 = ROOT.TPad()
ROOT.SetOwnership(pad1, False)
pad1.SetPad('pad1','pad1',0.0,0.5,1.0,1.0, 10)
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
pad2.SetPad('pad2','pad2',0.0,0.,1.0,0.5, 10)
pad2.Draw()
pad2.cd()
pad2.SetTopMargin(0.)
pad2.SetLeftMargin(0.1)
pad2.SetRightMargin(0.05)
pad2.SetBottomMargin(0.20)
pad2.SetTicks()

#var = "pt_eta"
var = "phi_eta"
xmin = Vars[var][0]
xmax = Vars[var][1]
ymin = Vars[var][2]
ymax = Vars[var][3]
            
leg = ROOT.TLegend(.75, .80, .98, .98)
leg.SetFillColor(0)
leg.SetTextSize(0.04)

leg2 = ROOT.TLegend(.75, .80, .98, .98)
leg2.SetFillColor(0)
leg2.SetTextSize(0.04)

           
all_hist = dict([ (s, ROOT.TFile.Open("NoNorm_"+ch+"/" + s+"EdmNtp.root" if ch !="All" else "NoNorm_El/" + s+"EdmNtp.root").Get(var)) for s in samples ])
print "NoNorm_",ch,"/", s,"EdmNtp.root"
print all_hist

[h.RebinX(rebinX) for h in all_hist.values()]
[h.RebinY(rebinY) for h in all_hist.values()]


var2 = "pt_phi"
#var2 = "zllmass_eta"
#xmin = Vars[var2][0]
#xmax = Vars[var2][1]
ymin2 = Vars[var2][2]
ymax2 = Vars[var2][3]

           
all_hist2 = dict([ (s, ROOT.TFile.Open("NoNorm_"+ch+"/" + s+"EdmNtp.root" if ch !="All" else "NoNorm_El/" + s+"EdmNtp.root").Get(var2)) for s in samples ])
print "NoNorm_",ch,"/", s,"EdmNtp.root"
print all_hist2

[h2.RebinX(rebinX) for h2 in all_hist2.values()]
[h2.RebinY(rebinY2) for h2 in all_hist2.values()]

# all_hist[samples[1]].Scale(0.01)
#all_hist[samples[0]].Scale(0.04)

i = 0
s = samples[i]
pad1.cd()
hist = all_hist[s]
hist.SetName(hist.GetName()+'_'+ch+'_'+s)
hist.SetStats(ROOT.kFALSE)
print hist.Integral()
#        hist.SetFillColor(l[i])
if(xmin!=xmax): hist.GetXaxis().SetRangeUser(xmin, xmax)
if(ymin!=ymax): hist.GetYaxis().SetRangeUser(ymin, ymax)
hist.SetMarkerColor(l[i])
hist.SetFillColor(l[i])
hist.SetMarkerSize(0.3)
hist.SetMarkerStyle(20)
hist.GetXaxis().SetTitle(var)
hist.Draw("box" if i==0 else "box same") 
leg.AddEntry(hist, s)
leg.Draw("SAME")

pad2.cd()
hist2 = all_hist2[s]
hist2.SetName(hist2.GetName()+'_'+ch+'_'+s)
hist2.SetStats(ROOT.kFALSE)
print hist.Integral()
#        hist.SetFillColor(l[i])
if(xmin!=xmax): hist2.GetXaxis().SetRangeUser(xmin, xmax)
if(ymin2!=ymax2): hist2.GetYaxis().SetRangeUser(ymin2, ymax2)
hist2.SetMarkerColor(l[i])
hist2.SetFillColor(l[i])
hist2.SetMarkerSize(0.3)
hist2.SetMarkerStyle(20)
hist2.GetXaxis().SetTitle(var2)
hist2.GetXaxis().SetLabelSize(0.07)
hist2.GetXaxis().SetTitleSize(0.08)
hist2.GetXaxis().SetTitle("#eta")
hist2.Draw("box" if i==0 else "box same") 
leg2.AddEntry(hist2, s)
leg2.Draw("SAME")


    # histogram style setting     
    
c1.Print(folder + name+"_"+ch+"_"+samples[0]+".eps")
c1.Print(folder + name+"_"+ch+"_"+samples[0]+".png")





# make plots for all variables in vars
#map(makePlot, Vars.iterkeys())
