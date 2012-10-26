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


ch= "Mu"

inputpath = "/data3/scratch/users/decosa/Higgs/Summer12/"
ntpdir = "/edmntp12Jul12Cleaned/"

if(ch=="Mu"): dirnames = ["MuRun2012B", "MuRun2012A", "MuRun2012B_Ext"]
else: dirnames = ["ElRun2012B", "ElRun2012A", "ElRun2012B_Ext"]
files = []

for a in dirnames:
    cmd = "ls "+inputpath + a +ntpdir
    print cmd
    status,ls_la = commands.getstatusoutput( cmd )
    list_ = ls_la.split(os.linesep)
    for c in list_:
        files.append(inputpath + a + ntpdir + c)
#print files        
chain = ROOT.TChain("Events")
for f in files:
    chain.Add(f)


if(ch=="Mu"):
    scan = "muHiggsMass:muHiggsEventNumber:muHiggsLeptDau1Pt:muHiggsJetDau1puBeta:muHiggsJetDau2puBeta"
    events = ["417132403", "1219306590", "494765722"]
    for evt in events:
        chain.Scan("muHiggsRefitMass", "muHiggsEventNumber=="+str(evt))
else:
    scan = "elHiggsMass:elHiggsEventNumber:elHiggsLeptDau1Pt:elHiggsJetDau1puBeta:elHiggsJetDau2puBeta:elHiggsLeptDau1EtaSC:elHiggsLeptDau2EtaSC"
#    scan2 = "elHiggsLeptDau1TrkIso:elHiggsLeptDau2TrkIso:elHiggsLeptDau1EcalIso:elHiggsLeptDau2EcalIso:elHiggsLeptDau1HcalIso:elHiggsLeptDau2HcalIso"
    scan2 = "elHiggsEleDau1TrkIso03:elHiggssEleDau1EcalIso03:elHiggssEleDau1HcalIso03:elHiggsLeptDau1Pt"
    events = ["807387348", "27833017", "1161557958", "43723161"]
    
    for evt in events:
        cut = "elHiggsEventNumber=="+str(evt)+" && elHiggsRefitMass>260 && elHiggsRefitMass<261"
        chain.Scan(scan2, cut)
