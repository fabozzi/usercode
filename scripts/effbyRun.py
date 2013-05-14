#!/usr/bin/python

import sys
#sys.argv.append('-b')
import os, commands
import math
import ROOT

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()        # don't pop up canvases


# open the official input run list from a file and make a vector from it
rpcRuns = []
runListFile = open("runslist.txt", "r")
for eachRun in runListFile:
    rpcRuns.append(int(eachRun.rstrip()))
ntotruns = len(rpcRuns) 

aveEffRuns = []
aveEffErrors = []
runErrors = []
for y in range(0, ntotruns):
    aveEffRuns.append(0.0)
    aveEffErrors.append(0.0)
    runErrors.append(0.0)
print rpcRuns
#print aveEffRuns

#print ntotruns

# This is the run number under study (actually one should loop over runs
runNum = 163269
aveEffInRun = 0.0 
aveEffErrorInRun = 0.0 

runNumIndex = rpcRuns.index( runNum )
#print runNumIndex

# plot eff for 1 or more chambers

# Specify the list of chambers to consider

# wheels = ['W+2', 'W+1', 'W+0', 'W-1', 'W-2', 'RE-3', 'RE-2', 'RE-1', 'RE+1', 'RE+2', 'RE+3']

# stations = ['RB1', 'RB2', 'RB3', 'RB4', 'R1', 'R2', 'R3']

# sectors = ['S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12',
#            'CH01', 'CH02', 'CH03', 'CH04', 'CH05', 'CH06', 'CH07', 'CH08', 'CH09', 'CH10', 'CH11', 'CH12',
#            'CH13', 'CH14', 'CH15', 'CH16', 'CH17', 'CH18', 'CH19', 'CH20', 'CH21', 'CH22', 'CH23', 'CH24',
#            'CH25', 'CH26', 'CH27', 'CH28', 'CH29', 'CH30', 'CH31', 'CH32', 'CH33', 'CH34', 'CH35', 'CH36']

# rolls = ['Forward', 'Middle', 'Backward', 'A', 'B', 'C']

wheels =['W-2']
stations =['RB1']
sectors =['S12']
rolls =['Forward']

# flags to bypass the above list and consider all wheels/stations/sectors/rolls
allwh = False
allstat = False
allsect = False
allrolls = True

# fill the list of bad chambers from input blacklist
badChambers = []
blackListFileName = "blacklist.txt"
blackListFile = open(blackListFileName, "r")
for chamber in blackListFile:
    chamber = (chamber.rstrip().split(' '))[1]
#    print chamber
    badChambers.append(chamber)
#print badChambers

#while 1:
#    chamber = blackListFile.readline().rstrip()
#    print chamber

# dictionaries to load eff/errors for selected chambers 
eff_map = {}
err_map = {}

# output file with eff histogram
rfile_out = ROOT.TFile.Open("output_eff.root","RECREATE")
#                wheelN = wheel.replace("+","p")
#                wheelNm = wheelN.replace("-","m")
#                stationN= station.replace("+","p") 
#                stationNm= stationN.replace("-","m") 
#hname = "eff_"+wheelNm+"_"+stationNm+"_"+sector+"_"+roll
hname = "eff"

h_effdistr = ROOT.TH1F(hname,hname,101, -0.5, 100.5) 
h_effruns = ROOT.TH1F("effruns","effruns",ntotruns, 0.5, ntotruns+0.5) 

#input fileName with chambers performances
inputEffFile = str(runNum)+"_rollYeff.txt"
infile = open(inputEffFile, "r")

# counter of good chambers considered for the efficiencies
ngoodChambers = 0
# looping on chambers in input file
while 1:
    line = infile.readline()
    newline = line.rstrip()    
    if newline.startswith("eff"):
        break;
#    print newline
    subpieces = newline.split(' ')
    chamberID = subpieces[0]
# check if it is a bad chamber: if so, skip the chamber in the loop    
    isBadChamber = False
    for badcham in badChambers :
        if chamberID == badcham :
            isBadChamber = True
            break
    if isBadChamber == True :
        continue
######
    idparts = chamberID.split('_')
    whID = idparts[0]
    statID = idparts[1]
    sectID = idparts[2]
    rollID = idparts[3]
    statID = statID.rstrip('inout+-')
# checking if that chamber in input file matches one of the chambers to study
    for wheel in wheels:
        for station in stations:
            for sector in sectors:
                for roll in rolls:
                    matchwh = allwh or (whID==wheel) 
                    matchstat = allstat or (statID==station) 
                    matchsect = allsect or (sectID==sector) 
                    matchroll = allrolls or (rollID==roll) 
#    print matchwh, matchstat, matchsect, matchroll
                    if (matchwh and matchstat and matchsect and matchroll):
                        effcham = float(subpieces[1])
                        errcham = float(subpieces[2])
                        eff_map[chamberID]=effcham
                        err_map[chamberID]=errcham
# filling histogram with efficiency for that chamber
                        h_effdistr.Fill(effcham)
                        aveEffInRun = aveEffInRun + effcham 
                        aveEffErrorInRun = aveEffErrorInRun + errcham * errcham 
                        ngoodChambers = ngoodChambers+1
    
print ngoodChambers
aveEffInRun = aveEffInRun / float(ngoodChambers)
aveEffErrorInRun = math.sqrt(aveEffErrorInRun) / float(ngoodChambers)
print aveEffInRun
print aveEffErrorInRun

infile.close()

aveEffRuns[runNumIndex] = aveEffInRun
aveEffErrors[runNumIndex] = aveEffErrorInRun

#print aveEffRuns
#print aveEffErrors
#print runErrors
# printout dictionaries
print eff_map
print err_map


# define graph with average eff. vs. run
#t_rpcRuns = ROOT.TVector(ntotruns)
#t_aveEffRuns = ROOT.TVector(ntotruns)
#t_runErrors = ROOT.TVector(ntotruns)
#t_aveEffErrors = ROOT.TVector(ntotruns)

for m in range(0, ntotruns) :
    h_effruns.SetBinContent(m+1,aveEffRuns[m])
    h_effruns.SetBinError(m+1,aveEffErrors[m])
    if m%50 == 0 :
        h_effruns.GetXaxis().SetBinLabel(m+1,str(rpcRuns[m]))
#    t_rpcRuns[m] = rpcRuns[m]
#    t_aveEffRuns[m] = aveEffRuns[m]
#    t_runErrors[m] = 0.0
#    t_aveEffErrors[m] = aveEffErrors[m]


#c_eff = ROOT.TCanvas("effruns", "effruns")
#myheff = c_eff.DrawFrame(1632680 0, 209160, 110)
#g_effbyrun = ROOT.TGraphErrors(t_rpcRuns, t_aveEffRuns, t_runErrors, t_aveEffErrors)
#g_effbyrun.SetName("effruns")
#g_effbyrun.SetTitle("effruns")


h_effdistr.Write()
h_effruns.Write()
#g_effbyrun.Write()
rfile_out.Close()    
