#!/usr/bin/python

import sys
#sys.argv.append('-b')
import os, commands
import math
import ROOT

from ROOT import *
gROOT.ProcessLine(
    "struct RpcEffStruct {\
    double run_number;\
    double raw_id;\
    double eff_seg;\
    double eff_seg_error;\
    double n_extrap;\
    double cluster_size;\
    double clus_size_bin01;\
    double clus_size_bin02;\
    double clus_size_bin03;\
    double clus_size_bin04;\
    double clus_size_bin05;\
    double clus_size_bin06;\
    double clus_size_bin07;\
    double clus_size_bin08;\
    double clus_size_bin09;\
    double clus_size_bin10;\
} ;");

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()        # don't pop up canvases


# open chamberID definition file and create a dictionary
chambDict = {}
chambCounter = 0
chambDictFile = open("mapRoll.ascii", "r")
for chambEntry in chambDictFile:
#    print chambEntry
    chambCounter = chambCounter + 1
    chambID = chambEntry.rstrip().split(' ')
    chambDict[ chambID[1] ] = chambID[0]

#print chambCounter
#print len(chambDict)
#print chambDict
############################################################

# flag to remove badchambers 
removeBadRolls = True
############################################################

# fill the list of bad chambers from input blacklist
badChambers = []
blackListFileName = "blacklist.txt"
blackListFile = open(blackListFileName, "r")
for chamber in blackListFile:
    chamber = (chamber.rstrip().split(' '))[1]
#    print chamber
    badChambers.append(chamber)
#print badChambers
########################################################

# open the official input run list from a file and make a vector from it
rpcRuns = []
runListFile = open("runslist.txt", "r")
for eachRun in runListFile:
    rpcRuns.append(int(eachRun.rstrip()))
ntotruns = len(rpcRuns) 
############################################################

# define vectors to store average efficiencies and errors
aveEffRuns = []
aveEffErrors = []

for y in range(0, ntotruns):
    aveEffRuns.append(0.0)
    aveEffErrors.append(0.0)
############################################################
#print rpcRuns
#print aveEffRuns
#print ntotruns

firstRun = 190679
lastRun = 194076

runNumIndexFirst = rpcRuns.index( firstRun )
runNumIndexLast = rpcRuns.index( lastRun )




# Specify the list of chambers to consider

# wheels = ['W+2', 'W+1', 'W+0', 'W-1', 'W-2', 'RE-3', 'RE-2', 'RE-1', 'RE+1', 'RE+2', 'RE+3']

# stations = ['RB1', 'RB2', 'RB3', 'RB4', 'R1', 'R2', 'R3']

# sectors = ['S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12',
#            'CH01', 'CH02', 'CH03', 'CH04', 'CH05', 'CH06', 'CH07', 'CH08', 'CH09', 'CH10', 'CH11', 'CH12',
#            'CH13', 'CH14', 'CH15', 'CH16', 'CH17', 'CH18', 'CH19', 'CH20', 'CH21', 'CH22', 'CH23', 'CH24',
#            'CH25', 'CH26', 'CH27', 'CH28', 'CH29', 'CH30', 'CH31', 'CH32', 'CH33', 'CH34', 'CH35', 'CH36']

# rolls = ['Forward', 'Middle', 'Backward', 'A', 'B', 'C']

# all barrel
wheels = ['W+2', 'W+1', 'W+0', 'W-1', 'W-2']
stations = ['RB1', 'RB2', 'RB3', 'RB4']
stat_postfixes = ['', 'in', 'out', '+', '-', '++', '--']
sectors = ['S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12']
rolls = ['Forward','Middle','Backward']

fileoutpostfix = 'barrel'

# all endcap
#wheels = ['RE-3', 'RE-2', 'RE-1', 'RE+1', 'RE+2', 'RE+3']
#stations = ['R1', 'R2', 'R3']
#stat_postfixes = ['', 'in', 'out', '+', '-', '++', '--']
#sectors = ['CH01', 'CH02', 'CH03', 'CH04', 'CH05', 'CH06', 'CH07', 'CH08', 'CH09', 'CH10', 'CH11', 'CH12',
#           'CH13', 'CH14', 'CH15', 'CH16', 'CH17', 'CH18', 'CH19', 'CH20', 'CH21', 'CH22', 'CH23', 'CH24',
#           'CH25', 'CH26', 'CH27', 'CH28', 'CH29', 'CH30', 'CH31', 'CH32', 'CH33', 'CH34', 'CH35', 'CH36']
#rolls = ['A', 'B', 'C']
#fileoutpostfix = 'endcap'


# vector to store rawID of chambers under study
# it will be used to get efficiencies
rawIds = []

# loop on the requested chambers
for wheel in wheels:
    for station in stations:
        for stat_postfix in stat_postfixes:
            for sector in sectors:
                for roll in rolls:
                    chambName = wheel+'_'+station+stat_postfix+'_'+sector+'_'+roll
#                    print chambName
                    # check if the chamber name exists                     
                    if (chambName in chambDict):
#                        print 'Chamber exists!'
                        # check if it is a bad chamber: if so, skip the chamber in the loop    
                        isBadChamber = False
                        if removeBadRolls :
                            for badcham in badChambers :
                                if chambName == badcham :
                                    isBadChamber = True
#                                    print 'BAD CHAMBER!!!'
                                    break
                        if not(isBadChamber) :
                            # get the Raw ID of the chamber
                            chambRawID = chambDict[chambName]
#                            print chambRawID
                            rawIds.append(chambRawID)


#print rawIds

outfilename = 'output_eff_'+fileoutpostfix+'.root'


rfile_in = ROOT.TFile.Open("TestFile.root")

dbTree = rfile_in.Get("DBTree")
# get the branch variables
rpcefficiency = RpcEffStruct()
dbTree.SetBranchAddress("rpcefficiency", AddressOf(rpcefficiency, "run_number") )

rfile_out = ROOT.TFile.Open(outfilename,"RECREATE")
h_effruns = ROOT.TH1F("effruns","effruns",ntotruns, 0.5, ntotruns+0.5) 



for myRunInd in range(runNumIndexFirst, runNumIndexLast+1) :
    runNum = rpcRuns[ myRunInd ]
    runNumIndex = myRunInd
    print 'analyzing run ', runNum
    h_effdistr = ROOT.TH1F("eff"+str(runNum),"eff"+str(runNum),1600, 0, 100) 

# This is the run number under study (actually one should loop over runs
#runNum = 163269
#runNum = 190679
#runNumIndex = rpcRuns.index( runNum )
#print runNumIndex

# temporary variables
    aveEffInRun = 0.0 
    aveEffErrorInRun = 0.0 

# plot eff for 1 or more chambers

# now loop on the TTree to manipulate efficiencies

    ngoodChambers = 0

    for i in xrange(dbTree.GetEntries()):
        dbTree.GetEntry(i)
        #    print runNum
        #    print int(rpcefficiency.run_number)
        if ( int(rpcefficiency.run_number) == runNum ) :
            for id in rawIds :
                if (id == str(int(rpcefficiency.raw_id))) :
                    #                print id, int(rpcefficiency.raw_id) 
                    #                print rpcefficiency.raw_id, rpcefficiency.eff_seg
                    effcham = rpcefficiency.eff_seg
                    errcham = rpcefficiency.eff_seg_error
                    if errcham > 0 :
                        h_effdistr.Fill(effcham) 
                        aveEffInRun = aveEffInRun + effcham 
                        aveEffErrorInRun = aveEffErrorInRun + errcham * errcham 
                        ngoodChambers = ngoodChambers+1                
                        
                        
                        #    

    aveEffInRun = aveEffInRun / float(ngoodChambers)
    aveEffErrorInRun = math.sqrt(aveEffErrorInRun) / float(ngoodChambers)
    print aveEffInRun
    print aveEffErrorInRun
    
    
    aveEffRuns[runNumIndex] = aveEffInRun
    aveEffErrors[runNumIndex] = aveEffErrorInRun
    h_effdistr.Write()
    

print aveEffRuns
print aveEffErrors
                        
for m in range(0, ntotruns) :
    h_effruns.SetBinContent(m+1,aveEffRuns[m])
    h_effruns.SetBinError(m+1,aveEffErrors[m])
    if m%50 == 0 :
        h_effruns.GetXaxis().SetBinLabel(m+1,str(rpcRuns[m]))


#h_effdistr.Write()
h_effruns.Write()
#

rfile_out.Close()    





