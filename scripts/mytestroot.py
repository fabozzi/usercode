#!/usr/bin/python

import sys
#sys.argv.append('-b')
import os, commands
import math
import ROOT

#gROOT.ProcessLine(
#    "struct MyStruct {\
#    double run_number;\
#    double rpc_events;\
#    double bx_wm2;\
#    double bx_rms_wm2;\ 
#    double bx_wm1;\ 
#    double bx_rms_wm1;\ 
#    double bx_w0;\ 
#    double bx_rms_w0;\ 
#    double bx_wp1;\ 
#    double bx_rms_wp1;\ 
#    double bx_wp2;\ 
#    double bx_rms_wp2;\ 
#    double bx_dm3;\ 
#    double bx_rms_dm3;\ 
#    double bx_dm2;\ 
#    double bx_rms_dm2;\ 
#    double bx_dm1;\ 
#    double bx_rms_dm1;\ 
#    double bx_dp1;\ 
#    double bx_rms_dp1;\ 
#    double bx_dp2;\ 
#    double bx_rms_dp2;\ 
#    double bx_dp3;\ 
#    double bx_rms_dp3;\ 
#    double n_digis_em;\ 
#    double n_clusters_em;\ 
#    double cluster_size_em;\ 
#    double n_digis_barrel;\ 
#    double n_clusters_barrel;\ 
#    double cluster_size_barrel;\ 
#    double n_digis_ep;\ 
#    double n_clusters_ep;\ 
#    double cluster_size_ep;\ 
#} ;");

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


rfile_in = ROOT.TFile.Open("TestFile.root")

dbTree = rfile_in.Get("DBTree")

# get the branch variables
rpcefficiency = RpcEffStruct()
dbTree.SetBranchAddress("rpcefficiency", AddressOf(rpcefficiency, "run_number") )
#dbTree.SetBranchAddress("rpcefficiency", AddressOf(rpcefficiency, "eff_seg") )
for i in xrange(dbTree.GetEntries()):
    dbTree.GetEntry(i)
    if ( int(rpcefficiency.run_number) == 190679 ) :
        print rpcefficiency.raw_id, rpcefficiency.eff_seg

        
