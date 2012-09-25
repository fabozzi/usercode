import os, commands
#from scaleFactors import scale

from HiggsAnalysis.Higgs2l2b.scaleFactors import scale

### ***************************
### *** BEGIN CONFIGURATION ***

### Indicate :
###     - channel on which you want to run the analysis -> channels
###     - the run period: 2011A, 2011B, 2011All, 2012Al -> runPeriod
###     - input edmNtuple folder -> ntupleFolder
###     - samples on which you want to run -> dirnames

channels = ['Mu', 'El']
#channels = ['El']
#channels = ['Mu']

runPeriod = "2012All"

ntupleFolder = "edmntp12Jul12Cleaned"

signalMap = {
   'GluGluToHToZZTo2L2Q_M-200_8TeV': 200,
   'GluGluToHToZZTo2L2Q_M-210_8TeV': 210,
   'GluGluToHToZZTo2L2Q_M-250_8TeV': 250,
   'GluGluToHToZZTo2L2Q_M-275_8TeV': 275,
   'GluGluToHToZZTo2L2Q_M-300_8TeV': 300,
   'GluGluToHToZZTo2L2Q_M-325_8TeV': 325,
   'GluGluToHToZZTo2L2Q_M-350_8TeV': 350,
   'GluGluToHToZZTo2L2Q_M-400_8TeV': 400,
   'GluGluToHToZZTo2L2Q_M-425_8TeV': 425,
   'GluGluToHToZZTo2L2Q_M-500_8TeV': 500,
   'GluGluToHToZZTo2L2Q_M-600_8TeV': 600,
   }

#ntupleFolder = "edmntpTestCleaned"
#ntupleFolder2 = "edmntp26Jun12"
#ntupleFolder2 = "edmntpTest"
#dirname = ['GluGluToHToZZTo2L2Q_M-300_8TeV', 'GluGluToHToZZTo2L2Q_M-200_8TeV','GluGluToHToZZTo2L2Q_M-525_8TeV', 'GluGluToHToZZTo2L2Q_M-600_8TeV','GluGluToHToZZTo2L2Q_M-700_8TeV', 'GluGluToHToZZTo2L2Q_M-800_8TeV', 'GluGluToHToZZTo2L2Q_M-900_8TeV']

#dirname = [  'GluGluToHToZZTo2L2Q_M-300_8TeV', 'GluGluToHToZZTo2L2Q_M-200_8TeV',  'GluGluToHToZZTo2L2Q_M-400_8TeV', 'GluGluToHToZZTo2L2Q_M-600_8TeV']
dirname = [  'GluGluToHToZZTo2L2Q_M-300_8TeV', 'GluGluToHToZZTo2L2Q_M-200_8TeV', 'GluGluToHToZZTo2L2Q_M-400_8TeV', 'GluGluToHToZZTo2L2Q_M-600_8TeV', 'GluGluToHToZZTo2L2Q_M-210_8TeV', 'GluGluToHToZZTo2L2Q_M-250_8TeV', 'GluGluToHToZZTo2L2Q_M-275_8TeV', 'GluGluToHToZZTo2L2Q_M-325_8TeV', 'GluGluToHToZZTo2L2Q_M-350_8TeV', 'GluGluToHToZZTo2L2Q_M-425_8TeV']
#dirname = [ 'DYJetsToLL_M-50']
#dirname = [  'MuRun2012A', 'MuRun2012A_23May2012', 'MuRun2012B', 'MuRun2012B_Ext', 'ElRun2012B', 'ElRun2012B_Ext','ElRun2012A', 'ElRun2012A_23May2012']
#dirname = ['MuRun2012B_Ext']
#dirname = ['TT', 'ZZ', 'WZ']
#dirname = [ 'ElRun2012B', 'ElRun2012B_Ext','ElRun2012A', 'ElRun2012A_23May2012']

#dirname = ['ElRun2012A']
### *** END CONFIGURATION ***
### *************************
   
pu_per = ""
   
if (runPeriod == "2011A") :
   pu_per = "2011A"
if (runPeriod == "2011B") :
   pu_per = "2011B"
if(runPeriod == "2012A") :
   pu_per = "2012A"
if(runPeriod == "2012All") :
   pu_per = "2012All"

   
print "PU period: "
print pu_per
sep = " "

#for selection in selections:
for channel in channels:
   cmdDir = "rm SelectedEvents"+channel+".txt"
   os.system(cmdDir)
   cmdDir = "mkdir NoNorm_" + channel
   os.system(cmdDir)
   cmdDir = "mkdir tree_NoNorm_" + channel
   os.system(cmdDir)
   
   
   for a in dirname:
      data = "MC"
      sf = "SF"
      wt = "tree"
      txt = "txt"
      path = "/data3/scratch/users/decosa/Higgs/Summer12/"+a+"/"+ ntupleFolder+ "/"
#      applyfix = "NOfixMu"
      scaleFact = 1
      applyLR = "noLR"
      mH = 600

      if( a.startswith("GluGlu") ):
         mH = signalMap[ a ]
         if( mH >= 400 ):
            applyLR = "LR"
      
      if( (not a.startswith("ElRun")) and (not  a.startswith("MuRun"))): scaleFact = scale[a]
      
      if( (channel == "Mu") and (a.startswith("ElRun")) ) :
         continue
      if( (channel == "El") and (a.startswith("MuRun")) ) :
         continue      
      if ( a.startswith("MuRun") or a.startswith("ElRun") ) :
         #path = "/data3/scratch/users/decosa/Higgs/Summer12/"+a+"/"+ ntupleFolder2+ "/"
         data= "DATA"
         sf = "noSF"
         wt = "tree"
         txt = "txt"
         scaleFact = 1

      print path
      cfgfilename = "test.py"
      
      cmd = "ls"+" "+path
      status,ls_la = commands.getstatusoutput( cmd )
      
      list = ls_la.split(os.linesep)
      dir = []
      
      merge = "mergeTFileServiceHistograms -i "
      cmd = "rm tree_NoNorm_"+ channel + "/"+a+"_tree.root"
      os.system(cmd)
      mergeTree = "hadd "+"tree_NoNorm_" + channel + "/"+a+"_tree.root"
      cmd = "rm "+a+"_tree.root"
      os.system(cmd)
      for d in list:
         if d.endswith('.root'):
            
            cmd="2l2b_2012 "+path+d+" "+"NoNorm_" + channel + "/"+ a + str(list.index(d))+".root" +sep+data+sep+sf+sep+channel+sep+wt+sep+txt+sep+pu_per+sep+str(scaleFact)+sep+applyLR+sep+str(mH)
            ### To run macro for ZZ analysis, please, uncomment the following line and comment the previous one.
            ### Electron ID in ZZanalysis is not up to date to 2l2j analysis - 16 April 2012 
            #cmd="zz "+path+d+" "+"NoNorm_" + channel + "/"+ a + str(list.index(d))+".root" +sep+data+sep+sf+sep+channel+sep+wt+sep+txt+sep+pu_per+sep+applyfix
            print cmd
            merge = merge + " "+"NoNorm_" + channel+ "/"+a+str(list.index(d))+".root"
            mergeTree = mergeTree + " "+"tree_NoNorm_" + channel+ "/"+a+str(list.index(d))+".root"
            os.system(cmd)

      merge = merge +" -o "+"NoNorm_"+ channel+ "/"+a+"EdmNtp.root"
      
      os.system(merge)
   
      if wt == "tree":
         print "MERGE TREE COMMAND "
         print mergeTree
         os.system(mergeTree)
      for d in list:
         if d.endswith('.root'):
            cmd ="rm "+ "NoNorm_"+ channel + "/"+a+str(list.index(d))+".root"
            os.system(cmd)
            if wt=="tree" :
               cmd ="rm "+ "tree_NoNorm_"+ channel + "/"+a+str(list.index(d))+".root"
               os.system(cmd)
            
if (runPeriod == "2011All") :
   if ( ("MuRun2011A" in dirname) and ("MuRun2011B" in dirname) ):
      cmd = "rm NoNorm_Mu/MuRun2011AllEdmNtp.root"
      os.system(cmd)
      cmd = "mergeTFileServiceHistograms -i NoNorm_Mu/MuRun2011*EdmNtp.root -o NoNorm_Mu/MuRun2011AllEdmNtp.root"
      os.system(cmd)
   if ( ("ElRun2011A" in dirname) and ("ElRun2011B" in dirname) ):
      cmd = "rm NoNorm_El/ElRun2011AllEdmNtp.root"
      os.system(cmd)
      cmd = "mergeTFileServiceHistograms -i NoNorm_El/ElRun2011*EdmNtp.root -o NoNorm_El/ElRun2011AllEdmNtp.root"
      os.system(cmd)


if (runPeriod == "2012All") :
   if ( ("MuRun2012A" in dirname) and ("MuRun2012B" in dirname) ):
      cmd = "rm NoNorm_Mu/MuRun2012AllEdmNtp.root"
      os.system(cmd)
      cmd = "mergeTFileServiceHistograms -i NoNorm_Mu/MuRun2012*EdmNtp.root -o NoNorm_Mu/MuRun2012AllEdmNtp.root"
      os.system(cmd)
      cmd = "rm tree_NoNorm_Mu/MuRun2012All_tree.root"
      os.system(cmd)
      cmd = "hadd tree_NoNorm_Mu/MuRun2012All_tree.root tree_NoNorm_Mu/MuRun2012*_tree.root"
      os.system(cmd)
   if ( ("ElRun2012A" in dirname) and ("ElRun2012B" in dirname) ):
      cmd = "rm NoNorm_El/ElRun2012AllEdmNtp.root"
      os.system(cmd)
      cmd = "mergeTFileServiceHistograms -i NoNorm_El/ElRun2012*EdmNtp.root -o NoNorm_El/ElRun2012AllEdmNtp.root"
      os.system(cmd)
      cmd = "rm tree_NoNorm_El/ElRun2012All_tree.root"
      os.system(cmd)
      cmd = "hadd tree_NoNorm_El/ElRun2012All_tree.root tree_NoNorm_El/ElRun2012*_tree.root"
      os.system(cmd)
