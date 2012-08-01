
### Run this script after  "source ~fabozzi/setdpm_host.csh cmsse02"

import os, commands
import string
### Dirnames is the list of directories from which copying files
### indirnames is the list of directories to which copying files


inputpath = '/data3/scratch/users/fabozzi/Higgs/Summer12/'
ntpdir = '/edmntp26Jun12/'

outputpath = '/data3/scratch/users/decosa/Higgs/Summer12/'
ntpdirout = '/edmntp11Jul12Cleaned/'

#ntpdir = '/edmntpTest/'
#ntpdirout = '/edmntpCleanedTest/'
#outpath = '/data3/scratch/users/decosa/Higgs/TestSummer12/'
#ntpdir = '/edmntp23May12/'

#dirnames =  [ 'TT', 'WZ', 'ZZ']
#dirnames =  [ 'TEST']
#dirnames =  ['ElRun2012A', 'ElRun2012B']
dirnames = [ 'MuRun2012A','MuRun2012A_23May2012', 'MuRun2012B', 'MuRun2012B_Ext' ]
#dirnames = ['ElRun2012A','ElRun2012A_23May2012', 'ElRun2012B', 'ElRun2012B_Ext', 'MuRun2012A','MuRun2012A_23May2012', 'MuRun2012B', 'MuRun2012B_Ext' ]
#dirnames = ['GluGluToHToZZTo2L2Q_M-200_8TeV', 'GluGluToHToZZTo2L2Q_M-300_8TeV', 'GluGluToHToZZTo2L2Q_M-525_8TeV', 'GluGluToHToZZTo2L2Q_M-600_8TeV', 'GluGluToHToZZTo2L2Q_M-700_8TeV', 'GluGluToHToZZTo2L2Q_M-800_8TeV', 'GluGluToHToZZTo2L2Q_M-900_8TeV']
#dirnames =  ['MuRun2012A']
#dirnames =  ['DYJetsToLL_M-50', 'TT', 'WZ', 'ZZ']

#####dirnames =  ['DYJetsToLL_M-50']

#dirnames =  ['GluGluToHToZZTo2L2Q_M-300', 'GluGluToHToZZTo2L2Q_M-200']

for a in dirnames:

    ### Look at the original folder and check if there are more copies of the same ntuple
    # cmd = "rfdir /dpm/na.infn.it/home/cms/store/user/decosa/" +indirnames[dirnames.index(a)]
    #cmd = "lcg-ls -D srmv2 srm://srm.ciemat.es:8443/pnfs/ciemat.es/data/cms" +inputpaths[dirnames.index(a)]

    cmd = "ls "+inputpath + a +ntpdir
    print cmd
    status,ls_la = commands.getstatusoutput( cmd )    
    files = []
    #cmd = "mkdir " + a
    #os.system(cmd)
    list_ = ls_la.split(os.linesep)
    #print list
        
    for c in list_:
        #print c
        if (c.startswith('h2l2q_ntuple') and c.endswith('.root')):
            #print c
            #cmd =  "lcg-cp -D srmv2 srm://srm.ciemat.es:8443/pnfs/ciemat.es/data/cms" +inputpaths[dirnames.index(a)] + c + " " +outpath + a + ntpdir + c
            #print cmd
            #os.system(cmd)
            splitc = c.split("_")
            ntuple = "_".join(splitc[:3])
            check = [h.startswith(ntuple) for h in files]
            if (True not in check):
                files.append(c)
                    

    cfgfilename = "higgs2l2qNtupleFilter.py"    
    if(a.startswith("MuRun")): cfgfilename = "higgs2l2qNtupleFilter_Muons.py"
    elif(a.startswith("ElRun")): cfgfilename = "higgs2l2qNtupleFilter_Electrons.py"
    else : cfgfilename = "higgs2l2qNtupleFilter_MC.py"

    print "cfgTocopy"
    print  cfgfilename
    n = 20
    numfiles = len(files)/n + 1
    if(len(files)%n == 0 ): numfiles = len(files)/n

    print len(files)
    print n
    print numfiles
    outpath = outputpath + a
    outdir = outpath + ntpdirout

    cmdLocalDisk = "mkdir "+ outpath 
    os.system(cmdLocalDisk)
    cmdLocalDisk = "mkdir "+ outdir 
    os.system(cmdLocalDisk)
    cmdPermission = "chmod 775 " + outdir
    os.system(cmdPermission)
    cmdClean = "rm "+ outdir + "*"
    os.system(cmdClean)

    
    for i in range(numfiles):     
        num = str(i)
        file = open(cfgfilename, "r")
        newfile = open('higgs2l2qNtupleFilter_'+a+'_'+num+'.py', "w")
        print "new file: ", 'higgs2l2qNtupleFilter_'+a+'_'+num+'.py'
        selectedfiles = files[i*n:(n+i*n)]
        if (numfiles == 1):  selectedfiles = files[:]
        if ( (numfiles > 1) and  (i == (numfiles -1)) and (len(files)%n != 0)): selectedfiles = files[i*n:]
        
        while 1:
            line = file.readline()
            newLine = line
            if line == "":
                break;
            if line.startswith("process.source.fileNames=cms.untracked.vstring('file:h2l2q_ntuple.root')"):
                sources = [('file:'+inputpath + a +ntpdir + c) for c in selectedfiles]
                poolsource = '\' , \''.join(sources)
                #print poolsource
                newLine = string.replace(line, 'file:h2l2q_ntuple.root', poolsource );
            if (('h2l2q_ntuple_clean.root') in line):
                newLine = string.replace(line, 'h2l2q_ntuple_clean.root', outputpath + a +ntpdirout +'h2l2q_ntuple_clean_'+num+'.root' );
                #print newLine
            newfile.write(newLine)

    
        newfile.close()

        file.close()
                
    #print "FILES"            
    #print files
                   



