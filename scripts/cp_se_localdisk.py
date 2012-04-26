import os, commands


#To run after  "source ~fabozzi/setdpm_host.csh cmsse02"

#dirnames =  ['MuRun2012A', 'DYJetsToLL_M-50']
#indirnames =  ['Run2012A_v3', 'DYJetsToLL_M-50_v1']

dirnames =  ['ElRun2012A']
indirnames =  ['ElRun2012A_v3']

#outpath = '/data3/scratch/users/decosa/Higgs/Summer12/'
outpath = '/data3/scratch/users/decosa/Higgs/TestSummer12/'
ntpdir = '/edmntp26Apr12/'


for a in dirnames:
    cmdLocalDisk = "mkdir "+ outpath 
    os.system(cmdLocalDisk)
    cmdLocalDisk = "mkdir "+ outpath + a 
    os.system(cmdLocalDisk)
    cmdLocalDisk = "mkdir "+ outpath + a +ntpdir
    os.system(cmdLocalDisk)   
    cmdPermission = "chmod 775 " + outpath + a +ntpdir
    os.system(cmdPermission)
    cmdClean = "rm "+ outpath + a +ntpdir + "*"
    os.system(cmdClean)
    print cmdLocalDisk
    print cmdPermission
    print  cmdClean
    cmd = "rfdir /dpm/na.infn.it/home/cms/store/user/decosa/" +indirnames[dirnames.index(a)]
    status,ls_la = commands.getstatusoutput( cmd )    
    files = []
    cmd = "mkdir " + a
    os.system(cmd)
    list = ls_la.split(os.linesep)
    
    for l in list:
        b = l.split()
        for c in b:
            if c.endswith('.root'):
                print "FILE"
                print c
                cmd = "rfcp /dpm/na.infn.it/home/cms/store/user/decosa/"+indirnames[dirnames.index(a)]+"/"+ c + " " +outpath +a+ntpdir
                os.system(cmd)
                print cmd
                files.append(c)
                
    print files
                   



