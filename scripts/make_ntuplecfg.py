import os, commands
import string

# script which divides collections listed in a samples.txt file obtained
# with a rfdir command like this:
# rfdir /dpm/na.infn.it/home/cms/store/user/fabozzi/Summer11MCSkim/WJetsToLNu > samples.txt
# the collections are divided in sets of size = 100
# for each collection set, a cfg file is created from a template
# and a job is sent
# ADAPT THE SCRIPT TO YOUR NEED

path_storage = "\'rfio:/dpm/na.infn.it/home/cms/store/user/fabozzi/Summer11MCSkim/WJetsToLNu/"

input_samplefile = open("samples.txt", "r")
samplebunches = [ ]
samplelist = [ ]

listedfiles = input_samplefile.readlines()
numsamples = len(listedfiles)
i = 0
j = 0

for tempfile in listedfiles: 
    dd = tempfile.split()
    for ddd in dd:
        if ddd.endswith('.root') :
            samplelist.append(path_storage+ddd+"\',\n")
#            print i, j
            j += 1
            i += 1
            if (j==100) or (i==numsamples) :
#                print len(samplelist)
                tempsample = samplelist[:]
                samplebunches.append(tempsample)
                del samplelist[0:j]
                j=0


#i = 0
#for sample in samplebunches :
#    print i
#    print sample
#    i += 1

i = 0
for samples in samplebunches :    
    cfgtemplate = open("higgs2l2bNtupleMC_42X.py", "r")
    cfgfile = open("higgs2l2bNtupleMC_42X_"+str(i)+".py", "w")
    writelines = cfgtemplate.readlines()
    for commandline in writelines :
        if commandline.endswith('h2l2bData.root\'\n') :
            for sample in samples :
 #               print sample
                cfgfile.write(sample)
        elif commandline.endswith('h2l2b450_histo.root\') )\n') :
            cfgfile.write(commandline.replace('h2l2b450_histo','h2l2b450_histo_'+str(i) ) )            
        elif commandline.endswith('h2l2b_ntuple.root\')\n') :
            cfgfile.write(commandline.replace('h2l2b_ntuple','h2l2b_ntuple_'+str(i) ) )                        
        else:
            cfgfile.write(commandline)
    cfgtemplate.close()
    cfgfile.close()
    os.system("cmsRun higgs2l2bNtupleMC_42X_"+str(i)+".py >& log_"+str(i)+".txt&");
    i += 1

