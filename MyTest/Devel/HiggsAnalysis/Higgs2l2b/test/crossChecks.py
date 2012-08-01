
ch = "El"
filename = "SelectedEvents"+ch+".txt"
f = open(filename, "r")
files = f.readlines()

h = [s.split() for s in files]

runNums_0btag = [k[1] for k in h if k[len(k)-1]=='0b']
runNums_1btag = [k[1] for k in h if k[len(k)-1]=='1b']
runNums_2btag = [k[1] for k in h if k[len(k)-1]=='2b']


f_0b = open("SelectedEvts"+ch+"_0btag.txt", "w")
[f_0b.write(newLine+"\n") for newLine in runNums_0btag]
f_1b = open("SelectedEvts"+ch+"_1btag.txt", "w")
[f_1b.write(newLine+"\n") for newLine in runNums_1btag]
f_2b = open("SelectedEvts"+ch+"_2btag.txt", "w")
[f_2b.write(newLine+"\n") for newLine in runNums_2btag]
