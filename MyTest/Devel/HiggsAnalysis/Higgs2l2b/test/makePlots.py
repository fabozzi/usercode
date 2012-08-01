import os, commands

os.system("mkdir plots")
os.system("mkdir histos_Mu")
os.system("mkdir histos_El")
os.system("mkdir reports")

# plot.C parameters
###################################

# plot(string runPeriod= "A", string var = "zllmass", string label="m_{ll}", int rebin_fact = 1, float ymax = 50., bool set_logScale=false, bool custom = false, bool drawSig = false, string masspoint="500", bool drawSigStack = false, float lowcut = 70., float upcut = 110.)
# runPeriod = runPeriod for 2011: A, B, All
# var = name of the histo variable to plot
# label = name of the variable label for stack plot
# rebin_fact = rebinning factor of histogram to plot
# ymax = offset to ymax for stack plot
# set_logScale = to draw in log scale
# custom = setting custom limits to zllmass histogram
# drawSig = superimpose signal distribution (eventually amplified by a factor)
# masspoint = mass of signal sample to use in the plot
# drawSigStack = add signal to stack plot
# lowcut = lower limit for integration (apply for variables != lljjmass)
# upcut = upper limit for integration (apply for variables != lljjmass)

drawplotsCmd = [
# after kine/id lepton cuts
    "root -l -b \'mergedPlot.C(\"All\", \"npv\", \"nPV\", 1, 500., false, true, false, \"500\")\'  ",
    "root -l -b \'mergedPlot.C(\"All\", \"npv1\", \"nPV\", 1, 500., false, false, false, \"500\")\'  ",
    "root -l -b \'mergedPlot.C(\"All\", \"npv2\", \"nPV\", 1, 500., false, false, false, \"500\")\'  ",
    "root -l -b \'mergedPlot.C(\"All\", \"npvfin\", \"nPV\", 1, 10., false, false, false, \"500\")\'  ",
    "root -l -b \'plot.C(\"All\", \"npvfin\", \"nPV\", 1, 10., false, false, false, \"500\", false, 0, 30)\'  ",
    "root -l -b \'plot.C(\"All\", \"npv\", \"nPV\", 1, 500., false, true, false, \"500\", false, 0, 30)\'  ",

    "root -l -b \'plot.C(\"All\", \"zllmass\", \"m_{ll}\", 1, 200., false, true, false, \"500\", false, 70, 110 )\'  ",
    "root -l -b \'plot.C(\"All\", \"zjjmass\", \"m_{jj}\", 4, 100., false, false, false, \"500\", false, 75, 105 )\'  ",
    "root -l -b \'plot.C(\"All\", \"zllmass\", \"m_{ll}\", 1, 200., true, true, false, \"500\", false, 70, 110 )\'  ",
    "root -l -b \'plot.C(\"All\", \"zjjmass\", \"m_{jj}\", 4, 100., true, false, false, \"500\", false, 75, 105 )\'  ",
    "root -l -b \'plot.C(\"All\", \"leptpt1\", \"p_{T}\", 4, 400., true, false, false, \"500\", false)\'  ",
    "root -l -b \'plot.C(\"All\", \"leptpt2\", \"p_{T}\", 4, 400., true, false, false, \"500\", false)\'  ",
    "root -l -b \'plot.C(\"All\", \"jetpt1\", \"p_{T}\", 4, 400., true, false, false, \"500\", false)\'  ",
    "root -l -b \'plot.C(\"All\", \"jetpt2\", \"p_{T}\", 4, 400., true, false, false, \"500\", false)\'  ",

# after m_ll, m_jj inv mass cuts 
    "root -l -b \'plot.C(\"All\", \"TCHEJet\", \"TCHE (all jets)\", 1, 50., true, false, true, \"500\", false, 1.7, 30)\'  ",
    "root -l -b \'mergedPlot.C(\"All\", \"TCHEJet\", \"TCHE (all jets)\", 1, 50., true, false, true, \"500\")\'  ",
    "root -l -b \'plot.C(\"All\", \"metSig\", \"metSignificance\", 1, 5., true, false, true, \"500\", false, 0, 10)\'  ",
    "root -l -b \'mergedPlot.C(\"All\", \"metSig\", \"metSignificance\", 1, 5., false, false, true, \"500\")\'  ",
    "root -l -b \'mergedPlot.C(\"All\", \"CSVJet\", \"CSV (all jets)\", 1, 50., false, false, false, \"500\")\'  ",

# after 2b-tag selection
  "root -l -b \'mergedPlot.C(\"All\", \"metSignif\", \"metSignificance\", 2, 60., false, false, false)\'  ",

    "root -l -b \'plot.C(\"All\", \"zllpt\", \"p_{T}\", 10, 10., false, false, true, \"500\", false)\'  ",
    "root -l -b \'plot.C(\"All\", \"DRjj\", \"#Delta R_{jj}\", 5, 10., false, false, true, \"500\", false )\'  ",
    "root -l -b \'plot.C(\"All\", \"metSignif\", \"metSignificance\", 4, 10., false, false, true, \"500\", false, 0, 10 )\'  ",
    "root -l -b \'plot.C(\"All\", \"lljjmass2tags\", \"m_{lljj}\", 20, 10., false, false, false,\"500\", false, 400, 600 )\'  "
    "root -l -b \'plot.C(\"All\", \"HelyLDRefit2tags\", \"HelyLD\", 5, 20., false, false, true, \"500\", false)\'  "

# after metSig selection
    "root -l -b \'plot.C(\"All\", \"zllptBefLD\", \"p_{T}\", 5, 10., false, false, true, \"500\", false)\'  ",
    "root -l -b \'plot.C(\"All\", \"DRjjBefLD\", \"#Delta(R)\", 5, 10., false, false, true, \"500\", false )\'  ",
    "root -l -b \'plot.C(\"All\", \"lljjmassBefLD\", \"m_{lljj}\", 20, 10., false, false, false,\"500\", false )\'  "
    "root -l -b \'mergedPlot.C(\"All\", \"lljjmassBefLD\", \"m_{lljj}\", 20, 100., false, false, false,\"500\")\'  "

######### ANGULAR VARIABLES #############
    "root -l -b \'mergedPlot.C(\"All\", \"HelyLDRefit\", \"HelyLD\", 4, 50., false, false, true, \"500\")\'  ",
    "root -l -b \'mergedPlot.C(\"All\", \"cosTheta1StarRefit\", \"cos#Theta^{*}\", 4, 50., false, false, true, \"500\")\'  ",
    "root -l -b \'mergedPlot.C(\"All\", \"cosTheta1Refit\", \"cos#Theta_{1}\", 4, 50., false, false, true, \"500\")\'  ",
    "root -l -b \'mergedPlot.C(\"All\", \"cosTheta2Refit\", \"cos#Theta_{2}\", 5, 50., false, false, true, \"500\")\'  ",
    "root -l -b \'mergedPlot.C(\"All\", \"phiRefit\", \"#Phi\", 4, 50, false, false, true, \"500\")\'  ",
    "root -l -b \'mergedPlot.C(\"All\", \"phiStarRefit\", \"#Phi^{*}\", 4, 50, false, false, true, \"500\")\'  ",


######### FINAL PLOTS  #############

   "root -l -b \'plot.C(\"All\", \"lljjmass\", \"m_{lljj}\", 20, 10., false, false, false,\"250\")\'  ",
   "root -l -b \'plot.C(\"All\", \"lljjmass\", \"m_{lljj}\", 20, 10., false, false, false,\"300\")\'  ",
   "root -l -b \'plot.C(\"All\", \"lljjmass\", \"m_{lljj}\", 20, 10., false, false, false,\"350\")\'  ",
   "root -l -b \'plot.C(\"All\", \"lljjmass\", \"m_{lljj}\", 20, 10., false, false, false,\"400\")\'  ",
   "root -l -b \'plot.C(\"All\", \"lljjmass\", \"m_{lljj}\", 20, 10., false, false, false,\"450\")\'  ",
   "root -l -b \'plot.C(\"All\", \"lljjmass\", \"m_{lljj}\", 20, 10., false, false, false,\"500\", true)\'  ",
   "root -l -b \'mergedPlot.C(\"All\", \"lljjmass\", \"m_{lljj}\", 20, 50., false, false, false,\"500\")\'  "

    ]

for cmd in drawplotsCmd:
    os.system(cmd)
