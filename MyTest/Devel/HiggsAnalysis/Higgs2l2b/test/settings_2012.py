import ROOT


#data = False
data = True
run = 'AB'
#ch = "El"
ch = 'Mu'


#ch = 'El'
#ch = 'All'

#normalizeToData = False
normalizeToData = True



folder = "26Oct"+ch+"Channel/"


samples = [ "WW", "WZ", "ZZ", "TT", "DYJetsToLL_M-50"]
#samples = [ "WW", "WZ", "ZZ", "TT", "DY0", "DY1", "DY2", "DY3", "DY4"]
#samples = [ "WW", "WZ",  'T-tW', 'Tbar-tW',"ZZ", "TT", "DYJets"]


## lumiA = 2132.6
## lumiB = 2508.0

if(ch=="Mu"):
    lumiA = 709
    lumiB = 4404
if(ch =="El"):
    lumiA = 705
    lumiB = 4404

    
zerobtag = 216736
onebtag = 21344
twobtag = 2136
final = 300

Vars = {}


#### Vars[variableName] = (rebin, xmin, xmax)


## Vars["zllmass"] = (1,0.,200.,False)
if (ch == 'Mu'): Vars["stdCandle"] = (1,70.,110.,False, 'm_{#mu#mu}')
if (ch == 'El'): Vars["stdCandle"] = (1,70.,110.,False, 'm_{ee}')
if (ch == 'All'): Vars["stdCandle"] = (1,70.,110.,False, 'm_{ll}')

#Vars["zllmass"] = (1,0.,200.,True)
#Vars["stdCandle"] = (1,70.,110.,True)

Vars["zjjmass"] = (5,0.,0.,False,'m_{jj}')
Vars["jetpt1"] = (4,0.,0., False, 'P_{T} (leading jet)')
if (ch == 'Mu'):Vars["leptpt1"] = (4,0.,0., False,'P_{T} (leading #mu)')
if (ch == 'El'):Vars["leptpt1"] = (4,0.,0., True,'P_{T} (leading e)')
if (ch == 'All'):Vars["leptpt1"] = (4,0.,0., True,'P_{T} (leading lepton)')
Vars["jetpt2"] = (4,0.,0., True, 'P_{T} (sub-leading jet)')
if (ch == 'Mu'):Vars["leptpt2"] = (4,0.,200., True,'P_{T} (sub-leading #mu)')
if (ch == 'El'):Vars["leptpt2"] = (4,0.,200., True,'P_{T} (sub-leading e)')
if (ch == 'All'):Vars["leptpt2"] = (4,0.,200., True,'P_{T} (sub-leading lepton)')

if (ch == 'Mu'):Vars["lepteta1"] = (1,-3.,3., False, '#eta (leading #mu)')
if (ch == 'Mu'):Vars["lepteta2"] = (1,-3.,3., False, '#eta (sub-leading #mu)')
if (ch == 'El'):Vars["lepteta1"] = (1,-3.,3., False, '#eta (leading e)')
if (ch == 'El'):Vars["lepteta2"] = (1,-3.,3., False, '#eta (sub-leading e)')
if (ch == 'All'):Vars["lepteta1"] = (1,-3.,3., False, '#eta (leading lepton)')
if (ch == 'All'):Vars["lepteta2"] = (1,-3.,3., False, '#eta (sub-leading lepton)')


#Vars["lepteta1"] = (1,-3.,3., True)
#Vars["lepteta2"] = (1,-3.,3., True)
## Vars["cosTheta1Refit"] = (10,0.,0.,True)
## Vars["cosTheta1StarRefit"] = (10,0.,0.,True)
## Vars["cosTheta2Refit"] = (10,0.,0.,True)
## Vars["phiRefit"] = (10,0.,0.,True)
## Vars["phiStarRefit"] = (10,0.,0.,True)
## Vars["HelyLDRefit"] = (10,0.,0.,True)

Vars["cosTheta1Refit"] = (5,0.,0.,False, 'cos#Theta_{1}')
Vars["cosTheta1StarRefit"] = (5,0.,0.,False, 'cos#Theta_{1}^{*}')
Vars["cosTheta2Refit"] = (5,0.,0.,False, 'cos#Theta_{2}')
Vars["phiRefit"] = (5,0.,0.,False, '#Phi')
Vars["phiStarRefit"] = (5,0.,0.,False, '#Phi^{*}')
Vars["HelyLDRefit"] = (5,0.,0.,False, 'Helicity LD')

Vars["TCHEJet"] = (1,0.,0., True, 'TCHE')
Vars["CSVJet"] = (1,0.,0., True, 'CSV')
Vars["JPJet"] = (2,0.,0., True, 'JP')
Vars["metSignif"] = (2,0.,0.,True, 'MET Significance')

if (ch == 'Mu'):Vars["npv"] = (1,0.,0.,False,'# of primary vertexes - #mu channel')
if (ch == 'El'):Vars["npv"] = (1,0.,0.,False,'# of primary vertexes - e channel')
Vars["npv1"] = (1,0.,0.,True, '# of primary vertexes')
Vars["npv1"] = (1,0.,0.,False, '# of primary vertexes')

if (ch == 'Mu'):Vars["npv_woReweight"] = (1,0.,0.,False, '# of primary vertexes (without reweighting) - #mu channel')
if (ch == 'El'):Vars["npv_woReweight"] = (1,0.,0.,False, '# of primary vertexes (without reweighting) - e channel')

#Vars["njets"] = (1,0.,0.,True,'# of jets')

#Vars["leptphi1"] = (2,0.,0.,True)
#Vars["qgd"] = (2,0.,0.,True)
Vars["qgd"] = (2,0.,0.,False, 'Quark-Gluon Discriminant')

#Vars["npv2"] = (1,0.,0.,False)
if (ch == 'Mu'):Vars["lljjmass_0btagSB"]=(10, 0., 0., True, "m_{#mu#mujj} sidebands - 0btag")
if (ch == 'El'):Vars["lljjmass_0btagSB"]=(10, 0., 0., True, "m_{eejj} sidebands - 0btag")
if (ch == 'Mu'):Vars["lljjmass_1btagSB"]=(10, 0., 0., True, "m_{#mu#mujj} sidebands - 1btag")
if (ch == 'El'):Vars["lljjmass_1btagSB"]=(10, 0., 0., True, "m_{eejj} sidebands - 1btag")
if (ch == 'Mu'):Vars["lljjmass_2btagSB"]=(20, 0., 0., True, "m_{#mu#mujj} sidebands - 2btag")
if (ch == 'El'):Vars["lljjmass_2btagSB"]=(20, 0., 0., True, "m_{eejj} sidebands - 2btag")

## if (ch == 'Mu'):Vars["lljjmass_0btag"]=(10, 0., 0., True, "m_{#mu#mujj} - 0btag")
## if (ch == 'El'):Vars["lljjmass_0btag"]=(10, 0., 0., True, "m_{eejj} - 0btag")
## if (ch == 'Mu'):Vars["lljjmass_1btag"]=(10, 0., 0., True, "m_{#mu#mujj} - 1btag")
## if (ch == 'El'):Vars["lljjmass_1btag"]=(10, 0., 0., True, "m_{eejj} - 1btag")
## if (ch == 'Mu'):Vars["lljjmass"]=(20, 0., 0., True, "m_{#mu#mujj} - 2btag")
## if (ch == 'El'):Vars["lljjmass"]=(20, 0., 0., True, "m_{eejj} - 2btag")


if (ch == 'Mu'):Vars["lljjmass_0btag"]=(10, 0., 0., False, "m_{#mu#mujj} - 0btag")
if (ch == 'El'):Vars["lljjmass_0btag"]=(10, 0., 0., False, "m_{eejj} - 0btag")
if (ch == 'Mu'):Vars["lljjmass_1btag"]=(20, 0., 0., False, "m_{#mu#mujj} - 1btag")
if (ch == 'El'):Vars["lljjmass_1btag"]=(20, 0., 0., False, "m_{eejj} - 1btag")
if (ch == 'Mu'):Vars["lljjmass"]=(40, 0., 0., False, "m_{#mu#mujj} - 2btag")
if (ch == 'El'):Vars["lljjmass"]=(40, 0., 0., False, "m_{eejj} - 2btag")


Vars["lljjmass_sb_all"]=(10, 0., 0., True, "m{lljj} sidebands - allCategories")

#Vars["beta"] = (10,0.,0.,True, '#beta')
Vars["beta"] = (10,0.,0.,False, '#beta')

#Vars["jjm_0btag"] = (2,0.,0.,True, final)
#Vars["jjm_1btag"] = (2,0.,0.,False, final)
#Vars["jjm_2btag"] = (2,0.,0.,False, final)

#Vars["stdCandle_0btag"] = (2,60.,120.,False)
#Vars["stdCandle_1btag"] = (2,60.,120.,False)
#Vars["stdCandle_2btag"] = (2,60.,120.,False)
#Vars["stdCandle_final_0"] = (2,60.,120.,False)


#Vars["zjjmass_0btag"] = (10,0.,0.,False, zerobtag)
#Vars["zjjmass_1btag"] = (10,0.,0.,False, onebtag)
#Vars["zjjmass_2btag"] = (10,0.,0.,False, twobtag)

#Vars["met_0btag"] = (2,0.,0.,True, zerobtag)
#Vars["met_1btag"] = (2,0.,0.,True, onebtag)
#Vars["met_2btag"] = (2,0.,0.,True, twobtag)

#Vars["t1corrMET_0btag"] = (2,0.,0.,True)
#Vars["t1corrMET_1btag"] = (2,0.,0.,True )
#Vars["t1corrMET_2btag"] = (2,0.,0.,True )

#Vars["trkMET_0btag"] = (2,0.,0.,True)
#Vars["trkMET_1btag"] = (2,0.,0.,True)
#Vars["trkMET_2btag"] = (2,0.,0.,True)

#Vars["jzb_0btag"] = (5,-200.,200.,True, zerobtag)
#Vars["jzb_1btag"] = (5,-400.,400.,True, onebtag)
#Vars["jzb_2btag"] = (5,-400.,400.,True, twobtag)

#Vars["drJ1L_0btag"] = (2,0.,0.,True, zerobtag)
#Vars["drJ1L_1btag"] = (2,0.,0.,True, onebtag)
#Vars["drJ1L_2btag"] = (2,0.,0.,True, twobtag)

#Vars["drJ2L_0btag"] = (2,0.,0.,True, zerobtag)
#Vars["drJ2L_1btag"] = (2,0.,0.,True, onebtag)
#Vars["drJ2L_2btag"] = (4,0.,0.,True, twobtag)






#Vars["jzbmjj"] = (2,0.,0.,True)

## Vars["lljjmass"]=(5, 0., 0., True)
## Vars["lljjmass_btag"]=(5, 0., 0., True)
## Vars["lljjpt"]=(5, 0., 0., True)
## Vars["lljjY"]=(5, 0.,0., True)
## Vars["jetpt1"] = (4,0.,0., True)
## Vars["jetpt2"] = (4,0.,0., True)

## Vars["leptpt2"] = (4,0.,0., True)
## Vars["DRjj"] = (5,0.,0., True)
## Vars["DRjj_btag"] = (5,0.,0., True)
## Vars["zzdr"] = (5,0.,0., True)
## Vars["zzdeta"] = (5,0.,0., True)
## Vars["zzdphi"] = (5,0.,0., True)
#Vars["zzdpt"] = (5,0.,0., False)
## Vars["jjdeta"] = (5,0.,0., True)
## Vars["jjdphi"] = (5,0.,0., True)
## Vars["zzdr_btag"] = (5,0.,0., True)
## Vars["lldr"] = (5,0.,0., True)

## Vars["zllpt"] = (10,0.,0., True)
## Vars["zjjpt"] = (10,0.,0., True)

## Vars["zjjmass"] = (10,0.,0.,True)

