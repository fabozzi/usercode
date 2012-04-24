## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# turn on when running on MC
runOnMC = False

# AK5 sequence with pileup substraction is the default
# the other sequences can be turned off with the following flags.
## True -> run also sequence without PU subtraction
runAK5NoPUSub = True 

### enable PU correction ##########
doJetPileUpCorrection = True
##################################

#add the L2L3Residual corrections only for data
if runOnMC:#MC
    jetCorrections=['L1FastJet','L2Relative','L3Absolute']
else:#Data
    jetCorrections=['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

############ general options ####################
process.options.wantSummary = True
process.maxEvents.input = 50
process.MessageLogger.cerr.FwkReport.reportEvery = 10
########### gloabl tag ############################
from CMGTools.Common.Tools.getGlobalTag import getGlobalTag
process.GlobalTag.globaltag = cms.string(getGlobalTag(runOnMC))
##################################################

############ PRINTOUT ###################

sep_line = "-" * 50
print sep_line
print 'running the following PFBRECO+PAT sequences:'
print '\tAK5'
if runAK5NoPUSub: print '\tAK5NoPUSub'
#print 'embedding in taus: ', doEmbedPFCandidatesInTaus
#print 'HPS taus         : ', hpsTaus
#print 'produce CMG tuple: ', runCMG
print 'run on MC        : ', runOnMC
print sep_line
print 'Global tag       : ', process.GlobalTag.globaltag
print sep_line

### INPUT COLLECTIONS ##########

process.source.fileNames = [
#    'file:/data3/scratch/cms/mc/Summer12/PU_S7_START52_V5-v2/DYJetsToLL_M-50/FE123555-F27A-E111-8E40-003048D46046.root'
#    'file:ACAEB147-ED80-E111-A73F-0025901D6268.root'
    'file:/data3/scratch/cms/data/Run2012A/DoubleMu/190704/A08772C8-4F83-E111-83B9-003048D2BA82.root'
]

### DEFINITION OF THE PFBRECO+PAT SEQUENCES ##########
# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.out.fileName = cms.untracked.string('patTuple_myTest.root')

# Configure PAT to use PFBRECO instead of AOD sources
# this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *

# ---------------- rho calculation for JEC ----------------------

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets

process.kt6PFJets = kt4PFJets.clone(
    rParam = cms.double(0.6),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
)

#compute rho correction for lepton isolation
############ to be specified Ghost_EtaMax???? ##############
#process.kt6PFJetsForIso = process.kt6PFJets.clone( Rho_EtaMax = cms.double(2.5), Ghost_EtaMax = cms.double(2.5) )
process.kt6PFJetsForIso = process.kt6PFJets.clone( Rho_EtaMax = cms.double(2.5) )
############################################################

# ---------------- Sequence AK5 ----------------------

# PFBRECO+PAT sequence 1:
# no lepton cleaning, AK5PFJets

postfixAK5 = "AK5"
jetAlgoAK5="AK5"

usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgoAK5, runOnMC=runOnMC, postfix=postfixAK5,
          jetCorrections=('AK5PFchs', jetCorrections))

############### remove useless modules ####################
def removeUseless( modName ):
    getattr(process,"patDefaultSequence"+postfixAK5).remove(
        getattr(process, modName+postfixAK5)
        )

removeUseless( "produceCaloMETCorrections" )
removeUseless( "pfCandsNotInJet" )
removeUseless( "pfJetMETcorr" )
removeUseless( "pfCandMETcorr" )
removeUseless( "pfchsMETcorr" )
removeUseless( "pfType1CorrectedMet" )
removeUseless( "pfType1p2CorrectedMet" )
#########################################################

############# set some parametern for leptons ###########
# removing stupid useless stuff from our muons:
getattr(process,"patMuons"+postfixAK5).embedCaloMETMuonCorrs = False 
getattr(process,"patMuons"+postfixAK5).embedTcMETMuonCorrs = False
#but embed the tracker track for cutting on 
getattr(process,"patMuons"+postfixAK5).embedTrack = True
#Patrick: Make the embedded track available for electrons (curing a bug in PAT)
getattr(process,"patElectrons"+postfixAK5).embedTrack = True

# removing default cuts on muons 
getattr(process,"pfMuonsFromVertexAK5").dzCut = 99
getattr(process,"pfMuonsFromVertexAK5").d0Cut = 99
getattr(process,"pfSelectedMuonsAK5").cut="pt()>3"

# removing default cuts on electrons 
getattr(process,"pfElectronsFromVertexAK5").dzCut = 99
getattr(process,"pfElectronsFromVertexAK5").d0Cut = 99
getattr(process,"pfSelectedElectronsAK5").cut="pt()>5"

############ recipe for PU correction ##########################

if doJetPileUpCorrection:
    from CommonTools.ParticleFlow.Tools.enablePileUpCorrection import enablePileUpCorrection
    enablePileUpCorrection( process, postfix=postfixAK5 )
################################################################

# NOT INTERESTED IN TAUS
#from CMGTools.Common.PAT.tauTools import *
#if doEmbedPFCandidatesInTaus:
#    embedPFCandidatesInTaus( process, postfix=postfixAK5, enable=True )
#if hpsTaus:
#    #Lucie : HPS is the default in 52X & V08-08-11-03, see PhysicsTools/PatAlgos/python/Tools/pfTools.py. Following line not needed (crash) in 5X.
#    import os
#    if os.environ['CMSSW_VERSION'] < "CMSSW_5_0":
#        adaptPFTaus(process,"hpsPFTau",postfix=postfixAK5) 
#    #  note that the following disables the tau cleaning in patJets
#    adaptSelectedPFJetForHPSTau(process,jetSelection="pt()>15.0",postfix=postfixAK5)
#    # currently (Sept 27,2011) there are three sets of tau isolation discriminators better to choose in CMG tuples.
#    removeHPSTauIsolation(process,postfix=postfixAK5)
###################################################################
   
# curing a weird bug in PAT..
from CMGTools.Common.PAT.removePhotonMatching import removePhotonMatching
removePhotonMatching( process, postfixAK5 )

###################################################################


# use non pileup substracted rho as in the Jan2012 JEC set
getattr(process,"patJetCorrFactors"+postfixAK5).rho = cms.InputTag("kt6PFJets","rho")
getattr(process,"patJets"+postfixAK5).addTagInfos = True
getattr(process,"patJets"+postfixAK5).tagInfoSources  = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAODAK5"),
    cms.InputTag("impactParameterTagInfosAODAK5")
    )
### non default embedding of AOD items for default patJets
getattr(process,"patJets"+postfixAK5).embedCaloTowers = False
getattr(process,"patJets"+postfixAK5).embedPFCandidates = True
###

getattr(process,"pfNoMuon"+postfixAK5).enable = False 
getattr(process,"pfNoElectron"+postfixAK5).enable = False 
getattr(process,"pfNoTau"+postfixAK5).enable = False 
getattr(process,"pfNoJet"+postfixAK5).enable = True
getattr(process,"pfIsolatedMuons"+postfixAK5).isolationCut = 999999
getattr(process,"pfIsolatedElectrons"+postfixAK5).isolationCut = 999999

########## insert the PFMET significance calculation #############
# insert the PFMET sifnificance calculation
from CMGTools.Common.PAT.addMETSignificance_cff import addMETSig
addMETSig( process, postfixAK5 )
####################################################################

######################################

#add user variables to PAT-jets 
process.customPFJets = cms.EDProducer(
    'PFJetUserData',
    JetInputCollection=cms.untracked.InputTag("selectedPatJetsAK5"),
    Verbosity=cms.untracked.bool(False)
    )

############## "Classic" PAT Muons and Electrons ########################
# (made from all reco muons, and all gsf electrons, respectively)
process.patMuons.embedTcMETMuonCorrs = False
process.patMuons.embedCaloMETMuonCorrs = False
process.patMuons.embedTrack = True

process.patElectrons.embedTrack = True
process.patElectrons.pfElectronSource = 'particleFlow'
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons', 'PFIso')
process.muIsoSequence = setupPFMuonIso(process, 'muons', 'PFIso')
adaptPFIsoMuons( process, applyPostfix(process,"patMuons",""), 'PFIso')
adaptPFIsoElectrons( process, applyPostfix(process,"patElectrons",""), 'PFIso')
process.stdMuonSeq = cms.Sequence(
    process.pfParticleSelectionSequence +
    process.muIsoSequence +
    process.makePatMuons +
    process.selectedPatMuons
    )
process.stdElectronSeq = cms.Sequence(
    process.pfParticleSelectionSequence +
    process.eleIsoSequence +
    process.makePatElectrons +
    process.selectedPatElectrons
    )

if not runOnMC:
    process.stdMuonSeq.remove( process.muonMatch )
    process.stdElectronSeq.remove( process.electronMatch )
    process.patMuons.embedGenMatch = False
    process.patElectrons.embedGenMatch = False

# Modules for the Cut-based Electron ID in the VBTF prescription
#import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi as vbtfid
# Switch to the official Electron VBTF Selection for 2011 Data (relax H/E cut in the Endcap):
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/VbtfEleID2011

import HiggsAnalysis.Higgs2l2b.simpleCutBasedElectronIDSummer11_cfi as vbtfid
process.eidVBTFRel95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95relIso' )
process.eidVBTFRel80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80relIso' )
process.eidVBTFCom95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95cIso'   )
process.eidVBTFCom80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80cIso'   )

process.eidSequence = cms.Sequence(
        process.eidVBTFRel95 +
        process.eidVBTFRel80 +
        process.eidVBTFCom95 +
        process.eidVBTFCom80
#        mvaTrigV0 +
#        mvaNonTrigV0 +
)

process.patElectronsAK5.electronIDSources = cms.PSet(
        eidVBTFRel95 = cms.InputTag("eidVBTFRel95"),
        eidVBTFRel80 = cms.InputTag("eidVBTFRel80"),
        eidVBTFCom95 = cms.InputTag("eidVBTFCom95"),
        eidVBTFCom80 = cms.InputTag("eidVBTFCom80"),
        #MVA (to be added)
#        mvaTrigV0 = cms.InputTag("mvaTrigV0"),
#        mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
)

process.patElectrons.electronIDSources = cms.PSet(
        eidVBTFRel95 = cms.InputTag("eidVBTFRel95"),
        eidVBTFRel80 = cms.InputTag("eidVBTFRel80"),
        eidVBTFCom95 = cms.InputTag("eidVBTFCom95"),
        eidVBTFCom80 = cms.InputTag("eidVBTFCom80"),
        #MVA (to be added)
#        mvaTrigV0 = cms.InputTag("mvaTrigV0"),
#        mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
)

getattr(process, 'patDefaultSequence' + postfixAK5).replace(
    getattr(process, "patElectrons" + postfixAK5),
    process.eidSequence  + getattr(process, "patElectrons" + postfixAK5) 
    )

# # adding custom detector based iso deposit ---> !!! this works only on V4 event content !!!
from RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi import *
from CMGTools.Common.PAT.addLeptCustomIsoDeposit_cff import addMuonCustomIsoDeposit
from CMGTools.Common.PAT.addLeptCustomIsoDeposit_cff import addElectronCustomIsoDeposit

# COLIN REMOVED CANNOT RUN
# uncomment -> segv 
# addMuonCustomIsoDeposit( process, 'patDefaultSequence', postfixAK5)
# uncomment -> segv
# track deposit alone -> segv
# ecal and hcal deposits alone -> segv
# addMuonCustomIsoDeposit( process, 'stdMuonSeq', '')
# addElectronCustomIsoDeposit( process, 'patDefaultSequence', postfixAK5)
# did not try to uncomment the following, as std electrons can't be made (bad refcore problem)
# addElectronCustomIsoDeposit( process, 'stdElectronSeq', '')



process.stdLeptonSequence = cms.Sequence(
    process.stdMuonSeq +
    process.eidSequence +
    process.stdElectronSeq 
    )

# Classic Electrons with UserData

process.userDataSelectedElectrons = cms.EDProducer(
    "Higgs2l2bElectronUserData",
    src = cms.InputTag("selectedPatElectrons"),
    rho = cms.InputTag("kt6PFJetsForIso:rho")
)

# ID-selected electrons (only ID and conversion, no isolation)
process.selectedIDElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("userDataSelectedElectrons"),
    cut = cms.string("(electronID('eidVBTFCom95') == 7) ||"               +
                     " (electronID('eidVBTFCom95') == 5) "
                     )
)

# Isolated electrons: standard isolation
process.selectedIsoElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("selectedIDElectrons"),
    cut = cms.string("electronID('eidVBTFCom95') == 7")
)

# Classic Muons with UserData
process.userDataSelectedMuons = cms.EDProducer(
    "Higgs2l2bMuonUserData",
    src = cms.InputTag("selectedPatMuons"),
    rho = cms.InputTag("kt6PFJetsForIso:rho")
)

# ID-selected muons
process.selectedIDMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("userDataSelectedMuons"),
    cut = cms.string(
            "pt > 10 && isGlobalMuon && isTrackerMuon && globalTrack().normalizedChi2 < 10 &&" +
            "globalTrack().hitPattern().numberOfValidTrackerHits > 10 && "                      +
            "globalTrack().hitPattern().numberOfValidPixelHits > 0 && "                         +
            "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "                         +
            "dB < 0.2 && numberOfMatches > 1 && abs(eta) < 2.4" )
)

# Isolated muons: standard isolation
process.selectedIsoMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedIDMuons"),
    cut = cms.string("trackIso + caloIso < 0.15 * pt")
)

process.userDataStandardLeptonSequence = cms.Sequence(
    process.userDataSelectedMuons *
    process.userDataSelectedElectrons *
    process.selectedIDMuons *
    process.selectedIDElectrons *
    process.selectedIsoMuons *
    process.selectedIsoElectrons 
    )


# Jet cleaning for patJets
process.cleanPatJetsIsoLept = cms.EDProducer(
    "PATJetCleaner",
    src = cms.InputTag("customPFJets"),
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(
    muons = cms.PSet(
    src = cms.InputTag("selectedIsoMuons"),
    algorithm = cms.string("byDeltaR"),
    preselection = cms.string(""),
    deltaR = cms.double(0.5),
    checkRecoComponents = cms.bool(False),
    pairCut = cms.string(""),
    requireNoOverlaps = cms.bool(True),
    ),
    electrons = cms.PSet(
    src = cms.InputTag("selectedIsoElectrons"),
    algorithm = cms.string("byDeltaR"),
    preselection = cms.string(""),
    deltaR = cms.double(0.5),
    checkRecoComponents = cms.bool(False),
    pairCut = cms.string(""),
    requireNoOverlaps = cms.bool(True),
    )
    ),
    finalCut = cms.string('')
    )


# ---------------- Sequence AK5NoPUSub, pfNoPileUp switched off ---------------

# PFBRECO+PAT sequence 2:
# pfNoPileUp switched off, AK5PFJets. This sequence is a clone of the AK5 sequence defined previously.

if runAK5NoPUSub:
    print 'Preparing AK5NoPUSub sequence...',

    postfixNoPUSub = 'NoPUSub'
    postfixAK5NoPUSub = postfixAK5+postfixNoPUSub

    from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
    cloneProcessingSnippet(process, getattr(process, 'patPF2PATSequence'+postfixAK5), postfixNoPUSub)

    getattr(process,"pfNoPileUp"+postfixAK5NoPUSub).enable = False
    getattr(process,"patJetCorrFactors"+postfixAK5NoPUSub).payload = "AK5PF"

    print 'Done'

    process.customPFJetsNoPUSub = cms.EDProducer(
        'PFJetUserData',
        JetInputCollection=cms.untracked.InputTag("selectedPatJetsAK5NoPUSub"),
        Verbosity=cms.untracked.bool(False)
        )

# Jet cleaning for patJets NoPUSub
    process.cleanPatJetsNoPUIsoLept = cms.EDProducer(
        "PATJetCleaner",
        src = cms.InputTag("customPFJetsNoPUSub"),
        preselection = cms.string(''),
        checkOverlaps = cms.PSet(
        muons = cms.PSet(
        src = cms.InputTag("selectedIsoMuons"),
        algorithm = cms.string("byDeltaR"),
        preselection = cms.string(""),
        deltaR = cms.double(0.5),
        checkRecoComponents = cms.bool(False),
        pairCut = cms.string(""),
        requireNoOverlaps = cms.bool(True),
        ),
        electrons = cms.PSet(
        src = cms.InputTag("selectedIsoElectrons"),
        algorithm = cms.string("byDeltaR"),
        preselection = cms.string(""),
        deltaR = cms.double(0.5),
        checkRecoComponents = cms.bool(False),
        pairCut = cms.string(""),
        requireNoOverlaps = cms.bool(True),
        )
        ),
        finalCut = cms.string('')
        )
    

# ---------------- Common stuff ---------------

process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
process.patTrigger.processName = cms.string('*')

### PATH DEFINITION #############################################

# counters that can be used at analysis level to know the processed events
process.prePathCounter = cms.EDProducer("EventCountProducer")
process.postPathCounter = cms.EDProducer("EventCountProducer")

# trigger information (no selection)

process.p = cms.Path( process.prePathCounter +
                      process.patTriggerDefaultSequence )
process.p += process.kt6PFJets
process.p += process.kt6PFJetsForIso

# PFBRECO+PAT ---

process.p += getattr(process,"patPF2PATSequence"+postfixAK5)

process.p += process.customPFJets

process.p += process.stdLeptonSequence
process.p += process.userDataStandardLeptonSequence
process.p += process.cleanPatJetsIsoLept



if runAK5NoPUSub:
    process.p += getattr(process,"patPF2PATSequence"+postfixAK5NoPUSub)
    process.p += process.customPFJetsNoPUSub
    process.p += process.cleanPatJetsNoPUIsoLept


# Select leptons
process.selectedPatMuons.cut = (
            "pt > 10 && abs(eta) < 2.4"
        )
process.selectedPatMuonsAK5.cut = (
            "pt > 10 && abs(eta) < 2.4"
        )
process.selectedPatElectronsAK5.cut = (
    "pt > 10.0 && abs(eta) < 2.5"
    )
process.selectedPatElectrons.cut = (
    "pt > 10.0 && abs(eta) < 2.5"
    )
# Select jets
process.selectedPatJetsAK5.cut = cms.string('pt > 25.0 && abs(eta) < 2.4')
process.selectedPatJetsAK5NoPUSub.cut = cms.string('pt > 25.0 && abs(eta) < 2.4')

################# COMBINATORIAL ANALYSIS ###########################

process.zee = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("userDataSelectedElectrons@+ userDataSelectedElectrons@-")
                             )

process.zmm = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("userDataSelectedMuons@+ userDataSelectedMuons@-")
                             )

process.zem = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("userDataSelectedElectrons@+ userDataSelectedMuons@-")
                             )

process.zjj = cms.EDProducer("CandViewShallowCloneCombiner",
                             checkCharge = cms.bool(False),
                             checkOverlap = cms.bool(False),
                             cut = cms.string(''),
                             decay = cms.string("cleanPatJetsIsoLept cleanPatJetsIsoLept")
)

process.hzzeejjBaseColl = cms.EDProducer("CandViewCombiner",
                                             checkCharge = cms.bool(False),
                                             cut = cms.string(''),
                                             decay = cms.string("zee zjj")
                                         )

process.hzzmmjjBaseColl = cms.EDProducer("CandViewCombiner",
                                             checkCharge = cms.bool(False),
                                             cut = cms.string(''),
                                             decay = cms.string("zmm zjj")
                                         )

process.hzzemjjBaseColl = cms.EDProducer("CandViewCombiner",
                                             checkCharge = cms.bool(False),
                                             cut = cms.string(''),
                                             decay = cms.string("zem zjj")
                                         )

process.hzzeejj = cms.EDProducer("Higgs2l2bUserDataNoMC",
                                     higgs = cms.InputTag("hzzeejjBaseColl"),
                                     gensTag = cms.InputTag("genParticles"),
                                     PFCandidates = cms.InputTag("particleFlow"),
                                     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                     dzCut = cms.double(0.1)
                                 )

process.hzzmmjj = cms.EDProducer("Higgs2l2bUserDataNoMC",
                                     higgs = cms.InputTag("hzzmmjjBaseColl"),
                                     gensTag = cms.InputTag("genParticles"),
                                     PFCandidates = cms.InputTag("particleFlow"),
                                     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                     dzCut = cms.double(0.1)
                                     )

process.hzzemjj = cms.EDProducer("Higgs2l2bUserDataNoMC",
                                     higgs = cms.InputTag("hzzemjjBaseColl"),
                                     gensTag = cms.InputTag("genParticles"),
                                     PFCandidates = cms.InputTag("particleFlow"),
                                     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                     dzCut = cms.double(0.1)
                                     )

####### saving also standard candles with PF leptons

process.zeePF = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("selectedPatElectronsAK5@+ selectedPatElectronsAK5@-")
                             )

process.zmmPF= cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("selectedPatMuonsAK5@+ selectedPatMuonsAK5@-")
                             )

process.zemPF = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("selectedPatElectronsAK5@+ selectedPatMuonsAK5@-")
                             )


process.combinatorialSequence = cms.Sequence(
    process.zee +
    process.zmm +
    process.zem +
    process.zjj +
    process.hzzeejjBaseColl +
    process.hzzmmjjBaseColl +
    process.hzzemjjBaseColl +
    process.hzzeejj +
    process.hzzmmjj +
    process.hzzemjj + 
    process.zeePF +
    process.zmmPF +
    process.zemPF 
)

process.p += process.combinatorialSequence

# event cleaning (in tagging mode, no event rejected)

process.load('CMGTools.Common.eventCleaning.eventCleaning_cff')

process.p += process.eventCleaningSequence

process.p += getattr(process,"postPathCounter") 
 
### Ntuplization ###
process.load("HiggsAnalysis.Higgs2l2b.Higgs2l2qedmNtuples_52_cff")

process.PUInfoNtuple = cms.EDProducer(
        "GenPUNtupleDump",
        isData = cms.bool(not runOnMC)
        )

# Event rho dumper
process.rhoDumper = cms.EDProducer(
    "EventRhoDumper",
    rho = cms.InputTag("kt6PFJets:rho"),
    restrictedRho = cms.InputTag("kt6PFJetsForIso:rho")
    )


# Met variables producer
process.metInfoProducer = cms.EDProducer(
    "MetVariablesProducer",
    metTag = cms.InputTag("patMETsAK5"),
    t1CorrMetTag = cms.InputTag("patType1CorrectedPFMetAK5")
    )


process.analysisPath = cms.Sequence(
    process.eventVtxInfoNtuple+
    process.PUInfoNtuple+
    process.rhoDumper+
    process.metInfoProducer+
    process.Higgs2e2bEdmNtuple+
    process.Higgs2mu2bEdmNtuple+
    process.Higgsemu2bEdmNtuple+
    process.jetinfos
    )

process.p += process.analysisPath

# Setup for a basic filtering
process.zll = cms.EDProducer("CandViewMerger",
                             src = cms.VInputTag("zee", "zmm", "zem")
)

process.zllFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("zll"),
                                 minNumber = cms.uint32(1),
)

process.jetFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("customPFJets"),
                                 minNumber = cms.uint32(2),
)

process.filterPath = cms.Path(
    process.zll *
    process.zllFilter *
    process.jetFilter
)

### OUTPUT DEFINITION #############################################
process.edmNtuplesOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('h2l2q_ntuple_52.root'),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_eventVtxInfoNtuple_*_*",
        "keep *_PUInfoNtuple_*_*",
        "keep *_rhoDumper_*_*",
        "keep *_metInfoProducer_*_*",
        "keep *_kt6PFJets_rho_PAT",
        "keep *_kt6PFJetsForIso_rho_*",
        "keep *_Higgs2e2bEdmNtuple_*_*",
        "keep *_Higgs2mu2bEdmNtuple_*_*",
        "keep *_Higgsemu2bEdmNtuple_*_*",
        "keep *_jetinfos_*_*"
        ),
    dropMetaData = cms.untracked.string('ALL'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring("filterPath")
        )
    )

process.edmNtuplesOut.outputCommands.extend([
    'keep *_ecalDeadCellTPfilter_*_*',
    'keep *_HBHENoiseFilterResultProducer*_*_*',
    'keep *_BeamHaloSummary_*_*',
    'keep *_recovRecHitFilter_*_*',
    'keep *_eeNoiseFilter_*_*',
    'keep *_trackingFailureFilter_*_*',
    'keep *_goodPrimaryVertexFilter_*_*',
    'keep *_scrapingFilter_*_*',
    'keep *_totalKinematicsFilterCMG_*_*'])

process.edmNtuplesOut.outputCommands.extend(['keep edmMergeableCounter_*_*_*'])

process.endPath = cms.EndPath(process.edmNtuplesOut)

