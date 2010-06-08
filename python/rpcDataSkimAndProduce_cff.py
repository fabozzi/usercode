import FWCore.ParameterSet.Config as cms

from L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff import *
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import *
hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)

# bsc minbias and veto on beam halo
hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

# this is for filtering on HLT MinBiasBSC bit (BPTX coincidence within 2 bunch crossings)
hltMinBiasBSC = cms.EDFilter("HLTHighLevel",
                                     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                     HLTPaths = cms.vstring("HLT_MinBiasBSC"),
                                     eventSetupPathsKey = cms.string(''),
                                     andOr = cms.bool(True),
                                     throw = cms.bool(True)
                                     )
# filter on good vertex
primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(15),	
                                           maxd0 = cms.double(2)	
                                           )

# filter to remove scraping ("monster") events (require at least 25% of high purity tracks)
scrapingFilter = cms.EDFilter("FilterOutScraping",
                                      applyfilter = cms.untracked.bool(True),
                                      debugOn = cms.untracked.bool(False),
                                      numtrack = cms.untracked.uint32(10),
                                      thresh = cms.untracked.double(0.25)
                                      )

# filter on l1 mu paths and trigger given by dt or csc
from MuonAnalysis.Fabozzi.rpcSkimComplDetTrigger_cff import *

# Tracker muon filters
goodTrackerMuons = cms.EDFilter("CandViewSelector",
  src = cms.InputTag("muons"),
  cut = cms.string('isTrackerMuon = 1 & abs(eta) < 1.6'),
  filter = cms.bool(True)                                
)

idTrackerMuons = cms.EDFilter("TrackerMuonFilter",
  src = cms.InputTag("goodTrackerMuons"),
  algo = cms.string("TMLastStationAngTight")
) 

# Global muon filters
goodGlobalMuons = cms.EDFilter("CandViewSelector",
  src = cms.InputTag("muons"),
  cut = cms.string('isGlobalMuon = 1 & abs(eta) < 1.6'),
  filter = cms.bool(True)                                
)

validHitsSelector = cms.EDFilter("MuonValidHitsSelector",
  src = cms.InputTag("goodGlobalMuons"),
)

###############################################
# Re-make muon reco without RPC hits
# standalone muons
from RecoMuon.StandAloneMuonProducer.standAloneMuons_cff import *
standAloneMuonsNoRPC = standAloneMuons.clone()
standAloneMuonsNoRPC.STATrajBuilderParameters.FilterParameters.EnableRPCMeasurement = cms.bool(False)
standAloneMuonsNoRPC.STATrajBuilderParameters.BWFilterParameters.EnableRPCMeasurement = cms.bool(False)
# global muons
from RecoMuon.GlobalMuonProducer.GlobalMuonProducer_cff import *
globalMuonsNoRPC = globalMuons.clone()
globalMuonsNoRPC.MuonCollectionLabel = cms.InputTag("standAloneMuonsNoRPC","UpdatedAtVtx")
muontrackingNoRPC = cms.Sequence(
    standAloneMuonsNoRPC *
    globalMuonsNoRPC
)

# muons
from RecoMuon.MuonIdentification.muons_cfi import *
muonsNoRPC = muons.clone()
muonsNoRPC.inputCollectionLabels = cms.VInputTag(cms.InputTag("generalTracks"),cms.InputTag("globalMuonsNoRPC"),cms.InputTag("standAloneMuonsNoRPC","UpdatedAtVtx"))
muonsNoRPC.fillGlobalTrackQuality = False

muonIdProducerSequenceNoRPC = cms.Sequence(
    muonsNoRPC
)

###############################################
### sequences with no trigger selection

rpcSkimTrackerMuonPath = cms.Sequence(
  hltLevel1GTSeed +
  hltMinBiasBSC +
  primaryVertexFilter +
  scrapingFilter +
  goodTrackerMuons +
  idTrackerMuons +
  muontrackingNoRPC +
  muonIdProducerSequenceNoRPC
)

rpcSkimGlobalMuonPath1 = cms.Sequence(
  hltLevel1GTSeed +
  hltMinBiasBSC +
  primaryVertexFilter +
  scrapingFilter +
  ~goodTrackerMuons +
  goodGlobalMuons +
  validHitsSelector +
  muontrackingNoRPC +
  muonIdProducerSequenceNoRPC
)

rpcSkimGlobalMuonPath2 = cms.Sequence(
  hltLevel1GTSeed +
  hltMinBiasBSC +
  primaryVertexFilter +
  scrapingFilter +
  goodTrackerMuons +
  ~idTrackerMuons +
  goodGlobalMuons +
  validHitsSelector + 
  muontrackingNoRPC +
  muonIdProducerSequenceNoRPC
)

###############################################
### sequences with trigger selection

rpcSkimTrigAndTrackerMuonPath = cms.Sequence(
  hltLevel1GTSeed +
  hltMinBiasBSC +
  primaryVertexFilter +
  scrapingFilter +
  complDetTrigger +
  goodTrackerMuons +
  idTrackerMuons +
  muontrackingNoRPC +
  muonIdProducerSequenceNoRPC
)

rpcSkimTrigAndGlobalMuonPath1 = cms.Sequence(
  hltLevel1GTSeed +
  hltMinBiasBSC +
  primaryVertexFilter +
  scrapingFilter +
  complDetTrigger +
  ~goodTrackerMuons +
  goodGlobalMuons +
  validHitsSelector +
  muontrackingNoRPC +
  muonIdProducerSequenceNoRPC
)

rpcSkimTrigAndGlobalMuonPath2 = cms.Sequence(
  hltLevel1GTSeed +
  hltMinBiasBSC +
  primaryVertexFilter +
  scrapingFilter +
  complDetTrigger +
  goodTrackerMuons +
  ~idTrackerMuons +
  goodGlobalMuons +
  validHitsSelector + 
  muontrackingNoRPC +
  muonIdProducerSequenceNoRPC
)
