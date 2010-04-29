import FWCore.ParameterSet.Config as cms

process = cms.Process("RPCSkim")

### standard includes
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

### global tag
process.GlobalTag.globaltag = "GR_R_35X_V5::All"

### source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/601/EC82F50B-F040-DF11-8386-00E08179172F.root'
   ),
)

### number of events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) 


#############REAL DATA######################
# this is for filtering on L1 technical trigger bit
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
# bsc minbias and veto on beam halo
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

# this is for filtering on HLT MinBiasBSC bit (BPTX coincidence within 2 bunch crossings)
process.hltMinBiasBSC = cms.EDFilter("HLTHighLevel",
                                     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                     HLTPaths = cms.vstring("HLT_MinBiasBSC"
                                                            ),
                                     eventSetupPathsKey = cms.string(''),
                                     andOr = cms.bool(True),
                                     throw = cms.bool(True)
                                     )

# filter on good vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(15),	
                                           maxd0 = cms.double(2)	
                                           )

# filter to remove scraping ("monster") events (require at least 25% of high purity tracks)
process.scrapingFilter = cms.EDFilter("FilterOutScraping",
                                      applyfilter = cms.untracked.bool(True),
                                      debugOn = cms.untracked.bool(False),
                                      numtrack = cms.untracked.uint32(10),
                                      thresh = cms.untracked.double(0.25)
                                      )

###############################################

# Tracker muon filters
process.goodTrackerMuons = cms.EDFilter("CandViewSelector",
  src = cms.InputTag("muons"),
  cut = cms.string('isTrackerMuon = 1 & abs(eta) < 1.6'),
  filter = cms.bool(True)                                
)

process.idTrackerMuons = cms.EDFilter("TrackerMuonFilter",
  src = cms.InputTag("goodTrackerMuons"),
  algo = cms.string("TMLastStationAngTight")
) 

# Global muon filters
process.goodGlobalMuons = cms.EDFilter("CandViewSelector",
  src = cms.InputTag("muons"),
  cut = cms.string('isGlobalMuon = 1 & abs(eta) < 1.6'),
  filter = cms.bool(True)                                
)

process.validHitsSelector = cms.EDFilter("MuonValidHitsSelector",
  src = cms.InputTag("goodGlobalMuons"),
)
                                         
### paths

process.rpcSkimTrackerMuonPath = cms.Path(
  process.hltLevel1GTSeed +
  process.hltMinBiasBSC +
  process.primaryVertexFilter +
  process.scrapingFilter +
  process.goodTrackerMuons +
  process.idTrackerMuons 
)

process.rpcSkimGlobalMuonPath1 = cms.Path(
  process.hltLevel1GTSeed +
  process.hltMinBiasBSC +
  process.primaryVertexFilter +
  process.scrapingFilter +
  ~process.goodTrackerMuons +
  process.goodGlobalMuons +
  process.validHitsSelector  
)

process.rpcSkimGlobalMuonPath2 = cms.Path(
  process.hltLevel1GTSeed +
  process.hltMinBiasBSC +
  process.primaryVertexFilter +
  process.scrapingFilter +
  process.goodTrackerMuons +
  ~process.idTrackerMuons +
  process.goodGlobalMuons +
  process.validHitsSelector  
)

# Output module configuration
from Configuration.EventContent.EventContent_cff import *

rpcSkimEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring()
)
rpcSkimEventContent.outputCommands.extend(RECOEventContent.outputCommands)

rpcSkimEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
           'rpcSkimTrackerMuonPath','rpcSkimGlobalMuonPath1',
           'rpcSkimGlobalMuonPath2')
    )
)

process.rpcSkimOutputModule = cms.OutputModule("PoolOutputModule",
    rpcSkimEventContent,
    rpcSkimEventSelection,
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('RPCSkim'),
        dataTier = cms.untracked.string('USER')
   ),
   fileName = cms.untracked.string('rpcSkim_data.root')
)

process.outpath = cms.EndPath(process.rpcSkimOutputModule)


