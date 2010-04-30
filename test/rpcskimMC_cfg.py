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
process.GlobalTag.globaltag = "START3X_V26A::All"

### source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V26A_356ReReco-v1/0009/FEFC70B6-F53D-DF11-B57E-003048679150.root'
   ),
)

### number of events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(4000) ) 


#############MC######################

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
                                         
###############################################
# Re-make muon reco without RPC hits
# standalone muons
process.load("RecoMuon.StandAloneMuonProducer.standAloneMuons_cff")
process.standAloneMuonsNoRPC = process.standAloneMuons.clone()
process.standAloneMuonsNoRPC.STATrajBuilderParameters.FilterParameters.EnableRPCMeasurement = cms.bool(False)
process.standAloneMuonsNoRPC.STATrajBuilderParameters.BWFilterParameters.EnableRPCMeasurement = cms.bool(False)
# global muons
process.load("RecoMuon.GlobalMuonProducer.GlobalMuonProducer_cff")
process.globalMuonsNoRPC = process.globalMuons.clone()
process.globalMuonsNoRPC.MuonCollectionLabel = cms.InputTag("standAloneMuonsNoRPC","UpdatedAtVtx")
process.muontrackingNoRPC = cms.Sequence(
    process.standAloneMuonsNoRPC *
    process.globalMuonsNoRPC
)

# muons
process.load("RecoMuon.MuonIdentification.muons_cfi")
process.muonsNoRPC = process.muons.clone()
process.muonsNoRPC.inputCollectionLabels = cms.VInputTag(cms.InputTag("generalTracks"),cms.InputTag("globalMuonsNoRPC"),cms.InputTag("standAloneMuonsNoRPC","UpdatedAtVtx"))
process.muonsNoRPC.fillGlobalTrackQuality = False

process.muonIdProducerSequenceNoRPC = cms.Sequence(
    process.muonsNoRPC
)

###############################################
# Unpack L1Gt info
process.load("EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi")

process.rawProd = cms.Sequence(
   process.l1GtUnpack
)


###############################################
### paths

process.rpcSkimTrackerMuonPath = cms.Path(
  process.primaryVertexFilter +
  process.scrapingFilter +
  process.goodTrackerMuons +
  process.idTrackerMuons +
  process.muontrackingNoRPC +
  process.muonIdProducerSequenceNoRPC +
  process.rawProd
)

process.rpcSkimGlobalMuonPath1 = cms.Path(
  process.primaryVertexFilter +
  process.scrapingFilter +
  ~process.goodTrackerMuons +
  process.goodGlobalMuons +
  process.validHitsSelector +
  process.muontrackingNoRPC +
  process.muonIdProducerSequenceNoRPC +
  process.rawProd
)

process.rpcSkimGlobalMuonPath2 = cms.Path(
  process.primaryVertexFilter +
  process.scrapingFilter +
  process.goodTrackerMuons +
  ~process.idTrackerMuons +
  process.goodGlobalMuons +
  process.validHitsSelector +
  process.muontrackingNoRPC +
  process.muonIdProducerSequenceNoRPC +
  process.rawProd  
)

# Output module configuration
from Configuration.EventContent.EventContent_cff import *

rpcSkimEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

noRPCMuonsContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
      'keep *_*NoRPC_*_*',
      'keep L1MuRegionalCand*_*_*_*',
      'keep *_simMuonRPCDigis_*_*',
      'keep L1GlobalTriggerObjectMapRecord*_*_*_*',
      'keep *_g4SimHits_*_*'
      )
    )

rpcSkimEventContent.outputCommands.extend(RECOEventContent.outputCommands)
rpcSkimEventContent.outputCommands.extend(noRPCMuonsContent.outputCommands)

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
   fileName = cms.untracked.string('file:/tmp/fabozzi/rpcSkimNew_MinBias.root')
)

process.outpath = cms.EndPath(process.rpcSkimOutputModule)


