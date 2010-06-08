import FWCore.ParameterSet.Config as cms
from L1Trigger.Skimmer.l1Filter_cfi import *

l1Filter.algorithms = cms.vstring("L1_SingleMuOpen","L1_SingleMu0","L1_SingleMu3","L1_SingleMu5","L1_SingleMu5",
                                  "L1_SingleMu7","L1_SingleMu10","L1_SingleMu14","L1_SingleMu20","L1_SingleMu",
                                  "L1_DoubleMuTopBottom","L1_DoubleMuOpen","L1_DoubleMu3","L1_Mu3QE8_Jet6",
                                  "L1_Mu5QE8_Jet6")

detTrigger = cms.EDFilter("DetTriggerFilter",
                          gtTag=cms.InputTag("gtDigis"),
                          absMaxBX=cms.int32(0)
                          )

complDetTrigger = cms.Sequence(l1Filter +
                               detTrigger
                               )
