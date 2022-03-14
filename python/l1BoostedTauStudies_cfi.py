import FWCore.ParameterSet.Config as cms

l1NtupleProducer = cms.EDAnalyzer("BoostedTauStudies",                                  
                                  genParticles            = cms.InputTag("genParticles", "", "HLT")
)
