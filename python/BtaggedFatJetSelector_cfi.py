import FWCore.ParameterSet.Config as cms

from BpbH.BprimeTobH.FatJetSelector_cfi import * 

defaultBTaggedFatJetSelectionParameters = defaultFatJetSelectionParameters.clone(
    jettype             = cms.string('BTAGGEDFATJET'), 
    fatJetCSVDiscMin    = cms.double(0.244), 
    fatJetCSVDiscMax    = cms.double(1.000), 
    subjet1CSVDiscMin   = cms.double(0.000),
    subjet1CSVDiscMax   = cms.double(1.000),
    subjet2CSVDiscMin   = cms.double(0.000),
    subjet2CSVDiscMax   = cms.double(1.000),
    )

