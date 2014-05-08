import FWCore.ParameterSet.Config as cms

from BpbH.BprimeTobH.JetSelector_cfi import * 

defaultFatJetSelectionParameters = cms.PSet(
    jettype             = cms.string('CA8JET'), 
    JetSelParams        = defaultJetSelectionParameters.clone(), 
    fatJetPtMin         = cms.double(150),
    fatJetPtMax         = cms.double(100000),
    fatJetAbsEtaMax     = cms.double(2.4),
    fatJetMassMin       = cms.double(0),
    fatJetMassMax       = cms.double(100000),
    fatJetPrunedMassMin = cms.double(0),
    fatJetPrunedMassMax = cms.double(100000),
    fatJetTau2ByTau1Max = cms.double(0.5),
    dRSubjetsMin        = cms.double(0.0), 
    dRSubjetsMax        = cms.double(999), 
    fatJetCSVDiscMin    = cms.double(-10000.), 
    fatJetCSVDiscMax    = cms.double(1.000), 
    subjet1CSVDiscMin   = cms.double(-10000),
    subjet1CSVDiscMax   = cms.double(1.000),
    subjet2CSVDiscMin   = cms.double(-10000.), 
    subjet2CSVDiscMax   = cms.double(1.000),
    )

