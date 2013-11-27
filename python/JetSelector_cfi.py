import FWCore.ParameterSet.Config as cms

defaultJetSelectionParameters = cms.PSet(
    jettype             = cms.string('AK5JET'),
    jetPtMin            = cms.double(50),
    jetPtMax            = cms.double(100000),
    jetAbsEtaMax        = cms.double(2.4),
    jetCSVDiscMin       = cms.double(0.000),   
    jetCSVDiscMax       = cms.double(1.000),   
    IsJetIDLoose        = cms.bool(True), 
    IsJetIDTight        = cms.bool(False), 
    )
