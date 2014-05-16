import FWCore.ParameterSet.Config as cms

HLTPathNames = [
    'HLT_HT750'
    ]

HLTMuPathNames = [
    'HLT_Mu22'
    ]
    
defaultTriggerSelectionParameters = cms.PSet(
    HLTPaths = cms.vstring(HLTPathNames) 
    )
    
defaultMuTriggerSelectionParameters = cms.PSet(
    HLTPaths = cms.vstring(HLTMuPathNames) 
    )
