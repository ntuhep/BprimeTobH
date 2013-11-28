import FWCore.ParameterSet.Config as cms

HLTPathNames = [
    'HLT_HT750'
    ]

defaultTriggerSelectionParameters = cms.PSet(
    HLTPaths = cms.vstring(HLTPathNames) 
    )
