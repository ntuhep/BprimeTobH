import FWCore.ParameterSet.Config as cms

PreselHLTPathNames = [
    'HLT_HT650'
    ]

HLTPathNames = [
    'HLT_HT750'
    ]

defaultTriggerSelectionParameters = cms.PSet(
    HLTPaths = cms.vstring(HLTPathNames) 
    )
