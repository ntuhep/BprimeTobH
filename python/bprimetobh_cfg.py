import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('BprimeTobH.root')
)


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'file:2011May10ReReco_numEvent1000.root'
    )
)

process.ntuple = cms.EDAnalyzer(
    'BprimeTobH',
    BeamSpotLabel = cms.InputTag('offlineBeamSpot')
)


process.p = cms.Path(process.ntuple)
