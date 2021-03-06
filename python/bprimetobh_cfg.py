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
    BeamSpotLabel = cms.InputTag('offlineBeamSpot'), 
    VertexLabel = cms.InputTag('offlinePrimaryVertices'),
    VertexBSLabel = cms.InputTag('offlinePrimaryVerticesWithBS'), 
    muonlabel = cms.VInputTag('selectedPatMuonsPFlowLoose','selectedPatMuons'),
    electronlabel = cms.VInputTag('selectedPatElectronsPFlowLoose', 
                              'selectedPatElectrons'),
    jetlabel = cms.VInputTag('goodPatJetsCA8PrunedPacked', 'AF5jets'),
    
    LepCollections = cms.vstring('PFLepInfo', 'LepInfo'),
    JetCollections = cms.vstring('CA8PrunedInfo', 'CA8PrunedSubjet1Info'),
    #JetType = cms.vint32(0,2), # 0: pfjet, 1: calojet, 2: fatjet
    JetTypes = cms.vstring('fatjet', 'subjet'),
)


process.p = cms.Path(process.ntuple)
