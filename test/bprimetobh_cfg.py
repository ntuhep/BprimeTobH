import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
import copy

options = VarParsing ('python')

options.register('outFilename',
		'BprimeTobH.root', 
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"Output file name"
		)

options.register('useData',
		False,
		VarParsing.multiplicity.singleton,
		VarParsing.varType.int,
		'Run this on real data')

options.setDefault('maxEvents', -1) 

options.parseArguments()

print options

process = cms.Process("Ntuple")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff') 
process.GlobalTag.globaltag = cms.string( 'START53_V21::All' )

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.MessageLogger = cms.Service("MessageLogger",
#		destinations = cms.untracked.vstring( 
#			'detailedInfo',
#			),
#		detailedInfo = cms.untracked.PSet(
#			threshold = cms.untracked.string('INFO'), 
#			), 
#		suppressInfo = cms.untracked.vstring("BprimeTobH"), 
#		)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.options = cms.untracked.PSet(
		SkipEvent = cms.untracked.vstring('ProductNotFound')
#		IgnoreCompletely = cms.untracked.vstring('ProductNotFound')
		)

process.TFileService = cms.Service(
		"TFileService",
		fileName = cms.string(options.outFilename) 
		)

process.source = cms.Source(
		"PoolSource",
		fileNames = cms.untracked.vstring(
			#'file:/tmp/petrakou/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3.root'
			'root://eoscms//eos/cms/store/user/devdatta/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3.root'
			#'file:tlbsm_53x_v3_mc_14_1_fZB.root'
			)
		)

process.ntuple = cms.EDAnalyzer(
    'BprimeTobH',
    IncludeL7 = cms.untracked.bool(False), 
    #DoSubjets = cms.untracked.bool(True), 
    BeamSpotLabel = cms.InputTag('offlineBeamSpot'), 
    VertexLabel = cms.InputTag('goodOfflinePrimaryVertices'),
    VertexBSLabel = cms.InputTag('goodOfflinePrimaryVertices'), 
    muonlabel = cms.VInputTag('selectedPatMuonsPFlow','selectedPatMuons'),
    electronlabel = cms.VInputTag('selectedPatElectrons::PAT', 'selectedPatElectronsPFlow::PAT'),
    jetlabel = cms.InputTag('goodPatJetsPFlow'),  
    fatjetlabel = cms.InputTag('goodPatJetsCA8PF'), 
    prunedfatjetlabel = cms.InputTag('goodPatJetsCA8PrunedPFPacked'),  
    subjetlabel = cms.InputTag('selectedPatJetsCA8PrunedSubjetsPF'),  
    genjetlabel = cms.InputTag('selectedPatJetsPFlow'),  
    hltlabel  = cms.VInputTag("TriggerResults::HLT"),
    gtdigilabel = cms.VInputTag("gtDigis"),
    genlabel = cms.InputTag("prunedGenParticles"), 
    LepCollections = cms.vstring('PFLepInfo', 'LepInfo'),

    # the jet branch names appear in the ntuple
    JetCollections = cms.vstring('FatJetInfo', 'SubJetInfo', 'JetInfo'), 
    # the types for the jets, must be correspond to the above types
    JetTypes = cms.vstring('fatjet', 'subjet', 'jet'),
    DoGenJets = cms.untracked.bool(False), 
    DoGenInfo = cms.untracked.bool(True), 
    JetMinPt = cms.untracked.double(20), # [GeV]
)


#------------------------------------------------------------
#  Gluon Tagger 
#------------------------------------------------------------
process.load('QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff')  
process.goodOfflinePrimaryVerticesQG.src  = cms.InputTag('goodOfflinePrimaryVertices')
process.QGTagger.srcJets = cms.InputTag("selectedPatJetsPFlow")
process.QGTagger.isPatJet  = cms.untracked.bool(True) 
process.QGTagger.useCHS  = cms.untracked.bool(True)
process.QGTagger.srcRho  = cms.InputTag('kt6PFJets','rho') 
process.QGTagger.srcRhoIso  = cms.InputTag('kt6PFJets','rho') 

#process.p = cms.Path(process.QuarkGluonTagger*process.ntuple)
process.p = cms.Path(process.ntuple)
