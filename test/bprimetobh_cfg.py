import FWCore.ParameterSet.Config as cms

#from FWCore.ParameterSet.VarParsing import VarParsing
import copy

# options = VarParsing ('python')

# options.register('outFilename',
# 		'BprimeTobH.root', 
# 		VarParsing.multiplicity.singleton,
# 		VarParsing.varType.string,
# 		"Output file name"
# 		)

# options.register('useData',
# 		False,
# 		VarParsing.multiplicity.singleton,
# 		VarParsing.varType.int,
# 		'Run this on real data')

# options.setDefault('maxEvents', -1) 

# options.parseArguments()

#print options

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

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet(
		SkipEvent = cms.untracked.vstring('ProductNotFound')
#		IgnoreCompletely = cms.untracked.vstring('ProductNotFound')
		)

process.TFileService = cms.Service(
		"TFileService",
    #fileName = cms.string('root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/SingleBprimeTobH_v1.root')
    fileName = cms.string('bprimeTobH.root')
)

from inputFiles_cfi import * 

process.source = cms.Source(
		"PoolSource",
		fileNames = cms.untracked.vstring(
      FileNames
      ), 
		)

process.ntuple = cms.EDAnalyzer(
    'BprimeTobH',
    IncludeL7          = cms.untracked.bool(False), 
    DoElectrons        = cms.untracked.bool(False), 
    BeamSpotLabel      = cms.InputTag('offlineBeamSpot'), 
    VertexLabel        = cms.InputTag('goodOfflinePrimaryVertices'),
    muonlabel          = cms.VInputTag('selectedPatMuonsPFlow','selectedPatMuons'),
    electronlabel      = cms.VInputTag('selectedPatElectrons::PAT', 'selectedPatElectronsPFlow::PAT'),
    jetlabel           = cms.InputTag('goodPatJetsPFlow'),  
    fatjetlabel        = cms.InputTag('goodPatJetsCA8PF'), 
    prunedfatjetlabel  = cms.InputTag('goodPatJetsCA8PrunedPFPacked'),  
    subjetlabel        = cms.InputTag('selectedPatJetsCA8PrunedSubjetsPF'),  
    genjetlabel        = cms.InputTag('ak5GenJetsNoNu'),  
    hltlabel           = cms.VInputTag("TriggerResults::HLT"),
    genevtlabel        = cms.VInputTag("generator"),
    gtdigilabel        = cms.VInputTag("gtDigis"),
    rhocorrectionlabel = cms.VInputTag("kt6PFJetsForIsolation:rho","kt6PFJetsForIsolation:rho"),  # [electron,muon]
    sigmaLabel         = cms.VInputTag("kt6PFJetsForIsolation:sigma","kt6PFJetsForIsolation:rho"),    # [electron,muon]
    puInfoLabel        = cms.VInputTag("addPileupInfo"),
    genlabel           = cms.InputTag("prunedGenParticles"), 
    LepCollections     = cms.vstring('PFLepInfo', 'LepInfo'),
    JetCollections     = cms.vstring('FatJetInfo', 'SubJetInfo', 'JetInfo', 'GenJetInfo'), 
    JetTypes           = cms.vstring('fatjet', 'subjet', 'jet', 'genjet'),
    DoGenInfo          = cms.untracked.bool(True), 
    DoGenJets          = cms.untracked.bool(True), 
    DoJESUncert        = cms.untracked.bool(False), 
    DoJERUncert        = cms.untracked.bool(False), 
    DobtagUncert       = cms.untracked.bool(False), 
    DoHtagUncert       = cms.untracked.bool(False), 
    JetPtMin           = cms.untracked.double(20), 
    JetYMax            = cms.untracked.double(2.5),  
    FatJetPtMin        = cms.untracked.double(150), 
)


#------------------------------------------------------------
#  Gluon Tagger 
#------------------------------------------------------------
# process.load('QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff')  
# process.goodOfflinePrimaryVerticesQG.src  = cms.InputTag('goodOfflinePrimaryVertices')
# process.QGTagger.srcJets = cms.InputTag("selectedPatJetsPFlow")
# process.QGTagger.isPatJet  = cms.untracked.bool(True) 
# process.QGTagger.useCHS  = cms.untracked.bool(True)
# process.QGTagger.srcRho  = cms.InputTag('kt6PFJets','rho') 
# process.QGTagger.srcRhoIso  = cms.InputTag('kt6PFJets','rho') 

#process.p = cms.Path(process.QuarkGluonTagger*process.ntuple)
process.p = cms.Path(process.ntuple)

