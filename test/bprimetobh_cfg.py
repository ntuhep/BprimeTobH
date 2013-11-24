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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.options = cms.untracked.PSet(
		SkipEvent = cms.untracked.vstring('ProductNotFound')
#		IgnoreCompletely = cms.untracked.vstring('ProductNotFound')
		)

process.TFileService = cms.Service(
		"TFileService",
    fileName = cms.string('BprimeTobH_v1.root')
)

process.source = cms.Source(
		"PoolSource",
		fileNames = cms.untracked.vstring(
      'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/yxin/Jet/Run2012A-22Jan2013-v1_TLBSM_53x_v3/45cbb6c27540456f7aaf244304c73a89/tlbsm_53x_v3_data_13_1_7EQ.root' 
      #### BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_10_1_Tyi.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_11_1_nn1.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_12_1_5De.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_13_1_mgq.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_1_1_HK5.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_2_1_y9R.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_4_1_CUr.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_5_1_T7R.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_6_1_uVu.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_7_1_Hzg.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_8_1_0VX.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-850_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_9_1_CFC.root', 
      #### BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v2/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_10_1_Nf4.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v2/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_11_1_DV5.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v2/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_12_1_I6g.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v2/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_13_1_V02.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v2/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_4_1_FN6.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v2/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_5_1_7NR.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v2/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_6_1_zQ9.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v2/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_7_1_bDw.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v2/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_9_1_okB.root', 
      #### BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_10_1_EcS.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_11_1_5wf.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_12_1_pwH.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_13_1_6k0.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_14_1_4wq.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_1_1_PiR.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_2_1_fWu.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_3_1_76H.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_4_1_BkD.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_5_1_Hjy.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_6_1_t0M.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_7_1_YlO.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_8_1_rQQ.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3_bugfix_v1/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_fat_bTagInfos_9_1_p5T.root', 
      #### TTJets_HadronicMGDecays_8TeV-madgraph 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_988_2_4Xr.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_989_3_rnE.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_98_1_XT8.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_990_2_HWP.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_991_2_Hjn.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_992_3_UY0.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_993_2_uA6.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_994_2_LGb.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_995_3_iew.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_996_2_8yt.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_997_2_ATn.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_998_2_aMo.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_999_1_7Dn.root', 
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/devdatta/TTJets_HadronicMGDecays_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_99_1_zBb.root', 
      #### BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph
      #'dcache:/pnfs/cms/WAX/11/store/user/lpctlbsm/ferencek/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1_TLBSM_53x_v3/99bd99199697666ff01397dad5652e9e/tlbsm_53x_v3_mc_10_1_FZV.root' 
			)
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
    DoGenJets          = cms.untracked.bool(False), 
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

