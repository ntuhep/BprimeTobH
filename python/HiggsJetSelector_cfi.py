from BpbH.BprimeTobH.FatJetSelector_cfi import * 

defaultHiggsJetSelectionParameters = defaultFatJetSelectionParameters.clone(
    jettype             = cms.string('HIGGSJET'), 
    subjet1CSVDiscMin   = cms.double(0.244),
    subjet1CSVDiscMax   = cms.double(1.000),
    subjet2CSVDiscMin   = cms.double(0.244),
    subjet2CSVDiscMax   = cms.double(1.000),
    )

