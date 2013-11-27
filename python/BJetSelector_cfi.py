import FWCore.ParameterSet.Config as cms

from  BpbH.BprimeTobH.JetSelector_cfi import * 

defaultBJetSelectionParameters = defaultJetSelectionParameters.clone( 
    jettype             = cms.string('BTAGGEDAK5JET'), 
    jetCSVDiscMin       = cms.double(0.679),   
    jetCSVDiscMax       = cms.double(1.000),   
    )
