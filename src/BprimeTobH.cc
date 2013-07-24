// -*- C++ -*-
//
// Package:BprimeTobH
// Class:  BprimeTobH
// 
/**\class BprimeTobH BprimeTobH.cc HbbAna/BprimeTobH/src/BprimeTobH.cc

Description: [one line class summary]

Implementation:
Based on bprimeKit 
*/
//
// Original Author:  Xin Shi <Xin.Shi@cern.ch>
// Created:  Sat Apr  6 09:49:36 CEST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "../interface/format.h"
#include "../interface/TriggerBooking.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/Utilities/interface/EDMException.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


// For JEC
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h" // make isetup work

#include "fastjet/PseudoJet.hh"
#include "../interface/Njettiness.hh" 

#include "QuarkGluonTagger/EightTeV/interface/QGTagger.h"

using namespace std;

//
// constants, enums and typedefs
//

const unsigned int MAX_LEPCOLLECTIONS=3; 
const unsigned int MAX_JETCOLLECTIONS=5; 

typedef vector<pat::Jet> PatJetCollection;
typedef vector<reco::GenJet> GenJetCollection;
typedef map<const pat::Jet* ,const pat::Jet*> JetToJetMap;

//
// class declaration
//

class BprimeTobH : public edm::EDAnalyzer {
  public:
    explicit BprimeTobH(const edm::ParameterSet&);
    ~BprimeTobH();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    bool hasBeamSpot(const edm::Event&);
    void clearVariables(); 
    bool hasPrimaryVertex(const edm::Event &); 
    bool hasPrimaryVertexBS(const edm::Event &); 
    bool hasMuons(const edm::Event &); 
    bool hasElectrons(const edm::Event &); 
    bool hasJets(const edm::Event &, const edm::EventSetup&); 
    void saveGenInfo(const edm::Event &); 
    void saveHLT(const edm::Event&);
    void saveL1T(const edm::Event&);
    void processJets(const edm::Handle<PatJetCollection>&, const edm::Handle<PatJetCollection>&,
		     const edm::Event&, const edm::EventSetup&, const JetToJetMap&, const unsigned int); 
    void processGenJets(const edm::Handle<GenJetCollection>&,
			const edm::Event&, const edm::EventSetup&, const unsigned int); 

    // ----------member data ---------------------------
    TTree* tree_;  

    bool includeL7_;
    bool doElectrons_; 

    edm::InputTag BeamSpotLabel_;
    edm::InputTag VertexLabel_;
    edm::InputTag VertexBSLabel_;
    vector<edm::InputTag> muonlabel_;
    vector<edm::InputTag> electronlabel_;
    edm::InputTag jetlabel_;
    edm::InputTag fatjetlabel_;
    edm::InputTag prunedfatjetlabel_;
    edm::InputTag subjetlabel_;
    edm::InputTag genjetlabel_;
    vector<edm::InputTag> hltlabel_;
    vector<edm::InputTag> gtdigilabel_;
    edm::InputTag genlabel_; 

    EvtInfoBranches EvtInfo;
    GenInfoBranches GenInfo;
    VertexInfoBranches VertexInfo;
    LepInfoBranches LepInfo[MAX_LEPCOLLECTIONS];
    JetInfoBranches JetInfo[MAX_JETCOLLECTIONS];

    // Across the event 
    reco::BeamSpot beamSpot_;  
    reco::Vertex primaryVertex_;
    reco::Vertex primaryVertexBS_;

    vector<std::string> lepcollections_;
    vector<std::string> jetcollections_;
    vector<std::string> jettypes_;

    bool doGenJets_ ; 
    bool doGenInfo_;

    Njettiness nsubjettinessCalculator;
    double JetMinPt_; 
};

//
// static data member definitions
//

//
// constructors and destructor
//
BprimeTobH::BprimeTobH(const edm::ParameterSet& iConfig):
  tree_(0), 
  includeL7_(iConfig.getUntrackedParameter<bool>("IncludeL7",false)), 
  doElectrons_(iConfig.getUntrackedParameter<bool>("DoElectrons",false)), 
  BeamSpotLabel_(iConfig.getParameter<edm::InputTag>("BeamSpotLabel")),
  VertexLabel_(iConfig.getParameter<edm::InputTag>("VertexLabel")), 
  VertexBSLabel_(iConfig.getParameter<edm::InputTag>("VertexBSLabel")), 
  muonlabel_(iConfig.getParameter<vector<edm::InputTag> >("muonlabel")), 
  electronlabel_(iConfig.getParameter<vector<edm::InputTag> >("electronlabel")),  
  jetlabel_(iConfig.getParameter<edm::InputTag>("jetlabel")),  
  fatjetlabel_(iConfig.getParameter<edm::InputTag>("fatjetlabel")),  
  prunedfatjetlabel_(iConfig.getParameter<edm::InputTag>("prunedfatjetlabel")),  
  subjetlabel_(iConfig.getParameter<edm::InputTag>("subjetlabel")),  
  genjetlabel_(iConfig.getParameter<edm::InputTag>("genjetlabel")),  
  hltlabel_(iConfig.getParameter<vector<edm::InputTag> >("hltlabel")),  
  gtdigilabel_(iConfig.getParameter<vector<edm::InputTag> >("gtdigilabel")), 
  genlabel_(iConfig.getParameter<edm::InputTag>("genlabel")),
  lepcollections_(iConfig.getParameter<std::vector<std::string> >("LepCollections")),
  jetcollections_(iConfig.getParameter<std::vector<std::string> >("JetCollections")),
  jettypes_(iConfig.getParameter<std::vector<std::string> >("JetTypes")),
  doGenJets_(iConfig.getUntrackedParameter<bool>("DoGenJets")),
  doGenInfo_(iConfig.getUntrackedParameter<bool>("DoGenInfo")),
  nsubjettinessCalculator(Njettiness::onepass_kt_axes, NsubParameters(1.0, 0.8, 0.8)), 
  JetMinPt_(iConfig.getUntrackedParameter<double>("JetMinPt"))
{

  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("results") );

}

BprimeTobH::~BprimeTobH() { 

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void BprimeTobH::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { 
  clearVariables(); 

  EvtInfo.RunNo = iEvent.id().run();
  EvtInfo.EvtNo = iEvent.id().event();

  if ( hasBeamSpot(iEvent) 
      && hasPrimaryVertex(iEvent)
      && hasPrimaryVertexBS(iEvent)
     ) {

    hasMuons(iEvent); 
    if (doElectrons_) hasElectrons(iEvent); 
    hasJets(iEvent, iSetup); 
    saveHLT(iEvent); 
    saveL1T(iEvent); 

    if ( doGenInfo_ && !iEvent.isRealData() ) saveGenInfo(iEvent); 

    tree_->Fill();
  }

  clearVariables(); 
}

// ------------ method called once each job just before starting event loop  ------------
void BprimeTobH::beginJob() { 
  tree_ = new TTree ("tree", "BprimeTobH");
  EvtInfo.RegisterTree(tree_);  
  VertexInfo.RegisterTree(tree_);

  if(lepcollections_.size() > MAX_LEPCOLLECTIONS) 
    cout << "WARNING: Too many lep collections, using first " 
      << MAX_LEPCOLLECTIONS << endl;

  for(unsigned i=0; i<lepcollections_.size(); i++) {
    if(i >= MAX_LEPCOLLECTIONS) break;
    LepInfo[i].RegisterTree(tree_,lepcollections_[i]);
  }

  for(unsigned i=0; i<jetcollections_.size(); i++) {
    if(i >= MAX_JETCOLLECTIONS) break;

    JetInfo[i].RegisterTree(tree_,jetcollections_[i]);
  }

  if(doGenInfo_) GenInfo.RegisterTree(tree_);  

}

// ------------ method called once each job just after ending the event loop  ------------
void BprimeTobH::endJob() { 
}

// ------------ method called when starting to processes a run  ------------
void BprimeTobH::beginRun(edm::Run const&, edm::EventSetup const&) { 
}

// ------------ method called when ending the processing of a run  ------------
void BprimeTobH::endRun(edm::Run const&, edm::EventSetup const&) { 
}

// ------------ method called when starting to processes a luminosity block  ------------
void BprimeTobH::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { 
}

// ------------ method called when ending the processing of a luminosity block  ------------
void BprimeTobH::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BprimeTobH::fillDescriptions(edm::ConfigurationDescriptions& descriptions) { 
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool BprimeTobH::hasBeamSpot(const edm::Event& iEvent) { 

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(BeamSpotLabel_, beamSpotHandle);

  if ( ! beamSpotHandle.isValid() ) {
    edm::LogInfo("MyAnalyzer") << "No beam spot available from EventSetup" ; 
    return false; 
  }

  beamSpot_ = *beamSpotHandle; 

  EvtInfo.BeamSpotX = beamSpot_.position().x();
  EvtInfo.BeamSpotY = beamSpot_.position().y();
  EvtInfo.BeamSpotZ = beamSpot_.position().z();

  return true; 
}

void BprimeTobH::clearVariables() {  
  memset(&EvtInfo,0x00,sizeof(EvtInfo));
}

bool BprimeTobH::hasPrimaryVertex(const edm::Event& iEvent) {
  memset(&VertexInfo,0x00,sizeof(VertexInfo));

  edm::Handle<reco::VertexCollection> recoVertexHandle;
  iEvent.getByLabel(VertexLabel_, recoVertexHandle);

  if ( ! recoVertexHandle.isValid() 
      or recoVertexHandle.failedToGet() 
      or recoVertexHandle->size() <= 0 )
    return false;

  bool gotPrimVtx = false;

  for (std::vector<reco::Vertex>::const_iterator iVertex = recoVertexHandle->begin();
      iVertex != recoVertexHandle->end(); iVertex++) { 

    if (VertexInfo.Size>=MAX_VERTICES) {
      cout << " PV " << VertexInfo.Size << endl;
      fprintf(stderr,"ERROR: number of  Tracks exceeds the size of array.\n");
      break; 
    }

    VertexInfo.Type[VertexInfo.Size] = 0; //Vertices WITHOUT the Beam Spot constraint
    VertexInfo.isValid[VertexInfo.Size] = iVertex->isValid();
    VertexInfo.isFake[VertexInfo.Size] = iVertex->isFake(); //Uly 2011-05-16
    VertexInfo.Ndof[VertexInfo.Size] = iVertex->ndof();
    VertexInfo.NormalizedChi2[VertexInfo.Size] = iVertex->normalizedChi2();
    VertexInfo.x[VertexInfo.Size] = iVertex->x();
    VertexInfo.y[VertexInfo.Size] = iVertex->y();
    VertexInfo.z[VertexInfo.Size] = iVertex->z();
    VertexInfo.Rho[VertexInfo.Size] = iVertex->position().Rho();
    VertexInfo.Pt_Sum[VertexInfo.Size] = 0.;
    VertexInfo.Pt_Sum2[VertexInfo.Size] = 0.;

    if (!gotPrimVtx && (!iVertex->isFake() && iVertex->ndof()>=4. && iVertex->z() <=24.
          && iVertex->position().Rho()<=2.)) {
      primaryVertex_ = *(iVertex); 
      gotPrimVtx=true;
    }

    VertexInfo.Size++;
  }

  if (!primaryVertex_.isValid()) {
    return false; 
  }

  return true; 
}


bool BprimeTobH::hasPrimaryVertexBS(const edm::Event& iEvent) { 
  edm::Handle<reco::VertexCollection>  VertexHandleBS; 
  double PVBS_Pt_Max = -100.;

  iEvent.getByLabel(VertexBSLabel_, VertexHandleBS);

  if( ! VertexHandleBS.isValid() or VertexHandleBS.failedToGet() 
      or VertexHandleBS->size() <= 0 )
    return false;

  const vector<reco::Vertex> VerticesBS = *VertexHandleBS;
  for(vector<reco::Vertex>::const_iterator it_vtx = VerticesBS.begin();
      it_vtx != VerticesBS.end(); it_vtx++ ) {
    if (VertexInfo.Size>=MAX_VERTICES) {
      cout << "PVBS " << VertexInfo.Size << endl;
      fprintf(stderr,"ERROR: number of  Vertices exceeds the size of array.\n");
      break;
    }
    VertexInfo.Type[VertexInfo.Size] = 1;  //Vertices WITH the Beam Spot constraint
    VertexInfo.isValid[VertexInfo.Size] = it_vtx->isValid();
    VertexInfo.isFake[VertexInfo.Size] = it_vtx->isFake(); //Uly 2011-05-16
    VertexInfo.Ndof[VertexInfo.Size] = it_vtx->ndof();
    VertexInfo.NormalizedChi2[VertexInfo.Size] = it_vtx->normalizedChi2();
    VertexInfo.x[VertexInfo.Size] = it_vtx->x();
    VertexInfo.y[VertexInfo.Size] = it_vtx->y();
    VertexInfo.z[VertexInfo.Size] = it_vtx->z();
    VertexInfo.Rho[VertexInfo.Size] = it_vtx->position().Rho();
    VertexInfo.Pt_Sum[VertexInfo.Size] = 0.;
    VertexInfo.Pt_Sum2[VertexInfo.Size] = 0.;

    if( VertexInfo.Pt_Sum[VertexInfo.Size] >= PVBS_Pt_Max ){
      PVBS_Pt_Max = VertexInfo.Pt_Sum[VertexInfo.Size];
      primaryVertexBS_ = *it_vtx;
    }			
    VertexInfo.Size++;
  }

  if (!primaryVertexBS_.isValid()) 
    return false; 

  return true; 
}


bool BprimeTobH::hasMuons(const edm::Event& iEvent) { 
  vector<edm::Handle<vector<pat::Muon> > > MuonHandle;
  for(unsigned il=0; il<muonlabel_.size(); il++) {
    MuonHandle.push_back(edm::Handle<vector<pat::Muon> >());
    iEvent.getByLabel( muonlabel_[il], MuonHandle[il]);
  }

  for(unsigned icoll=0; icoll<lepcollections_.size(); icoll++) {
    //loop over collections
    if(icoll >= MAX_LEPCOLLECTIONS) break;

    memset(&LepInfo[icoll],0x00,sizeof(LepInfo[icoll]));

    if (MuonHandle.size() <= icoll) continue;  

    //loop over muons in collection
    for(vector<pat::Muon>::const_iterator it_mu = MuonHandle[icoll]->begin(); 
        it_mu != MuonHandle[icoll]->end(); it_mu++ ) { 

      if (LepInfo[icoll].Size>=MAX_LEPTONS) {
        fprintf(stderr,"ERROR: number of leptons exceeds the size of array.\n");
        break;//exit(0);
      }

      LepInfo[icoll].Index[LepInfo[icoll].Size] = LepInfo[icoll].Size;
    }
  }
  return true; 
}

bool BprimeTobH::hasElectrons(const edm::Event& iEvent) { 
  vector<edm::Handle<vector<pat::Electron> > > ElectronHandle;
  for(unsigned il=0; il<electronlabel_.size(); il++) {
    do { {	
      ElectronHandle.push_back(edm::Handle<vector<pat::Electron> >());
      iEvent.getByLabel( electronlabel_[il], ElectronHandle[il]);
    } {
      throw edm::Exception(edm::errors::NotFound) << " Object " << electronlabel_[il] << " not found" ; 
    } } while (false) ; 
  }

  for(unsigned icoll=0; icoll<lepcollections_.size(); icoll++) {
    //loop over collections
    if(icoll >= MAX_LEPCOLLECTIONS) break;

    memset(&LepInfo[icoll],0x00,sizeof(LepInfo[icoll]));

    if (ElectronHandle.size() <= icoll) continue;  

    for( vector<pat::Electron>::const_iterator it_el = ElectronHandle[icoll]->begin(); 
        it_el != ElectronHandle[icoll]->end(); it_el++ ) {

      if (LepInfo[icoll].Size>=MAX_LEPTONS) {
        fprintf(stderr,"ERROR: number of leptons exceeds the size of array.\n");
        break;//exit(0);
      }

      LepInfo[icoll].Index[LepInfo[icoll].Size] = LepInfo[icoll].Size;
    }

  }

  return true; 
}

bool BprimeTobH::hasJets(const edm::Event& iEvent, const edm::EventSetup& iSetup) { 
  // vector<edm::Handle<vector<pat::Jet> > > JetHandle;
  // for(unsigned il=0; il<jetlabel_.size(); il++) {
  //   JetHandle.push_back(edm::Handle<vector<pat::Jet> >());
  //   iEvent.getByLabel( jetlabel_[il], JetHandle[il]);
  // }

  // currently, just deal with these three jetlabels as defined in the py: 
  // jetlabel = cms.VInputTag(
  //   'selectedPatJets',
  //   'selectedPatJetsCA8PrunedPFPacked', 
  //   'selectedPatJetsCA8PrunedSubJetsPF'), 

  edm::Handle <PatJetCollection> jetsColl; //  = JetHandle[0];
  edm::Handle <PatJetCollection> fatjetsColl; //  = JetHandle[0];
  edm::Handle <PatJetCollection> prunedfatjetsColl; //  = JetHandle[1];
  edm::Handle <PatJetCollection> subjetsColl; //  = JetHandle[2];
  edm::Handle <GenJetCollection> genjetsColl; 

  // edm::Handle<pat::Jet> FatJetHandle;
  // iEvent.getByLabel( fatjetlabel_, FatJetHandle);
  iEvent.getByLabel( jetlabel_, jetsColl);
  iEvent.getByLabel( fatjetlabel_, fatjetsColl);
  iEvent.getByLabel( prunedfatjetlabel_, prunedfatjetsColl);
  iEvent.getByLabel( subjetlabel_, subjetsColl);
  iEvent.getByLabel( genjetlabel_, genjetsColl);
  // edm::Handle <PatJetCollection> fatjetsColl = FatJetHandle;

  // get the fatjet to prunedfat jet match map

  JetToJetMap fatJetToPrunedFatJetMap;

  for(PatJetCollection::const_iterator it = fatjetsColl->begin(); it != fatjetsColl->end(); ++it)
  {
    PatJetCollection::const_iterator prunedJetMatch;
    bool prunedJetMatchFound = false;
    float dR = 100.;//0.8; // hard coded for now. 
    for(PatJetCollection::const_iterator pjIt = prunedfatjetsColl->begin();
        pjIt != prunedfatjetsColl->end(); ++pjIt)
    {
      //cout << ". ";
      float dR_temp = reco::deltaR( it->p4(), pjIt->p4() );
      if( dR_temp < dR )
      {
        prunedJetMatchFound = true;
        dR = dR_temp;
        prunedJetMatch = pjIt;
      }
    }
    if( !prunedJetMatchFound )
      edm::LogError("NoMatchingGroomedJet") << " Matching pruned jet not found"; // This should never happen but just in case
    fatJetToPrunedFatJetMap[&(*it)] = &(*prunedJetMatch);
  }

  // Now process 'FatJetInfo', 'SubJetInfo', 'JetInfo', 'GenJetInfo':
  unsigned int iJetColl = 0 ; // FatJetInfo 
  processJets(fatjetsColl, subjetsColl, iEvent, iSetup, fatJetToPrunedFatJetMap, iJetColl) ;

  iJetColl = 1; // SubJetInfo 
  processJets(subjetsColl, fatjetsColl, iEvent, iSetup, fatJetToPrunedFatJetMap, iJetColl) ;

  iJetColl = 2 ; // JetInfo 
  processJets(jetsColl, subjetsColl, iEvent, iSetup, fatJetToPrunedFatJetMap, iJetColl) ;

  iJetColl = 3 ; // GenJetInfo 
  processGenJets(genjetsColl, iEvent, iSetup, iJetColl) ;

  return true; 
}

void BprimeTobH::processJets(const edm::Handle<PatJetCollection>& jetsColl, 
    const edm::Handle<PatJetCollection>& jetsColl2,
    const edm::Event& iEvent, 
    const edm::EventSetup& iSetup, 
    const JetToJetMap& fatJetToPrunedFatJetMap, 
    const unsigned int icoll) { 

  if(icoll >= MAX_JETCOLLECTIONS) return;

  memset(&JetInfo[icoll],0x00,sizeof(JetInfo[icoll]));

  //// Jet ID 
  PFJetIDSelectionFunctor pfjetIDLoose( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE ) ;
  PFJetIDSelectionFunctor pfjetIDTight( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT );
  pat::strbitset ret = pfjetIDLoose.getBitTemplate() ;

  // For Jet Uncertainty

  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  if (icoll == 0)  iSetup.get<JetCorrectionsRecord>().get("AK7PFchs", JetCorParColl); 
  else if (icoll == 1 || icoll == 2)  iSetup.get<JetCorrectionsRecord>().get("AK5PFchs", JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

  for( vector<pat::Jet>::const_iterator it_jet = jetsColl->begin();
      it_jet != jetsColl->end(); it_jet++ ) { 

    if (iEvent.isRealData() && it_jet->pt() < JetMinPt_) continue ;

    JetInfo[icoll].Index   [JetInfo[icoll].Size] = JetInfo[icoll].Size;
    JetInfo[icoll].NTracks [JetInfo[icoll].Size] = it_jet->associatedTracks().size();

    JetInfo[icoll].Et      [JetInfo[icoll].Size] = it_jet->et();
    JetInfo[icoll].Pt      [JetInfo[icoll].Size] = it_jet->pt();
    JetInfo[icoll].Eta     [JetInfo[icoll].Size] = it_jet->eta();
    JetInfo[icoll].Phi     [JetInfo[icoll].Size] = it_jet->phi();
    JetInfo[icoll].Energy  [JetInfo[icoll].Size] = it_jet->energy(); 
    JetInfo[icoll].Px      [JetInfo[icoll].Size] = it_jet->px(); 
    JetInfo[icoll].Py      [JetInfo[icoll].Size] = it_jet->py(); 
    JetInfo[icoll].Pz      [JetInfo[icoll].Size] = it_jet->pz(); 
    JetInfo[icoll].Mass    [JetInfo[icoll].Size] = it_jet->mass();
    JetInfo[icoll].Area    [JetInfo[icoll].Size] = it_jet->jetArea();

    jecUnc->setJetEta(it_jet->eta());
    jecUnc->setJetPt(it_jet->pt()); // here you must use the CORRECTED jet pt
    if(fabs(it_jet->eta())<=5.0) JetInfo[icoll].Unc  [JetInfo[icoll].Size] = jecUnc->getUncertainty(true);

    JetInfo[icoll].JetCharge   [JetInfo[icoll].Size] = it_jet->jetCharge();
    JetInfo[icoll].NConstituents[JetInfo[icoll].Size] = it_jet->numberOfDaughters();

    if (it_jet->isPFJet())  {
      JetInfo[icoll].NCH[JetInfo[icoll].Size] = it_jet->chargedMultiplicity();
      JetInfo[icoll].CEF[JetInfo[icoll].Size] = it_jet->chargedEmEnergyFraction();
      JetInfo[icoll].NHF[JetInfo[icoll].Size] = it_jet->neutralHadronEnergyFraction();
      JetInfo[icoll].NEF[JetInfo[icoll].Size] = it_jet->neutralEmEnergyFraction();
      JetInfo[icoll].CHF[JetInfo[icoll].Size] = it_jet->chargedHadronEnergyFraction();
    }

    bool JetIDLoose(false);
    bool JetIDTight(false);

    JetInfo[icoll].QGTagsMLP       [JetInfo[icoll].Size] = -999;
    JetInfo[icoll].QGTagsLikelihood       [JetInfo[icoll].Size] = -1;

    //// Gluon tagger
    edm::Handle<edm::ValueMap<float> >  QGTagsHandleMLP;
    edm::Handle<edm::ValueMap<float> >  QGTagsHandleLikelihood;
    iEvent.getByLabel("QGTagger","qgMLP", QGTagsHandleMLP);
    //iEvent.getByLabel("QGTagger","qgLikelihood", QGTagsHandleLikelihood);
    iEvent.getByLabel("QGTagger", QGTagsHandleLikelihood); 

    int ijet = it_jet - jetsColl->begin();
    edm::RefToBase<reco::Jet> jetRef(edm::Ref<std::vector <pat::Jet> >(jetsColl,ijet));

    if (QGTagsHandleMLP.isValid()){
      // std::cout << "QGTagsHandleMLP is Valid\n" ; 
      JetInfo[icoll].QGTagsMLP       [JetInfo[icoll].Size] = (*QGTagsHandleMLP)[jetRef];
    }
    // else std::cout << "QGTagsHandleMLP is not Valid\n" ;
    if (QGTagsHandleLikelihood.isValid()){
      // std::cout << "QGTagsHandleLikelihood is Valid\n" ; 
      JetInfo[icoll].QGTagsLikelihood       [JetInfo[icoll].Size] = (*QGTagsHandleLikelihood)[jetRef];
    } 
    // else std::cout << "QGTagsHandleLikelihood is not Valid\n" ; 

    if(it_jet->isPFJet()) {
      //Jet ID for PFJet
      ret.set(false);
      JetIDLoose = pfjetIDLoose(*it_jet, ret);
      ret.set(false);
      JetIDTight = pfjetIDTight(*it_jet, ret); 

    } //// If isPFJet 
    else { 
      JetIDLoose = false; 
      JetIDTight = false; 
    }

    JetInfo[icoll].JetIDLOOSE[JetInfo[icoll].Size] = (JetIDLoose) ?  1 : 0;
    JetInfo[icoll].JetIDTIGHT[JetInfo[icoll].Size] = (JetIDTight) ?  1 : 0;

    // Jet corrections, B-tagging, and Jet ID information
    // now we just fill everything (regardless of availability)
    JetInfo[icoll].PtCorrRaw   [JetInfo[icoll].Size] = it_jet->correctedJet("Uncorrected").pt();	   
    JetInfo[icoll].PtCorrL2[JetInfo[icoll].Size] = it_jet->correctedJet("L2Relative" ).pt(); // L2(rel) 
    JetInfo[icoll].PtCorrL3[JetInfo[icoll].Size] = it_jet->correctedJet("L3Absolute" ).pt(); // L3(abs) 
    if(includeL7_) {
      JetInfo[icoll].PtCorrL7g   [JetInfo[icoll].Size] = it_jet->correctedJet("L7Parton", "gluon" ).pt(); // L7(gluon)
      JetInfo[icoll].PtCorrL7uds [JetInfo[icoll].Size] = it_jet->correctedJet("L7Parton", "uds"   ).pt(); // L7(uds-jet) 
      JetInfo[icoll].PtCorrL7c   [JetInfo[icoll].Size] = it_jet->correctedJet("L7Parton", "charm" ).pt(); // L7(c-jet)  
      JetInfo[icoll].PtCorrL7b   [JetInfo[icoll].Size] = it_jet->correctedJet("L7Parton", "bottom").pt(); // L7(b-jet) 
    }

    JetInfo[icoll].JetBProbBJetTags[JetInfo[icoll].Size] = it_jet->bDiscriminator("jetBProbabilityBJetTags");
    JetInfo[icoll].JetProbBJetTags [JetInfo[icoll].Size] = it_jet->bDiscriminator("jetProbabilityBJetTags");
    JetInfo[icoll].TrackCountHiPurBJetTags [JetInfo[icoll].Size] = it_jet->bDiscriminator("trackCountingHighPurBJetTags"); 
    JetInfo[icoll].CombinedSVBJetTags  [JetInfo[icoll].Size] = it_jet->bDiscriminator("combinedSecondaryVertexBJetTags");
    JetInfo[icoll].CombinedSVMVABJetTags   [JetInfo[icoll].Size] = it_jet->bDiscriminator("combinedSecondaryVertexMVABJetTags");
    JetInfo[icoll].SoftElecByIP3dBJetTags  [JetInfo[icoll].Size] = it_jet->bDiscriminator("softElectronByIP3dBJetTags");
    JetInfo[icoll].SoftElecByPtBJetTags[JetInfo[icoll].Size] = it_jet->bDiscriminator("softElectronByPtBJetTags");
    JetInfo[icoll].SoftMuonBJetTags[JetInfo[icoll].Size] = it_jet->bDiscriminator("softMuonBJetTags");
    JetInfo[icoll].SoftMuonByIP3dBJetTags  [JetInfo[icoll].Size] = it_jet->bDiscriminator("softMuonByIP3dBJetTags");
    JetInfo[icoll].SoftMuonByPtBJetTags[JetInfo[icoll].Size] = it_jet->bDiscriminator("softMuonByPtBJetTags");
    JetInfo[icoll].DoubleSVHighEffBJetTags [JetInfo[icoll].Size] = it_jet->bDiscriminator("doubleSecondaryVertexHighEffBJetTags"); 

    if (!iEvent.isRealData()) {
      JetInfo[icoll].GenFlavor[JetInfo[icoll].Size] = it_jet->partonFlavour(); 
      if (doGenJets_) {
        const reco::GenJet * genjet = it_jet->genJet();
        JetInfo[icoll].GenJetPt [JetInfo[icoll].Size] = genjet->pt(); 
        JetInfo[icoll].GenJetEta[JetInfo[icoll].Size] = genjet->eta();
        JetInfo[icoll].GenJetPhi[JetInfo[icoll].Size] = genjet->eta(); 
        const reco::GenParticle* parton = it_jet->genParton();
        const reco::Candidate* genCand = parton;
        int qpTag(0) ; 
        while (genCand!=NULL && genCand->numberOfMothers()==1) {
          genCand = genCand->mother(0);
          if (abs(genCand->pdgId())==7 ) qpTag = 7;	// check if it's bprime
          else if (abs(genCand->pdgId())==8 ) qpTag = 8;	// check if it's tprime.
          else qpTag = 0 ; 
          if (abs(genCand->pdgId())==23) JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] = 23;
          if (abs(genCand->pdgId())==24) JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] = 24;
          if (abs(genCand->pdgId())==25) JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] = 25;
        }
        JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] += qpTag*100 ; 
      }
    } //// GenJet for MC 

    // in the case of fatjet, store the two subjets index 
    if (jettypes_[icoll] == "fatjet")  {

      int subjet1Idx=-1, subjet2Idx=-1;
      for( PatJetCollection::const_iterator jIt = jetsColl2->begin();
          jIt != jetsColl2->end(); ++jIt )
      {
        if( &(*jIt) == fatJetToPrunedFatJetMap.find(&(*it_jet))->second->daughter(0) )
          subjet1Idx = int( jIt - jetsColl2->begin() );
        if( &(*jIt) == fatJetToPrunedFatJetMap.find(&(*it_jet))->second->daughter(1) )
          subjet2Idx = int( jIt - jetsColl2->begin() );
        if( subjet1Idx>=0 && subjet2Idx>=0 ) break;
      }

      JetInfo[icoll].Jet_SubJet1Idx[JetInfo[icoll].Size] = subjet1Idx;
      JetInfo[icoll].Jet_SubJet2Idx[JetInfo[icoll].Size] = subjet2Idx;

      std::vector<fastjet::PseudoJet> fjConstituents;
      std::vector<edm::Ptr<reco::PFCandidate> > constituents = it_jet->getPFConstituents();
      std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator m;
      for ( m = constituents.begin(); m != constituents.end(); ++m ) {
        reco::PFCandidatePtr constit = *m;
        if (constit->pt() == 0) {
          edm::LogWarning("NullTransverseMomentum") << "dropping input candidate with pt=0";
          continue;
        } 
        fjConstituents.push_back(fastjet::PseudoJet(constit->px(),constit->py(),constit->pz(),constit->energy()));
        fjConstituents.back().set_user_index(m - constituents.begin());
      } //// Loop over fat jet constituents 
      JetInfo[icoll].tau1[JetInfo[icoll].Size] = nsubjettinessCalculator.getTau(1,fjConstituents);
      JetInfo[icoll].tau2[JetInfo[icoll].Size] = nsubjettinessCalculator.getTau(2,fjConstituents);
      JetInfo[icoll].tau3[JetInfo[icoll].Size] = nsubjettinessCalculator.getTau(3,fjConstituents);

      JetInfo[icoll].EtPruned      [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(&(*it_jet))->second->et();
      JetInfo[icoll].PtPruned      [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(&(*it_jet))->second->pt();
      JetInfo[icoll].EtaPruned     [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(&(*it_jet))->second->eta();
      JetInfo[icoll].PhiPruned     [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(&(*it_jet))->second->phi();
      JetInfo[icoll].EnergyPruned  [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(&(*it_jet))->second->energy(); 
      JetInfo[icoll].PxPruned      [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(&(*it_jet))->second->px(); 
      JetInfo[icoll].PyPruned      [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(&(*it_jet))->second->py(); 
      JetInfo[icoll].PzPruned      [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(&(*it_jet))->second->pz(); 
      JetInfo[icoll].MassPruned    [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(&(*it_jet))->second->mass();
      JetInfo[icoll].AreaPruned    [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(&(*it_jet))->second->jetArea();

    } //// If fat jets 

    // in the case of subjet, store the fatjet index
    if (jettypes_[icoll] == "subjet")  {
      int fatjetIdx=-1;
      for( PatJetCollection::const_iterator jIt = jetsColl2->begin(); jIt != jetsColl2->end(); ++jIt )
      {
        if( &(*it_jet) == fatJetToPrunedFatJetMap.find(&(*jIt))->second->daughter(0) ||
            &(*it_jet) == fatJetToPrunedFatJetMap.find(&(*jIt))->second->daughter(1) )
        {
          fatjetIdx = int( jIt - jetsColl2->begin() );
          break;
        }
      }

      JetInfo[icoll].Jet_FatJetIdx[JetInfo[icoll].Size] = fatjetIdx;

    } //// If subjets 

    JetInfo[icoll].Size++;

  } //loop over jets in collection

}


void BprimeTobH::saveHLT(const edm::Event& iEvent) {  
  // HLT: Booking trigger bits  
  edm::Handle<edm::TriggerResults> TrgResultsHandle;
  bool with_TriggerResults = (hltlabel_.size()>0) ? iEvent.getByLabel(hltlabel_[0],TrgResultsHandle) : false;
  if (with_TriggerResults) {

    const edm::TriggerNames & TrgNames = iEvent.triggerNames(*TrgResultsHandle);
    EvtInfo.TrgCount = 0;  
    for(int i=0; i<N_TRIGGER_BOOKINGS; i++) {

      unsigned int TrgIndex = TrgNames.triggerIndex( TriggerBooking[i] );

      if (TrgIndex == TrgNames.size()) {
        EvtInfo.TrgBook[i] = -4; // The trigger path is not known in this event.
      }else if ( !TrgResultsHandle->wasrun( TrgIndex ) ) {
        EvtInfo.TrgBook[i] = -3; // The trigger path was not included in this event.
      }else if ( !TrgResultsHandle->accept( TrgIndex ) ) {
        EvtInfo.TrgBook[i] = -2; // The trigger path was not accepted in this event.
      }else if (  TrgResultsHandle->error ( TrgIndex ) ) {
        EvtInfo.TrgBook[i] = -1; // The trigger path has an error in this event.
      }else {
        EvtInfo.TrgBook[i] = +1; // It's triggered.
        EvtInfo.TrgCount++; 
      }
    }
    EvtInfo.nHLT = TrgNames.size();
    for(unsigned int i=0; i<TrgNames.size();i++){
      EvtInfo.HLTbits[i] = (TrgResultsHandle->accept(i) == true) ? 1:0;
      // Print out Trigger table
      //std::cout << "trigger path= " << TrgNames.triggerName(i) << std::endl;
    }
  }
}

void BprimeTobH::saveL1T(const edm::Event& iEvent) { 
  //   L1 trigger and techincal trigger bits
  edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
  if(gtdigilabel_.size() > 0) iEvent.getByLabel( gtdigilabel_[0], gtRecord);

  if(gtRecord.isValid()) {
    DecisionWord dWord = gtRecord->decisionWord();
    if ( ! dWord.empty() ) { // if board not there this is zero
      // loop over dec. bit to get total rate (no overlap)
      for ( int i = 0; i < 128; ++i ) {
        //	if(dWord[i]!=0 && debug)cout << i << " " << dWord[i] << ": ";
        EvtInfo.L1[i]=dWord[i];
      }
    }
    TechnicalTriggerWord tw = gtRecord->technicalTriggerWord();
    if ( ! tw.empty() ) {
      // loop over dec. bit to get total rate (no overlap)
      for ( int i = 0; i < 64; ++i ) {
        //	if(tw[i]!=0 && debug) cout << i << "  " << tw[i] << ": ";
        EvtInfo.TT[i]=tw[i];
      }
    }
  }

}

void BprimeTobH::saveGenInfo(const edm::Event& iEvent) {

  memset(&GenInfo,0x00,sizeof(GenInfo));

  edm::Handle<reco::GenParticleCollection> GenHandle;  
  iEvent.getByLabel(genlabel_, GenHandle);

  vector<const reco::Candidate *> cands;
  vector<const reco::Candidate *>::const_iterator found = cands.begin();

  for( std::vector<reco::GenParticle>::const_iterator it_gen = GenHandle->begin(); it_gen != GenHandle->end(); it_gen++ ) { 
    if(it_gen->status() == 3){	
      int iMo1 = -1;
      int iMo2 = -1;
      int iDa1 = -1;
      int iDa2 = -1;
      int NMo = it_gen->numberOfMothers();
      int NDa = it_gen->numberOfDaughters();
      found = find(cands.begin(), cands.end(), it_gen->mother(0));
      if(found != cands.end()) iMo1 = found - cands.begin() ;

      found = find(cands.begin(), cands.end(), it_gen->mother(NMo-1));
      if(found != cands.end()) iMo2 = found - cands.begin() ;

      found = find(cands.begin(), cands.end(), it_gen->daughter(0));
      if(found != cands.end()) iDa1 = found - cands.begin() ;

      found = find(cands.begin(), cands.end(), it_gen->daughter(NDa-1));
      if(found != cands.end()) iDa2 = found - cands.begin() ;

      GenInfo.Pt[GenInfo.Size] 		= it_gen->pt();
      GenInfo.Eta[GenInfo.Size]	 	= it_gen->eta();
      GenInfo.Phi[GenInfo.Size]	 	= it_gen->phi();
      GenInfo.Mass[GenInfo.Size]		= it_gen->mass();
      GenInfo.PdgID[GenInfo.Size]		= it_gen->pdgId();
      GenInfo.Status[GenInfo.Size]	        = it_gen->status();

      GenInfo.nMo[GenInfo.Size]		= NMo; 
      GenInfo.nDa[GenInfo.Size]		= NDa; 
      GenInfo.Mo1[GenInfo.Size]		= iMo1; 
      GenInfo.Mo2[GenInfo.Size]		= iMo2; 
      GenInfo.Da1[GenInfo.Size]		= iDa1; 
      GenInfo.Da2[GenInfo.Size]		= iDa2; 

      ++GenInfo.Size ; 
    } //// Storing status == 3 particles only 
  } //// Looping over GenParticles

}


void BprimeTobH::processGenJets(const edm::Handle<GenJetCollection>& jetsColl, 
				const edm::Event& iEvent, 
				const edm::EventSetup& iSetup, 
				const unsigned int icoll) { 

  if(icoll >= MAX_JETCOLLECTIONS) return;

  memset(&JetInfo[icoll],0x00,sizeof(JetInfo[icoll]));

  for( vector<reco::GenJet>::const_iterator it_jet = jetsColl->begin();
       it_jet != jetsColl->end(); it_jet++ ) { 
    
    if (iEvent.isRealData() && it_jet->pt() < JetMinPt_) continue ;
     JetInfo[icoll].Index   [JetInfo[icoll].Size] = JetInfo[icoll].Size;
     // JetInfo[icoll].NTracks [JetInfo[icoll].Size] = it_jet->associatedTracks().size();

     JetInfo[icoll].Et      [JetInfo[icoll].Size] = it_jet->et();
     JetInfo[icoll].Pt      [JetInfo[icoll].Size] = it_jet->pt();
     JetInfo[icoll].Eta     [JetInfo[icoll].Size] = it_jet->eta();
     JetInfo[icoll].Phi     [JetInfo[icoll].Size] = it_jet->phi();
     JetInfo[icoll].Energy  [JetInfo[icoll].Size] = it_jet->energy(); 
     JetInfo[icoll].Px      [JetInfo[icoll].Size] = it_jet->px(); 
     JetInfo[icoll].Py      [JetInfo[icoll].Size] = it_jet->py(); 
     JetInfo[icoll].Pz      [JetInfo[icoll].Size] = it_jet->pz(); 
     JetInfo[icoll].Mass    [JetInfo[icoll].Size] = it_jet->mass();
     JetInfo[icoll].Area    [JetInfo[icoll].Size] = it_jet->jetArea();
      
     // JetInfo[icoll].JetCharge   [JetInfo[icoll].Size] = it_jet->jetCharge();
     JetInfo[icoll].NConstituents[JetInfo[icoll].Size] = it_jet->numberOfDaughters();

     JetInfo[icoll].Size++;
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(BprimeTobH);
