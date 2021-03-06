//-*- C++ -*-
//
// Package:BprimeTobH
// Class: BprimeTobH
//
/**\class BprimeTobH BprimeTobH.cc HbbAna/BprimeTobH/src/BprimeTobH.cc

Description: [one line class summary]

Implementation:
Based on bprimeKit
*/
//
// Original Author: Xin Shi <Xin.Shi@cern.ch>
// Created: Sat Apr 6 09:49:36 CEST 2013
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

//// PileupSummaryInfo
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//// For JEC
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h" 

#include "fastjet/PseudoJet.hh"
#include "../interface/Njettiness.hh"

// #include "QuarkGluonTagger/EightTeV/interface/QGTagger.h"

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
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, 
        edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, 
        edm::EventSetup const&);
    bool hasBeamSpot(const edm::Event&);
    void clearVariables();
    bool hasPrimaryVertex(const edm::Event &);
    bool hasMuons(const edm::Event &);
    bool hasElectrons(const edm::Event &);
    bool hasJets(const edm::Event &, const edm::EventSetup&);
    void saveGenInfo(const edm::Event &);
    void saveHLT(const edm::Event&);
    void saveL1T(const edm::Event&);
    void processJets(const edm::Handle<PatJetCollection>&, 
        const edm::Handle<PatJetCollection>&,
        const edm::Event&,
        const edm::EventSetup&, const JetToJetMap&,
        const unsigned int);
    void processGenJets(const edm::Handle<GenJetCollection>&,
        const edm::Event&, const edm::EventSetup&,
        const unsigned int);
    double calcJetY(pat::Jet) ;

    // ----------member data ---------------------------
    TTree* tree_;

    TH1F* h_events_;

    bool includeL7_;
    bool doElectrons_;

    edm::InputTag BeamSpotLabel_;
    edm::InputTag VertexLabel_;
    vector<edm::InputTag> muonlabel_;
    vector<edm::InputTag> electronlabel_;
    edm::InputTag jetlabel_;
    edm::InputTag fatjetlabel_;
    edm::InputTag prunedfatjetlabel_;
    edm::InputTag subjetlabel_;
    edm::InputTag genjetlabel_;
    vector<edm::InputTag> hltlabel_;
    vector<edm::InputTag> genevtlabel_;
    vector<edm::InputTag> gtdigilabel_;
    vector<edm::InputTag> rhocorrectionlabel_;
    vector<edm::InputTag> sigmaLabel_;
    vector<edm::InputTag> puInfoLabel_;
    edm::InputTag genlabel_;

    EvtInfoBranches EvtInfo;
    GenInfoBranches GenInfo;
    VertexInfoBranches VertexInfo;
    LepInfoBranches LepInfo             [MAX_LEPCOLLECTIONS];
    JetInfoBranches JetInfo             [MAX_JETCOLLECTIONS];
    JetInfoBranches JetInfo_JEShigh     [MAX_JETCOLLECTIONS];
    JetInfoBranches JetInfo_JESlow      [MAX_JETCOLLECTIONS];
    JetInfoBranches JetInfo_JERhigh     [MAX_JETCOLLECTIONS];
    JetInfoBranches JetInfo_JERlow      [MAX_JETCOLLECTIONS];
    JetInfoBranches JetInfo_btag_bc_high;
    JetInfoBranches JetInfo_btag_bc_low ;
    JetInfoBranches JetInfo_btag_l_high ;
    JetInfoBranches JetInfo_btag_l_low  ;
    JetInfoBranches JetInfo_Htag_bc_high;
    JetInfoBranches JetInfo_Htag_bc_low ;
    JetInfoBranches JetInfo_Htag_l_high ;
    JetInfoBranches JetInfo_Htag_l_low  ;

    // Across the event
    reco::BeamSpot beamSpot_;
    reco::Vertex primaryVertex_;

    vector<std::string> lepcollections_;
    vector<std::string> jetcollections_;
    vector<std::string> jettypes_;

    bool doGenInfo_;
    bool doGenJets_ ;

    bool doJESUncert_ ; 
    bool doJERUncert_ ; 
    bool dobtagUncert_ ; 
    bool doHtagUncert_ ; 

    double JetPtMin_;
    double JetYMax_ ; 
    double FatJetPtMin_;

    Njettiness nsubjettinessCalculator;

};

//
// static data member definitions
//

//
// constructors and destructor
//
BprimeTobH::BprimeTobH(const edm::ParameterSet& iConfig):
  tree_                  (0),
  h_events_              (0), 
  includeL7_             (iConfig.getUntrackedParameter<bool>("IncludeL7")),
  doElectrons_           (iConfig.getUntrackedParameter<bool>("DoElectrons")),
  BeamSpotLabel_         (iConfig.getParameter<edm::InputTag>("BeamSpotLabel")),
  VertexLabel_           (iConfig.getParameter<edm::InputTag>("VertexLabel")),
  muonlabel_             (iConfig.getParameter<vector<edm::InputTag> >("muonlabel")),
  electronlabel_         (iConfig.getParameter<vector<edm::InputTag> >("electronlabel")),
  jetlabel_              (iConfig.getParameter<edm::InputTag>("jetlabel")),
  fatjetlabel_           (iConfig.getParameter<edm::InputTag>("fatjetlabel")),
  prunedfatjetlabel_     (iConfig.getParameter<edm::InputTag>("prunedfatjetlabel")),
  subjetlabel_           (iConfig.getParameter<edm::InputTag>("subjetlabel")),
  genjetlabel_           (iConfig.getParameter<edm::InputTag>("genjetlabel")),
  hltlabel_              (iConfig.getParameter<vector<edm::InputTag> >("hltlabel")),
  genevtlabel_           (iConfig.getParameter<vector<edm::InputTag> >("genevtlabel")),
  gtdigilabel_           (iConfig.getParameter<vector<edm::InputTag> >("gtdigilabel")),
  rhocorrectionlabel_    (iConfig.getParameter<vector<edm::InputTag>>("rhocorrectionlabel")), 
  sigmaLabel_            (iConfig.getParameter<vector<edm::InputTag>>("sigmaLabel")),
  puInfoLabel_           (iConfig.getParameter<vector<edm::InputTag>>("puInfoLabel")),
  genlabel_              (iConfig.getParameter<edm::InputTag>("genlabel")),
  lepcollections_        (iConfig.getParameter<std::vector<std::string> >("LepCollections")),
  jetcollections_        (iConfig.getParameter<std::vector<std::string> >("JetCollections")),
  jettypes_              (iConfig.getParameter<std::vector<std::string> >("JetTypes")),
  doGenInfo_             (iConfig.getUntrackedParameter<bool>("DoGenInfo")),
  doGenJets_             (iConfig.getUntrackedParameter<bool>("DoGenJets")),
  doJESUncert_           (iConfig.getUntrackedParameter<bool>("DoJESUncert")),
  doJERUncert_           (iConfig.getUntrackedParameter<bool>("DoJERUncert")),
  dobtagUncert_          (iConfig.getUntrackedParameter<bool>("DobtagUncert")),
  doHtagUncert_          (iConfig.getUntrackedParameter<bool>("DoHtagUncert")),
  JetPtMin_              (iConfig.getUntrackedParameter<double>("JetPtMin")), 
  JetYMax_               (iConfig.getUntrackedParameter<double>("JetYMax")), 
  FatJetPtMin_           (iConfig.getUntrackedParameter<double>("FatJetPtMin")), 
  nsubjettinessCalculator(Njettiness::onepass_kt_axes, NsubParameters(1.0, 0.8, 0.8)) 
{

}

BprimeTobH::~BprimeTobH() {

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event ------------
void BprimeTobH::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  clearVariables();

  h_events_->Fill(0);   

  EvtInfo.RunNo  = iEvent.id().run();
  EvtInfo.EvtNo  = iEvent.id().event();
  EvtInfo.BxNo   = iEvent.bunchCrossing();
  EvtInfo.LumiNo = iEvent.luminosityBlock();
  EvtInfo.Orbit  = iEvent.orbitNumber();
  EvtInfo.McFlag = iEvent.isRealData()? 0: 1;  

  EvtInfo.nTrgBook = N_TRIGGER_BOOKINGS;

  std::vector<edm::Handle<double> > rhoHandles;
  std::vector<edm::Handle<double> > sigmaHandles;
  for(unsigned il=0; il<rhocorrectionlabel_.size(); il++) {
    rhoHandles.push_back(edm::Handle<double> ());
    iEvent.getByLabel( rhocorrectionlabel_[il],rhoHandles[il]);
    sigmaHandles.push_back(edm::Handle<double> ());
    iEvent.getByLabel( sigmaLabel_[il],sigmaHandles[il]);
  }

  for(unsigned int ii = 0; ii < 2; ii++) {  
    if(rhoHandles[ii].isValid())   EvtInfo.RhoPU[ii]   = *(rhoHandles[ii].product());
    if(sigmaHandles[ii].isValid()) EvtInfo.SigmaPU[ii] = *(sigmaHandles[ii].product());
  }

  if( iEvent.isRealData() == false ) { 

    //// Direct filling of PU weights: reserved for later use 
    //// edm::EventBase* iEventB = dynamic_cast<edm::EventBase*>(&iEvent);
    //// EvtInfo.WeightPU = LumiWeights_.weight( (*iEventB) );

    edm::Handle<std::vector<PileupSummaryInfo> >  puInfoHandles;
    if(puInfoLabel_.size() > 0) iEvent.getByLabel(puInfoLabel_[0], puInfoHandles); 
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = puInfoHandles->begin(); PVI != puInfoHandles->end(); ++PVI) {
      EvtInfo.nPU[EvtInfo.nBX]    = PVI->getPU_NumInteractions();
      EvtInfo.BXPU[EvtInfo.nBX]   = PVI->getBunchCrossing();
      EvtInfo.TrueIT[EvtInfo.nBX] = PVI->getTrueNumInteractions();
      EvtInfo.nBX += 1;
    }

    edm::Handle<GenEventInfoProduct> genEventInfoHandle;	 
    bool with_GenEventInfo = (genevtlabel_.size() > 0) ? iEvent.getByLabel( genevtlabel_[0], genEventInfoHandle ) : false ;  

    if (with_GenEventInfo) {
      if (genEventInfoHandle->hasPDF()) { 
        EvtInfo.PDFid1   = genEventInfoHandle->pdf()->id.first;
        EvtInfo.PDFid2   = genEventInfoHandle->pdf()->id.second;
        EvtInfo.PDFx1    = genEventInfoHandle->pdf()->x.first;
        EvtInfo.PDFx2    = genEventInfoHandle->pdf()->x.second;
        EvtInfo.PDFscale = genEventInfoHandle->pdf()->scalePDF;
        EvtInfo.PDFv1    = genEventInfoHandle->pdf()->xPDF.first;
        EvtInfo.PDFv2    = genEventInfoHandle->pdf()->xPDF.second;
      }
      EvtInfo.qScale   = genEventInfoHandle->qScale();
      EvtInfo.alphaQCD = genEventInfoHandle->alphaQCD();
      EvtInfo.alphaQED = genEventInfoHandle->alphaQED();
      EvtInfo.Weight   = genEventInfoHandle->weight();
    }

  }

  if ( hasBeamSpot(iEvent) && hasPrimaryVertex(iEvent) ) { 

    saveL1T(iEvent);
    saveHLT(iEvent);
    hasMuons(iEvent);
    if (doElectrons_) hasElectrons(iEvent);
    hasJets(iEvent, iSetup);
    if ( doGenInfo_ && !iEvent.isRealData() ) saveGenInfo(iEvent);

    tree_->Fill();

    h_events_->Fill(1);   

  }

  clearVariables();
}

// ------------ method called once each job just before starting event loop ------------
void BprimeTobH::beginJob() {

  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("results") );

  h_events_ = fs->make<TH1F>( "h_events", "Processed events", 2,  0, 2); 
  h_events_->Sumw2() ; 

  tree_ = new TTree ("tree", "BprimeTobH");
  EvtInfo.RegisterTree(tree_);
  VertexInfo.RegisterTree(tree_);

  if(lepcollections_.size() > MAX_LEPCOLLECTIONS) 
    edm::LogWarning("LepCollection") << "WARNING: Too many lep collections, using first " << MAX_LEPCOLLECTIONS ; 

  for(unsigned i=0; i<lepcollections_.size(); i++) { 
    if(i >= MAX_LEPCOLLECTIONS) break;
    LepInfo[i].RegisterTree(tree_,lepcollections_[i]);
  }

  if(jetcollections_.size() > MAX_JETCOLLECTIONS) 
    edm::LogWarning("LepCollection") << "WARNING: Too many jet collections, using first " << MAX_JETCOLLECTIONS ; 

  for(unsigned i=0; i<jetcollections_.size(); i++) {
    if(i >= MAX_JETCOLLECTIONS) break;
    JetInfo[i].RegisterTree(tree_,jetcollections_[i]);
  }

  if(doGenInfo_) GenInfo.RegisterTree(tree_);

}

// ------------ method called once each job just after ending the event loop ------------
void BprimeTobH::endJob() {
}

// ------------ method called when starting to processes a run ------------
void BprimeTobH::beginRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run ------------
void BprimeTobH::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block ------------
void BprimeTobH::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block ------------
void BprimeTobH::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module ------------
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

  if ( !recoVertexHandle.isValid()
      or recoVertexHandle.failedToGet()
      or recoVertexHandle->size() <= 0 )
    return false;

  bool gotPrimVtx = false;

  for (std::vector<reco::Vertex>::const_iterator iVertex = recoVertexHandle->begin();
      iVertex != recoVertexHandle->end(); iVertex++) {

    if (VertexInfo.Size>=MAX_VERTICES) {
      cout << " PV " << VertexInfo.Size << endl;
      fprintf(stderr,"ERROR: number of Tracks exceeds the size of array.\n");
      break;
    }

    VertexInfo.Type[VertexInfo.Size] = 1; //Vertices WITH the Beam Spot constraint
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

    if (!gotPrimVtx && (!iVertex->isFake() && iVertex->ndof()>=4. 
          && iVertex->z() <=24.
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
  edm::Handle <PatJetCollection> jetsColl; // = JetHandle[0];
  edm::Handle <PatJetCollection> fatjetsColl; // = JetHandle[0];
  edm::Handle <PatJetCollection> prunedfatjetsColl; // = JetHandle[1];
  edm::Handle <PatJetCollection> subjetsColl; // = JetHandle[2];
  edm::Handle <GenJetCollection> genjetsColl;

  iEvent.getByLabel( jetlabel_         , jetsColl);
  iEvent.getByLabel( fatjetlabel_      , fatjetsColl);
  iEvent.getByLabel( prunedfatjetlabel_, prunedfatjetsColl);
  iEvent.getByLabel( subjetlabel_      , subjetsColl);
  if ( iEvent.isRealData() == false ) iEvent.getByLabel( genjetlabel_      , genjetsColl);

  // get the fatjet to prunedfat jet match map
  JetToJetMap fatJetToPrunedFatJetMap;

  for(PatJetCollection::const_iterator it = fatjetsColl->begin(); 
      it != fatjetsColl->end(); ++it) {

    double fatjety = calcJetY(*it) ; 
    if (it->pt() < FatJetPtMin_ || fabs(fatjety) > JetYMax_) continue ; 

    PatJetCollection::const_iterator prunedJetMatch;
    bool prunedJetMatchFound = false;
    float dR = 0.8; //// hard coded for now 
    for(PatJetCollection::const_iterator pjIt = prunedfatjetsColl->begin(); 
        pjIt != prunedfatjetsColl->end(); ++pjIt) { 
      float dR_temp = reco::deltaR( it->p4(), pjIt->p4() ); 
      if( dR_temp < dR ) { 
        prunedJetMatchFound = true;
        dR = dR_temp;
        prunedJetMatch = pjIt;
      }
    }
    if( !prunedJetMatchFound )
      edm::LogError("NoMatchingGroomedJet") << 
        " Matching pruned jet not found" << endl; 
    else fatJetToPrunedFatJetMap[&(*it)] = &(*prunedJetMatch); 

  } // Loop over fat jets 

  // Now process 'FatJetInfo', 'SubJetInfo', 'JetInfo', 'GenJetInfo':
  unsigned int iJetColl ;

  iJetColl = 0 ; // FatJetInfo
  processJets(fatjetsColl, subjetsColl, iEvent, iSetup, fatJetToPrunedFatJetMap, iJetColl) ;

  iJetColl = 1; //// SubJetInfo
  processJets(subjetsColl, fatjetsColl, iEvent, iSetup, fatJetToPrunedFatJetMap, iJetColl) ;

  iJetColl = 2 ; //// JetInfo
  processJets(jetsColl, subjetsColl, iEvent, iSetup, fatJetToPrunedFatJetMap, iJetColl) ;

  if ( iEvent.isRealData() == false ) {
    iJetColl = 3 ; // GenJetInfo
    processGenJets(genjetsColl, iEvent, iSetup, iJetColl) ;
  }

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
  if (icoll == 0) iSetup.get<JetCorrectionsRecord>().get("AK7PFchs", JetCorParColl);
  else if (icoll == 1 || icoll == 2) iSetup.get<JetCorrectionsRecord>().get("AK5PFchs", JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

  for( vector<pat::Jet>::const_iterator it_jet = jetsColl->begin();
      it_jet != jetsColl->end(); it_jet++ ) {

    double jety = calcJetY(*it_jet) ; 
    if (it_jet->pt() < JetPtMin_  || fabs(jety) > JetYMax_) continue ;

    // avoid no match found jets
    const pat::Jet* p_jet = &(*it_jet);     
    if (jettypes_[icoll] == "fatjet" && 
        fatJetToPrunedFatJetMap.find(p_jet) == fatJetToPrunedFatJetMap.end() )
      continue; 

    JetInfo[icoll].Index [JetInfo[icoll].Size] = JetInfo[icoll].Size;
    JetInfo[icoll].NTracks [JetInfo[icoll].Size] = it_jet->associatedTracks().size();

    JetInfo[icoll].Et [JetInfo[icoll].Size] = it_jet->et();
    JetInfo[icoll].Pt [JetInfo[icoll].Size] = it_jet->pt();
    JetInfo[icoll].Eta [JetInfo[icoll].Size] = it_jet->eta();
    JetInfo[icoll].Phi [JetInfo[icoll].Size] = it_jet->phi();
    JetInfo[icoll].Energy [JetInfo[icoll].Size] = it_jet->energy();
    JetInfo[icoll].Px [JetInfo[icoll].Size] = it_jet->px();
    JetInfo[icoll].Py [JetInfo[icoll].Size] = it_jet->py();
    JetInfo[icoll].Pz [JetInfo[icoll].Size] = it_jet->pz();
    JetInfo[icoll].Mass [JetInfo[icoll].Size] = it_jet->mass();
    JetInfo[icoll].Area [JetInfo[icoll].Size] = it_jet->jetArea();

    jecUnc->setJetEta(it_jet->eta());
    jecUnc->setJetPt(it_jet->pt()); // here you must use the CORRECTED jet pt
    if(fabs(it_jet->eta())<=5.0) JetInfo[icoll].Unc [JetInfo[icoll].Size] = jecUnc->getUncertainty(true);

    JetInfo[icoll].JetCharge [JetInfo[icoll].Size] = it_jet->jetCharge();
    JetInfo[icoll].NConstituents[JetInfo[icoll].Size] = it_jet->numberOfDaughters();

    if (it_jet->isPFJet()) {
      JetInfo[icoll].NCH[JetInfo[icoll].Size] = it_jet->chargedMultiplicity();
      JetInfo[icoll].CEF[JetInfo[icoll].Size] = it_jet->chargedEmEnergyFraction();
      JetInfo[icoll].NHF[JetInfo[icoll].Size] = it_jet->neutralHadronEnergyFraction();
      JetInfo[icoll].NEF[JetInfo[icoll].Size] = it_jet->neutralEmEnergyFraction();
      JetInfo[icoll].CHF[JetInfo[icoll].Size] = it_jet->chargedHadronEnergyFraction();
    }

    bool JetIDLoose(false);
    bool JetIDTight(false);

    JetInfo[icoll].QGTagsMLP [JetInfo[icoll].Size] = -999;
    JetInfo[icoll].QGTagsLikelihood [JetInfo[icoll].Size] = -1;

    //// Gluon tagger
    // edm::Handle<edm::ValueMap<float> > QGTagsHandleMLP;
    // edm::Handle<edm::ValueMap<float> > QGTagsHandleLikelihood;
    // iEvent.getByLabel("QGTagger","qgMLP", QGTagsHandleMLP);
    // //iEvent.getByLabel("QGTagger","qgLikelihood", QGTagsHandleLikelihood);
    // iEvent.getByLabel("QGTagger", QGTagsHandleLikelihood);

    int ijet = it_jet - jetsColl->begin();
    edm::RefToBase<reco::Jet> jetRef(edm::Ref<std::vector <pat::Jet> >(jetsColl,ijet));

    // if (QGTagsHandleMLP.isValid()){
    //   // std::cout << "QGTagsHandleMLP is Valid\n" ;
    //   JetInfo[icoll].QGTagsMLP [JetInfo[icoll].Size] = (*QGTagsHandleMLP)[jetRef];
    // }
    // // else std::cout << "QGTagsHandleMLP is not Valid\n" ;
    // if (QGTagsHandleLikelihood.isValid()){
    //   // std::cout << "QGTagsHandleLikelihood is Valid\n" ;
    //   JetInfo[icoll].QGTagsLikelihood [JetInfo[icoll].Size] = (*QGTagsHandleLikelihood)[jetRef];
    // }
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

    JetInfo[icoll].JetIDLOOSE[JetInfo[icoll].Size] = (JetIDLoose) ? 1 : 0;
    JetInfo[icoll].JetIDTIGHT[JetInfo[icoll].Size] = (JetIDTight) ? 1 : 0;

    // Jet corrections, B-tagging, and Jet ID information
    // now we just fill everything (regardless of availability)
    JetInfo[icoll].PtCorrRaw [JetInfo[icoll].Size] = it_jet->correctedJet("Uncorrected").pt();  
    JetInfo[icoll].PtCorrL2[JetInfo[icoll].Size] = it_jet->correctedJet("L2Relative" ).pt(); // L2(rel)
    JetInfo[icoll].PtCorrL3[JetInfo[icoll].Size] = it_jet->correctedJet("L3Absolute" ).pt(); // L3(abs)
    if(includeL7_) {
      JetInfo[icoll].PtCorrL7g [JetInfo[icoll].Size] = it_jet->correctedJet("L7Parton", "gluon" ).pt(); // L7(gluon)
      JetInfo[icoll].PtCorrL7uds [JetInfo[icoll].Size] = it_jet->correctedJet("L7Parton", "uds" ).pt(); // L7(uds-jet)
      JetInfo[icoll].PtCorrL7c [JetInfo[icoll].Size] = it_jet->correctedJet("L7Parton", "charm" ).pt(); // L7(c-jet)
      JetInfo[icoll].PtCorrL7b [JetInfo[icoll].Size] = it_jet->correctedJet("L7Parton", "bottom").pt(); // L7(b-jet)
    }

    JetInfo[icoll].JetBProbBJetTags[JetInfo[icoll].Size] = it_jet->bDiscriminator("jetBProbabilityBJetTags");
    JetInfo[icoll].JetProbBJetTags [JetInfo[icoll].Size] = it_jet->bDiscriminator("jetProbabilityBJetTags");
    JetInfo[icoll].TrackCountHiPurBJetTags [JetInfo[icoll].Size] = it_jet->bDiscriminator("trackCountingHighPurBJetTags");
    JetInfo[icoll].CombinedSVBJetTags [JetInfo[icoll].Size] = it_jet->bDiscriminator("combinedSecondaryVertexBJetTags");
    JetInfo[icoll].CombinedSVMVABJetTags [JetInfo[icoll].Size] = it_jet->bDiscriminator("combinedSecondaryVertexMVABJetTags");
    JetInfo[icoll].SoftElecByIP3dBJetTags [JetInfo[icoll].Size] = it_jet->bDiscriminator("softElectronByIP3dBJetTags");
    JetInfo[icoll].SoftElecByPtBJetTags[JetInfo[icoll].Size] = it_jet->bDiscriminator("softElectronByPtBJetTags");
    JetInfo[icoll].SoftMuonBJetTags[JetInfo[icoll].Size] = it_jet->bDiscriminator("softMuonBJetTags");
    JetInfo[icoll].SoftMuonByIP3dBJetTags [JetInfo[icoll].Size] = it_jet->bDiscriminator("softMuonByIP3dBJetTags");
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
          if (abs(genCand->pdgId())==7 ) qpTag = 7; // check if it's bprime
          else if (abs(genCand->pdgId())==8 ) qpTag = 8;  // check if it's tprime.
          else qpTag = 0 ;
          if (abs(genCand->pdgId())==23) JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] = 23;
          if (abs(genCand->pdgId())==24) JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] = 24;
          if (abs(genCand->pdgId())==25) JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] = 25;
        }
        JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] += qpTag*100 ;

        // const reco::GenParticle* parton = it_jet->genParton();
        if (parton!=NULL) 
        {
          JetInfo[icoll].GenPt     [JetInfo[icoll].Size] = parton->pt();
          JetInfo[icoll].GenEta    [JetInfo[icoll].Size] = parton->eta();
          JetInfo[icoll].GenPhi    [JetInfo[icoll].Size] = parton->phi();
          JetInfo[icoll].GenPdgID  [JetInfo[icoll].Size] = parton->pdgId();
          JetInfo[icoll].GenFlavor [JetInfo[icoll].Size] = it_jet->partonFlavour();

          const reco::Candidate* genCand = parton;

          int bprime_tag = 0;	// 0: not b' or t'; 1: b'; 2:t'
          while(genCand!=NULL && genCand->numberOfMothers()==1) 
          {
            genCand = genCand->mother(0);
            if (abs(genCand->pdgId())==7 ) bprime_tag = 1;	// check if it's bprime
            if (abs(genCand->pdgId())==8 ) bprime_tag = 2;	// check if it's tprime.
            if (abs(genCand->pdgId())==23) JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] = 2;
            if (abs(genCand->pdgId())==24) JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] = 1;
          }
          if (bprime_tag==1) JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] += 10;
          if (bprime_tag==2) JetInfo[icoll].GenMCTag[JetInfo[icoll].Size] += 20;
        }
      }
    } //// GenJet for MC



    // in the case of fatjet, store the two subjets index
    if (jettypes_[icoll] == "fatjet") {

      // if ( fatJetToPrunedFatJetMap.find(p_jet) != fatJetToPrunedFatJetMap.end() ) { 
      JetInfo[icoll].EtPruned [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(p_jet)->second->et();
      JetInfo[icoll].PtPruned [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(p_jet)->second->pt();
      JetInfo[icoll].EtaPruned [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(p_jet)->second->eta();
      JetInfo[icoll].PhiPruned [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(p_jet)->second->phi();
      JetInfo[icoll].EnergyPruned [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(p_jet)->second->energy();
      JetInfo[icoll].PxPruned [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(p_jet)->second->px();
      JetInfo[icoll].PyPruned [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(p_jet)->second->py();
      JetInfo[icoll].PzPruned [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(p_jet)->second->pz();
      JetInfo[icoll].MassPruned [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(p_jet)->second->mass();
      JetInfo[icoll].AreaPruned [JetInfo[icoll].Size] = fatJetToPrunedFatJetMap.find(p_jet)->second->jetArea();

      int subjet1Idx=-1, subjet2Idx=-1;
      for( PatJetCollection::const_iterator jIt = jetsColl2->begin(); jIt != jetsColl2->end(); ++jIt ) { 
        if( &(*jIt) == fatJetToPrunedFatJetMap.find(p_jet)->second->daughter(0) )
          subjet1Idx = int( jIt - jetsColl2->begin() );
        if( &(*jIt) == fatJetToPrunedFatJetMap.find(p_jet)->second->daughter(1) )
          subjet2Idx = int( jIt - jetsColl2->begin() );
        if( subjet1Idx>=0 && subjet2Idx>=0 ) break;
      }

      JetInfo[icoll].Jet_SubJet1Idx[JetInfo[icoll].Size] = subjet1Idx;
      JetInfo[icoll].Jet_SubJet2Idx[JetInfo[icoll].Size] = subjet2Idx;

      // } 

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

    } //// If fat jets

    // in the case of subjet, store the fatjet index

    if (jettypes_[icoll] == "subjet") {
      int fatjetIdx=-1;
      for( PatJetCollection::const_iterator jIt = jetsColl2->begin(); jIt != jetsColl2->end(); ++jIt )
      {
        if ( fatJetToPrunedFatJetMap.find(&(*jIt)) == fatJetToPrunedFatJetMap.end() ) 
          continue; 

        if( p_jet == fatJetToPrunedFatJetMap.find(&(*jIt))->second->daughter(0) ||
            p_jet == fatJetToPrunedFatJetMap.find(&(*jIt))->second->daughter(1) )
        {
          fatjetIdx = int( jIt - jetsColl2->begin() );
          break;
        }
        JetInfo[icoll].Jet_FatJetIdx[JetInfo[icoll].Size] = fatjetIdx;
      }
    } //// If subjets

    JetInfo[icoll].Size++;

  } //loop over jets in collection

  delete jecUnc; // thanks to Jacky!  

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
      }else if ( TrgResultsHandle->error ( TrgIndex ) ) {
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
  // L1 trigger and techincal trigger bits
  edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
  if(gtdigilabel_.size() > 0) iEvent.getByLabel( gtdigilabel_[0], gtRecord);

  if(gtRecord.isValid()) {
    DecisionWord dWord = gtRecord->decisionWord();
    if ( ! dWord.empty() ) { // if board not there this is zero
      // loop over dec. bit to get total rate (no overlap)
      for ( int i = 0; i < 128; ++i ) {
        // if(dWord[i]!=0 && debug)cout << i << " " << dWord[i] << ": ";
        EvtInfo.L1[i]=dWord[i];
      }
    }
    TechnicalTriggerWord tw = gtRecord->technicalTriggerWord();
    if ( ! tw.empty() ) {
      // loop over dec. bit to get total rate (no overlap)
      for ( int i = 0; i < 64; ++i ) {
        // if(tw[i]!=0 && debug) cout << i << " " << tw[i] << ": ";
        EvtInfo.TT[i]=tw[i];
      }
    }
  }

}

void BprimeTobH::saveGenInfo(const edm::Event& iEvent) {

  memset(&GenInfo,0x00,sizeof(GenInfo));

  edm::Handle<reco::GenParticleCollection> GenHandle;
  iEvent.getByLabel(genlabel_, GenHandle);

  for( std::vector<reco::GenParticle>::const_iterator it_gen = GenHandle->begin(); it_gen != GenHandle->end(); it_gen++ ) {
    if(it_gen->status() == 3){  

      GenInfo.Pt[GenInfo.Size] = it_gen->pt();
      GenInfo.Eta[GenInfo.Size] = it_gen->eta();
      GenInfo.Phi[GenInfo.Size] = it_gen->phi();
      GenInfo.Mass[GenInfo.Size]  = it_gen->mass();
      GenInfo.PdgID[GenInfo.Size] = it_gen->pdgId();
      GenInfo.Status[GenInfo.Size]  = it_gen->status();

      GenInfo.nDa[GenInfo.Size] = it_gen->numberOfDaughters(); 
      GenInfo.nMo[GenInfo.Size] = it_gen->numberOfMothers();

      int NDa = it_gen->numberOfDaughters();
      int NMo = it_gen->numberOfMothers();

      if (NDa == 1) {
        const reco::Candidate* dau = it_gen->daughter(0);
        GenInfo.Da0Pt    [GenInfo.Size] = dau->pt();
        GenInfo.Da0Eta   [GenInfo.Size] = dau->eta();
        GenInfo.Da0Phi   [GenInfo.Size] = dau->phi();
        GenInfo.Da0Mass  [GenInfo.Size] = dau->mass();
        GenInfo.Da0PdgID [GenInfo.Size] = dau->pdgId();
        GenInfo.Da0Status[GenInfo.Size] = dau->status();

        GenInfo.Da1Pt    [GenInfo.Size] = -1 ; 
        GenInfo.Da1Eta   [GenInfo.Size] = -1 ; 
        GenInfo.Da1Phi   [GenInfo.Size] = -1 ; 
        GenInfo.Da1Mass  [GenInfo.Size] = -1 ; 
        GenInfo.Da1PdgID [GenInfo.Size] = -1 ; 
        GenInfo.Da1Status[GenInfo.Size] = -1 ; 
      }
      else if (NDa >= 2) {
        const reco::Candidate* daus[2] = { it_gen->daughter(0), it_gen->daughter(1) };
        GenInfo.Da0Pt    [GenInfo.Size] = daus[0]->pt();
        GenInfo.Da0Eta   [GenInfo.Size] = daus[0]->eta();
        GenInfo.Da0Phi   [GenInfo.Size] = daus[0]->phi();
        GenInfo.Da0Mass  [GenInfo.Size] = daus[0]->mass();
        GenInfo.Da0PdgID [GenInfo.Size] = daus[0]->pdgId();
        GenInfo.Da0Status[GenInfo.Size] = daus[0]->status();

        GenInfo.Da1Pt    [GenInfo.Size] = daus[1]->pt();
        GenInfo.Da1Eta   [GenInfo.Size] = daus[1]->eta();
        GenInfo.Da1Phi   [GenInfo.Size] = daus[1]->phi();
        GenInfo.Da1Mass  [GenInfo.Size] = daus[1]->mass();
        GenInfo.Da1PdgID [GenInfo.Size] = daus[1]->pdgId();
        GenInfo.Da1Status[GenInfo.Size] = daus[1]->status();
      }

      if (NMo == 1) {
        const reco::Candidate* mom = it_gen->mother(0) ; 
        GenInfo.Mo0Pt    [GenInfo.Size] = mom->pt() ; 
        GenInfo.Mo0Eta   [GenInfo.Size] = mom->eta();
        GenInfo.Mo0Phi   [GenInfo.Size] = mom->phi();
        GenInfo.Mo0Mass  [GenInfo.Size] = mom->mass();
        GenInfo.Mo0PdgID [GenInfo.Size] = mom->pdgId();
        GenInfo.Mo0Status[GenInfo.Size] = mom->status();

        GenInfo.Mo1Pt    [GenInfo.Size] = -1 ; 
        GenInfo.Mo1Eta   [GenInfo.Size] = -1 ; 
        GenInfo.Mo1Phi   [GenInfo.Size] = -1 ; 
        GenInfo.Mo1Mass  [GenInfo.Size] = -1 ; 
        GenInfo.Mo1PdgID [GenInfo.Size] = -1 ; 
        GenInfo.Mo1Status[GenInfo.Size] = -1 ; 
      }
      else if (NMo >= 2) { 
        const reco::Candidate* moms[2] = { it_gen->mother(0), it_gen->mother(1) }; 
        GenInfo.Mo0Eta   [GenInfo.Size] = moms[0]->eta();
        GenInfo.Mo0Phi   [GenInfo.Size] = moms[0]->phi();
        GenInfo.Mo0Mass  [GenInfo.Size] = moms[0]->mass();
        GenInfo.Mo0PdgID [GenInfo.Size] = moms[0]->pdgId();
        GenInfo.Mo0Status[GenInfo.Size] = moms[0]->status();

        GenInfo.Mo1Pt    [GenInfo.Size] = moms[1]->pt();
        GenInfo.Mo1Eta   [GenInfo.Size] = moms[1]->eta();
        GenInfo.Mo1Phi   [GenInfo.Size] = moms[1]->phi();
        GenInfo.Mo1Mass  [GenInfo.Size] = moms[1]->mass();
        GenInfo.Mo1PdgID [GenInfo.Size] = moms[1]->pdgId();
        GenInfo.Mo1Status[GenInfo.Size] = moms[1]->status();
      }

      ++GenInfo.Size; 

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

    if (iEvent.isRealData() && it_jet->pt() < JetPtMin_) continue ;
    JetInfo[icoll].Index [JetInfo[icoll].Size] = JetInfo[icoll].Size;
    // JetInfo[icoll].NTracks [JetInfo[icoll].Size] = it_jet->associatedTracks().size();

    JetInfo[icoll].Et [JetInfo[icoll].Size] = it_jet->et();
    JetInfo[icoll].Pt [JetInfo[icoll].Size] = it_jet->pt();
    JetInfo[icoll].Eta [JetInfo[icoll].Size] = it_jet->eta();
    JetInfo[icoll].Phi [JetInfo[icoll].Size] = it_jet->phi();
    JetInfo[icoll].Energy [JetInfo[icoll].Size] = it_jet->energy();
    JetInfo[icoll].Px [JetInfo[icoll].Size] = it_jet->px();
    JetInfo[icoll].Py [JetInfo[icoll].Size] = it_jet->py();
    JetInfo[icoll].Pz [JetInfo[icoll].Size] = it_jet->pz();
    JetInfo[icoll].Mass [JetInfo[icoll].Size] = it_jet->mass();
    JetInfo[icoll].Area [JetInfo[icoll].Size] = it_jet->jetArea();

    // JetInfo[icoll].JetCharge [JetInfo[icoll].Size] = it_jet->jetCharge();
    JetInfo[icoll].NConstituents[JetInfo[icoll].Size] = it_jet->numberOfDaughters();

    JetInfo[icoll].Size++;
  }
}

double BprimeTobH::calcJetY(pat::Jet it) {
  return TMath::Log( (TMath::Sqrt( (it.mass()*it.mass()) + (it.pt()*it.pt())*
          (TMath::CosH(it.eta())*TMath::CosH(it.eta())) )
        + (it.pt()*TMath::SinH(it.eta()) ) ) / 
      (TMath::Sqrt( (it.mass()*it.mass())
                    + (it.pt()*it.pt()) )) ) ; 

}

//define this as a plug-in
DEFINE_FWK_MODULE(BprimeTobH);

