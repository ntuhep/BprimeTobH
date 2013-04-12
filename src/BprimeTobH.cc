// -*- C++ -*-
//
// Package:    BprimeTobH
// Class:      BprimeTobH
// 
/**\class BprimeTobH BprimeTobH.cc HbbAna/BprimeTobH/src/BprimeTobH.cc

 Description: [one line class summary]

 Implementation:
     Based on bprimeKit 
*/
//
// Original Author:  Xin Shi <Xin.Shi@cern.ch>
//         Created:  Sat Apr  6 09:49:36 CEST 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "../interface/format.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"


using namespace std;

const unsigned int MAX_LEPCOLLECTIONS=3; 
const unsigned int MAX_JETCOLLECTIONS=3; 

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
  bool hasJets(const edm::Event &); 

  
  // ----------member data ---------------------------
  TTree* tree_;  
  edm::InputTag BeamSpotLabel_;
  edm::InputTag VertexLabel_;
  vector<edm::InputTag> muonlabel_;
  vector<edm::InputTag> electronlabel_;
  vector<edm::InputTag> jetlabel_;

  EvtInfoBranches EvtInfo;
  VertexInfoBranches VertexInfo;
  LepInfoBranches LepInfo[MAX_LEPCOLLECTIONS];
  JetInfoBranches JetInfo[MAX_JETCOLLECTIONS];

  // Across the event 
  reco::BeamSpot beamSpot_;  
  reco::Vertex primaryVertex_;
  reco::Vertex primaryVertexBS_;
  
  vector<std::string> lepcollections_;
  vector<std::string> jetcollections_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BprimeTobH::BprimeTobH(const edm::ParameterSet& iConfig):
  tree_(0), 
  BeamSpotLabel_(iConfig.getParameter<edm::InputTag>("BeamSpotLabel")),
  VertexLabel_(iConfig.getParameter<edm::InputTag>("VertexLabel")), 
  muonlabel_(iConfig.getParameter<vector<edm::InputTag> >("muonlabel")), 
  electronlabel_(iConfig.getParameter<vector<edm::InputTag> >("electronlabel")),  
  jetlabel_(iConfig.getParameter<vector<edm::InputTag> >("jetlabel")),  
  lepcollections_(iConfig.getParameter<std::vector<std::string> >("LepCollections"))
{
  edm::Service<TFileService> fs;
  TFileDirectory results = TFileDirectory( fs->mkdir("results") );


}


BprimeTobH::~BprimeTobH()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BprimeTobH::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  clearVariables(); 
  
  EvtInfo.RunNo = iEvent.id().run();
  EvtInfo.EvtNo = iEvent.id().event();

  if ( hasBeamSpot(iEvent) 
       && hasPrimaryVertex(iEvent)
       && hasPrimaryVertexBS(iEvent)
              ) {

    hasMuons(iEvent); 
    hasElectrons(iEvent); 
    hasJets(iEvent); 
      
    tree_->Fill();
  }

  clearVariables(); 
}

// ------------ method called once each job just before starting event loop  ------------
void 
BprimeTobH::beginJob()
{
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

}

// ------------ method called once each job just after ending the event loop  ------------
void 
BprimeTobH::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
BprimeTobH::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
BprimeTobH::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
BprimeTobH::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
BprimeTobH::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BprimeTobH::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool
BprimeTobH::hasBeamSpot(const edm::Event& iEvent)
{
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

void 
BprimeTobH::clearVariables(){
  memset(&EvtInfo,0x00,sizeof(EvtInfo));
  memset(&VertexInfo,0x00,sizeof(VertexInfo));
}

bool
BprimeTobH::hasPrimaryVertex(const edm::Event& iEvent)
{
  edm::Handle<reco::VertexCollection> recoVertexHandle;
  iEvent.getByLabel(VertexLabel_, recoVertexHandle);

  if ( ! recoVertexHandle.isValid() 
       or recoVertexHandle.failedToGet() 
       or recoVertexHandle->size() <= 0 )
    return false;

  // nprivtx = recoVertexHandle->size(); 
  bool gotPrimVtx = false;

  for (std::vector<reco::Vertex>::const_iterator iVertex = recoVertexHandle->begin();
       iVertex != recoVertexHandle->end(); iVertex++) { 
    
    if (VertexInfo.Size>=MAX_VERTICES) {
      cout << "PV " << VertexInfo.Size << endl;
      fprintf(stderr,"ERROR: number of  Tracks exceeds the size of array.\n");
      break;//exit(0);
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
    
    for (reco::Vertex::trackRef_iterator it = iVertex->tracks_begin(); 
	 it != iVertex->tracks_end(); it++) {
      VertexInfo.Pt_Sum[VertexInfo.Size] += (*it)->pt();
      VertexInfo.Pt_Sum2[VertexInfo.Size] += ((*it)->pt() * (*it)->pt());
    }

    if (!gotPrimVtx && (!iVertex->isFake() && iVertex->ndof()>=4. && iVertex->z() <=24. && iVertex->position().Rho()<=2.)) {
      primaryVertex_ = *(iVertex); 
      gotPrimVtx=true;
    }

    VertexInfo.Size++;
  }

  if (!primaryVertex_.isValid()) {
    // edm::LogError("myVertex") << "No valid vertex found!." ;
    return false; 
  }
 
  // edm::LogInfo("myVertex") << "Found number of vertex: " << recVtxs->size(); 
  
  return true; 
}


bool
BprimeTobH::hasPrimaryVertexBS(const edm::Event& iEvent)
{
  edm::Handle<reco::VertexCollection>  VertexHandleBS; 
  double PVBS_Pt_Max = -100.;

  if( ! VertexHandleBS.isValid() or VertexHandleBS.failedToGet() 
      or VertexHandleBS->size() <= 0 )
    return false;

  const vector<reco::Vertex> VerticesBS = *VertexHandleBS;
  for(vector<reco::Vertex>::const_iterator it_vtx = VerticesBS.begin();
      it_vtx != VerticesBS.end(); it_vtx++ ) {
    if (VertexInfo.Size>=MAX_VERTICES) {
      cout << "PVBS " << VertexInfo.Size << endl;
      fprintf(stderr,"ERROR: number of  Vertices exceeds the size of array.\n");
      break;//exit(0);
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
    for (reco::Vertex::trackRef_iterator it = it_vtx->tracks_begin(); it != it_vtx->tracks_end(); it++) {
      VertexInfo.Pt_Sum[VertexInfo.Size] += (*it)->pt();
      VertexInfo.Pt_Sum2[VertexInfo.Size] += ((*it)->pt() * (*it)->pt());
    }
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


bool
BprimeTobH::hasMuons(const edm::Event& iEvent)
{
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

bool
BprimeTobH::hasElectrons(const edm::Event& iEvent)
{
  vector<edm::Handle<vector<pat::Electron> > > ElectronHandle;
  for(unsigned il=0; il<electronlabel_.size(); il++) {
    ElectronHandle.push_back(edm::Handle<vector<pat::Electron> >());
    iEvent.getByLabel( electronlabel_[il], ElectronHandle[il]);
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

bool
BprimeTobH::hasJets(const edm::Event& iEvent)
{
  vector<edm::Handle<vector<pat::Jet> > > JetHandle;
  for(unsigned il=0; il<jetlabel_.size(); il++) {
    JetHandle.push_back(edm::Handle<vector<pat::Jet> >());
    iEvent.getByLabel( jetlabel_[il], JetHandle[il]);
  }
  
  for(unsigned icoll=0; icoll<jetcollections_.size(); icoll++) {
    //loop over collections
    if(icoll >= MAX_JETCOLLECTIONS) break;
    
    memset(&JetInfo[icoll],0x00,sizeof(JetInfo[icoll]));
    
    if (JetHandle.size() <= icoll) continue;  

    for( vector<pat::Jet>::const_iterator it_el = JetHandle[icoll]->begin(); 
	 it_el != JetHandle[icoll]->end(); it_el++ ) {
      
      JetInfo[icoll].Index[JetInfo[icoll].Size] = JetInfo[icoll].Size;
    }

  }

  return true; 
}



//define this as a plug-in
DEFINE_FWK_MODULE(BprimeTobH);
