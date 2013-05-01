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
		bool hasJets(const edm::Event &, const edm::EventSetup&); 
		void saveHLT(const edm::Event&);
		void saveL1T(const edm::Event&);


		// ----------member data ---------------------------
		TTree* tree_;  

		bool includeL7_;
		bool doSubjets_; 

		edm::InputTag BeamSpotLabel_;
		edm::InputTag VertexLabel_;
		edm::InputTag VertexBSLabel_;
		vector<edm::InputTag> muonlabel_;
		vector<edm::InputTag> electronlabel_;
		vector<edm::InputTag> jetlabel_;
		vector<edm::InputTag> hltlabel_;
		vector<edm::InputTag> gtdigilabel_;

		EvtInfoBranches EvtInfo;
		VertexInfoBranches VertexInfo;
		LepInfoBranches LepInfo[MAX_LEPCOLLECTIONS];
		JetInfoBranches JetInfo[MAX_JETCOLLECTIONS];
		JetInfoBranches SubJetInfo[MAX_JETCOLLECTIONS];

		// Across the event 
		reco::BeamSpot beamSpot_;  
		reco::Vertex primaryVertex_;
		reco::Vertex primaryVertexBS_;

		vector<std::string> lepcollections_;
		vector<std::string> jetcollections_;
		vector<std::string> jettypes_;

		bool doGenJets_ ; 


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
	includeL7_(iConfig.getUntrackedParameter<bool>("IncludeL7",false)), 
	doSubjets_(iConfig.getUntrackedParameter<bool>("DoSubjets",true)), 
	BeamSpotLabel_(iConfig.getParameter<edm::InputTag>("BeamSpotLabel")),
	VertexLabel_(iConfig.getParameter<edm::InputTag>("VertexLabel")), 
	VertexBSLabel_(iConfig.getParameter<edm::InputTag>("VertexBSLabel")), 
	muonlabel_(iConfig.getParameter<vector<edm::InputTag> >("muonlabel")), 
	electronlabel_(iConfig.getParameter<vector<edm::InputTag> >("electronlabel")),  
	jetlabel_(iConfig.getParameter<vector<edm::InputTag> >("jetlabel")),  
	hltlabel_(iConfig.getParameter<vector<edm::InputTag> >("hltlabel")),  
	gtdigilabel_(iConfig.getParameter<vector<edm::InputTag> >("gtdigilabel")),  
	lepcollections_(iConfig.getParameter<std::vector<std::string> >("LepCollections")),
	jetcollections_(iConfig.getParameter<std::vector<std::string> >("JetCollections")),
	jettypes_(iConfig.getParameter<std::vector<std::string> >("JetTypes")),
        doGenJets_(iConfig.getUntrackedParameter<bool>("DoGenJets")) 
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
		//hasElectrons(iEvent); 
		hasJets(iEvent, iSetup); 
		saveHLT(iEvent); 
		saveL1T(iEvent); 

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

	for(unsigned i=0; i<jetcollections_.size(); i++) {
		if(i >= MAX_JETCOLLECTIONS) break;

		JetInfo[i].RegisterTree(tree_,jetcollections_[i]);
		SubJetInfo[i].RegisterTree(tree_,jetcollections_[i].append("_subjets"));
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
	//   memset(&VertexInfo,0x00,sizeof(VertexInfo));
}

	bool
BprimeTobH::hasPrimaryVertex(const edm::Event& iEvent)
{
	memset(&VertexInfo,0x00,sizeof(VertexInfo));
	// Signal_Vz = 0.; 

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

		//for (reco::Vertex::trackRef_iterator it = iVertex->tracks_begin(); 
		// it != iVertex->tracks_end(); it++) {
		//  VertexInfo.Pt_Sum[VertexInfo.Size] += (*it)->pt();
		//  VertexInfo.Pt_Sum2[VertexInfo.Size] += ((*it)->pt() * (*it)->pt());
		//}

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
		//for (reco::Vertex::trackRef_iterator it = it_vtx->tracks_begin(); it != it_vtx->tracks_end(); it++) {
		//  VertexInfo.Pt_Sum[VertexInfo.Size] += (*it)->pt();
		//  VertexInfo.Pt_Sum2[VertexInfo.Size] += ((*it)->pt() * (*it)->pt());
		//}
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

	bool
BprimeTobH::hasJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	vector<edm::Handle<vector<pat::Jet> > > JetHandle;
	for(unsigned il=0; il<jetlabel_.size(); il++) {
		JetHandle.push_back(edm::Handle<vector<pat::Jet> >());

		// if ( jetlabel_[il].encode() == "goodPatJetsCA8PrunedPackedSubjet1")
		//   iEvent.getByLabel( "goodPatJetsCA8PrunedPacked", JetHandle[il]);
		// else
		iEvent.getByLabel( jetlabel_[il], JetHandle[il]);
	}

	for(unsigned icoll=0; icoll<jetcollections_.size(); icoll++) {
		//loop over collections
		if(icoll >= MAX_JETCOLLECTIONS) break;

		memset(&JetInfo[icoll],0x00,sizeof(JetInfo[icoll]));
		memset(&SubJetInfo[icoll],0x00,sizeof(SubJetInfo[icoll]));

		// For Jet Uncertainty

		edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
		iSetup.get<JetCorrectionsRecord>().get("AK5PF", JetCorParColl); //?? Hardcode ??
		JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
		JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

		if (JetHandle.size() <= icoll) continue;  

		for( vector<pat::Jet>::const_iterator it_jet = JetHandle[icoll]->begin(); it_jet != JetHandle[icoll]->end(); it_jet++ ) { 

			JetInfo[icoll].Index   [JetInfo[icoll].Size] = JetInfo[icoll].Size;
			JetInfo[icoll].NTracks [JetInfo[icoll].Size] = it_jet->associatedTracks().size();
			JetInfo[icoll].Et  [JetInfo[icoll].Size] = it_jet->et();
			JetInfo[icoll].Pt  [JetInfo[icoll].Size] = it_jet->pt();
			JetInfo[icoll].Eta [JetInfo[icoll].Size] = it_jet->eta();
			JetInfo[icoll].Phi [JetInfo[icoll].Size] = it_jet->phi();

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

			bool JetID(false);

			if (jettypes_[icoll] == "fatjet")  {
				JetID = true; //Apply JetID in PAT level
			}
			else if(it_jet->isPFJet()) {
				//Jet ID for PFJet
				edm::ParameterSet PS_pf;
				PS_pf.addParameter<std::string>("version", "FIRSTDATA");
				PS_pf.addParameter<std::string>("quality", "LOOSE");
				PFJetIDSelectionFunctor pfjetIDLOOSE(PS_pf) ;
				pat::strbitset ret = pfjetIDLOOSE.getBitTemplate() ;
				ret.set(false);
				JetID = pfjetIDLOOSE(*it_jet, ret);
			}
			else JetID = false; 

			JetInfo[icoll].JetIDLOOSE[JetInfo[icoll].Size] = (JetID) ?  1 : 0;

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

			JetInfo[icoll].Px  [JetInfo[icoll].Size] = it_jet->px(); 
			JetInfo[icoll].Py  [JetInfo[icoll].Size] = it_jet->py(); 
			JetInfo[icoll].Pz  [JetInfo[icoll].Size] = it_jet->pz(); 
			JetInfo[icoll].Energy  [JetInfo[icoll].Size] = it_jet->energy(); 

			JetInfo[icoll].Mass[JetInfo[icoll].Size] = it_jet->mass();
			JetInfo[icoll].Area[JetInfo[icoll].Size] = it_jet->jetArea();

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
			}

			// loop for subjets
			if (doSubjets_) { 

				for (unsigned ndau = 0; ndau < it_jet->numberOfDaughters(); ++ndau) {
					reco::CandidatePtr dptr = it_jet->daughterPtr(ndau);
					if(dptr.isAvailable()){
						pat::Jet const * subjet = dynamic_cast<pat::Jet const *>( &(*dptr));
						if (subjet) { 
							if (subjet->isPFJet()) { 

								//edm::LogInfo("Subjets") << " Coll = " << jetcollections_[icoll] << " Jet no. =  " << JetInfo[icoll].Size << " Exist for jetCollection " << jetcollections_[icoll] << " index = " << SubJetInfo[icoll].Size ; 

								SubJetInfo[icoll].Index   [SubJetInfo[icoll].Size] = SubJetInfo[icoll].Size;
								SubJetInfo[icoll].ParentIndex [SubJetInfo[icoll].Size] = JetInfo[icoll].Size;

								SubJetInfo[icoll].NTracks [SubJetInfo[icoll].Size] = subjet->associatedTracks().size();
								SubJetInfo[icoll].Et  [SubJetInfo[icoll].Size] = subjet->et();
								SubJetInfo[icoll].Pt  [SubJetInfo[icoll].Size] = subjet->pt();
								SubJetInfo[icoll].Eta [SubJetInfo[icoll].Size] = subjet->eta();
								SubJetInfo[icoll].Phi [SubJetInfo[icoll].Size] = subjet->phi();

								SubJetInfo[icoll].JetCharge[SubJetInfo[icoll].Size] = subjet->jetCharge();

								if (subjet->isPFJet())  {
									SubJetInfo[icoll].NCH[JetInfo[icoll].Size] = subjet->chargedMultiplicity();
									SubJetInfo[icoll].CEF[JetInfo[icoll].Size] = subjet->chargedEmEnergyFraction();
									SubJetInfo[icoll].NHF[JetInfo[icoll].Size] = subjet->neutralHadronEnergyFraction();
									SubJetInfo[icoll].NEF[JetInfo[icoll].Size] = subjet->neutralEmEnergyFraction();
									SubJetInfo[icoll].CHF[JetInfo[icoll].Size] = subjet->chargedHadronEnergyFraction();
								}

								SubJetInfo[icoll].JetBProbBJetTags[SubJetInfo[icoll].Size] = subjet->bDiscriminator("jetBProbabilityBJetTags");
								SubJetInfo[icoll].JetProbBJetTags [SubJetInfo[icoll].Size] = subjet->bDiscriminator("jetProbabilityBJetTags");
								SubJetInfo[icoll].TrackCountHiPurBJetTags [SubJetInfo[icoll].Size] = subjet->bDiscriminator("trackCountingHighPurBJetTags"); 
								SubJetInfo[icoll].CombinedSVBJetTags  [SubJetInfo[icoll].Size] = subjet->bDiscriminator("combinedSecondaryVertexBJetTags");
								SubJetInfo[icoll].CombinedSVMVABJetTags   [SubJetInfo[icoll].Size] = subjet->bDiscriminator("combinedSecondaryVertexMVABJetTags");
								SubJetInfo[icoll].SoftElecByIP3dBJetTags  [SubJetInfo[icoll].Size] = subjet->bDiscriminator("softElectronByIP3dBJetTags");
								SubJetInfo[icoll].SoftElecByPtBJetTags[SubJetInfo[icoll].Size] = subjet->bDiscriminator("softElectronByPtBJetTags");
								SubJetInfo[icoll].SoftMuonBJetTags[SubJetInfo[icoll].Size] = subjet->bDiscriminator("softMuonBJetTags");
								SubJetInfo[icoll].SoftMuonByIP3dBJetTags  [SubJetInfo[icoll].Size] = subjet->bDiscriminator("softMuonByIP3dBJetTags");
								SubJetInfo[icoll].SoftMuonByPtBJetTags[SubJetInfo[icoll].Size] = subjet->bDiscriminator("softMuonByPtBJetTags");
								SubJetInfo[icoll].DoubleSVHighEffBJetTags [SubJetInfo[icoll].Size] = subjet->bDiscriminator("doubleSecondaryVertexHighEffBJetTags"); 

								SubJetInfo[icoll].Px  [SubJetInfo[icoll].Size] = subjet->px(); 
								SubJetInfo[icoll].Py  [SubJetInfo[icoll].Size] = subjet->py(); 
								SubJetInfo[icoll].Pz  [SubJetInfo[icoll].Size] = subjet->pz(); 
								SubJetInfo[icoll].Energy  [SubJetInfo[icoll].Size] = subjet->energy(); 
								SubJetInfo[icoll].Mass[SubJetInfo[icoll].Size] = subjet->mass();
								SubJetInfo[icoll].Area[SubJetInfo[icoll].Size] = subjet->jetArea();

								if (!iEvent.isRealData()) {
									SubJetInfo[icoll].GenFlavor[SubJetInfo[icoll].Size] = subjet->partonFlavour(); 
									if (doGenJets_) {
										const reco::GenJet * genjetSubJet = subjet->genJet();
										SubJetInfo[icoll].GenJetPt [SubJetInfo[icoll].Size] = genjetSubJet->pt(); 
										SubJetInfo[icoll].GenJetEta[SubJetInfo[icoll].Size] = genjetSubJet->eta();
										SubJetInfo[icoll].GenJetPhi[SubJetInfo[icoll].Size] = genjetSubJet->eta(); 
										const reco::GenParticle* partonSubJet = subjet->genParton();
										const reco::Candidate*   genCandSubJet = partonSubJet;
										int qpTagSubJet(0) ; 
										while (genCandSubJet!=NULL && genCandSubJet->numberOfMothers()==1) {
											genCandSubJet = genCandSubJet->mother(0);
											if (abs(genCandSubJet->pdgId())==7 ) qpTagSubJet = 7;	// check if it's bprime
											else if (abs(genCandSubJet->pdgId())==8 ) qpTagSubJet = 8;	// check if it's tprime.
											else qpTagSubJet = 0 ; 
											if (abs(genCandSubJet->pdgId())==23) SubJetInfo[icoll].GenMCTag[JetInfo[icoll].Size] = 23;
											if (abs(genCandSubJet->pdgId())==24) SubJetInfo[icoll].GenMCTag[JetInfo[icoll].Size] = 24;
											if (abs(genCandSubJet->pdgId())==25) SubJetInfo[icoll].GenMCTag[JetInfo[icoll].Size] = 25;
										}
										SubJetInfo[icoll].GenMCTag[JetInfo[icoll].Size] += qpTagSubJet*100 ; 
									}
								}

								SubJetInfo[icoll].Size++;

							}
							else {
								edm::LogInfo("Subjets") << " Coll = " << jetcollections_[icoll] << " Jet no. =  " << JetInfo[icoll].Size << " Does not for jetCollection " << jetcollections_[icoll] ; 
							}
						} // subjet pointer exists
					}

				} // loop over subjets 
			} // if doSubJets 

			JetInfo[icoll].Size++;

		}//loop over jets in collection
	}//have jet collection

	return true; 
}


	void
BprimeTobH::saveHLT(const edm::Event& iEvent)
{
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

	void
BprimeTobH::saveL1T(const edm::Event& iEvent)
{
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


//define this as a plug-in
DEFINE_FWK_MODULE(BprimeTobH);
