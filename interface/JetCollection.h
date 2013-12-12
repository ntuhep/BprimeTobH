#ifndef BPRIMETOBH_INTERFACE_JETCOLLECTION_H
#define BPRIMETOBH_INTERFACE_JETCOLLECTION_H

#include <vector>
#include "BpbH/BprimeTobH/interface/Jet.h"

typedef std::vector<Jet> JetCollection ; 

class MakeJetInfoBranches {

  public: 

    MakeJetInfoBranches (JetCollection& jetcoll, TTree* tree, std::string jetinfoname) { 

    for (JetInfo_.Size = 0; JetInfo_.Size < (int)jetcoll.size(); ++JetInfo_.Size) { 

      JetInfo_.Et[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Et() ;
      JetInfo_.Pt[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Pt() ;
      JetInfo_.Unc[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Unc() ;
      JetInfo_.Eta[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Eta() ;
      JetInfo_.Phi[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Phi() ;

      JetInfo_.Energy[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Energy() ; 
      JetInfo_.Px[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Px() ;     
      JetInfo_.Py[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Py() ;     
      JetInfo_.Pz[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Pz() ;     
      JetInfo_.Mass[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Mass() ;
      JetInfo_.Area[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Area() ;

      JetInfo_.EtPruned[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).EtPruned() ;            
      JetInfo_.PtPruned[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PtPruned() ;            
      JetInfo_.UncPruned[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).UncPruned() ;            
      JetInfo_.EtaPruned[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).EtaPruned() ;            
      JetInfo_.PhiPruned[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PhiPruned() ;             

      JetInfo_.EnergyPruned[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).EnergyPruned() ; 
      JetInfo_.PxPruned[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PxPruned() ; 
      JetInfo_.PyPruned[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PyPruned() ; 
      JetInfo_.PzPruned[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PzPruned() ; 
      JetInfo_.MassPruned[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).MassPruned() ; 
      JetInfo_.AreaPruned[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).AreaPruned() ; 

      JetInfo_.tau1[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).tau1() ; 
      JetInfo_.tau2[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).tau2() ; 
      JetInfo_.tau3[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).tau3() ; 

      JetInfo_.GenJetPt[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).GenJetPt() ;
      JetInfo_.GenJetEta[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).GenJetEta() ;
      JetInfo_.GenJetPhi[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).GenJetPhi() ;
      JetInfo_.GenPt[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).GenPt() ;
      JetInfo_.GenEta[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).GenEta() ;
      JetInfo_.GenPhi[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).GenPhi() ;
      JetInfo_.GenPdgID[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).GenPdgID() ;
      JetInfo_.GenFlavor[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).GenFlavor() ;
      JetInfo_.GenMCTag[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).GenMCTag() ; 

      JetInfo_.JetIDLOOSE[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).JetIDLOOSE() ; 
      JetInfo_.JetIDTIGHT[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).JetIDTIGHT() ; 
      JetInfo_.JetCharge[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).JetCharge() ;
      JetInfo_.QGTagsMLP[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).QGTagsMLP() ;
      JetInfo_.QGTagsLikelihood[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).QGTagsLikelihood() ;
      JetInfo_.NConstituents[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).NConstituents() ;
      JetInfo_.NTracks[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).NTracks() ;
      JetInfo_.NCH[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).NCH() ;
      JetInfo_.CEF[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).CEF() ;
      JetInfo_.NHF[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).NHF() ;
      JetInfo_.NEF[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).NEF() ;
      JetInfo_.CHF[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).CHF() ;

      JetInfo_.PtCorrRaw[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PtCorrRaw() ;  
      JetInfo_.PtCorrL2[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PtCorrL2() ;  
      JetInfo_.PtCorrL3[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PtCorrL3() ;  
      JetInfo_.PtCorrL7g[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PtCorrL7g() ;
      JetInfo_.PtCorrL7uds[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PtCorrL7uds() ;
      JetInfo_.PtCorrL7c[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PtCorrL7c() ;  
      JetInfo_.PtCorrL7b[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).PtCorrL7b() ;  

      JetInfo_.JetBProbBJetTags[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).JetBProbBJetTags() ;
      JetInfo_.JetProbBJetTags[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).JetProbBJetTags() ;
      JetInfo_.TrackCountHiPurBJetTags[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).TrackCountHiPurBJetTags() ;  
      JetInfo_.CombinedSVBJetTags[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).CombinedSVBJetTags() ;
      JetInfo_.CombinedSVMVABJetTags[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).CombinedSVMVABJetTags() ;
      JetInfo_.SoftElecByIP3dBJetTags[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).SoftElecByIP3dBJetTags() ;
      JetInfo_.SoftElecByPtBJetTags[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).SoftElecByPtBJetTags() ;  
      JetInfo_.SoftMuonBJetTags[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).SoftMuonBJetTags() ;      
      JetInfo_.SoftMuonByIP3dBJetTags[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).SoftMuonByIP3dBJetTags() ;
      JetInfo_.SoftMuonByPtBJetTags[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).SoftMuonByPtBJetTags() ;  
      JetInfo_.DoubleSVHighEffBJetTags[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).DoubleSVHighEffBJetTags() ; 

      JetInfo_.Jet_FatJetIdx[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Jet_FatJetIdx() ;
      JetInfo_.Jet_SubJet1Idx[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Jet_SubJet1Idx() ;
      JetInfo_.Jet_SubJet2Idx[JetInfo_.Size] = (jetcoll.at(JetInfo_.Size)).Jet_SubJet2Idx() ;

    }

    JetInfo_.RegisterTree(tree, jetinfoname) ; 

  }

    ~MakeJetInfoBranches () {} 
  private:
    JetInfoBranches JetInfo_ ; 

};

#endif 

