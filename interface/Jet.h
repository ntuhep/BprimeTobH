#ifndef BPRIMETOBH_INTERFACE_JET_H
#define BPRIMETOBH_INTERFACE_JET_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "BpbH/BprimeTobH/interface/format.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <boost/algorithm/string.hpp>

class Jet {

  public:

    Jet (const Jet& jet) { 

      Et_ = jet.Et() ;
      Pt_ = jet.Pt() ;
      Unc_ = jet.Unc() ;
      Eta_ = jet.Eta() ;
      Phi_ = jet.Phi() ;

      Energy_ = jet.Energy() ; 
      Px_ = jet.Px() ;     
      Py_ = jet.Py() ;     
      Pz_ = jet.Pz() ;     
      Mass_ = jet.Mass() ;
      Area_ = jet.Area() ;

      EtPruned_ = jet.EtPruned() ;            
      PtPruned_ = jet.PtPruned() ;            
      UncPruned_ = jet.UncPruned() ;            
      EtaPruned_ = jet.EtaPruned() ;            
      PhiPruned_ = jet.PhiPruned() ;             

      EnergyPruned_ = jet.EnergyPruned() ; 
      PxPruned_ = jet.PxPruned() ; 
      PyPruned_ = jet.PyPruned() ; 
      PzPruned_ = jet.PzPruned() ; 
      MassPruned_ = jet.MassPruned() ; 
      AreaPruned_ = jet.AreaPruned() ; 

      tau1_ = jet.tau1() ; 
      tau2_ = jet.tau2() ; 
      tau3_ = jet.tau3() ; 

      GenJetPt_ = jet.GenJetPt() ;
      GenJetEta_ = jet.GenJetEta() ;
      GenJetPhi_ = jet.GenJetPhi() ;
      GenPt_ = jet.GenPt() ;
      GenEta_ = jet.GenEta() ;
      GenPhi_ = jet.GenPhi() ;
      GenPdgID_ = jet.GenPdgID() ;
      GenFlavor_ = jet.GenFlavor() ;
      GenMCTag_ = jet.GenMCTag() ; 

      JetIDLOOSE_ = jet.JetIDLOOSE() ; 
      JetIDTIGHT_ = jet.JetIDTIGHT() ; 
      JetCharge_ = jet.JetCharge() ;
      QGTagsMLP_ = jet.QGTagsMLP() ;
      QGTagsLikelihood_ = jet.QGTagsLikelihood() ;
      NConstituents_ = jet.NConstituents() ;
      NTracks_ = jet.NTracks() ;
      NCH_ = jet.NCH() ;
      CEF_ = jet.CEF() ;
      NHF_ = jet.NHF() ;
      NEF_ = jet.NEF() ;
      CHF_ = jet.CHF() ;

      PtCorrRaw_ = jet.PtCorrRaw() ;  
      PtCorrL2_ = jet.PtCorrL2() ;  
      PtCorrL3_ = jet.PtCorrL3() ;  
      PtCorrL7g_ = jet.PtCorrL7g() ;
      PtCorrL7uds_ = jet.PtCorrL7uds() ;
      PtCorrL7c_ = jet.PtCorrL7c() ;  
      PtCorrL7b_ = jet.PtCorrL7b() ;  

      JetBProbBJetTags_ = jet.JetBProbBJetTags() ;
      JetProbBJetTags_ = jet.JetProbBJetTags() ;
      TrackCountHiPurBJetTags_ = jet.TrackCountHiPurBJetTags() ;  
      CombinedSVBJetTags_ = jet.CombinedSVBJetTags() ;
      CombinedSVMVABJetTags_ = jet.CombinedSVMVABJetTags() ;
      SoftElecByIP3dBJetTags_ = jet.SoftElecByIP3dBJetTags() ;
      SoftElecByPtBJetTags_ = jet.SoftElecByPtBJetTags() ;  
      SoftMuonBJetTags_ = jet.SoftMuonBJetTags() ;      
      SoftMuonByIP3dBJetTags_ = jet.SoftMuonByIP3dBJetTags() ;
      SoftMuonByPtBJetTags_ = jet.SoftMuonByPtBJetTags() ;  
      DoubleSVHighEffBJetTags_ = jet.DoubleSVHighEffBJetTags() ; 

      Jet_FatJetIdx_ = jet.Jet_FatJetIdx() ;
      Jet_SubJet1Idx_ = jet.Jet_SubJet1Idx() ;
      Jet_SubJet2Idx_ = jet.Jet_SubJet2Idx() ;

      IsBtaggedCSVL_ = false ;  
      IsBtaggedCSVM_ = false ;  
      IsBtaggedCSVT_ = false ;  
    } 

    Jet(JetInfoBranches& JetInfo, int& JetIndex) { 

      Et_ = JetInfo.Et[JetIndex];
      Pt_ = JetInfo.Pt[JetIndex];
      Unc_ = JetInfo.Unc[JetIndex];
      Eta_ = JetInfo.Eta[JetIndex];
      Phi_ = JetInfo.Phi[JetIndex];

      Energy_ = JetInfo.Energy[JetIndex]; 
      Px_ = JetInfo.Px[JetIndex];     
      Py_ = JetInfo.Py[JetIndex];     
      Pz_ = JetInfo.Pz[JetIndex];     
      Mass_ = JetInfo.Mass[JetIndex];
      Area_ = JetInfo.Area[JetIndex];

      EtPruned_ = JetInfo.EtPruned[JetIndex];            
      PtPruned_ = JetInfo.PtPruned[JetIndex];            
      UncPruned_ = JetInfo.UncPruned[JetIndex];            
      EtaPruned_ = JetInfo.EtaPruned[JetIndex];            
      PhiPruned_ = JetInfo.PhiPruned[JetIndex];             

      EnergyPruned_ = JetInfo.EnergyPruned[JetIndex]; 
      PxPruned_ = JetInfo.PxPruned[JetIndex]; 
      PyPruned_ = JetInfo.PyPruned[JetIndex]; 
      PzPruned_ = JetInfo.PzPruned[JetIndex]; 
      MassPruned_ = JetInfo.MassPruned[JetIndex]; 
      AreaPruned_ = JetInfo.AreaPruned[JetIndex]; 

      tau1_ = JetInfo.tau1[JetIndex]; 
      tau2_ = JetInfo.tau2[JetIndex]; 
      tau3_ = JetInfo.tau3[JetIndex]; 

      GenJetPt_ = JetInfo.GenJetPt[JetIndex];
      GenJetEta_ = JetInfo.GenJetEta[JetIndex];
      GenJetPhi_ = JetInfo.GenJetPhi[JetIndex];
      GenPt_ = JetInfo.GenPt[JetIndex];
      GenEta_ = JetInfo.GenEta[JetIndex];
      GenPhi_ = JetInfo.GenPhi[JetIndex];
      GenPdgID_ = JetInfo.GenPdgID[JetIndex];
      GenFlavor_ = JetInfo.GenFlavor[JetIndex];
      GenMCTag_ = JetInfo.GenMCTag[JetIndex]; 

      JetIDLOOSE_ = JetInfo.JetIDLOOSE[JetIndex]; 
      JetIDTIGHT_ = JetInfo.JetIDTIGHT[JetIndex]; 
      JetCharge_ = JetInfo.JetCharge[JetIndex];
      QGTagsMLP_ = JetInfo.QGTagsMLP[JetIndex];
      QGTagsLikelihood_ = JetInfo.QGTagsLikelihood[JetIndex];
      NConstituents_ = JetInfo.NConstituents[JetIndex];
      NTracks_ = JetInfo.NTracks[JetIndex];
      NCH_ = JetInfo.NCH[JetIndex];
      CEF_ = JetInfo.CEF[JetIndex];
      NHF_ = JetInfo.NHF[JetIndex];
      NEF_ = JetInfo.NEF[JetIndex];
      CHF_ = JetInfo.CHF[JetIndex];

      PtCorrRaw_ = JetInfo.PtCorrRaw[JetIndex];  
      PtCorrL2_ = JetInfo.PtCorrL2[JetIndex];  
      PtCorrL3_ = JetInfo.PtCorrL3[JetIndex];  
      PtCorrL7g_ = JetInfo.PtCorrL7g[JetIndex];
      PtCorrL7uds_ = JetInfo.PtCorrL7uds[JetIndex];
      PtCorrL7c_ = JetInfo.PtCorrL7c[JetIndex];  
      PtCorrL7b_ = JetInfo.PtCorrL7b[JetIndex];  

      JetBProbBJetTags_ = JetInfo.JetBProbBJetTags[JetIndex];
      JetProbBJetTags_ = JetInfo.JetProbBJetTags[JetIndex];
      TrackCountHiPurBJetTags_ = JetInfo.TrackCountHiPurBJetTags[JetIndex];  
      CombinedSVBJetTags_ = JetInfo.CombinedSVBJetTags[JetIndex];
      CombinedSVMVABJetTags_ = JetInfo.CombinedSVMVABJetTags[JetIndex];
      SoftElecByIP3dBJetTags_ = JetInfo.SoftElecByIP3dBJetTags[JetIndex];
      SoftElecByPtBJetTags_ = JetInfo.SoftElecByPtBJetTags[JetIndex];  
      SoftMuonBJetTags_ = JetInfo.SoftMuonBJetTags[JetIndex];      
      SoftMuonByIP3dBJetTags_ = JetInfo.SoftMuonByIP3dBJetTags[JetIndex];
      SoftMuonByPtBJetTags_ = JetInfo.SoftMuonByPtBJetTags[JetIndex];  
      DoubleSVHighEffBJetTags_ = JetInfo.DoubleSVHighEffBJetTags[JetIndex]; 

      Jet_FatJetIdx_ = JetInfo.Jet_FatJetIdx[JetIndex];
      Jet_SubJet1Idx_ = JetInfo.Jet_SubJet1Idx[JetIndex];
      Jet_SubJet2Idx_ = JetInfo.Jet_SubJet2Idx[JetIndex];

      IsBtaggedCSVL_ = false ;  
      IsBtaggedCSVM_ = false ;  
      IsBtaggedCSVT_ = false ;  
    }

    ~Jet () {} ;

    void Set_Et( float Et ) { Et_ = Et ; } 
    void Set_Pt( float Pt ) { Pt_ = Pt ; } 
    void Set_Unc( float Unc ) { Unc_ = Unc ; } 
    void Set_Eta( float Eta ) { Eta_ = Eta ; } 
    void Set_Phi( float Phi ) { Phi_ = Phi ; } 

    void Set_Energy( float Energy ) { Energy_ = Energy ; }  
    void Set_Px( float Px ) { Px_ = Px ; }      
    void Set_Py( float Py ) { Py_ = Py ; }      
    void Set_Pz( float Pz ) { Pz_ = Pz ; }      
    void Set_Mass( float Mass ) { Mass_ = Mass ; } 
    void Set_Area( float Area ) { Area_ = Area ; } 

    void Set_EtPruned( float EtPruned ) { EtPruned_ = EtPruned ; }             
    void Set_PtPruned( float PtPruned ) { PtPruned_ = PtPruned ; }             
    void Set_UncPruned( float UncPruned ) { UncPruned_ = UncPruned ; }             
    void Set_EtaPruned( float EtaPruned ) { EtaPruned_ = EtaPruned ; }             
    void Set_PhiPruned( float PhiPruned ) { PhiPruned_ = PhiPruned ; }              

    void Set_EnergyPruned( float EnergyPruned ) { EnergyPruned_ = EnergyPruned ; }  
    void Set_PxPruned( float PxPruned ) { PxPruned_ = PxPruned ; }  
    void Set_PyPruned( float PyPruned ) { PyPruned_ = PyPruned ; }  
    void Set_PzPruned( float PzPruned ) { PzPruned_ = PzPruned ; }  
    void Set_MassPruned( float MassPruned ) { MassPruned_ = MassPruned ; }  
    void Set_AreaPruned( float AreaPruned ) { AreaPruned_ = AreaPruned ; }  

    void Set_tau1( float tau1 ) { tau1_ = tau1 ; }  
    void Set_tau2( float tau2 ) { tau2_ = tau2 ; }  
    void Set_tau3( float tau3 ) { tau3_ = tau3 ; }  

    void Set_GenJetPt( float GenJetPt ) { GenJetPt_ = GenJetPt ; } 
    void Set_GenJetEta( float GenJetEta ) { GenJetEta_ = GenJetEta ; } 
    void Set_GenJetPhi( float GenJetPhi ) { GenJetPhi_ = GenJetPhi ; } 
    void Set_GenPt( float GenPt ) { GenPt_ = GenPt ; } 
    void Set_GenEta( float GenEta ) { GenEta_ = GenEta ; } 
    void Set_GenPhi( float GenPhi ) { GenPhi_ = GenPhi ; } 
    void Set_GenPdgID( int     GenPdgID ) {   GenPdgID_ =   GenPdgID ; } 
    void Set_GenFlavor( int     GenFlavor ) {   GenFlavor_ =   GenFlavor ; } 
    void Set_GenMCTag( int     GenMCTag ) {   GenMCTag_ =   GenMCTag ; }  

    void Set_JetIDLOOSE( bool   JetIDLOOSE ) {  JetIDLOOSE_ =  JetIDLOOSE ; }  
    void Set_JetIDTIGHT( bool   JetIDTIGHT ) {  JetIDTIGHT_ =  JetIDTIGHT ; }  
    void Set_JetCharge( float JetCharge ) { JetCharge_ = JetCharge ; } 
    void Set_QGTagsMLP( float QGTagsMLP ) { QGTagsMLP_ = QGTagsMLP ; } 
    void Set_QGTagsLikelihood( float QGTagsLikelihood ) { QGTagsLikelihood_ = QGTagsLikelihood ; } 
    void Set_NConstituents( int     NConstituents ) {   NConstituents_ =   NConstituents ; } 
    void Set_NTracks( int     NTracks ) {   NTracks_ =   NTracks ; } 
    void Set_NCH( int     NCH ) {   NCH_ =   NCH ; } 
    void Set_CEF( float CEF ) { CEF_ = CEF ; } 
    void Set_NHF( float NHF ) { NHF_ = NHF ; } 
    void Set_NEF( float NEF ) { NEF_ = NEF ; } 
    void Set_CHF( float CHF ) { CHF_ = CHF ; } 

    void Set_PtCorrRaw( float PtCorrRaw ) { PtCorrRaw_ = PtCorrRaw ; }   
    void Set_PtCorrL2( float PtCorrL2 ) { PtCorrL2_ = PtCorrL2 ; }   
    void Set_PtCorrL3( float PtCorrL3 ) { PtCorrL3_ = PtCorrL3 ; }   
    void Set_PtCorrL7g( float PtCorrL7g ) { PtCorrL7g_ = PtCorrL7g ; } 
    void Set_PtCorrL7uds( float PtCorrL7uds ) { PtCorrL7uds_ = PtCorrL7uds ; } 
    void Set_PtCorrL7c( float PtCorrL7c ) { PtCorrL7c_ = PtCorrL7c ; }   
    void Set_PtCorrL7b( float PtCorrL7b ) { PtCorrL7b_ = PtCorrL7b ; }   

    void Set_JetBProbBJetTags( float JetBProbBJetTags ) { JetBProbBJetTags_ = JetBProbBJetTags ; } 
    void Set_JetProbBJetTags( float JetProbBJetTags ) { JetProbBJetTags_ = JetProbBJetTags ; } 
    void Set_TrackCountHiPurBJetTags( float TrackCountHiPurBJetTags ) { TrackCountHiPurBJetTags_ = TrackCountHiPurBJetTags ; }   
    void Set_CombinedSVBJetTags( float CombinedSVBJetTags ) { CombinedSVBJetTags_ = CombinedSVBJetTags ; } 
    void Set_CombinedSVMVABJetTags( float CombinedSVMVABJetTags ) { CombinedSVMVABJetTags_ = CombinedSVMVABJetTags ; } 
    void Set_SoftElecByIP3dBJetTags( float SoftElecByIP3dBJetTags ) { SoftElecByIP3dBJetTags_ = SoftElecByIP3dBJetTags ; } 
    void Set_SoftElecByPtBJetTags( float SoftElecByPtBJetTags ) { SoftElecByPtBJetTags_ = SoftElecByPtBJetTags ; }   
    void Set_SoftMuonBJetTags( float SoftMuonBJetTags ) { SoftMuonBJetTags_ = SoftMuonBJetTags ; }       
    void Set_SoftMuonByIP3dBJetTags( float SoftMuonByIP3dBJetTags ) { SoftMuonByIP3dBJetTags_ = SoftMuonByIP3dBJetTags ; } 
    void Set_SoftMuonByPtBJetTags( float SoftMuonByPtBJetTags ) { SoftMuonByPtBJetTags_ = SoftMuonByPtBJetTags ; }   
    void Set_DoubleSVHighEffBJetTags( float DoubleSVHighEffBJetTags ) { DoubleSVHighEffBJetTags_ = DoubleSVHighEffBJetTags ; }  

    void Set_Jet_FatJetIdx( int Jet_FatJetIdx ) { Jet_FatJetIdx_ = Jet_FatJetIdx ; } 
    void Set_Jet_SubJet1Idx( int Jet_SubJet1Idx ) { Jet_SubJet1Idx_ = Jet_SubJet1Idx ; } 
    void Set_Jet_SubJet2Idx( int Jet_SubJet2Idx ) { Jet_SubJet2Idx_ = Jet_SubJet2Idx ; } 

    void Set_IsBtaggedCSVL( bool IsBtaggedCSVL ) { IsBtaggedCSVL_ = IsBtaggedCSVL ; }  
    void Set_IsBtaggedCSVM( bool IsBtaggedCSVM ) { IsBtaggedCSVM_ = IsBtaggedCSVM ; }  
    void Set_IsBtaggedCSVT( bool IsBtaggedCSVT ) { IsBtaggedCSVT_ = IsBtaggedCSVT ; }  

    void setIsBTagged (std::string algo, bool pass) { 
      if ( boost::iequals(algo, "CSVL") ) pass ? IsBtaggedCSVL_ = true : IsBtaggedCSVL_ = false ; 
      else if ( boost::iequals(algo, "CSVM") ) { pass ? IsBtaggedCSVL_ = true , IsBtaggedCSVM_ = true : IsBtaggedCSVM_ = false ; }
      else if ( boost::iequals(algo, "CSVT") ) { pass ? IsBtaggedCSVL_ = true , IsBtaggedCSVM_ = true , IsBtaggedCSVT_ = true : IsBtaggedCSVT_ = false ; } 
      else edm::LogError("JetBtag") << ">>>>> Wrong algo name given! " ; 
    }

    float Et() const { return Et_ ; } 
    float Pt() const { return Pt_ ; } 
    float Unc() const { return Unc_ ; } 
    float Eta() const { return Eta_ ; } 
    float Phi() const { return Phi_ ; } 

    float Energy() const { return Energy_ ; }  
    float Px() const { return Px_ ; }      
    float Py() const { return Py_ ; }      
    float Pz() const { return Pz_ ; }      
    float Mass() const { return Mass_ ; } 
    float Area() const { return Area_ ; } 

    float EtPruned() const { return EtPruned_ ; }             
    float PtPruned() const { return PtPruned_ ; }             
    float UncPruned() const { return UncPruned_ ; }             
    float EtaPruned() const { return EtaPruned_ ; }             
    float PhiPruned() const { return PhiPruned_ ; }              

    float EnergyPruned() const { return EnergyPruned_ ; }  
    float PxPruned() const { return PxPruned_ ; }  
    float PyPruned() const { return PyPruned_ ; }  
    float PzPruned() const { return PzPruned_ ; }  
    float MassPruned() const { return MassPruned_ ; }  
    float AreaPruned() const { return AreaPruned_ ; }  

    float tau1() const { return tau1_ ; }  
    float tau2() const { return tau2_ ; }  
    float tau3() const { return tau3_ ; }  

    float GenJetPt() const { return GenJetPt_ ; } 
    float GenJetEta() const { return GenJetEta_ ; } 
    float GenJetPhi() const { return GenJetPhi_ ; } 
    float GenPt() const { return GenPt_ ; } 
    float GenEta() const { return GenEta_ ; } 
    float GenPhi() const { return GenPhi_ ; } 
    int   GenPdgID() const { return   GenPdgID_ ; } 
    int   GenFlavor() const { return   GenFlavor_ ; } 
    int   GenMCTag() const { return   GenMCTag_ ; }  

    bool  JetIDLOOSE() const { return  JetIDLOOSE_ ; }  
    bool  JetIDTIGHT() const { return  JetIDTIGHT_ ; }  
    float JetCharge() const { return JetCharge_ ; } 
    float QGTagsMLP() const { return QGTagsMLP_ ; } 
    float QGTagsLikelihood() const { return QGTagsLikelihood_ ; } 
    int   NConstituents() const { return   NConstituents_ ; } 
    int   NTracks() const { return   NTracks_ ; } 
    int   NCH() const { return   NCH_ ; } 
    float CEF() const { return CEF_ ; } 
    float NHF() const { return NHF_ ; } 
    float NEF() const { return NEF_ ; } 
    float CHF() const { return CHF_ ; } 

    float PtCorrRaw() const { return PtCorrRaw_ ; }   
    float PtCorrL2() const { return PtCorrL2_ ; }   
    float PtCorrL3() const { return PtCorrL3_ ; }   
    float PtCorrL7g() const { return PtCorrL7g_ ; } 
    float PtCorrL7uds() const { return PtCorrL7uds_ ; } 
    float PtCorrL7c() const { return PtCorrL7c_ ; }   
    float PtCorrL7b() const { return PtCorrL7b_ ; }   

    float JetBProbBJetTags() const { return JetBProbBJetTags_ ; } 
    float JetProbBJetTags() const { return JetProbBJetTags_ ; } 
    float TrackCountHiPurBJetTags() const { return TrackCountHiPurBJetTags_ ; }   
    float CombinedSVBJetTags() const { return CombinedSVBJetTags_ ; } 
    float CombinedSVMVABJetTags() const { return CombinedSVMVABJetTags_ ; } 
    float SoftElecByIP3dBJetTags() const { return SoftElecByIP3dBJetTags_ ; } 
    float SoftElecByPtBJetTags() const { return SoftElecByPtBJetTags_ ; }   
    float SoftMuonBJetTags() const { return SoftMuonBJetTags_ ; }       
    float SoftMuonByIP3dBJetTags() const { return SoftMuonByIP3dBJetTags_ ; } 
    float SoftMuonByPtBJetTags() const { return SoftMuonByPtBJetTags_ ; }   
    float DoubleSVHighEffBJetTags() const { return DoubleSVHighEffBJetTags_ ; }  

    int Jet_FatJetIdx() const { return Jet_FatJetIdx_ ; } 
    int Jet_SubJet1Idx() const { return Jet_SubJet1Idx_ ; } 
    int Jet_SubJet2Idx() const { return Jet_SubJet2Idx_ ; } 

    bool IsBtaggedCSVL() const { return IsBtaggedCSVL_ ; }  
    bool IsBtaggedCSVM() const { return IsBtaggedCSVM_ ; }  
    bool IsBtaggedCSVT() const { return IsBtaggedCSVT_ ; }  

  private:

    float Et_ ;
    float Pt_ ;
    float Unc_ ;
    float Eta_ ;
    float Phi_ ;

    float Energy_ ; 
    float Px_ ;     
    float Py_ ;     
    float Pz_ ;     
    float Mass_ ;
    float Area_ ;

    float EtPruned_ ;            
    float PtPruned_ ;            
    float UncPruned_ ;            
    float EtaPruned_ ;            
    float PhiPruned_ ;             

    float EnergyPruned_ ; 
    float PxPruned_ ; 
    float PyPruned_ ; 
    float PzPruned_ ; 
    float MassPruned_ ; 
    float AreaPruned_ ; 

    float tau1_ ; 
    float tau2_ ; 
    float tau3_ ; 

    float GenJetPt_ ;
    float GenJetEta_ ;
    float GenJetPhi_ ;
    float GenPt_ ;
    float GenEta_ ;
    float GenPhi_ ;
    int   GenPdgID_ ;
    int   GenFlavor_ ;
    int   GenMCTag_ ; 

    bool  JetIDLOOSE_ ; 
    bool  JetIDTIGHT_ ; 
    float JetCharge_ ;
    float QGTagsMLP_ ;
    float QGTagsLikelihood_ ;
    int   NConstituents_ ;
    int   NTracks_ ;
    int   NCH_ ;
    float CEF_ ;
    float NHF_ ;
    float NEF_ ;
    float CHF_ ;

    float PtCorrRaw_ ;  
    float PtCorrL2_ ;  
    float PtCorrL3_ ;  
    float PtCorrL7g_ ;
    float PtCorrL7uds_ ;
    float PtCorrL7c_ ;  
    float PtCorrL7b_ ;  

    float JetBProbBJetTags_ ;
    float JetProbBJetTags_ ;
    float TrackCountHiPurBJetTags_ ;  
    float CombinedSVBJetTags_ ;
    float CombinedSVMVABJetTags_ ;
    float SoftElecByIP3dBJetTags_ ;
    float SoftElecByPtBJetTags_ ;  
    float SoftMuonBJetTags_ ;      
    float SoftMuonByIP3dBJetTags_ ;
    float SoftMuonByPtBJetTags_ ;  
    float DoubleSVHighEffBJetTags_ ; 

    int Jet_FatJetIdx_ ;
    int Jet_SubJet1Idx_ ;
    int Jet_SubJet2Idx_ ;

    bool IsBtaggedCSVL_ ;  
    bool IsBtaggedCSVM_ ;  
    bool IsBtaggedCSVT_ ;  

};

#endif 
