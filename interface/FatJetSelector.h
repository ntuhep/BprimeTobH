#ifndef BPRIMETOBH_INTERFACE_FATJETSELECTOR_H
#define BPRIMETOBH_INTERFACE_FATJETSELECTOR_H

#include "BpbH/BprimeTobH/interface/format.h"
#include "BpbH/BprimeTobH/interface/JetSelector.h"

#include <TLorentzVector.h>

class FatJetSelector : public Selector<int> { 
  public:
    enum JETTYPES_t { CA8JET, PRUNEDCA8JET, BTAGGEDCA8JET, HIGGSJET, N_JETTYPES} ; 

    FatJetSelector () {}

    FatJetSelector (edm::ParameterSet const& params) {

      std::string jettypeStr = params.getParameter<std::string>("jettype") ;
      if (jettypeStr == "CA8JET") jettype_ = CA8JET; 
      else if (jettypeStr == "PRUNEDCA8JET") jettype_ = PRUNEDCA8JET;
      else if (jettypeStr == "BTAGGEDCA8JET") jettype_ = BTAGGEDCA8JET;
      else if (jettypeStr == "HIGGSJET") jettype_ = HIGGSJET;  
      else edm::LogError("WrongFatJetType") << " Check jet type!!!!" ; 

      jetSelParams_ = params.getParameter<edm::ParameterSet>("JetSelParams") ; 

      push_back("fatJetPtMin") ;
      push_back("fatJetPtMax") ;
      push_back("fatJetAbsEtaMax") ;
      push_back("fatJetMassMin") ;
      push_back("fatJetMassMax") ;
      push_back("fatJetPrunedMassMin") ;
      push_back("fatJetPrunedMassMax") ;
      push_back("fatJetTau2ByTau1Max") ;
      push_back("dRSubjetsMin") ; 
      push_back("dRSubjetsMax") ; 
      push_back("fatJetCSVDiscMin") ;
      push_back("fatJetCSVDiscMax") ;
      push_back("subjet1CSVDiscMin") ;
      push_back("subjet1CSVDiscMax") ;
      push_back("subjet2CSVDiscMin") ;
      push_back("subjet2CSVDiscMax") ;

      set("fatJetPtMin"          ,params.getParameter<double>("fatJetPtMin") ) ;
      set("fatJetPtMax"          ,params.getParameter<double>("fatJetPtMax") ) ;
      set("fatJetAbsEtaMax"      ,params.getParameter<double>("fatJetAbsEtaMax") ) ;
      set("fatJetMassMin"        ,params.getParameter<double>("fatJetMassMin") ) ;
      set("fatJetMassMax"        ,params.getParameter<double>("fatJetMassMax") ) ;
      set("fatJetPrunedMassMin"  ,params.getParameter<double>("fatJetPrunedMassMin") ) ;
      set("fatJetPrunedMassMax"  ,params.getParameter<double>("fatJetPrunedMassMax") ) ;
      set("fatJetTau2ByTau1Max"  ,params.getParameter<double>("fatJetTau2ByTau1Max") ) ;
      set("dRSubjetsMin"         ,params.getParameter<double>("dRSubjetsMin") ) ;
      set("dRSubjetsMax"         ,params.getParameter<double>("dRSubjetsMax") ) ;
      if (jettype_ == BTAGGEDCA8JET) {
        set("fatJetCSVDiscMin"     ,params.getParameter<double>("fatJetCSVDiscMin") ,true) ;
        set("fatJetCSVDiscMax"     ,params.getParameter<double>("fatJetCSVDiscMax") ,true) ;
      }
      else {
        set("fatJetCSVDiscMin"     ,params.getParameter<double>("fatJetCSVDiscMin") ,false) ;
        set("fatJetCSVDiscMax"     ,params.getParameter<double>("fatJetCSVDiscMax") ,false) ;
      }
      if (jettype_ == HIGGSJET) {
        set("subjet1CSVDiscMin"    ,params.getParameter<double>("subjet1CSVDiscMin") ,true) ;
        set("subjet1CSVDiscMax"    ,params.getParameter<double>("subjet1CSVDiscMax") ,true) ;
        set("subjet2CSVDiscMin"    ,params.getParameter<double>("subjet2CSVDiscMin") ,true) ;
        set("subjet2CSVDiscMax"    ,params.getParameter<double>("subjet2CSVDiscMax") ,true) ;
      }
      else {
        set("subjet1CSVDiscMin"    ,params.getParameter<double>("subjet1CSVDiscMin") ,false) ;
        set("subjet1CSVDiscMax"    ,params.getParameter<double>("subjet1CSVDiscMax") ,false) ;
        set("subjet2CSVDiscMin"    ,params.getParameter<double>("subjet2CSVDiscMin") ,false) ;
        set("subjet2CSVDiscMax"    ,params.getParameter<double>("subjet2CSVDiscMax") ,false) ;
      }

      indexfatJetPtMin_         = index_type(&bits_ ,"fatJetPtMin") ;
      indexfatJetPtMax_         = index_type(&bits_ ,"fatJetPtMax") ;
      indexfatJetAbsEtaMax_     = index_type(&bits_ ,"fatJetAbsEtaMax") ;
      indexfatJetMassMin_       = index_type(&bits_ ,"fatJetMassMin") ;
      indexfatJetMassMax_       = index_type(&bits_ ,"fatJetMassMax") ;
      indexfatJetPrunedMassMin_ = index_type(&bits_ ,"fatJetPrunedMassMin") ;
      indexfatJetPrunedMassMax_ = index_type(&bits_ ,"fatJetPrunedMassMax") ;
      indexfatJetTau2ByTau1Max_ = index_type(&bits_ ,"fatJetTau2ByTau1Max") ;
      indexdRSubjetsMin_        = index_type(&bits_ ,"dRSubjetsMin") ;
      indexdRSubjetsMax_        = index_type(&bits_ ,"dRSubjetsMax") ;
      indexfatJetCSVDiscMin_    = index_type(&bits_ ,"fatJetCSVDiscMin") ;
      indexfatJetCSVDiscMax_    = index_type(&bits_ ,"fatJetCSVDiscMax") ;
      indexsubjet1CSVDiscMin_   = index_type(&bits_ ,"subjet1CSVDiscMin") ;
      indexsubjet1CSVDiscMax_   = index_type(&bits_ ,"subjet1CSVDiscMax") ;
      indexsubjet2CSVDiscMin_   = index_type(&bits_ ,"subjet2CSVDiscMin") ;
      indexsubjet2CSVDiscMax_   = index_type(&bits_ ,"subjet2CSVDiscMax") ;

      retInternal_ = getBitTemplate();   

    }

    bool operator()( int & jet, pat::strbitset & ret ) {
      return true; 
    }

    bool operator()( int const & jet, pat::strbitset & ret ) {
      return true; 
    }

    using Selector<int>::operator();

    bool operator()(JetInfoBranches& jetInfo,  int const & jet, JetInfoBranches& SubJetInfo, pat::strbitset & ret ) {
      return firstDataCuts (jetInfo, jet, SubJetInfo, ret ) ; 
    }

    bool firstDataCuts (JetInfoBranches& jetInfo, const int & jet, JetInfoBranches& SubJetInfo, pat::strbitset & ret ) { 

      bool isJetSel(false) ;
      ret.set(false);
      JetSelector jetSel(jetSelParams_) ;
      pat::strbitset retjetid = jetSel.getBitTemplate() ;
      retjetid.set(false) ;
      if (jetSel(jetInfo, jet, retjetid) != 0) isJetSel = true ;

      if (isJetSel == false) return false ; 

      double jetPt         = jetInfo.Pt[jet];
      double jetAbsEta     = std::abs(jetInfo.Eta[jet]); 
      double jetMass       = jetInfo.Mass[jet];
      double jetMassPruned = jetInfo.MassPruned[jet];
      double jetTau1       = jetInfo.tau1[jet];
      double jetTau2       = jetInfo.tau2[jet];
      double jetCSVDisc    = jetInfo.CombinedSVBJetTags[jet]; 
      int    iSubJet1      = jetInfo.Jet_SubJet1Idx[jet]; 
      int    iSubJet2      = jetInfo.Jet_SubJet2Idx[jet]; 
      TLorentzVector subjet1_p4, subjet2_p4;
      subjet1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subjet2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]);
      double subjet_dy     = subjet1_p4.Rapidity() - subjet2_p4.Rapidity() ;
      double subjet_dphi   = subjet1_p4.DeltaPhi(subjet2_p4); ;
      double subjet_dyphi  = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;
      double subjet1CSVDisc = SubJetInfo.CombinedSVBJetTags[iSubJet1] ; 
      double subjet2CSVDisc = SubJetInfo.CombinedSVBJetTags[iSubJet2] ; 

      if( SubJetInfo.Pt[iSubJet1]==0. || SubJetInfo.Pt[iSubJet2]==0. ) return false; //// skip fat jets for which one of the subjets has pT=0
      if ( subjet_dyphi < jetMass/jetPt ) return false ; 

      if ( ignoreCut(indexfatJetPtMin_)         || jetPt > cut(indexfatJetPtMin_, double() ) ) passCut( ret ,indexfatJetPtMin_) ; 
      if ( ignoreCut(indexfatJetPtMax_)         || jetPt < cut(indexfatJetPtMax_, double() ) ) passCut( ret ,indexfatJetPtMax_) ; 
      if ( ignoreCut(indexfatJetAbsEtaMax_)     || jetAbsEta < cut(indexfatJetAbsEtaMax_, double() ) ) passCut( ret ,indexfatJetAbsEtaMax_) ; 
      if ( ignoreCut(indexfatJetMassMin_)       || jetMass > cut(indexfatJetMassMin_, double() ) ) passCut( ret ,indexfatJetMassMin_) ;
      if ( ignoreCut(indexfatJetMassMax_)       || jetMass < cut(indexfatJetMassMax_, double() ) ) passCut( ret ,indexfatJetMassMax_) ;
      if ( ignoreCut(indexfatJetPrunedMassMin_) || jetMassPruned > cut(indexfatJetPrunedMassMin_, double() ) ) passCut( ret ,indexfatJetPrunedMassMin_) ;
      if ( ignoreCut(indexfatJetPrunedMassMax_) || jetMassPruned < cut(indexfatJetPrunedMassMax_, double() ) ) passCut( ret ,indexfatJetPrunedMassMax_) ;
      if ( ignoreCut(indexfatJetTau2ByTau1Max_) || jetTau2/jetTau1 < cut(indexfatJetTau2ByTau1Max_, double() ) ) passCut( ret ,indexfatJetTau2ByTau1Max_) ;
      if ( ignoreCut(indexdRSubjetsMin_)        || subjet_dyphi > cut(indexdRSubjetsMin_, double() ) ) passCut( ret, indexdRSubjetsMin_) ; 
      if ( ignoreCut(indexdRSubjetsMax_)        || subjet_dyphi < cut(indexdRSubjetsMax_, double() ) ) passCut( ret, indexdRSubjetsMax_) ; 

      if ( ignoreCut(indexfatJetCSVDiscMin_)  || jetCSVDisc > cut(indexfatJetCSVDiscMin_, double() ) ) passCut( ret ,indexfatJetCSVDiscMin_) ;
      if ( ignoreCut(indexfatJetCSVDiscMax_)  || jetCSVDisc < cut(indexfatJetCSVDiscMax_, double() ) ) passCut( ret ,indexfatJetCSVDiscMax_) ;

      if ( ignoreCut(indexsubjet1CSVDiscMin_) || subjet1CSVDisc > cut(indexsubjet1CSVDiscMin_, double() ) ) passCut( ret ,indexsubjet1CSVDiscMin_) ;
      if ( ignoreCut(indexsubjet1CSVDiscMax_) || subjet1CSVDisc < cut(indexsubjet1CSVDiscMax_, double() ) ) passCut( ret ,indexsubjet1CSVDiscMax_) ;
      if ( ignoreCut(indexsubjet2CSVDiscMin_) || subjet2CSVDisc > cut(indexsubjet2CSVDiscMin_, double() ) ) passCut( ret ,indexsubjet2CSVDiscMin_) ;
      if ( ignoreCut(indexsubjet2CSVDiscMax_) || subjet2CSVDisc < cut(indexsubjet2CSVDiscMax_, double() ) ) passCut( ret ,indexsubjet2CSVDiscMax_) ;

      setIgnored( ret ) ; 
      return (bool)ret ; 
    }

  private:

    JETTYPES_t jettype_ ; 

    edm::ParameterSet jetSelParams_ ; 
    index_type indexfatJetPtMin_ ; 
    index_type indexfatJetPtMax_ ; 
    index_type indexfatJetAbsEtaMax_ ;
    index_type indexfatJetMassMin_ ;
    index_type indexfatJetMassMax_ ;
    index_type indexfatJetPrunedMassMin_ ;
    index_type indexfatJetPrunedMassMax_ ;
    index_type indexfatJetTau2ByTau1Max_ ;
    index_type indexdRSubjetsMin_; 
    index_type indexdRSubjetsMax_; 
    index_type indexfatJetCSVDiscMin_; 
    index_type indexfatJetCSVDiscMax_; 
    index_type indexsubjet1CSVDiscMin_ ;
    index_type indexsubjet1CSVDiscMax_ ;
    index_type indexsubjet2CSVDiscMin_ ;
    index_type indexsubjet2CSVDiscMax_ ;

} ;

#endif 


