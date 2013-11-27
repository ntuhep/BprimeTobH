#ifndef BPRIMETOBH_INTERFACE_JETSELECTOR_H
#define BPRIMETOBH_INTERFACE_JETSELECTOR_H

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "BpbH/BprimeTobH/interface/format.h"
#include "BpbH/BprimeTobH/interface/JetID.h"
#include <boost/algorithm/string.hpp>
#include <string>

class JetSelector : public Selector<int> {
  public:

    enum JETTYPES_t { AK5JET, BTAGGEDAK5JET, PRUNEDSUBJET, N_JETTYPES} ; 

    JetSelector () {} 

    JetSelector (edm::ParameterSet const& params) {

      std::string jettypeStr = params.getParameter<std::string>("jettype") ;
      if (jettypeStr == "AK5JET") jettype_ = AK5JET; 
      else if (jettypeStr == "BTAGGEDAK5JET") jettype_ = BTAGGEDAK5JET;
      else if (jettypeStr == "PRUNEDSUBJET") jettype_ = PRUNEDSUBJET; 
      else edm::LogError("WrongJetType") << " >>>> Check jet type !!!" ; 

      push_back("jetPtMin") ;
      push_back("jetPtMax") ;
      push_back("jetAbsEtaMax") ;
      push_back("jetCSVDiscMin") ;
      push_back("jetCSVDiscMax") ;
      push_back("IsJetIDLoose") ; 
      push_back("IsJetIDTight") ; 

      set("jetPtMin"             ,params.getParameter<double>("jetPtMin") ) ;
      set("jetPtMax"             ,params.getParameter<double>("jetPtMax") ) ;
      set("jetAbsEtaMax"         ,params.getParameter<double>("jetAbsEtaMax") ) ;
      set("jetCSVDiscMin"        ,params.getParameter<double>("jetCSVDiscMin") ) ;
      set("jetCSVDiscMax"        ,params.getParameter<double>("jetCSVDiscMax") ) ;
      set("IsJetIDLoose"         ,params.getParameter<bool>("IsJetIDLoose") ) ; 
      set("IsJetIDTight"         ,params.getParameter<bool>("IsJetIDTight") ) ; 

      indexjetPtMin_            = index_type(&bits_ ,"jetPtMin") ;
      indexjetPtMax_            = index_type(&bits_ ,"jetPtMax") ;
      indexjetAbsEtaMax_        = index_type(&bits_ ,"jetAbsEtaMax") ;
      indexjetCSVDiscMin_       = index_type(&bits_ ,"jetCSVDiscMin") ;
      indexjetCSVDiscMax_       = index_type(&bits_ ,"jetCSVDiscMax") ;
      indexIsJetIDLoose_        = index_type(&bits_ ,"IsJetIDLoose") ; 
      indexIsJetIDTight_        = index_type(&bits_ ,"IsJetIDTight") ; 

      retInternal_ = getBitTemplate();   

      if (params.getParameter<bool>("IsJetIDLoose") == true && params.getParameter<bool>("IsJetIDTight") == false) 
        jetID_ = new JetID(JetID::FIRSTDATA,JetID::LOOSE, jetInfo_) ;
      else if (params.getParameter<bool>("IsJetIDTight") == true && params.getParameter<bool>("IsJetIDLoose") == false) 
        jetID_ = new JetID(JetID::FIRSTDATA,JetID::TIGHT, jetInfo_) ;
      else edm::LogError("JetID") << "Ambiguous JetID: Please select only one (LOOSE or TIGHT) as True!!!" ; 
      retjetid_ = jetID_->getBitTemplate() ;

    }

    // 
    // Accessor from PAT jets
    // 
    bool operator()( int & jet, pat::strbitset & ret ) {
      return true;
    }
    using Selector<int>::operator();

    // 
    // Accessor from *CORRECTED* 4-vector, EMF, and Jet ID. 
    // This can be used with reco quantities. 
    // 
    bool operator()( int const & jet, pat::strbitset & ret ) {
      return true;
    }

    bool firstDataCuts( const int & jet, pat::strbitset & ret) {
      ret.set(false);

      double jetPt         = jetInfo_.Pt[jet];
      double jetEta        = jetInfo_.Eta[jet];
      double jetCSVDisc    = jetInfo_.CombinedSVBJetTags[jet]; 
      JetID id(*jetID_) ;
      retjetid_.set(false);
      bool   jetid         = id(jet, retjetid_) ; 

      if ( ignoreCut(indexjetPtMin_)            || jetPt > cut(indexjetPtMin_, int() ) ) passCut( ret ,indexjetPtMin_) ;
      if ( ignoreCut(indexjetPtMax_)            || jetPt < cut(indexjetPtMax_, int() ) ) passCut( ret ,indexjetPtMax_) ;
      if ( ignoreCut(indexjetAbsEtaMax_)        || jetEta < cut(indexjetAbsEtaMax_, int() ) ) passCut( ret ,indexjetAbsEtaMax_) ;
      if ( ignoreCut(indexIsJetIDLoose_)        || jetid == cut(indexIsJetIDLoose_, int() ) ) passCut( ret ,indexIsJetIDLoose_) ; 
      if ( ignoreCut(indexIsJetIDTight_)        || jetid == cut(indexIsJetIDTight_, int() ) ) passCut( ret ,indexIsJetIDTight_) ; 

      if (jettype_ == BTAGGEDAK5JET) {
        if ( ignoreCut(indexjetCSVDiscMin_)       || jetCSVDisc > cut(indexjetCSVDiscMin_, int() ) ) passCut( ret ,indexjetCSVDiscMin_) ;
        if ( ignoreCut(indexjetCSVDiscMax_)       || jetCSVDisc < cut(indexjetCSVDiscMax_, int() ) ) passCut( ret ,indexjetCSVDiscMax_) ;
      }

      setIgnored( ret ) ; 
      return (bool)ret ; 
    }

  private:

    JETTYPES_t jettype_ ; 
    JetInfoBranches jetInfo_ ; 

    JetID* jetID_ ; 
    pat::strbitset retjetid_ ; 

    index_type indexjetPtMin_ ;
    index_type indexjetPtMax_ ;
    index_type indexjetAbsEtaMax_ ;
    index_type indexjetCSVDiscMin_; 
    index_type indexjetCSVDiscMax_; 
    index_type indexIsJetIDLoose_ ; 
    index_type indexIsJetIDTight_ ; 

};

#endif 
