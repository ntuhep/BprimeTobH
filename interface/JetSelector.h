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

      set("jetPtMin"             ,params.getParameter<double>("jetPtMin") ) ;
      set("jetPtMax"             ,params.getParameter<double>("jetPtMax") ) ;
      set("jetAbsEtaMax"         ,params.getParameter<double>("jetAbsEtaMax") ) ;
      if (jettype_ == BTAGGEDAK5JET) {
        set("jetCSVDiscMin"        ,params.getParameter<double>("jetCSVDiscMin") ,true) ;
        set("jetCSVDiscMax"        ,params.getParameter<double>("jetCSVDiscMax") ,true) ;
      }
      else {
        set("jetCSVDiscMin"        ,params.getParameter<double>("jetCSVDiscMin") ,false) ;
        set("jetCSVDiscMax"        ,params.getParameter<double>("jetCSVDiscMax") ,false) ;
      }

      indexjetPtMin_            = index_type(&bits_ ,"jetPtMin") ;
      indexjetPtMax_            = index_type(&bits_ ,"jetPtMax") ;
      indexjetAbsEtaMax_        = index_type(&bits_ ,"jetAbsEtaMax") ;
      indexjetCSVDiscMin_       = index_type(&bits_ ,"jetCSVDiscMin") ;
      indexjetCSVDiscMax_       = index_type(&bits_ ,"jetCSVDiscMax") ;

      retInternal_ = getBitTemplate();   

      if (params.getParameter<bool>("IsJetIDLoose") == true && params.getParameter<bool>("IsJetIDTight") == false) quality_ = JetID::LOOSE ; 
      else if (params.getParameter<bool>("IsJetIDTight") == true && params.getParameter<bool>("IsJetIDLoose") == false) quality_ = JetID::TIGHT ; 
      else edm::LogError("JetID") << "Ambiguous JetID: Please select only one (LOOSE or TIGHT) as True!!!" ; 

    }

    bool operator()(int & jet, pat::strbitset & ret ) {
      return true ;  
    }

    bool operator()(int const  & jet, pat::strbitset & ret ) {
      return true ;  
    }

    using Selector<int>::operator();

    bool operator()(JetInfoBranches& jetInfo, int const & jet, pat::strbitset & ret ) {
      return  firstDataCuts(jetInfo, jet, ret ) ; 
    }

    bool firstDataCuts(JetInfoBranches& jetInfo,  const int & jet, pat::strbitset & ret) {

      ret.set(false);
      bool isJetID(false) ; 
      JetID jetID(JetID::FIRSTDATA,quality_, jetInfo) ; 
      pat::strbitset retjetid = jetID.getBitTemplate() ;
      retjetid.set(false) ;
      if (jetID(jetInfo, jet,retjetid) != 0) isJetID = true ;

      if (isJetID == false) return false ; 

      double jetPt         = jetInfo.Pt[jet];
      double jetAbsEta     = std::abs(jetInfo.Eta[jet]);
      double jetCSVDisc    = jetInfo.CombinedSVBJetTags[jet]; 

      if ( ignoreCut(indexjetPtMin_)        || jetPt > cut(indexjetPtMin_, double() ) ) passCut( ret ,indexjetPtMin_) ; 
      if ( ignoreCut(indexjetPtMax_)        || jetPt < cut(indexjetPtMax_, double() ) ) passCut( ret ,indexjetPtMax_) ; 
      if ( ignoreCut(indexjetAbsEtaMax_)    || jetAbsEta < cut(indexjetAbsEtaMax_, double() ) ) passCut( ret ,indexjetAbsEtaMax_) ; 

      if ( ignoreCut(indexjetCSVDiscMin_) || jetCSVDisc > cut(indexjetCSVDiscMin_, double() ) ) passCut( ret ,indexjetCSVDiscMin_) ;
      if ( ignoreCut(indexjetCSVDiscMax_) || jetCSVDisc < cut(indexjetCSVDiscMax_, double() ) ) passCut( ret ,indexjetCSVDiscMax_) ;

      setIgnored( ret ) ; 
      return (bool)ret ; 
    }

  private:

    JETTYPES_t jettype_ ; 

    JetID* jetID_ ; 
    pat::strbitset retjetid_ ; 
    JetID::Quality_t quality_; 

    index_type indexjetPtMin_ ;
    index_type indexjetPtMax_ ;
    index_type indexjetAbsEtaMax_ ;
    index_type indexjetCSVDiscMin_; 
    index_type indexjetCSVDiscMax_; 

};

#endif 
