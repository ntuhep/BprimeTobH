#ifndef BPRIMETOBH_INTERFACE_JETSELECTOR_H
#define BPRIMETOBH_INTERFACE_JETSELECTOR_H

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "../../BprimeTobH/interface/format.h"
#include "../../BprimeTobH/interface/JetID.h"
#include <boost/algorithm/string.hpp>

class JetSelector : public Selector<int> {
  public:

    enum JETTYPES_t { AK5JET, CA8JET, PRUNEDSUBJET, N_JETTYPES} ; 

    JetSelector () {} 

    JetSelector (edm::ParameterSet const& params) {
      std::string jettypeStr = params.getParameter<std::string>("jettype") ;

      if (jettypeStr == "AK5JET") jettype_ = AK5JET; 
      else if (jettypeStr == "CA8JET") jettype_ = CA8JET;
      else if (jettypeStr == "PRUNEDSUBJET") jettype_ = PRUNEDSUBJET; 
      else edm::LogError("WrongJetType") << " Check jet type!!!!" ; 

      push_back("jetPtMin") ;
      push_back("jetPtMax") ;
      push_back("jetAbsEtaMax") ;
      push_back("bjetPtMin") ;
      push_back("fatJetPtMin") ;
      push_back("fatJetPtMax") ;
      push_back("fatJetAbsEtaMax") ;
      push_back("fatJetMassMin") ;
      push_back("fatJetMassMax") ;
      push_back("fatJetPrunedMassMin") ;
      push_back("fatJetPrunedMassMax") ;
      push_back("fatJetTau2ByTau1Max") ;
      push_back("subjet1CSVDiscMin") ;
      push_back("subjet1CSVDiscMax") ;
      push_back("subjet2CSVDiscMin") ;
      push_back("subjet2CSVDiscMax") ;
      push_back("JetID") ; 

      set("jetPtMin"             ,params.getParameter<double>("jetPtMin") ) ;
      set("jetPtMax"             ,params.getParameter<double>("jetPtMax") ) ;
      set("jetAbsEtaMax"         ,params.getParameter<double>("jetAbsEtaMax") ) ;
      set("bjetPtMin"            ,params.getParameter<double>("bjetPtMin") ) ;
      set("fatJetPtMin"          ,params.getParameter<double>("fatJetPtMin") ) ;
      set("fatJetPtMax"          ,params.getParameter<double>("fatJetPtMax") ) ;
      set("fatJetAbsEtaMax"      ,params.getParameter<double>("fatJetAbsEtaMax") ) ;
      set("fatJetMassMin"        ,params.getParameter<double>("fatJetMassMin") ) ;
      set("fatJetMassMax"        ,params.getParameter<double>("fatJetMassMax") ) ;
      set("fatJetPrunedMassMin"  ,params.getParameter<double>("fatJetPrunedMassMin") ) ;
      set("fatJetPrunedMassMax"  ,params.getParameter<double>("fatJetPrunedMassMax") ) ;
      set("fatJetTau2ByTau1Max"  ,params.getParameter<double>("fatJetTau2ByTau1Max") ) ;
      set("subjet1CSVDiscMin"    ,params.getParameter<double>("subjet1CSVDiscMin") ) ;
      set("subjet1CSVDiscMax"    ,params.getParameter<double>("subjet1CSVDiscMax") ) ;
      set("subjet2CSVDiscMin"    ,params.getParameter<double>("subjet2CSVDiscMin") ) ;
      set("subjet2CSVDiscMax"    ,params.getParameter<double>("subjet2CSVDiscMax") ) ;
      set("JetID"                ,params.getParameter<double>("JetID") ) ; 

      indexjetPtMin_            = index_type(&bits_ ,"jetPtMin") ;
      indexjetPtMax_            = index_type(&bits_ ,"jetPtMax") ;
      indexjetAbsEtaMax_        = index_type(&bits_ ,"jetAbsEtaMax") ;
      indexbjetPtMin_           = index_type(&bits_ ,"bjetPtMin") ;
      indexfatJetPtMin_         = index_type(&bits_ ,"fatJetPtMin") ;
      indexfatJetPtMax_         = index_type(&bits_ ,"fatJetPtMax") ;
      indexfatJetAbsEtaMax_     = index_type(&bits_ ,"fatJetAbsEtaMax") ;
      indexfatJetMassMin_       = index_type(&bits_ ,"fatJetMassMin") ;
      indexfatJetMassMax_       = index_type(&bits_ ,"fatJetMassMax") ;
      indexfatJetPrunedMassMin_ = index_type(&bits_ ,"fatJetPrunedMassMin") ;
      indexfatJetPrunedMassMax_ = index_type(&bits_ ,"fatJetPrunedMassMax") ;
      indexfatJetTau2ByTau1Max_ = index_type(&bits_ ,"fatJetTau2ByTau1Max") ;
      indexsubjet1CSVDiscMin_   = index_type(&bits_ ,"subjet1CSVDiscMin") ;
      indexsubjet1CSVDiscMax_   = index_type(&bits_ ,"subjet1CSVDiscMax") ;
      indexsubjet2CSVDiscMin_   = index_type(&bits_ ,"subjet2CSVDiscMin") ;
      indexsubjet2CSVDiscMax_   = index_type(&bits_ ,"subjet2CSVDiscMax") ;
      indexJetID_               = index_type(&bits_ ,"JetID") ; 

      retInternal_ = getBitTemplate();   

      if (boost::iequals(indexJetID_.str(),"LOOSE"))      jetID_ = new JetID(JetID::FIRSTDATA,JetID::LOOSE, jetInfo_) ;
      else if (boost::iequals(indexJetID_.str(),"TIGHT")) jetID_ = new JetID(JetID::FIRSTDATA,JetID::TIGHT, jetInfo_) ;
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
      double jetMass       = jetInfo_.Mass[jet];
      double jetMassPruned = jetInfo_.MassPruned[jet];
      double jetTau1       = jetInfo_.tau1[jet];
      double jetTau2       = jetInfo_.tau2[jet];
      double jetCSVDisc    = jetInfo_.CombinedSVBJetTags[jet];
      JetID id(*jetID_) ;
      retjetid.set(false);
      bool   jetid         = id(jet, retjetid_) ; 

      if ( ignoreCut(indexjetPtMin_)            || jetPt > cut(indexjetPtMin_, int() ) ) passCut( ret ,indexjetPtMin_) ;
      if ( ignoreCut(indexjetPtMax_)            || jetPt < cut(indexjetPtMax_, int() ) ) passCut( ret ,indexjetPtMax_) ;
      if ( ignoreCut(indexjetAbsEtaMax_)        || jetEta < cut(indexjetAbsEtaMax_, int() ) ) passCut( ret ,indexjetAbsEtaMax_) ;
      if ( ignoreCut(indexfatJetMassMin_)       || jetMass > cut(indexfatJetMassMin_, int() ) ) passCut( ret ,indexfatJetMassMin_) ;
      if ( ignoreCut(indexfatJetMassMax_)       || jetMass < cut(indexfatJetMassMax_, int() ) ) passCut( ret ,indexfatJetMassMax_) ;
      if ( ignoreCut(indexfatJetPrunedMassMin_) || jetMassPruned > cut(indexfatJetPrunedMassMin_, int() ) ) passCut( ret ,indexfatJetPrunedMassMin_) ;
      if ( ignoreCut(indexfatJetPrunedMassMax_) || jetMassPruned < cut(indexfatJetPrunedMassMax_, int() ) ) passCut( ret ,indexfatJetPrunedMassMax_) ;
      if ( ignoreCut(indexfatJetTau2ByTau1Max_) || jetTau2/jetTau1 < cut(indexfatJetTau2ByTau1Max_, int() ) ) passCut( ret ,indexfatJetTau2ByTau1Max_) ;
      if ( ignoreCut(indexsubjet1CSVDiscMin_)   || jetCSVDisc > cut(indexsubjet1CSVDiscMin_, int() ) ) passCut( ret ,indexsubjet1CSVDiscMin_) ;
      if ( ignoreCut(indexsubjet1CSVDiscMax_)   || jetCSVDisc < cut(indexsubjet1CSVDiscMax_, int() ) ) passCut( ret ,indexsubjet1CSVDiscMax_) ;
      if ( ignoreCut(indexJetID_)               || jetid == cut(indexJetID_, int() ) ) passCut( ret ,indexJetID_) ; 

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
    index_type indexbjetPtMin_ ;
    index_type indexfatJetPtMin_ ;
    index_type indexfatJetPtMax_ ;
    index_type indexfatJetAbsEtaMax_ ;
    index_type indexfatJetMassMin_ ;
    index_type indexfatJetMassMax_ ;
    index_type indexfatJetPrunedMassMin_ ;
    index_type indexfatJetPrunedMassMax_ ;
    index_type indexfatJetTau2ByTau1Max_ ;
    index_type indexsubjet1CSVDiscMin_ ;
    index_type indexsubjet1CSVDiscMax_ ;
    index_type indexsubjet2CSVDiscMin_ ;
    index_type indexsubjet2CSVDiscMax_ ;
    index_type indexJetID_ ; 

};

#endif 
