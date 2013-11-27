#ifndef BPRIMETOBH_INTERFACE_HTSELECTOR_H
#define BPRIMETOBH_INTERFACE_HTSELECTOR_H

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "BpbH/BprimeTobH/interface/HT.h"

class HTSelector : public Selector<HT> {
  public:
    HTSelector () {}

    ~HTSelector () {} 

    HTSelector (edm::ParameterSet const& params) {
      push_back("HTMin") ; 
      push_back("HTMax") ; 

      set("HTMin" ,params.getParameter<double>("HTMin") ) ;
      set("HTMax" ,params.getParameter<double>("HTMax") ) ;

      indexHTMin_            = index_type(&bits_ ,"HTMin") ;
      indexHTMax_            = index_type(&bits_ ,"HTMax") ;

      retInternal_ = getBitTemplate();   

    }

    // 
    // Accessor from PAT jets
    // 
    bool operator()( HT & HTIn, pat::strbitset & ret ) {
      return firstDataCuts( HTIn, ret );
    }

    // 
    // Accessor from *CORRECTED* 4-vector, EMF, and Jet ID. 
    // This can be used with reco quantities. 
    // 
    bool operator()( HT const & HTIn, pat::strbitset & ret ) {
      return firstDataCuts( HTIn, ret );
    }

    using Selector<HT>::operator();

    bool firstDataCuts( const HT & HTIn, pat::strbitset & ret) {
      ret.set(false);

      double thisHT = HTIn.getHT() ; 

      if ( ignoreCut(indexHTMin_) || thisHT > cut(indexHTMin_, double() ) ) passCut( ret ,indexHTMin_) ; 
      if ( ignoreCut(indexHTMax_) || thisHT < cut(indexHTMax_, double() ) ) passCut( ret ,indexHTMax_) ; 

      setIgnored( ret ) ; 
      return (bool)ret ; 
    }

  private:
    index_type indexHTMin_ ; 
    index_type indexHTMax_ ; 
}; 

#endif 

