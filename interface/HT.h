#ifndef BPRIMETOBH_INTERFACE_HT_H
#define BPRIMETOBH_INTERFACE_HT_H

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "BpbH/BprimeTobH/interface/format.h"
#include "BpbH/BprimeTobH/interface/JetCollection.h"

#include <vector>

class HT {

  public :

    HT () : HT_(0) { ; } 
    HT (HT&) { ; }  
    ~HT() { ; } 

    void setJetCollection (JetCollection& JetColl) { vJets_.push_back(JetColl) ; } 

    void buildHT() { 

      if ( vJets_.size() == 0 ) HT_ = 0; 

      for (std::vector<JetCollection>::const_iterator itjetcolls =  vJets_.begin(); itjetcolls != vJets_.end(); ++itjetcolls) { 
        for (JetCollection::const_iterator itjet = itjetcolls->begin(); itjet != itjetcolls->end(); ++itjet) { 
          HT_ += itjet->Pt() ; 
        }
      }

    } 

    double getHT () const { return HT_ ; } 

  private: 

    double HT_ ; 
    std::vector<JetCollection> vJets_ ; 

};


#endif 

