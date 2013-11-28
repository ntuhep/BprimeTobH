#ifndef BPRIMETOBH_INTERFACE_TRIGGERSELECTOR_H
#define BPRIMETOBH_INTERFACE_TRIGGERSELECTOR_H

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "BpbH/BprimeTobH/interface/format.h"
#include "BpbH/BprimeTobH/interface/TriggerBooking.h"
#include <string>

class TriggerSelector {

  public:
    TriggerSelector (edm::ParameterSet const& pset) : 
      hltPaths_(pset.getParameter<std::vector<std::string> >("HLTPaths")),
      trigDecision_(false) { 
      } 

    bool getTrigDecision (EvtInfoBranches& evtinfo) {

      trigDecision_ = false;
      if(hltPaths_.empty()) {
        edm::LogWarning("TriggerSelector") << " >>>> Warning: No HLT paths specified. Event trigDecision_es trigger. " ; 
        trigDecision_ = true ; 
        return trigDecision_ ; 
      }

      for ( std::vector<std::string>::const_iterator ittrig = hltPaths_.begin(); ittrig != hltPaths_.end(); ++ittrig ) {
        if (trigDecision_) break ; 
        for( int ii = 0; ii < N_TRIGGER_BOOKINGS; ++ii ) {
          if(trigDecision_) break;
          std::string trigname = TriggerBooking[ii] ; 
          if( std::string::npos != trigname.find(*ittrig) ) {//found trigger 
            if( evtinfo.TrgBook[ii] == 1 ) {
              trigDecision_ = true ; 
            }
          }
        }
      }

      return trigDecision_ ; 

    }

  private:
    std::vector<std::string> hltPaths_ ; 
    bool trigDecision_ ; 

};

#endif 
