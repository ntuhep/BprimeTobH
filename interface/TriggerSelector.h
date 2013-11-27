#ifndef BPRIMETOBH_INTERFACE_TRIGGERSELECTOR_H
#define BPRIMETOBH_INTERFACE_TRIGGERSELECTOR_H

#include "BpbH/BprimeTobH/interface/format.h"
#include "BpbH/BprimeTobH/interface/TriggerBooking.h"

class TriggerSelector {

  public:
    TriggerSelector (EventInfoBramches& evtinfo, edm::ParamaterSet const& pset) : 
      EvtInfo_(evtinfo), 
      hltPaths_(pset.getParameter<std::vector<int> >("HLTPaths")) { 
      } 

  private:
    EventInfoBranches EvtInfo_ ; 
    std::vector<std::int> hltPaths_ ; 

};

#endif 
