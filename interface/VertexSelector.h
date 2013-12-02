#ifndef BPRIMETOBH_INTERFACE_VERTEXSELECTOR_H
#define BPRIMETOBH_INTERFACE_VERTEXSELECTOR_H

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "BpbH/BprimeTobH/interface/VertexCollection.h"

class VertexSelector {
  public:

    VertexSelector (const VertexInfoBranches& VtxInfo) {
      for (int ivtx = 0; ivtx < VtxInfo.Size; ++ivtx) { 
        Vertex vtx(VtxInfo, ivtx) ; 
        vertices_.push_back(vtx) ; 
      }
    }

    int NGoodVtxs () const { 
      int nGoodVtxs(0) ; 
      for (VertexCollection::const_iterator itvtx = vertices_.begin(); itvtx != vertices_.end(); ++itvtx) { 
        if ( itvtx->IsGoodOfflinePrimaryVtx() ) ++nGoodVtxs ; 
      }
      return nGoodVtxs ; 
    }

  private:
    VertexCollection vertices_ ; 
    
}; 

#endif 
