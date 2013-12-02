#ifndef BPRIMETOBH_INTERFACE_VERTEX_H
#define BPRIMETOBH_INTERFACE_VERTEX_H

#include "BpbH/BprimeTobH/interface/format.h"

class Vertex {

  public:

    Vertex(const VertexInfoBranches& VtxInfo, const int& VtxIndex) :
      isValid_(VtxInfo.isValid[VtxIndex]), 
      isFake_(VtxInfo.isFake[VtxIndex]),  
      Type_(VtxInfo.Type[VtxIndex]),    //0 - Offline Primary Vertices, 1 - Offline Primary Vertices with beam spot constraint, 2 - Pixel Vertices
      Ndof_(VtxInfo.Ndof[VtxIndex]), 
      NormalizedChi2_(VtxInfo.NormalizedChi2[VtxIndex]), 
      Pt_Sum_(VtxInfo.Pt_Sum[VtxIndex]), 
      Pt_Sum2_(VtxInfo.Pt_Sum2[VtxIndex]), 
      x_(VtxInfo.x[VtxIndex]), 
      y_(VtxInfo.y[VtxIndex]), 
      z_(VtxInfo.z[VtxIndex]), 
      Rho_(VtxInfo.Rho[VtxIndex]),
      isGoodOfflinePrimaryVtx_(false) 
  { 
    isGoodOfflinePrimaryVtx() ; 
  }

    ~Vertex () {} ;

    int     isValid() const { return isValid_ ; } ;
    bool    isFake() const { return  isFake_ ; } ; 
    int     Type() const { return Type_ ; } ;   //0 - Offline Primary Vertices, 1 - Offline Primary Vertices with beam spot constraint, 2 - Pixel Vertices
    float   Ndof() const { return Ndof_ ; } ;
    float   NormalizedChi2() const { return NormalizedChi2_ ; } ;
    float   Pt_Sum() const { return  Pt_Sum_ ; } ;
    float   Pt_Sum2() const { return Pt_Sum2_ ; } ;
    float   x() const { return x_ ; } ;
    float   y() const { return y_ ; } ;
    float   z() const { return z_ ; } ;
    float   Rho() const { return Rho_ ; } ;
    bool    IsGoodOfflinePrimaryVtx() const { return isGoodOfflinePrimaryVtx_ ; }

  private:

    void isGoodOfflinePrimaryVtx () { 
      if (   Type_ == 1 
          && isFake_ == false
          && Ndof_ > 4
          && Rho_ < 2.
          && z_ < 24.) isGoodOfflinePrimaryVtx_ = true ; 
    }

    int     isValid_ ;
    bool    isFake_ ; 
    int     Type_ ;   //0 - Offline Primary Vertices, 1 - Offline Primary Vertices with beam spot constraint, 2 - Pixel Vertices
    float   Ndof_ ;
    float   NormalizedChi2_ ;
    float   Pt_Sum_ ;
    float   Pt_Sum2_ ;
    float   x_ ;
    float   y_ ;
    float   z_ ;
    float   Rho_ ;
    bool    isGoodOfflinePrimaryVtx_ ; 
};

#endif 
