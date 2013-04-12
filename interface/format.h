#ifndef FORMAT_H
#define FORMAT_H


#include <TTree.h>

static const int MAX_VERTICES=256; 
static const int MAX_LEPTONS=256; 
static const int MAX_JETS=256; 

class EvtInfoBranches {
 public:
  int      RunNo;
  long long int EvtNo;
  
  //BeamSpot
  float BeamSpotX;
  float BeamSpotY;
  float BeamSpotZ;

  
  void RegisterTree(TTree *tree) {
    tree->Branch("EvtInfo.RunNo" , &RunNo, "EvtInfo.RunNo/I");
    tree->Branch("EvtInfo.EvtNo" , &EvtNo, "EvtInfo.EvtNo/L");
    tree->Branch("EvtInfo.BeamSpotX", &BeamSpotX, "EvtInfo.BeamSpotX/F");
    tree->Branch("EvtInfo.BeamSpotY", &BeamSpotY, "EvtInfo.BeamSpotY/F");
    tree->Branch("EvtInfo.BeamSpotZ", &BeamSpotZ, "EvtInfo.BeamSpotZ/F");
  }
  

  void Register(TTree *tree) {
    tree->SetBranchAddress("EvtInfo.RunNo" , &RunNo);
    tree->SetBranchAddress("EvtInfo.EvtNo" , &EvtNo);
    tree->SetBranchAddress("EvtInfo.BeamSpotX"    , &BeamSpotX	 );
    tree->SetBranchAddress("EvtInfo.BeamSpotY"    , &BeamSpotY	 );
    tree->SetBranchAddress("EvtInfo.BeamSpotZ"    , &BeamSpotZ	 );
  }
  
}; 

  
class VertexInfoBranches {
	public:
		int     Size;
		int     isValid[MAX_VERTICES];
		bool    isFake[MAX_VERTICES]; 
		int     Type[MAX_VERTICES];   //0 - Offline Primary Vertices, 1 - Offline Primary Vertices with beam spot constraint, 2 - Pixel Vertices
		float   Ndof[MAX_VERTICES];
		float   NormalizedChi2[MAX_VERTICES];
		float   Pt_Sum[MAX_VERTICES];
		float   Pt_Sum2[MAX_VERTICES];
		float   x[MAX_VERTICES];
		float   y[MAX_VERTICES];
		float   z[MAX_VERTICES];
		float   Rho[MAX_VERTICES];

		void RegisterTree(TTree *tree) {
		  tree->Branch("VertexInfo.Size", &Size, "VertexInfo.Size/I");
		  tree->Branch("VertexInfo.isValid", &isValid[0], "VertexInfo.isValid[VertexInfo.Size]/I");
		  tree->Branch("VertexInfo.isFake", &isFake[0], "VertexInfo.isFake[VertexInfo.Size]/O"); 
		  tree->Branch("VertexInfo.Type", &Type[0], "VertexInfo.Type[VertexInfo.Size]/I");
		  tree->Branch("VertexInfo.Ndof", &Ndof[0], "VertexInfo.Ndof[VertexInfo.Size]/F");
		  tree->Branch("VertexInfo.NormalizedChi2", &NormalizedChi2[0], "VertexInfo.NormalizedChi2[VertexInfo.Size]/F");
		  tree->Branch("VertexInfo.Pt_Sum", &Pt_Sum[0], "VertexInfo.Pt_Sum[VertexInfo.Size]/F");
		  tree->Branch("VertexInfo.Pt_Sum2", &Pt_Sum2[0], "VertexInfo.Pt_Sum2[VertexInfo.Size]/F");
		  tree->Branch("VertexInfo.x", &x[0], "VertexInfo.x[VertexInfo.Size]/F");
		  tree->Branch("VertexInfo.y", &y[0], "VertexInfo.y[VertexInfo.Size]/F"	);
		  tree->Branch("VertexInfo.z", &z[0], "VertexInfo.z[VertexInfo.Size]/F"	);
		  tree->Branch("VertexInfo.Rho", &Rho[0] , "VertexInfo.Rho[VertexInfo.Size]/F");
		}										    
		void Register(TTree *tree) {
		  tree->SetBranchAddress("VertexInfo.Size", &Size);
		  tree->SetBranchAddress("VertexInfo.isValid", &isValid[0]);
		  tree->SetBranchAddress("VertexInfo.isFake", &isFake[0]);
		  tree->SetBranchAddress("VertexInfo.Type", &Type[0]);
		  tree->SetBranchAddress("VertexInfo.Ndof", &Ndof[0]);
		  tree->SetBranchAddress("VertexInfo.NormalizedChi2", &NormalizedChi2[0]);
		  tree->SetBranchAddress("VertexInfo.Pt_Sum", &Pt_Sum[0]);
		  tree->SetBranchAddress("VertexInfo.Pt_Sum2", &Pt_Sum2[0]);
		  tree->SetBranchAddress("VertexInfo.x" , &x[0]);
		  tree->SetBranchAddress("VertexInfo.y"        , &y[0]  	 );
		  tree->SetBranchAddress("VertexInfo.z"        , &z[0]  	 );
		  tree->SetBranchAddress("VertexInfo.Rho"        , &Rho[0]  	 );
		}  	
		
}; 


class LepInfoBranches {
 public:
  int	Size; 
  int	Index[MAX_LEPTONS];
  
  void RegisterTree(TTree *tree, std::string name="LepInfo") {
    tree->Branch((name+".Size").c_str()	, &Size , (name+".Size/I").c_str() );
    tree->Branch((name+".Index").c_str(), &Index[0], (name+".Index["+name+".Size]/I").c_str());
  }
    
  void Register(TTree *tree, std::string name="LepInfo") {
      tree->SetBranchAddress((name+".Size").c_str(), &Size);
      tree->SetBranchAddress((name+".Index").c_str(), &Index[0] );
  }

      
};


class JetInfoBranches {
 public:
  int   Size; 
  int   Index[MAX_JETS];

  void RegisterTree(TTree *root, std::string name="JetInfo") {
    root->Branch((name+".Size").c_str(), &Size, (name+".Size/I").c_str() );
    root->Branch((name+".Index").c_str(), &Index[0], 
		 (name+".Index["+name+".Size]/I").c_str());

  }
  
  void Register(TTree *root, std::string name="JetInfo") {
    root->SetBranchAddress((name+".Size").c_str() , &Size);
    root->SetBranchAddress((name+".Index").c_str(), &Index[0]);
  }

};





#endif
    
