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
  int   NTracks[MAX_JETS];
  float Et[MAX_JETS];
  float Pt[MAX_JETS];
  float Unc[MAX_JETS];
  float Eta[MAX_JETS];
  float Phi[MAX_JETS];

  float GenJetPt[MAX_JETS];
  float GenJetEta[MAX_JETS];
  float GenJetPhi[MAX_JETS];
  float GenPt[MAX_JETS];
  float GenEta[MAX_JETS];
  float GenPhi[MAX_JETS];
  int   GenPdgID[MAX_JETS];
  int   GenFlavor[MAX_JETS];
  int   GenMCTag[MAX_JETS]; // 0: unknown, 1: decay from W, 2: decay from Z, (+10) from b', (+20) from t'

  bool  JetIDLOOSE[MAX_JETS]; //Add by Chiyi
  float JetCharge[MAX_JETS];
  int   NConstituents[MAX_JETS];
  int   NCH[MAX_JETS];
  float CEF[MAX_JETS];
  float NHF[MAX_JETS];
  float NEF[MAX_JETS];
  float CHF[MAX_JETS];
  float JVAlpha[MAX_JETS];
  float JVBeta[MAX_JETS];
  float PtCorrRaw[MAX_JETS];  
  float PtCorrL2[MAX_JETS];  
  float PtCorrL3[MAX_JETS];  
  float PtCorrL7g[MAX_JETS];
  float PtCorrL7uds[MAX_JETS];
  float PtCorrL7c[MAX_JETS];  
  float PtCorrL7b[MAX_JETS];  

  float JetBProbBJetTags[MAX_JETS];
  float JetProbBJetTags[MAX_JETS];
  float TrackCountHiPurBJetTags[MAX_JETS];  
  float TrackCountHiEffBJetTags[MAX_JETS]; 
  float SimpleSVBJetTags[MAX_JETS];  //for 35X sample //Add by Chiyi
  float SimpleSVHEBJetTags[MAX_JETS];  //for 36X sample //Add by Chiyi
  float SimpleSVHPBJetTags[MAX_JETS];  //for 36X sample //Add by Chiyi
  float CombinedSVBJetTags[MAX_JETS];
  float CombinedSVMVABJetTags[MAX_JETS];
  float SoftElecByIP3dBJetTags[MAX_JETS];
  float SoftElecByPtBJetTags[MAX_JETS];  
  float SoftMuonBJetTags[MAX_JETS];      
  float SoftMuonByIP3dBJetTags[MAX_JETS];
  float SoftMuonByPtBJetTags[MAX_JETS];  
  float DoubleSVHighEffBJetTags[MAX_JETS]; //// Added by DM 

  float Px[MAX_JETS]; //Uly 2011-04-04
  float Py[MAX_JETS]; //Uly 2011-04-04
  float Pz[MAX_JETS]; //Uly 2011-04-04
  float Energy[MAX_JETS]; //Uly 2011-04-04

  float Mass[MAX_JETS];
  float Area[MAX_JETS];

  float MassD1[MAX_JETS];
  float MassD2[MAX_JETS];
  float PtD1[MAX_JETS];
  float PtD2[MAX_JETS];
  float EtD1[MAX_JETS];
  float EtD2[MAX_JETS];
  float EtaD1[MAX_JETS];
  float EtaD2[MAX_JETS];
  float PhiD1[MAX_JETS];
  float PhiD2[MAX_JETS];

  void RegisterTree(TTree *root, std::string name="JetInfo") {
    root->Branch((name+".Size").c_str(), &Size, (name+".Size/I").c_str() );
    root->Branch((name+".Index").c_str(), &Index[0], 
		 (name+".Index["+name+".Size]/I").c_str());
    root->Branch((name+".NTracks").c_str()  	       , &NTracks[0]		     , (name+".NTracks["+name+".Size]/I").c_str()		);
    root->Branch((name+".Et").c_str()		       , &Et[0] 		     , (name+".Et["+name+".Size]/F").c_str()			);
    root->Branch((name+".Pt").c_str()		       , &Pt[0] 		     , (name+".Pt["+name+".Size]/F").c_str()			);
    root->Branch((name+".Unc").c_str()		       , &Unc[0] 		     , (name+".Unc["+name+".Size]/F").c_str()			);
    root->Branch((name+".Eta").c_str()		       , &Eta[0]		     , (name+".Eta["+name+".Size]/F").c_str()			);
    root->Branch((name+".Phi").c_str()		       , &Phi[0]		     , (name+".Phi["+name+".Size]/F").c_str()			);
    root->Branch((name+".JetIDLOOSE").c_str()	       ,&JetIDLOOSE[0]	     , (name+".JetIDLOOSE["+name+".Size]/I").c_str()		); //Add by Chiyi
    root->Branch((name+".JetCharge").c_str()	       , &JetCharge[0]  	     , (name+".JetCharge["+name+".Size]/F").c_str()		);
    root->Branch((name+".NConstituents").c_str()	       , &NConstituents[0]	     , (name+".NConstituents["+name+".Size]/I").c_str()		);
    root->Branch((name+".NCH").c_str()	       , &NCH[0]	     , (name+".NCH["+name+".Size]/I").c_str()		);
    root->Branch((name+".CEF").c_str() 	       	       , &CEF[0]		     , (name+".CEF["+name+".Size]/F").c_str()		);
    root->Branch((name+".NHF").c_str() 	       	       , &NHF[0]		     , (name+".NHF["+name+".Size]/F").c_str()		);
    root->Branch((name+".NEF").c_str() 	       	       , &NEF[0]		     , (name+".NEF["+name+".Size]/F").c_str()		);
    root->Branch((name+".CHF").c_str() 	       	       , &CHF[0]		     , (name+".CHF["+name+".Size]/F").c_str()		);
    root->Branch((name+".JVAlpha").c_str() 	       	       , &JVAlpha[0]		     , (name+".JVAlpha["+name+".Size]/F").c_str()		);
    root->Branch((name+".JVBeta").c_str()	       	       , &JVBeta[0]	     	     , (name+".JVBeta["+name+".Size]/F").c_str()  		); 
    root->Branch((name+".PtCorrRaw").c_str() 	       , &PtCorrRaw[0]		     , (name+".PtCorrRaw["+name+".Size]/F").c_str()		);
    root->Branch((name+".PtCorrL2").c_str() 	       , &PtCorrL2[0]		     , (name+".PtCorrL2["+name+".Size]/F").c_str()		);
    root->Branch((name+".PtCorrL3").c_str() 	       , &PtCorrL3[0]		     , (name+".PtCorrL3["+name+".Size]/F").c_str()		);
    root->Branch((name+".PtCorrL7g").c_str() 	       , &PtCorrL7g[0]	  	     , (name+".PtCorrL7g["+name+".Size]/F").c_str()		);
    root->Branch((name+".PtCorrL7uds").c_str() 	       , &PtCorrL7uds[0]	     , (name+".PtCorrL7uds["+name+".Size]/F").c_str()		);	
    root->Branch((name+".PtCorrL7c").c_str() 	       , &PtCorrL7c[0]	  	     , (name+".PtCorrL7c["+name+".Size]/F").c_str()		);	
    root->Branch((name+".PtCorrL7b").c_str() 	       , &PtCorrL7b[0]	  	     , (name+".PtCorrL7b["+name+".Size]/F").c_str()		);	

    root->Branch((name+".JetBProbBJetTags").c_str()        , &JetBProbBJetTags[0]	     , (name+".JetBProbBJetTags["+name+".Size]/F").c_str()	);
    root->Branch((name+".JetProbBJetTags").c_str()	       , &JetProbBJetTags[0]	     , (name+".JetProbBJetTags["+name+".Size]/F").c_str()	);
    root->Branch((name+".TrackCountHiPurBJetTags").c_str() , &TrackCountHiPurBJetTags[0] , (name+".TrackCountHiPurBJetTags["+name+".Size]/F").c_str());	    
    root->Branch((name+".TrackCountHiEffBJetTags").c_str() , &TrackCountHiEffBJetTags[0] , (name+".TrackCountHiEffBJetTags["+name+".Size]/F").c_str());
    root->Branch((name+".SimpleSVBJetTags").c_str()        , &SimpleSVBJetTags[0]	     , (name+".SimpleSVBJetTags["+name+".Size]/F").c_str()	); 
    root->Branch((name+".SimpleSVHEBJetTags").c_str()        , &SimpleSVHEBJetTags[0]	     , (name+".SimpleSVHEBJetTags["+name+".Size]/F").c_str()	); 
    root->Branch((name+".SimpleSVHPBJetTags").c_str()        , &SimpleSVHPBJetTags[0]	     , (name+".SimpleSVHPBJetTags["+name+".Size]/F").c_str()	); 
    root->Branch((name+".CombinedSVBJetTags").c_str()      , &CombinedSVBJetTags[0]      , (name+".CombinedSVBJetTags["+name+".Size]/F").c_str()	);
    root->Branch((name+".CombinedSVMVABJetTags").c_str()   , &CombinedSVMVABJetTags[0]   , (name+".CombinedSVMVABJetTags["+name+".Size]/F").c_str()  );
    root->Branch((name+".SoftElecByIP3dBJetTags").c_str()  , &SoftElecByIP3dBJetTags[0]  , (name+".SoftElecByIP3dBJetTags["+name+".Size]/F").c_str()	);
    root->Branch((name+".SoftElecByPtBJetTags").c_str()    , &SoftElecByPtBJetTags[0]    , (name+".SoftElecByPtBJetTags["+name+".Size]/F").c_str()	);
    root->Branch((name+".SoftMuonBJetTags").c_str()        , &SoftMuonBJetTags[0]	     , (name+".SoftMuonBJetTags["+name+".Size]/F").c_str()	);
    root->Branch((name+".SoftMuonByIP3dBJetTags").c_str()  , &SoftMuonByIP3dBJetTags[0]  , (name+".SoftMuonByIP3dBJetTags["+name+".Size]/F").c_str()	);	
    root->Branch((name+".SoftMuonByPtBJetTags").c_str()    , &SoftMuonByPtBJetTags[0]    , (name+".SoftMuonByPtBJetTags["+name+".Size]/F").c_str()	);	
    root->Branch((name+".DoubleSVHighEffBJetTags").c_str() , &DoubleSVHighEffBJetTags[0] , (name+".DoubleSVHighEffBJetTags["+name+".Size]/F").c_str()); 



  }
  
  void Register(TTree *root, std::string name="JetInfo") {
    root->SetBranchAddress((name+".Size").c_str() , &Size);
    root->SetBranchAddress((name+".Index").c_str(), &Index[0]);
    	root->SetBranchAddress((name+".NTracks").c_str()		 , &NTracks[0]  	       );
			root->SetBranchAddress((name+".Et").c_str()			 , &Et[0]		       );
			root->SetBranchAddress((name+".Pt").c_str()			 , &Pt[0]		       );
			root->SetBranchAddress((name+".Unc").c_str()			 , &Unc[0]		       );
			root->SetBranchAddress((name+".Eta").c_str()			 , &Eta[0]		       );
			root->SetBranchAddress((name+".Phi").c_str()			 , &Phi[0]		       );
			root->SetBranchAddress((name+".JetIDLOOSE").c_str()		 , &JetIDLOOSE[0]	       ); //Add by Chiyi
			root->SetBranchAddress((name+".JetCharge").c_str()		 , &JetCharge[0]	       );
			root->SetBranchAddress((name+".NConstituents").c_str()		 , &NConstituents[0]	       );
			root->SetBranchAddress((name+".NCH").c_str()		 , &NCH[0]	       );
			root->SetBranchAddress((name+".CEF").c_str()		 , &CEF[0]  	       );
			root->SetBranchAddress((name+".NHF").c_str()		 , &NHF[0]  	       );
			root->SetBranchAddress((name+".NEF").c_str()		 , &NEF[0]  	       );
			root->SetBranchAddress((name+".CHF").c_str()		 , &CHF[0]  	       );
			root->SetBranchAddress((name+".JVAlpha").c_str()		 , &JVAlpha[0]  	       );
			root->SetBranchAddress((name+".JVBeta").c_str() 		 , &JVBeta[0]		       ); 
			root->SetBranchAddress((name+".PtCorrRaw").c_str()		 , &PtCorrRaw[0]	       );
			root->SetBranchAddress((name+".PtCorrL2").c_str()		 , &PtCorrL2[0] 	       );
			root->SetBranchAddress((name+".PtCorrL3").c_str()		 , &PtCorrL3[0] 	       );
			root->SetBranchAddress((name+".PtCorrL7g").c_str()		 , &PtCorrL7g[0]	       );
			root->SetBranchAddress((name+".PtCorrL7uds").c_str()		 , &PtCorrL7uds[0]	       );      
			root->SetBranchAddress((name+".PtCorrL7c").c_str()		 , &PtCorrL7c[0]	       );      
			root->SetBranchAddress((name+".PtCorrL7b").c_str()		 , &PtCorrL7b[0]	       );      

			root->SetBranchAddress((name+".JetBProbBJetTags").c_str()	 , &JetBProbBJetTags[0]        );
			root->SetBranchAddress((name+".JetProbBJetTags").c_str()	 , &JetProbBJetTags[0]         );
			root->SetBranchAddress((name+".TrackCountHiPurBJetTags").c_str() , &TrackCountHiPurBJetTags[0] );	   
			root->SetBranchAddress((name+".TrackCountHiEffBJetTags").c_str() , &TrackCountHiEffBJetTags[0] );
			root->SetBranchAddress((name+".SimpleSVBJetTags").c_str()	 , &SimpleSVBJetTags[0]        ); 
			root->SetBranchAddress((name+".SimpleSVHEBJetTags").c_str()	 , &SimpleSVHEBJetTags[0]      ); 
			root->SetBranchAddress((name+".SimpleSVHPBJetTags").c_str()	 , &SimpleSVHPBJetTags[0]      ); 
			root->SetBranchAddress((name+".CombinedSVBJetTags").c_str()	 , &CombinedSVBJetTags[0]      );
			root->SetBranchAddress((name+".CombinedSVMVABJetTags").c_str()   , &CombinedSVMVABJetTags[0]   );
			root->SetBranchAddress((name+".SoftElecByIP3dBJetTags").c_str()  , &SoftElecByIP3dBJetTags[0]  );
			root->SetBranchAddress((name+".SoftElecByPtBJetTags").c_str()    , &SoftElecByPtBJetTags[0]    );
			root->SetBranchAddress((name+".SoftMuonBJetTags").c_str()        , &SoftMuonBJetTags[0]	       );
			root->SetBranchAddress((name+".SoftMuonByIP3dBJetTags").c_str()  , &SoftMuonByIP3dBJetTags[0]  );
			root->SetBranchAddress((name+".SoftMuonByPtBJetTags").c_str()    , &SoftMuonByPtBJetTags[0]    ); 
			root->SetBranchAddress((name+".DoubleSVHighEffBJetTags").c_str() , &DoubleSVHighEffBJetTags[0] ); 

  }

};





#endif
    
