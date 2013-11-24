#ifndef FORMAT_H
#define FORMAT_H


#include <TTree.h>


static const int MAX_LEPTONS 	      = 256 ;
static const int MAX_TRACKS 	      = 256 ;
static const int MAX_JETS 	        = 128 ;
static const int MAX_PAIRS 	        = 512 ;
static const int MAX_PHOTONS	      = 128 ;
static const int MAX_GENS	          = 128 ;
static const int MAX_VERTICES       = 256 ;
static const int MAX_BX		          = 128 ;
static const int N_TRIGGER_BOOKINGS = 5842;  

class EvtInfoBranches {
 public:
  int      RunNo;
  long int EvtNo;
  int      BxNo;
  int      LumiNo;
  int      Orbit;
  int      McFlag;	        // MC or not MC, that's the question
  int      McSigTag;        // MC Signature tag	  - 0: others, 1: 2L (opposite-sign), 2: 2L (same-sign), 3: 3L, 4: 4L
  int      McbprimeMode[2]; // b'(bar) decay mode   - 0: others, 1: tW, 2: cW, 3: bZ, 4: bH
  int      MctprimeMode[2]; // t'(bar) decay mode   - 0: others, 1: bW, 2: tZ, 3: tH, 4: tgamma
  int      McWMode[4];      // W- from b'(t'bar)/W+ from t/W+ from b'bar(t')/W- from tbar - 0: others, 1: enu, 2: munu, 3: taunu, 4: jj
  int      McZMode[2];      // Z from b'(bar)/t'(t'bar)	 - 0: others, 1: ee, 2: mumu, 3: tautau, 4: nunu, 5: bb, 6: jj
  float    McbprimeMass[2]; // mass: b'(bar)  
  float    MctprimeMass[2]; // mass: t'(bar)  
  float    MctopMass[2];    // mass: top(bar)
  float    McWMass[4];      // mass: W- from b'(t'bar)/W+ from t/W+ from b'bar(t')/W- from tbar 
  float    McZMass[2];      // mass: Z from b'(bar) or t'(bar)
  float    McDauPt[14];     // Generic MC daughter information
  float    McDauEta[14];    // MC daughters: 0-1: hard jet from b'bar/t'bar, 2-9: W daughters, 10-13: Z daughters 
  float    McDauPhi[14];    // 
  float    RhoPU[2];        // [electron,muon]	 
  float    SigmaPU[2];      // [electron,muon] 
  int      McDauPdgID[14];	 
  int      PDFid1;
  int      PDFid2;
  float    PDFx1;
  float    PDFx2;
  float    PDFscale;
  float    PDFv1;
  float    PDFv2;
  float    qScale ; 
  float    alphaQCD ; 
  float    alphaQED ; 
  float    Weight;
  float    WeightPU; 

  //BeamSpot
  float BeamSpotX;
  float BeamSpotY;
  float BeamSpotZ;

  // PU
  int   nBX;
  int   nPU[MAX_BX];
  int   BXPU[MAX_BX];
  float TrueIT[MAX_BX];

  float PFMET;
  float PFMETPhi;
  float PFRawMET;
  float PFRawMETPhi;
  float PFSumEt;
  float PFMETSig;
  float PFGenMET;
  float PFGenMETPhi;

  float PFMETx; //Uly 2011-04-04
  float PFMETy; //Uly 2011-04-04

  int   TrgCount;	// No. of fired booking bits
  int   nTrgBook;
  char  TrgBook[N_TRIGGER_BOOKINGS];	// Trigger bits
  int   nHLT;
  bool  HLTbits[N_TRIGGER_BOOKINGS];
  int   L1[128]; // L1 trigger bits
  int   TT[64]; 	// Techical trigger bits
  float HighPurityFraction; //Added by Dmitry to help filter out bad events
  int   NofTracks; //Added by Dmitry to help filter out bad events
  
  void RegisterTree(TTree *root) {
    root->Branch("EvtInfo.RunNo"	    , &RunNo	       , "EvtInfo.RunNo/I"	    );
    root->Branch("EvtInfo.EvtNo"	    , &EvtNo	       , "EvtInfo.EvtNo/L"	    );
    root->Branch("EvtInfo.BxNo"	      , &BxNo	       , "EvtInfo.BxNo/I"	        );
    root->Branch("EvtInfo.LumiNo"	    , &LumiNo	       , "EvtInfo.LumiNo/I"	    );
    root->Branch("EvtInfo.Orbit"	    , &Orbit	       , "EvtInfo.Orbit/I"	    );
    root->Branch("EvtInfo.McFlag"	    , &McFlag	       , "EvtInfo.McFlag/I"	    );
    root->Branch("EvtInfo.McSigTag"	    , &McSigTag        , "EvtInfo.McSigTag/I"	    );
    root->Branch("EvtInfo.McbprimeMode" , &McbprimeMode[0] , "EvtInfo.McbprimeMode[2]/I");
    root->Branch("EvtInfo.MctprimeMode" , &MctprimeMode[0] , "EvtInfo.MctprimeMode[2]/I");
    root->Branch("EvtInfo.McWMode"      , &McWMode[0]      , "EvtInfo.McWMode[4]/I"     );
    root->Branch("EvtInfo.McZMode"      , &McZMode[0]      , "EvtInfo.McZMode[2]/I"     );
    root->Branch("EvtInfo.McbprimeMass" , &McbprimeMass[0] , "EvtInfo.McbprimeMass[2]/F");
    root->Branch("EvtInfo.MctprimeMass" , &MctprimeMass[0] , "EvtInfo.MctprimeMass[2]/F");
    root->Branch("EvtInfo.MctopMass"    , &MctopMass[0]    , "EvtInfo.MctopMass[2]/F"   );
    root->Branch("EvtInfo.McWMass"      , &McWMass[0]      , "EvtInfo.McWMass[4]/F"     );
    root->Branch("EvtInfo.McZMass"      , &McZMass[0]      , "EvtInfo.McZMass[2]/F"     );
    root->Branch("EvtInfo.McDauPt"      , &McDauPt[0]      , "EvtInfo.McDauPt[14]/F"    );
    root->Branch("EvtInfo.McDauEta"     , &McDauEta[0]     , "EvtInfo.McDauEta[14]/F"   );
    root->Branch("EvtInfo.McDauPhi"     , &McDauPhi[0]     , "EvtInfo.McDauPhi[14]/F"   );  	 
    root->Branch("EvtInfo.McDauPdgID"   , &McDauPdgID[0]   , "EvtInfo.McDauPdgID[14]/I" );  		 
    root->Branch("EvtInfo.PDFid1"	    , &PDFid1	       , "EvtInfo.PDFid1/I"	    );
    root->Branch("EvtInfo.PDFid2"	    , &PDFid2	       , "EvtInfo.PDFid2/I"	    );
    root->Branch("EvtInfo.PDFx1"	    , &PDFx1	       , "EvtInfo.PDFx1/F"	    );
    root->Branch("EvtInfo.RhoPU"	    , &RhoPU[0]	     , "EvtInfo.RhoPU[2]/F"	  );
    root->Branch("EvtInfo.SigmaPU"	  , &SigmaPU[0]    , "EvtInfo.SigmaPU[2]/F"	);
    root->Branch("EvtInfo.PDFx2"	    , &PDFx2	       , "EvtInfo.PDFx2/F"	    );
    root->Branch("EvtInfo.PDFscale"	  , &PDFscale	     , "EvtInfo.PDFscale/F"	  );
    root->Branch("EvtInfo.PDFv1"	    , &PDFv1	       , "EvtInfo.PDFv1/F"	    );
    root->Branch("EvtInfo.PDFv2"	    , &PDFv2	       , "EvtInfo.PDFv2/F"	    );		 
    root->Branch("EvtInfo.qScale"     , &qScale  	     , "EvtInfo.qScale/F"			);
    root->Branch("EvtInfo.alphaQCD"   , &alphaQCD	     , "EvtInfo.alphaQCD/F"		);
    root->Branch("EvtInfo.alphaQED"   , &alphaQED	     , "EvtInfo.alphaQED/F"		);
    root->Branch("EvtInfo.Weight"     , &Weight	       , "EvtInfo.Weight/F"			);
    root->Branch("EvtInfo.WeightPU"   , &WeightPU      , "EvtInfo.Weight/F"			);

    root->Branch("EvtInfo.BeamSpotX"    , &BeamSpotX       , "EvtInfo.BeamSpotX/F"      );
    root->Branch("EvtInfo.BeamSpotY"    , &BeamSpotY       , "EvtInfo.BeamSpotY/F"      );
    root->Branch("EvtInfo.BeamSpotZ"    , &BeamSpotZ       , "EvtInfo.BeamSpotZ/F"      );
    root->Branch("EvtInfo.PFMET"        , &PFMET           , "EvtInfo.PFMET/F"          );
    root->Branch("EvtInfo.PFMETPhi"     , &PFMETPhi        , "EvtInfo.PFMETPhi/F"       );
    root->Branch("EvtInfo.PFRawMET"     , &PFRawMET        , "EvtInfo.PFRawMET/F"       );
    root->Branch("EvtInfo.PFRawMETPhi"  , &PFRawMETPhi     , "EvtInfo.PFRawMETPhi/F"    );
    root->Branch("EvtInfo.PFSumEt"      , &PFSumEt         , "EvtInfo.PFSumEt/F"        );
    root->Branch("EvtInfo.PFMETSig"     , &PFMETSig        , "EvtInfo.PFMETSig/F"       );
    root->Branch("EvtInfo.PFGenMET"     , &PFGenMET        , "EvtInfo.PFGenMET/F"       );
    root->Branch("EvtInfo.PFGenMETPhi"  , &PFGenMETPhi     , "EvtInfo.PFGenMETPhi/F"    );
    root->Branch("EvtInfo.PFMETx"       , &PFMETx          , "EvtInfo.PFMETx/F"         ); //Uly 2011-04-04
    root->Branch("EvtInfo.PFMETy"       , &PFMETy          , "EvtInfo.PFMETy/F"         ); //Uly 2011-04-04
    root->Branch("EvtInfo.TrgCount"     , &TrgCount        , "EvtInfo.TrgCount/I"	    );
    root->Branch("EvtInfo.nTrgBook"     , &nTrgBook	       , "EvtInfo.nTrgBook/I"       );
    root->Branch("EvtInfo.TrgBook"      , &TrgBook[0]      , "EvtInfo.TrgBook[EvtInfo.nTrgBook]/B"    );
    root->Branch("EvtInfo.L1"	    , &L1[0]	       , "EvtInfo.L1[128]/I"	    );
    root->Branch("EvtInfo.TT"	    , &TT[0]	       , "EvtInfo.TT[64]/I"	    );
    root->Branch("EvtInfo.HighPurityFraction"    , &HighPurityFraction  , "EvtInfo.HighPurityFraction/F"   );
    root->Branch("EvtInfo.NofTracks"             , &NofTracks           , "EvtInfo.NofTracks/I"            );
    root->Branch("EvtInfo.nHLT"		     , &nHLT	 	    , "EvtInfo.nHLT/I"                 );
    root->Branch("EvtInfo.HLTbits"		     , HLTbits		    , "EvtInfo.HLTbits[EvtInfo.nHLT]/O");
    root->Branch("EvtInfo.nBX"		     , &nBX	 	    , "EvtInfo.nBX/I"                  );
    root->Branch("EvtInfo.nPU"		     , &nPU[0]		    , "EvtInfo.nPU[EvtInfo.nBX]/I"     );
    root->Branch("EvtInfo.BXPU"		     , &BXPU[0]		    , "EvtInfo.BXPU[EvtInfo.nBX]/I"    );
    root->Branch("EvtInfo.TrueIT"		     , &TrueIT[0]	    , "EvtInfo.TrueIT[EvtInfo.nBX]/F"  ); 
  }										    

  void Register(TTree *root) {
    root->SetBranchAddress("EvtInfo.RunNo"        , &RunNo  	 );
    root->SetBranchAddress("EvtInfo.EvtNo"        , &EvtNo  	 );
    root->SetBranchAddress("EvtInfo.BxNo"         , &BxNo   	 );
    root->SetBranchAddress("EvtInfo.LumiNo"       , &LumiNo  	 );
    root->SetBranchAddress("EvtInfo.Orbit"        , &Orbit  	 );
    root->SetBranchAddress("EvtInfo.McFlag"       , &McFlag 	 );
    root->SetBranchAddress("EvtInfo.McSigTag"     , &McSigTag	 );
    root->SetBranchAddress("EvtInfo.McbprimeMode" , &McbprimeMode[0] );
    root->SetBranchAddress("EvtInfo.MctprimeMode" , &MctprimeMode[0] );
    root->SetBranchAddress("EvtInfo.McWMode"      , &McWMode[0]	 );
    root->SetBranchAddress("EvtInfo.McZMode"      , &McZMode[0]	 );
    root->SetBranchAddress("EvtInfo.McbprimeMass" , &McbprimeMass[0] );
    root->SetBranchAddress("EvtInfo.MctprimeMass" , &MctprimeMass[0] );
    root->SetBranchAddress("EvtInfo.MctopMass"    , &MctopMass[0]    );
    root->SetBranchAddress("EvtInfo.McWMass"      , &McWMass[0]	 );
    root->SetBranchAddress("EvtInfo.McZMass"      , &McZMass[0]	 );
    root->SetBranchAddress("EvtInfo.McDauPt"      , &McDauPt[0]	 );
    root->SetBranchAddress("EvtInfo.McDauEta"     , &McDauEta[0]	 );
    root->SetBranchAddress("EvtInfo.McDauPhi"     , &McDauPhi[0]	 );	      
    root->SetBranchAddress("EvtInfo.McDauPdgID"   , &McDauPdgID[0]   );		      
    root->SetBranchAddress("EvtInfo.PDFid1"       , &PDFid1 	 );
    root->SetBranchAddress("EvtInfo.PDFid2"       , &PDFid2 	 );
    root->SetBranchAddress("EvtInfo.PDFx1"        , &PDFx1  	 );
    root->SetBranchAddress("EvtInfo.RhoPU"        , &RhoPU[0]  );
    root->SetBranchAddress("EvtInfo.SigmaPU"      , &SigmaPU[0]);
    root->SetBranchAddress("EvtInfo.PDFx2"        , &PDFx2  	 );
    root->SetBranchAddress("EvtInfo.PDFscale"     , &PDFscale	 );
    root->SetBranchAddress("EvtInfo.PDFv1"        , &PDFv1  	 );
    root->SetBranchAddress("EvtInfo.PDFv2"        , &PDFv2  	 );	      
    root->SetBranchAddress("EvtInfo.qScale"       , &qScale  	 );
    root->SetBranchAddress("EvtInfo.alphaQCD"     , &alphaQCD	 );
    root->SetBranchAddress("EvtInfo.alphaQED"     , &alphaQED	 );
    root->SetBranchAddress("EvtInfo.Weight"       , &Weight		 ); 
    root->SetBranchAddress("EvtInfo.WeightPU"     , &WeightPU  ); 

    root->SetBranchAddress("EvtInfo.BeamSpotX"    , &BeamSpotX );
    root->SetBranchAddress("EvtInfo.BeamSpotY"    , &BeamSpotY );
    root->SetBranchAddress("EvtInfo.BeamSpotZ"    , &BeamSpotZ );

    root->SetBranchAddress("EvtInfo.PFMET"        , &PFMET           );
    root->SetBranchAddress("EvtInfo.PFMETPhi"     , &PFMETPhi        );
    root->SetBranchAddress("EvtInfo.PFRawMET"     , &PFRawMET        );
    root->SetBranchAddress("EvtInfo.PFRawMETPhi"  , &PFRawMETPhi     );
    root->SetBranchAddress("EvtInfo.PFSumEt"      , &PFSumEt         );
    root->SetBranchAddress("EvtInfo.PFMETSig"     , &PFMETSig        );
    root->SetBranchAddress("EvtInfo.PFGenMET"     , &PFGenMET        );
    root->SetBranchAddress("EvtInfo.PFGenMETPhi"  , &PFGenMETPhi     );
    root->SetBranchAddress("EvtInfo.PFMETx"       , &PFMETx          ); //Uly 2011-04-04
    root->SetBranchAddress("EvtInfo.PFMETy"       , &PFMETy          ); //Uly 2011-04-04

    root->SetBranchAddress("EvtInfo.TrgCount"     , &TrgCount	 );
    root->SetBranchAddress("EvtInfo.nTrgBook"     , &nTrgBook	 );
    root->SetBranchAddress("EvtInfo.TrgBook"      , TrgBook	         );
    root->SetBranchAddress("EvtInfo.L1"           , &L1[0]  	 );
    root->SetBranchAddress("EvtInfo.TT"           , &TT[0]  	 );
    root->SetBranchAddress("EvtInfo.HighPurityFraction"    , &HighPurityFraction	 );
    root->SetBranchAddress("EvtInfo.NofTracks"    , &NofTracks	 );
    root->SetBranchAddress("EvtInfo.nHLT"         , &nHLT	         );
    root->SetBranchAddress("EvtInfo.HLTbits"      , HLTbits          );
    root->SetBranchAddress("EvtInfo.nBX"	      , &nBX             );
    root->SetBranchAddress("EvtInfo.nPU"	      , &nPU[0]          );
    root->SetBranchAddress("EvtInfo.BXPU"	      , &BXPU[0]         );
    root->SetBranchAddress("EvtInfo.TrueIT"	      , &TrueIT[0]       );
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

  void RegisterTree(TTree *root) {
    root->Branch("VertexInfo.Size"     , &Size	          , "VertexInfo.Size/I"	    );
    root->Branch("VertexInfo.isValid"  , &isValid[0]      , "VertexInfo.isValid[VertexInfo.Size]/I"	    );
    root->Branch("VertexInfo.isFake"   , &isFake[0]       , "VertexInfo.isFake[VertexInfo.Size]/O"	    ); //Uly 2011-04-04
    root->Branch("VertexInfo.Type"	    , &Type[0]	       , "VertexInfo.Type[VertexInfo.Size]/I"	    );
    root->Branch("VertexInfo.Ndof"	    , &Ndof[0]	       , "VertexInfo.Ndof[VertexInfo.Size]/F"	    );
    root->Branch("VertexInfo.NormalizedChi2"	    , &NormalizedChi2[0]	       , "VertexInfo.NormalizedChi2[VertexInfo.Size]/F"	    );
    root->Branch("VertexInfo.Pt_Sum"	    , &Pt_Sum[0]	       , "VertexInfo.Pt_Sum[VertexInfo.Size]/F"	    );
    root->Branch("VertexInfo.Pt_Sum2"	    , &Pt_Sum2[0]	       , "VertexInfo.Pt_Sum2[VertexInfo.Size]/F"	    );
    root->Branch("VertexInfo.x"	    , &x[0]	       , "VertexInfo.x[VertexInfo.Size]/F"	    );
    root->Branch("VertexInfo.y"	    , &y[0]	       , "VertexInfo.y[VertexInfo.Size]/F"	    );
    root->Branch("VertexInfo.z"	    , &z[0]	       , "VertexInfo.z[VertexInfo.Size]/F"	    );
    root->Branch("VertexInfo.Rho"	    , &Rho[0]	       , "VertexInfo.Rho[VertexInfo.Size]/F"	    );
  }										    
  void Register(TTree *root) {
    root->SetBranchAddress("VertexInfo.Size"        , &Size  	 );
    root->SetBranchAddress("VertexInfo.isValid"     , &isValid[0]  	 );
    root->SetBranchAddress("VertexInfo.isFake"      , &isFake[0]  	 ); //Uly 2011-04-04
    root->SetBranchAddress("VertexInfo.Type"        , &Type[0]  	 );
    root->SetBranchAddress("VertexInfo.Ndof"        , &Ndof[0]  	 );
    root->SetBranchAddress("VertexInfo.NormalizedChi2"        , &NormalizedChi2[0]  	 );
    root->SetBranchAddress("VertexInfo.Pt_Sum"        , &Pt_Sum[0]  	 );
    root->SetBranchAddress("VertexInfo.Pt_Sum2"        , &Pt_Sum2[0]  	 );
    root->SetBranchAddress("VertexInfo.x"        , &x[0]  	 );
    root->SetBranchAddress("VertexInfo.y"        , &y[0]  	 );
    root->SetBranchAddress("VertexInfo.z"        , &z[0]  	 );
    root->SetBranchAddress("VertexInfo.Rho"        , &Rho[0]  	 );
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


class GenInfoBranches {
 public:
  int Size;

  float Pt[MAX_GENS];
  float Eta[MAX_GENS];
  float Phi[MAX_GENS];
  float Mass[MAX_GENS];
  int PdgID[MAX_GENS];
  int Status[MAX_GENS];

  int nMo[MAX_GENS];
  int nDa[MAX_GENS];

  float Mo0Pt[MAX_GENS];
  float Mo0Eta[MAX_GENS];
  float Mo0Phi[MAX_GENS];
  float Mo0Mass[MAX_GENS];
  int   Mo0PdgID[MAX_GENS];
  int   Mo0Status[MAX_GENS];

  float Mo1Pt[MAX_GENS];
  float Mo1Eta[MAX_GENS];
  float Mo1Phi[MAX_GENS];
  float Mo1Mass[MAX_GENS];
  int   Mo1PdgID[MAX_GENS];
  int   Mo1Status[MAX_GENS];

  float Da0Pt[MAX_GENS];
  float Da0Eta[MAX_GENS];
  float Da0Phi[MAX_GENS];
  float Da0Mass[MAX_GENS];
  int   Da0PdgID[MAX_GENS];
  int   Da0Status[MAX_GENS];

  float Da1Pt[MAX_GENS];
  float Da1Eta[MAX_GENS];
  float Da1Phi[MAX_GENS];
  float Da1Mass[MAX_GENS];
  int   Da1PdgID[MAX_GENS];
  int   Da1Status[MAX_GENS];

  int Mo1[MAX_GENS];
  int Mo2[MAX_GENS];
  int Da1[MAX_GENS];
  int Da2[MAX_GENS];

  int   ncQuarks               ;
  float cQuark_pT     [MAX_GENS];
  float cQuark_eta    [MAX_GENS];
  float cQuark_phi    [MAX_GENS];
  int   cQuark_pdgID  [MAX_GENS];
  int   cQuark_status [MAX_GENS];
  int   cQuark_fromGSP[MAX_GENS];

  int   nbQuarks               ;
  float bQuark_pT     [MAX_GENS];
  float bQuark_eta    [MAX_GENS];
  float bQuark_phi    [MAX_GENS];
  int   bQuark_pdgID  [MAX_GENS];
  int   bQuark_status [MAX_GENS];
  int   bQuark_fromGSP[MAX_GENS];

  int   nBHadrons                    ;
  float BHadron_pT          [MAX_GENS];
  float BHadron_eta         [MAX_GENS];
  float BHadron_phi         [MAX_GENS];
  float BHadron_mass        [MAX_GENS];
  int   BHadron_pdgID       [MAX_GENS];
  int   BHadron_status      [MAX_GENS];
  int   BHadron_mother      [MAX_GENS];
  int   BHadron_hasBdaughter[MAX_GENS];

  int   nDHadrons                    ;
  float DHadron_pT          [MAX_GENS];
  float DHadron_eta         [MAX_GENS];
  float DHadron_phi         [MAX_GENS];
  float DHadron_mass        [MAX_GENS];
  int   DHadron_pdgID       [MAX_GENS];
  int   DHadron_status      [MAX_GENS];
  int   DHadron_mother      [MAX_GENS];
  int   DHadron_hasBdaughter[MAX_GENS];

  void RegisterTree(TTree *root) {
    root->Branch("GenInfo.Size"	, &Size		, "GenInfo.Size/I"			);
    root->Branch("GenInfo.Pt"	, &Pt[0]	, "GenInfo.Pt[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Eta"	, &Eta[0]	, "GenInfo.Eta[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Phi"	, &Phi[0]	, "GenInfo.Phi[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Mass"	, &Mass[0]	, "GenInfo.Mass[GenInfo.Size]/F"	);
    root->Branch("GenInfo.PdgID" , &PdgID[0]	, "GenInfo.PdgID[GenInfo.Size]/I"	);
    root->Branch("GenInfo.Status", &Status[0]	, "GenInfo.Status[GenInfo.Size]/I"	);
    root->Branch("GenInfo.Mo0Pt"	  , &Mo0Pt[0]	    , "GenInfo.Mo0Pt[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Mo0Eta"	  , &Mo0Eta[0]	  , "GenInfo.Mo0Eta[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Mo0Phi"	  , &Mo0Phi[0]	  , "GenInfo.Mo0Phi[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Mo0Mass"	, &Mo0Mass[0]	  , "GenInfo.Mo0Mass[GenInfo.Size]/F"	);
    root->Branch("GenInfo.Mo0PdgID" , &Mo0PdgID[0]	, "GenInfo.Mo0PdgID[GenInfo.Size]/I"	);
    root->Branch("GenInfo.Mo0Status", &Mo0Status[0]	, "GenInfo.Mo0Status[GenInfo.Size]/I"	);
    root->Branch("GenInfo.Mo1Pt"	  , &Mo1Pt[0]	    , "GenInfo.Mo1Pt[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Mo1Eta"	  , &Mo1Eta[0]	  , "GenInfo.Mo1Eta[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Mo1Phi"	  , &Mo1Phi[0]	  , "GenInfo.Mo1Phi[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Mo1Mass"	, &Mo1Mass[0]	  , "GenInfo.Mo1Mass[GenInfo.Size]/F"	);
    root->Branch("GenInfo.Mo1PdgID" , &Mo1PdgID[0]	, "GenInfo.Mo1PdgID[GenInfo.Size]/I"	);
    root->Branch("GenInfo.Mo1Status", &Mo1Status[0]	, "GenInfo.Mo1Status[GenInfo.Size]/I"	);
    root->Branch("GenInfo.Da0Pt"	  , &Da0Pt[0]	    , "GenInfo.Da0Pt[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Da0Eta"	  , &Da0Eta[0]	  , "GenInfo.Da0Eta[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Da0Phi"	  , &Da0Phi[0]	  , "GenInfo.Da0Phi[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Da0Mass"	, &Da0Mass[0]	  , "GenInfo.Da0Mass[GenInfo.Size]/F"	);
    root->Branch("GenInfo.Da0PdgID" , &Da0PdgID[0]	, "GenInfo.Da0PdgID[GenInfo.Size]/I"	);
    root->Branch("GenInfo.Da0Status", &Da0Status[0]	, "GenInfo.Da0Status[GenInfo.Size]/I"	);
    root->Branch("GenInfo.Da1Pt"	  , &Da1Pt[0]	    , "GenInfo.Da1Pt[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Da1Eta"	  , &Da1Eta[0]	  , "GenInfo.Da1Eta[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Da1Phi"	  , &Da1Phi[0]	  , "GenInfo.Da1Phi[GenInfo.Size]/F"		);
    root->Branch("GenInfo.Da1Mass"	, &Da1Mass[0]	  , "GenInfo.Da1Mass[GenInfo.Size]/F"	);
    root->Branch("GenInfo.Da1PdgID" , &Da1PdgID[0]	, "GenInfo.Da1PdgID[GenInfo.Size]/I"	);
    root->Branch("GenInfo.Da1Status", &Da1Status[0]	, "GenInfo.Da1Status[GenInfo.Size]/I"	);
    root->Branch("GenInfo.nMo"	, &nMo[0]	, "GenInfo.nMo[GenInfo.Size]/I"		);
    root->Branch("GenInfo.nDa"	, &nDa[0]	, "GenInfo.nDa[GenInfo.Size]/I"		);
    root->Branch("GenInfo.Mo1"	, &Mo1[0]	, "GenInfo.Mo1[GenInfo.Size]/I"		);
    root->Branch("GenInfo.Mo2"	, &Mo2[0]	, "GenInfo.Mo2[GenInfo.Size]/I"		);
    root->Branch("GenInfo.Da1"	, &Da1[0]	, "GenInfo.Da1[GenInfo.Size]/I"		);
    root->Branch("GenInfo.Da2"	, &Da2[0]	, "GenInfo.Da2[GenInfo.Size]/I"		);
    root->Branch("GenInfo.ncQuarks"      	, &ncQuarks         	, "GenInfo.ncQuarks/I"		                  );
    root->Branch("GenInfo.cQuark_pT"     	, &cQuark_pT     [0]	, "GenInfo.cQuark_pT[GenInfo.Size]/F"		    );
    root->Branch("GenInfo.cQuark_eta"    	, &cQuark_eta    [0]	, "GenInfo.cQuark_eta[GenInfo.Size]/F"		  );
    root->Branch("GenInfo.cQuark_phi"    	, &cQuark_phi    [0]	, "GenInfo.cQuark_phi[GenInfo.Size]/F"		  );
    root->Branch("GenInfo.cQuark_pdgID"  	, &cQuark_pdgID  [0]	, "GenInfo.cQuark_pdgID[GenInfo.Size]/I"		);
    root->Branch("GenInfo.cQuark_status" 	, &cQuark_status [0]	, "GenInfo.cQuark_status[GenInfo.Size]/I"		);
    root->Branch("GenInfo.cQuark_fromGSP"	, &cQuark_fromGSP[0]	, "GenInfo.cQuark_fromGSP[GenInfo.Size]/I"	);
  }

  void Register(TTree *root) {
    root->SetBranchAddress("GenInfo.Size"       , &Size		);
    root->SetBranchAddress("GenInfo.Pt"	        , &Pt[0]	);
    root->SetBranchAddress("GenInfo.Eta"	, &Eta[0]	);
    root->SetBranchAddress("GenInfo.Phi"	, &Phi[0]	);
    root->SetBranchAddress("GenInfo.Mass"	, &Mass[0]	);
    root->SetBranchAddress("GenInfo.PdgID"	, &PdgID[0]	);
    root->SetBranchAddress("GenInfo.Status"	, &Status[0]	);
    root->SetBranchAddress("GenInfo.Mo0Pt"	  , &Mo0Pt[0]	    ) ; 
    root->SetBranchAddress("GenInfo.Mo0Eta"	  , &Mo0Eta[0]	  ) ; 
    root->SetBranchAddress("GenInfo.Mo0Phi"	  , &Mo0Phi[0]	  ) ; 
    root->SetBranchAddress("GenInfo.Mo0Mass"	, &Mo0Mass[0]	  ) ; 
    root->SetBranchAddress("GenInfo.Mo0PdgID" , &Mo0PdgID[0]	) ; 
    root->SetBranchAddress("GenInfo.Mo0Status", &Mo0Status[0]	) ; 
    root->SetBranchAddress("GenInfo.Mo1Pt"	  , &Mo1Pt[0]	    ) ; 
    root->SetBranchAddress("GenInfo.Mo1Eta"	  , &Mo1Eta[0]	  ) ; 
    root->SetBranchAddress("GenInfo.Mo1Phi"	  , &Mo1Phi[0]	  ) ; 
    root->SetBranchAddress("GenInfo.Mo1Mass"	, &Mo1Mass[0]	  ) ; 
    root->SetBranchAddress("GenInfo.Mo1PdgID" , &Mo1PdgID[0]	) ; 
    root->SetBranchAddress("GenInfo.Mo1Status", &Mo1Status[0]	) ; 
    root->SetBranchAddress("GenInfo.nMo"	, &nMo[0]	);
    root->SetBranchAddress("GenInfo.nDa"	, &nDa[0]	);
    root->SetBranchAddress("GenInfo.Mo1"	, &Mo1[0]	);
    root->SetBranchAddress("GenInfo.Mo2"	, &Mo2[0]	);
    root->SetBranchAddress("GenInfo.Da1"	, &Da1[0]	);
    root->SetBranchAddress("GenInfo.Da2"	, &Da2[0]	);
    root->SetBranchAddress("GenInfo.Da0Pt"    , &Da0Pt[0]       );
    root->SetBranchAddress("GenInfo.Da0Eta"   , &Da0Eta[0]    );
    root->SetBranchAddress("GenInfo.Da0Phi"   , &Da0Phi[0]    );
    root->SetBranchAddress("GenInfo.Da0Mass"  , &Da0Mass[0]     );
    root->SetBranchAddress("GenInfo.Da0PdgID" , &Da0PdgID[0]       );
    root->SetBranchAddress("GenInfo.Da0Status", &Da0Status[0]      );
    root->SetBranchAddress("GenInfo.Da1Pt"    , &Da1Pt[0]       );
    root->SetBranchAddress("GenInfo.Da1Eta"   , &Da1Eta[0]    );
    root->SetBranchAddress("GenInfo.Da1Phi"   , &Da1Phi[0]    );
    root->SetBranchAddress("GenInfo.Da1Mass"  , &Da1Mass[0]     );
    root->SetBranchAddress("GenInfo.Da1PdgID" , &Da1PdgID[0]       );
    root->SetBranchAddress("GenInfo.Da1Status", &Da1Status[0]      );
    root->SetBranchAddress("GenInfo.ncQuarks"               , &ncQuarks             );
    root->SetBranchAddress("GenInfo.cQuark_pT"              , &cQuark_pT     [0]    );
    root->SetBranchAddress("GenInfo.cQuark_eta"             , &cQuark_eta    [0]    );
    root->SetBranchAddress("GenInfo.cQuark_phi"             , &cQuark_phi    [0]    );
    root->SetBranchAddress("GenInfo.cQuark_pdgID"           , &cQuark_pdgID  [0]    );
    root->SetBranchAddress("GenInfo.cQuark_status"          , &cQuark_status [0]    );
    root->SetBranchAddress("GenInfo.cQuark_fromGSP"         , &cQuark_fromGSP[0]    );
  }

};


class JetInfoBranches {

  public:
    int   Size; 
    int   Index[MAX_JETS];
    int   ParentIndex[MAX_JETS]; // used for subjets

    float Et[MAX_JETS];
    float Pt[MAX_JETS];
    float Unc[MAX_JETS];
    float Eta[MAX_JETS];
    float Phi[MAX_JETS];

    float Energy[MAX_JETS]; 
    float Px[MAX_JETS];     
    float Py[MAX_JETS];     
    float Pz[MAX_JETS];     
    float Mass[MAX_JETS];
    float Area[MAX_JETS];

    float EtPruned[MAX_JETS]; //// To add
    float PtPruned[MAX_JETS]; //// To add
    float UncPruned[MAX_JETS]; //// To add
    float EtaPruned[MAX_JETS]; //// To add
    float PhiPruned[MAX_JETS]; //// To add 

    float EnergyPruned[MAX_JETS]; //// To add 
    float PxPruned[MAX_JETS]; //// To add     
    float PyPruned[MAX_JETS]; //// To add     
    float PzPruned[MAX_JETS]; //// To add     
    float MassPruned[MAX_JETS]; //// To add
    float AreaPruned[MAX_JETS]; //// To add

    float tau1[MAX_JETS]; //// To add 
    float tau2[MAX_JETS]; //// To add 
    float tau3[MAX_JETS]; //// To add 

    float GenJetPt[MAX_JETS];
    float GenJetEta[MAX_JETS];
    float GenJetPhi[MAX_JETS];
    float GenPt[MAX_JETS];
    float GenEta[MAX_JETS];
    float GenPhi[MAX_JETS];
    int   GenPdgID[MAX_JETS];
    int   GenFlavor[MAX_JETS];
    int   GenMCTag[MAX_JETS]; 

    bool  JetIDLOOSE[MAX_JETS]; 
    bool  JetIDTIGHT[MAX_JETS]; //// To add 
    float JetCharge[MAX_JETS];
    float QGTagsMLP[MAX_JETS];
    float QGTagsLikelihood[MAX_JETS];
    int   NConstituents[MAX_JETS];
    int   NTracks[MAX_JETS];
    int   NCH[MAX_JETS];
    float CEF[MAX_JETS];
    float NHF[MAX_JETS];
    float NEF[MAX_JETS];
    float CHF[MAX_JETS];

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
    float CombinedSVBJetTags[MAX_JETS];
    float CombinedSVMVABJetTags[MAX_JETS];
    float SoftElecByIP3dBJetTags[MAX_JETS];
    float SoftElecByPtBJetTags[MAX_JETS];  
    float SoftMuonBJetTags[MAX_JETS];      
    float SoftMuonByIP3dBJetTags[MAX_JETS];
    float SoftMuonByPtBJetTags[MAX_JETS];  
    float DoubleSVHighEffBJetTags[MAX_JETS]; 

    int Jet_FatJetIdx[MAX_JETS];
    int Jet_SubJet1Idx[MAX_JETS];
    int Jet_SubJet2Idx[MAX_JETS];

    void RegisterTree(TTree *root, std::string name="JetInfo") {
      root->Branch((name+".Size").c_str(), &Size, (name+".Size/I").c_str() );
      root->Branch((name+".Index").c_str()               , &Index[0]                   , (name+".Index["+name+".Size]/I").c_str());
      root->Branch((name+".NTracks").c_str()  	       , &NTracks[0]		     , (name+".NTracks["+name+".Size]/I").c_str()		);
      root->Branch((name+".Et").c_str()		       , &Et[0] 		     , (name+".Et["+name+".Size]/F").c_str()			);
      root->Branch((name+".Pt").c_str()		       , &Pt[0] 		     , (name+".Pt["+name+".Size]/F").c_str()			);
      root->Branch((name+".Unc").c_str()		       , &Unc[0] 		     , (name+".Unc["+name+".Size]/F").c_str()			);
      root->Branch((name+".Eta").c_str()		       , &Eta[0]		     , (name+".Eta["+name+".Size]/F").c_str()			);
      root->Branch((name+".Phi").c_str()		       , &Phi[0]		     , (name+".Phi["+name+".Size]/F").c_str()			);
      root->Branch((name+".Energy").c_str()	       , &Energy[0] 		     , (name+".Energy["+name+".Size]/F").c_str()		); 
      root->Branch((name+".Px").c_str()		       , &Px[0] 		     , (name+".Px["+name+".Size]/F").c_str()			); 
      root->Branch((name+".Py").c_str()		       , &Py[0] 		     , (name+".Py["+name+".Size]/F").c_str()			); 
      root->Branch((name+".Pz").c_str()		       , &Pz[0] 		     , (name+".Pz["+name+".Size]/F").c_str()			); 
      root->Branch((name+".Mass").c_str()	       , &Mass[0] 		     , (name+".Mass["+name+".Size]/F").c_str()			); 
      root->Branch((name+".Area").c_str()	       , &Area[0] 		     , (name+".Area["+name+".Size]/F").c_str()			); 
      root->Branch((name+".EtPruned").c_str()		       , &EtPruned[0] 		     , (name+".EtPruned["+name+".Size]/F").c_str()			);
      root->Branch((name+".PtPruned").c_str()		       , &PtPruned[0] 		     , (name+".PtPruned["+name+".Size]/F").c_str()			);
      root->Branch((name+".UncPruned").c_str()		       , &UncPruned[0] 		     , (name+".UncPruned["+name+".Size]/F").c_str()			);
      root->Branch((name+".EtaPruned").c_str()		       , &EtaPruned[0]		     , (name+".EtaPruned["+name+".Size]/F").c_str()			);
      root->Branch((name+".PhiPruned").c_str()		       , &PhiPruned[0]		     , (name+".PhiPruned["+name+".Size]/F").c_str()			);
      root->Branch((name+".EnergyPruned").c_str()	       , &EnergyPruned[0] 		     , (name+".EnergyPruned["+name+".Size]/F").c_str()		); 
      root->Branch((name+".PxPruned").c_str()		       , &PxPruned[0] 		     , (name+".PxPruned["+name+".Size]/F").c_str()			); 
      root->Branch((name+".PyPruned").c_str()		       , &PyPruned[0] 		     , (name+".PyPruned["+name+".Size]/F").c_str()			); 
      root->Branch((name+".PzPruned").c_str()		       , &PzPruned[0] 		     , (name+".PzPruned["+name+".Size]/F").c_str()			); 
      root->Branch((name+".MassPruned").c_str()	       , &MassPruned[0] 		     , (name+".MassPruned["+name+".Size]/F").c_str()			); 
      root->Branch((name+".AreaPruned").c_str()	       , &AreaPruned[0] 		     , (name+".AreaPruned["+name+".Size]/F").c_str()			); 
      root->Branch((name+".tau1").c_str()	             , &tau1[0] 		       , (name+".tau1["+name+".Size]/F").c_str()			); 
      root->Branch((name+".tau2").c_str()	             , &tau2[0] 		       , (name+".tau2["+name+".Size]/F").c_str()			); 
      root->Branch((name+".tau3").c_str()	             , &tau3[0] 		       , (name+".tau3["+name+".Size]/F").c_str()			); 
      root->Branch((name+".JetIDLOOSE").c_str()	       ,&JetIDLOOSE[0]	     , (name+".JetIDLOOSE["+name+".Size]/O").c_str()		); 
      root->Branch((name+".JetIDTIGHT").c_str()	       ,&JetIDTIGHT[0]	     , (name+".JetIDTIGHT["+name+".Size]/O").c_str()		); 
      root->Branch((name+".JetCharge").c_str()	       , &JetCharge[0]  	   , (name+".JetCharge["+name+".Size]/F").c_str()		);
      root->Branch((name+".QGTagsMLP").c_str()	       , &QGTagsMLP[0]  	   , (name+".QGTagsMLP["+name+".Size]/F").c_str()		);
      root->Branch((name+".QGTagsLikelihood").c_str()	       , &QGTagsLikelihood[0]  	     , (name+".QGTagsLikelihood["+name+".Size]/F").c_str()		);
      root->Branch((name+".NConstituents").c_str()	       , &NConstituents[0]	     , (name+".NConstituents["+name+".Size]/I").c_str()		);
      root->Branch((name+".NCH").c_str()	       , &NCH[0]	     , (name+".NCH["+name+".Size]/I").c_str()		);
      root->Branch((name+".CEF").c_str() 	       	       , &CEF[0]		     , (name+".CEF["+name+".Size]/F").c_str()		);
      root->Branch((name+".NHF").c_str() 	       	       , &NHF[0]		     , (name+".NHF["+name+".Size]/F").c_str()		);
      root->Branch((name+".NEF").c_str() 	       	       , &NEF[0]		     , (name+".NEF["+name+".Size]/F").c_str()		);
      root->Branch((name+".CHF").c_str() 	       	       , &CHF[0]		     , (name+".CHF["+name+".Size]/F").c_str()		);
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
      root->Branch((name+".CombinedSVBJetTags").c_str()      , &CombinedSVBJetTags[0]      , (name+".CombinedSVBJetTags["+name+".Size]/F").c_str()	);
      root->Branch((name+".CombinedSVMVABJetTags").c_str()   , &CombinedSVMVABJetTags[0]   , (name+".CombinedSVMVABJetTags["+name+".Size]/F").c_str()  );
      root->Branch((name+".SoftElecByIP3dBJetTags").c_str()  , &SoftElecByIP3dBJetTags[0]  , (name+".SoftElecByIP3dBJetTags["+name+".Size]/F").c_str()	);
      root->Branch((name+".SoftElecByPtBJetTags").c_str()    , &SoftElecByPtBJetTags[0]    , (name+".SoftElecByPtBJetTags["+name+".Size]/F").c_str()	);
      root->Branch((name+".SoftMuonBJetTags").c_str()        , &SoftMuonBJetTags[0]	     , (name+".SoftMuonBJetTags["+name+".Size]/F").c_str()	);
      root->Branch((name+".SoftMuonByIP3dBJetTags").c_str()  , &SoftMuonByIP3dBJetTags[0]  , (name+".SoftMuonByIP3dBJetTags["+name+".Size]/F").c_str()	);	
      root->Branch((name+".SoftMuonByPtBJetTags").c_str()    , &SoftMuonByPtBJetTags[0]    , (name+".SoftMuonByPtBJetTags["+name+".Size]/F").c_str()	);	
      root->Branch((name+".DoubleSVHighEffBJetTags").c_str() , &DoubleSVHighEffBJetTags[0] , (name+".DoubleSVHighEffBJetTags["+name+".Size]/F").c_str()); 

      root->Branch((name+".GenJetPt").c_str() 	       , &GenJetPt[0]		     , (name+".GenJetPt["+name+".Size]/F").c_str()		);
      root->Branch((name+".GenJetEta").c_str()	       , &GenJetEta[0]  	     , (name+".GenJetEta["+name+".Size]/F").c_str()		);
      root->Branch((name+".GenJetPhi").c_str()	       , &GenJetPhi[0]  	     , (name+".GenJetPhi["+name+".Size]/F").c_str()		);
      root->Branch((name+".GenPt").c_str()		       , &GenPt[0]		     , (name+".GenPt["+name+".Size]/F").c_str()  		);
      root->Branch((name+".GenEta").c_str()		       , &GenEta[0]		     , (name+".GenEta["+name+".Size]/F").c_str() 		);
      root->Branch((name+".GenPhi").c_str()		       , &GenPhi[0]		     , (name+".GenPhi["+name+".Size]/F").c_str() 		);
      root->Branch((name+".GenPdgID").c_str() 	       , &GenPdgID[0]		     , (name+".GenPdgID["+name+".Size]/I").c_str()		);
      root->Branch((name+".GenFlavor").c_str()	       , &GenFlavor[0]  	     , (name+".GenFlavor["+name+".Size]/I").c_str()		);
      root->Branch((name+".GenMCTag").c_str()	       	       , &GenMCTag[0]  	             , (name+".GenMCTag["+name+".Size]/I").c_str()		);

      root->Branch((name+".Jet_FatJetIdx").c_str(),  Jet_FatJetIdx  ,(name+".Jet_FatJetIdx["+name+".Size]/I").c_str());
      root->Branch((name+".Jet_SubJet1Idx").c_str(),  Jet_SubJet1Idx  ,(name+".Jet_SubJet1Idx["+name+".Size]/I").c_str());
      root->Branch((name+".Jet_SubJet2Idx").c_str(),  Jet_SubJet2Idx  ,(name+".Jet_SubJet2Idx["+name+".Size]/I").c_str());

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
      root->SetBranchAddress((name+".EtPruned").c_str()		       , &EtPruned[0] 		    ) ; 
      root->SetBranchAddress((name+".PtPruned").c_str()		       , &PtPruned[0] 		    ) ; 
      root->SetBranchAddress((name+".UncPruned").c_str()		       , &UncPruned[0] 		    ) ; 
      root->SetBranchAddress((name+".EtaPruned").c_str()		       , &EtaPruned[0]		    ) ; 
      root->SetBranchAddress((name+".PhiPruned").c_str()		       , &PhiPruned[0]		    ) ; 
      root->SetBranchAddress((name+".EnergyPruned").c_str()	       , &EnergyPruned[0] 	) ; 
      root->SetBranchAddress((name+".PxPruned").c_str()		       , &PxPruned[0] 		    ) ; 
      root->SetBranchAddress((name+".PyPruned").c_str()		       , &PyPruned[0] 		    ) ; 
      root->SetBranchAddress((name+".PzPruned").c_str()		       , &PzPruned[0] 		    ) ; 
      root->SetBranchAddress((name+".MassPruned").c_str()	       , &MassPruned[0] 		  ) ; 
      root->SetBranchAddress((name+".AreaPruned").c_str()	       , &AreaPruned[0] 		  ) ; 
      root->SetBranchAddress((name+".tau1").c_str()	       , &tau1[0]			); 
      root->SetBranchAddress((name+".tau2").c_str()	       , &tau2[0]			); 
      root->SetBranchAddress((name+".tau3").c_str()	       , &tau3[0]			); 
      root->SetBranchAddress((name+".JetIDLOOSE").c_str()		 , &JetIDLOOSE[0]	       ); 
      root->SetBranchAddress((name+".JetIDTIGHT").c_str()		 , &JetIDTIGHT[0]	       ); 
      root->SetBranchAddress((name+".JetCharge").c_str()		 , &JetCharge[0]	       );
      root->SetBranchAddress((name+".QGTagsMLP").c_str()		 , &QGTagsMLP[0]	       );
      root->SetBranchAddress((name+".QGTagsLikelihood").c_str()		 , &QGTagsLikelihood[0]	       );
      root->SetBranchAddress((name+".NConstituents").c_str()		 , &NConstituents[0]	       );
      root->SetBranchAddress((name+".NCH").c_str()		 , &NCH[0]	       );
      root->SetBranchAddress((name+".CEF").c_str()		 , &CEF[0]  	       );
      root->SetBranchAddress((name+".NHF").c_str()		 , &NHF[0]  	       );
      root->SetBranchAddress((name+".NEF").c_str()		 , &NEF[0]  	       );
      root->SetBranchAddress((name+".CHF").c_str()		 , &CHF[0]  	       );
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
      root->SetBranchAddress((name+".CombinedSVBJetTags").c_str()	 , &CombinedSVBJetTags[0]      );
      root->SetBranchAddress((name+".CombinedSVMVABJetTags").c_str()   , &CombinedSVMVABJetTags[0]   );
      root->SetBranchAddress((name+".SoftElecByIP3dBJetTags").c_str()  , &SoftElecByIP3dBJetTags[0]  );
      root->SetBranchAddress((name+".SoftElecByPtBJetTags").c_str()    , &SoftElecByPtBJetTags[0]    );
      root->SetBranchAddress((name+".SoftMuonBJetTags").c_str()        , &SoftMuonBJetTags[0]	       );
      root->SetBranchAddress((name+".SoftMuonByIP3dBJetTags").c_str()  , &SoftMuonByIP3dBJetTags[0]  );
      root->SetBranchAddress((name+".SoftMuonByPtBJetTags").c_str()    , &SoftMuonByPtBJetTags[0]    ); 
      root->SetBranchAddress((name+".DoubleSVHighEffBJetTags").c_str() , &DoubleSVHighEffBJetTags[0] ); 

      root->SetBranchAddress((name+".GenJetPt").c_str()		 , &GenJetPt[0] 	       );
      root->SetBranchAddress((name+".GenJetPt").c_str()		 , &GenJetPt[0] 	       );
      root->SetBranchAddress((name+".GenJetEta").c_str()		 , &GenJetEta[0]	       );
      root->SetBranchAddress((name+".GenJetPhi").c_str()		 , &GenJetPhi[0]	       );
      root->SetBranchAddress((name+".GenPt").c_str()  		 , &GenPt[0]		       );
      root->SetBranchAddress((name+".GenEta").c_str() 		 , &GenEta[0]		       );
      root->SetBranchAddress((name+".GenPhi").c_str() 		 , &GenPhi[0]		       );
      root->SetBranchAddress((name+".GenPdgID").c_str()		 , &GenPdgID[0] 	       );
      root->SetBranchAddress((name+".GenFlavor").c_str()		 , &GenFlavor[0]	       );
      root->SetBranchAddress((name+".GenMCTag").c_str()		 , &GenMCTag[0] 	       );
      root->SetBranchAddress((name+".Px").c_str()			 , &Px[0]		       ); 
      root->SetBranchAddress((name+".Py").c_str()			 , &Py[0]		       ); 
      root->SetBranchAddress((name+".Pz").c_str()			 , &Pz[0]		       ); 
      root->SetBranchAddress((name+".Energy").c_str()		 , &Energy[0]		       ); 

      root->SetBranchAddress((name+".Mass").c_str()		 , &Mass[0]		       );
      root->SetBranchAddress((name+".Area").c_str()		 , &Area[0]		       );

      root->SetBranchAddress((name+".Jet_FatJetIdx").c_str(),  Jet_FatJetIdx );
      root->SetBranchAddress((name+".Jet_SubJet1Idx").c_str(),  Jet_SubJet1Idx );
      root->SetBranchAddress((name+".Jet_SubJet2Idx").c_str(),  Jet_SubJet2Idx );

    }

};





#endif

