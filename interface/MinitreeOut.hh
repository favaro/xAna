#ifndef __minitree_h
#define __minitree_h

#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

class MiniTree {
 public:
  /// TFile
  TFile *mtree_file;
  /// main tree with di-photon system variable
  TChain *mtree_mainTree;
  TChain *mtree_miniTree;
  TChain *mtree_mcTrue;

  /// dedicated tree to specific channel analysis
  TChain *mtree_elecTree;
  TChain *mtree_muonTree;
  TChain *mtree_dijetTree;

  bool mtree_fillLepTree;
  bool mtree_fillDijetTree;
  
  bool addSyncVariables;      

  int mtree_ievt; /// pointer to event in main tree

  /// for minitree everything is done in the .hh file so declare everything inline
  inline MiniTree( const char * filename, bool read = false );
  inline ~MiniTree(void);
  inline void createBranches(void);
  inline void setBranchesAddresses(void);
  inline void initEvent(void);
  inline void initRead(void);
  inline void fill(void);
  inline void fillMCtrueOnly(void);
  inline void end(void);

  ///////////////////////////////////////////////////////////
  /////////////// Declare all variables /////////////////////
  ///////////////////////////////////////////////////////////

  /// MC weights
  float mc_wei;
  float mc_wBSz;
  float mc_wPhoEffi;
  float mc_wTrigEffi;
  float mc_wPU, mc_wVtxId;
  float mc_wHQT;
  float mc_wXsec, mc_wNgen;

  /// truth variables
  float mc_mH, mc_cThetaStar_CS;
  float mc_H_Pt;
  float mc_H_Eta;
  float mc_H_Phi;
  float mc_H_E;
  float mc_partonOut_Pt1;
  float mc_partonOut_Pt2;
  float mc_partonOut_Phi1;
  float mc_partonOut_Phi2;
  float mc_partonOut_Eta1;
  float mc_partonOut_Eta2;
  float mc_partonOut_E1;
  float mc_partonOut_E2;

  float mc_partonOut_deltaPhi;
  float mc_partonOut_deltaPhiRF;

  float mc_genPhoEMatched[2];

  /// event id
  Int_t    mtree_runNum;
  Long64_t mtree_evtNum;
  Int_t    mtree_lumiSec;
  Float_t  weight;
  
  /// general variables
  Float_t mtree_rho , mtree_rho25;
  Float_t mtree_zVtx;
  Int_t   mtree_nVtx;
  Int_t   mtree_nVtxNoBS; // CF
  Float_t nVtxFloat;  // CF
  Float_t nVtxNoBSFloat;   // CF
  Float_t catFloat; //JM
  Int_t   mtree_ivtx1, mtree_ivtx2, mtree_ivtx3;
  Float_t mtree_mass, mtree_pt, mtree_piT, mtree_massDefVtx;
  Float_t mtree_y;
  Float_t mtree_massResoTot,mtree_massResoAng,mtree_massResoEng,mtree_massResoEngPosErr,mtree_massResoEngNegErr;
  Float_t mtree_massResoRightVtx, mtree_massResoRandVtx;
  Float_t mtree_vtxProb, mtree_vtxMva;
  Float_t mtree_diphoMva;
  Float_t mtree_rawMet, mtree_rawMetPhi;
  Float_t mtree_corMet, mtree_corMetPhi;
  Int_t   mtree_catMva;
  // **** CF *****
  Int_t   mtree_tagCat;
  Int_t   mtree_untagCat;
  // *************
  Int_t   mtree_catBase;
  Int_t   mtree_vbfTag;
  Int_t   mtree_vbfCat;
  Int_t   mtree_metTag;
  Int_t   mtree_lepTag;
  Int_t   mtree_lepCat;
  Int_t   mtree_hvjjTag;  
  

  /// for photon specific use array[2]
  /// + pho1, pho2 only for variables going in the fitter (roofit)
  float mtree_eg[2],mtree_ptg[2], mtree_eta[2], mtree_phi[2],mtree_sceta[2];
  float mtree_relResOverE[2], mtree_relSmearing[2];
  int   mtree_isElec[2];
  int   mtree_cat[2];  

  /// spin specific
  float mtree_cThetaLead_heli, mtree_cThetaTrail_heli,mtree_cThetaStar_CS;
  float mtree_pt1, mtree_pt2, mtree_mvaid1, mtree_mvaid2, mtree_eta1, mtree_eta2, mtree_r91, mtree_r92;
  float mtree_pit1, mtree_pit2;
  float mtree_minR9, mtree_minPhoIdEB, mtree_minPhoIdEE;
  float mtree_minSCEta, mtree_maxSCEta;

  /// phoid variables
  float mtree_mvaid[2], mtree_s4Ratio[2];
  float mtree_phoIso1[2], mtree_phoIso2[2], mtree_phoIso3[2]; 
  float mtree_coviEiE[2], mtree_coviEiP[2], mtree_r9[2], mtree_drTrk[2];
  float mtree_etCorrEcalIso[2], mtree_etCorrHcalIso[2], mtree_etCorrTrkIso[2];
  float mtree_pfChgIso02_selVtx[2];
  float mtree_esEffSigmaRR[2];
  float mtree_HoE[2];
    
  /// photon PF isolation
  float mtree_pfPhoIso03[2], mtree_pfPhoIso04[2];
  float mtree_pfChgIso03_selVtx[2], mtree_pfChgIso04_selVtx[2];
  float mtree_pfChgIso03_badVtx[2], mtree_pfChgIso03_defVtx[2];
  float mtree_pfChgIso04_badVtx[2], mtree_pfChgIso04_defVtx[2];

  
  /// sync only variables (float only var for simpl)
  float convInd[2], convNtk[2], convValidVtx[2];
  float convChi2Prob[2], convPt[2];

  float phoID_covIEtaIPhi[2], phoID_S4[2], phoID_ESEffSigmaRR[2],  phoID_pfPhoIso03[2];
  float phoID_pfChIso03[2], phoID_pfChIso03worst[2],phoID_E2x2[2], phoID_E5x5[2];
  float mtree_enScale[2];    // CF
  float mtree_phoSCRawE[2];  // CF
  float mtree_phoSmear[2];   // CF
  float mtree_phoRegrE[2];     
      
 
  float vtxId[3];
  float vtxMva[3]; 
  float vtxDz[3];
  float vtxPtBal, vtxAsym, vtxLogSumPt2, vtxPullToConv, vtxNConv;

  /// dijet selection variables
  float dijet_nJet;
  float dijet_dEtaJJ, dijet_dPhiJJ, dijet_zeppenfeld, dijet_mJJ, dijet_dphiGGJJ;
  float dijet_mvaVbf, dijet_mvaVbf_sig, dijet_mvaVbf_dipho, dijet_mvaVbf_gluglu, dijet_mvaVbfCombDipho;
  float dijet_ptj1, dijet_ptj2, dijet_etaj1, dijet_etaj2, dijet_jec1, dijet_jec2; 
  float dijet_jecUncrtPosj1, dijet_jecUncrtNegj1, dijet_jecUncrtPosj2, dijet_jecUncrtNegj2;
  float dijet_ptj[2], dijet_jec[2], dijet_etj[2], dijet_jecUncrtPos[2], dijet_jecUncrtNeg[2]; 
  float dijet_eta[2], dijet_phi[2], dijet_en[2];
  float dijet_betaStar[2], dijet_rms[2];


  /// MELA for two-jet selection
  // matrix elements:
  float mela_dXsec_HJJ;
  float mela_dXsec_HJJVBF;
  float mela_dXsec_HJJ_PS;
  float mela_dXsec_HJJVBF_PS;
  // discriminants
  float mela_VBFvsgg;
  float mela_SMvsPS_VBF;
  float mela_SMvsPS_gg;
  
  // dilepton variables
  float lept_nMu, lept_nElec;
  float ele_pt[2], ele_eta[2], ele_phi[2];
  float mu_pt[2],  mu_eta[2],  mu_phi[2];

  // JM Add ttH Jets to tree for Sync *******************************************************************
  // This is awful but it seems impossible to read and write vectors with only pointers or only objects....

  int jets_nbtag;
  int jets_njets;
  Float_t jets_nbtagFloat; 
  Float_t jets_njetsFloat; 

  vector<float> *pjets_bjet_csv;
  vector<float> *pjets_pt;
  vector<int> *pjets_index;
  /*
  vector<float> *pjets_et;
  vector<float> *pjets_eta;
  vector<float> *pjets_phi;
  vector<float> *pjets_en;
  vector<float> *pjets_betaStar;
  vector<float> *pjets_rms;
  vector<float> *pjets_z;
  */
  TBranch *brjets_bjet_csv;
  TBranch *brjets_pt;
  TBranch *brjets_index;
  
  vector<float> jets_bjet_csv;
  vector<float> jets_pt;
  vector<int> jets_index;
  vector<float> jets_et;
  vector<float> jets_eta;
  vector<float> jets_phi;
  vector<float> jets_en;
  vector<float> jets_betaStar;
  vector<float> jets_rms;
  vector<float> jets_z;
  

  //  ****************************************************************************


  /// vector to do synchronisation printout
  vector<string> sync_iName, sync_lName,sync_fName, sync_mtreeName;
  vector<int*> sync_iVal;
  vector<Long64_t*> sync_lVal;
  vector<float*> sync_fVal;
  //    virtual void loopSync(void);
  inline void setSynchVariables( string );

};


MiniTree::MiniTree( const char * filenames, bool read ) {
  addSyncVariables = true;
  
  if( !read ) {
    mtree_ievt = 0;
    mtree_file = new TFile( filenames,"RECREATE", "H->gg input tree for unbinned maximum-likelihood fit");
    
    mtree_mainTree  = (TChain*) new TTree( "HToGG","main minitree for unbinned fit" );
    mtree_elecTree  = (TChain*) new TTree( "HToGG_elevar","electron variables ONLY" );
    mtree_muonTree  = (TChain*) new TTree( "HToGG_muovar","muon variables ONLY"     );
    mtree_dijetTree = (TChain*) new TTree( "HToGG_dijet" ,"dijet variables ONLY"    );
    mtree_mcTrue    = (TChain*) new TTree( "HToGG_mcTrue" ,"dijet variables ONLY"    );
	
    createBranches();
    initEvent();
  } else {
    /// reread already existing minitree
    mtree_mainTree  = new TChain( "HToGG" );
    mtree_elecTree  = new TChain( "HToGG_elevar");
    mtree_muonTree  = new TChain( "HToGG_muovar");
    mtree_dijetTree = new TChain( "HToGG_dijet" );
    mtree_mcTrue    = new TChain( "HToGG_mcTrue" );

    /// adding files in filenames may need a loop if file name too complex
    mtree_mainTree ->Add(filenames);
    mtree_elecTree ->Add(filenames);
    mtree_muonTree ->Add(filenames);
    mtree_dijetTree->Add(filenames);
    mtree_mcTrue   ->Add(filenames);

    setBranchesAddresses();
  }
}

void MiniTree::fillMCtrueOnly(void) {
   mtree_mcTrue->Fill();
}

void MiniTree::fill(void) {
  /// roofti doesnot really like arrays
  mtree_pt1  = mtree_ptg[0];
  mtree_pt2  = mtree_ptg[1];
  mtree_pit1 = mtree_ptg[0]/mtree_mass;
  mtree_pit2 = mtree_ptg[1]/mtree_mass;

  mtree_mainTree->Fill();
  if( mtree_fillLepTree && mtree_lepCat  == 0 ) mtree_elecTree->Fill();
  if( mtree_fillLepTree && mtree_lepCat  == 1 ) mtree_muonTree->Fill();
  if( mtree_fillDijetTree                     ) mtree_dijetTree->Fill();
  
  // increase internal ievt pointer
  mtree_ievt++;
}


void MiniTree::initEvent(void) {

  mtree_fillLepTree = false;
  mtree_fillDijetTree = false;

  /// MC weights
  mc_wei = mc_wBSz = 1;
  mc_wPhoEffi = mc_wTrigEffi = 1;
  mc_wPU   = 1;
  mc_wHQT  = 1; 
  mc_wXsec = 1;
  mc_wNgen = 1;
  mc_wVtxId = 1;
  
  /// MC truth variables
  mc_mH = -1;
  mc_cThetaStar_CS = -999;
  mc_H_Pt = -999;
  mc_H_Eta = -999;
  mc_H_Phi = -999;
  mc_H_E = -999;
  mc_partonOut_Pt1 = -999;
  mc_partonOut_Pt2 = -999;
  mc_partonOut_Phi1 = -999;
  mc_partonOut_Phi2 = -999;
  mc_partonOut_Eta1 = -999;
  mc_partonOut_Eta2 = -999;
  mc_partonOut_E1 = -999;
  mc_partonOut_E2 = -999;

  /// event id
  mtree_runNum = mtree_evtNum = mtree_lumiSec = -1;
  
  /// weight section  
  weight = 1;
  
  /// general variables
  mtree_rho = mtree_rho25 = -1;
  mtree_zVtx = -999;
  mtree_nVtx = -1;
  mtree_ivtx1 = mtree_ivtx2 = mtree_ivtx3 = -1;
  mtree_mass  = mtree_pt = mtree_piT = mtree_massDefVtx = -1;
  mtree_massResoTot = mtree_massResoAng = mtree_massResoEng = mtree_massResoEngPosErr = mtree_massResoEngNegErr = -1;
  mtree_massResoRightVtx = mtree_massResoRandVtx = -1;
  mtree_y = -999;
  mtree_vtxProb  = -1;
  mtree_vtxMva   = -1;
  mtree_diphoMva = -1;
  mtree_rawMet = mtree_rawMetPhi = -1;
  mtree_corMet = mtree_corMetPhi = -1;
  mtree_catMva   = -1;
  mtree_tagCat   = -1;  // CF
  mtree_untagCat = -1;  // CF
  mtree_catBase=-1;
  /// set exclusive tagging to zero
  mtree_vbfTag = mtree_metTag = mtree_lepTag = 0;
  mtree_hvjjTag = 0;  
  mtree_vbfCat = -1;
  mtree_lepCat = -1;
   
 mtree_cThetaLead_heli = mtree_cThetaTrail_heli = mtree_cThetaStar_CS = -999;

 lept_nElec = lept_nMu = -1;

  /// photon specific use array[2]
  for( int i = 0; i < 2; i++ ) {
    
    mc_genPhoEMatched[i] = -999;
    
    mtree_eg[i] = mtree_ptg[i] = mtree_eta[i] = mtree_phi[i] = mtree_sceta[i] = -999;
    mtree_relResOverE[i] = mtree_relSmearing[i] = -999;
    mtree_cat[i]    = -1;
    mtree_isElec[i] = -1;
    
    mtree_pt1 = mtree_pt2 = -999;
    mtree_pit1 = mtree_pit2 = -999;
    mtree_mvaid1 = mtree_mvaid2 = -999;
    mtree_eta1 = mtree_eta2 = -999;
    mtree_r91 = mtree_r92 = -999;
    
    /// phoid variables
    mtree_mvaid[i] = mtree_s4Ratio[i] = -999;
    mtree_phoIso1[i] = mtree_phoIso2[i] = mtree_phoIso3[i] = -999; 
    mtree_coviEiE[i] = mtree_coviEiP[i] = mtree_r9[i] = mtree_drTrk[i] = -999;
    mtree_etCorrEcalIso[i] = mtree_etCorrHcalIso[i] = mtree_etCorrTrkIso[i] = -999;
    mtree_pfChgIso02_selVtx[i] = -999;
    mtree_esEffSigmaRR[i] = -999;
    mtree_HoE[i] = -999;
    
    /// photon PF isolation
    mtree_pfPhoIso03[i] = mtree_pfPhoIso04[i] = -999;
    mtree_pfChgIso03_selVtx[i] = mtree_pfChgIso04_selVtx[i] = -999;
    mtree_pfChgIso03_badVtx[i] = mtree_pfChgIso03_defVtx[i] = -999;
    mtree_pfChgIso04_badVtx[i] = mtree_pfChgIso04_defVtx[i] = -999;
    
    /// dijet selection variables
    dijet_dEtaJJ = dijet_dPhiJJ = dijet_zeppenfeld = dijet_mJJ = dijet_dphiGGJJ = dijet_ptj1 = dijet_ptj2 = -999;
    dijet_etaj1 = dijet_etaj2 = dijet_jec1 = dijet_jec2 = -999;
    dijet_jecUncrtPosj1 = dijet_jecUncrtPosj2 = dijet_jecUncrtNegj1 = dijet_jecUncrtNegj2 = -999;
    dijet_ptj[i] = dijet_jec[i] = dijet_etj[i] = dijet_jecUncrtPos[i] = dijet_jecUncrtNeg[i] = -999;
    dijet_eta[i] = dijet_phi[i] = dijet_en[i] = -999;
    dijet_betaStar[i] = dijet_rms[i] = -999;

    ele_eta[i] = ele_pt[i] = ele_phi[i] = -999;
    mu_eta[i] = mu_pt[i] = mu_phi[i] = -999;

    if( addSyncVariables ) {
      convInd[i] = convNtk[i] = convValidVtx[i] = -999;
      convChi2Prob[i] = convPt[i] = -999;      
      
      phoID_covIEtaIPhi[i] = phoID_S4[i] = phoID_ESEffSigmaRR[i] = phoID_pfPhoIso03[i] = -999;
      phoID_pfChIso03[i]   = phoID_pfChIso03worst[i] = phoID_E2x2[i] = phoID_E5x5[i] = -999;

    }
  }

  // MELA discriminants:
  mela_dXsec_HJJ = mela_dXsec_HJJVBF =  mela_dXsec_HJJ_PS =  mela_dXsec_HJJVBF_PS = -999;      
  mela_VBFvsgg    = -999;
  mela_SMvsPS_VBF = -999;
  mela_SMvsPS_gg  = -999;
   
  if( addSyncVariables ) {
    vtxPtBal =vtxAsym = vtxLogSumPt2 =  vtxPullToConv = -999;
    for( int iv = 0 ; iv < 3; iv ++ ) vtxId[iv] = vtxMva[iv] = vtxDz[iv] = -999;
  }

  dijet_mvaVbf        = -1;
  dijet_mvaVbf_sig    = -1;
  dijet_mvaVbf_dipho  = -1;
  dijet_mvaVbf_gluglu = -1;
  dijet_mvaVbfCombDipho = -1;
    
  // JM Add ttH Jets to tree for Sync ******************************************

  jets_nbtag=0;
  jets_njets=0;

  
  
  jets_pt.clear();
  jets_et.clear();
  jets_en.clear();
  jets_eta.clear();
  jets_phi.clear();
  jets_betaStar.clear();
  jets_rms.clear();
  jets_z.clear();
  jets_bjet_csv.clear();
  jets_index.clear();

  //if( avocado > -999 ) cout << " avocado end = " <<  avocado << endl;

  // *******************************************************************
}

void MiniTree::end(void) {

  
    mtree_file->Write() ;
    mtree_file->Close(); 
}

MiniTree::~MiniTree(void) {  

}



void MiniTree::createBranches(void) {
  cout << " ------------- MiniTree set output branches ------------" << endl;
  mtree_mcTrue->Branch("mH"           , &mc_mH          , "mH/F");
  mtree_mcTrue->Branch("cThetaStar_CS", &mc_cThetaStar_CS, "cThetaStar_CS/F");
    
  //----------------------- MC weights
  mtree_mainTree->Branch("mH"       , &mc_mH      , "mH/F");
  mtree_mainTree->Branch("cThetaStar_CS_truth", &mc_cThetaStar_CS, "cThetaStar_CS_truth/F");

  mtree_mainTree->Branch("wei"      , &mc_wei      , "wei/F");
  mtree_mainTree->Branch("wPU"      , &mc_wPU      , "wPU/F");
  mtree_mainTree->Branch("wHQT"     , &mc_wHQT     , "wHQT/F");
  mtree_mainTree->Branch("wXsec"    , &mc_wXsec    , "wXsec/F");
  mtree_mainTree->Branch("wNgen"    , &mc_wNgen    , "wNgen/F");
  mtree_mainTree->Branch("wBSz"     , &mc_wBSz     , "wBSz/F");
  mtree_mainTree->Branch("wVtxId"   , &mc_wVtxId   , "wVtxId/F");
  mtree_mainTree->Branch("wPhoEffi" , &mc_wPhoEffi , "wPhoEffi/F");
  mtree_mainTree->Branch("wTrigEffi", &mc_wTrigEffi, "wTrigEffi/F");

  mtree_mainTree->Branch("genPhoEMatched", mc_genPhoEMatched, "genPhoEMatched[2]/F");

  // MC information for study of correlations at gen level 
  mtree_mcTrue->Branch("mc_H_Pt"        , &mc_H_Pt        , "H_Pt/F");
  mtree_mcTrue->Branch("mc_H_Eta"       , &mc_H_Eta       , "H_Eta/F");
  mtree_mcTrue->Branch("mc_H_Phi"       , &mc_H_Phi       , "H_Phi/F");  
  mtree_mcTrue->Branch("mc_H_E"         , &mc_H_E         , "H_E/F");

  mtree_mcTrue->Branch("mc_partonOut_Pt1"   , &mc_partonOut_Pt1     , "partonOut_Pt1/F"); 
  mtree_mcTrue->Branch("mc_partonOut_Pt2"   , &mc_partonOut_Pt2     , "partonOut_Pt2/F");  
  mtree_mcTrue->Branch("mc_partonOut_Eta1"  , &mc_partonOut_Eta1    , "partonOut_Eta1/F"); 
  mtree_mcTrue->Branch("mc_partonOut_Eta2"  , &mc_partonOut_Eta2    , "partonOut_Eta2/F");
  mtree_mcTrue->Branch("mc_partonOut_Phi1"  , &mc_partonOut_Phi1    , "partonOut_Phi1/F"); 
  mtree_mcTrue->Branch("mc_partonOut_Phi2"  , &mc_partonOut_Phi2    , "partonOut_Phi2/F");
  mtree_mcTrue->Branch("mc_partonOut_E1"    , &mc_partonOut_E1      , "partonOut_E1/F"); 
  mtree_mcTrue->Branch("mc_partonOut_E2"    , &mc_partonOut_E2      , "partonOut_E2/F");
  
  mtree_mcTrue->Branch("mc_partonOut_deltaPhi"   , &mc_partonOut_deltaPhi     , "partonOut_deltaPhi/F"); 
  mtree_mcTrue->Branch("mc_partonOut_deltaPhiRF" , &mc_partonOut_deltaPhiRF   , "partonOut_deltaPhiRF/F"); 

  //-----------------------
  mtree_mainTree->Branch("iMainEvt", &mtree_ievt    , "iMainEvt/I" ); ///pointer to main tree
  mtree_mainTree->Branch("runNum"  , &mtree_runNum  , "runNum/I"   );
  mtree_mainTree->Branch("evtNum"  , &mtree_evtNum  , "evtNum/L" );
  mtree_mainTree->Branch("lumiSec" , &mtree_lumiSec , "lumiSec/I"  );
  mtree_mainTree->Branch("rho"     , &mtree_rho     , "rho/F");
  mtree_mainTree->Branch("rho25"   , &mtree_rho25   , "rho25/F");

  // primary vertex variables *************************************************************

  mtree_mainTree->Branch("zVtx"    , &mtree_zVtx  , "zVtx/F");
  mtree_mainTree->Branch("nVtx"    , &mtree_nVtx  , "nVtx/I");
  mtree_mainTree->Branch("nVtxNoBS"  , &mtree_nVtxNoBS, "nVtxNoBS/I");  // CF
  mtree_mainTree->Branch("iVtx1"   , &mtree_ivtx1 , "iVtx1/I");
  mtree_mainTree->Branch("iVtx2"   , &mtree_ivtx2 , "iVtx2/I");
  mtree_mainTree->Branch("iVtx3"   , &mtree_ivtx3 , "iVtx3/I");
  mtree_mainTree->Branch("vtxProb" , &mtree_vtxProb, "vtxProb/F");
  mtree_mainTree->Branch("vtxMva"  , &mtree_vtxMva , "vtxMva/F");
  
  // ************************************************************************************

  // skeleton
  //  mtree_mainTree->Branch("x", &mtree_x ,"x/F"  );
  //--------------- diphoton event variables

  mtree_mainTree->Branch("mass"      , &mtree_mass       , "mass/F");
  mtree_mainTree->Branch("massDefVtx", &mtree_massDefVtx , "massDefVtx/F");
  mtree_mainTree->Branch("pt"        , &mtree_pt         , "pt/F"   );
  mtree_mainTree->Branch("y"         , &mtree_y          , "y/F"    );
  mtree_mainTree->Branch("piT"       , &mtree_piT        , "piT/F"  );
  mtree_mainTree->Branch("pit1"      , &mtree_pit1       , "pit1/F"  );
  mtree_mainTree->Branch("pit2"      , &mtree_pit2       , "pit2/F"  );
  mtree_mainTree->Branch("minR9"     , &mtree_minR9      , "minR9/F"  );
  mtree_mainTree->Branch("maxSCEta"  , &mtree_maxSCEta   , "maxSCEta/F"  );
  mtree_mainTree->Branch("minSCEta"  , &mtree_minSCEta   , "minSCEta/F"  );
  mtree_mainTree->Branch("minPhoIdEB", &mtree_minPhoIdEB , "minPhoIdEB/F"  );
  mtree_mainTree->Branch("minPhoIdEE", &mtree_minPhoIdEE , "minPhoIdEE/F"  );
  mtree_mainTree->Branch("massResoTot"     , &mtree_massResoTot      , "massResoTot/F"  );
  mtree_mainTree->Branch("massResoAng"     , &mtree_massResoAng      , "massResoAng/F"  );
  mtree_mainTree->Branch("massResoEng"     , &mtree_massResoEng      , "massResoEng/F"  );
  mtree_mainTree->Branch("massResoEngPosErr"     , &mtree_massResoEngPosErr      , "massResoEngPosErr/F"  );
  mtree_mainTree->Branch("massResoEngNegErr"     , &mtree_massResoEngNegErr      , "massResoEngNegErr/F"  );
  mtree_mainTree->Branch("massResoRightVtx", &mtree_massResoRightVtx , "massResoRightVtx/F" );
  mtree_mainTree->Branch("massResoRandVtx" , &mtree_massResoRandVtx  , "massResoRandVtx/F"  );
  mtree_mainTree->Branch("diphoMva" , &mtree_diphoMva  , "diphoMva/F"  );
  mtree_mainTree->Branch("rawMet"   , &mtree_rawMet    , "rawMet/F"  );
  mtree_mainTree->Branch("corMet"   , &mtree_corMet    , "corMet/F"  );
  mtree_mainTree->Branch("rawMetPhi", &mtree_rawMetPhi , "rawMetPhi/F"  );
  mtree_mainTree->Branch("corMetPhi", &mtree_corMetPhi , "corMetPhi/F"  );

  mtree_mainTree->Branch("catMva" , &mtree_catMva  , "catMva/I"  );
  mtree_mainTree->Branch("tagCat" ,   &mtree_tagCat  ,   "tagCat/I" );  // CF
  mtree_mainTree->Branch("untagCat" , &mtree_untagCat  , "untagCat/I" );  // CF
  mtree_mainTree->Branch("catBase", &mtree_catBase , "catBase/I" );
  mtree_mainTree->Branch("lepTag" , &mtree_lepTag  , "lepTag/I"  );
  mtree_mainTree->Branch("vbfTag" , &mtree_vbfTag  , "vbfTag/I"  );
  mtree_mainTree->Branch("metTag" , &mtree_metTag  , "metTag/I"  );
  mtree_mainTree->Branch("hvjjTag", &mtree_hvjjTag , "hvjjTag/I" );
  mtree_mainTree->Branch("lepCat" , &mtree_lepCat  , "lepCat/I"  );
  mtree_mainTree->Branch("vbfCat" , &mtree_vbfCat  , "vbfCat/I"  );
  mtree_mainTree->Branch("cThetaLead_heli" ,  &mtree_cThetaLead_heli , "cThetaLead_heli/F"  );
  mtree_mainTree->Branch("cThetaTrail_heli",  &mtree_cThetaTrail_heli, "cThetaTrail_heli/F"  );
  mtree_mainTree->Branch("cThetaStar_CS"   ,  &mtree_cThetaStar_CS   , "cThetaStar_CS/F"  );


  // CF: add all jet variables to main tree for synch.  *************************************************

  mtree_mainTree->Branch("nJet", &dijet_nJet , "nJet/F" );

  mtree_mainTree->Branch("jet1Pt", &dijet_ptj1 , "jet1Pt/F" );
  mtree_mainTree->Branch("jet2Pt", &dijet_ptj2 , "jet2Pt/F" );
  mtree_mainTree->Branch("jet1Eta", &dijet_etaj1 , "jet1Eta/F" );
  mtree_mainTree->Branch("jet2Eta", &dijet_etaj2 , "jet2Eta/F" );
  mtree_mainTree->Branch("jec1", &dijet_jec1 , "jec1/F" );
  mtree_mainTree->Branch("jec2", &dijet_jec2 , "jec2/F" );
  mtree_mainTree->Branch("jecUncrtPosj1", &dijet_jecUncrtPosj1 , "jecUncrtPosj1/F" );
  mtree_mainTree->Branch("jecUncrtPosj2", &dijet_jecUncrtPosj2 , "jecUncrtPosj2/F" );
  mtree_mainTree->Branch("jecUncrtNegj1", &dijet_jecUncrtNegj1 , "jecUncrtNegj1/F" );
  mtree_mainTree->Branch("jecUncrtNegj2", &dijet_jecUncrtNegj2 , "jecUncrtNegj2/F" );

  mtree_mainTree->Branch("etj", dijet_etj , "etj[2]/F" );
  mtree_mainTree->Branch("ptj", dijet_ptj , "ptj[2]/F" );
  mtree_mainTree->Branch("jec", dijet_jec , "jec[2]/F" );
  mtree_mainTree->Branch("jecUncrtPos", dijet_jecUncrtPos , "jecUncrtPos[2]/F" );
  mtree_mainTree->Branch("jecUncrtNeg", dijet_jecUncrtNeg , "jecUncrtNeg[2]/F" );
  mtree_mainTree->Branch("enj" , dijet_en  , "enj[2]/F"  );
  mtree_mainTree->Branch("etaj", dijet_eta , "etaj[2]/F" );
  mtree_mainTree->Branch("phij", dijet_phi , "phij[2]/F" );
  mtree_mainTree->Branch("rmsj", dijet_rms , "rmsj[2]/F" );
  mtree_mainTree->Branch("betaStar", dijet_betaStar , "betaStar[2]/F"  );

  mtree_mainTree->Branch("dEtaJJ"    , &dijet_dEtaJJ     , "dEtaJJ/F"    );
  mtree_mainTree->Branch("dPhiJJ"    , &dijet_dPhiJJ     , "dPhiJJ/F"    );
  mtree_mainTree->Branch("zeppenfeld", &dijet_zeppenfeld , "zeppenfeld/F");
  mtree_mainTree->Branch("mJJ"       , &dijet_mJJ        , "mJJ/F"       );
  mtree_mainTree->Branch("dphiGGJJ"  , &dijet_dphiGGJJ   , "dphiGGJJ/F"  );
 
  mtree_mainTree->Branch( "mvaVbf",          &dijet_mvaVbf,         "mvaVbf/F");
  mtree_mainTree->Branch( "mvaVbf_sig"   ,   &dijet_mvaVbf_sig,     "mvaVbf_sig/F"   );
  mtree_mainTree->Branch( "mvaVbf_dipho" ,   &dijet_mvaVbf_dipho,   "mvaVbf_dipho/F" );
  mtree_mainTree->Branch( "mvaVbf_gluglu",   &dijet_mvaVbf_gluglu , "mvaVbf_gluglu/F");
  mtree_mainTree->Branch( "mvaVbfCombDipho", &dijet_mvaVbfCombDipho,"mvaVbfCombDipho/F");

  mtree_mainTree->Branch( "dXsec_HJJ",       &mela_dXsec_HJJ,       "dXsec_HJJ/F");
  mtree_mainTree->Branch( "dXsec_HJJVBF",    &mela_dXsec_HJJVBF,    "dXsec_HJJVBF/F");
  mtree_mainTree->Branch( "dXsec_HJJ_PS",    &mela_dXsec_HJJ_PS,     "dXsec_HJJ_PS/F");
  mtree_mainTree->Branch( "dXsec_HJJVBF_PS", &mela_dXsec_HJJVBF_PS, "dXsec_HJJVBF_PS/F");

  mtree_mainTree->Branch( "mela_VBFvsgg",    &mela_VBFvsgg,    "mela_VBFvsgg/F");
  mtree_mainTree->Branch( "mela_SMvsPS_VBF", &mela_SMvsPS_VBF, "mela_SMvsPS_VBF/F");
  mtree_mainTree->Branch( "mela_SMvsPS_gg",  &mela_SMvsPS_gg,  "mela_SMvsPS_gg/F");

  //***********************************************************************************************

  mtree_mainTree->Branch("eg" , mtree_eg , "eg[2]/F"  );
  mtree_mainTree->Branch("ptg", mtree_ptg, "ptg[2]/F" );
  mtree_mainTree->Branch("eta", mtree_eta, "eta[2]/F" );
  mtree_mainTree->Branch("phi", mtree_phi, "phi[2]/F" );
  mtree_mainTree->Branch("cat", mtree_cat, "cat[2]/I" );
  mtree_mainTree->Branch("sceta", mtree_sceta ,"sceta[2]/F"  );
  mtree_mainTree->Branch("relResOverE", mtree_relResOverE ,"relResOverE[2]/F"  );
  mtree_mainTree->Branch("relSmearing" , mtree_relSmearing  ,"relSmearing[2]/F"  );
  mtree_mainTree->Branch("isElec", mtree_isElec ,"isElec[2]/I"  );
  
/// phoid variables
  mtree_mainTree->Branch("mvaid"  , mtree_mvaid   , "mvaid[2]/F"    );
  mtree_mainTree->Branch("s4Ratio", mtree_s4Ratio , "s4Ratio[2]/F"  );
  mtree_mainTree->Branch("phoIso1", mtree_phoIso1 , "phoIso1[2]/F"  );
  mtree_mainTree->Branch("phoIso2", mtree_phoIso2 , "phoIso2[2]/F"  );
  mtree_mainTree->Branch("phoIso3", mtree_phoIso3 , "phoIso3[2]/F"  );
  mtree_mainTree->Branch("coviEiE", mtree_coviEiE , "coviEiE[2]/F"  );
  mtree_mainTree->Branch("coviEiP", mtree_coviEiP , "coviEiP[2]/F"  );
  mtree_mainTree->Branch("r9"     , mtree_r9      , "r9[2]/F"  );
  mtree_mainTree->Branch("HoE"    , mtree_HoE     , "HoE[2]/F"  );
  mtree_mainTree->Branch("drTrk"  , mtree_drTrk   , "drTrk[2]/F"  );
  mtree_mainTree->Branch("mvaid1" , &mtree_mvaid1 , "mvaid1/F"    );
  mtree_mainTree->Branch("mvaid2" , &mtree_mvaid2 , "mvaid2/F"    );
  mtree_mainTree->Branch("eta1"   , &mtree_eta1   , "eta1/F"    );
  mtree_mainTree->Branch("eta2"   , &mtree_eta2   , "eta2/F"    );
  mtree_mainTree->Branch("r91"    , &mtree_r91    , "r91/F"    );
  mtree_mainTree->Branch("r92"    , &mtree_r92    , "r92/F"    );

  // JM Add ttH Jets to tree for Sync ****************************************  
 
  mtree_mainTree->Branch("jets_njets",    &jets_njets,"jets_njets/I");
  mtree_mainTree->Branch("jets_nbtag", &jets_nbtag, "jets_nbtag/I");
  mtree_mainTree->Branch("jets_bjet_csv", &jets_bjet_csv );  
  mtree_mainTree->Branch("jets_index",&jets_index );
  mtree_mainTree->Branch("jets_pt",    &jets_pt );
  mtree_mainTree->Branch("jets_et",    &jets_et );
  mtree_mainTree->Branch("jets_en",    &jets_en );
  mtree_mainTree->Branch("jets_eta",   &jets_eta );
  mtree_mainTree->Branch("jets_phi",   &jets_phi );
  mtree_mainTree->Branch("jets_betaStar", &jets_betaStar );
  mtree_mainTree->Branch("jets_rms", &jets_rms );
  mtree_mainTree->Branch("jets_z", &jets_z );
 
  // *************************************************************************  

  /// sync only variables
  if( addSyncVariables ) {
    cout << "    ******* adding sync variables to minitree " << endl;
    mtree_mainTree->Branch("etCorrEcalIso", mtree_etCorrEcalIso, "etCorrEcalIso[2]/F");
    mtree_mainTree->Branch("etCorrHcalIso", mtree_etCorrHcalIso, "etCorrHcalIso[2]/F");
    mtree_mainTree->Branch("etCorrTrkIso" , mtree_etCorrTrkIso , "etCorrTrkIso[2]/F" );
    mtree_mainTree->Branch("esEffSigmaRR" , mtree_esEffSigmaRR , "esEffSigmaRR[2]/F"  );
    mtree_mainTree->Branch("enScale"      , mtree_enScale      , "enScale[2]/F" );   // CF
    mtree_mainTree->Branch("phoSCRawE"    , mtree_phoSCRawE    , "phoSCRawE[2]/F" ); // CF
    mtree_mainTree->Branch("phoSmear"     , mtree_phoSmear     , "phoSmear[2]/F" ); // CF
    mtree_mainTree->Branch("phoRegrE"     , mtree_phoRegrE     , "phoRegrE[2]/F" ); // CF

    
    /// photon PF isolation
    mtree_mainTree->Branch("pfChgIso02_selVtx", mtree_pfChgIso02_selVtx , "pfChgIso02_selVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso03_selVtx", mtree_pfChgIso03_selVtx , "pfChgIso03_selVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso03_badVtx", mtree_pfChgIso03_badVtx , "pfChgIso03_badVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso03_defVtx", mtree_pfChgIso03_defVtx , "pfChgIso03_defVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso04_selVtx", mtree_pfChgIso04_selVtx , "pfChgIso04_selVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso04_badVtx", mtree_pfChgIso04_badVtx , "pfChgIso04_badVtx[2]/F"  );
    mtree_mainTree->Branch("pfChgIso04_defVtx", mtree_pfChgIso04_defVtx , "pfChgIso04_defVtx[2]/F"  );
    mtree_mainTree->Branch("pfPhoIso03", mtree_pfPhoIso03 ,"pfPhoIso03[2]/F"  );
    mtree_mainTree->Branch("pfPhoIso04", mtree_pfPhoIso04 ,"pfPhoIso04[2]/F"  );
    
    
    mtree_mainTree->Branch("convInd", convInd, "convInd[2]/F");
    mtree_mainTree->Branch("convNtk", convNtk, "convNtk[2]/F");
    mtree_mainTree->Branch("convPt" , convPt ,"convPt[2]/F" );
    mtree_mainTree->Branch("convValidVtx", convValidVtx,"convValidVtx[2]/F" );
    mtree_mainTree->Branch("convChi2Prob", convChi2Prob, "convChi2Prob[2]/F");
    mtree_mainTree->Branch("vId"  , vtxId   , "vId[3]/F");
    mtree_mainTree->Branch("vDz"  , vtxDz   , "vDz[3]/F");
    mtree_mainTree->Branch("vMva" , vtxMva  , "vMva[3]/F");
    mtree_mainTree->Branch("vAsym"      , &vtxAsym      ,"vAsym/F"      );
    mtree_mainTree->Branch("vPtBal"     , &vtxPtBal     ,"vPtBal/F"     );
    mtree_mainTree->Branch("vLogSumPt2" , &vtxLogSumPt2 ,"vLogSumPt2/F" );
    mtree_mainTree->Branch("vPullToConv", &vtxPullToConv,"vPullToConv/F");
    mtree_mainTree->Branch("vNConv"     , &vtxNConv     ,"vNConv/F"     );
    mtree_mainTree->Branch("phoID_covIEtaIPhi" , phoID_covIEtaIPhi ,"phoID_covIEtaIPhi[2]/F");
    mtree_mainTree->Branch("phoID_S4"          , phoID_S4          ,"phoID_S4[2]/F"         );
    mtree_mainTree->Branch("phoID_ESEffSigmaRR", phoID_ESEffSigmaRR,"phoID_ESEffSigmaRR[2]/F");
    mtree_mainTree->Branch("phoID_pfPhoIso03"  , phoID_pfPhoIso03  ,"phoID_pfPhoIso03[2]/F");
    mtree_mainTree->Branch("phoID_pfChIso03"   , phoID_pfChIso03   ,"phoID_pfChIso03[2]/F");
    mtree_mainTree->Branch("phoID_pfChIso03worst", phoID_pfChIso03worst,"phoID_pfChIso03worst[2]/F");
    mtree_mainTree->Branch("phoID_E2x2", phoID_E2x2,"phoID_E2x2[2]/F");
    mtree_mainTree->Branch("phoID_E5x5", phoID_E5x5,"phoID_E5x5[2]/F");
  }

 

  /// lepton selection variables
  // skeleton

  mtree_mainTree->Branch("nMu",   &lept_nMu ,   "nMu/F" );
  mtree_mainTree->Branch("nElec", &lept_nElec , "nElec/F" );

  mtree_mainTree->Branch("elec_pt",  ele_pt ,  "ele_pt[2]/F" );
  mtree_mainTree->Branch("elec_eta", ele_eta , "ele_eta[2]/F" );
  mtree_mainTree->Branch("elec_phi", ele_phi , "ele_phi[2]/F" );
  mtree_mainTree->Branch("muon_pt",  mu_pt ,   "mu_pt[2]/F" );
  mtree_mainTree->Branch("muon_eta", mu_eta ,  "mu_eta[2]/F" );
  mtree_mainTree->Branch("muon_phi", mu_phi ,  "mu_phi[2]/F" );
  

  mtree_elecTree->Branch("iMainEvt"  , &mtree_ievt       , "iMainEvt/I"  ); ///pointer to main tree
  mtree_elecTree->Branch("mGG"       , &mtree_mass       , "mgg/F"       ); 
  mtree_elecTree->Branch("mvaGG"     , &mtree_diphoMva   , "mvaGG/F"     ); 
  mtree_elecTree->Branch("catBaseGG" , &mtree_catBase    , "catBaseGG/I" ); 

  mtree_muonTree->Branch("iMainEvt"  , &mtree_ievt       , "iMainEvt/I"  ); ///pointer to main tree
  mtree_muonTree->Branch("mGG"       , &mtree_mass       , "mgg/F"       ); 
  mtree_muonTree->Branch("mvaGG"     , &mtree_mvaid      , "mvaGG/F"     ); 
  mtree_muonTree->Branch("catBaseGG" , &mtree_catBase    , "catBaseGG/I" ); 

}







//////////////////////////////////////////////////////////////////////////////////
/////////////////////// MiniTree reading specific functions //////////////////////
//////////////////////////////////////////////////////////////////////////////////

void MiniTree::setSynchVariables( string object ) {
 

  sync_iName.push_back("run"  ); sync_iVal.push_back( &mtree_runNum ); sync_mtreeName.push_back("runNum");
  sync_iName.push_back("lumi" ); sync_iVal.push_back( &mtree_lumiSec); sync_mtreeName.push_back("lumiSec");  
  sync_lName.push_back("event"); sync_lVal.push_back( &mtree_evtNum ); sync_mtreeName.push_back("evtNum");
  //sync_iName.push_back("cat");   sync_iVal.push_back( &mtree_catMva ); sync_mtreeName.push_back("catMva");
  //sync_fName.push_back("rho"  ); sync_fVal.push_back( &mtree_rho    ); sync_mtreeName.push_back("rho");

  sync_iName.push_back("tcat"  );        sync_iVal.push_back( &mtree_tagCat );          sync_mtreeName.push_back("tagCat");
  sync_iName.push_back("ucat"  );        sync_iVal.push_back( &mtree_untagCat );          sync_mtreeName.push_back("untagCat");

  if( object == "vbf" ) {
    sync_fName.push_back("vbf_massJJ"  ); sync_fVal.push_back( &dijet_mJJ );             sync_mtreeName.push_back("vbf_massjj");
    sync_fName.push_back("combiMVA"  );   sync_fVal.push_back( &dijet_mvaVbfCombDipho ); sync_mtreeName.push_back("vbf_ggjjmva");
    sync_fName.push_back("dijetMVA"  );   sync_fVal.push_back( &dijet_mvaVbf ); sync_mtreeName.push_back("vbf_jjmva");
    
    sync_fName.push_back("vbf_dPhiJJGG"  );   sync_fVal.push_back( &dijet_dphiGGJJ ); sync_mtreeName.push_back("vbf_dphijjgg");
    
    sync_fName.push_back("ptom1");    sync_fVal.push_back( &mtree_pit1 ); sync_mtreeName.push_back("ptom1");
    sync_fName.push_back("ptom2");    sync_fVal.push_back( &mtree_pit2 ); sync_mtreeName.push_back("ptom2");

    sync_fName.push_back("numJets");      sync_fVal.push_back( &dijet_nJet );     sync_mtreeName.push_back("njet");
    
    string Ind[2] = {"1","2"};
    for( int i = 0; i<2; i++) {
      sync_fName.push_back("pho" + Ind[i] + "_eta" );  sync_fVal.push_back( &mtree_eta[i]  );         if(i==0 ) sync_mtreeName.push_back("eta");   
      sync_fName.push_back("pho" + Ind[i] + "_phi" );  sync_fVal.push_back( &mtree_phi[i] );          if(i==0 ) sync_mtreeName.push_back("phi");    
      sync_fName.push_back("ele" + Ind[i] + "_eta"  );   sync_fVal.push_back( &ele_eta[i] );       if(i==0 ) sync_mtreeName.push_back("eleeta");
      sync_fName.push_back("ele" + Ind[i] + "_phi"  );   sync_fVal.push_back( &ele_phi[i] );       if(i==0 ) sync_mtreeName.push_back("elephi");
      sync_fName.push_back("jet" + Ind[i] + "_eta"  );  sync_fVal.push_back( &dijet_eta[i] );       if(i==0 ) sync_mtreeName.push_back("jeteta");
    }
  }
  
  //-------------------------------------------
  
  if( object == "PV" ) {

    sync_fName.push_back( "vtxprob" );       sync_fVal.push_back( &mtree_vtxProb ); sync_mtreeName.push_back("vprob");  

    sync_fName.push_back( "vertexId0" );     sync_fVal.push_back( &vtxId[0] );      sync_mtreeName.push_back("vId");
    sync_fName.push_back( "vertexId1" );     sync_fVal.push_back( &vtxId[1] );      sync_mtreeName.push_back("vId1");
    sync_fName.push_back( "vertexId2" );     sync_fVal.push_back( &vtxId[2] );      sync_mtreeName.push_back("vId2");
    
    sync_fName.push_back( "vertex_z" );      sync_fVal.push_back( &mtree_zVtx );    sync_mtreeName.push_back("vz");
    
    sync_fName.push_back( "nVtx" );          sync_fVal.push_back( &nVtxFloat );    sync_mtreeName.push_back("nVtx");
    sync_fName.push_back( "nVtxNoBS" );      sync_fVal.push_back( &nVtxNoBSFloat );    sync_mtreeName.push_back("nVtxNoBS");
    sync_fName.push_back( "vtxNConv" );      sync_fVal.push_back( &vtxNConv );    sync_mtreeName.push_back("vtxNConv");
    
    sync_fName.push_back( "vertexMva0" );       sync_fVal.push_back( &vtxMva[0] );        sync_mtreeName.push_back("vmva");  
    sync_fName.push_back( "vertexMva1" );       sync_fVal.push_back( &vtxMva[1] );        sync_mtreeName.push_back("vmva1");
    sync_fName.push_back( "vertexMva2" );       sync_fVal.push_back( &vtxMva[2] );        sync_mtreeName.push_back("vmva2");  
    
    string cvInd[2] = {"1","2"};
    for( int i = 0; i < 2; i++ ) {
      sync_fName.push_back("convindex"   + cvInd[i] ); sync_fVal.push_back( &convInd[i] );       if(i==0 ) sync_mtreeName.push_back("convInd");
      sync_fName.push_back("convNtrk"    + cvInd[i] ); sync_fVal.push_back( &convNtk[i] );       if(i==0 ) sync_mtreeName.push_back("convNtk");
      sync_fName.push_back("convpt"  +  cvInd[i] );    sync_fVal.push_back( &convPt[i]  );       if(i==0 ) sync_mtreeName.push_back("convPt" );
      sync_fName.push_back("convChiProb" + cvInd[i] ); sync_fVal.push_back( &convChi2Prob[i]  ); if(i==0 ) sync_mtreeName.push_back("convChi2Prob" );
    } 

    sync_fName.push_back("numJets");      sync_fVal.push_back( &dijet_nJet );     sync_mtreeName.push_back("njet");
    sync_fName.push_back("numElec");      sync_fVal.push_back( &lept_nElec );   sync_mtreeName.push_back("nelec");
    sync_fName.push_back("numMu");        sync_fVal.push_back( &lept_nMu );     sync_mtreeName.push_back("nmu");   
  }

  //-------------------------------------------

  if( object == "photons" ) {

    
    sync_fName.push_back("sigmamom"  ); sync_fVal.push_back( &mtree_massResoRightVtx); sync_mtreeName.push_back("massResoRightVtx");
    sync_fName.push_back("sigmamom_wrong_vtx"  ); sync_fVal.push_back( &mtree_massResoRandVtx ); sync_mtreeName.push_back("massResoRandVtx");
    sync_fName.push_back( "vtxprob" );       sync_fVal.push_back( &mtree_vtxProb ); sync_mtreeName.push_back("vprob");  
    // ptom1
    // ptom2
    // dphi
    
    sync_fName.push_back("ptom1");    sync_fVal.push_back( &mtree_pit1 ); sync_mtreeName.push_back("ptom1");
    sync_fName.push_back("ptom2");    sync_fVal.push_back( &mtree_pit2 ); sync_mtreeName.push_back("ptom2");
    
    string phoInd[2] = {"1","2"};
    for( int i = 0; i < 2; i++ ) {
      
      sync_fName.push_back("pho" + phoInd[i] + "_e" );    sync_fVal.push_back( &mtree_eg[i]  ); if(i==0 ) sync_mtreeName.push_back("eg"); 
      // pho_EnScale
      sync_fName.push_back("pho" + phoInd[i] + "_enScale" );    sync_fVal.push_back( &mtree_enScale[i]  ); if(i==0 ) sync_mtreeName.push_back("enScale"); 
      sync_fName.push_back("pho" + phoInd[i] + "_eErr"); sync_fVal.push_back( &mtree_relResOverE[i]  ); if(i==0 ) sync_mtreeName.push_back("relResOverE"); 
      sync_fName.push_back("pho" + phoInd[i] + "_eta" );  sync_fVal.push_back( &mtree_eta[i]  );         if(i==0 ) sync_mtreeName.push_back("eta");
      sync_fName.push_back("pho" + phoInd[i] + "_SCeta" );  sync_fVal.push_back( &mtree_sceta[i]  );         if(i==0 ) sync_mtreeName.push_back("SCeta");
      sync_fName.push_back("pho" + phoInd[i] + "_phi" );  sync_fVal.push_back( &mtree_phi[i] );          if(i==0 ) sync_mtreeName.push_back("phi");    
      sync_fName.push_back("pho" + phoInd[i] + "_r9"  );  sync_fVal.push_back( &mtree_r9[i] );           if(i==0 ) sync_mtreeName.push_back("r9");
      sync_fName.push_back("pho" + phoInd[i] + "_idMVA"); sync_fVal.push_back( &mtree_mvaid[i]  );       if(i==0 ) sync_mtreeName.push_back("mvaid");
      sync_fName.push_back("pho" + phoInd[i] + "_phoSCRawE" );    sync_fVal.push_back( &mtree_phoSCRawE[i]  );   if(i==0 ) sync_mtreeName.push_back("phoSCRawE"); 
      sync_fName.push_back("pho" + phoInd[i] + "_phoSmear" );     sync_fVal.push_back( &mtree_relSmearing[i]  ); if(i==0 ) sync_mtreeName.push_back("phoSmear");
      sync_fName.push_back("pho" + phoInd[i] + "_phoRegrE" );     sync_fVal.push_back( &mtree_phoRegrE[i]  );    if(i==0 ) sync_mtreeName.push_back("phoRegrE");
      sync_fName.push_back("pho" + phoInd[i] + "_S4Ratio" ); sync_fVal.push_back( &mtree_s4Ratio[i] ); if(i==0 ) sync_mtreeName.push_back("s4Ratio");
      
    }
    /*
      sync_fName.push_back("sceta"   + phoInd[i] ); sync_fVal.push_back( &mtree_sceta[i] );  if(i==0 ) sync_mtreeName.push_back("sceta");
      sync_fName.push_back("hoe"     + phoInd[i] ); sync_fVal.push_back( &mtree_HoE[i] );     if(i==0 ) sync_mtreeName.push_back("HoE");
      sync_fName.push_back("sigieie" + phoInd[i] ); sync_fVal.push_back( &mtree_coviEiE[i] );  if(i==0 ) sync_mtreeName.push_back("coviEiE");
      sync_fName.push_back("ecaliso" + phoInd[i] ); sync_fVal.push_back( &mtree_etCorrEcalIso[i] );  if(i==0 ) sync_mtreeName.push_back("etCorrEcalIso");
      sync_fName.push_back("hcaliso" + phoInd[i] ); sync_fVal.push_back( &mtree_etCorrHcalIso[i] ); if(i==0 ) sync_mtreeName.push_back("etCorrHcalIso");
      sync_fName.push_back("trckiso" + phoInd[i] ); sync_fVal.push_back( &mtree_etCorrTrkIso[i]  ); if(i==0 ) sync_mtreeName.push_back("etCorrTrkIso");
      sync_fName.push_back("chpfiso2"+ phoInd[i] ); sync_fVal.push_back( &mtree_pfChgIso02_selVtx[i]  ); if(i==0 ) sync_mtreeName.push_back("pfChgIso02_selVtx");
      sync_fName.push_back("phoid"   + phoInd[i] ); sync_fVal.push_back( &mtree_mvaid[i]  ); if(i==0 ) sync_mtreeName.push_back("mvaid");
      //sync_fName.push_back("phoeta"  + phoInd[i] ); sync_fVal.push_back( &mtree_eta[i]  ); if(i==0 ) sync_mtreeName.push_back("eta");
      sync_fName.push_back("pt"      + phoInd[i] ); sync_fVal.push_back( &mtree_ptg[i]  ); if(i==0 ) sync_mtreeName.push_back("ptg");
      
      //sync_fName.push_back("eerr"    + phoInd[i] ); sync_fVal.push_back( &mtree_relResOverE[i]  ); if(i==0 ) sync_mtreeName.push_back("relResOverE");
      sync_fName.push_back("eerrsmeared" + phoInd[i] ); sync_fVal.push_back( &mtree_relSmearing[i]  ); if(i==0 ) sync_mtreeName.push_back("relSmearing");
    */    
  }
  
  // ESEffSigmaRR
  sync_fName.push_back("mass");      sync_fVal.push_back( &mtree_mass );     sync_mtreeName.push_back("mass");
  sync_fName.push_back("diphoMVA" ); sync_fVal.push_back( &mtree_diphoMva ); sync_mtreeName.push_back("diphoMva");

  //-------------------------------------------
    
  if( object == "jets" ) {
    
    sync_fName.push_back("cat");        sync_fVal.push_back( &catFloat);          sync_mtreeName.push_back("cat"); //JM

    string jetInd[2] = {"1","2"};
    for( int ij = 0; ij < 2; ij++ ) {

      sync_fName.push_back("jet" + jetInd[ij] + "_pt"  );   sync_fVal.push_back( &dijet_ptj[ij] );       if(ij==0 ) sync_mtreeName.push_back("jetpt");
      sync_fName.push_back("jet" + jetInd[ij] + "_eta"  );  sync_fVal.push_back( &dijet_eta[ij] );       if(ij==0 ) sync_mtreeName.push_back("jeteta");
      sync_fName.push_back("jet" + jetInd[ij] + "_phi"  );  sync_fVal.push_back( &dijet_phi[ij] );       if(ij==0 ) sync_mtreeName.push_back("jetphi");
    }

    // JM Add ttH Jets to tree for Sync
    
    sync_fName.push_back("numJets");      sync_fVal.push_back( &jets_njetsFloat );     sync_mtreeName.push_back("njetstth");    
    sync_fName.push_back("numBJets");      sync_fVal.push_back( &jets_nbtagFloat );     sync_mtreeName.push_back("numBJets");
    sync_fName.push_back("met");      sync_fVal.push_back( &mtree_corMet );     sync_mtreeName.push_back("met");
    sync_fName.push_back("met_phi");      sync_fVal.push_back( &mtree_corMetPhi );     sync_mtreeName.push_back("met_phi");
    sync_fName.push_back("uncorrMet");      sync_fVal.push_back( &mtree_rawMet );     sync_mtreeName.push_back("uncorrMet");
    sync_fName.push_back("uncorrMet_phi");      sync_fVal.push_back( &mtree_rawMetPhi );     sync_mtreeName.push_back("uncorrMet_phi");
  }
  
  //-------------------------------------------

  if( object == "leptons" ) {

    string lInd[2] = {"1","2"};
    for( int il = 0; il < 2; il++ ) {

      sync_fName.push_back("ele" + lInd[il] + "_pt"  );    sync_fVal.push_back( &ele_pt[il] );        if(il==0 ) sync_mtreeName.push_back("elept");
      sync_fName.push_back("ele" + lInd[il] + "_eta"  );   sync_fVal.push_back( &ele_eta[il] );       if(il==0 ) sync_mtreeName.push_back("eleeta");
      sync_fName.push_back("ele" + lInd[il] + "_phi"  );   sync_fVal.push_back( &ele_phi[il] );       if(il==0 ) sync_mtreeName.push_back("elephi");

      sync_fName.push_back("mu" + lInd[il] + "_pt"  );     sync_fVal.push_back( &mu_pt[il] );         if(il==0 ) sync_mtreeName.push_back("mupt");
      sync_fName.push_back("mu" + lInd[il] + "_eta"  );    sync_fVal.push_back( &mu_eta[il] );        if(il==0 ) sync_mtreeName.push_back("mueta");
      sync_fName.push_back("mu" + lInd[il] + "_phi"  );    sync_fVal.push_back( &mu_phi[il] );        if(il==0 ) sync_mtreeName.push_back("muphi");
    }
  }

}





void MiniTree::setBranchesAddresses(void) {
  cout << " ------------- MiniTree set input branches ------------" << endl;
    //----------------------- MC true
  mtree_mcTrue->Branch("mH"       , &mc_wei      );
  mtree_mcTrue->Branch("cThetaStarTruth_CS", &mc_cThetaStar_CS); 

  mtree_mcTrue->SetBranchAddress("mc_H_Pt"        , &mc_H_Pt   );
  mtree_mcTrue->SetBranchAddress("mc_H_Eta"       , &mc_H_Eta  );
  mtree_mcTrue->SetBranchAddress("mc_H_Phi"       , &mc_H_Phi  );
  mtree_mcTrue->SetBranchAddress("mc_H_E"         , &mc_H_E    );

  mtree_mcTrue->SetBranchAddress("mc_partonOut_Pt1"   , &mc_partonOut_Pt1  );
  mtree_mcTrue->SetBranchAddress("mc_partonOut_Pt2"   , &mc_partonOut_Pt2  );
  mtree_mcTrue->SetBranchAddress("mc_partonOut_Eta1"  , &mc_partonOut_Eta1 );
  mtree_mcTrue->SetBranchAddress("mc_partonOut_Eta2"  , &mc_partonOut_Eta2 );
  mtree_mcTrue->SetBranchAddress("mc_partonOut_Phi1"  , &mc_partonOut_Phi1 );
  mtree_mcTrue->SetBranchAddress("mc_partonOut_Phi2"  , &mc_partonOut_Phi2 );
  mtree_mcTrue->SetBranchAddress("mc_partonOut_E1"    , &mc_partonOut_E1   );
  mtree_mcTrue->SetBranchAddress("mc_partonOut_E2"    , &mc_partonOut_E2   );

  mtree_mcTrue->SetBranchAddress("mc_partonOut_deltaPhi"   , &mc_partonOut_deltaPhi ); 
  mtree_mcTrue->SetBranchAddress("mc_partonOut_deltaPhiRF" , &mc_partonOut_deltaPhiRF ); 
  
 

  //----------------------- MC weights
  mtree_mainTree->SetBranchAddress("wei"      , &mc_wei      );
  mtree_mainTree->SetBranchAddress("wPU"      , &mc_wPU      );
  mtree_mainTree->SetBranchAddress("wHQT"     , &mc_wHQT     );
  mtree_mainTree->SetBranchAddress("wXsec"    , &mc_wXsec    );
  mtree_mainTree->SetBranchAddress("wNgen"    , &mc_wNgen    );
  mtree_mainTree->SetBranchAddress("wBSz"     , &mc_wBSz     );
  mtree_mainTree->SetBranchAddress("wVtxId"   , &mc_wVtxId   );
  mtree_mainTree->SetBranchAddress("wPhoEffi" , &mc_wPhoEffi );
  mtree_mainTree->SetBranchAddress("wTrigEffi", &mc_wTrigEffi);

  mtree_mainTree->SetBranchAddress("genPhoEMatched" , mc_genPhoEMatched   );

  //-----------------------
  mtree_mainTree->SetBranchAddress("iMainEvt", &mtree_ievt    );
  mtree_mainTree->SetBranchAddress("runNum"  , &mtree_runNum  );
  mtree_mainTree->SetBranchAddress("evtNum"  , &mtree_evtNum  );
  mtree_mainTree->SetBranchAddress("lumiSec" , &mtree_lumiSec );
  mtree_mainTree->SetBranchAddress("rho"     , &mtree_rho     );
  mtree_mainTree->SetBranchAddress("rho25"   , &mtree_rho25   );

  mtree_mainTree->SetBranchAddress("zVtx"    , &mtree_zVtx   );
  mtree_mainTree->SetBranchAddress("nVtx"    , &mtree_nVtx   );
  mtree_mainTree->SetBranchAddress("nVtxNoBS"  , &mtree_nVtxNoBS );
  mtree_mainTree->SetBranchAddress("iVtx1"   , &mtree_ivtx1  );
  mtree_mainTree->SetBranchAddress("iVtx2"   , &mtree_ivtx2  );
  mtree_mainTree->SetBranchAddress("iVtx3"   , &mtree_ivtx3  );
  mtree_mainTree->SetBranchAddress("vtxProb" , &mtree_vtxProb);
  mtree_mainTree->SetBranchAddress("vtxMva"  , &mtree_vtxMva );
  
  // skeleton
  //  mtree_mainTree->SetBranchAddress("x", &mtree_x ,"x/F"  );
  //--------------- diphoton event variables
  mtree_mainTree->SetBranchAddress("mass"      , &mtree_mass      );
  mtree_mainTree->SetBranchAddress("massDefVtx", &mtree_massDefVtx);
  mtree_mainTree->SetBranchAddress("pt"        , &mtree_pt        );
  mtree_mainTree->SetBranchAddress("y"         , &mtree_y         );
  mtree_mainTree->SetBranchAddress("piT"       , &mtree_piT       );
  mtree_mainTree->SetBranchAddress("pit1"      , &mtree_pit1      );
  mtree_mainTree->SetBranchAddress("pit2"      , &mtree_pit2      );
  mtree_mainTree->SetBranchAddress("massResoTot"     , &mtree_massResoTot      );
  mtree_mainTree->SetBranchAddress("massResoAng"     , &mtree_massResoAng      );
  mtree_mainTree->SetBranchAddress("massResoEng"     , &mtree_massResoEng      );
  mtree_mainTree->SetBranchAddress("massResoEngPosErr"     , &mtree_massResoEngPosErr      );
  mtree_mainTree->SetBranchAddress("massResoEngNegErr"     , &mtree_massResoEngNegErr      );
  mtree_mainTree->SetBranchAddress("massResoRightVtx", &mtree_massResoRightVtx );
  mtree_mainTree->SetBranchAddress("massResoRandVtx" , &mtree_massResoRandVtx  );
  mtree_mainTree->SetBranchAddress("diphoMva" , &mtree_diphoMva  );
  mtree_mainTree->SetBranchAddress("rawMet"   , &mtree_rawMet    );
  mtree_mainTree->SetBranchAddress("corMet"   , &mtree_corMet    );
  mtree_mainTree->SetBranchAddress("rawMetPhi", &mtree_rawMetPhi );
  mtree_mainTree->SetBranchAddress("corMetPhi", &mtree_corMetPhi );

  mtree_mainTree->SetBranchAddress( "catMva"  , &mtree_catMva );
  mtree_mainTree->SetBranchAddress( "tagCat"  , &mtree_tagCat );
  mtree_mainTree->SetBranchAddress( "untagCat", &mtree_untagCat );
  mtree_mainTree->SetBranchAddress("catBase", &mtree_catBase );
  mtree_mainTree->SetBranchAddress("lepTag" , &mtree_lepTag  );
  mtree_mainTree->SetBranchAddress("vbfTag" , &mtree_vbfTag  );
  mtree_mainTree->SetBranchAddress("metTag" , &mtree_metTag  );
  mtree_mainTree->SetBranchAddress("hvjjTag", &mtree_hvjjTag );
  mtree_mainTree->SetBranchAddress("lepCat" , &mtree_lepCat  );
  mtree_mainTree->SetBranchAddress("vbfCat" , &mtree_vbfCat  );
  mtree_mainTree->SetBranchAddress("cThetaLead_heli" ,  &mtree_cThetaLead_heli );
  mtree_mainTree->SetBranchAddress("cThetaTrail_heli",  &mtree_cThetaTrail_heli);
  mtree_mainTree->SetBranchAddress("cThetaStar_CS"   ,  &mtree_cThetaStar_CS   );


  // CF: add jet variables to main tree ****************************************************

  mtree_mainTree->SetBranchAddress( "nJet", &dijet_nJet );

  mtree_mainTree->SetBranchAddress( "jet1Pt"  , &dijet_ptj1 );
  mtree_mainTree->SetBranchAddress( "jet2Pt"  , &dijet_ptj2 );

  mtree_mainTree->SetBranchAddress( "jet1Eta", &dijet_etaj1 );
  mtree_mainTree->SetBranchAddress( "jet2Eta", &dijet_etaj2 );
  mtree_mainTree->SetBranchAddress( "jec1",    &dijet_jec1 );
  mtree_mainTree->SetBranchAddress( "jec2",    &dijet_jec2 );
  mtree_mainTree->SetBranchAddress( "jecUncrtPosj1", &dijet_jecUncrtPosj1 );
  mtree_mainTree->SetBranchAddress( "jecUncrtPosj2", &dijet_jecUncrtPosj2 );
  mtree_mainTree->SetBranchAddress( "jecUncrtNegj1", &dijet_jecUncrtNegj1 );
  mtree_mainTree->SetBranchAddress( "jecUncrtNegj2", &dijet_jecUncrtNegj2 );

  mtree_mainTree->SetBranchAddress( "etj", dijet_etj );
  mtree_mainTree->SetBranchAddress( "ptj", dijet_ptj );
  mtree_mainTree->SetBranchAddress( "jec", dijet_jec );
  mtree_mainTree->SetBranchAddress( "jecUncrtPos", dijet_jecUncrtPos );
  mtree_mainTree->SetBranchAddress( "jecUncrtNeg", dijet_jecUncrtNeg );
  mtree_mainTree->SetBranchAddress( "enj" , dijet_en  );
  mtree_mainTree->SetBranchAddress( "etaj", dijet_eta );
  mtree_mainTree->SetBranchAddress( "phij", dijet_phi );
  mtree_mainTree->SetBranchAddress( "rmsj", dijet_rms );
  mtree_mainTree->SetBranchAddress( "betaStar", dijet_betaStar );

  mtree_mainTree->SetBranchAddress( "dEtaJJ"    , &dijet_dEtaJJ     );
  mtree_mainTree->SetBranchAddress( "dPhiJJ"    , &dijet_dPhiJJ     );
  mtree_mainTree->SetBranchAddress( "zeppenfeld", &dijet_zeppenfeld );
  mtree_mainTree->SetBranchAddress( "mJJ"       , &dijet_mJJ        );
  mtree_mainTree->SetBranchAddress( "dphiGGJJ"  , &dijet_dphiGGJJ   );
 
  // ******************************************************************************************

  mtree_mainTree->SetBranchAddress("eg" , mtree_eg );
  mtree_mainTree->SetBranchAddress("ptg", mtree_ptg);
  mtree_mainTree->SetBranchAddress("eta", mtree_eta);
  mtree_mainTree->SetBranchAddress("phi", mtree_phi);
  mtree_mainTree->SetBranchAddress("cat", mtree_cat);
  mtree_mainTree->SetBranchAddress("sceta", mtree_sceta);
  mtree_mainTree->SetBranchAddress("relResOverE", mtree_relResOverE);
  mtree_mainTree->SetBranchAddress("relSmearing" , mtree_relSmearing);
  mtree_mainTree->SetBranchAddress("isElec", mtree_isElec);

  
/// phoid variables
  mtree_mainTree->SetBranchAddress("mvaid"  , mtree_mvaid   );
  mtree_mainTree->SetBranchAddress("s4Ratio", mtree_s4Ratio );
  mtree_mainTree->SetBranchAddress("phoIso1", mtree_phoIso1 );
  mtree_mainTree->SetBranchAddress("phoIso2", mtree_phoIso2 );
  mtree_mainTree->SetBranchAddress("phoIso3", mtree_phoIso3 );
  mtree_mainTree->SetBranchAddress("coviEiE", mtree_coviEiE );
  mtree_mainTree->SetBranchAddress("coviEiP", mtree_coviEiP );
  mtree_mainTree->SetBranchAddress("r9", mtree_r9 );
  mtree_mainTree->SetBranchAddress("HoE", mtree_HoE );
  mtree_mainTree->SetBranchAddress("drTrk", mtree_drTrk );
  mtree_mainTree->SetBranchAddress("etCorrEcalIso", mtree_etCorrEcalIso);
  mtree_mainTree->SetBranchAddress("etCorrHcalIso", mtree_etCorrHcalIso);
  mtree_mainTree->SetBranchAddress("etCorrTrkIso" , mtree_etCorrTrkIso );
  mtree_mainTree->SetBranchAddress("esEffSigmaRR" , mtree_esEffSigmaRR );
  mtree_mainTree->SetBranchAddress("mvaid1"  , &mtree_mvaid1   );
  mtree_mainTree->SetBranchAddress("mvaid2"  , &mtree_mvaid2   );
  mtree_mainTree->SetBranchAddress("eta1"    , &mtree_eta1   );
  mtree_mainTree->SetBranchAddress("eta2"    , &mtree_eta2   );
  mtree_mainTree->SetBranchAddress("r91"     , &mtree_r91   );
  mtree_mainTree->SetBranchAddress("r92"     , &mtree_r92   );

  /// photon PF isolation
  mtree_mainTree->SetBranchAddress("pfChgIso02_selVtx", mtree_pfChgIso02_selVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso03_selVtx", mtree_pfChgIso03_selVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso03_badVtx", mtree_pfChgIso03_badVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso03_defVtx", mtree_pfChgIso03_defVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso04_selVtx", mtree_pfChgIso04_selVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso04_badVtx", mtree_pfChgIso04_badVtx );
  mtree_mainTree->SetBranchAddress("pfChgIso04_defVtx", mtree_pfChgIso04_defVtx );
  mtree_mainTree->SetBranchAddress("pfPhoIso03", mtree_pfPhoIso03 );
  mtree_mainTree->SetBranchAddress("pfPhoIso04", mtree_pfPhoIso04 );

  mtree_mainTree->SetBranchAddress("mvaVbf", &dijet_mvaVbf );
  mtree_mainTree->SetBranchAddress("mvaVbf_sig", &dijet_mvaVbf_sig ); 
  mtree_mainTree->SetBranchAddress("mvaVbf_dipho", &dijet_mvaVbf_dipho );
  mtree_mainTree->SetBranchAddress("mvaVbf_gluglu", &dijet_mvaVbf_gluglu );
  mtree_mainTree->SetBranchAddress("mvaVbfCombDipho", &dijet_mvaVbfCombDipho );

  mtree_mainTree->SetBranchAddress( "dXsec_HJJ",       &mela_dXsec_HJJ );
  mtree_mainTree->SetBranchAddress( "dXsec_HJJVBF",    &mela_dXsec_HJJVBF );
  mtree_mainTree->SetBranchAddress( "dXsec_HJJ_PS",    &mela_dXsec_HJJ_PS );
  mtree_mainTree->SetBranchAddress( "dXsec_HJJVBF_PS", &mela_dXsec_HJJVBF_PS );

  mtree_mainTree->SetBranchAddress( "mela_VBFvsgg",    &mela_VBFvsgg );
  mtree_mainTree->SetBranchAddress( "mela_SMvsPS_VBF", &mela_SMvsPS_VBF );
  mtree_mainTree->SetBranchAddress( "mela_SMvsPS_gg",  &mela_SMvsPS_gg );

  // JM Add ttH Jets to tree for Sync ******************************************

  // initialize
  pjets_index=0;
  brjets_index=0;
  pjets_bjet_csv=0;
  brjets_bjet_csv=0;
  pjets_pt=0;
  brjets_pt=0;

  mtree_mainTree->SetBranchAddress("jets_nbtag", &jets_nbtag );
  mtree_mainTree->SetBranchAddress("jets_njets" , &jets_njets  );
  mtree_mainTree->SetBranchAddress("jets_bjet_csv" , &pjets_bjet_csv, &brjets_bjet_csv ); 
  mtree_mainTree->SetBranchAddress("jets_index" , &pjets_index, &brjets_index     );
  mtree_mainTree->SetBranchAddress("jets_pt" , &pjets_pt, &brjets_pt ); 
 
  
  // ********************************************************************************

  /// add sync variables
  /// sync only variables
  if( addSyncVariables ) {
    mtree_mainTree->SetBranchAddress("convInd", convInd );
    mtree_mainTree->SetBranchAddress("convNtk", convNtk );
    mtree_mainTree->SetBranchAddress("convPt" , convPt  );
    mtree_mainTree->SetBranchAddress("convValidVtx", convValidVtx );
    mtree_mainTree->SetBranchAddress("convChi2Prob", convChi2Prob );
    mtree_mainTree->SetBranchAddress("vId"  , vtxId    );
    mtree_mainTree->SetBranchAddress("vDz"  , vtxDz    );
    mtree_mainTree->SetBranchAddress("vMva" , vtxMva   );
    mtree_mainTree->SetBranchAddress("vAsym", &vtxAsym );
    mtree_mainTree->SetBranchAddress("vPtBal"     , &vtxPtBal     );
    mtree_mainTree->SetBranchAddress("vLogSumPt2" , &vtxLogSumPt2 );
    mtree_mainTree->SetBranchAddress("vPullToConv", &vtxPullToConv);
    mtree_mainTree->SetBranchAddress("vNConv"     , &vtxNConv     );
    mtree_mainTree->SetBranchAddress("enScale"    , mtree_enScale );    // CF
    mtree_mainTree->SetBranchAddress("phoSCRawE"  , mtree_phoSCRawE );  // CF
    mtree_mainTree->SetBranchAddress("phoRegrE"   , mtree_phoRegrE );   // CF

    mtree_mainTree->Branch("phoID_covIEtaIPhi"   , phoID_covIEtaIPhi );
    mtree_mainTree->Branch("phoID_S4"            , phoID_S4          );
    mtree_mainTree->Branch("phoID_ESEffSigmaRR"  , phoID_ESEffSigmaRR);
    mtree_mainTree->Branch("phoID_pfPhoIso03"    , phoID_pfPhoIso03  );
    mtree_mainTree->Branch("phoID_pfChIso03"     , phoID_pfChIso03   );
    mtree_mainTree->Branch("phoID_pfChIso03worst", phoID_pfChIso03worst);
    mtree_mainTree->Branch("phoID_E2x2", phoID_E2x2);
    mtree_mainTree->Branch("phoID_E5x5", phoID_E5x5);
  }
    


  // dilepton variables

  mtree_mainTree->SetBranchAddress( "nElec", &lept_nElec );
  mtree_mainTree->SetBranchAddress( "nMu",   &lept_nMu );

  mtree_mainTree->SetBranchAddress("elec_pt",  ele_pt );
  mtree_mainTree->SetBranchAddress("elec_eta", ele_eta );
  mtree_mainTree->SetBranchAddress("elec_phi", ele_phi );
  mtree_mainTree->SetBranchAddress("muon_pt",  mu_pt );
  mtree_mainTree->SetBranchAddress("muon_eta", mu_eta );
  mtree_mainTree->SetBranchAddress("muon_phi", mu_phi );
  
  
  

}



#endif
