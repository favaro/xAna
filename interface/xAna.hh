#ifndef xAna_hh__
#define xAna_hh__

#include "interface/category.hh"
#include "interface/weightManager.hh"
#include "interface/configReader.hh"
#include "interface/massResolutionCalculator.hh"
#include "interface/photonScaleAndSmearing.hh"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooAbsPdf.h>
#include <RooArgList.h>
#include <RooConstVar.h>

#include <TMVA/Reader.h>

//// external libraries
#include "VertexAnalysis/HggVertexAnalyzer.h"
#include "VertexAnalysis/HggVertexFromConversions.h"
#include "VertexAnalysis/VertexAlgoParameters.h"

#include "GBRLikelihood/HybridGBRForest.h"
#include "GBRLikelihood/RooHybridBDTAutoPdf.h"
#include "GBRLikelihood/RooDoubleCBFast.h"

#include "JetMETObjects/FactorizedJetCorrector.h"
#include "JetMETObjects/JetCorrectionUncertainty.h"
#include "JetMETObjects/JetCorrectorParameters.h"
#include "JetMETObjects/JetResolution.h"
#include "JetMETObjects/SimpleJetCorrectionUncertainty.h"
#include "JetMETObjects/SimpleJetCorrector.h"

// ME for MELA estimator:

#include "JHUGenMELA/TVar.hh"
#include "JHUGenMELA/TEvtProb.hh"
//#include "JHUGenMELA/math.h"

#include <fstream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>

class MiniTree;

const float vcicSTmva[17][4] = {
  { 3.8,     2.2,     1.77,     1.29},    // baseline rel combIso (good vtx)
  {11.7,     3.4,     3.90,     1.84},    // baseline rel combIso (bad vtx)
  { 3.5,     2.2,     2.30,     1.45},    // presel trkIso (good vtx)
  { 0.014,    0.014,  0.034,    0.034},   // presel sigma_ie_ie 
  { 0.082,   0.075,   0.075,    0.075},   // presel H/E       
  { 0.94,    0.36,    0.94,     0.32},    // baseline R9
  { 1,       0.062,   0.97,     0.97},    // baseline dR to trk
  { 1.5,     1.5,     1.5,      1.5},     // baseline pixel     
  { 50,      4,       50,       4},       // presel EtCorrEcalIso
  { 50,      4,       50,       4},       // presel EtCorrHcalIso
  { 50,      4,       50,       4},       // presel EtCorrTrkIso
  { 3,       3,       3,        3},       // presel PUCorrHcalEcal
  { 2.8,     2.8,     2.8,      2.8},     // presel AbsTrkIsoCIC
  { 4,       4,       4,        4},       // presel HollowConeTrkIsoDr03
  { 0.0106,  0.0097,  0.028,    0.027},   // baseline sigma_ie_ie
  { 0.082,   0.062,   0.065,    0.048},   // baseline H/E     
  { 4,       4,       4,        4}        // presel psChargedIso02
};



// ============= CiC4PF ==================

const float vcicST[7][4] = {
  { 6,     4.7,     5.6,     3.6},    //                             
  {10,     6.5,     5.6,     4.4},    //                                                 
  { 3.8,     2.5,    3.1,   2.2},    // Charged PF iso                                                          
  { 0.0108,  0.0102,  0.028,    0.028},   // sigma_ie_ie                                                                  
  { 0.124,   0.092,   0.142,    0.063},   // H/E       
  { 0.94,    0.298,    0.94,     0.24},    // R9
  //  { 1,       0.062,   0.97,     0.97},    // dR to trk
  { 1,       1,     1,     1},    // dR to trk
};





// ********************************************************************

class xAna {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
   
   
   // TH1F *hMCpu;
  float phoID_HoE,phoID_covIEtaIEta,
    phoID_tIso1abs,phoID_tIso2abs,phoID_tIso3abs,
    phoID_R9,phoID_absIsoEcal,
    phoID_absIsoHcal,phoID_NVertexes,phoID_ScEta,phoID_EtaWidth,phoID_PhiWidth;
  float phoID_covIEtaIPhi, phoID_S4, phoID_rho, 
    phoID_ESEffSigmaRR, phoID_pfPhoIso03, phoID_pfChIso03, phoID_pfChIso03worst ;
  float phoID_SCRawE;
  float Ddipho_masserr,Ddipho_masserrwrongvtx,Ddipho_vtxprob;
  float Ddipho_piT1,Ddipho_piT2,Ddipho_eta1,Ddipho_eta2,Ddipho_cosdPhi,Ddipho_id1,Ddipho_id2;
  
  float myVBFLeadJPt,myVBFSubJPt,myVBFdEta,myVBF_Mjj,myVBFZep,myVBF_Mgg,myVBFLeadJEta,myVBFSubJEta;
  float myVBFdPhi,myVBFDiPhoPtOverM,myVBFLeadPhoPtOverM,myVBFSubPhoPtOverM,myVBFdPhiTrunc;
  float myVBF_MVA,myVBFDIPHObdt;
  
  TMVA::Reader *phoID_2011[2];
  TMVA::Reader *phoID_2012[2];
  TMVA::Reader *phoID_mva[2]; 
  
  TMVA::Reader *DiscriDiPho_2011;
  TMVA::Reader *DiscriDiPho_2012;
  TMVA::Reader *Ddipho_mva;
  
  JetCorrectorParameters *ResJetPar;
  JetCorrectorParameters *L3JetPar ;
  JetCorrectorParameters *L2JetPar ;
  JetCorrectorParameters *L1JetPar ;  
  FactorizedJetCorrector *JetCorrector;
  bool setupJEC(void); 

  HybridGBRForest *_forestebPho;
  HybridGBRForest *_foresteePho;
  HybridGBRForest *_forestebEle;
  HybridGBRForest *_foresteeEle;
   
  RooRealVar *_mean;
  RooRealVar *_tgt;
  RooRealVar *_sigma;
  RooRealVar *_n1;
  RooRealVar *_n2;
  RooAbsReal *_meanlim;
  RooAbsReal *_sigmalim;
  RooAbsReal *_n1lim;
  RooAbsReal *_n2lim;
  RooAbsPdf *_pdf;
  std::vector<float> _GBRInvals;
   
  RooArgList _args;
  Bool_t _isInitialized;
   
   
  TMVA::Reader *DiscriVBF;
  TMVA::Reader *DiscriVBFCombDipho;
  bool DiscriVBF_UseDiPhoPt;
  bool DiscriVBF_UsePhoPt;
  bool DiscriVBF_useMvaSel;
  bool DiscriVBF_useCombMvaSel;

  bool vetoElec[2];
  
  std::string DiscriVBF_Method;
  std::vector<float> DiscriVBF_cat;
  std::vector<float> DiscriVBFCombDipho_cat;
  std::vector<float> DiscriVBF_sig_cat   ; /// not used for now
  std::vector<float> DiscriVBF_dipho_cat ; /// not used for now
  std::vector<float> DiscriVBF_gluglu_cat; /// not used for now
  
  
  float sumTrackPtInCone(int,int,float,float,float,float,float,float);
  float worstSumTrackPtInCone(int,int&,float,float,float,float,float,float);
  
  void SetConfigReader( const ConfigReader &config );  

   //// variables compute off
   float vertProb;
   float vertMVA;
  
   std::vector<float>  *phoRegrErr;
   std::vector<float>  *phoS4ratio;
   std::vector<float>  *phoRegrSmear;
   std::vector<int>    *eleNClus;   //<==== not in the Ntuples 
   std::vector<float>  *jetRawPtSmeared;
   std::vector<float>  *jetRawEnSmeared;

   std::vector<float>  *phoStdE;
   std::vector<float>  *phoEnergyScale;
   

  // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nHLT;
   Int_t           HLT[2500];   //[nHLT]
   Int_t           HLTIndex[70];
   Float_t         bspotPos[3];
   Int_t           nVtx;
   std::vector<float>   *vtx_x;
   std::vector<float>   *vtx_y;
   std::vector<float>   *vtx_z;
   Int_t           IsVtxGood;
   Int_t           nGoodVtx;
   Int_t           nVtxBS;
   std::vector<float>   *vtxbs_x;
   std::vector<float>   *vtxbs_y;
   std::vector<float>   *vtxbs_z;
   std::vector<float>   *vtxbsPtMod;
   std::vector<float>   *vtxbsSumPt2;
   std::vector<std::vector<int> > *vtxbsTkIndex;
   std::vector<std::vector<float> > *vtxbsTkWeight;
   Int_t           nTrk;
   std::vector<float>   *trkP_x;
   std::vector<float>   *trkP_y;
   std::vector<float>   *trkP_z;
   std::vector<float>   *trkVtx_x;
   std::vector<float>   *trkVtx_y;
   std::vector<float>   *trkVtx_z;
   std::vector<float>   *trkd0;
   std::vector<float>   *trkd0Err;
   std::vector<float>   *trkdz;
   std::vector<float>   *trkdzErr;
   std::vector<float>   *trkPtErr;
   std::vector<int>     *trkQuality;
   Int_t           nGoodTrk;
   Int_t           IsTracksGood;
   Float_t         pdf[7];
   Float_t         pthat;
   Float_t         processID;
   Int_t           nMC;
   std::vector<int>     *mcPID;
   std::vector<float>   *mcVtx_x;
   std::vector<float>   *mcVtx_y;
   std::vector<float>   *mcVtx_z;
   std::vector<float>   *mcPt;
   std::vector<float>   *mcMass;
   std::vector<float>   *mcEta;
   std::vector<float>   *mcPhi;
   std::vector<float>   *mcE;
   std::vector<float>   *mcEt;
   std::vector<int>     *mcGMomPID;
   std::vector<int>     *mcMomPID;
   std::vector<float>   *mcMomPt;
   std::vector<float>   *mcMomMass;
   std::vector<float>   *mcMomEta;
   std::vector<float>   *mcMomPhi;
   std::vector<int>     *mcIndex;
   std::vector<int>     *mcDecayType;
   std::vector<int>     *mcParentage;
   std::vector<int>     *mcStatus;
   Float_t         genMET;
   Float_t         genMETPhi;
   Int_t           nPUInfo;
   std::vector<int>     *nPU;
   std::vector<int>     *puBX;
   std::vector<float>   *puTrue;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         pfType01MET;
   Float_t         pfType01METPhi;
   Float_t         pfType01METsumEt;
   Float_t         pfType01METmEtSig;
   Float_t         pfType01METSig;
   Float_t         recoPfMET;
   Float_t         recoPfMETPhi;
   Float_t         recoPfMETsumEt;
   Float_t         recoPfMETmEtSig;
   Float_t         recoPfMETSig;
   Float_t         trkMETxPV;
   Float_t         trkMETyPV;
   Float_t         trkMETPhiPV;
   Float_t         trkMETPV;
   std::vector<float>   *trkMETx;
   std::vector<float>   *trkMETy;
   std::vector<float>   *trkMETPhi;
   std::vector<float>   *trkMET;
   Int_t           metFilters[10];
   Int_t           nEle;
   std::vector<unsigned long> *eleTrg;
   std::vector<int>     *eleClass;
   std::vector<int>     *eleIsEcalDriven;
   std::vector<int>     *eleCharge;
   std::vector<int>     *eleChargeConsistent;
   std::vector<float>   *eleEn;
   std::vector<float>   *eleEcalEn;
   std::vector<float>   *eleSCRawEn;
   std::vector<float>   *eleSCEn;
   std::vector<float>   *eleESEn;
   std::vector<float>   *elePt;
   std::vector<float>   *eleEta;
   std::vector<float>   *elePhi;
   std::vector<float>   *eleR9;
   std::vector<std::vector<float> > *eleEtaVtx;
   std::vector<std::vector<float> > *elePhiVtx;
   std::vector<std::vector<float> > *eleEtVtx;
   std::vector<float>   *eleSCEta;
   std::vector<float>   *eleSCPhi;
   std::vector<float>   *eleSCEtaWidth;
   std::vector<float>   *eleSCPhiWidth;
   std::vector<float>   *eleVtx_x;
   std::vector<float>   *eleVtx_y;
   std::vector<float>   *eleVtx_z;
   std::vector<float>   *eleD0;
   std::vector<float>   *eleDz;
   std::vector<float>   *eleD0GV;
   std::vector<float>   *eleDzGV;
   std::vector<std::vector<float> > *eleD0Vtx;
   std::vector<std::vector<float> > *eleDzVtx;
   std::vector<float>   *eleHoverE;
   std::vector<float>   *eleHoverE12;
   std::vector<float>   *eleEoverP;
   std::vector<float>   *elePin;
   std::vector<float>   *elePout;
   std::vector<float>   *eleTrkMomErr;
   std::vector<float>   *eleBrem;
   std::vector<float>   *eledEtaAtVtx;
   std::vector<float>   *eledPhiAtVtx;
   std::vector<float>   *eleSigmaIEtaIEta;
   std::vector<float>   *eleSigmaIEtaIPhi;
   std::vector<float>   *eleSigmaIPhiIPhi;
   std::vector<float>   *eleEmax;
   std::vector<float>   *eleE2ndMax;
   std::vector<float>   *eleETop;
   std::vector<float>   *eleEBottom;
   std::vector<float>   *eleELeft;
   std::vector<float>   *eleERight;
   std::vector<float>   *eleE1x5;
   std::vector<float>   *eleE3x3;
   std::vector<float>   *eleE5x5;
   std::vector<float>   *eleE2x5Max;
   std::vector<float>   *eleE2x5Top;
   std::vector<float>   *eleE2x5Bottom;
   std::vector<float>   *eleE2x5Left;
   std::vector<float>   *eleE2x5Right;
   std::vector<float>   *eleSeedEta;
   std::vector<float>   *eleSeedE;
   std::vector<float>   *eleSeedPhi;
   std::vector<float>   *eleCrysEta;
   std::vector<float>   *eleCrysPhi;
   std::vector<int>     *eleCrysIEta;
   std::vector<int>     *eleCrysIPhi;
   std::vector<float>   *eleRegrE;
   std::vector<float>   *eleRegrEerr;
   std::vector<float>   *elePhoRegrE;
   std::vector<float>   *elePhoRegrEerr;
   std::vector<float>   *eleSeedTime;
   std::vector<int>     *eleRecoFlag;
   std::vector<int>     *elePos;
   std::vector<int>     *eleGenIndex;
   std::vector<int>     *eleGenGMomPID;
   std::vector<int>     *eleGenMomPID;
   std::vector<float>   *eleGenMomPt;
   std::vector<float>   *eleIsoTrkDR03;
   std::vector<float>   *eleIsoEcalDR03;
   std::vector<float>   *eleIsoHcalDR03;
   std::vector<float>   *eleIsoHcalDR0312;
   std::vector<float>   *eleIsoTrkDR04;
   std::vector<float>   *eleIsoEcalDR04;
   std::vector<float>   *eleIsoHcalDR04;
   std::vector<float>   *eleIsoHcalDR0412;
   std::vector<float>   *eleModIsoTrk;
   std::vector<float>   *eleModIsoEcal;
   std::vector<float>   *eleModIsoHcal;
   std::vector<int>     *eleMissHits;
   std::vector<float>   *eleConvDist;
   std::vector<float>   *eleConvDcot;
   std::vector<int>     *eleConvVtxFit;
   std::vector<float>   *eleIP3D;
   std::vector<float>   *eleIP3DErr;
   std::vector<float>   *eleIDMVANonTrig;
   std::vector<float>   *eleIDMVATrig;
   std::vector<float>   *elePFChIso03;
   std::vector<float>   *elePFPhoIso03;
   std::vector<float>   *elePFNeuIso03;
   std::vector<float>   *elePFChIso04;
   std::vector<float>   *elePFPhoIso04;
   std::vector<float>   *elePFNeuIso04;
   std::vector<float>   *eleESEffSigmaRR_x;
   std::vector<float>   *eleESEffSigmaRR_y;
   std::vector<float>   *eleESEffSigmaRR_z;
   Int_t           nPho;
   std::vector<unsigned long> *phoTrg;
   std::vector<unsigned long> *phoTrgFilter;
   std::vector<bool>    *phoIsPhoton;
   std::vector<float>   *phoSCPos_x;
   std::vector<float>   *phoSCPos_y;
   std::vector<float>   *phoSCPos_z;
   std::vector<float>   *phoCaloPos_x;
   std::vector<float>   *phoCaloPos_y;
   std::vector<float>   *phoCaloPos_z;
   std::vector<float>   *phoE;
   std::vector<float>   *phoEt;
   std::vector<float>   *phoEta;
   std::vector<float>   *phoVtx_x;
   std::vector<float>   *phoVtx_y;
   std::vector<float>   *phoVtx_z;
   std::vector<float>   *phoPhi;
   std::vector<std::vector<float> > *phoEtVtx;
   std::vector<std::vector<float> > *phoEtaVtx;
   std::vector<std::vector<float> > *phoPhiVtx;
   std::vector<float>   *phoR9;
   std::vector<int>     *phoNClus;
   std::vector<float>   *phoTrkIsoHollowDR03;
   std::vector<float>   *phoEcalIsoDR03;
   std::vector<float>   *phoHcalIsoDR03;
   std::vector<float>   *phoHcalIsoDR0312;
   std::vector<float>   *phoTrkIsoHollowDR04;
   std::vector<float>   *phoCiCdRtoTrk;
   std::vector<float>   *phoEcalIsoDR04;
   std::vector<float>   *phoHcalIsoDR04;
   std::vector<float>   *phoHcalIsoDR0412;
   std::vector<float>   *phoHoverE;
   std::vector<float>   *phoHoverE12;
   std::vector<int>     *phoEleVeto;
   std::vector<float>   *phoSigmaIEtaIEta;
   std::vector<float>   *phoSigmaIEtaIPhi;
   std::vector<float>   *phoSigmaIPhiIPhi;
   std::vector<float>   *phoCiCPF4phopfIso03;
   std::vector<float>   *phoCiCPF4phopfIso04;
   std::vector<std::vector<float> > *phoCiCPF4chgpfIso02;
   std::vector<std::vector<float> > *phoCiCPF4chgpfIso03;
   std::vector<std::vector<float> > *phoCiCPF4chgpfIso04;
   std::vector<float>   *phoEmax;
   std::vector<float>   *phoETop;
   std::vector<float>   *phoEBottom;
   std::vector<float>   *phoELeft;
   std::vector<float>   *phoERight;
   std::vector<float>   *phoE2ndMax;
   std::vector<float>   *phoE3x3;
   std::vector<float>   *phoE3x1;
   std::vector<float>   *phoE1x3;
   std::vector<float>   *phoE5x5;
   std::vector<float>   *phoE1x5;
   std::vector<float>   *phoE2x2;
   std::vector<float>   *phoE2x5Max;
   std::vector<float>   *phoE2x5Top;
   std::vector<float>   *phoE2x5Bottom;
   std::vector<float>   *phoE2x5Left;
   std::vector<float>   *phoE2x5Right;
   std::vector<float>   *phoSeedE;
   std::vector<float>   *phoSeedEta;
   std::vector<float>   *phoSeedPhi;
   std::vector<float>   *phoCrysEta;
   std::vector<float>   *phoCrysPhi;
   std::vector<int>     *phoCrysIEta;
   std::vector<int>     *phoCrysIPhi;
   std::vector<float>   *phoPFChIso;
   std::vector<float>   *phoPFPhoIso;
   std::vector<float>   *phoPFNeuIso;
   std::vector<float>   *phoSCRChIso;
   std::vector<float>   *phoSCRPhoIso;
   std::vector<float>   *phoSCRNeuIso;
   std::vector<float>   *phoSCRChIso04;
   std::vector<float>   *phoSCRPhoIso04;
   std::vector<float>   *phoSCRNeuIso04;
   std::vector<float>   *phoRandConeChIso;
   std::vector<float>   *phoRandConePhoIso;
   std::vector<float>   *phoRandConeNeuIso;
   std::vector<float>   *phoRandConeChIso04;
   std::vector<float>   *phoRandConePhoIso04;
   std::vector<float>   *phoRandConeNeuIso04;
   std::vector<float>   *phoRegrE;
   std::vector<float>   *phoRegrEerr;
   std::vector<float>   *phoSeedTime;
   std::vector<int>     *phoSeedDetId1;
   std::vector<int>     *phoSeedDetId2;
   std::vector<float>   *phoLICTD;
   std::vector<int>     *phoRecoFlag;
   std::vector<int>     *phoPos;
   std::vector<int>     *phoGenIndex;
   std::vector<int>     *phoGenGMomPID;
   std::vector<int>     *phoGenMomPID;
   std::vector<float>   *phoGenMomPt;
   std::vector<float>   *phoSCE;
   std::vector<float>   *phoSCRawE;
   std::vector<float>   *phoESEn;
   std::vector<float>   *phoSCEt;
   std::vector<float>   *phoSCEta;
   std::vector<float>   *phoSCPhi;
   std::vector<float>   *phoSCEtaWidth;
   std::vector<float>   *phoSCPhiWidth;
   std::vector<float>   *phoSCBrem;
   std::vector<int>     *phoOverlap;
   std::vector<int>     *phohasPixelSeed;
   std::vector<int>     *pho_hasConvPf;
   std::vector<int>     *pho_hasSLConvPf;
   std::vector<float>   *pho_pfconvVtxZ;
   std::vector<float>   *pho_pfconvVtxZErr;
   std::vector<int>     *pho_nSLConv;
   std::vector<std::vector<float> > *pho_pfSLConvPos_x;
   std::vector<std::vector<float> > *pho_pfSLConvPos_y;
   std::vector<std::vector<float> > *pho_pfSLConvPos_z;
   std::vector<std::vector<float> > *pho_pfSLConvVtxZ;
   std::vector<int>     *phoIsConv;
   std::vector<int>     *phoNConv;
   std::vector<float>   *phoConvInvMass;
   std::vector<float>   *phoConvCotTheta;
   std::vector<float>   *phoConvEoverP;
   std::vector<float>   *phoConvZofPVfromTrks;
   std::vector<float>   *phoConvMinDist;
   std::vector<float>   *phoConvdPhiAtVtx;
   std::vector<float>   *phoConvdPhiAtCalo;
   std::vector<float>   *phoConvdEtaAtCalo;
   std::vector<float>   *phoConvTrkd0_x;
   std::vector<float>   *phoConvTrkd0_y;
   std::vector<float>   *phoConvTrkPin_x;
   std::vector<float>   *phoConvTrkPin_y;
   std::vector<float>   *phoConvTrkPout_x;
   std::vector<float>   *phoConvTrkPout_y;
   std::vector<float>   *phoConvTrkdz_x;
   std::vector<float>   *phoConvTrkdz_y;
   std::vector<float>   *phoConvTrkdzErr_x;
   std::vector<float>   *phoConvTrkdzErr_y;
   std::vector<float>   *phoConvChi2;
   std::vector<float>   *phoConvChi2Prob;
   std::vector<int>     *phoConvNTrks;
   std::vector<float>   *phoConvCharge1;
   std::vector<float>   *phoConvCharge2;
   std::vector<int>     *phoConvValidVtx;
   std::vector<float>   *phoConvLikeLihood;
   std::vector<float>   *phoConvP4_0;
   std::vector<float>   *phoConvP4_1;
   std::vector<float>   *phoConvP4_2;
   std::vector<float>   *phoConvP4_3;
   std::vector<float>   *phoConvVtx_x;
   std::vector<float>   *phoConvVtx_y;
   std::vector<float>   *phoConvVtx_z;
   std::vector<float>   *phoConvVtxErr_x;
   std::vector<float>   *phoConvVtxErr_y;
   std::vector<float>   *phoConvVtxErr_z;
   std::vector<float>   *phoConvPairMomentum_x;
   std::vector<float>   *phoConvPairMomentum_y;
   std::vector<float>   *phoConvPairMomentum_z;
   std::vector<float>   *phoConvRefittedMomentum_x;
   std::vector<float>   *phoConvRefittedMomentum_y;
   std::vector<float>   *phoConvRefittedMomentum_z;
   std::vector<int>     *SingleLegConv;
   std::vector<std::vector<float> > *phoPFConvVtx_x;
   std::vector<std::vector<float> > *phoPFConvVtx_y;
   std::vector<std::vector<float> > *phoPFConvVtx_z;
   std::vector<std::vector<float> > *phoPFConvMom_x;
   std::vector<std::vector<float> > *phoPFConvMom_y;
   std::vector<std::vector<float> > *phoPFConvMom_z;
   std::vector<float>   *phoESEffSigmaRR_x;
   std::vector<float>   *phoESEffSigmaRR_y;
   std::vector<float>   *phoESEffSigmaRR_z;
   Int_t           nMu;
   std::vector<unsigned long> *muTrg;
   std::vector<float>   *muEta;
   std::vector<float>   *muPhi;
   std::vector<int>     *muCharge;
   std::vector<float>   *muPt;
   std::vector<float>   *muPz;
   std::vector<float>   *muVtx_x;
   std::vector<float>   *muVtx_y;
   std::vector<float>   *muVtx_z;
   std::vector<float>   *muVtxGlb_x;
   std::vector<float>   *muVtxGlb_y;
   std::vector<float>   *muVtxGlb_z;
   std::vector<int>     *muGenIndex;
   std::vector<float>   *mucktPt;
   std::vector<float>   *mucktPtErr;
   std::vector<float>   *mucktEta;
   std::vector<float>   *mucktPhi;
   std::vector<float>   *mucktdxy;
   std::vector<float>   *mucktdz;
   std::vector<float>   *muIsoTrk;
   std::vector<float>   *muIsoCalo;
   std::vector<float>   *muIsoEcal;
   std::vector<float>   *muIsoHcal;
   std::vector<float>   *muChi2NDF;
   std::vector<float>   *muInnerChi2NDF;
   std::vector<float>   *muPFIsoR04_CH;
   std::vector<float>   *muPFIsoR04_NH;
   std::vector<float>   *muPFIsoR04_Pho;
   std::vector<float>   *muPFIsoR04_PU;
   std::vector<float>   *muPFIsoR04_CPart;
   std::vector<float>   *muPFIsoR04_NHHT;
   std::vector<float>   *muPFIsoR04_PhoHT;
   std::vector<float>   *muPFIsoR03_CH;
   std::vector<float>   *muPFIsoR03_NH;
   std::vector<float>   *muPFIsoR03_Pho;
   std::vector<float>   *muPFIsoR03_PU;
   std::vector<float>   *muPFIsoR03_CPart;
   std::vector<float>   *muPFIsoR03_NHHT;
   std::vector<float>   *muPFIsoR03_PhoHT;
   std::vector<int>     *muType;
   std::vector<float>   *muD0;
   std::vector<float>   *muDz;
   std::vector<float>   *muD0GV;
   std::vector<float>   *muDzGV;
   std::vector<std::vector<float> > *muD0Vtx;
   std::vector<std::vector<float> > *muDzVtx;
   std::vector<float>   *muInnerD0;
   std::vector<float>   *muInnerDz;
   std::vector<float>   *muInnerD0GV;
   std::vector<float>   *muInnerDzGV;
   std::vector<float>   *muInnerPt;
   std::vector<float>   *muInnerPtErr;
   std::vector<int>     *muNumberOfValidTrkLayers;
   std::vector<int>     *muNumberOfValidTrkHits;
   std::vector<int>     *muNumberOfValidPixelLayers;
   std::vector<int>     *muNumberOfValidPixelHits;
   std::vector<int>     *muNumberOfValidMuonHits;
   std::vector<int>     *muStations;
   std::vector<int>     *muChambers;
   std::vector<float>   *muIP3D;
   std::vector<float>   *muIP3DErr;
   Int_t           nTau;
   std::vector<bool>    *tauDecayModeFinding;
   std::vector<bool>    *tauAgainstElectronLooseMVA3;
   std::vector<bool>    *tauAgainstElectronMediumMVA3;
   std::vector<bool>    *tauAgainstElectronTightMVA3;
   std::vector<bool>    *tauAgainstElectronVTightMVA3;
   std::vector<bool>    *tauAgainstElectronDeadECAL;
   std::vector<bool>    *tauAgainstMuonLoose2;
   std::vector<bool>    *tauAgainstMuonMedium2;
   std::vector<bool>    *tauAgainstMuonTight2;
   std::vector<bool>    *tauCombinedIsolationDeltaBetaCorrRaw3Hits;
   std::vector<bool>    *tauLooseCombinedIsolationDeltaBetaCorr3Hits;
   std::vector<bool>    *tauMediumCombinedIsolationDeltaBetaCorr3Hits;
   std::vector<bool>    *tauTightCombinedIsolationDeltaBetaCorr3Hits;
   std::vector<float>   *tauEta;
   std::vector<float>   *tauPhi;
   std::vector<float>   *tauPt;
   std::vector<float>   *tauEt;
   std::vector<float>   *tauCharge;
   std::vector<int>     *tauDecayMode;
   std::vector<float>   *tauEMFraction;
   std::vector<float>   *tauHCAL3x3OverPLead;
   std::vector<float>   *tauHCALMaxOverPLead;
   std::vector<float>   *tauHCALTotOverPLead;
   std::vector<float>   *tauIsolationPFChargedHadrCandsPtSum;
   std::vector<float>   *tauIsolationPFGammaCandsEtSum;
   std::vector<float>   *tauLeadPFChargedHadrCandsignedSipt;
   std::vector<bool>    *tauLeadChargedHadronExists;
   std::vector<float>   *tauLeadChargedHadronEta;
   std::vector<float>   *tauLeadChargedHadronPhi;
   std::vector<float>   *tauLeadChargedHadronPt;
   Float_t         rho25;
   Float_t         rho25_neu;
   Float_t         rho25_muPFiso;
   Float_t         rho25_elePFiso;
   Float_t         rho2011;
   Float_t         rho2012;
   Float_t         QGTag_MLP;
   Float_t         QGTag_likelihood;
   Int_t           nCA8Jet;
   std::vector<float>   *CA8JetPt;
   std::vector<float>   *CA8JetEta;
   std::vector<float>   *CA8JetPhi;
   std::vector<float>   *CA8JetMass;
   std::vector<float>   *CA8JetArea;
   std::vector<float>   *CA8Jet_tau1;
   std::vector<float>   *CA8Jet_tau2;
   std::vector<float>   *CA8Jet_tau3;
   std::vector<float>   *CA8JetCHF;
   std::vector<float>   *CA8JetNHF;
   std::vector<float>   *CA8JetCEF;
   std::vector<float>   *CA8JetNEF;
   std::vector<int>     *CA8JetNCH;
   std::vector<int>     *CA8Jetnconstituents;
   std::vector<float>   *CA8prunedJetMass;
   std::vector<int>     *CA8prunedJet_nSubJets;
   std::vector<std::vector<float> > *CA8prunedJet_SubjetPt;
   std::vector<std::vector<float> > *CA8prunedJet_SubjetEta;
   std::vector<std::vector<float> > *CA8prunedJet_SubjetPhi;
   std::vector<std::vector<float> > *CA8prunedJet_SubjetMass;
   Int_t           nJet;
   std::vector<unsigned long> *jetTrg;
   std::vector<float>   *jetEn;
   std::vector<float>   *jetPt;
   std::vector<float>   *jetEta;
   std::vector<float>   *jetPhi;
   std::vector<float>   *jetCharge;
   std::vector<float>   *jetEt;
   std::vector<float>   *jetRawPt;
   std::vector<float>   *jetRawEn;
   std::vector<float>   *jetArea;
   std::vector<float>   *jetCHF;
   std::vector<float>   *jetNHF;
   std::vector<float>   *jetCEF;
   std::vector<float>   *jetNEF;
   std::vector<int>     *jetNCH;
   std::vector<float>   *jetHFHAE;
   std::vector<float>   *jetHFEME;
   std::vector<int>     *jetNConstituents;
   std::vector<float>   *jetCombinedSecondaryVtxBJetTags;
   std::vector<float>   *jetCombinedSecondaryVtxMVABJetTags;
   std::vector<float>   *jetJetProbabilityBJetTags;
   std::vector<float>   *jetJetBProbabilityBJetTags;
   std::vector<std::vector<float> > *jetBetaStar;
   std::vector<bool>    *jetPFLooseId;
   std::vector<float>   *jetDRMean;
   std::vector<float>   *jetDR2Mean;
   std::vector<float>   *jetDZ;
   std::vector<float>   *jetFrac01;
   std::vector<float>   *jetFrac02;
   std::vector<float>   *jetFrac03;
   std::vector<float>   *jetFrac04;
   std::vector<float>   *jetFrac05;
   std::vector<float>   *jetFrac06;
   std::vector<float>   *jetFrac07;
   std::vector<float>   *jetBeta;
   std::vector<float>   *jetBetaStarCMG;
   std::vector<float>   *jetBetaStarClassic;
   std::vector<std::vector<float> > *jetBetaExt;
   std::vector<std::vector<float> > *jetBetaStarCMGExt;
   std::vector<std::vector<float> > *jetBetaStarClassicExt;
   std::vector<float>   *jetNNeutrals;
   std::vector<float>   *jetNCharged;
   std::vector<std::vector<float> > *jetMVAs;
   std::vector<std::vector<int> > *jetWPLevels;
   std::vector<std::vector<float> > *jetMVAsExt_simple;
   std::vector<std::vector<int> > *jetWPLevelsExt_simple;
   std::vector<std::vector<float> > *jetMVAsExt_full;
   std::vector<std::vector<int> > *jetWPLevelsExt_full;
   std::vector<std::vector<float> > *jetMVAsExt_cutBased;
   std::vector<std::vector<int> > *jetWPLevelsExt_cutBased;
   std::vector<std::vector<float> > *jetMVAsExt_philv1;
   std::vector<std::vector<int> > *jetWPLevelsExt_philv1;
   std::vector<float>   *jetMt;
   std::vector<float>   *jetJECUnc;
   std::vector<float>   *jetLeadTrackPt;
   std::vector<float>   *jetVtxPt;
   std::vector<float>   *jetVtxMass;
   std::vector<float>   *jetVtx3dL;
   std::vector<float>   *jetVtx3deL;
   std::vector<float>   *jetSoftLeptPt;
   std::vector<float>   *jetSoftLeptPtRel;
   std::vector<float>   *jetSoftLeptdR;
   std::vector<float>   *jetSoftLeptIdlooseMu;
   std::vector<float>   *jetSoftLeptIdEle95;
   std::vector<float>   *jetDPhiMETJet;
   std::vector<float>   *jetPuJetIdL;
   std::vector<float>   *jetPuJetIdM;
   std::vector<float>   *jetPuJetIdT;
   std::vector<int>     *jetPartonID;
   std::vector<int>     *jetGenJetIndex;
   std::vector<float>   *jetGenJetEn;
   std::vector<float>   *jetGenJetPt;
   std::vector<float>   *jetGenJetEta;
   std::vector<float>   *jetGenJetPhi;
   std::vector<int>     *jetGenPartonID;
   std::vector<float>   *jetGenEn;
   std::vector<float>   *jetGenPt;
   std::vector<float>   *jetGenEta;
   std::vector<float>   *jetGenPhi;
   std::vector<int>     *jetGenPartonMomID;
   Int_t           nLowPtJet;
   std::vector<float>   *jetLowPtEn;
   std::vector<float>   *jetLowPtPt;
   std::vector<float>   *jetLowPtEta;
   std::vector<float>   *jetLowPtPhi;
   std::vector<float>   *jetLowPtCharge;
   std::vector<float>   *jetLowPtEt;
   std::vector<float>   *jetLowPtRawPt;
   std::vector<float>   *jetLowPtRawEn;
   std::vector<float>   *jetLowPtArea;
   std::vector<int>     *jetLowPtPartonID;
   std::vector<float>   *jetLowPtGenJetEn;
   std::vector<float>   *jetLowPtGenJetPt;
   std::vector<float>   *jetLowPtGenJetEta;
   std::vector<float>   *jetLowPtGenJetPhi;
   std::vector<int>     *jetLowPtGenPartonID;
   std::vector<float>   *jetLowPtGenEn;
   std::vector<float>   *jetLowPtGenPt;
   std::vector<float>   *jetLowPtGenEta;
   std::vector<float>   *jetLowPtGenPhi;
   Int_t           nConv;
   std::vector<float>   *convP4_x;
   std::vector<float>   *convP4_y;
   std::vector<float>   *convP4_z;
   std::vector<float>   *convP4_E;
   std::vector<float>   *convVtx_x;
   std::vector<float>   *convVtx_y;
   std::vector<float>   *convVtx_z;
   std::vector<float>   *convVtxErr_x;
   std::vector<float>   *convVtxErr_y;
   std::vector<float>   *convVtxErr_z;
   std::vector<float>   *convPairMomentum_x;
   std::vector<float>   *convPairMomentum_y;
   std::vector<float>   *convPairMomentum_z;
   std::vector<float>   *convRefittedMomentum_x;
   std::vector<float>   *convRefittedMomentum_y;
   std::vector<float>   *convRefittedMomentum_z;
   std::vector<int>     *convNTracks;
   std::vector<float>   *convPairInvMass;
   std::vector<float>   *convPairCotThetaSep;
   std::vector<float>   *convEoverP;
   std::vector<float>   *convDistOfMinApproach;
   std::vector<float>   *convDPhiTrksAtVtx;
   std::vector<float>   *convDPhiTrksAtEcal;
   std::vector<float>   *convDEtaTrksAtEcal;
   std::vector<float>   *convDxy;
   std::vector<float>   *convDz;
   std::vector<float>   *convLxy;
   std::vector<float>   *convLz;
   std::vector<float>   *convZofPrimVtxFromTrks;
   std::vector<int>     *convNHitsBeforeVtx_0;
   std::vector<int>     *convNHitsBeforeVtx_1;
   std::vector<int>     *convNSharedHits;
   std::vector<int>     *convValidVtx;
   std::vector<float>   *convMVALikelihood;
   std::vector<float>   *convChi2;
   std::vector<float>   *convChi2Probability;
   std::vector<float>   *convTk1Dz;
   std::vector<float>   *convTk2Dz;
   std::vector<float>   *convTk1DzErr;
   std::vector<float>   *convTk2DzErr;
   std::vector<int>     *convCh1Ch2;
   std::vector<float>   *convTk1D0;
   std::vector<float>   *convTk1Pout;
   std::vector<float>   *convTk1Pin;
   std::vector<float>   *convTk2D0;
   std::vector<float>   *convTk2Pout;
   std::vector<float>   *convTk2Pin;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_HLT;   //!
   TBranch        *b_HLTIndex;   //!
   TBranch        *b_bspotPos;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_IsVtxGood;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_nVtxBS;   //!
   TBranch        *b_vtxbs_x;   //!
   TBranch        *b_vtxbs_y;   //!
   TBranch        *b_vtxbs_z;   //!
   TBranch        *b_vtxbsPtMod;   //!
   TBranch        *b_vtxbsSumPt2;   //!
   TBranch        *b_vtxbsTkIndex;   //!
   TBranch        *b_vtxbsTkWeight;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkP_x;   //!
   TBranch        *b_trkP_y;   //!
   TBranch        *b_trkP_z;   //!
   TBranch        *b_trkVtx_x;   //!
   TBranch        *b_trkVtx_y;   //!
   TBranch        *b_trkVtx_z;   //!
   TBranch        *b_trkd0;   //!
   TBranch        *b_trkd0Err;   //!
   TBranch        *b_trkdz;   //!
   TBranch        *b_trkdzErr;   //!
   TBranch        *b_trkPtErr;   //!
   TBranch        *b_trkQuality;   //!
   TBranch        *b_nGoodTrk;   //!
   TBranch        *b_IsTracksGood;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx_x;   //!
   TBranch        *b_mcVtx_y;   //!
   TBranch        *b_mcVtx_z;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_mcDecayType;   //!
   TBranch        *b_mcParentage;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_pfType01MET;   //!
   TBranch        *b_pfType01METPhi;   //!
   TBranch        *b_pfType01METsumEt;   //!
   TBranch        *b_pfType01METmEtSig;   //!
   TBranch        *b_pfType01METSig;   //!
   TBranch        *b_recoPfMET;   //!
   TBranch        *b_recoPfMETPhi;   //!
   TBranch        *b_recoPfMETsumEt;   //!
   TBranch        *b_recoPfMETmEtSig;   //!
   TBranch        *b_recoPfMETSig;   //!
   TBranch        *b_trkMETxPV;   //!
   TBranch        *b_trkMETyPV;   //!
   TBranch        *b_trkMETPhiPV;   //!
   TBranch        *b_trkMETPV;   //!
   TBranch        *b_trkMETx;   //!
   TBranch        *b_trkMETy;   //!
   TBranch        *b_trkMETPhi;   //!
   TBranch        *b_trkMET;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleTrg;   //!
   TBranch        *b_eleClass;   //!
   TBranch        *b_eleIsEcalDriven;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleEcalEn;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleESEn;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleEtaVtx;   //!
   TBranch        *b_elePhiVtx;   //!
   TBranch        *b_eleEtVtx;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleVtx_x;   //!
   TBranch        *b_eleVtx_y;   //!
   TBranch        *b_eleVtx_z;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleD0GV;   //!
   TBranch        *b_eleDzGV;   //!
   TBranch        *b_eleD0Vtx;   //!
   TBranch        *b_eleDzVtx;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleHoverE12;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_elePin;   //!
   TBranch        *b_elePout;   //!
   TBranch        *b_eleTrkMomErr;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIEtaIPhi;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleEmax;   //!
   TBranch        *b_eleE2ndMax;   //!
   TBranch        *b_eleETop;   //!
   TBranch        *b_eleEBottom;   //!
   TBranch        *b_eleELeft;   //!
   TBranch        *b_eleERight;   //!
   TBranch        *b_eleE1x5;   //!
   TBranch        *b_eleE3x3;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleE2x5Max;   //!
   TBranch        *b_eleE2x5Top;   //!
   TBranch        *b_eleE2x5Bottom;   //!
   TBranch        *b_eleE2x5Left;   //!
   TBranch        *b_eleE2x5Right;   //!
   TBranch        *b_eleSeedEta;   //!
   TBranch        *b_eleSeedE;   //!
   TBranch        *b_eleSeedPhi;   //!
   TBranch        *b_eleCrysEta;   //!
   TBranch        *b_eleCrysPhi;   //!
   TBranch        *b_eleCrysIEta;   //!
   TBranch        *b_eleCrysIPhi;   //!
   TBranch        *b_eleRegrE;   //!
   TBranch        *b_eleRegrEerr;   //!
   TBranch        *b_elePhoRegrE;   //!
   TBranch        *b_elePhoRegrEerr;   //!
   TBranch        *b_eleSeedTime;   //!
   TBranch        *b_eleRecoFlag;   //!
   TBranch        *b_elePos;   //!
   TBranch        *b_eleGenIndex;   //!
   TBranch        *b_eleGenGMomPID;   //!
   TBranch        *b_eleGenMomPID;   //!
   TBranch        *b_eleGenMomPt;   //!
   TBranch        *b_eleIsoTrkDR03;   //!
   TBranch        *b_eleIsoEcalDR03;   //!
   TBranch        *b_eleIsoHcalDR03;   //!
   TBranch        *b_eleIsoHcalDR0312;   //!
   TBranch        *b_eleIsoTrkDR04;   //!
   TBranch        *b_eleIsoEcalDR04;   //!
   TBranch        *b_eleIsoHcalDR04;   //!
   TBranch        *b_eleIsoHcalDR0412;   //!
   TBranch        *b_eleModIsoTrk;   //!
   TBranch        *b_eleModIsoEcal;   //!
   TBranch        *b_eleModIsoHcal;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleConvDist;   //!
   TBranch        *b_eleConvDcot;   //!
   TBranch        *b_eleConvVtxFit;   //!
   TBranch        *b_eleIP3D;   //!
   TBranch        *b_eleIP3DErr;   //!
   TBranch        *b_eleIDMVANonTrig;   //!
   TBranch        *b_eleIDMVATrig;   //!
   TBranch        *b_elePFChIso03;   //!
   TBranch        *b_elePFPhoIso03;   //!
   TBranch        *b_elePFNeuIso03;   //!
   TBranch        *b_elePFChIso04;   //!
   TBranch        *b_elePFPhoIso04;   //!
   TBranch        *b_elePFNeuIso04;   //!
   TBranch        *b_eleESEffSigmaRR_x;   //!
   TBranch        *b_eleESEffSigmaRR_y;   //!
   TBranch        *b_eleESEffSigmaRR_z;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoTrg;   //!
   TBranch        *b_phoTrgFilter;   //!
   TBranch        *b_phoIsPhoton;   //!
   TBranch        *b_phoSCPos_x;   //!
   TBranch        *b_phoSCPos_y;   //!
   TBranch        *b_phoSCPos_z;   //!
   TBranch        *b_phoCaloPos_x;   //!
   TBranch        *b_phoCaloPos_y;   //!
   TBranch        *b_phoCaloPos_z;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoVtx_x;   //!
   TBranch        *b_phoVtx_y;   //!
   TBranch        *b_phoVtx_z;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEtVtx;   //!
   TBranch        *b_phoEtaVtx;   //!
   TBranch        *b_phoPhiVtx;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoNClus;   //!
   TBranch        *b_phoTrkIsoHollowDR03;   //!
   TBranch        *b_phoEcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR0312;   //!
   TBranch        *b_phoTrkIsoHollowDR04;   //!
   TBranch        *b_phoCiCdRtoTrk;   //!
   TBranch        *b_phoEcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoDR0412;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoHoverE12;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoSigmaIEtaIPhi;   //!
   TBranch        *b_phoSigmaIPhiIPhi;   //!
   TBranch        *b_phoCiCPF4phopfIso03;   //!
   TBranch        *b_phoCiCPF4phopfIso04;   //!
   TBranch        *b_phoCiCPF4chgpfIso02;   //!
   TBranch        *b_phoCiCPF4chgpfIso03;   //!
   TBranch        *b_phoCiCPF4chgpfIso04;   //!
   TBranch        *b_phoEmax;   //!
   TBranch        *b_phoETop;   //!
   TBranch        *b_phoEBottom;   //!
   TBranch        *b_phoELeft;   //!
   TBranch        *b_phoERight;   //!
   TBranch        *b_phoE2ndMax;   //!
   TBranch        *b_phoE3x3;   //!
   TBranch        *b_phoE3x1;   //!
   TBranch        *b_phoE1x3;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoE1x5;   //!
   TBranch        *b_phoE2x2;   //!
   TBranch        *b_phoE2x5Max;   //!
   TBranch        *b_phoE2x5Top;   //!
   TBranch        *b_phoE2x5Bottom;   //!
   TBranch        *b_phoE2x5Left;   //!
   TBranch        *b_phoE2x5Right;   //!
   TBranch        *b_phoSeedE;   //!
   TBranch        *b_phoSeedEta;   //!
   TBranch        *b_phoSeedPhi;   //!
   TBranch        *b_phoCrysEta;   //!
   TBranch        *b_phoCrysPhi;   //!
   TBranch        *b_phoCrysIEta;   //!
   TBranch        *b_phoCrysIPhi;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoSCRChIso;   //!
   TBranch        *b_phoSCRPhoIso;   //!
   TBranch        *b_phoSCRNeuIso;   //!
   TBranch        *b_phoSCRChIso04;   //!
   TBranch        *b_phoSCRPhoIso04;   //!
   TBranch        *b_phoSCRNeuIso04;   //!
   TBranch        *b_phoRandConeChIso;   //!
   TBranch        *b_phoRandConePhoIso;   //!
   TBranch        *b_phoRandConeNeuIso;   //!
   TBranch        *b_phoRandConeChIso04;   //!
   TBranch        *b_phoRandConePhoIso04;   //!
   TBranch        *b_phoRandConeNeuIso04;   //!
   TBranch        *b_phoRegrE;   //!
   TBranch        *b_phoRegrEerr;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedDetId1;   //!
   TBranch        *b_phoSeedDetId2;   //!
   TBranch        *b_phoLICTD;   //!
   TBranch        *b_phoRecoFlag;   //!
   TBranch        *b_phoPos;   //!
   TBranch        *b_phoGenIndex;   //!
   TBranch        *b_phoGenGMomPID;   //!
   TBranch        *b_phoGenMomPID;   //!
   TBranch        *b_phoGenMomPt;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoSCEt;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phoOverlap;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_pho_hasConvPf;   //!
   TBranch        *b_pho_hasSLConvPf;   //!
   TBranch        *b_pho_pfconvVtxZ;   //!
   TBranch        *b_pho_pfconvVtxZErr;   //!
   TBranch        *b_pho_nSLConv;   //!
   TBranch        *b_pho_pfSLConvPos_x;   //!
   TBranch        *b_pho_pfSLConvPos_y;   //!
   TBranch        *b_pho_pfSLConvPos_z;   //!
   TBranch        *b_pho_pfSLConvVtxZ;   //!
   TBranch        *b_phoIsConv;   //!
   TBranch        *b_phoNConv;   //!
   TBranch        *b_phoConvInvMass;   //!
   TBranch        *b_phoConvCotTheta;   //!
   TBranch        *b_phoConvEoverP;   //!
   TBranch        *b_phoConvZofPVfromTrks;   //!
   TBranch        *b_phoConvMinDist;   //!
   TBranch        *b_phoConvdPhiAtVtx;   //!
   TBranch        *b_phoConvdPhiAtCalo;   //!
   TBranch        *b_phoConvdEtaAtCalo;   //!
   TBranch        *b_phoConvTrkd0_x;   //!
   TBranch        *b_phoConvTrkd0_y;   //!
   TBranch        *b_phoConvTrkPin_x;   //!
   TBranch        *b_phoConvTrkPin_y;   //!
   TBranch        *b_phoConvTrkPout_x;   //!
   TBranch        *b_phoConvTrkPout_y;   //!
   TBranch        *b_phoConvTrkdz_x;   //!
   TBranch        *b_phoConvTrkdz_y;   //!
   TBranch        *b_phoConvTrkdzErr_x;   //!
   TBranch        *b_phoConvTrkdzErr_y;   //!
   TBranch        *b_phoConvChi2;   //!
   TBranch        *b_phoConvChi2Prob;   //!
   TBranch        *b_phoConvNTrks;   //!
   TBranch        *b_phoConvCharge1;   //!
   TBranch        *b_phoConvCharge2;   //!
   TBranch        *b_phoConvValidVtx;   //!
   TBranch        *b_phoConvLikeLihood;   //!
   TBranch        *b_phoConvP4_0;   //!
   TBranch        *b_phoConvP4_1;   //!
   TBranch        *b_phoConvP4_2;   //!
   TBranch        *b_phoConvP4_3;   //!
   TBranch        *b_phoConvVtx_x;   //!
   TBranch        *b_phoConvVtx_y;   //!
   TBranch        *b_phoConvVtx_z;   //!
   TBranch        *b_phoConvVtxErr_x;   //!
   TBranch        *b_phoConvVtxErr_y;   //!
   TBranch        *b_phoConvVtxErr_z;   //!
   TBranch        *b_phoConvPairMomentum_x;   //!
   TBranch        *b_phoConvPairMomentum_y;   //!
   TBranch        *b_phoConvPairMomentum_z;   //!
   TBranch        *b_phoConvRefittedMomentum_x;   //!
   TBranch        *b_phoConvRefittedMomentum_y;   //!
   TBranch        *b_phoConvRefittedMomentum_z;   //!
   TBranch        *b_SingleLegConv;   //!
   TBranch        *b_phoPFConvVtx_x;   //!
   TBranch        *b_phoPFConvVtx_y;   //!
   TBranch        *b_phoPFConvVtx_z;   //!
   TBranch        *b_phoPFConvMom_x;   //!
   TBranch        *b_phoPFConvMom_y;   //!
   TBranch        *b_phoPFConvMom_z;   //!
   TBranch        *b_phoESEffSigmaRR_x;   //!
   TBranch        *b_phoESEffSigmaRR_y;   //!
   TBranch        *b_phoESEffSigmaRR_z;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muTrg;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muPz;   //!
   TBranch        *b_muVtx_x;   //!
   TBranch        *b_muVtx_y;   //!
   TBranch        *b_muVtx_z;   //!
   TBranch        *b_muVtxGlb_x;   //!
   TBranch        *b_muVtxGlb_y;   //!
   TBranch        *b_muVtxGlb_z;   //!
   TBranch        *b_muGenIndex;   //!
   TBranch        *b_mucktPt;   //!
   TBranch        *b_mucktPtErr;   //!
   TBranch        *b_mucktEta;   //!
   TBranch        *b_mucktPhi;   //!
   TBranch        *b_mucktdxy;   //!
   TBranch        *b_mucktdz;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muIsoCalo;   //!
   TBranch        *b_muIsoEcal;   //!
   TBranch        *b_muIsoHcal;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerChi2NDF;   //!
   TBranch        *b_muPFIsoR04_CH;   //!
   TBranch        *b_muPFIsoR04_NH;   //!
   TBranch        *b_muPFIsoR04_Pho;   //!
   TBranch        *b_muPFIsoR04_PU;   //!
   TBranch        *b_muPFIsoR04_CPart;   //!
   TBranch        *b_muPFIsoR04_NHHT;   //!
   TBranch        *b_muPFIsoR04_PhoHT;   //!
   TBranch        *b_muPFIsoR03_CH;   //!
   TBranch        *b_muPFIsoR03_NH;   //!
   TBranch        *b_muPFIsoR03_Pho;   //!
   TBranch        *b_muPFIsoR03_PU;   //!
   TBranch        *b_muPFIsoR03_CPart;   //!
   TBranch        *b_muPFIsoR03_NHHT;   //!
   TBranch        *b_muPFIsoR03_PhoHT;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muD0GV;   //!
   TBranch        *b_muDzGV;   //!
   TBranch        *b_muD0Vtx;   //!
   TBranch        *b_muDzVtx;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muInnerD0GV;   //!
   TBranch        *b_muInnerDzGV;   //!
   TBranch        *b_muInnerPt;   //!
   TBranch        *b_muInnerPtErr;   //!
   TBranch        *b_muNumberOfValidTrkLayers;   //!
   TBranch        *b_muNumberOfValidTrkHits;   //!
   TBranch        *b_muNumberOfValidPixelLayers;   //!
   TBranch        *b_muNumberOfValidPixelHits;   //!
   TBranch        *b_muNumberOfValidMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muChambers;   //!
   TBranch        *b_muIP3D;   //!
   TBranch        *b_muIP3DErr;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_tauDecayModeFinding;   //!
   TBranch        *b_tauAgainstElectronLooseMVA3;   //!
   TBranch        *b_tauAgainstElectronMediumMVA3;   //!
   TBranch        *b_tauAgainstElectronTightMVA3;   //!
   TBranch        *b_tauAgainstElectronVTightMVA3;   //!
   TBranch        *b_tauAgainstElectronDeadECAL;   //!
   TBranch        *b_tauAgainstMuonLoose2;   //!
   TBranch        *b_tauAgainstMuonMedium2;   //!
   TBranch        *b_tauAgainstMuonTight2;   //!
   TBranch        *b_tauCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_tauLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauEta;   //!
   TBranch        *b_tauPhi;   //!
   TBranch        *b_tauPt;   //!
   TBranch        *b_tauEt;   //!
   TBranch        *b_tauCharge;   //!
   TBranch        *b_tauDecayMode;   //!
   TBranch        *b_tauEMFraction;   //!
   TBranch        *b_tauHCAL3x3OverPLead;   //!
   TBranch        *b_tauHCALMaxOverPLead;   //!
   TBranch        *b_tauHCALTotOverPLead;   //!
   TBranch        *b_tauIsolationPFChargedHadrCandsPtSum;   //!
   TBranch        *b_tauIsolationPFGammaCandsEtSum;   //!
   TBranch        *b_tauLeadPFChargedHadrCandsignedSipt;   //!
   TBranch        *b_tauLeadChargedHadronExists;   //!
   TBranch        *b_tauLeadChargedHadronEta;   //!
   TBranch        *b_tauLeadChargedHadronPhi;   //!
   TBranch        *b_tauLeadChargedHadronPt;   //!
   TBranch        *b_rho25;   //!
   TBranch        *b_rho25_neu;   //!
   TBranch        *b_rho25_muPFiso;   //!
   TBranch        *b_rho25_elePFiso;   //!
   TBranch        *b_rho2011;   //!
   TBranch        *b_rho2012;   //!
   TBranch        *b_QGTag_MLP;   //!
   TBranch        *b_QGTag_likelihood;   //!
   TBranch        *b_nCA8Jet;   //!
   TBranch        *b_CA8JetPt;   //!
   TBranch        *b_CA8JetEta;   //!
   TBranch        *b_CA8JetPhi;   //!
   TBranch        *b_CA8JetMass;   //!
   TBranch        *b_CA8JetArea;   //!
   TBranch        *b_CA8Jet_tau1;   //!
   TBranch        *b_CA8Jet_tau2;   //!
   TBranch        *b_CA8Jet_tau3;   //!
   TBranch        *b_CA8JetCHF;   //!
   TBranch        *b_CA8JetNHF;   //!
   TBranch        *b_CA8JetCEF;   //!
   TBranch        *b_CA8JetNEF;   //!
   TBranch        *b_CA8JetNCH;   //!
   TBranch        *b_CA8Jetnconstituents;   //!
   TBranch        *b_CA8prunedJetMass;   //!
   TBranch        *b_CA8prunedJet_nSubJets;   //!
   TBranch        *b_CA8prunedJet_SubjetPt;   //!
   TBranch        *b_CA8prunedJet_SubjetEta;   //!
   TBranch        *b_CA8prunedJet_SubjetPhi;   //!
   TBranch        *b_CA8prunedJet_SubjetMass;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetTrg;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetCharge;   //!
   TBranch        *b_jetEt;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetHFHAE;   //!
   TBranch        *b_jetHFEME;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_jetCombinedSecondaryVtxBJetTags;   //!
   TBranch        *b_jetCombinedSecondaryVtxMVABJetTags;   //!
   TBranch        *b_jetJetProbabilityBJetTags;   //!
   TBranch        *b_jetJetBProbabilityBJetTags;   //!
   TBranch        *b_jetBetaStar;   //!
   TBranch        *b_jetPFLooseId;   //!
   TBranch        *b_jetDRMean;   //!
   TBranch        *b_jetDR2Mean;   //!
   TBranch        *b_jetDZ;   //!
   TBranch        *b_jetFrac01;   //!
   TBranch        *b_jetFrac02;   //!
   TBranch        *b_jetFrac03;   //!
   TBranch        *b_jetFrac04;   //!
   TBranch        *b_jetFrac05;   //!
   TBranch        *b_jetFrac06;   //!
   TBranch        *b_jetFrac07;   //!
   TBranch        *b_jetBeta;   //!
   TBranch        *b_jetBetaStarCMG;   //!
   TBranch        *b_jetBetaStarClassic;   //!
   TBranch        *b_jetBetaExt;   //!
   TBranch        *b_jetBetaStarCMGExt;   //!
   TBranch        *b_jetBetaStarClassicExt;   //!
   TBranch        *b_jetNNeutrals;   //!
   TBranch        *b_jetNCharged;   //!
   TBranch        *b_jetMVAs;   //!
   TBranch        *b_jetWPLevels;   //!
   TBranch        *b_jetMVAsExt_simple;   //!
   TBranch        *b_jetWPLevelsExt_simple;   //!
   TBranch        *b_jetMVAsExt_full;   //!
   TBranch        *b_jetWPLevelsExt_full;   //!
   TBranch        *b_jetMVAsExt_cutBased;   //!
   TBranch        *b_jetWPLevelsExt_cutBased;   //!
   TBranch        *b_jetMVAsExt_philv1;   //!
   TBranch        *b_jetWPLevelsExt_philv1;   //!
   TBranch        *b_jetMt;   //!
   TBranch        *b_jetJECUnc;   //!
   TBranch        *b_jetLeadTrackPt;   //!
   TBranch        *b_jetVtxPt;   //!
   TBranch        *b_jetVtxMass;   //!
   TBranch        *b_jetVtx3dL;   //!
   TBranch        *b_jetVtx3deL;   //!
   TBranch        *b_jetSoftLeptPt;   //!
   TBranch        *b_jetSoftLeptPtRel;   //!
   TBranch        *b_jetSoftLeptdR;   //!
   TBranch        *b_jetSoftLeptIdlooseMu;   //!
   TBranch        *b_jetSoftLeptIdEle95;   //!
   TBranch        *b_jetDPhiMETJet;   //!
   TBranch        *b_jetPuJetIdL;   //!
   TBranch        *b_jetPuJetIdM;   //!
   TBranch        *b_jetPuJetIdT;   //!
   TBranch        *b_jetPartonID;   //!
   TBranch        *b_jetGenJetIndex;   //!
   TBranch        *b_jetGenJetEn;   //!
   TBranch        *b_jetGenJetPt;   //!
   TBranch        *b_jetGenJetEta;   //!
   TBranch        *b_jetGenJetPhi;   //!
   TBranch        *b_jetGenPartonID;   //!
   TBranch        *b_jetGenEn;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenPhi;   //!
   TBranch        *b_jetGenPartonMomID;   //!
   TBranch        *b_nLowPtJet;   //!
   TBranch        *b_jetLowPtEn;   //!
   TBranch        *b_jetLowPtPt;   //!
   TBranch        *b_jetLowPtEta;   //!
   TBranch        *b_jetLowPtPhi;   //!
   TBranch        *b_jetLowPtCharge;   //!
   TBranch        *b_jetLowPtEt;   //!
   TBranch        *b_jetLowPtRawPt;   //!
   TBranch        *b_jetLowPtRawEn;   //!
   TBranch        *b_jetLowPtArea;   //!
   TBranch        *b_jetLowPtPartonID;   //!
   TBranch        *b_jetLowPtGenJetEn;   //!
   TBranch        *b_jetLowPtGenJetPt;   //!
   TBranch        *b_jetLowPtGenJetEta;   //!
   TBranch        *b_jetLowPtGenJetPhi;   //!
   TBranch        *b_jetLowPtGenPartonID;   //!
   TBranch        *b_jetLowPtGenEn;   //!
   TBranch        *b_jetLowPtGenPt;   //!
   TBranch        *b_jetLowPtGenEta;   //!
   TBranch        *b_jetLowPtGenPhi;   //!
   TBranch        *b_nConv;   //!
   TBranch        *b_convP4_x;   //!
   TBranch        *b_convP4_y;   //!
   TBranch        *b_convP4_z;   //!
   TBranch        *b_convP4_E;   //!
   TBranch        *b_convVtx_x;   //!
   TBranch        *b_convVtx_y;   //!
   TBranch        *b_convVtx_z;   //!
   TBranch        *b_convVtxErr_x;   //!
   TBranch        *b_convVtxErr_y;   //!
   TBranch        *b_convVtxErr_z;   //!
   TBranch        *b_convPairMomentum_x;   //!
   TBranch        *b_convPairMomentum_y;   //!
   TBranch        *b_convPairMomentum_z;   //!
   TBranch        *b_convRefittedMomentum_x;   //!
   TBranch        *b_convRefittedMomentum_y;   //!
   TBranch        *b_convRefittedMomentum_z;   //!
   TBranch        *b_convNTracks;   //!
   TBranch        *b_convPairInvMass;   //!
   TBranch        *b_convPairCotThetaSep;   //!
   TBranch        *b_convEoverP;   //!
   TBranch        *b_convDistOfMinApproach;   //!
   TBranch        *b_convDPhiTrksAtVtx;   //!
   TBranch        *b_convDPhiTrksAtEcal;   //!
   TBranch        *b_convDEtaTrksAtEcal;   //!
   TBranch        *b_convDxy;   //!
   TBranch        *b_convDz;   //!
   TBranch        *b_convLxy;   //!
   TBranch        *b_convLz;   //!
   TBranch        *b_convZofPrimVtxFromTrks;   //!
   TBranch        *b_convNHitsBeforeVtx_0;   //!
   TBranch        *b_convNHitsBeforeVtx_1;   //!
   TBranch        *b_convNSharedHits;   //!
   TBranch        *b_convValidVtx;   //!
   TBranch        *b_convMVALikelihood;   //!
   TBranch        *b_convChi2;   //!
   TBranch        *b_convChi2Probability;   //!
   TBranch        *b_convTk1Dz;   //!
   TBranch        *b_convTk2Dz;   //!
   TBranch        *b_convTk1DzErr;   //!
   TBranch        *b_convTk2DzErr;   //!
   TBranch        *b_convCh1Ch2;   //!
   TBranch        *b_convTk1D0;   //!
   TBranch        *b_convTk1Pout;   //!
   TBranch        *b_convTk1Pin;   //!
   TBranch        *b_convTk2D0;   //!
   TBranch        *b_convTk2Pout;   //!
   TBranch        *b_convTk2Pin;   //!


   //// adding variables for MC Truth
   std::vector<float>   *mcHardPt;
   std::vector<float>   *mcHardEta;
   std::vector<float>   *mcHardPhi;
   std::vector<float>   *mcHardM;
   std::vector<int>     *mcHardPID;
   std::vector<int>     *mcHardFun;
  
  //CF 
   std::vector<float>   *mcHardOutPt;
   std::vector<float>   *mcHardOutP;
   std::vector<float>   *mcHardOutE; 
   std::vector<float>   *mcHardOutEta;
   std::vector<float>   *mcHardOutPhi;
   std::vector<float>   *mcHardOutPdgId;

   TBranch        *b_mcHardPt;   //!
   TBranch        *b_mcHardEta;   //!
   TBranch        *b_mcHardPhi;   //!
   TBranch        *b_mcHardM;   //!
   TBranch        *b_mcHardPID;   //!
   TBranch        *b_mcHardFun;   //!
  
  //CF
   TBranch        *b_mcHardOutPt;
   TBranch        *b_mcHardOutP;
   TBranch        *b_mcHardOutE; 
   TBranch        *b_mcHardOutEta;
   TBranch        *b_mcHardOutPhi;
   TBranch        *b_mcHardOutPdgId;

   xAna(TTree *tree=0, bool isItData=false);
   ~xAna();
   //Int_t    Cut(Long64_t entry);
   Int_t    GetEntry(Long64_t entry);
   Long64_t LoadTree(Long64_t entry);
   void     Init(TTree *tree, bool isItData);


   /// FC: add some function for reweighting
  void SetWeightManager( WeightManager *weight_manager) { _weight_manager = weight_manager ;}
  bool PassTriggerSelection2011(void);
  bool PassTriggerSelection2012(void);

   //CF: move photon energy corrections to separate function
   //void ApplyEnergyCorrections( bool DoOverSmearing );
   void  ApplyEnergyCorrections( bool DoOverSmearing, int i ); // CF: per-photon function

   /// FC: MVA function
   void CorrLikelihood(int ipho, float &ecorr, float &esigma , bool forPho);
   float PhoID_MVA( int i, int ivtx ) ;
   float DiPhoID_MVA( const TLorentzVector &glead, const TLorentzVector &gtrail, const TLorentzVector &hcand, 
		      const MassResolution & massResoCalc, float vtxProb,
		      const float &idlead, const float &idtrail );

   void rescalePhoIdVar( int i );
   void rescalePhoIdVarBaseline( int i );

   /// FC: electron and conversion matching
   bool isMatchedToAPromptElec( int ipho, float &dRclosestElec );
   bool isAGoodConv( int ipho );
   bool isInGAP_EB( int ipho);
   
   void Setup_MVA(void);
   
   /// FC: met selection
   bool MetTagSelection( const TLorentzVector &g1, const TLorentzVector &g2, const std::vector<int> &jindex );

   /// FC: MC truth functions
   bool isZgamma(void);
   void fillMCtruthInfo_HiggsSignal(void);
   void recoToTruthMatching( HggEvtCandidate& event );
   
   int     HggTreeWriteLoop(const char* filename, int ijob );
   //void     ZggTreeWriteLoop(const char* filename, Int_t mode = 0);
   Bool_t   Notify();
   void     Show(Long64_t entry = -1);

   Int_t photonID(int ipho, int ivtx, float corEt_ipho = -1, bool invEleVeto = false );
   Int_t cicPhotonID(int ipho, int ivtx, float corEt_ipho = -1 );
   Int_t cic4PFPhotonID(int ipho, int ivtx, float corEt_ipho = -1, bool invEleVeto = false );
   Int_t mvaPhotonID(int ipho, int ivtx, float corEt_ipho = -1, bool invEleVeto = false );
   void  fillPhotonVariablesToMiniTree(int ipho, int ivtx, int index );

//   Int_t    romePhotonID(int ipho);
   Int_t    phocat(int i);
   Int_t    phocatCIC(int i);
   Double_t deltaPhi(Double_t phi1, Double_t phi2);
   Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
   void     InitHists();

   Int_t   matchPhotonToConversion( int lpho);
   Float_t   etaTransformation(  float EtaParticle , float vertex);
   PhotonInfo  fillPhotonInfo(  int iPho1 );
   //int    getSelectedVertex(int iPho1, int iPho2, bool useAllConvs=true);
   //  int    getSelectedVertex(int iPho1, int iPho2, bool useAllConvs=true);
   std::vector<int>    getSelectedVertex(int iPho1, int iPho2,  bool useMva=true);

   //Int_t photonCat(int i) ;
   //bool preselectPhoton(int i) ;
  bool preselectPhoton(int i, float Et) ;
  TVector3 getCorPhotonTVector3(int i, int selVtx);
 
  float computeMCweights(const HggEvtCandidate &hcand, const bool isprompt);

   Double_t getHqTWeight(Float_t HiggsMCMass, Float_t HiggsMCPt);
   
   // CF
   std::vector<HggEvtCandidate> buildHggCandidates( bool DoPreselection );

  bool selectTwoPhotons( int &i, int &j, int selVtx, 
			 TLorentzVector &g1, TLorentzVector &g2 );



   bool selectTwoPhotonsAwayFromLeptons( int i, int j,  
					 std::vector<int> elecIndex ,  std::vector<int> muonIndex);//JM

   Int_t*   selectElectrons2011(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd);
   Int_t*   selectMuons2011(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd);//JM
 
   std::vector<int> selectElectrons(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd);
   std::vector<int> selectElectronsAwayFromTwoPhotons(int i, int j, std::vector<int> elecIndex ); //JM
   std::vector<int> selectMuons(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd);
   std::vector<int> selectElectronsHCP2012(  int &iLepVtx, double ptcut=20 );
   std::vector<int> selectMuonsHCP2012( int &iLepVtx, double ptcut=20 );

  
   std::vector<int>   selectJets(Int_t iPho1, Int_t iPho2, Int_t nVtxJetID, Double_t wei, Int_t selVtx);   
   std::vector<int>   selectJetsJEC(Int_t iPho1, Int_t iPho2, Int_t nVtxJetID, Int_t selVtx,std::vector<int>& allGoodJets);   
   std::vector<int>   selectJetsTTH( std::vector<int> goodJetsIndex ,  std::vector<int> elecIndex, std::vector<int> muonIndex, int& nBTag ); //JM
   std::vector<int>   selectJetsVHLep( std::vector<int> goodJetsIndex ,  std::vector<int> elecIndex, std::vector<int> muonIndex ); //JM
   std::vector<int>   selectJetsVHHad( std::vector<int> goodJetsIndex  ); //JM

   void dijetSelection( const TLorentzVector &corPho1,const TLorentzVector &corPho2,					 
			const std::vector<int> &goodJetsIndex,
			int ivtx,
			int &vbftag, int &hstratag, int &catjet );


   void           applyJEC4SoftJets();
   void           METSmearCorrection(Int_t iPho1, Int_t iPho2);
   double         getMETJetSmearFactor(Float_t metJetEta);
   double         ErrEt( double Et, double Eta);
   void           METPhiShiftCorrection();
   void           METScaleCorrection(Int_t iPho1, Int_t iPho2);
   TLorentzVector corMET( Int_t iPho1, Int_t iPho2);
   double         getMETJetScaleFactor(Float_t metJetEta);

   double         evaluateDiPhotonMVA( int ilead, int itrail,  const TLorentzVector & gKlead, const TLorentzVector & gKtrail, MassResolution & massResoCalc );

   void           fillJetVariablesToMinitree( const HggEvtCandidate hggevent );
   float          combineMVAforVBF( const HggEvtCandidate hggevent );

   void           evaluateMELA( const HggEvtCandidate hggevent );

   //------------------------------------------------

   // CF FOR CATEGORIES:

   
   int                   convertToOfficialCat( const TString catname, const int subcat );
   void                  createCategories( std::vector<category>& tag, std::vector<category>& untag );
   bool                  isEvtInCategory( HggEvtCandidate hggevent, category cat );
   void                  analyzeHggCandidates( HggEvtCandidate& hggevt, std::vector<HggEvtCandidate>& hggevents, std::vector<category> allCategories, MassResolution massresoc );
   HggEvtCandidate       sortHggCandidates( std::vector<HggEvtCandidate>& hggevents, std::vector<category> tagCategories, std::vector<category> untagCategories, MassResolution massresoc );
   float                 stringToFloat( const TString varname, HggEvtCandidate& hggevent );

   //------------------------------------------------
   
   std::ofstream _xcheckTextFile;
   std::ofstream debugfile;

   bool DoSetRho0_;
   bool DoDebugEvent;
   Int_t mode_;

   float lumi_,xsec_;
   int nEvts_;

   WeightManager *_weight_manager;
   
   HggVertexAnalyzer        *_vtxAna;
   HggVertexFromConversions *_vtxAnaFromConversions;
   VertexAlgoParameters     _vtxAlgoParams;
   TMVA::Reader *_perVtxReader;
   std::string   _perVtxMvaWeights;
   std::string   _perVtxMvaMethod;
   TMVA::Reader *_perEvtReader;
   std::string   _perEvtMvaWeights;
   std::string   _perEvtMvaMethod;

   ConfigReader *_config;
   MiniTree *_minitree;     

   //CF:
   photonEnergyScale    enScaleSkimEOS; /// skim EOS bugged so need to undo the energy scale in the skim
   photonEnergyScale    enScale;
   photonEnergySmearing overSmear_;
   int                  overSmearSyst_;
   
   //ofstream outTextFile;
   TFile outFile;
   TFile *fout;
   TH1F* hPURunAB;
   TH1F* hPURunC;
   TH1F* hPURunD;
   TH1F *hTotEvents;
   TH1D *data_npu_estimated;

   TH1F *hITPU_wei;
   TH1F *hITPU;

   TH1F *hITPUTrue_wei;
   TH1F *hITPUTrue;

   TH1F *h_pho1Et;
   TH1F *h_pho2Et;
   TH1F *h_pho3Et;

   TH1F *hPhoCorrEt1;
   TH1F *hPhoCorrEt2;

   TH1F *hPhoCorrEt1_WH;
   TH1F *hPhoCorrEt2_WH;
 
   TH1F *hPhoCorrEt1_ZH;
   TH1F *hPhoCorrEt2_ZH;
   
   //------------------ lepton tag----------
   
   TH1F *hDMZmin;
   TH1F *hDPhi_gele_g;

   // TH2F *hDPhi_geleg_gele;


   TH2F *hDPhi_geleg_gele;
   
   TH1F *hcat4;
   TH1F *hggMass_eleVeto;
   TH1F *hggMass_eleVeto_jetVeto;
   
   TH1F *hDRtoTrk;
   TH1F *hggMass_muTag;
   TH1F *hggMass_eleTag;
   TH1F *hggMassWH_muTag;
   TH1F *hggMassWH_eleTag;
   TH1F *hggMassZH_muTag;
   TH1F *hggMassZH_eleTag;
   TH1F *hggMass_eletag_Zveto;
   TH1F *hggMassWH_eleTag_Zveto;
   TH1F *hggMassZH_eleTag_Zveto;
   TH1F *hDR_Pho1Ele;
   TH1F *hDR_Pho2Ele;
   TH1F *hDR_SCPho1Ele;
   TH1F *hDR_SCPho2Ele;
   TH1F *hDRmin_PhoEle;
   TH1F *hDRmin_SCPhoEle;
   TH1F *hDRmin_PhoMu;

   TH1F *hAllElePt;
   TH1F *hAllSCEleEta;
   TH1F *hEleCombRelIso_EB;
   TH1F *hEleCombRelIso_EE;
   TH1F *hSelEle1Pt;
   TH1F *hSelSCEle1Eta; 
   TH1F *hSelEle2Pt; 
   TH1F *hSelSCEle2Eta;

   TH1F *helePFRelIso_BID_B;
   TH1F *helePFRelIso_AID_B;
   TH1F *helePFRelIso_BID_EC;
   TH1F *helePFRelIso_AID_EC;

   TH1F *hmuPFRelIso_BID;
   TH1F *hmuPFRelIso_AID;

   TH1F *hMuCombRelIso;
   TH1F *hAllMuPt;
   TH1F *hAllMuEta;
  

   TH1F *hSelMu1Pt;
   TH1F *hSelMu1Eta; 
   TH1F *hSelMu2Pt; 
   TH1F *hSelMu2Eta;
  
   TH1F *h_pfMET;
   TH1F *h_tcMET;
   TH1F *h_CaloMET;

   TH1F *h_pfMET_lepTag;
   TH1F *h_tcMET_lepTag;
   TH1F *h_CaloMET_lepTag;

   TH1F *hMassEleGamGam;
   TH1F *hPtEleGamGam;

   TH1F *hMassMuGamGam;
   TH1F *hPtMuGamGam;


   TH1F *hMassEleEleGam1;
   TH1F *hMassEleEleGam2;
   TH1F *hPtEleEleGam1;
   TH1F *hPtEleEleGam2;
   TH1F *hMassEleGam1;
   TH1F *hMassEleGam2;
   TH1F *hPtEleGam1;
   TH1F *hPtEleGam2;

   TH1F *hMassEleFakePho;
   TH1F *hMassEleRealPho;



  //------------------two jet tag----------


   TH1F *hJetNHF;
   TH1F *hJetNEF;
   TH1F *hJetNConst;

   TH1F *hJetCHF;
   TH1F *hJetNCH;
   TH1F *hJetCEF;

   TH1F *hDRAllJetg1;
   TH1F *hDRAllJetg2;
   TH1F *hAllJetPt;
   TH1F *hAllJetEta;

   TH1F *hDRjet1g1;
   TH1F *hDRjet1g2;
   TH1F *hDRjet2g1;
   TH1F *hDRjet2g2;

   TH1F *hDPhijet1g1;
   TH1F *hDPhijet1g2;
   TH1F *hDPhijet2g1;
   TH1F *hDPhijet2g2;

   TH1F *hDPhijet1jet2;


   TH1F *hggMass_twoJetTag;
   TH1F *hggMass;
   TH1F *hggPt_twoJetTag;

   TH1F *hMassJetJet;   
   TH1F *hJet1Pt;
   TH1F *hJet2Pt;
   TH1F *hJet1Eta;
   TH1F *hJet2Eta;
   TH1F *hDeltaEtaJets;
   TH1F *hZeppenfeld;
   TH1F *hDeltaPhiDiJetDiPho;
  
   TH1F *hZeppenfeld_bc;
   TH1F *hDeltaEtaJets_bc;
   TH1F *hJetEtaProd_bc;
   TH1F *hMassJetJet_bc;   
   TH1F *hDeltaPhiDiJetDiPho_bc;
   TH1F *hMgg_twoJetTag_jet1Pt;
   TH1F *hMgg_twoJetTag_deltaEta;
   TH1F *hMgg_twoJetTag_zep;
   TH1F *hMgg_twoJetTag_Mjets;
   TH1F *hMgg_twoJetTag_dPhi;
  
   TH1F *hJet1Pt_HSTRA_bc;
   TH1F *hDeltaEtaJets_HSTRA_bc;
   TH1F *hZeppenfeld_HSTRA_bc;
   TH1F *hMassJetJet_HSTRA_bc;

   TH1F *hMgg_HSTRA_twoJetTag;
   TH1F *hMgg_HSTRA_twoJetTag_jet1Pt;
   TH1F *hMgg_HSTRA_twoJetTag_deltaEta;
   TH1F *hMgg_HSTRA_twoJetTag_zep;
   TH1F *hMgg_HSTRA_twoJetTag_Mjets;


   //----------------------------------------

};



#endif // #ifdef xAna_cxx
