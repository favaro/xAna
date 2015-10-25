#include "interface/xAna.hh"
#include "interface/vertexUtils.hh"
#include "interface/MinitreeOut.hh"

#include <TMVA/Tools.h>
#include <TMVA/Reader.h>

#include <string>
using namespace std;

string mvaDiscriDir = "etc/inputs/mvaDiscriminants/";

//////////////////////
//Setup PhotonID MVA//
//////////////////////
void xAna::Setup_MVA( ) {
  phoID_2011[0] = new TMVA::Reader("!Color:Silent");
  phoID_2011[1] = new TMVA::Reader("!Color:Silent");
  phoID_2012[0] = new TMVA::Reader("!Color:Silent");
  phoID_2012[1] = new TMVA::Reader("!Color:Silent");
  DiscriDiPho_2011 = new TMVA::Reader("!Color:Silent");
  DiscriDiPho_2012 = new TMVA::Reader("!Color:Silent");

  /// ------ setup 2011 mva phoID, EB ------ ///
  cout << " Booking 2011 photon ID: EB" << endl;
  string mvamethod = "BDT";
  phoID_2011[0]->AddVariable("ph.scrawe", &phoID_SCRawE);
  phoID_2011[0]->AddVariable("ph.r9",                &phoID_R9 );
  phoID_2011[0]->AddVariable("ph.sigietaieta",       &phoID_covIEtaIEta );
  phoID_2011[0]->AddVariable("ph.scetawidth",        &phoID_EtaWidth);
  phoID_2011[0]->AddVariable("ph.scphiwidth",        &phoID_PhiWidth);
  phoID_2011[0]->AddVariable("ph.idmva_CoviEtaiPhi", &phoID_covIEtaIPhi);
  phoID_2011[0]->AddVariable("ph.idmva_s4ratio",     &phoID_S4);
  phoID_2011[0]->AddVariable("ph.idmva_GammaIso",    &phoID_pfPhoIso03 );
  phoID_2011[0]->AddVariable("ph.idmva_ChargedIso_selvtx",   &phoID_pfChIso03 );
  phoID_2011[0]->AddVariable("ph.idmva_ChargedIso_worstvtx", &phoID_pfChIso03worst );
  phoID_2011[0]->AddVariable("ph.sceta",             &phoID_ScEta );
  phoID_2011[0]->AddVariable("rho",                  &phoID_rho); 
  phoID_2011[0]->BookMVA(mvamethod.c_str(), mvaDiscriDir + "2013FinalPaper_PhotonID_Barrel_BDT_TrainRangePT15_7TeV.weights.xml");

  /// ------ setup 2011 mva phoID,EE ------ ///
  cout << " Booking 2011 photon ID: EE" << endl;
  phoID_2011[1]->AddVariable("ph.scrawe", &phoID_SCRawE);
  phoID_2011[1]->AddVariable("ph.r9",                &phoID_R9 );
  phoID_2011[1]->AddVariable("ph.sigietaieta",       &phoID_covIEtaIEta );
  phoID_2011[1]->AddVariable("ph.scetawidth",        &phoID_EtaWidth);
  phoID_2011[1]->AddVariable("ph.scphiwidth",        &phoID_PhiWidth);
  phoID_2011[1]->AddVariable("ph.idmva_CoviEtaiPhi", &phoID_covIEtaIPhi);
  phoID_2011[1]->AddVariable("ph.idmva_s4ratio",     &phoID_S4);
  phoID_2011[1]->AddVariable("ph.idmva_GammaIso",    &phoID_pfPhoIso03 );
  phoID_2011[1]->AddVariable("ph.idmva_ChargedIso_selvtx",   &phoID_pfChIso03 );
  phoID_2011[1]->AddVariable("ph.idmva_ChargedIso_worstvtx", &phoID_pfChIso03worst );
  phoID_2011[1]->AddVariable("ph.sceta",             &phoID_ScEta );
  phoID_2011[1]->AddVariable("rho",                  &phoID_rho);
  phoID_2011[1]->AddVariable("ph.idmva_PsEffWidthSigmaRR",   &phoID_ESEffSigmaRR );
  phoID_2011[1]->BookMVA(mvamethod.c_str(),mvaDiscriDir + "2013FinalPaper_PhotonID_Endcap_BDT_TrainRangePT15_7TeV.weights.xml");
  
  /// ------ setup 2012 mva phoID, EB ------ ///
  cout << " Booking 2012 photon ID: EB" << endl;
  phoID_2012[0]->AddVariable("ph.scrawe", &phoID_SCRawE);
  phoID_2012[0]->AddVariable("ph.r9",                &phoID_R9 );
  phoID_2012[0]->AddVariable("ph.sigietaieta",       &phoID_covIEtaIEta );
  phoID_2012[0]->AddVariable("ph.scetawidth",        &phoID_EtaWidth);
  phoID_2012[0]->AddVariable("ph.scphiwidth",        &phoID_PhiWidth);
  phoID_2012[0]->AddVariable("ph.idmva_CoviEtaiPhi", &phoID_covIEtaIPhi);
  phoID_2012[0]->AddVariable("ph.idmva_s4ratio",     &phoID_S4);
  phoID_2012[0]->AddVariable("ph.idmva_GammaIso",    &phoID_pfPhoIso03 );
  phoID_2012[0]->AddVariable("ph.idmva_ChargedIso_selvtx",   &phoID_pfChIso03 );
  phoID_2012[0]->AddVariable("ph.idmva_ChargedIso_worstvtx", &phoID_pfChIso03worst );
  phoID_2012[0]->AddVariable("ph.sceta",             &phoID_ScEta );
  phoID_2012[0]->AddVariable("rho",                  &phoID_rho);
  phoID_2012[0]->BookMVA(mvamethod.c_str(), mvaDiscriDir + "2013FinalPaper_PhotonID_Barrel_BDT_TrainRangePT15_8TeV.weights.xml");

  /// ------ setup 2012 mva phoID, EE ------ ///
  cout << " Booking 2011 photon ID: EE" << endl;
  phoID_2012[1]->AddVariable("ph.scrawe", &phoID_SCRawE);
  phoID_2012[1]->AddVariable("ph.r9",                &phoID_R9 );
  phoID_2012[1]->AddVariable("ph.sigietaieta",       &phoID_covIEtaIEta );
  phoID_2012[1]->AddVariable("ph.scetawidth",        &phoID_EtaWidth);
  phoID_2012[1]->AddVariable("ph.scphiwidth",        &phoID_PhiWidth);
  phoID_2012[1]->AddVariable("ph.idmva_CoviEtaiPhi", &phoID_covIEtaIPhi);
  phoID_2012[1]->AddVariable("ph.idmva_s4ratio",     &phoID_S4);
  phoID_2012[1]->AddVariable("ph.idmva_GammaIso",    &phoID_pfPhoIso03 );
  phoID_2012[1]->AddVariable("ph.idmva_ChargedIso_selvtx",   &phoID_pfChIso03 );
  phoID_2012[1]->AddVariable("ph.idmva_ChargedIso_worstvtx", &phoID_pfChIso03worst );
  phoID_2012[1]->AddVariable("ph.sceta",             &phoID_ScEta );
  phoID_2012[1]->AddVariable("rho",                  &phoID_rho);
  phoID_2012[1]->AddVariable("ph.idmva_PsEffWidthSigmaRR",   &phoID_ESEffSigmaRR );
  phoID_2012[1]->BookMVA(mvamethod.c_str(), mvaDiscriDir + "2013FinalPaper_PhotonID_Endcap_BDT_TrainRangePT15_8TeV.weights.xml");
  
  /// ------ setup 2011 diphoton ------ ///
  cout << " Booking 2011 MVA diphoton discriminant" << endl;
  DiscriDiPho_2011->AddVariable("masserrsmeared/mass",&Ddipho_masserr);
  DiscriDiPho_2011->AddVariable("masserrsmearedwrongvtx/mass",&Ddipho_masserrwrongvtx);
  DiscriDiPho_2011->AddVariable("vtxprob",&Ddipho_vtxprob);
  DiscriDiPho_2011->AddVariable("ph1.pt/mass",&Ddipho_piT1);
  DiscriDiPho_2011->AddVariable("ph2.pt/mass",&Ddipho_piT2);
  DiscriDiPho_2011->AddVariable("ph1.eta",&Ddipho_eta1);
  DiscriDiPho_2011->AddVariable("ph2.eta",&Ddipho_eta2);
  DiscriDiPho_2011->AddVariable("TMath::Cos(ph1.phi-ph2.phi)",&Ddipho_cosdPhi);
  DiscriDiPho_2011->AddVariable("ph1.idmva",&Ddipho_id1);
  DiscriDiPho_2011->AddVariable("ph2.idmva",&Ddipho_id2);
  //DiscriDiPho_2011->BookMVA("BDTG","mvaDiscriInputs/HggBambu_SMDipho_Jan16_BDTG.weights.xml");
  DiscriDiPho_2011->BookMVA("BDTG", mvaDiscriDir + "HggBambu_SMDipho_Oct29_allsigevenbkg7TeV_BDTG.weights.xml");

  /// ------ setup 2012 diphoton ------ ///
  cout << " Booking 2012 MVA diphoton discriminant" << endl;
  DiscriDiPho_2012->AddVariable("masserrsmeared/mass",&Ddipho_masserr);
  DiscriDiPho_2012->AddVariable("masserrsmearedwrongvtx/mass",&Ddipho_masserrwrongvtx);
  DiscriDiPho_2012->AddVariable("vtxprob",&Ddipho_vtxprob);
  DiscriDiPho_2012->AddVariable("ph1.pt/mass",&Ddipho_piT1);
  DiscriDiPho_2012->AddVariable("ph2.pt/mass",&Ddipho_piT2);
  DiscriDiPho_2012->AddVariable("ph1.eta",&Ddipho_eta1);
  DiscriDiPho_2012->AddVariable("ph2.eta",&Ddipho_eta2);
  DiscriDiPho_2012->AddVariable("TMath::Cos(ph1.phi-ph2.phi)",&Ddipho_cosdPhi);
  DiscriDiPho_2012->AddVariable("ph1.idmva",&Ddipho_id1);
  DiscriDiPho_2012->AddVariable("ph2.idmva",&Ddipho_id2);
  //  DiscriDiPho_2012->BookMVA("BDTG","mvaDiscriInputs/HggBambu_SMDipho_Jun19_BDTG.weights.xml");
  DiscriDiPho_2012->BookMVA("BDTG", mvaDiscriDir + "HggBambu_SMDipho_Oct01_redqcdweightallsigevenbkg_BDTG.weights.xml");




  if( _config->setup() == "Legacy2013_7TeV" )  for( int i = 0 ; i < 2; i++ ) phoID_mva[i] = phoID_2011[i];
  else                                         for( int i = 0 ; i < 2; i++ ) phoID_mva[i] = phoID_2012[i];
  if( _config->setup() == "Legacy2013_7TeV" )  Ddipho_mva = DiscriDiPho_2011;
  else                                         Ddipho_mva = DiscriDiPho_2012;


  /// ------------------------------///
//  /// ----- setup 2012 MVA VBF -----///
//  DiscriVBF = new TMVA::Reader("!Color:Silent");
//  string tmvaVbfMvaWeights = "mvaDiscriInputs/TMVA_vbf_6var_mjj100_diphopt_phopt_BDTG.weights.xml";
//  DiscriVBF_Method = "BDTG";
//
//  cout << " Booking 2012 MVA vbf discriminant" << endl;
//  DiscriVBF = new TMVA::Reader( "!Color:!Silent" );
//  DiscriVBF->AddVariable("jet1pt", &myVBFLeadJPt);
//  DiscriVBF->AddVariable("jet2pt", &myVBFSubJPt);
//  DiscriVBF->AddVariable("abs(jet1eta-jet2eta)", &myVBFdEta);
//  DiscriVBF->AddVariable("mj1j2", &myVBF_Mjj);
//  DiscriVBF->AddVariable("zepp" , &myVBFZep);
//  DiscriVBF->AddVariable("dphi" , &myVBFdPhi);
//  if( DiscriVBF_UseDiPhoPt )  DiscriVBF->AddVariable("diphopt/diphoM", &myVBFDiPhoPtOverM); 
//  if( DiscriVBF_UsePhoPt   ) {
//      DiscriVBF->AddVariable("pho1pt/diphoM", &myVBFLeadPhoPtOverM);
//      DiscriVBF->AddVariable("pho2pt/diphoM", &myVBFSubPhoPtOverM);
//    }
//  DiscriVBF->BookMVA( DiscriVBF_Method, tmvaVbfMvaWeights );
//  // end of MVA VBF 
  /// ----- setup 2012 MVA VBF -----///
  string tmvaVbfMvaWeights = mvaDiscriDir + "TMVA_dijet_sherpa_scalewt50_2evenb_powheg200_maxdPhi_oct9_Gradient.weights.xml";
  DiscriVBF_Method = "BDTG";
  cout << " Booking 2012 MVA vbf discriminant" << endl;
  DiscriVBF = new TMVA::Reader( "!Color:!Silent" );
  DiscriVBF->AddVariable("dijet_leadEta", &myVBFLeadJEta);
  DiscriVBF->AddVariable("dijet_subleadEta", &myVBFSubJEta);
  DiscriVBF->AddVariable("dijet_LeadJPt", &myVBFLeadJPt);
  DiscriVBF->AddVariable("dijet_SubJPt", &myVBFSubJPt);
  DiscriVBF->AddVariable("dijet_Zep" , &myVBFZep);
  DiscriVBF->AddVariable("min(dijet_dPhi,2.916)", &myVBFdPhiTrunc);
  DiscriVBF->AddVariable("dijet_Mjj", &myVBF_Mjj);
  DiscriVBF->AddVariable("dipho_pt/mass", &myVBFDiPhoPtOverM); 
  DiscriVBF->BookMVA( DiscriVBF_Method, tmvaVbfMvaWeights );

  DiscriVBFCombDipho = new TMVA::Reader("!Color:!Silent"); 
  string tmvaVbfCombDiphoMvaWeights = mvaDiscriDir + "TMVA_vbf_dijet_dipho_evenbkg_scaledwt50_maxdPhi_Gradient.weights.xml";
  cout << " Booking 2012 MVA combined vbf-dipho discriminant" << endl;
  DiscriVBFCombDipho->AddVariable("dipho_mva",                       &myVBFDIPHObdt);
  DiscriVBFCombDipho->AddVariable("bdt_dijet_maxdPhi",               &myVBF_MVA);
  DiscriVBFCombDipho->AddVariable("dipho_pt/mass",                  &myVBFDiPhoPtOverM);
  DiscriVBFCombDipho->BookMVA( DiscriVBF_Method, tmvaVbfCombDiphoMvaWeights );

  // end of MVA VBF 


 

  /// -------------------------------///
  /// ----- setup vertex finder -----///
  vector<string> tmvaPerVtxVariables;
  tmvaPerVtxVariables.push_back("ptbal");
  tmvaPerVtxVariables.push_back("ptasym");
  tmvaPerVtxVariables.push_back("logsumpt2");  
  tmvaPerVtxVariables.push_back("limPullToConv");
  tmvaPerVtxVariables.push_back("nConv");
  string perVtxMvaWeights = mvaDiscriDir + "TMVAClassification_BDTCat_conversions_tmva_407.weights.xml";
  _perVtxMvaMethod  = "BDTCat_conversions";
  string perEvtMvaWeights = mvaDiscriDir + "TMVAClassification_evtBDTG_conversions_tmva_407.weights.xml";
  _perEvtMvaMethod  = "evtBDTG";
  string vertexsetup = mvaDiscriDir + "vertex_selection_7TeV.dat";
  if( _config->setup().find("2012") != string::npos || _config->setup().find("2013") != string::npos) {
    perEvtMvaWeights = mvaDiscriDir + "TMVAClassification_BDTvtxprob2012.weights.xml";
    vertexsetup      = mvaDiscriDir + "vertex_selection_hcp2012.dat";
    _perEvtMvaMethod = "BDTvtxprob2012";
  }
  _vtxAlgoParams = VertexAlgoParametersReader(vertexsetup);
  _vtxAna  = new HggVertexAnalyzer(_vtxAlgoParams);

  _perVtxReader = 0;
  _perEvtReader = 0;
  if(  perVtxMvaWeights != "" ) { 
    _perVtxReader = new TMVA::Reader( "!Color:!Silent" );
    HggVertexAnalyzer::bookVariables( *_perVtxReader, tmvaPerVtxVariables );
    _perVtxReader->BookMVA( _perVtxMvaMethod, perVtxMvaWeights );
  }
  if( perEvtMvaWeights != "" ) {
    _perEvtReader = new TMVA::Reader( "!Color:!Silent" );
    HggVertexAnalyzer::bookPerEventVariables( *_perEvtReader );
    _perEvtReader->BookMVA( _perEvtMvaMethod, perEvtMvaWeights );
  }
  cout << " Booking Vertex discriminants" << endl;
  cout << "       - use conversions     : " << _vtxAlgoParams.useAllConversions << endl;
  cout << "       - single leg sigma1Pix: " << _vtxAlgoParams.singlelegsigma1Pix << endl;
  cout << "       - check that if use single leg conversion (option=3) then sigma1Pix has meaningful value "
       << endl;

  // end of Vertex MVA 


//BOOKING energy corrections:
  cout << " Booking GBRLiklihood ECorr" << endl;
  TFile* fgbrPho=new TFile("etc/external/GBRLikelihood/data/regweights_v5_forest_ph.root");
  TFile* fgbrEle=new TFile("etc/external/GBRLikelihood/data/regweights_v5_forest_ele.root");  
  fgbrPho->GetObject("EGRegressionForest_EB", _forestebPho);
  fgbrPho->GetObject("EGRegressionForest_EE", _foresteePho);
  fgbrEle->GetObject("EGRegressionForest_EB", _forestebEle);
  fgbrEle->GetObject("EGRegressionForest_EE", _foresteeEle);

  _tgt   = new RooRealVar("tgt","",1.);
  _mean  = new RooRealVar("mean","",1.);
  _sigma = new RooRealVar("sigma","",1.);
  _n1    = new RooRealVar("n1","",2.);
  _n2    = new RooRealVar("n2","",2.);
  
  _sigmalim = new RooRealConstraint("sigmalim","",*_sigma,0.0002,0.5);
  _meanlim  = new RooRealConstraint("meanlim","",*_mean,0.2,2.0);
  _n1lim    = new RooRealConstraint("n1lim","",*_n1,1.01,110.);
  _n2lim    = new RooRealConstraint("n2lim","",*_n2,1.01,110.);
  
  RooConstVar *cbmean = new RooConstVar("cbmean","",1.0);
  RooConstVar *alpha1 = new RooConstVar("alpha1","",2.0);
  RooConstVar *alpha2 = new RooConstVar("alpha2","",1.0);
  
  
  //  _pdf = new RooDoubleCBFast("sigpdf","",*_tgt,RooFit::RooConst(1.),*_sigmalim,RooFit::RooConst(2.0),*_n1lim,RooFit::RooConst(1.0),*_n2lim);
  //add to RooArgList for proper garbage collection
  
  _args.addOwned(*_tgt);
  _args.addOwned(*_mean);
  _args.addOwned(*_sigma);
  _args.addOwned(*_n1);
  _args.addOwned(*_n2);
  _args.addOwned(*cbmean);
  _args.addOwned(*alpha1);
  _args.addOwned(*alpha2);
  _args.addOwned(*_sigmalim);
  _args.addOwned(*_meanlim);
  _args.addOwned(*_n1lim);
  _args.addOwned(*_n2lim);
  //   _args.addOwned(*_pdf);
  
  delete fgbrPho;
  delete fgbrEle;
}


float xAna::PhoID_MVA( int i, int selVtx ) {
  ////////////////////
  //do photon id mva//
  ////////////////////
  phoID_HoE         = (*phoHoverE)[i];
  phoID_covIEtaIEta = (*phoSigmaIEtaIEta)[i];
  phoID_R9          = (*phoR9)[i];
  phoID_absIsoEcal  = (*phoEcalIsoDR03)[i] - 0.17*rho25;
  phoID_absIsoHcal  = (*phoHcalIsoDR04)[i] - 0.17*rho25;
  phoID_NVertexes   = nVtx;
  phoID_ScEta       = (*phoSCEta)[i];
  phoID_EtaWidth    = (*phoSCEtaWidth)[i];
  phoID_PhiWidth    = (*phoSCPhiWidth)[i];
  phoID_SCRawE      =(*phoSCRawE)[i];

  if( _config->setup() == "Legacy2013_7TeV" || _config->setup() == "Prompt2012_ichep" || _config->setup() == "Legacy2013_8TeV" ) {
    phoID_covIEtaIPhi =  (*phoSigmaIEtaIPhi)[i];
    phoID_S4 = (*phoS4ratio)[i];
    phoID_ESEffSigmaRR = (*phoESEffSigmaRR_x)[i];
    phoID_rho = rho2012;
    phoID_pfPhoIso03 = (*phoCiCPF4phopfIso03)[i];
    phoID_pfChIso03  = (*phoCiCPF4chgpfIso03)[i][selVtx];
    phoID_pfChIso03worst = -99;
    for( int iv=0; iv<nVtxBS ; iv++) 
      if( (*phoCiCPF4chgpfIso03)[i][iv] >  phoID_pfChIso03worst )
	phoID_pfChIso03worst = (*phoCiCPF4chgpfIso03)[i][iv];
  } else {
    int badvtx;
    float pho_tkiso_badvtx = worstSumTrackPtInCone(i,badvtx,0,0.4,0.02,0.0,1.0,0.1);
    float pho_tkiso_goodvtx = sumTrackPtInCone(i,selVtx,0,0.3,0.02,0.0,1.0,0.1);
    phoID_tIso1abs = (pho_tkiso_goodvtx + (*phoEcalIsoDR03)[i] + (*phoHcalIsoDR04)[i] - 0.17*rho25);
    phoID_tIso2abs = (pho_tkiso_badvtx + (*phoEcalIsoDR04)[i] + (*phoHcalIsoDR04)[i] - 0.52*rho25);
    phoID_tIso3abs = pho_tkiso_goodvtx;
  }

  return -1;
}


float xAna::DiPhoID_MVA( const TLorentzVector &glead, const TLorentzVector &gtrail, const TLorentzVector &hcand, 
				       const MassResolution & massResoCalc, float vtxProb,
				       const float &idlead, const float &idtrail) {
  //Ddipho_masserr         = massResoCalc.relativeMassResolutionCorrVertex();
  Ddipho_masserr         = massResoCalc.relativeMassResolutionFab_energy();
  Ddipho_masserrwrongvtx = massResoCalc.relativeMassResolutionWrongVertex();

  Ddipho_vtxprob = vtxProb;
  Ddipho_piT1 = glead.Pt()/hcand.M();
  Ddipho_piT2 = gtrail.Pt()/hcand.M();
  Ddipho_eta1 = glead.Eta();
  Ddipho_eta2 = gtrail.Eta();
  Ddipho_cosdPhi = TMath::Cos(glead.Phi()-gtrail.Phi());
  Ddipho_id1 = idlead;
  Ddipho_id2 = idtrail;
  
  float diphomva = -1;
  if(  Ddipho_id1 > -0.3 &&  Ddipho_id2 > -0.3  ) 
    diphomva = Ddipho_mva->EvaluateMVA("BDTG");
  
  
  /// fill minitree with dipho inputs
  _minitree->mtree_massResoRightVtx = Ddipho_masserr;
  _minitree->mtree_massResoRandVtx  = Ddipho_masserrwrongvtx;
  _minitree->mtree_vtxProb          = Ddipho_vtxprob;
  
  _minitree->mtree_diphoMva = diphomva;

  return diphomva;
}

void xAna::CorrLikelihood(int ipho, float &ecorr, float &esigma, bool forPho){
	
  std::vector<float> GBRInvals;
  GBRInvals.resize(37);
  if(forPho){
    GBRInvals[0]=(*phoSCRawE)[ipho];
    GBRInvals[1]=(*phoSCEta)[ipho];
    GBRInvals[2]=(*phoR9)[ipho];
    GBRInvals[3]=(*phoSCEtaWidth)[ipho];
    GBRInvals[4]=(*phoSCPhiWidth)[ipho];
    GBRInvals[5]=(*phoNClus)[ipho];
    GBRInvals[6]=(*phoHoverE12)[ipho];
    GBRInvals[7]= rho2012;
    //GBRInvals[8]=nVtx; //CF
    GBRInvals[8]  = nVtxBS;
    GBRInvals[9]  = (*phoSeedEta)[ipho]- (*phoSCEta)[ipho];
    GBRInvals[10] = asin(sin((*phoSeedPhi)[ipho]-(*phoSCPhi)[ipho]));
    GBRInvals[11] = (*phoSeedE)[ipho]/ (*phoSCRawE)[ipho];
    GBRInvals[12] = (*phoE3x3)[ipho]/ (*phoE5x5)[ipho];
    GBRInvals[13] = ((*phoSigmaIEtaIEta)[ipho]);
    GBRInvals[14] = sqrt((*phoSigmaIPhiIPhi)[ipho]);
    GBRInvals[15] = ((*phoSigmaIEtaIPhi)[ipho]);
    GBRInvals[16] = (*phoEmax)[ipho]/ (*phoE5x5)[ipho];
    GBRInvals[17]=(*phoE2ndMax)[ipho]/ (*phoE5x5)[ipho];
    GBRInvals[18]=(*phoETop)[ipho]/ (*phoE5x5)[ipho];
    GBRInvals[19]=(*phoEBottom)[ipho]/ (*phoE5x5)[ipho];
    GBRInvals[20]=(*phoELeft)[ipho]/ (*phoE5x5)[ipho];
    GBRInvals[21]=(*phoERight)[ipho]/(*phoE5x5)[ipho];
    GBRInvals[22]=(*phoE2x5Max)[ipho]/(*phoE5x5)[ipho];
    GBRInvals[23]=(*phoE2x5Top)[ipho]/ (*phoE5x5)[ipho];
    GBRInvals[24]=(*phoE2x5Bottom)[ipho]/(*phoE5x5)[ipho];
    GBRInvals[25]=(*phoE2x5Left)[ipho]/(*phoE5x5)[ipho];
    GBRInvals[26]=(*phoE2x5Right)[ipho]/ (*phoE5x5)[ipho];
    if(fabs((*phoSCEta)[ipho]) <1.479){
      GBRInvals[27]= (*phoE5x5)[ipho]/(*phoSeedE)[ipho];
      GBRInvals[28]=(*phoCrysIEta)[ipho];
      GBRInvals[29]=(*phoCrysIPhi)[ipho]%18;
      GBRInvals[30]=((*phoCrysIEta)[ipho]%5);
      GBRInvals[31]=((*phoCrysIPhi)[ipho]%2);
      int ieta=(int)(*phoCrysIEta)[ipho];
      int iphi=(int)(*phoCrysIPhi)[ipho];  
      if(ieta!=0)GBRInvals[32]=(TMath::Abs(ieta)<=25)*(ieta%25) + (TMath::Abs(ieta)>25)*((ieta-25*TMath::Abs(ieta)/ieta)%20);
      else GBRInvals[32]=(TMath::Abs(ieta)<=25)*(ieta%25) + (TMath::Abs(ieta)>25)*((ieta-25*0)%20);
      GBRInvals[33]=(iphi%20);
      GBRInvals[34]=(*phoCrysEta)[ipho];
      GBRInvals[35]=(*phoCrysPhi)[ipho];
    }
    else{
      GBRInvals[27]=(*phoESEn)[ipho]/(*phoSCRawE)[ipho];
    }
  }
  else{//ecorr for electrons
    GBRInvals[0]=(*eleSCRawEn)[ipho];
    GBRInvals[1]=(*eleSCEta)[ipho];
    GBRInvals[2]=(*eleR9)[ipho];
    GBRInvals[3]=(*eleSCEtaWidth)[ipho];
    GBRInvals[4]=(*eleSCPhiWidth)[ipho];
    GBRInvals[5]=(*eleNClus)[ipho];
    GBRInvals[6]=(*eleHoverE12)[ipho];
    GBRInvals[7]=rho2012;
    GBRInvals[8]=nVtx;
    GBRInvals[9]=(*eleSeedEta)[ipho]- (*eleSCEta)[ipho];
    GBRInvals[10]=asin(sin((*eleSeedPhi)[ipho]-(*eleSCPhi)[ipho]));
    GBRInvals[11]=(*eleSeedE)[ipho]/ (*eleSCRawEn)[ipho];
    GBRInvals[12]=(*eleE3x3)[ipho]/ (*eleE5x5)[ipho];
    GBRInvals[13] = ((*eleSigmaIEtaIEta)[ipho]);
    GBRInvals[14]=sqrt((*eleSigmaIPhiIPhi)[ipho]);
    GBRInvals[15]=((*eleSigmaIEtaIPhi)[ipho]);
    GBRInvals[16]=(*eleEmax)[ipho]/ (*eleE5x5)[ipho];
    GBRInvals[17]=(*eleE2ndMax)[ipho]/ (*eleE5x5)[ipho];
    GBRInvals[18]=(*eleETop)[ipho]/ (*eleE5x5)[ipho];
    GBRInvals[19]=(*eleEBottom)[ipho]/ (*eleE5x5)[ipho];
    GBRInvals[20]=(*eleELeft)[ipho]/ (*eleE5x5)[ipho];
    GBRInvals[21]=(*eleERight)[ipho]/(*eleE5x5)[ipho];
    GBRInvals[22]=(*eleE2x5Max)[ipho]/(*eleE5x5)[ipho];
    GBRInvals[23]=(*eleE2x5Top)[ipho]/ (*eleE5x5)[ipho];
    GBRInvals[24]=(*eleE2x5Bottom)[ipho]/(*eleE5x5)[ipho];
    GBRInvals[25]=(*eleE2x5Left)[ipho]/(*eleE5x5)[ipho];
    GBRInvals[26]=(*eleE2x5Right)[ipho]/ (*eleE5x5)[ipho];
    if(fabs((*eleSCEta)[ipho]) <1.479){
      GBRInvals[27]= (*eleE5x5)[ipho]/(*eleSeedE)[ipho];
      GBRInvals[28]=(*eleCrysIEta)[ipho];
      GBRInvals[29]=(*eleCrysIPhi)[ipho]%18;
      GBRInvals[30]=((*eleCrysIEta)[ipho]%5);
      GBRInvals[31]=((*eleCrysIPhi)[ipho]%2);
      int ieta=(int)(*eleCrysIEta)[ipho];
      int iphi=(int)(*eleCrysIPhi)[ipho];  
      if(ieta!=0)GBRInvals[32]=(TMath::Abs(ieta)<=25)*(ieta%25) + (TMath::Abs(ieta)>25)*((ieta-25*TMath::Abs(ieta)/ieta)%20);
      else GBRInvals[32]=(TMath::Abs(ieta)<=25)*(ieta%25) + (TMath::Abs(ieta)>25)*((ieta-25*0)%20);
      GBRInvals[33]=(iphi%20);
      GBRInvals[34]=(*eleCrysEta)[ipho];
      GBRInvals[35]=(*eleCrysPhi)[ipho];
    }
    else{
      GBRInvals[27]=(*eleESEn)[ipho]/(*eleSCRawEn)[ipho];
    }
  }
     double den;
  HybridGBRForest *forest;
 if(forPho){
  	if(fabs((*phoSCEta)[ipho]) <1.479){
    	den = (*phoSCRawE)[ipho];
    	forest = _forestebPho;
  	}
  	else {
   	 den = (*phoSCRawE)[ipho] + (*phoESEn)[ipho];
    	forest = _foresteePho;
  	}
  }
  else{
    if(fabs((*eleSCEta)[ipho]) <1.479){
      den = (*eleSCRawEn)[ipho];
      forest = _forestebEle;
    }
    else {
      den = (*eleSCRawEn)[ipho] + (*eleESEn)[ipho];
      forest = _foresteeEle;
    }
    
  }

 _tgt->setVal(1.0); 
 _sigma->setVal(forest->GetResponse(&GBRInvals[0],0));
 _mean->setVal(forest->GetResponse(&GBRInvals[0],1));
 _n1->setVal(forest->GetResponse(&GBRInvals[0],2));
 _n2->setVal(forest->GetResponse(&GBRInvals[0],3));
 ecorr = den/_meanlim->getVal();
 esigma = _sigmalim->getVal()*ecorr;

 }




///-------------------------------------------------------------------------///
/// function not really used anymore: compute trk iso for old MVA photon id ///

float xAna::sumTrackPtInCone(int phoID, int vtxID, float minPt, 
					   float outerConeRadius, float innerConeRadius,
					   float etaStripHalfWidth, float dzMax, float d0Max) {
  if(vtxID < 0)return -999;
  TVector3 vtxpos((*vtx_x)[vtxID],(*vtx_x)[vtxID],(*vtx_z)[vtxID]);
  float sum = 0;
  for(int i = 0; i < nTrk; i++){
    TVector3 trackp((*trkP_x)[i],(*trkP_y)[i],(*trkP_z)[i]);
    if(trackp.Pt() < minPt)continue;
    TVector3 trkVtxVec((*trkVtx_x)[i],(*trkVtx_y)[i],(*trkVtx_z)[i]);
    double deltaz = fabs((vtxpos.Z() - trkVtxVec.Z()) - 
			 ( (trkVtxVec.X()-vtxpos.X())*trackp.Px() + (trkVtxVec.Y()-vtxpos.Y())*trackp.Py() )/trackp.Pt() * trackp.Pz()/trackp.Pt() );
    if(deltaz > dzMax)continue;
    double dxy = (-(trkVtxVec.X() - vtxpos.X())*trackp.Py() + (trkVtxVec.Y() - vtxpos.Y())*trackp.Px())/trackp.Pt();
    if(fabs(dxy) > d0Max)continue;
    double deta = fabs((*phoEtaVtx)[phoID][vtxID] - trackp.Eta());
    double dphi = fabs((*phoPhiVtx)[phoID][vtxID] - trackp.Phi());
    if(dphi > TMath::Pi())dphi = TMath::TwoPi() - dphi;
    double dR = sqrt(deta*deta + dphi*dphi);
    if(dR < outerConeRadius && dR >= innerConeRadius && deta >= etaStripHalfWidth)sum += trackp.Pt();
  }
  return sum;
}

float xAna::worstSumTrackPtInCone(int phoID, int &vtxID, float minPt, 
				  float outerConeRadius, float innerConeRadius, 
				  float etaStripHalfWidth, float dzMax, float d0Max) {
  int worstvtxID = -1;
  float maxisosum = -1;
  for(int i = 0; i < nVtx; i++){
    float isosum = sumTrackPtInCone(phoID,i,minPt,outerConeRadius,innerConeRadius,etaStripHalfWidth,dzMax,d0Max);
    if(isosum > maxisosum){
      maxisosum = isosum;
      worstvtxID = i;
    }
  }
  vtxID = worstvtxID;
  return maxisosum;
}


void xAna::evaluateMELA( const HggEvtCandidate hggevent )
{

  TVar::VerbosityLevel verbosity = TVar::ERROR; 
  TEvtProb Xcal2;  

  TLorentzVector kinematics[3];
  if( hggevent.njets < 2 ) cout << " ****** less than two jets! I'm about to crash! ****** " << endl; 

  size_t jl = hggevent.jetIdx[0];
  size_t jt = hggevent.jetIdx[1]; 
  size_t j1, j2;
  if( (*jetEta)[jl] > (*jetEta)[jt] ) { j1 = jl; j2 = jt;}
  else { j1 = jt; j2 = jl;}

  kinematics[0].SetPtEtaPhiE((*jetPt)[j1], (*jetEta)[j1], (*jetPhi)[j1], (*jetEn)[j1]);
  kinematics[1].SetPtEtaPhiE((*jetPt)[j2], (*jetEta)[j2], (*jetPhi)[j2], (*jetEn)[j2]);
  kinematics[2] = hggevent.gammaK_lead + hggevent.gammaK_trail;

  _minitree->mela_dXsec_HJJ       = Xcal2.XsecCalcXJJ(TVar::HJJNONVBF, kinematics, verbosity);
  _minitree->mela_dXsec_HJJVBF    = Xcal2.XsecCalcXJJ(TVar::HJJVBF, kinematics, verbosity);

  _minitree->mela_dXsec_HJJ_PS    = Xcal2.XsecCalcXJJ(TVar::HJJNONVBF_PS, kinematics, verbosity);
  _minitree->mela_dXsec_HJJVBF_PS = Xcal2.XsecCalcXJJ(TVar::HJJVBF_PS, kinematics, verbosity);
  
  _minitree->mela_VBFvsgg    = ( _minitree->mela_dXsec_HJJVBF/(_minitree->mela_dXsec_HJJVBF + 3*_minitree->mela_dXsec_HJJ) );   
  _minitree->mela_SMvsPS_VBF = ( _minitree->mela_dXsec_HJJVBF/(_minitree->mela_dXsec_HJJVBF + 0.1*_minitree->mela_dXsec_HJJVBF_PS) );
  _minitree->mela_SMvsPS_gg  = ( _minitree->mela_dXsec_HJJ/(_minitree->mela_dXsec_HJJ+_minitree->mela_dXsec_HJJ_PS) ); 

}

