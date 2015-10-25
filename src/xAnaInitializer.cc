#include "interface/xAna.hh"

using namespace std;

void xAna::SetConfigReader( const ConfigReader &config ) {
  if( _config != 0 ) delete _config;
  _config = new ConfigReader(config);
}

//JM: add isItData to disable MC branches when data
xAna::xAna(TTree *tree, bool isItData ) {
   /// create missing vectors 
   /// FIX ME some of them should be there alredy!!!!
   phoRegrErr      = new std::vector<float>();
   phoS4ratio      = new std::vector<float>();
   phoRegrSmear    = new std::vector<float>();
   eleNClus        = new std::vector<int>();
   jetRawPtSmeared = new std::vector<float>();
   jetRawEnSmeared = new std::vector<float>();
   phoStdE         = new std::vector<float>();

   phoEnergyScale  = new std::vector<float>(); //CF
   
  
  for( int i = 0 ; i < 2 ; i++ ) {
    phoID_2011[i] = 0;
    phoID_2012[i] = 0;
  }
  _config = 0;
  
  ResJetPar = 0;
  L3JetPar  = 0;
  L2JetPar  = 0;
  L1JetPar  = 0;  
  JetCorrector = 0;
  Init(tree,isItData);
  
}


Double_t xAna::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {

  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);

  return sqrt(dEta*dEta+dPhi*dPhi);
}

Double_t xAna::deltaPhi(Double_t phi1, Double_t phi2) {

  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi())  dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();

  return dPhi;
}


Double_t xAna::getHqTWeight(Float_t HiggsMCMass, Float_t HiggsMCPt) {
  // HiggsMCMass=120;
  Double_t wei = 1;
  Double_t x = (Double_t) HiggsMCPt;
  

  Double_t par[6];
  if (HiggsMCMass == 90) {
    par[0] = 1.5512; par[1] = -0.0415444; par[2] = 0.00110316; par[3] = -1.51299e-05; par[4] = 7.51189e-08; par[5] = 0.574011;
  } else if (HiggsMCMass == 95) {
    par[0] = 1.50275; par[1] = -0.0342925; par[2] = 0.000796335; par[3] = -9.8689e-06; par[4] = 4.48063e-08; par[5] = 0.570917;
  } else if (HiggsMCMass == 100) {
    par[0] = 1.45047; par[1] = -0.0284948; par[2] = 0.000624926; par[3] = -7.66006e-06; par[4] = 3.42959e-08; par[5] = 0.574828;
  } else if (HiggsMCMass == 105) {
    par[0] = 1.29301; par[1] = -0.018701; par[2] = 0.000466271; par[3] = -6.21643e-06; par[4] = 2.84922e-08; par[5] = 0.659637;
  } else if (HiggsMCMass == 110) {
    par[0] = 1.41094; par[1] = -0.0216905; par[2] = 0.000379204; par[3] = -4.22212e-06; par[4] = 1.79075e-08; par[5] = 0.569545;
  } else if (HiggsMCMass == 115) {
    par[0] = 1.38536; par[1] = -0.0196079; par[2] = 0.000340012; par[3] = -3.85773e-06; par[4] = 1.64619e-08; par[5] = 0.569461;
  } else if (HiggsMCMass == 120) {
    par[0] = 1.40129; par[1] = -0.0206659; par[2] = 0.000339413; par[3] = -3.33799e-06; par[4] = 1.23947e-08; par[5] = 0.575397;
  } else if (HiggsMCMass == 121) {
    par[0] = 1.41674; par[1] = -0.0216618; par[2] = 0.000382855; par[3] = -4.06627e-06; par[4] = 1.59015e-08; par[5] = 0.563047;
  } else if (HiggsMCMass == 123) {
    par[0] = 1.38972; par[1] = -0.0189709; par[2] = 0.000287065; par[3] = -2.7141e-06; par[4] = 9.80252e-09; par[5] = 0.575012;
  } else if (HiggsMCMass == 125) {
    par[0] = 1.36709; par[1] = -0.0164974; par[2] = 0.000211968; par[3] = -1.82729e-06; par[4] = 6.32798e-09; par[5] = 0.588449;
  } else if (HiggsMCMass == 130) {
    par[0] = 1.40301; par[1] = -0.0204618; par[2] = 0.00034626; par[3] = -3.42286e-06; par[4] = 1.24308e-08; par[5] = 0.57312;
  } else if (HiggsMCMass == 135) {
    par[0] = 1.36409; par[1] = -0.0175548; par[2] = 0.000286312; par[3] = -2.86484e-06; par[4] = 1.05347e-08; par[5] = 0.560269;
  } else if (HiggsMCMass == 140) {
    par[0] = 1.32729; par[1] = -0.0133144; par[2] = 0.000156537; par[3] = -1.31663e-06; par[4] = 4.32718e-09; par[5] = 0.590839;
  } else if (HiggsMCMass == 150) {
    par[0] = 1.3531; par[1] = -0.0152357; par[2] = 0.000211272; par[3] = -1.84944e-06; par[4] = 6.00008e-09; par[5] = 0.582344;
  } else if (HiggsMCMass == 155) {
    par[0] = 1.33302; par[1] = -0.0132161; par[2] = 0.000154738; par[3] = -1.20301e-06; par[4] = 3.58784e-09; par[5] = 0.569991;
  } else if (HiggsMCMass == 160) {
    par[0] = 1.30238; par[1] = -0.0113894; par[2] = 0.000130351; par[3] = -1.05532e-06; par[4] = 3.16211e-09; par[5] = 0.598551;
  } else {
    return wei; 
  }

  wei = par[0] + par[1] * x + par[2] * pow(x, 2) + par[3] * pow(x, 3) + par[4] * pow(x, 4);
  if (x > HiggsMCMass) wei = par[5]; 

  return wei;
}



xAna::~xAna() {
   if (!fChain) return;

   /// create missing vectors 
   /// FIX ME some of them should be there alredy!!!!
     
   delete ResJetPar;
   delete L3JetPar ;
   delete L2JetPar ;
   delete L1JetPar ;  
   delete JetCorrector;
   delete phoRegrErr  ;
   delete phoS4ratio  ;
   delete phoRegrSmear;
   delete eleNClus    ;
   delete jetRawPtSmeared;
   delete jetRawEnSmeared;
   delete phoStdE;
   delete phoEnergyScale;
   
   
   ///   delete _pdf;
   delete _forestebPho;
   delete _foresteePho;
   
   
   //  delete fChain->GetCurrentFile();
}

Int_t xAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   int res = fChain->GetEntry(entry);
   /// resize additional vectors
   phoRegrErr  ->resize( phoE->size() );
   phoS4ratio  ->resize( phoE->size() ); 
   phoRegrSmear->resize( phoE->size() );
   eleNClus    ->resize( eleEn->size() );
   jetRawPtSmeared->resize( jetEn->size() );
   jetRawEnSmeared->resize( jetEn->size() );

   phoStdE        ->resize( phoE->size() );
   phoEnergyScale -> resize( phoE->size() ); // CF
   phoRegrE       -> resize( phoE->size() ); // CF
   return res;
}


Long64_t xAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void xAna::InitHists() {

  /* float xpt[11] = {10, 15, 20, 30, 40, 50, 60, 80, 100, 140, 200}; */
  /* float xvt[8]  = { 0,  2,  4,  6,  8, 10, 12,  20}; */

  hTotEvents =   new TH1F("hTotEvents",   "hTotEvents",3,0,3);
  //hITPU_wei =  new TH1F("hITPU_wei",   "hITPU_wei",60,0,60);
  // hITPU =  new TH1F("hITPU",   "hITPU",60,0,60);

  
  hITPU_wei =  new TH1F("hITPU_wei",   "hITPU_wei",100,0,100);
  hITPU =  new TH1F("hITPU",   "hITPU",100,0,100);

  hITPUTrue_wei =  new TH1F("hITPUTrue_wei",   "hITPUTrue_wei",100,0,100);
  hITPUTrue =  new TH1F("hITPUTrue",   "hITPUTrue",100,0,100);

  h_pho1Et =  new TH1F("h_pho1Et",   "h_pho1Et",200,0,200);
  h_pho2Et =  new TH1F("h_pho2Et",   "h_pho2Et",200,0,200);
  h_pho3Et =  new TH1F("h_pho3Et",   "h_pho3Et",200,0,200);
  hPhoCorrEt1 =  new TH1F("hPhoCorrEt1",   "hPhoCorrEt1",200,0,200);
  hPhoCorrEt2 =  new TH1F("hPhoCorrEt2",   "hPhoCorrEt2",200,0,200);
  hPhoCorrEt1_WH =  new TH1F("hPhoCorrEt1_WH",   "hPhoCorrEt1_WH",200,0,200);
  hPhoCorrEt2_WH =  new TH1F("hPhoCorrEt2_WH",   "hPhoCorrEt2_WH",200,0,200);
  hPhoCorrEt1_ZH =  new TH1F("hPhoCorrEt1_ZH",   "hPhoCorrEt1_ZH",200,0,200);
  hPhoCorrEt2_ZH =  new TH1F("hPhoCorrEt2_ZH",   "hPhoCorrEt2_ZH",200,0,200);

  //=========================== lepton tag==============================
  hDPhi_gele_g = new TH1F("hDPhi_gele_g", "hDPhi_gele_g", 70, -3.5, 3.5);

  hDPhi_geleg_gele = new TH2F("hDPhi_geleg_gele", "hDPhi_geleg_gele", 70, -3.5, 3.5, 70,-3.5, 3.5);
  
  hcat4 = new TH1F("hcat4","hcat4",4,0,4);
  hDMZmin = new TH1F("hDMZmin","hDMZmin",500,0,50);
  hDRtoTrk = new TH1F("hDRtoTrk","hDRtoTrk", 300,0,30);
  hMassEleGamGam   = new TH1F("hMassEleGamGam", "M_{#gamma#gamma e}", 500, 0, 500);
  hMassMuGamGam   = new TH1F("hMassMuGamGam", "M_{#gamma#gamma #mu}", 500, 0, 500);
  hPtEleGamGam   = new TH1F("hPtEleGamGam", "p_{T}(#gamma#gamma e)", 500, 0, 500);
  hPtMuGamGam   = new TH1F("hPtMuGamGam", "p_{t}(#gamma#gamma #mu)", 500, 0, 500);


  hggMass_eleVeto = new TH1F("hggMass_eleVeto", "M_{#gamma#gamma}", 1000, 0, 500);
  hggMass_eleVeto_jetVeto = new TH1F("hggMass_eleVeto_jetVeto", "M_{#gamma#gamma}", 1000, 0, 500);

  hggMass_muTag     = new TH1F("hggMass_muTag", "M_{#gamma#gamma}", 1000, 0, 500);
  hggMass_eleTag     = new TH1F("hggMass_eleTag", "M_{#gamma#gamma}", 1000, 0, 500);
  hggMass_eletag_Zveto = new TH1F("hggMass_eleTag_Zveto", "M_{#gamma#gamma} Zveto", 1000, 0, 500);


  hggMassWH_eleTag_Zveto = new TH1F("hggMassWH_eleTag_Zveto", "M_{#gamma#gamma} Zveto", 1000, 0, 500);
  hggMassZH_eleTag_Zveto = new TH1F("hggMassZH_eleTag_Zveto", "M_{#gamma#gamma} Zveto", 1000, 0, 500);

  hggMassWH_muTag     = new TH1F("hggMassWH_muTag", "M_{#gamma#gamma}", 1000, 0,500);
  hggMassWH_eleTag     = new TH1F("hggMassWH_eleTag", "M_{#gamma#gamma}", 1000, 0, 500);
  hggMassZH_muTag     = new TH1F("hggMassZH_muTag", "M_{#gamma#gamma}", 1000, 0, 500);
  hggMassZH_eleTag     = new TH1F("hggMassZH_eleTag", "M_{#gamma#gamma}", 1000, 0, 500);
  hDR_Pho1Ele     = new TH1F("hDR_Pho1Ele", "#DeltaR(#gamma_{1},e)", 60, 0,6);
  hDR_Pho2Ele     = new TH1F("hDR_Pho2Ele", "#DeltaR(#gamma_{2},e)", 60, 0, 6);
  hDR_SCPho1Ele     = new TH1F("hDR_SCPho1Ele", "#DeltaR(#gamma_{1}(SC),e(SC))", 60, 0,6);
  hDR_SCPho2Ele     = new TH1F("hDR_SCPho2Ele", "#DeltaR(#gamma_{2}(SC),e(SC))", 60, 0, 6);
  hDRmin_PhoMu     = new TH1F("hDRmin_PhoMu", "#DeltaR_{min}(#gamma,#mu)", 60, 0, 6);
  hDRmin_PhoEle     = new TH1F("hDRmin_PhoEle", "#DeltaR_{min}(#gamma,e) respect to teh selected vertex", 60, 0, 6);
  hDRmin_SCPhoEle     = new TH1F("hDRmin_SCPhoEle", "#DeltaR_{min}(#gamma(SC),e(SC))", 60, 0, 6);
  hAllElePt  = new TH1F("hAllElePt", "p_{T}(e)", 200, 0, 200);
  hAllSCEleEta  = new TH1F("hAllSCEleEta", "#eta(e)", 60, -3, 3);
  hEleCombRelIso_EB =  new TH1F("hEleCombRelIso_EB", "CombRelIso(e) EB", 200, 0, 0.2);
  hEleCombRelIso_EE =  new TH1F("hEleCombRelIso_EE", "CombRelIso(e) EE", 200, 0, 0.2);

 /*  hPFRelIso_BID_B =  new TH1F("hPFRelIso_BID_B", "CombRelPFIso(e) EB", 200, 0, 0.2); */
/*   hPFRelIso_AID_B =  new TH1F("hPFRelIso_AID_B", "CombRelPFIso(e) EB", 200, 0, 0.2); */
/*   hPFRelIso_BID_EC =  new TH1F("hPFRelIso_BID_EC", "CombRelPFIso(e) EC", 200, 0, 0.2); */
/*   hPFRelIso_AID_EC =  new TH1F("hPFRelIso_AID_EC", "CombRelPFIso(e) EC", 200, 0, 0.2); */


  helePFRelIso_BID_B =  new TH1F("helePFRelIso_BID_B", "CombRelPFIso(e) EB", 200, 0, 0.2);
  helePFRelIso_AID_B =  new TH1F("helePFRelIso_AID_B", "CombRelPFIso(e) EB", 200, 0, 0.2);
  helePFRelIso_BID_EC =  new TH1F("helePFRelIso_BID_EC", "CombRelPFIso(e) EC", 200, 0, 0.2);
  helePFRelIso_AID_EC =  new TH1F("helePFRelIso_AID_EC", "CombRelPFIso(e) EC", 200, 0, 0.2);
  
  hmuPFRelIso_BID =  new TH1F("hmuPFRelIso_BID", "CombRelPFIso(#mu) EC", 200, 0, 0.2);
  hmuPFRelIso_AID =  new TH1F("hmuPFRelIso_AID", "CombRelPFIso(#mu) EC", 200, 0, 0.2);


  hSelEle1Pt  = new TH1F("hSelEle1Pt", "p_{T}(e)", 200, 0, 200);
  hSelSCEle1Eta  = new TH1F("hSelSCEle1Eta", "#eta(e)", 60, -3, 3);
  hSelEle2Pt  = new TH1F("hSelEle2Pt", "p_{T}(e)", 200, 0, 200);
  hSelSCEle2Eta  = new TH1F("hSelSCEle2Eta", "#eta(e)", 60, -3, 3);

  hAllMuPt  = new TH1F("hAllMuPt", "p_{T}(#mu)", 200, 0, 200);
  hAllMuEta  = new TH1F("hAllMuEta", "#eta(#mu)", 60, -3, 3);
  hMuCombRelIso =  new TH1F("hMuCombRelIso", "CombRelIso(#mu)" , 200, 0, 2);

  hSelMu1Pt  = new TH1F("hSelMu1Pt", "p_{T}(#mu)", 200, 0, 200);
  hSelMu1Eta  = new TH1F("hSelMu1Eta", "#eta(#mu)", 60, -3, 3);
  hSelMu2Pt  = new TH1F("hSelMu2Pt", "p_{T}(#mu)", 200, 0, 200);
  hSelMu2Eta  = new TH1F("hSelMu2Eta", "#eta(#mu)", 60, -3, 3);


  hMassEleEleGam1    = new TH1F("hMassEleEleGam1", "M(ee#gamma_{1})", 500, 0, 500);
  hMassEleEleGam2    = new TH1F("hMassEleEleGam2", "M(ee#gamma_{2})", 1000, 0, 500);
  hMassEleGam1    = new TH1F("hMassEleGam1", "M(e#gamma_{1})", 500, 0, 500);
  hMassEleGam2    = new TH1F("hMassEleGam2", "M(e#gamma_{2})", 500, 0, 500);


  hMassEleFakePho = new TH1F("hMassEleFakePho", "M(e gamma_{fake}) (GeV/c^{2})", 500, 0, 500)  ;
  hMassEleRealPho = new TH1F("hMassEleRealPho", "M(e gamma_{real}) (GeV/c^{2})", 500, 0, 500);

  hPtEleEleGam1    = new TH1F("hPtEleEleGam1", "p_{T}(ee#gamma_{1})", 500, 0, 500);
  hPtEleEleGam2    = new TH1F("hPtEleEleGam2", "p_{T}(ee#gamma_{2})", 1000, 0, 500);
  hPtEleGam1    = new TH1F("hPtEleGam1", "p_{T}(e#gamma_{1})", 500, 0, 500);
  hPtEleGam2    = new TH1F("hPtEleGam2", "p_{T}(e#gamma_{2})", 500, 0, 500);






  h_pfMET_lepTag = new TH1F("h_pfMET_lepTag", "PFMET lepTag", 200, 0, 200);
  h_tcMET_lepTag = new TH1F("h_tcMET_lepTag", "TCMET lepTag", 200, 0, 200);
  h_CaloMET_lepTag = new TH1F("h_CaloMET_lepTag", "Calo MET lepTag", 200, 0, 200);


  h_pfMET = new TH1F("h_pfMET", "PFMET ", 200, 0, 200);
  h_tcMET = new TH1F("h_tcMET", "TCMET ", 200, 0, 200);
  h_CaloMET = new TH1F("h_CaloMET", "Calo MET ", 200, 0, 200);
  //================================================
  hggMass     = new TH1F("hggMass", "M_{#gamma#gamma}", 1000, 0, 500);
  //==================================================

  hJetNHF =  new TH1F("hJetNHF", "hJetNHF", 100, 0, 1);
  hJetNEF =  new TH1F("hJetNEF", "hJetNEF", 100, 0, 1);
  hJetNConst =  new TH1F("hJetNConst", "hJetNConst", 10, 0, 10);

  hJetCHF =  new TH1F("hJetCHF", "hJetCHF", 100, 0, 1);
  hJetNCH =  new TH1F("hJetNCH", "hJetNCH", 100, 0, 10);
  hJetCEF =  new TH1F("hJetCEF", "hJetCEF", 100, 0, 1);



  hDRAllJetg1 = new TH1F("hDRAllJetg1", "hDRAllJetg1", 70, 0, 7);
  hDRAllJetg2 = new TH1F("hDRAllJetg2", "hDRAllJetg2", 70, 0, 7);
  hAllJetPt = new TH1F("hAllJetPt", "hAllJetPt", 200, 0, 200);
  hAllJetEta = new TH1F("hAllJetEta", "hAllJetEta", 100, -5, 5);

  hDRjet1g1   = new TH1F("hDRjet1g1", "hDRjet1g1", 70, 0, 7);
  hDRjet1g2   = new TH1F("hDRjet1g2", "hDRjet1g2", 70, 0, 7);
  hDRjet2g1   = new TH1F("hDRjet2g1", "hDRjet2g1", 70, 0, 7);
  hDRjet2g2   = new TH1F("hDRjet2g2", "hDRjet2g2", 70, 0, 7);

  hDPhijet1g1   = new TH1F("hDPhijet1g1", "hDPhijet1g1", 70, -3.5, 3.5);
  hDPhijet1g2   = new TH1F("hDPhijet1g2", "hDPhijet1g2", 70, -3.5, 3.5);
  hDPhijet2g1   = new TH1F("hDPhijet2g1", "hDPhijet2g1",70, -3.5, 3.5);
  hDPhijet2g2   = new TH1F("hDPhijet2g2", "hDPhijet2g2", 70, -3.5, 3.5);
  hDPhijet1jet2   = new TH1F("hDPhijet1jet2", "hDPhijet1jet2", 70, -3.5, 3.5);

  hggMass_twoJetTag     = new TH1F("hggMass_twoJetTag", "M_{#gamma#gamma}", 1000, 0, 1000);
  hggPt_twoJetTag    = new TH1F("hggPt_twoJetTag", "Pt_{#gamma#gamma}", 100, 0, 1000);
  hJet1Pt    = new TH1F("hJet1Pt", "p_{T}(jet1)", 300, 0, 300);
  hJet2Pt    = new TH1F("hJet2Pt", "p_{T}(jet2)", 300, 0, 300);
 
  hJet1Eta    = new TH1F("hJet1Eta", "#eta (jet1)", 100, -5, 5);
  hJet2Eta    = new TH1F("hJet2Eta", "#eta (jet2)", 100, -5, 5);


  hDeltaEtaJets    = new TH1F("hDeltaEtaJets", "#Delta #eta (jet1,jet2)", 80, 0, 8);
  hDeltaEtaJets_bc    = new TH1F("hDeltaEtaJets_bc", "#Delta #eta (jet1,jet2) bc", 80, 0, 8);

 
  hJetEtaProd_bc   = new TH1F("hJetEtaProd_bc", "hJetEtaProd_bc", 3, -1,2 );

  //hDeltaEtaJetsVsMgg    = new TH2F("hDeltaEtaJetsVsMgg", "#Delta #eta (jet1,jet2) vs M(#gamma #gamma)", 1000, 0, 1000, 80, 0, 8);
  //hMassJetJetVsMgg  = new TH2F("hMassJetJetVsMgg", " M(jet1,jet2) vs M(#gamma #gamma)", 1000, 0, 1000, 200, 0, 2000);

  hZeppenfeld = new TH1F("hZeppenfeld", "hZeppenfeld", 200, -10, 10);
  hZeppenfeld_bc = new TH1F("hZeppenfeld_bc", "hZeppenfeld bc", 200, -10, 10);

  hMassJetJet = new TH1F("hMassJetJet", "M(jet_{1},#gamma_{1}) [GeV]", 200, 0, 2000); 
  hMassJetJet_bc = new TH1F("hMassJetJet_bc", "M(jet_{1},#gamma_{1}) [GeV]  bc", 200, 0, 2000); 

  hDeltaPhiDiJetDiPho   = new TH1F("hDeltaPhiDiJetDiPho", "dPhi 2jet 2pho", 70, -3.5, 3.5);
  hDeltaPhiDiJetDiPho_bc = new TH1F("hDeltaPhiDiJetDiPho_bc", "dPhi 2jet 2pho bc", 70, -3.5, 3.5);



  //hMgg_twoJetTag_etaprod   = new TH1F("hMgg_twoJetTag_etaprod",  "M_{#gamma#gamma} #eta_{1} #times #eta_{2}<0", 1000, 0, 1000);
  hMgg_twoJetTag_jet1Pt    = new TH1F("hMgg_twoJetTag_jet1Pt",   "M_{#gamma#gamma} jet1 pt > 30",               1000, 0, 1000);
  hMgg_twoJetTag_deltaEta  = new TH1F("hMgg_twoJetTag_deltaEta", "M_{#gamma#gamma} #Delta #eta > 3.5",          1000, 0, 1000);
  hMgg_twoJetTag_zep       = new TH1F("hMgg_twoJetTag_zep",      "M_{#gamma#gamma} Zeppenfeld < 2.5",           1000, 0, 1000);
  hMgg_twoJetTag_Mjets     = new TH1F("hMgg_twoJetTag_Mjets",    "M_{#gamma#gamma} M(jet1,jet2) > 350",         1000, 0, 1000);
  hMgg_twoJetTag_dPhi      = new TH1F("hMgg_twoJetTag_dPhi",     "M_{#gamma#gamma} dPhi(jj,gg) > 2.6",          1000, 0, 1000);


  hJet1Pt_HSTRA_bc               = new TH1F("hJet1Pt_HSTRA_bc", "p_{T}(jet1) bc", 300, 0, 300);
  hDeltaEtaJets_HSTRA_bc         = new TH1F("hDeltaEtaJets_HSTRA_bc", "#Delta #eta (jet1,jet2) bc", 80, 0, 8);
  hZeppenfeld_HSTRA_bc           = new TH1F("hZeppenfeld_HSTRA_bc",          "hZeppenfeld bc",                            200, -10, 10);
  hMassJetJet_HSTRA_bc           = new TH1F("hMassJetJet_HSTRA_bc",          "M(jet_{1},#gamma_{1}) [GeV]  bc",           200, 0, 2000); 

  hMgg_HSTRA_twoJetTag           = new TH1F("hMgg_HSTRA_twoJetTag",          "M_{#gamma#gamma} pho1 pt > 60",            1000, 0, 1000);
  hMgg_HSTRA_twoJetTag_jet1Pt    = new TH1F("hMgg_HSTRA_twoJetTag_jet1Pt",   "M_{#gamma#gamma} jet1 pt > 30",            1000, 0, 1000);
  hMgg_HSTRA_twoJetTag_deltaEta  = new TH1F("hMgg_HSTRA_twoJetTag_deltaEta", "M_{#gamma#gamma} #Delta #eta < 2.5",       1000, 0, 1000);
  hMgg_HSTRA_twoJetTag_zep       = new TH1F("hMgg_HSTRA_twoJetTag_zep",      "M_{#gamma#gamma} Zeppenfeld < 1.5",        1000, 0, 1000);
  hMgg_HSTRA_twoJetTag_Mjets     = new TH1F("hMgg_HSTRA_twoJetTag_Mjets",    "M_{#gamma#gamma} 55 < M(jet1,jet2) < 115", 1000, 0, 1000);


}
// JM: add isItData to disable MC branches when data

void xAna::Init(TTree *tree, bool isItData) {
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   vtx_x = 0;
   vtx_y = 0;
   vtx_z = 0;
   vtxbs_x = 0;
   vtxbs_y = 0;
   vtxbs_z = 0;
   vtxbsPtMod = 0;
   vtxbsSumPt2 = 0;
   vtxbsTkIndex = 0;
   vtxbsTkWeight = 0;
   trkP_x = 0;
   trkP_y = 0;
   trkP_z = 0;
   trkVtx_x = 0;
   trkVtx_y = 0;
   trkVtx_z = 0;
   trkd0 = 0;
   trkd0Err = 0;
   trkdz = 0;
   trkdzErr = 0;
   trkPtErr = 0;
   trkQuality = 0;
   mcHardPt = 0;
   mcHardEta = 0;
   mcHardPhi = 0;
   mcHardM = 0;
   mcHardPID = 0;
   mcHardFun = 0;

   //CF
   mcHardOutPt = 0; 
   mcHardOutP = 0;
   mcHardOutEta = 0;
   mcHardOutPhi = 0;
   mcHardOutPdgId = 0;
   mcHardOutE = 0;

   mcPID = 0;
   mcVtx_x = 0;
   mcVtx_y = 0;
   mcVtx_z = 0;
   mcPt = 0;
   mcMass = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcGMomPID = 0;
   mcMomPID = 0;
   mcMomPt = 0;
   mcMomMass = 0;
   mcMomEta = 0;
   mcMomPhi = 0;
   mcIndex = 0;
   mcDecayType = 0;
   mcParentage = 0;
   mcStatus = 0;
   nPU = 0;
   puBX = 0;
   puTrue = 0;
   trkMETx = 0;
   trkMETy = 0;
   trkMETPhi = 0;
   trkMET = 0;
   eleTrg = 0;
   eleClass = 0;
   eleIsEcalDriven = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleEn = 0;
   eleEcalEn = 0;
   eleSCRawEn = 0;
   eleSCEn = 0;
   eleESEn = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9 = 0;
   eleEtaVtx = 0;
   elePhiVtx = 0;
   eleEtVtx = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleVtx_x = 0;
   eleVtx_y = 0;
   eleVtx_z = 0;
   eleD0 = 0;
   eleDz = 0;
   eleD0GV = 0;
   eleDzGV = 0;
   eleD0Vtx = 0;
   eleDzVtx = 0;
   eleHoverE = 0;
   eleHoverE12 = 0;
   eleEoverP = 0;
   elePin = 0;
   elePout = 0;
   eleTrkMomErr = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eleSigmaIEtaIEta = 0;
   eleSigmaIEtaIPhi = 0;
   eleSigmaIPhiIPhi = 0;
   eleEmax = 0;
   eleE2ndMax = 0;
   eleETop = 0;
   eleEBottom = 0;
   eleELeft = 0;
   eleERight = 0;
   eleE1x5 = 0;
   eleE3x3 = 0;
   eleE5x5 = 0;
   eleE2x5Max = 0;
   eleE2x5Top = 0;
   eleE2x5Bottom = 0;
   eleE2x5Left = 0;
   eleE2x5Right = 0;
   eleSeedEta = 0;
   eleSeedE = 0;
   eleSeedPhi = 0;
   eleCrysEta = 0;
   eleCrysPhi = 0;
   eleCrysIEta = 0;
   eleCrysIPhi = 0;
   eleRegrE = 0;
   eleRegrEerr = 0;
   elePhoRegrE = 0;
   elePhoRegrEerr = 0;
   eleSeedTime = 0;
   eleRecoFlag = 0;
   elePos = 0;
   eleGenIndex = 0;
   eleGenGMomPID = 0;
   eleGenMomPID = 0;
   eleGenMomPt = 0;
   eleIsoTrkDR03 = 0;
   eleIsoEcalDR03 = 0;
   eleIsoHcalDR03 = 0;
   eleIsoHcalDR0312 = 0;
   eleIsoTrkDR04 = 0;
   eleIsoEcalDR04 = 0;
   eleIsoHcalDR04 = 0;
   eleIsoHcalDR0412 = 0;
   eleModIsoTrk = 0;
   eleModIsoEcal = 0;
   eleModIsoHcal = 0;
   eleMissHits = 0;
   eleConvDist = 0;
   eleConvDcot = 0;
   eleConvVtxFit = 0;
   eleIP3D = 0;
   eleIP3DErr = 0;
   eleIDMVANonTrig = 0;
   eleIDMVATrig = 0;
   elePFChIso03 = 0;
   elePFPhoIso03 = 0;
   elePFNeuIso03 = 0;
   elePFChIso04 = 0;
   elePFPhoIso04 = 0;
   elePFNeuIso04 = 0;
   eleESEffSigmaRR_x = 0;
   eleESEffSigmaRR_y = 0;
   eleESEffSigmaRR_z = 0;
   phoTrg = 0;
   phoTrgFilter = 0;
   phoIsPhoton = 0;
   phoSCPos_x = 0;
   phoSCPos_y = 0;
   phoSCPos_z = 0;
   phoCaloPos_x = 0;
   phoCaloPos_y = 0;
   phoCaloPos_z = 0;
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoVtx_x = 0;
   phoVtx_y = 0;
   phoVtx_z = 0;
   phoPhi = 0;
   phoEtVtx = 0;
   phoEtaVtx = 0;
   phoPhiVtx = 0;
   phoR9 = 0;
   phoNClus = 0;
   phoTrkIsoHollowDR03 = 0;
   phoEcalIsoDR03 = 0;
   phoHcalIsoDR03 = 0;
   phoHcalIsoDR0312 = 0;
   phoTrkIsoHollowDR04 = 0;
   phoCiCdRtoTrk = 0;
   phoEcalIsoDR04 = 0;
   phoHcalIsoDR04 = 0;
   phoHcalIsoDR0412 = 0;
   phoHoverE = 0;
   phoHoverE12 = 0;
   phoEleVeto = 0;
   phoSigmaIEtaIEta = 0;
   phoSigmaIEtaIPhi = 0;
   phoSigmaIPhiIPhi = 0;
   phoCiCPF4phopfIso03 = 0;
   phoCiCPF4phopfIso04 = 0;
   phoCiCPF4chgpfIso02 = 0;
   phoCiCPF4chgpfIso03 = 0;
   phoCiCPF4chgpfIso04 = 0;
   phoEmax = 0;
   phoETop = 0;
   phoEBottom = 0;
   phoELeft = 0;
   phoERight = 0;
   phoE2ndMax = 0;
   phoE3x3 = 0;
   phoE3x1 = 0;
   phoE1x3 = 0;
   phoE5x5 = 0;
   phoE1x5 = 0;
   phoE2x2 = 0;
   phoE2x5Max= 0;
   phoE2x5Top = 0;
   phoE2x5Bottom = 0;
   phoE2x5Left = 0;
   phoE2x5Right = 0;
   phoSeedE = 0;
   phoSeedEta = 0;
   phoSeedPhi = 0;
   phoCrysEta = 0;
   phoCrysPhi = 0;
   phoCrysIEta = 0;
   phoCrysIPhi = 0;
   phoPFChIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoSCRChIso = 0;
   phoSCRPhoIso = 0;
   phoSCRNeuIso = 0;
   phoSCRChIso04 = 0;
   phoSCRPhoIso04 = 0;
   phoSCRNeuIso04 = 0;
   phoRandConeChIso = 0;
   phoRandConePhoIso = 0;
   phoRandConeNeuIso = 0;
   phoRandConeChIso04 = 0;
   phoRandConePhoIso04 = 0;
   phoRandConeNeuIso04 = 0;
   phoRegrE = 0;
   phoRegrEerr = 0;
   phoSeedTime = 0;
   phoSeedDetId1 = 0;
   phoSeedDetId2 = 0;
   phoLICTD = 0;
   phoRecoFlag = 0;
   phoPos = 0;
   phoGenIndex = 0;
   phoGenGMomPID = 0;
   phoGenMomPID = 0;
   phoGenMomPt = 0;
   phoSCE = 0;
   phoSCRawE = 0;
   phoESEn = 0;
   phoSCEt = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoSCBrem = 0;
   phoOverlap = 0;
   phohasPixelSeed = 0;
   pho_hasConvPf = 0;
   pho_hasSLConvPf = 0;
   pho_pfconvVtxZ = 0;
   pho_pfconvVtxZErr = 0;
   pho_nSLConv = 0;
   pho_pfSLConvPos_x = 0;
   pho_pfSLConvPos_y = 0;
   pho_pfSLConvPos_z = 0;
   pho_pfSLConvVtxZ = 0;
   phoIsConv = 0;
   phoNConv = 0;
   phoConvInvMass = 0;
   phoConvCotTheta = 0;
   phoConvEoverP = 0;
   phoConvZofPVfromTrks = 0;
   phoConvMinDist = 0;
   phoConvdPhiAtVtx = 0;
   phoConvdPhiAtCalo = 0;
   phoConvdEtaAtCalo = 0;
   phoConvTrkd0_x = 0;
   phoConvTrkd0_y = 0;
   phoConvTrkPin_x = 0;
   phoConvTrkPin_y = 0;
   phoConvTrkPout_x = 0;
   phoConvTrkPout_y = 0;
   phoConvTrkdz_x = 0;
   phoConvTrkdz_y = 0;
   phoConvTrkdzErr_x = 0;
   phoConvTrkdzErr_y = 0;
   phoConvChi2 = 0;
   phoConvChi2Prob = 0;
   phoConvNTrks = 0;
   phoConvCharge1 = 0;
   phoConvCharge2 = 0;
   phoConvValidVtx = 0;
   phoConvLikeLihood = 0;
   phoConvP4_0 = 0;
   phoConvP4_1 = 0;
   phoConvP4_2 = 0;
   phoConvP4_3 = 0;
   phoConvVtx_x = 0;
   phoConvVtx_y = 0;
   phoConvVtx_z = 0;
   phoConvVtxErr_x = 0;
   phoConvVtxErr_y = 0;
   phoConvVtxErr_z = 0;
   phoConvPairMomentum_x = 0;
   phoConvPairMomentum_y = 0;
   phoConvPairMomentum_z = 0;
   phoConvRefittedMomentum_x = 0;
   phoConvRefittedMomentum_y = 0;
   phoConvRefittedMomentum_z = 0;
   SingleLegConv = 0;
   phoPFConvVtx_x = 0;
   phoPFConvVtx_y = 0;
   phoPFConvVtx_z = 0;
   phoPFConvMom_x = 0;
   phoPFConvMom_y = 0;
   phoPFConvMom_z = 0;
   phoESEffSigmaRR_x = 0;
   phoESEffSigmaRR_y = 0;
   phoESEffSigmaRR_z = 0;
   muTrg = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muPt = 0;
   muPz = 0;
   muVtx_x = 0;
   muVtx_y = 0;
   muVtx_z = 0;
   muVtxGlb_x = 0;
   muVtxGlb_y = 0;
   muVtxGlb_z = 0;
   muGenIndex = 0;
   mucktPt = 0;
   mucktPtErr = 0;
   mucktEta = 0;
   mucktPhi = 0;
   mucktdxy = 0;
   mucktdz = 0;
   muIsoTrk = 0;
   muIsoCalo = 0;
   muIsoEcal = 0;
   muIsoHcal = 0;
   muChi2NDF = 0;
   muInnerChi2NDF = 0;
   muPFIsoR04_CH = 0;
   muPFIsoR04_NH = 0;
   muPFIsoR04_Pho = 0;
   muPFIsoR04_PU = 0;
   muPFIsoR04_CPart = 0;
   muPFIsoR04_NHHT = 0;
   muPFIsoR04_PhoHT = 0;
   muPFIsoR03_CH = 0;
   muPFIsoR03_NH = 0;
   muPFIsoR03_Pho = 0;
   muPFIsoR03_PU = 0;
   muPFIsoR03_CPart = 0;
   muPFIsoR03_NHHT = 0;
   muPFIsoR03_PhoHT = 0;
   muType = 0;
   muD0 = 0;
   muDz = 0;
   muD0GV = 0;
   muDzGV = 0;
   muD0Vtx = 0;
   muDzVtx = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muInnerD0GV = 0;
   muInnerDzGV = 0;
   muInnerPt = 0;
   muInnerPtErr = 0;
   muNumberOfValidTrkLayers = 0;
   muNumberOfValidTrkHits = 0;
   muNumberOfValidPixelLayers = 0;
   muNumberOfValidPixelHits = 0;
   muNumberOfValidMuonHits = 0;
   muStations = 0;
   muChambers = 0;
   muIP3D = 0;
   muIP3DErr = 0;
   tauDecayModeFinding = 0;
   tauAgainstElectronLooseMVA3 = 0;
   tauAgainstElectronMediumMVA3 = 0;
   tauAgainstElectronTightMVA3 = 0;
   tauAgainstElectronVTightMVA3 = 0;
   tauAgainstElectronDeadECAL = 0;
   tauAgainstMuonLoose2 = 0;
   tauAgainstMuonMedium2 = 0;
   tauAgainstMuonTight2 = 0;
   tauCombinedIsolationDeltaBetaCorrRaw3Hits = 0;
   tauLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   tauMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
   tauTightCombinedIsolationDeltaBetaCorr3Hits = 0;
   tauEta = 0;
   tauPhi = 0;
   tauPt = 0;
   tauEt = 0;
   tauCharge = 0;
   tauDecayMode = 0;
   tauEMFraction = 0;
   tauHCAL3x3OverPLead = 0;
   tauHCALMaxOverPLead = 0;
   tauHCALTotOverPLead = 0;
   tauIsolationPFChargedHadrCandsPtSum = 0;
   tauIsolationPFGammaCandsEtSum = 0;
   tauLeadPFChargedHadrCandsignedSipt = 0;
   tauLeadChargedHadronExists = 0;
   tauLeadChargedHadronEta = 0;
   tauLeadChargedHadronPhi = 0;
   tauLeadChargedHadronPt = 0;
   CA8JetPt = 0;
   CA8JetEta = 0;
   CA8JetPhi = 0;
   CA8JetMass = 0;
   CA8JetArea = 0;
   CA8Jet_tau1 = 0;
   CA8Jet_tau2 = 0;
   CA8Jet_tau3 = 0;
   CA8JetCHF = 0;
   CA8JetNHF = 0;
   CA8JetCEF = 0;
   CA8JetNEF = 0;
   CA8JetNCH = 0;
   CA8Jetnconstituents = 0;
   CA8prunedJetMass = 0;
   CA8prunedJet_nSubJets = 0;
   CA8prunedJet_SubjetPt = 0;
   CA8prunedJet_SubjetEta = 0;
   CA8prunedJet_SubjetPhi = 0;
   CA8prunedJet_SubjetMass = 0;
   jetTrg = 0;
   jetEn = 0;
   jetPt = 0;
   jetEta = 0;
   jetPhi = 0;
   jetCharge = 0;
   jetEt = 0;
   jetRawPt = 0;
   jetRawEn = 0;
   jetArea = 0;
   jetCHF = 0;
   jetNHF = 0;
   jetCEF = 0;
   jetNEF = 0;
   jetNCH = 0;
   jetHFHAE = 0;
   jetHFEME = 0;
   jetNConstituents = 0;
   jetCombinedSecondaryVtxBJetTags = 0;
   jetCombinedSecondaryVtxMVABJetTags = 0;
   jetJetProbabilityBJetTags = 0;
   jetJetBProbabilityBJetTags = 0;
   jetBetaStar = 0;
   jetPFLooseId = 0;
   jetDRMean = 0;
   jetDR2Mean = 0;
   jetDZ = 0;
   jetFrac01 = 0;
   jetFrac02 = 0;
   jetFrac03 = 0;
   jetFrac04 = 0;
   jetFrac05 = 0;
   jetFrac06 = 0;
   jetFrac07 = 0;
   jetBeta = 0;
   jetBetaStarCMG = 0;
   jetBetaStarClassic = 0;
   jetBetaExt = 0;
   jetBetaStarCMGExt = 0;
   jetBetaStarClassicExt = 0;
   jetNNeutrals = 0;
   jetNCharged = 0;
   jetMVAs = 0;
   jetWPLevels = 0;
   jetMVAsExt_simple = 0;
   jetWPLevelsExt_simple = 0;
   jetMVAsExt_full = 0;
   jetWPLevelsExt_full = 0;
   jetMVAsExt_cutBased = 0;
   jetWPLevelsExt_cutBased = 0;
   jetMVAsExt_philv1 = 0;
   jetWPLevelsExt_philv1 = 0;
   jetMt = 0;
   jetJECUnc = 0;
   jetLeadTrackPt = 0;
   jetVtxPt = 0;
   jetVtxMass = 0;
   jetVtx3dL = 0;
   jetVtx3deL = 0;
   jetSoftLeptPt = 0;
   jetSoftLeptPtRel = 0;
   jetSoftLeptdR = 0;
   jetSoftLeptIdlooseMu = 0;
   jetSoftLeptIdEle95 = 0;
   jetDPhiMETJet = 0;
   jetPuJetIdL = 0;
   jetPuJetIdM = 0;
   jetPuJetIdT = 0;
   jetPartonID = 0;
   jetGenJetIndex = 0;
   jetGenJetEn = 0;
   jetGenJetPt = 0;
   jetGenJetEta = 0;
   jetGenJetPhi = 0;
   jetGenPartonID = 0;
   jetGenEn = 0;
   jetGenPt = 0;
   jetGenEta = 0;
   jetGenPhi = 0;
   jetGenPartonMomID = 0;
   jetLowPtEn = 0;
   jetLowPtPt = 0;
   jetLowPtEta = 0;
   jetLowPtPhi = 0;
   jetLowPtCharge = 0;
   jetLowPtEt = 0;
   jetLowPtRawPt = 0;
   jetLowPtRawEn = 0;
   jetLowPtArea = 0;
   jetLowPtPartonID = 0;
   jetLowPtGenJetEn = 0;
   jetLowPtGenJetPt = 0;
   jetLowPtGenJetEta = 0;
   jetLowPtGenJetPhi = 0;
   jetLowPtGenPartonID = 0;
   jetLowPtGenEn = 0;
   jetLowPtGenPt = 0;
   jetLowPtGenEta = 0;
   jetLowPtGenPhi = 0;
   convP4_x = 0;
   convP4_y = 0;
   convP4_z = 0;
   convP4_E = 0;
   convVtx_x = 0;
   convVtx_y = 0;
   convVtx_z = 0;
   convVtxErr_x = 0;
   convVtxErr_y = 0;
   convVtxErr_z = 0;
   convPairMomentum_x = 0;
   convPairMomentum_y = 0;
   convPairMomentum_z = 0;
   convRefittedMomentum_x = 0;
   convRefittedMomentum_y = 0;
   convRefittedMomentum_z = 0;
   convNTracks = 0;
   convPairInvMass = 0;
   convPairCotThetaSep = 0;
   convEoverP = 0;
   convDistOfMinApproach = 0;
   convDPhiTrksAtVtx = 0;
   convDPhiTrksAtEcal = 0;
   convDEtaTrksAtEcal = 0;
   convDxy = 0;
   convDz = 0;
   convLxy = 0;
   convLz = 0;
   convZofPrimVtxFromTrks = 0;
   convNHitsBeforeVtx_0 = 0;
   convNHitsBeforeVtx_1 = 0;
   convNSharedHits = 0;
   convValidVtx = 0;
   convMVALikelihood = 0;
   convChi2 = 0;
   convChi2Probability = 0;
   convTk1Dz = 0;
   convTk2Dz = 0;
   convTk1DzErr = 0;
   convTk2DzErr = 0;
   convCh1Ch2 = 0;
   convTk1D0 = 0;
   convTk1Pout = 0;
   convTk1Pin = 0;
   convTk2D0 = 0;
   convTk2Pout = 0;
   convTk2Pin = 0;

   /// non automatic root variables
   DiscriVBF_UseDiPhoPt = true ;
   DiscriVBF_UsePhoPt   = true ;
   DiscriVBF_useMvaSel  = false;
   DiscriVBF_useCombMvaSel  = false;
   vtxbsTkIndex = 0;
   //vtxbsTkWeight = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("HLT", HLT, &b_HLT);
   fChain->SetBranchAddress("HLTIndex", HLTIndex, &b_HLTIndex);
   fChain->SetBranchAddress("bspotPos", bspotPos, &b_bspotPos);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("IsVtxGood", &IsVtxGood, &b_IsVtxGood);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("nVtxBS", &nVtxBS, &b_nVtxBS);
   fChain->SetBranchAddress("vtxbs_x", &vtxbs_x, &b_vtxbs_x);
   fChain->SetBranchAddress("vtxbs_y", &vtxbs_y, &b_vtxbs_y);
   fChain->SetBranchAddress("vtxbs_z", &vtxbs_z, &b_vtxbs_z);
   fChain->SetBranchAddress("vtxbsPtMod", &vtxbsPtMod, &b_vtxbsPtMod);
   fChain->SetBranchAddress("vtxbsSumPt2", &vtxbsSumPt2, &b_vtxbsSumPt2);
   fChain->SetBranchAddress("vtxbsTkIndex", &vtxbsTkIndex, &b_vtxbsTkIndex);
   fChain->SetBranchAddress("vtxbsTkWeight", &vtxbsTkWeight, &b_vtxbsTkWeight);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("trkP_x", &trkP_x, &b_trkP_x);
   fChain->SetBranchAddress("trkP_y", &trkP_y, &b_trkP_y);
   fChain->SetBranchAddress("trkP_z", &trkP_z, &b_trkP_z);
   fChain->SetBranchAddress("trkVtx_x", &trkVtx_x, &b_trkVtx_x);
   fChain->SetBranchAddress("trkVtx_y", &trkVtx_y, &b_trkVtx_y);
   fChain->SetBranchAddress("trkVtx_z", &trkVtx_z, &b_trkVtx_z);
   fChain->SetBranchAddress("trkd0", &trkd0, &b_trkd0);
   fChain->SetBranchAddress("trkd0Err", &trkd0Err, &b_trkd0Err);
   fChain->SetBranchAddress("trkdz", &trkdz, &b_trkdz);
   fChain->SetBranchAddress("trkdzErr", &trkdzErr, &b_trkdzErr);
   fChain->SetBranchAddress("trkPtErr", &trkPtErr, &b_trkPtErr);
   fChain->SetBranchAddress("trkQuality", &trkQuality, &b_trkQuality);
   fChain->SetBranchAddress("nGoodTrk", &nGoodTrk, &b_nGoodTrk);
   fChain->SetBranchAddress("IsTracksGood", &IsTracksGood, &b_IsTracksGood);
   //JM
   if(!isItData){
     fChain->SetBranchAddress("pdf", pdf, &b_pdf);
     fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
     fChain->SetBranchAddress("processID", &processID, &b_processID);   
     fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
     fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
     fChain->SetBranchAddress("mcVtx_x", &mcVtx_x, &b_mcVtx_x);
     fChain->SetBranchAddress("mcVtx_y", &mcVtx_y, &b_mcVtx_y);
     fChain->SetBranchAddress("mcVtx_z", &mcVtx_z, &b_mcVtx_z);
     fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
     fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
     fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
     fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
     fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
     fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
     fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
     fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
     fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);
     fChain->SetBranchAddress("mcMomMass", &mcMomMass, &b_mcMomMass);
     fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);
     fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);
     fChain->SetBranchAddress("mcIndex", &mcIndex, &b_mcIndex);
     fChain->SetBranchAddress("mcDecayType", &mcDecayType, &b_mcDecayType);
     fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);
     fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
     fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
     fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
     fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
     fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
     fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
     fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   }
 
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("pfType01MET", &pfType01MET, &b_pfType01MET);
   fChain->SetBranchAddress("pfType01METPhi", &pfType01METPhi, &b_pfType01METPhi);
   fChain->SetBranchAddress("pfType01METsumEt", &pfType01METsumEt, &b_pfType01METsumEt);
   fChain->SetBranchAddress("pfType01METmEtSig", &pfType01METmEtSig, &b_pfType01METmEtSig);
   fChain->SetBranchAddress("pfType01METSig", &pfType01METSig, &b_pfType01METSig);
   fChain->SetBranchAddress("recoPfMET", &recoPfMET, &b_recoPfMET);
   fChain->SetBranchAddress("recoPfMETPhi", &recoPfMETPhi, &b_recoPfMETPhi);
   fChain->SetBranchAddress("recoPfMETsumEt", &recoPfMETsumEt, &b_recoPfMETsumEt);
   fChain->SetBranchAddress("recoPfMETmEtSig", &recoPfMETmEtSig, &b_recoPfMETmEtSig);
   fChain->SetBranchAddress("recoPfMETSig", &recoPfMETSig, &b_recoPfMETSig);
   fChain->SetBranchAddress("trkMETxPV", &trkMETxPV, &b_trkMETxPV);
   fChain->SetBranchAddress("trkMETyPV", &trkMETyPV, &b_trkMETyPV);
   fChain->SetBranchAddress("trkMETPhiPV", &trkMETPhiPV, &b_trkMETPhiPV);
   fChain->SetBranchAddress("trkMETPV", &trkMETPV, &b_trkMETPV);
   fChain->SetBranchAddress("trkMETx", &trkMETx, &b_trkMETx);
   fChain->SetBranchAddress("trkMETy", &trkMETy, &b_trkMETy);
   fChain->SetBranchAddress("trkMETPhi", &trkMETPhi, &b_trkMETPhi);
   fChain->SetBranchAddress("trkMET", &trkMET, &b_trkMET);
   fChain->SetBranchAddress("metFilters", metFilters, &b_metFilters);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleTrg", &eleTrg, &b_eleTrg);
   fChain->SetBranchAddress("eleClass", &eleClass, &b_eleClass);
   fChain->SetBranchAddress("eleIsEcalDriven", &eleIsEcalDriven, &b_eleIsEcalDriven);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleEcalEn", &eleEcalEn, &b_eleEcalEn);
   fChain->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleESEn", &eleESEn, &b_eleESEn);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleEtaVtx", &eleEtaVtx, &b_eleEtaVtx);
   fChain->SetBranchAddress("elePhiVtx", &elePhiVtx, &b_elePhiVtx);
   fChain->SetBranchAddress("eleEtVtx", &eleEtVtx, &b_eleEtVtx);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleVtx_x", &eleVtx_x, &b_eleVtx_x);
   fChain->SetBranchAddress("eleVtx_y", &eleVtx_y, &b_eleVtx_y);
   fChain->SetBranchAddress("eleVtx_z", &eleVtx_z, &b_eleVtx_z);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleD0GV", &eleD0GV, &b_eleD0GV);
   fChain->SetBranchAddress("eleDzGV", &eleDzGV, &b_eleDzGV);
   fChain->SetBranchAddress("eleD0Vtx", &eleD0Vtx, &b_eleD0Vtx);
   fChain->SetBranchAddress("eleDzVtx", &eleDzVtx, &b_eleDzVtx);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleHoverE12", &eleHoverE12, &b_eleHoverE12);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("elePin", &elePin, &b_elePin);
   fChain->SetBranchAddress("elePout", &elePout, &b_elePout);
   fChain->SetBranchAddress("eleTrkMomErr", &eleTrkMomErr, &b_eleTrkMomErr);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIEtaIPhi", &eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   fChain->SetBranchAddress("eleEmax", &eleEmax, &b_eleEmax);
   fChain->SetBranchAddress("eleE2ndMax", &eleE2ndMax, &b_eleE2ndMax);
   fChain->SetBranchAddress("eleETop", &eleETop, &b_eleETop);
   fChain->SetBranchAddress("eleEBottom", &eleEBottom, &b_eleEBottom);
   fChain->SetBranchAddress("eleELeft", &eleELeft, &b_eleELeft);
   fChain->SetBranchAddress("eleERight", &eleERight, &b_eleERight);
   fChain->SetBranchAddress("eleE1x5", &eleE1x5, &b_eleE1x5);
   fChain->SetBranchAddress("eleE3x3", &eleE3x3, &b_eleE3x3);
   fChain->SetBranchAddress("eleE5x5", &eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleE2x5Max", &eleE2x5Max, &b_eleE2x5Max);
   fChain->SetBranchAddress("eleE2x5Top", &eleE2x5Top, &b_eleE2x5Top);
   fChain->SetBranchAddress("eleE2x5Bottom", &eleE2x5Bottom, &b_eleE2x5Bottom);
   fChain->SetBranchAddress("eleE2x5Left", &eleE2x5Left, &b_eleE2x5Left);
   fChain->SetBranchAddress("eleE2x5Right", &eleE2x5Right, &b_eleE2x5Right);
   fChain->SetBranchAddress("eleSeedEta", &eleSeedEta, &b_eleSeedEta);
   fChain->SetBranchAddress("eleSeedE", &eleSeedE, &b_eleSeedE);
   fChain->SetBranchAddress("eleSeedPhi", &eleSeedPhi, &b_eleSeedPhi);
   fChain->SetBranchAddress("eleCrysEta", &eleCrysEta, &b_eleCrysEta);
   fChain->SetBranchAddress("eleCrysPhi", &eleCrysPhi, &b_eleCrysPhi);
   fChain->SetBranchAddress("eleCrysIEta", &eleCrysIEta, &b_eleCrysIEta);
   fChain->SetBranchAddress("eleCrysIPhi", &eleCrysIPhi, &b_eleCrysIPhi);
   fChain->SetBranchAddress("eleRegrE", &eleRegrE, &b_eleRegrE);
   fChain->SetBranchAddress("eleRegrEerr", &eleRegrEerr, &b_eleRegrEerr);
   fChain->SetBranchAddress("elePhoRegrE", &elePhoRegrE, &b_elePhoRegrE);
   fChain->SetBranchAddress("elePhoRegrEerr", &elePhoRegrEerr, &b_elePhoRegrEerr);
   fChain->SetBranchAddress("eleSeedTime", &eleSeedTime, &b_eleSeedTime);
   fChain->SetBranchAddress("eleRecoFlag", &eleRecoFlag, &b_eleRecoFlag);
   fChain->SetBranchAddress("elePos", &elePos, &b_elePos);
   //JM
   if(!isItData){
     fChain->SetBranchAddress("eleGenIndex", &eleGenIndex, &b_eleGenIndex);
     fChain->SetBranchAddress("eleGenGMomPID", &eleGenGMomPID, &b_eleGenGMomPID);
     fChain->SetBranchAddress("eleGenMomPID", &eleGenMomPID, &b_eleGenMomPID);
     fChain->SetBranchAddress("eleGenMomPt", &eleGenMomPt, &b_eleGenMomPt);
   }
   fChain->SetBranchAddress("eleIsoTrkDR03", &eleIsoTrkDR03, &b_eleIsoTrkDR03);
   fChain->SetBranchAddress("eleIsoEcalDR03", &eleIsoEcalDR03, &b_eleIsoEcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR03", &eleIsoHcalDR03, &b_eleIsoHcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR0312", &eleIsoHcalDR0312, &b_eleIsoHcalDR0312);
   fChain->SetBranchAddress("eleIsoTrkDR04", &eleIsoTrkDR04, &b_eleIsoTrkDR04);
   fChain->SetBranchAddress("eleIsoEcalDR04", &eleIsoEcalDR04, &b_eleIsoEcalDR04);
   fChain->SetBranchAddress("eleIsoHcalDR04", &eleIsoHcalDR04, &b_eleIsoHcalDR04);
   fChain->SetBranchAddress("eleIsoHcalDR0412", &eleIsoHcalDR0412, &b_eleIsoHcalDR0412);
   fChain->SetBranchAddress("eleModIsoTrk", &eleModIsoTrk, &b_eleModIsoTrk);
   fChain->SetBranchAddress("eleModIsoEcal", &eleModIsoEcal, &b_eleModIsoEcal);
   fChain->SetBranchAddress("eleModIsoHcal", &eleModIsoHcal, &b_eleModIsoHcal);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleConvDist", &eleConvDist, &b_eleConvDist);
   fChain->SetBranchAddress("eleConvDcot", &eleConvDcot, &b_eleConvDcot);
   fChain->SetBranchAddress("eleConvVtxFit", &eleConvVtxFit, &b_eleConvVtxFit);
   fChain->SetBranchAddress("eleIP3D", &eleIP3D, &b_eleIP3D);
   fChain->SetBranchAddress("eleIP3DErr", &eleIP3DErr, &b_eleIP3DErr);
   fChain->SetBranchAddress("eleIDMVANonTrig", &eleIDMVANonTrig, &b_eleIDMVANonTrig);
   fChain->SetBranchAddress("eleIDMVATrig", &eleIDMVATrig, &b_eleIDMVATrig);
   fChain->SetBranchAddress("elePFChIso03", &elePFChIso03, &b_elePFChIso03);
   fChain->SetBranchAddress("elePFPhoIso03", &elePFPhoIso03, &b_elePFPhoIso03);
   fChain->SetBranchAddress("elePFNeuIso03", &elePFNeuIso03, &b_elePFNeuIso03);
   fChain->SetBranchAddress("elePFChIso04", &elePFChIso04, &b_elePFChIso04);
   fChain->SetBranchAddress("elePFPhoIso04", &elePFPhoIso04, &b_elePFPhoIso04);
   fChain->SetBranchAddress("elePFNeuIso04", &elePFNeuIso04, &b_elePFNeuIso04);
   fChain->SetBranchAddress("eleESEffSigmaRR_x", &eleESEffSigmaRR_x, &b_eleESEffSigmaRR_x);
   fChain->SetBranchAddress("eleESEffSigmaRR_y", &eleESEffSigmaRR_y, &b_eleESEffSigmaRR_y);
   fChain->SetBranchAddress("eleESEffSigmaRR_z", &eleESEffSigmaRR_z, &b_eleESEffSigmaRR_z);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoTrg", &phoTrg, &b_phoTrg);
   fChain->SetBranchAddress("phoTrgFilter", &phoTrgFilter, &b_phoTrgFilter);
   fChain->SetBranchAddress("phoIsPhoton", &phoIsPhoton, &b_phoIsPhoton);
   fChain->SetBranchAddress("phoSCPos_x", &phoSCPos_x, &b_phoSCPos_x);
   fChain->SetBranchAddress("phoSCPos_y", &phoSCPos_y, &b_phoSCPos_y);
   fChain->SetBranchAddress("phoSCPos_z", &phoSCPos_z, &b_phoSCPos_z);
   fChain->SetBranchAddress("phoCaloPos_x", &phoCaloPos_x, &b_phoCaloPos_x);
   fChain->SetBranchAddress("phoCaloPos_y", &phoCaloPos_y, &b_phoCaloPos_y);
   fChain->SetBranchAddress("phoCaloPos_z", &phoCaloPos_z, &b_phoCaloPos_z);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoVtx_x", &phoVtx_x, &b_phoVtx_x);
   fChain->SetBranchAddress("phoVtx_y", &phoVtx_y, &b_phoVtx_y);
   fChain->SetBranchAddress("phoVtx_z", &phoVtx_z, &b_phoVtx_z);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEtVtx", &phoEtVtx, &b_phoEtVtx);
   fChain->SetBranchAddress("phoEtaVtx", &phoEtaVtx, &b_phoEtaVtx);
   fChain->SetBranchAddress("phoPhiVtx", &phoPhiVtx, &b_phoPhiVtx);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoNClus", &phoNClus, &b_phoNClus);
   fChain->SetBranchAddress("phoTrkIsoHollowDR03", &phoTrkIsoHollowDR03, &b_phoTrkIsoHollowDR03);
   fChain->SetBranchAddress("phoEcalIsoDR03", &phoEcalIsoDR03, &b_phoEcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoDR03", &phoHcalIsoDR03, &b_phoHcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoDR0312", &phoHcalIsoDR0312, &b_phoHcalIsoDR0312);
   fChain->SetBranchAddress("phoTrkIsoHollowDR04", &phoTrkIsoHollowDR04, &b_phoTrkIsoHollowDR04);
   fChain->SetBranchAddress("phoCiCdRtoTrk", &phoCiCdRtoTrk, &b_phoCiCdRtoTrk);
   fChain->SetBranchAddress("phoEcalIsoDR04", &phoEcalIsoDR04, &b_phoEcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR04", &phoHcalIsoDR04, &b_phoHcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR0412", &phoHcalIsoDR0412, &b_phoHcalIsoDR0412);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoHoverE12", &phoHoverE12, &b_phoHoverE12);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi", &phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
   fChain->SetBranchAddress("phoCiCPF4phopfIso03", &phoCiCPF4phopfIso03, &b_phoCiCPF4phopfIso03);
   fChain->SetBranchAddress("phoCiCPF4phopfIso04", &phoCiCPF4phopfIso04, &b_phoCiCPF4phopfIso04);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso02", &phoCiCPF4chgpfIso02, &b_phoCiCPF4chgpfIso02);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso03", &phoCiCPF4chgpfIso03, &b_phoCiCPF4chgpfIso03);
   fChain->SetBranchAddress("phoCiCPF4chgpfIso04", &phoCiCPF4chgpfIso04, &b_phoCiCPF4chgpfIso04);
   fChain->SetBranchAddress("phoEmax", &phoEmax, &b_phoEmax);
   fChain->SetBranchAddress("phoETop", &phoETop, &b_phoETop);
   fChain->SetBranchAddress("phoEBottom", &phoEBottom, &b_phoEBottom);
   fChain->SetBranchAddress("phoELeft", &phoELeft, &b_phoELeft);
   fChain->SetBranchAddress("phoERight", &phoERight, &b_phoERight);
   fChain->SetBranchAddress("phoE2ndMax", &phoE2ndMax, &b_phoE2ndMax);
   fChain->SetBranchAddress("phoE3x3", &phoE3x3, &b_phoE3x3);
   fChain->SetBranchAddress("phoE3x1", &phoE3x1, &b_phoE3x1);
   fChain->SetBranchAddress("phoE1x3", &phoE1x3, &b_phoE1x3);
   fChain->SetBranchAddress("phoE5x5", &phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoE1x5", &phoE1x5, &b_phoE1x5);
   fChain->SetBranchAddress("phoE2x2", &phoE2x2, &b_phoE2x2);
   fChain->SetBranchAddress("phoE2x5Max", &phoE2x5Max, &b_phoE2x5Max);
   fChain->SetBranchAddress("phoE2x5Top", &phoE2x5Top, &b_phoE2x5Top);
   fChain->SetBranchAddress("phoE2x5Bottom", &phoE2x5Bottom, &b_phoE2x5Bottom);
   fChain->SetBranchAddress("phoE2x5Left", &phoE2x5Left, &b_phoE2x5Left);
   fChain->SetBranchAddress("phoE2x5Right", &phoE2x5Right, &b_phoE2x5Right);
   fChain->SetBranchAddress("phoSeedE", &phoSeedE, &b_phoSeedE);
   fChain->SetBranchAddress("phoSeedEta", &phoSeedEta, &b_phoSeedEta);
   fChain->SetBranchAddress("phoSeedPhi", &phoSeedPhi, &b_phoSeedPhi);
   fChain->SetBranchAddress("phoCrysEta", &phoCrysEta, &b_phoCrysEta);
   fChain->SetBranchAddress("phoCrysPhi", &phoCrysPhi, &b_phoCrysPhi);
   fChain->SetBranchAddress("phoCrysIEta", &phoCrysIEta, &b_phoCrysIEta);
   fChain->SetBranchAddress("phoCrysIPhi", &phoCrysIPhi, &b_phoCrysIPhi);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoSCRChIso", &phoSCRChIso, &b_phoSCRChIso);
   fChain->SetBranchAddress("phoSCRPhoIso", &phoSCRPhoIso, &b_phoSCRPhoIso);
   fChain->SetBranchAddress("phoSCRNeuIso", &phoSCRNeuIso, &b_phoSCRNeuIso);
   fChain->SetBranchAddress("phoSCRChIso04", &phoSCRChIso04, &b_phoSCRChIso04);
   fChain->SetBranchAddress("phoSCRPhoIso04", &phoSCRPhoIso04, &b_phoSCRPhoIso04);
   fChain->SetBranchAddress("phoSCRNeuIso04", &phoSCRNeuIso04, &b_phoSCRNeuIso04);
   fChain->SetBranchAddress("phoRandConeChIso", &phoRandConeChIso, &b_phoRandConeChIso);
   fChain->SetBranchAddress("phoRandConePhoIso", &phoRandConePhoIso, &b_phoRandConePhoIso);
   fChain->SetBranchAddress("phoRandConeNeuIso", &phoRandConeNeuIso, &b_phoRandConeNeuIso);
   fChain->SetBranchAddress("phoRandConeChIso04", &phoRandConeChIso04, &b_phoRandConeChIso04);
   fChain->SetBranchAddress("phoRandConePhoIso04", &phoRandConePhoIso04, &b_phoRandConePhoIso04);
   fChain->SetBranchAddress("phoRandConeNeuIso04", &phoRandConeNeuIso04, &b_phoRandConeNeuIso04);
   fChain->SetBranchAddress("phoRegrE", &phoRegrE, &b_phoRegrE);
   fChain->SetBranchAddress("phoRegrEerr", &phoRegrEerr, &b_phoRegrEerr);
   fChain->SetBranchAddress("phoSeedTime", &phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedDetId1", &phoSeedDetId1, &b_phoSeedDetId1);
   fChain->SetBranchAddress("phoSeedDetId2", &phoSeedDetId2, &b_phoSeedDetId2);
   fChain->SetBranchAddress("phoLICTD", &phoLICTD, &b_phoLICTD);
   fChain->SetBranchAddress("phoRecoFlag", &phoRecoFlag, &b_phoRecoFlag);
   fChain->SetBranchAddress("phoPos", &phoPos, &b_phoPos);
   //JM
   if(!isItData){
     fChain->SetBranchAddress("phoGenIndex", &phoGenIndex, &b_phoGenIndex);
     fChain->SetBranchAddress("phoGenGMomPID", &phoGenGMomPID, &b_phoGenGMomPID);
     fChain->SetBranchAddress("phoGenMomPID", &phoGenMomPID, &b_phoGenMomPID);
     fChain->SetBranchAddress("phoGenMomPt", &phoGenMomPt, &b_phoGenMomPt);
   }
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEn", &phoESEn, &b_phoESEn);
   fChain->SetBranchAddress("phoSCEt", &phoSCEt, &b_phoSCEt);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phoOverlap", &phoOverlap, &b_phoOverlap);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("pho_hasConvPf", &pho_hasConvPf, &b_pho_hasConvPf);
   fChain->SetBranchAddress("pho_hasSLConvPf", &pho_hasSLConvPf, &b_pho_hasSLConvPf);
   fChain->SetBranchAddress("pho_pfconvVtxZ", &pho_pfconvVtxZ, &b_pho_pfconvVtxZ);
   fChain->SetBranchAddress("pho_pfconvVtxZErr", &pho_pfconvVtxZErr, &b_pho_pfconvVtxZErr);
   fChain->SetBranchAddress("pho_nSLConv", &pho_nSLConv, &b_pho_nSLConv);
   fChain->SetBranchAddress("pho_pfSLConvPos_x", &pho_pfSLConvPos_x, &b_pho_pfSLConvPos_x);
   fChain->SetBranchAddress("pho_pfSLConvPos_y", &pho_pfSLConvPos_y, &b_pho_pfSLConvPos_y);
   fChain->SetBranchAddress("pho_pfSLConvPos_z", &pho_pfSLConvPos_z, &b_pho_pfSLConvPos_z);
   fChain->SetBranchAddress("pho_pfSLConvVtxZ", &pho_pfSLConvVtxZ, &b_pho_pfSLConvVtxZ);
   fChain->SetBranchAddress("phoIsConv", &phoIsConv, &b_phoIsConv);
   fChain->SetBranchAddress("phoNConv", &phoNConv, &b_phoNConv);
   fChain->SetBranchAddress("phoConvInvMass", &phoConvInvMass, &b_phoConvInvMass);
   fChain->SetBranchAddress("phoConvCotTheta", &phoConvCotTheta, &b_phoConvCotTheta);
   fChain->SetBranchAddress("phoConvEoverP", &phoConvEoverP, &b_phoConvEoverP);
   fChain->SetBranchAddress("phoConvZofPVfromTrks", &phoConvZofPVfromTrks, &b_phoConvZofPVfromTrks);
   fChain->SetBranchAddress("phoConvMinDist", &phoConvMinDist, &b_phoConvMinDist);
   fChain->SetBranchAddress("phoConvdPhiAtVtx", &phoConvdPhiAtVtx, &b_phoConvdPhiAtVtx);
   fChain->SetBranchAddress("phoConvdPhiAtCalo", &phoConvdPhiAtCalo, &b_phoConvdPhiAtCalo);
   fChain->SetBranchAddress("phoConvdEtaAtCalo", &phoConvdEtaAtCalo, &b_phoConvdEtaAtCalo);
   fChain->SetBranchAddress("phoConvTrkd0_x", &phoConvTrkd0_x, &b_phoConvTrkd0_x);
   fChain->SetBranchAddress("phoConvTrkd0_y", &phoConvTrkd0_y, &b_phoConvTrkd0_y);
   fChain->SetBranchAddress("phoConvTrkPin_x", &phoConvTrkPin_x, &b_phoConvTrkPin_x);
   fChain->SetBranchAddress("phoConvTrkPin_y", &phoConvTrkPin_y, &b_phoConvTrkPin_y);
   fChain->SetBranchAddress("phoConvTrkPout_x", &phoConvTrkPout_x, &b_phoConvTrkPout_x);
   fChain->SetBranchAddress("phoConvTrkPout_y", &phoConvTrkPout_y, &b_phoConvTrkPout_y);
   fChain->SetBranchAddress("phoConvTrkdz_x", &phoConvTrkdz_x, &b_phoConvTrkdz_x);
   fChain->SetBranchAddress("phoConvTrkdz_y", &phoConvTrkdz_y, &b_phoConvTrkdz_y);
   fChain->SetBranchAddress("phoConvTrkdzErr_x", &phoConvTrkdzErr_x, &b_phoConvTrkdzErr_x);
   fChain->SetBranchAddress("phoConvTrkdzErr_y", &phoConvTrkdzErr_y, &b_phoConvTrkdzErr_y);
   fChain->SetBranchAddress("phoConvChi2", &phoConvChi2, &b_phoConvChi2);
   fChain->SetBranchAddress("phoConvChi2Prob", &phoConvChi2Prob, &b_phoConvChi2Prob);
   fChain->SetBranchAddress("phoConvNTrks", &phoConvNTrks, &b_phoConvNTrks);
   fChain->SetBranchAddress("phoConvCharge1", &phoConvCharge1, &b_phoConvCharge1);
   fChain->SetBranchAddress("phoConvCharge2", &phoConvCharge2, &b_phoConvCharge2);
   fChain->SetBranchAddress("phoConvValidVtx", &phoConvValidVtx, &b_phoConvValidVtx);
   fChain->SetBranchAddress("phoConvLikeLihood", &phoConvLikeLihood, &b_phoConvLikeLihood);
   fChain->SetBranchAddress("phoConvP4_0", &phoConvP4_0, &b_phoConvP4_0);
   fChain->SetBranchAddress("phoConvP4_1", &phoConvP4_1, &b_phoConvP4_1);
   fChain->SetBranchAddress("phoConvP4_2", &phoConvP4_2, &b_phoConvP4_2);
   fChain->SetBranchAddress("phoConvP4_3", &phoConvP4_3, &b_phoConvP4_3);
   fChain->SetBranchAddress("phoConvVtx_x", &phoConvVtx_x, &b_phoConvVtx_x);
   fChain->SetBranchAddress("phoConvVtx_y", &phoConvVtx_y, &b_phoConvVtx_y);
   fChain->SetBranchAddress("phoConvVtx_z", &phoConvVtx_z, &b_phoConvVtx_z);
   fChain->SetBranchAddress("phoConvVtxErr_x", &phoConvVtxErr_x, &b_phoConvVtxErr_x);
   fChain->SetBranchAddress("phoConvVtxErr_y", &phoConvVtxErr_y, &b_phoConvVtxErr_y);
   fChain->SetBranchAddress("phoConvVtxErr_z", &phoConvVtxErr_z, &b_phoConvVtxErr_z);
   fChain->SetBranchAddress("phoConvPairMomentum_x", &phoConvPairMomentum_x, &b_phoConvPairMomentum_x);
   fChain->SetBranchAddress("phoConvPairMomentum_y", &phoConvPairMomentum_y, &b_phoConvPairMomentum_y);
   fChain->SetBranchAddress("phoConvPairMomentum_z", &phoConvPairMomentum_z, &b_phoConvPairMomentum_z);
   fChain->SetBranchAddress("phoConvRefittedMomentum_x", &phoConvRefittedMomentum_x, &b_phoConvRefittedMomentum_x);
   fChain->SetBranchAddress("phoConvRefittedMomentum_y", &phoConvRefittedMomentum_y, &b_phoConvRefittedMomentum_y);
   fChain->SetBranchAddress("phoConvRefittedMomentum_z", &phoConvRefittedMomentum_z, &b_phoConvRefittedMomentum_z);
   fChain->SetBranchAddress("SingleLegConv", &SingleLegConv, &b_SingleLegConv);
   fChain->SetBranchAddress("phoPFConvVtx_x", &phoPFConvVtx_x, &b_phoPFConvVtx_x);
   fChain->SetBranchAddress("phoPFConvVtx_y", &phoPFConvVtx_y, &b_phoPFConvVtx_y);
   fChain->SetBranchAddress("phoPFConvVtx_z", &phoPFConvVtx_z, &b_phoPFConvVtx_z);
   fChain->SetBranchAddress("phoPFConvMom_x", &phoPFConvMom_x, &b_phoPFConvMom_x);
   fChain->SetBranchAddress("phoPFConvMom_y", &phoPFConvMom_y, &b_phoPFConvMom_y);
   fChain->SetBranchAddress("phoPFConvMom_z", &phoPFConvMom_z, &b_phoPFConvMom_z);
   fChain->SetBranchAddress("phoESEffSigmaRR_x", &phoESEffSigmaRR_x, &b_phoESEffSigmaRR_x);
   fChain->SetBranchAddress("phoESEffSigmaRR_y", &phoESEffSigmaRR_y, &b_phoESEffSigmaRR_y);
   fChain->SetBranchAddress("phoESEffSigmaRR_z", &phoESEffSigmaRR_z, &b_phoESEffSigmaRR_z);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muTrg", &muTrg, &b_muTrg);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muPz", &muPz, &b_muPz);
   fChain->SetBranchAddress("muVtx_x", &muVtx_x, &b_muVtx_x);
   fChain->SetBranchAddress("muVtx_y", &muVtx_y, &b_muVtx_y);
   fChain->SetBranchAddress("muVtx_z", &muVtx_z, &b_muVtx_z);
   fChain->SetBranchAddress("muVtxGlb_x", &muVtxGlb_x, &b_muVtxGlb_x);
   fChain->SetBranchAddress("muVtxGlb_y", &muVtxGlb_y, &b_muVtxGlb_y);
   fChain->SetBranchAddress("muVtxGlb_z", &muVtxGlb_z, &b_muVtxGlb_z);
   //JM
   if(!isItData)
     fChain->SetBranchAddress("muGenIndex", &muGenIndex, &b_muGenIndex);
   fChain->SetBranchAddress("mucktPt", &mucktPt, &b_mucktPt);
   fChain->SetBranchAddress("mucktPtErr", &mucktPtErr, &b_mucktPtErr);
   fChain->SetBranchAddress("mucktEta", &mucktEta, &b_mucktEta);
   fChain->SetBranchAddress("mucktPhi", &mucktPhi, &b_mucktPhi);
   fChain->SetBranchAddress("mucktdxy", &mucktdxy, &b_mucktdxy);
   fChain->SetBranchAddress("mucktdz", &mucktdz, &b_mucktdz);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muIsoCalo", &muIsoCalo, &b_muIsoCalo);
   fChain->SetBranchAddress("muIsoEcal", &muIsoEcal, &b_muIsoEcal);
   fChain->SetBranchAddress("muIsoHcal", &muIsoHcal, &b_muIsoHcal);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerChi2NDF", &muInnerChi2NDF, &b_muInnerChi2NDF);
   fChain->SetBranchAddress("muPFIsoR04_CH", &muPFIsoR04_CH, &b_muPFIsoR04_CH);
   fChain->SetBranchAddress("muPFIsoR04_NH", &muPFIsoR04_NH, &b_muPFIsoR04_NH);
   fChain->SetBranchAddress("muPFIsoR04_Pho", &muPFIsoR04_Pho, &b_muPFIsoR04_Pho);
   fChain->SetBranchAddress("muPFIsoR04_PU", &muPFIsoR04_PU, &b_muPFIsoR04_PU);
   fChain->SetBranchAddress("muPFIsoR04_CPart", &muPFIsoR04_CPart, &b_muPFIsoR04_CPart);
   fChain->SetBranchAddress("muPFIsoR04_NHHT", &muPFIsoR04_NHHT, &b_muPFIsoR04_NHHT);
   fChain->SetBranchAddress("muPFIsoR04_PhoHT", &muPFIsoR04_PhoHT, &b_muPFIsoR04_PhoHT);
   fChain->SetBranchAddress("muPFIsoR03_CH", &muPFIsoR03_CH, &b_muPFIsoR03_CH);
   fChain->SetBranchAddress("muPFIsoR03_NH", &muPFIsoR03_NH, &b_muPFIsoR03_NH);
   fChain->SetBranchAddress("muPFIsoR03_Pho", &muPFIsoR03_Pho, &b_muPFIsoR03_Pho);
   fChain->SetBranchAddress("muPFIsoR03_PU", &muPFIsoR03_PU, &b_muPFIsoR03_PU);
   fChain->SetBranchAddress("muPFIsoR03_CPart", &muPFIsoR03_CPart, &b_muPFIsoR03_CPart);
   fChain->SetBranchAddress("muPFIsoR03_NHHT", &muPFIsoR03_NHHT, &b_muPFIsoR03_NHHT);
   fChain->SetBranchAddress("muPFIsoR03_PhoHT", &muPFIsoR03_PhoHT, &b_muPFIsoR03_PhoHT);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muD0GV", &muD0GV, &b_muD0GV);
   fChain->SetBranchAddress("muDzGV", &muDzGV, &b_muDzGV);
   fChain->SetBranchAddress("muD0Vtx", &muD0Vtx, &b_muD0Vtx);
   fChain->SetBranchAddress("muDzVtx", &muDzVtx, &b_muDzVtx);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muInnerD0GV", &muInnerD0GV, &b_muInnerD0GV);
   fChain->SetBranchAddress("muInnerDzGV", &muInnerDzGV, &b_muInnerDzGV);
   fChain->SetBranchAddress("muInnerPt", &muInnerPt, &b_muInnerPt);
   fChain->SetBranchAddress("muInnerPtErr", &muInnerPtErr, &b_muInnerPtErr);
   fChain->SetBranchAddress("muNumberOfValidTrkLayers", &muNumberOfValidTrkLayers, &b_muNumberOfValidTrkLayers);
   fChain->SetBranchAddress("muNumberOfValidTrkHits", &muNumberOfValidTrkHits, &b_muNumberOfValidTrkHits);
   fChain->SetBranchAddress("muNumberOfValidPixelLayers", &muNumberOfValidPixelLayers, &b_muNumberOfValidPixelLayers);
   fChain->SetBranchAddress("muNumberOfValidPixelHits", &muNumberOfValidPixelHits, &b_muNumberOfValidPixelHits);
   fChain->SetBranchAddress("muNumberOfValidMuonHits", &muNumberOfValidMuonHits, &b_muNumberOfValidMuonHits);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muChambers", &muChambers, &b_muChambers);
   fChain->SetBranchAddress("muIP3D", &muIP3D, &b_muIP3D);
   fChain->SetBranchAddress("muIP3DErr", &muIP3DErr, &b_muIP3DErr);
   /*
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("tauDecayModeFinding", &tauDecayModeFinding, &b_tauDecayModeFinding);
   fChain->SetBranchAddress("tauAgainstElectronLooseMVA3", &tauAgainstElectronLooseMVA3, &b_tauAgainstElectronLooseMVA3);
   fChain->SetBranchAddress("tauAgainstElectronMediumMVA3", &tauAgainstElectronMediumMVA3, &b_tauAgainstElectronMediumMVA3);
   fChain->SetBranchAddress("tauAgainstElectronTightMVA3", &tauAgainstElectronTightMVA3, &b_tauAgainstElectronTightMVA3);
   fChain->SetBranchAddress("tauAgainstElectronVTightMVA3", &tauAgainstElectronVTightMVA3, &b_tauAgainstElectronVTightMVA3);
   fChain->SetBranchAddress("tauAgainstElectronDeadECAL", &tauAgainstElectronDeadECAL, &b_tauAgainstElectronDeadECAL);
   fChain->SetBranchAddress("tauAgainstMuonLoose2", &tauAgainstMuonLoose2, &b_tauAgainstMuonLoose2);
   fChain->SetBranchAddress("tauAgainstMuonMedium2", &tauAgainstMuonMedium2, &b_tauAgainstMuonMedium2);
   fChain->SetBranchAddress("tauAgainstMuonTight2", &tauAgainstMuonTight2, &b_tauAgainstMuonTight2);
   fChain->SetBranchAddress("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits, &b_tauCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("tauLooseCombinedIsolationDeltaBetaCorr3Hits", &tauLooseCombinedIsolationDeltaBetaCorr3Hits, &b_tauLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauMediumCombinedIsolationDeltaBetaCorr3Hits", &tauMediumCombinedIsolationDeltaBetaCorr3Hits, &b_tauMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauTightCombinedIsolationDeltaBetaCorr3Hits", &tauTightCombinedIsolationDeltaBetaCorr3Hits, &b_tauTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauEta", &tauEta, &b_tauEta);
   fChain->SetBranchAddress("tauPhi", &tauPhi, &b_tauPhi);
   fChain->SetBranchAddress("tauPt", &tauPt, &b_tauPt);
   fChain->SetBranchAddress("tauEt", &tauEt, &b_tauEt);
   fChain->SetBranchAddress("tauCharge", &tauCharge, &b_tauCharge);
   fChain->SetBranchAddress("tauDecayMode", &tauDecayMode, &b_tauDecayMode);
   fChain->SetBranchAddress("tauEMFraction", &tauEMFraction, &b_tauEMFraction);
   fChain->SetBranchAddress("tauHCAL3x3OverPLead", &tauHCAL3x3OverPLead, &b_tauHCAL3x3OverPLead);
   fChain->SetBranchAddress("tauHCALMaxOverPLead", &tauHCALMaxOverPLead, &b_tauHCALMaxOverPLead);
   fChain->SetBranchAddress("tauHCALTotOverPLead", &tauHCALTotOverPLead, &b_tauHCALTotOverPLead);
   fChain->SetBranchAddress("tauIsolationPFChargedHadrCandsPtSum", &tauIsolationPFChargedHadrCandsPtSum, &b_tauIsolationPFChargedHadrCandsPtSum);
   fChain->SetBranchAddress("tauIsolationPFGammaCandsEtSum", &tauIsolationPFGammaCandsEtSum, &b_tauIsolationPFGammaCandsEtSum);
   fChain->SetBranchAddress("tauLeadPFChargedHadrCandsignedSipt", &tauLeadPFChargedHadrCandsignedSipt, &b_tauLeadPFChargedHadrCandsignedSipt);
   fChain->SetBranchAddress("tauLeadChargedHadronExists", &tauLeadChargedHadronExists, &b_tauLeadChargedHadronExists);
   fChain->SetBranchAddress("tauLeadChargedHadronEta", &tauLeadChargedHadronEta, &b_tauLeadChargedHadronEta);
   fChain->SetBranchAddress("tauLeadChargedHadronPhi", &tauLeadChargedHadronPhi, &b_tauLeadChargedHadronPhi);
   fChain->SetBranchAddress("tauLeadChargedHadronPt", &tauLeadChargedHadronPt, &b_tauLeadChargedHadronPt);
   */
   fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
   fChain->SetBranchAddress("rho25_neu", &rho25_neu, &b_rho25_neu);
   fChain->SetBranchAddress("rho25_muPFiso", &rho25_muPFiso, &b_rho25_muPFiso);
   fChain->SetBranchAddress("rho25_elePFiso", &rho25_elePFiso, &b_rho25_elePFiso);
   fChain->SetBranchAddress("rho2011", &rho2011, &b_rho2011);
   fChain->SetBranchAddress("rho2012", &rho2012, &b_rho2012);
   fChain->SetBranchAddress("QGTag_MLP", &QGTag_MLP, &b_QGTag_MLP);
   fChain->SetBranchAddress("QGTag_likelihood", &QGTag_likelihood, &b_QGTag_likelihood);
   fChain->SetBranchAddress("nCA8Jet", &nCA8Jet, &b_nCA8Jet);
   fChain->SetBranchAddress("CA8JetPt", &CA8JetPt, &b_CA8JetPt);
   fChain->SetBranchAddress("CA8JetEta", &CA8JetEta, &b_CA8JetEta);
   fChain->SetBranchAddress("CA8JetPhi", &CA8JetPhi, &b_CA8JetPhi);
   fChain->SetBranchAddress("CA8JetMass", &CA8JetMass, &b_CA8JetMass);
   fChain->SetBranchAddress("CA8JetArea", &CA8JetArea, &b_CA8JetArea);
   fChain->SetBranchAddress("CA8Jet_tau1", &CA8Jet_tau1, &b_CA8Jet_tau1);
   fChain->SetBranchAddress("CA8Jet_tau2", &CA8Jet_tau2, &b_CA8Jet_tau2);
   fChain->SetBranchAddress("CA8Jet_tau3", &CA8Jet_tau3, &b_CA8Jet_tau3);
   fChain->SetBranchAddress("CA8JetCHF", &CA8JetCHF, &b_CA8JetCHF);
   fChain->SetBranchAddress("CA8JetNHF", &CA8JetNHF, &b_CA8JetNHF);
   fChain->SetBranchAddress("CA8JetCEF", &CA8JetCEF, &b_CA8JetCEF);
   fChain->SetBranchAddress("CA8JetNEF", &CA8JetNEF, &b_CA8JetNEF);
   fChain->SetBranchAddress("CA8JetNCH", &CA8JetNCH, &b_CA8JetNCH);
   fChain->SetBranchAddress("CA8Jetnconstituents", &CA8Jetnconstituents, &b_CA8Jetnconstituents);
   fChain->SetBranchAddress("CA8prunedJetMass", &CA8prunedJetMass, &b_CA8prunedJetMass);
   fChain->SetBranchAddress("CA8prunedJet_nSubJets", &CA8prunedJet_nSubJets, &b_CA8prunedJet_nSubJets);
   fChain->SetBranchAddress("CA8prunedJet_SubjetPt", &CA8prunedJet_SubjetPt, &b_CA8prunedJet_SubjetPt);
   fChain->SetBranchAddress("CA8prunedJet_SubjetEta", &CA8prunedJet_SubjetEta, &b_CA8prunedJet_SubjetEta);
   fChain->SetBranchAddress("CA8prunedJet_SubjetPhi", &CA8prunedJet_SubjetPhi, &b_CA8prunedJet_SubjetPhi);
   fChain->SetBranchAddress("CA8prunedJet_SubjetMass", &CA8prunedJet_SubjetMass, &b_CA8prunedJet_SubjetMass);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetTrg", &jetTrg, &b_jetTrg);
   fChain->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetCharge", &jetCharge, &b_jetCharge);
   fChain->SetBranchAddress("jetEt", &jetEt, &b_jetEt);
   fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", &jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", &jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", &jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", &jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetHFHAE", &jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", &jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConstituents", &jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("jetCombinedSecondaryVtxBJetTags", &jetCombinedSecondaryVtxBJetTags, &b_jetCombinedSecondaryVtxBJetTags);
   fChain->SetBranchAddress("jetCombinedSecondaryVtxMVABJetTags", &jetCombinedSecondaryVtxMVABJetTags, &b_jetCombinedSecondaryVtxMVABJetTags);
   fChain->SetBranchAddress("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags, &b_jetJetProbabilityBJetTags);
   fChain->SetBranchAddress("jetJetBProbabilityBJetTags", &jetJetBProbabilityBJetTags, &b_jetJetBProbabilityBJetTags);
   fChain->SetBranchAddress("jetBetaStar", &jetBetaStar, &b_jetBetaStar);
   fChain->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetDRMean", &jetDRMean, &b_jetDRMean);
   fChain->SetBranchAddress("jetDR2Mean", &jetDR2Mean, &b_jetDR2Mean);
   fChain->SetBranchAddress("jetDZ", &jetDZ, &b_jetDZ);
   fChain->SetBranchAddress("jetFrac01", &jetFrac01, &b_jetFrac01);
   fChain->SetBranchAddress("jetFrac02", &jetFrac02, &b_jetFrac02);
   fChain->SetBranchAddress("jetFrac03", &jetFrac03, &b_jetFrac03);
   fChain->SetBranchAddress("jetFrac04", &jetFrac04, &b_jetFrac04);
   fChain->SetBranchAddress("jetFrac05", &jetFrac05, &b_jetFrac05);
   fChain->SetBranchAddress("jetFrac06", &jetFrac06, &b_jetFrac06);
   fChain->SetBranchAddress("jetFrac07", &jetFrac07, &b_jetFrac07);
   fChain->SetBranchAddress("jetBeta", &jetBeta, &b_jetBeta);
   fChain->SetBranchAddress("jetBetaStarCMG", &jetBetaStarCMG, &b_jetBetaStarCMG);
   fChain->SetBranchAddress("jetBetaStarClassic", &jetBetaStarClassic, &b_jetBetaStarClassic);
   fChain->SetBranchAddress("jetBetaExt", &jetBetaExt, &b_jetBetaExt);
   fChain->SetBranchAddress("jetBetaStarCMGExt", &jetBetaStarCMGExt, &b_jetBetaStarCMGExt);
   fChain->SetBranchAddress("jetBetaStarClassicExt", &jetBetaStarClassicExt, &b_jetBetaStarClassicExt);
   fChain->SetBranchAddress("jetNNeutrals", &jetNNeutrals, &b_jetNNeutrals);
   fChain->SetBranchAddress("jetNCharged", &jetNCharged, &b_jetNCharged);
   fChain->SetBranchAddress("jetMVAs", &jetMVAs, &b_jetMVAs);
   fChain->SetBranchAddress("jetWPLevels", &jetWPLevels, &b_jetWPLevels);
   fChain->SetBranchAddress("jetMVAsExt_simple", &jetMVAsExt_simple, &b_jetMVAsExt_simple);
   fChain->SetBranchAddress("jetWPLevelsExt_simple", &jetWPLevelsExt_simple, &b_jetWPLevelsExt_simple);
   fChain->SetBranchAddress("jetMVAsExt_full", &jetMVAsExt_full, &b_jetMVAsExt_full);
   fChain->SetBranchAddress("jetWPLevelsExt_full", &jetWPLevelsExt_full, &b_jetWPLevelsExt_full);
   fChain->SetBranchAddress("jetMVAsExt_cutBased", &jetMVAsExt_cutBased, &b_jetMVAsExt_cutBased);
   fChain->SetBranchAddress("jetWPLevelsExt_cutBased", &jetWPLevelsExt_cutBased, &b_jetWPLevelsExt_cutBased);
   fChain->SetBranchAddress("jetMVAsExt_philv1", &jetMVAsExt_philv1, &b_jetMVAsExt_philv1);
   fChain->SetBranchAddress("jetWPLevelsExt_philv1", &jetWPLevelsExt_philv1, &b_jetWPLevelsExt_philv1);
   fChain->SetBranchAddress("jetMt", &jetMt, &b_jetMt);
   fChain->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
   fChain->SetBranchAddress("jetLeadTrackPt", &jetLeadTrackPt, &b_jetLeadTrackPt);
   fChain->SetBranchAddress("jetVtxPt", &jetVtxPt, &b_jetVtxPt);
   fChain->SetBranchAddress("jetVtxMass", &jetVtxMass, &b_jetVtxMass);
   fChain->SetBranchAddress("jetVtx3dL", &jetVtx3dL, &b_jetVtx3dL);
   fChain->SetBranchAddress("jetVtx3deL", &jetVtx3deL, &b_jetVtx3deL);
   fChain->SetBranchAddress("jetSoftLeptPt", &jetSoftLeptPt, &b_jetSoftLeptPt);
   fChain->SetBranchAddress("jetSoftLeptPtRel", &jetSoftLeptPtRel, &b_jetSoftLeptPtRel);
   fChain->SetBranchAddress("jetSoftLeptdR", &jetSoftLeptdR, &b_jetSoftLeptdR);
   fChain->SetBranchAddress("jetSoftLeptIdlooseMu", &jetSoftLeptIdlooseMu, &b_jetSoftLeptIdlooseMu);
   fChain->SetBranchAddress("jetSoftLeptIdEle95", &jetSoftLeptIdEle95, &b_jetSoftLeptIdEle95);
   fChain->SetBranchAddress("jetDPhiMETJet", &jetDPhiMETJet, &b_jetDPhiMETJet);
   fChain->SetBranchAddress("jetPuJetIdL", &jetPuJetIdL, &b_jetPuJetIdL);
   fChain->SetBranchAddress("jetPuJetIdM", &jetPuJetIdM, &b_jetPuJetIdM);
   fChain->SetBranchAddress("jetPuJetIdT", &jetPuJetIdT, &b_jetPuJetIdT);
   //JM
   if(!isItData){
     fChain->SetBranchAddress("jetPartonID", &jetPartonID, &b_jetPartonID);
     fChain->SetBranchAddress("jetGenJetIndex", &jetGenJetIndex, &b_jetGenJetIndex);
     fChain->SetBranchAddress("jetGenJetEn", &jetGenJetEn, &b_jetGenJetEn);
     fChain->SetBranchAddress("jetGenJetPt", &jetGenJetPt, &b_jetGenJetPt);
     fChain->SetBranchAddress("jetGenJetEta", &jetGenJetEta, &b_jetGenJetEta);
     fChain->SetBranchAddress("jetGenJetPhi", &jetGenJetPhi, &b_jetGenJetPhi);
     fChain->SetBranchAddress("jetGenPartonID", &jetGenPartonID, &b_jetGenPartonID);
     fChain->SetBranchAddress("jetGenEn", &jetGenEn, &b_jetGenEn);
     fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
     fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
     fChain->SetBranchAddress("jetGenPhi", &jetGenPhi, &b_jetGenPhi);
     fChain->SetBranchAddress("jetGenPartonMomID", &jetGenPartonMomID, &b_jetGenPartonMomID);
   }
   fChain->SetBranchAddress("nLowPtJet", &nLowPtJet, &b_nLowPtJet);
   fChain->SetBranchAddress("jetLowPtEn", &jetLowPtEn, &b_jetLowPtEn);
   fChain->SetBranchAddress("jetLowPtPt", &jetLowPtPt, &b_jetLowPtPt);
   fChain->SetBranchAddress("jetLowPtEta", &jetLowPtEta, &b_jetLowPtEta);
   fChain->SetBranchAddress("jetLowPtPhi", &jetLowPtPhi, &b_jetLowPtPhi);
   fChain->SetBranchAddress("jetLowPtCharge", &jetLowPtCharge, &b_jetLowPtCharge);
   fChain->SetBranchAddress("jetLowPtEt", &jetLowPtEt, &b_jetLowPtEt);
   fChain->SetBranchAddress("jetLowPtRawPt", &jetLowPtRawPt, &b_jetLowPtRawPt);
   fChain->SetBranchAddress("jetLowPtRawEn", &jetLowPtRawEn, &b_jetLowPtRawEn);
   fChain->SetBranchAddress("jetLowPtArea", &jetLowPtArea, &b_jetLowPtArea);
   //JM
   if(!isItData){
     fChain->SetBranchAddress("jetLowPtPartonID", &jetLowPtPartonID, &b_jetLowPtPartonID);
     fChain->SetBranchAddress("jetLowPtGenJetEn", &jetLowPtGenJetEn, &b_jetLowPtGenJetEn);
     fChain->SetBranchAddress("jetLowPtGenJetPt", &jetLowPtGenJetPt, &b_jetLowPtGenJetPt);
     fChain->SetBranchAddress("jetLowPtGenJetEta", &jetLowPtGenJetEta, &b_jetLowPtGenJetEta);
     fChain->SetBranchAddress("jetLowPtGenJetPhi", &jetLowPtGenJetPhi, &b_jetLowPtGenJetPhi);
     fChain->SetBranchAddress("jetLowPtGenPartonID", &jetLowPtGenPartonID, &b_jetLowPtGenPartonID);
     fChain->SetBranchAddress("jetLowPtGenEn", &jetLowPtGenEn, &b_jetLowPtGenEn);
     fChain->SetBranchAddress("jetLowPtGenPt", &jetLowPtGenPt, &b_jetLowPtGenPt);
     fChain->SetBranchAddress("jetLowPtGenEta", &jetLowPtGenEta, &b_jetLowPtGenEta);
     fChain->SetBranchAddress("jetLowPtGenPhi", &jetLowPtGenPhi, &b_jetLowPtGenPhi);
   }
   fChain->SetBranchAddress("nConv", &nConv, &b_nConv);
   fChain->SetBranchAddress("convP4_x", &convP4_x, &b_convP4_x);
   fChain->SetBranchAddress("convP4_y", &convP4_y, &b_convP4_y);
   fChain->SetBranchAddress("convP4_z", &convP4_z, &b_convP4_z);
   fChain->SetBranchAddress("convP4_E", &convP4_E, &b_convP4_E);
   fChain->SetBranchAddress("convVtx_x", &convVtx_x, &b_convVtx_x);
   fChain->SetBranchAddress("convVtx_y", &convVtx_y, &b_convVtx_y);
   fChain->SetBranchAddress("convVtx_z", &convVtx_z, &b_convVtx_z);
   fChain->SetBranchAddress("convVtxErr_x", &convVtxErr_x, &b_convVtxErr_x);
   fChain->SetBranchAddress("convVtxErr_y", &convVtxErr_y, &b_convVtxErr_y);
   fChain->SetBranchAddress("convVtxErr_z", &convVtxErr_z, &b_convVtxErr_z);
   fChain->SetBranchAddress("convPairMomentum_x", &convPairMomentum_x, &b_convPairMomentum_x);
   fChain->SetBranchAddress("convPairMomentum_y", &convPairMomentum_y, &b_convPairMomentum_y);
   fChain->SetBranchAddress("convPairMomentum_z", &convPairMomentum_z, &b_convPairMomentum_z);
   fChain->SetBranchAddress("convRefittedMomentum_x", &convRefittedMomentum_x, &b_convRefittedMomentum_x);
   fChain->SetBranchAddress("convRefittedMomentum_y", &convRefittedMomentum_y, &b_convRefittedMomentum_y);
   fChain->SetBranchAddress("convRefittedMomentum_z", &convRefittedMomentum_z, &b_convRefittedMomentum_z);
   fChain->SetBranchAddress("convNTracks", &convNTracks, &b_convNTracks);
   fChain->SetBranchAddress("convPairInvMass", &convPairInvMass, &b_convPairInvMass);
   fChain->SetBranchAddress("convPairCotThetaSep", &convPairCotThetaSep, &b_convPairCotThetaSep);
   fChain->SetBranchAddress("convEoverP", &convEoverP, &b_convEoverP);
   fChain->SetBranchAddress("convDistOfMinApproach", &convDistOfMinApproach, &b_convDistOfMinApproach);
   fChain->SetBranchAddress("convDPhiTrksAtVtx", &convDPhiTrksAtVtx, &b_convDPhiTrksAtVtx);
   fChain->SetBranchAddress("convDPhiTrksAtEcal", &convDPhiTrksAtEcal, &b_convDPhiTrksAtEcal);
   fChain->SetBranchAddress("convDEtaTrksAtEcal", &convDEtaTrksAtEcal, &b_convDEtaTrksAtEcal);
   fChain->SetBranchAddress("convDxy", &convDxy, &b_convDxy);
   fChain->SetBranchAddress("convDz", &convDz, &b_convDz);
   fChain->SetBranchAddress("convLxy", &convLxy, &b_convLxy);
   fChain->SetBranchAddress("convLz", &convLz, &b_convLz);
   fChain->SetBranchAddress("convZofPrimVtxFromTrks", &convZofPrimVtxFromTrks, &b_convZofPrimVtxFromTrks);
   fChain->SetBranchAddress("convNHitsBeforeVtx_0", &convNHitsBeforeVtx_0, &b_convNHitsBeforeVtx_0);
   fChain->SetBranchAddress("convNHitsBeforeVtx_1", &convNHitsBeforeVtx_1, &b_convNHitsBeforeVtx_1);
   fChain->SetBranchAddress("convNSharedHits", &convNSharedHits, &b_convNSharedHits);
   fChain->SetBranchAddress("convValidVtx", &convValidVtx, &b_convValidVtx);
   fChain->SetBranchAddress("convMVALikelihood", &convMVALikelihood, &b_convMVALikelihood);
   fChain->SetBranchAddress("convChi2", &convChi2, &b_convChi2);
   fChain->SetBranchAddress("convChi2Probability", &convChi2Probability, &b_convChi2Probability);
   fChain->SetBranchAddress("convTk1Dz", &convTk1Dz, &b_convTk1Dz);
   fChain->SetBranchAddress("convTk2Dz", &convTk2Dz, &b_convTk2Dz);
   fChain->SetBranchAddress("convTk1DzErr", &convTk1DzErr, &b_convTk1DzErr);
   fChain->SetBranchAddress("convTk2DzErr", &convTk2DzErr, &b_convTk2DzErr);
   fChain->SetBranchAddress("convCh1Ch2", &convCh1Ch2, &b_convCh1Ch2);
   fChain->SetBranchAddress("convTk1D0", &convTk1D0, &b_convTk1D0);
   fChain->SetBranchAddress("convTk1Pout", &convTk1Pout, &b_convTk1Pout);
   fChain->SetBranchAddress("convTk1Pin", &convTk1Pin, &b_convTk1Pin);
   fChain->SetBranchAddress("convTk2D0", &convTk2D0, &b_convTk2D0);
   fChain->SetBranchAddress("convTk2Pout", &convTk2Pout, &b_convTk2Pout);
   fChain->SetBranchAddress("convTk2Pin", &convTk2Pin, &b_convTk2Pin);

   /// variables only in Fabrice ggNtuples for now.
   TBranch *brHard = 0;
   brHard = (TBranch*) fChain->GetListOfBranches()->FindObject("mcHardPt");

   if( brHard ) {

     fChain->SetBranchAddress("mcHardPt", &mcHardPt, &b_mcHardPt);
     fChain->SetBranchAddress("mcHardEta", &mcHardEta, &b_mcHardEta);
     fChain->SetBranchAddress("mcHardPhi", &mcHardPhi, &b_mcHardPhi);
     fChain->SetBranchAddress("mcHardM", &mcHardM, &b_mcHardM);
     fChain->SetBranchAddress("mcHardPID", &mcHardPID, &b_mcHardPID);
     fChain->SetBranchAddress("mcHardFun", &mcHardFun, &b_mcHardFun);

     fChain->SetBranchAddress("mcHardOutPt", &mcHardOutPt, &b_mcHardOutPt);
     fChain->SetBranchAddress("mcHardOutP", &mcHardOutP, &b_mcHardOutP);
     fChain->SetBranchAddress("mcHardOutEta", &mcHardOutEta, &b_mcHardOutEta);
     fChain->SetBranchAddress("mcHardOutPhi", &mcHardOutPhi, &b_mcHardOutPhi);
     fChain->SetBranchAddress("mcHardOutE", &mcHardOutE, &b_mcHardOutE);
     fChain->SetBranchAddress("mcHardOutPdgId", &mcHardOutPdgId, &b_mcHardOutPdgId);
   }

   //   Notify();

}

Bool_t xAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
     cout<<"Before Crash? "<<endl;
   return kTRUE;
}

void xAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

