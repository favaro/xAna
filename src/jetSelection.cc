#include "interface/xAna.hh"
#include "interface/MinitreeOut.hh"


#include <map>
#include <vector>

using namespace std;

string jecDir = "etc/inputs/jecAK5PF/";

vector<int> xAna::selectJets(Int_t iPho1, Int_t iPho2, Int_t nvtx, Double_t wei, Int_t selVtx) {
  vector<int> jetIndex; 
  for( int i = 0 ; i < nJet; ++i ) { 
    //  if(i>20)continue;
    hAllJetPt->Fill(  (*jetPt)[i],wei);
    hAllJetEta->Fill((*jetEta)[i],wei);
    
    float dR_jetg1 = deltaR((*phoEta)[iPho1], (*phoPhi)[iPho1],(*jetEta)[i],(*jetPhi)[i]);
    float dR_jetg2 = deltaR((*phoEta)[iPho2], (*phoPhi)[iPho2],(*jetEta)[i],(*jetPhi)[i]);
    //cout <<   " dR_jetg1  "<<  dR_jetg1  <<" dR_jetg2  "<<  dR_jetg2   <<endl;
    
    hDRAllJetg1->Fill(dR_jetg1,wei);
    hDRAllJetg2->Fill(dR_jetg2,wei);
    
    if( (*jetPt)[i] > 20 ){
      
      if( dR_jetg1 < 0.5 || dR_jetg2 < 0.5 ) continue;    
      
      hJetNHF->Fill((*jetNHF)[i],wei);
      hJetNEF->Fill((*jetNEF)[i],wei);
      hJetNConst->Fill((*jetNConstituents)[i],wei);
      
      hJetCHF->Fill((*jetCHF)[i],wei);
      hJetNCH->Fill((*jetNCH)[i],wei);
      hJetCEF->Fill((*jetCEF)[i],wei);

      // 2011
      if( nvtx <= 0 ) {     
        if ( fabs((*jetEta)[i]) > 4.7 ) continue;
      }
      // PU Jet ID 2012
      else if( nvtx > 0 ) {
        if ( fabs((*jetEta)[i]) > 4.7 ) continue;
        // if ( (jetWPLevels[i][2] & 4) == 0 ) continue; // [0] simple MVA, [1] full MVA, [2] cutbased, [3] PhilV1; 1(= 1<<0): tight, 2(= 1<<1): medium, 4(= 2<<2): loose /// !! for ntuples after May30
        if (      fabs((*jetEta)[i]) < 2.5  ) { if ((*jetBetaStarClassicExt)[i][selVtx]/log(nvtx-0.64) > 0.2 || (*jetDR2Mean)[i] > 0.06 ) continue; }
        else if ( fabs((*jetEta)[i]) < 2.75 ) { if ((*jetBetaStarClassicExt)[i][selVtx]/log(nvtx-0.64) > 0.3 || (*jetDR2Mean)[i] > 0.05 ) continue; }
        else if ( fabs((*jetEta)[i]) < 3.0  ) { if (                                                         (*jetDR2Mean)[i] > 0.05 ) continue; }
        else                               { if (                                                         (*jetDR2Mean)[i] > 0.055) continue; }
      }

      jetIndex.push_back( i );
      
    }	     
  }   
  return jetIndex;  


}


bool xAna::setupJEC(void) {
  if( JetCorrector != 0 ) return true;

  /// ------------------------------- JEC --------------------------------------------------//
  /// step 1
  // Create the JetCorrectorParameter objects, the order does not matter.
  // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
  //cout << " --- Setup ResJetPar" << endl;
  if( isData ) ResJetPar = new JetCorrectorParameters(jecDir + "FT_53_V21_AN4::All_L2L3Residual_AK5PF.txt"); 
  else         ResJetPar = new JetCorrectorParameters(jecDir + "START53_V23::All_L2L3Residual_AK5PF.txt");
  //cout <<  " --- Setup L1JetPar" << endl;   
  if( isData ) L1JetPar  = new JetCorrectorParameters(jecDir + "FT_53_V21_AN4::All_L1FastJet_AK5PF.txt");
  else         L1JetPar  = new JetCorrectorParameters(jecDir + "START53_V23::All_L1FastJet_AK5PF.txt");
  //cout << " --- Setup L2JetPar" << endl;   
  if( isData ) L2JetPar  = new JetCorrectorParameters(jecDir + "FT_53_V21_AN4::All_L2Relative_AK5PF.txt");
  else         L2JetPar  = new JetCorrectorParameters(jecDir + "START53_V23::All_L2Relative_AK5PF.txt");
  //cout << " --- Setup L3JetPar" << endl;   
  if( isData ) L3JetPar  = new JetCorrectorParameters(jecDir + "FT_53_V21_AN4::All_L3Absolute_AK5PF.txt");
  else         L3JetPar  = new JetCorrectorParameters(jecDir + "START53_V23::All_L3Absolute_AK5PF.txt");
  /// step 2
  //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
  vector<JetCorrectorParameters> vPar;
  vPar.push_back(*L1JetPar);
  vPar.push_back(*L2JetPar);
  vPar.push_back(*L3JetPar);
  if( isData ) vPar.push_back(*ResJetPar);

  /// step 3
  JetCorrector = new FactorizedJetCorrector(vPar);

  return true;
}

// original ------
vector<int> xAna::selectJetsJEC(Int_t iPho1, Int_t iPho2, Int_t nvtx, Int_t selVtx, vector<int>& jetIndex) {
  
  //vector<int> jetIndex;
  jetIndex.clear();
  
  setupJEC();

  vector<int> jetSortedIndex;
  for( int i = 0 ; i < nJet; ++i ) {
    /// step 4
    JetCorrector->setJetEta((*jetEta)[i]);
    JetCorrector->setJetPt((*jetRawPt)[i]);
    JetCorrector->setJetA((*jetArea)[i]);
    JetCorrector->setRho(rho2012); 
    /// step 5
    double correction = JetCorrector->getCorrection();
    // vector<double> factors = JetCorrector->getSubCorrections();

    // cout << "   " << (*jetPt)[i] << "   " << (*jetRawPt)[i]] << "   " << (*jetPt)[i]/(*jetRawPt)[i]] << "   " << correction << endl;
    // cout << "   " << (*jetEn)[i] << "   " << (*jetRawEn)[i] << "   " << (*jetEn)[i]/(*jetRawEn)[i] << "   " << correction << endl;

    /// step 6
    (*jetPt)[i] = (*jetRawPt)[i] * correction;
    (*jetEn)[i] = (*jetRawEn)[i] * correction;
    /// --------------------------------------------------------------------------------------//

    

    float dR_jetg1 = deltaR((*phoEta)[iPho1], (*phoPhi)[iPho1],(*jetEta)[i],(*jetPhi)[i]);
    float dR_jetg2 = deltaR((*phoEta)[iPho2], (*phoPhi)[iPho2],(*jetEta)[i],(*jetPhi)[i]);
    //cout <<   " dR_jetg1  "<<  dR_jetg1  <<" dR_jetg2  "<<  dR_jetg2   <<endl;


    if( (*jetPt)[i] > 20 ){


      if( dR_jetg1 < 0.5 || dR_jetg2 < 0.5 ) continue;


      // 2011
      if( nvtx <= 0 ) {
        if ( fabs((*jetEta)[i]) > 4.7 ) continue;
      }
      // PU Jet ID 2012
      else if( nvtx > 0 ) {
        if ( fabs((*jetEta)[i]) > 4.7 ) continue;
        // if ( (jetWPLevels[i][2] & 4) == 0 ) continue; // [0] simple MVA, [1] full MVA, [2] cutbased, [3] PhilV1; 1(= 1<<0): tight, 2(= 1<<1): medium, 4(= 2<<2): loose /// !! for ntuples after May30
        if (      fabs((*jetEta)[i]) < 2.5  ) { if ((*jetBetaStarClassicExt)[i][selVtx]/log(nvtx-0.64) > 0.2 || (*jetDR2Mean)[i] > 0.06 ) continue; }
        else if ( fabs((*jetEta)[i]) < 2.75 ) { if ((*jetBetaStarClassicExt)[i][selVtx]/log(nvtx-0.64) > 0.3 || (*jetDR2Mean)[i] > 0.05 ) continue; }
        else if ( fabs((*jetEta)[i]) < 3.0  ) { if (                                                         (*jetDR2Mean)[i] > 0.05 ) continue; }
        else                               { if (                                                         (*jetDR2Mean)[i] > 0.055) continue; }
      }

      jetIndex.push_back( i );

    }
  }
    
  int leadJetId    = -99;
  int subleadJetId = -99;
  Float_t leadJetPt    = -99.;
  Float_t subleadJetPt = -99.;

  _minitree->dijet_nJet = (float)jetIndex.size();



  for ( unsigned ijet =0; ijet < jetIndex.size(); ijet++) {
    if      ((*jetPt)[jetIndex[ijet]] > leadJetPt)    { subleadJetPt = leadJetPt;             subleadJetId=leadJetId; 
                                                     leadJetPt    = (*jetPt)[jetIndex[ijet]]; leadJetId=ijet; }
    else if ((*jetPt)[jetIndex[ijet]] > subleadJetPt) { subleadJetPt = (*jetPt)[jetIndex[ijet]]; subleadJetId=ijet;}
  }

  if (leadJetId    > -1) jetSortedIndex.push_back(jetIndex[leadJetId]);
  if (subleadJetId > -1) jetSortedIndex.push_back(jetIndex[subleadJetId]);

 
  return jetSortedIndex;
}


// JM TTH JETS SELECTION

vector<int> xAna::selectJetsTTH( vector<int> goodJetsIndex ,  vector<int> elecIndex, vector<int> muonIndex, int& nBTag ) {

  vector<int> ttHjetIndex; 
  ttHjetIndex.clear();
  nBTag=0;
  
  map<float,int,greater<float> > jetMap; 
  vector<int> ttHjetSortedIndex; 
  
  for( size_t ij = 0; ij < goodJetsIndex.size(); ij++ ) {
   
    float jetBtag = (*jetCombinedSecondaryVtxBJetTags)[goodJetsIndex[ij]];
    if( fabs( (*jetEta)[goodJetsIndex[ij]] ) > 2.4 || (*jetPt)[goodJetsIndex[ij]] < 25) continue; 

    // If leptons, select jets away from leptons  
    // FIX these hardcoded cuts ?
    
    bool pass_dRLep=true;
    
    for(unsigned ie=0;ie<elecIndex.size();ie++){
      double dR = deltaR((*eleEta)[elecIndex[ie]], (*elePhi)[elecIndex[ie]], (*jetEta)[goodJetsIndex[ij]],(*jetPhi)[goodJetsIndex[ij]]);
      if(dR<0.5) pass_dRLep=false;
    }
    
    for(unsigned im=0;im<muonIndex.size();im++){
      double dR = deltaR( (*muEta)[muonIndex[im]], (*muPhi)[muonIndex[im]],(*jetEta)[goodJetsIndex[ij]],(*jetPhi)[goodJetsIndex[ij]]);
      if(dR<0.5) pass_dRLep=false; 
    }
    if(!pass_dRLep) continue;
  
    ttHjetIndex.push_back(goodJetsIndex[ij]);
       
    if( jetBtag > 0.679 )  nBTag++;
    
    // if you want to sort ttH jets by CSV:
    //jetMap.insert( make_pair( jetBtag,  goodJetsIndex[ij]  ) );
    
    // if you want to sort ttH jets by pT:
    jetMap.insert( make_pair( (*jetPt)[goodJetsIndex[ij]], goodJetsIndex[ij] ) );
    
  }
  
  map<float,int>::iterator p; 
  for( p = jetMap.begin(); p != jetMap.end(); ++p ) ttHjetSortedIndex.push_back( p->second );

  //  return ttHjetIndex;
  return ttHjetSortedIndex;
}


// JM VH JETS SELECTION

vector<int> xAna::selectJetsVHLep( vector<int> goodJetsIndex ,  vector<int> elecIndex, vector<int> muonIndex ) {

  vector<int> VHjetIndex; 
  VHjetIndex.clear();
  
  map<float,int,greater<float> > jetMap; 
  vector<int> VHjetSortedIndex; 

  for( size_t ij = 0; ij < goodJetsIndex.size(); ij++ ) {

    if( fabs( (*jetEta)[goodJetsIndex[ij]] ) > 2.4 ) continue; 

    // If leptons, select jets away from leptons  
    
    bool pass_dRLep=true;
    
    for(unsigned ie=0;ie<elecIndex.size();ie++){
      double dR = deltaR((*eleEta)[elecIndex[ie]], (*elePhi)[elecIndex[ie]], (*jetEta)[goodJetsIndex[ij]],(*jetPhi)[goodJetsIndex[ij]]);
      if(dR<0.5) pass_dRLep=false;
    }
    
    for(unsigned im=0;im<muonIndex.size();im++){
      double dR = deltaR( (*muEta)[muonIndex[im]], (*muPhi)[muonIndex[im]],(*jetEta)[goodJetsIndex[ij]],(*jetPhi)[goodJetsIndex[ij]]);
      if(dR<0.5) pass_dRLep=false; 
    }
    if(!pass_dRLep) continue;
  
    VHjetIndex.push_back(goodJetsIndex[ij]);
    
    // sort VH jets by pT:
    jetMap.insert( make_pair( (*jetPt)[goodJetsIndex[ij]], goodJetsIndex[ij] ) );
    
  }
  
  map<float,int>::iterator p; 
  for( p = jetMap.begin(); p != jetMap.end(); ++p ) VHjetSortedIndex.push_back( p->second );

  return VHjetSortedIndex;
}


vector<int> xAna::selectJetsVHHad( vector<int> goodJetsIndex ) {

  vector<int> VHjetIndex; 
  VHjetIndex.clear();
  
  map<float,int,greater<float> > jetMap; 
  vector<int> VHjetSortedIndex; 
  
  for( size_t ij = 0; ij < goodJetsIndex.size(); ij++ ) {
   
    if( fabs( (*jetEta)[goodJetsIndex[ij]] ) > 2.4 || (*jetPt)[goodJetsIndex[ij]] < 40 ) continue; 
  
    VHjetIndex.push_back(goodJetsIndex[ij]);    

    // if you want to sort VH jets by pT:
    jetMap.insert( make_pair( (*jetPt)[goodJetsIndex[ij]], goodJetsIndex[ij] ) );    
  }
  
  map<float,int>::iterator p; 
  for( p = jetMap.begin(); p != jetMap.end(); ++p ) VHjetSortedIndex.push_back( p->second );

  return VHjetSortedIndex;
}



void xAna::dijetSelection( const TLorentzVector & corPho1, const TLorentzVector & corPho2,
					 const vector<int> & goodJetsIndex,
					  int ivtx,
					 int &vbftag, int &hstratag, int &catjet ) {
  vbftag   = 0;
  hstratag = 0;
  catjet   =-1;

  // CF: horrible...but just for quick synch...
  
  if( goodJetsIndex.size() > 0 ) {
    _minitree->dijet_ptj[0] = (*jetPt)[goodJetsIndex[0]];
    _minitree->dijet_jec[0] = (*jetPt)[goodJetsIndex[0]]/(*jetRawPt)[goodJetsIndex[0]];
    _minitree->dijet_etj[0] = (*jetEt)[goodJetsIndex[0]];
    _minitree->dijet_eta[0] = (*jetEta)[goodJetsIndex[0]];
    _minitree->dijet_phi[0] = (*jetPhi)[goodJetsIndex[0]];
    _minitree->dijet_en[0]  = (*jetEn)[goodJetsIndex[0]];
    _minitree->dijet_rms[0] = (*jetDR2Mean)[goodJetsIndex[0]];
    _minitree->dijet_betaStar[0]  = (*jetBetaStarClassicExt)[goodJetsIndex[0]][ivtx];

 
  }
  
  // -------------------------------------
  
  if( goodJetsIndex.size() < 2 ) return;
  int iJet1 = goodJetsIndex[0];
  int iJet2 = goodJetsIndex[1];

  TLorentzVector l_corr_gg;
  l_corr_gg = corPho1 + corPho2 ;

  TLorentzVector jet1;
  TLorentzVector jet2;
  TLorentzVector ljetjet;

  jet1.SetPtEtaPhiE((*jetPt)[iJet1], (*jetEta)[iJet1], (*jetPhi)[iJet1], (*jetEn)[iJet1]);
  jet2.SetPtEtaPhiE((*jetPt)[iJet2], (*jetEta)[iJet2], (*jetPhi)[iJet2], (*jetEn)[iJet2]);
  ljetjet = jet1 + jet2;
  

  float  deltaEta_jj = (*jetEta)[iJet1]- (*jetEta)[iJet2];
  float  deltaPhi_jj  = fabs(jet1.Phi()-jet2.Phi())<3.14259? fabs(jet1.Phi()-jet2.Phi()):6.28-fabs(jet1.Phi()-jet2.Phi());   
  
  float  massJetJet = ljetjet.M();

  //float dR_jet1g1 = deltaR(corPho1.Eta(), corPho1.Phi(),jet1.Eta(),jet1.Phi());
  //float dR_jet1g2 = deltaR(corPho1.Eta(), corPho1.Phi(),jet1.Eta(),jet1.Phi());
  //float dR_jet2g1 = deltaR(corPho1.Eta(), corPho1.Phi(),jet2.Eta(),jet2.Phi());
  //float dR_jet2g2 = deltaR(corPho1.Eta(), corPho1.Phi(),jet2.Eta(),jet2.Phi());

 
  //float dPhi_jet1jet2 = jet1.DeltaPhi(jet2);
  double zeppenfeld = l_corr_gg.Eta() - ((*jetEta)[iJet1]+ (*jetEta)[iJet2])/2.;
  double dPhiDiJetDiPho = ljetjet.DeltaPhi(l_corr_gg);

  ////ProdEtaJets=(*jetEta)[iJet1]*(*jetEta)[iJet2];
  ////dEtaJets=fabs(deltaEta_jj);
  ////ZapJets=zeppenfeld;
  ////MJets=massJetJet;

  // Di-Jet VBF Selection
 
  if ((*jetPt)[iJet1] > 30.) { // dijet VBF cut2
   
    if( fabs(deltaEta_jj) > 3.5){ // dijet VBF cut3
         
      if(fabs(zeppenfeld) < 2.5){ // dijet VBF cut4

	if(massJetJet > 350){ // dijet VBF cut5
	 	 	  
	  if (fabs(dPhiDiJetDiPho) > 2.6) { // dijet VBF cut6
	    
	  } // End of dijet VBF cut6 (dPhiDiJetDiPho > 2.6)
	} // End of dijet VBF cut5 (M_jj>350)
      } // End of dijet VBF cut4 (|Zeppenfeld|<2.5)
    } // End of dijet VBF cut3 (|deltaEta|>3.5)
  } // End of dijet VBF cut2 (lead Jet Pt > 30.)
  // } // End of dijet VBF cut1 (lead Pho Et > 55. * mass/120. && trail Pho Et > 25.)
  
  // Di-Jet VH Selection
  // if (corPho1.Et()>l_corr_gg.M()*60./120. && corPho2.Et()>25.) { // dijet VH cut1

  if ((*jetPt)[iJet1] > 30.) { // dijet VH cut2
    
    if( fabs(deltaEta_jj) < 2.5){ // dijet VH cut3
      
      if(fabs(zeppenfeld)<1.5){ // dijet VH cut4

	if(fabs(massJetJet-85.)<30.){ // dijet VH cut5

          if (corPho1.Pt()/l_corr_gg.M() > 0.5 && corPho2.Pt() > 25.) hstratag = 1;
	  
	} // End of dijet VH cut5 (|M_jj-85.|<30.)
      } // End of dijet VH cut4 (|Zeppenfeld|<1.5)
    } // End of dijet VH cut3 (|deltaEta|<2.5)
  } // End of dijet VH cut2 (lead Jet Pt > 30.)
  // } // End of dijet VH cut1 (lead Pho Et > 60. * mass/120. && trail Pho Et > 25.)

  // Daniele VBF jet category
  vbftag = 0;
  

  /*
  if( corPho1.Pt()/l_corr_gg.M() >60./120. && corPho2.Pt() > 25.) {

    if     ((corPho1.Pt()/l_corr_gg.M()>0.5)&&(corPho2.Pt()>25.)&&((*jetPt)[iJet1]>30.)&&((*jetPt)[iJet2]>30.)
	    &&(fabs(deltaEta_jj)>3.0)&&(fabs(zeppenfeld)< 2.5)&&(massJetJet>500)&&(fabs(dPhiDiJetDiPho)>2.6)) catjet = 0;
    else if((corPho1.Pt()/l_corr_gg.M()>0.5)&&(corPho2.Pt()>25.)&&((*jetPt)[iJet1]>30.)&&((*jetPt)[iJet2]>20.)
	    &&(fabs(deltaEta_jj)>3.0)&&(fabs(zeppenfeld)< 2.5)&&(massJetJet>250)&&(fabs(dPhiDiJetDiPho)>2.6)) catjet = 1;      // CF: ?????
  }
  */
  

  
  _minitree->dijet_dEtaJJ    = deltaEta_jj;
  _minitree->dijet_dPhiJJ    = deltaPhi_jj;
  _minitree->dijet_zeppenfeld = zeppenfeld;
  _minitree->dijet_mJJ        = massJetJet;
  _minitree->dijet_dphiGGJJ   = dPhiDiJetDiPho;

  for( int ij = 0 ; ij < 2; ij++ ) {   
    
    _minitree->dijet_ptj[ij] = (*jetPt)[goodJetsIndex[ij]];
    _minitree->dijet_jec[ij] = (*jetPt)[goodJetsIndex[ij]]/(*jetRawPt)[goodJetsIndex[ij]];
    _minitree->dijet_etj[ij] = (*jetEt)[goodJetsIndex[ij]];
    _minitree->dijet_eta[ij] = (*jetEta)[goodJetsIndex[ij]];
    _minitree->dijet_phi[ij] = (*jetPhi)[goodJetsIndex[ij]];
    _minitree->dijet_en[ij]  = (*jetEn)[goodJetsIndex[ij]];
    _minitree->dijet_rms[ij] = (*jetDR2Mean)[goodJetsIndex[ij]];
    _minitree->dijet_betaStar[ij]  = (*jetBetaStarClassicExt)[goodJetsIndex[ij]][ivtx];

    
  }
  _minitree->dijet_ptj1 = _minitree->dijet_ptj[0];
  _minitree->dijet_ptj2 = _minitree->dijet_ptj[1];
  

  /// compute VBF mva discriminant
  myVBFLeadJEta = (*jetEta)[goodJetsIndex[0]];
  myVBFSubJEta  = (*jetEta)[goodJetsIndex[1]];
  myVBFLeadJPt= (*jetPt)[goodJetsIndex[0]];
  myVBFSubJPt = (*jetPt)[goodJetsIndex[1]];
  myVBF_Mjj   = massJetJet;
  myVBFdEta   = fabs(deltaEta_jj);
  myVBFZep    = fabs(zeppenfeld);
  myVBFdPhi   = fabs(dPhiDiJetDiPho);
  myVBFdPhiTrunc   = min(fabs(dPhiDiJetDiPho),2.916);
  myVBF_Mgg   = l_corr_gg.M();
  myVBFDiPhoPtOverM   = l_corr_gg.Pt() / myVBF_Mgg;
  myVBFLeadPhoPtOverM = corPho1.Pt()   / myVBF_Mgg;
  myVBFSubPhoPtOverM  = corPho2.Pt()   / myVBF_Mgg;
  

  /// reset catjet for mva VBF selection
  if( DiscriVBF_useMvaSel || DiscriVBF_useCombMvaSel) catjet = -1;
  if( myVBFLeadPhoPtOverM > _config->diJetPtg1M() && 
      myVBFSubPhoPtOverM  > _config->diJetPtg2M() &&
      myVBFLeadJPt>30. && myVBFSubJPt>20. &&
      myVBF_Mjj           > _config->diJetMjj() ) {

    _minitree->mtree_fillDijetTree = true;
    _minitree->dijet_mvaVbf = DiscriVBF->EvaluateMVA(DiscriVBF_Method);
    _minitree->dijet_mvaVbf_sig    = 1-DiscriVBF->EvaluateMulticlass(DiscriVBF_Method)[0];
    _minitree->dijet_mvaVbf_dipho  = 1-DiscriVBF->EvaluateMulticlass(DiscriVBF_Method)[1];
    _minitree->dijet_mvaVbf_gluglu = 1-DiscriVBF->EvaluateMulticlass(DiscriVBF_Method)[2];

  

    myVBFDIPHObdt =  _minitree->mtree_diphoMva;  // diphoMva already evaluated in xAna_all
    myVBF_MVA = _minitree->dijet_mvaVbf;
    _minitree->dijet_mvaVbfCombDipho = DiscriVBFCombDipho->EvaluateMVA(DiscriVBF_Method);

    
    if( DiscriVBF_useMvaSel )  {      
      if( myVBFLeadPhoPtOverM>1./3. && myVBFSubPhoPtOverM>1./4. && myVBFLeadJPt>30. && myVBFSubJPt>20. && myVBF_Mjj > 250. ) 
      for( unsigned icat = 0; icat < DiscriVBF_cat.size(); icat++ ) 
	if(  _minitree->dijet_mvaVbf >= DiscriVBF_cat[icat] ) { catjet = icat; break;  }
    }
    if( DiscriVBF_useCombMvaSel )  {      
      if( myVBFLeadPhoPtOverM>1./3. && myVBFSubPhoPtOverM>1./4. && myVBFLeadJPt>30. && myVBFSubJPt>20. && myVBF_Mjj > 250. ) {
	for( unsigned icat = 0; icat < DiscriVBFCombDipho_cat.size(); icat++ ) {
	  catjet = icat; 
	  break;  
	 
	}      
      }
    }
   

  }

  if( catjet >= 0 ) vbftag = 1;
  if( vbftag == 1 ) _minitree->mtree_fillDijetTree = true;

  if( _minitree->mtree_fillDijetTree && catjet < 0 ) catjet = 999;
  

}


//**********************************************************************

float xAna::combineMVAforVBF( const HggEvtCandidate hggevent )
{
  
  float MvaValue = -999;

  if( hggevent.njets != 2 ) cout << " WARNING trying to apply VBF MVA with #jets != 2" << endl;
  
  TLorentzVector jet1, jet2;
  jet1.SetPtEtaPhiE((*jetPt)[hggevent.jetIdx[0]], (*jetEta)[hggevent.jetIdx[0]], (*jetPhi)[hggevent.jetIdx[0]], (*jetEn)[hggevent.jetIdx[0]]);
  jet2.SetPtEtaPhiE((*jetPt)[hggevent.jetIdx[1]], (*jetEta)[hggevent.jetIdx[1]], (*jetPhi)[hggevent.jetIdx[1]], (*jetEn)[hggevent.jetIdx[1]]);
  float deltaEta_jj     = (*jetEta)[hggevent.jetIdx[0]]- (*jetEta)[hggevent.jetIdx[1]];
  double dPhiDiJetDiPho = ( jet1+jet2 ).DeltaPhi( hggevent.gammaK_lead + hggevent.gammaK_trail );  
  float zeppenfeld      = ( hggevent.gammaK_lead + hggevent.gammaK_trail ).Eta() - ((*jetEta)[hggevent.jetIdx[0]]+ (*jetEta)[hggevent.jetIdx[1]])/2.;
  
 /// compute VBF mva discriminant   (ALL VARIABLES DEFINED AS GLOBAL IN MAIN .h)
  myVBFLeadJEta  = (*jetEta)[hggevent.jetIdx[0]];
  myVBFSubJEta   = (*jetEta)[hggevent.jetIdx[1]];
  myVBFLeadJPt   = (*jetPt)[hggevent.jetIdx[0]];
  myVBFSubJPt    = (*jetPt)[hggevent.jetIdx[1]];
  myVBF_Mjj      = ( jet1 + jet2 ).M();
  myVBFdEta      = fabs(deltaEta_jj);
  myVBFZep       = fabs(zeppenfeld);
  myVBFdPhi      = fabs(dPhiDiJetDiPho);
  myVBFdPhiTrunc = min(fabs(dPhiDiJetDiPho),2.916);
  myVBF_Mgg      = ( hggevent.gammaK_lead + hggevent.gammaK_trail ).M();
  myVBFDiPhoPtOverM   = ( hggevent.gammaK_lead + hggevent.gammaK_trail ).Pt() / myVBF_Mgg;
  myVBFLeadPhoPtOverM = hggevent.gammaK_lead.Pt()   / myVBF_Mgg;
  myVBFSubPhoPtOverM  = hggevent.gammaK_trail.Pt()   / myVBF_Mgg;



  if( myVBFLeadPhoPtOverM > _config->diJetPtg1M() && 
      myVBFSubPhoPtOverM  > _config->diJetPtg2M() &&
      myVBFLeadJPt>30. && myVBFSubJPt>20. &&
      myVBF_Mjj           > _config->diJetMjj() ) {

    _minitree->mtree_fillDijetTree = true;
    _minitree->dijet_mvaVbf = DiscriVBF->EvaluateMVA(DiscriVBF_Method);
    _minitree->dijet_mvaVbf_sig    = 1-DiscriVBF->EvaluateMulticlass(DiscriVBF_Method)[0];
    _minitree->dijet_mvaVbf_dipho  = 1-DiscriVBF->EvaluateMulticlass(DiscriVBF_Method)[1];
    _minitree->dijet_mvaVbf_gluglu = 1-DiscriVBF->EvaluateMulticlass(DiscriVBF_Method)[2];   
      
    myVBFDIPHObdt = _minitree->mtree_diphoMva; 
    myVBF_MVA     = _minitree->dijet_mvaVbf;
    _minitree->dijet_mvaVbfCombDipho = DiscriVBFCombDipho->EvaluateMVA(DiscriVBF_Method);
    
  }
  
  if( DiscriVBF_useMvaSel )           MvaValue = _minitree->dijet_mvaVbf;
  else if(  DiscriVBF_useCombMvaSel ) MvaValue = _minitree->dijet_mvaVbfCombDipho;

  return MvaValue;
}



//**********************************************************************

void xAna::fillJetVariablesToMinitree( const HggEvtCandidate hggevent )
{

  TLorentzVector corPho1 = hggevent.gammaK_lead;
  TLorentzVector corPho2 = hggevent.gammaK_trail;

  if( hggevent.njets > 0 ) {
    _minitree->dijet_ptj[0] = (*jetPt)[hggevent.jetIdx[0]];
    _minitree->dijet_jec[0] = (*jetPt)[hggevent.jetIdx[0]]/(*jetRawPt)[hggevent.jetIdx[0]];
    _minitree->dijet_etj[0] = (*jetEt)[hggevent.jetIdx[0]];
    _minitree->dijet_eta[0] = (*jetEta)[hggevent.jetIdx[0]];
    _minitree->dijet_phi[0] = (*jetPhi)[hggevent.jetIdx[0]];
    _minitree->dijet_en[0]  = (*jetEn)[hggevent.jetIdx[0]];
    _minitree->dijet_rms[0] = (*jetDR2Mean)[hggevent.jetIdx[0]];
    _minitree->dijet_betaStar[0]  = (*jetBetaStarClassicExt)[hggevent.jetIdx[0]][hggevent.vertexIdx];    
  }

  if( hggevent.njets == 2 ) {

    int iJet1 = hggevent.jetIdx[0];
    int iJet2 = hggevent.jetIdx[1];

    TLorentzVector l_corr_gg;
    l_corr_gg = corPho1 + corPho2;
    
    TLorentzVector jet1;
    TLorentzVector jet2;
    TLorentzVector ljetjet;
    
    jet1.SetPtEtaPhiE((*jetPt)[iJet1], (*jetEta)[iJet1], (*jetPhi)[iJet1], (*jetEn)[iJet1]);
    jet2.SetPtEtaPhiE((*jetPt)[iJet2], (*jetEta)[iJet2], (*jetPhi)[iJet2], (*jetEn)[iJet2]);
    ljetjet = jet1 + jet2;

    //cout << " from jetSel " << iJet1 << " " << iJet2 << " " << (*jetPt)[iJet1] << " " << (*jetPt)[iJet2] << endl;
    
    float  deltaEta_jj = (*jetEta)[iJet1]- (*jetEta)[iJet2];
    float  deltaPhi_jj  = fabs(jet1.Phi()-jet2.Phi())<3.14259? fabs(jet1.Phi()-jet2.Phi()):6.28-fabs(jet1.Phi()-jet2.Phi());   
    
    float  massJetJet = ljetjet.M();
    
    //float dR_jet1g1 = deltaR(corPho1.Eta(), corPho1.Phi(),jet1.Eta(),jet1.Phi());
    //float dR_jet1g2 = deltaR(corPho1.Eta(), corPho1.Phi(),jet1.Eta(),jet1.Phi());
    //float dR_jet2g1 = deltaR(corPho1.Eta(), corPho1.Phi(),jet2.Eta(),jet2.Phi());
    //float dR_jet2g2 = deltaR(corPho1.Eta(), corPho1.Phi(),jet2.Eta(),jet2.Phi());
    
    // float dPhi_jet1jet2 = jet1.DeltaPhi(jet2);
    
    
    double zeppenfeld = l_corr_gg.Eta() - ((*jetEta)[iJet1]+ (*jetEta)[iJet2])/2.;
    
    double dPhiDiJetDiPho = ljetjet.DeltaPhi(l_corr_gg);
    
    _minitree->dijet_dEtaJJ    = deltaEta_jj;
    _minitree->dijet_dPhiJJ    = deltaPhi_jj;
    _minitree->dijet_zeppenfeld = zeppenfeld;
    _minitree->dijet_mJJ        = massJetJet;
    _minitree->dijet_dphiGGJJ   = dPhiDiJetDiPho;

    for( size_t ij = 0 ; ij < 2; ij++ ) {   
      
      _minitree->dijet_ptj[ij] = (*jetPt)[hggevent.jetIdx[ij]];
      _minitree->dijet_jec[ij] = (*jetPt)[hggevent.jetIdx[ij]]/(*jetRawPt)[hggevent.jetIdx[ij]];
      _minitree->dijet_etj[ij] = (*jetEt)[hggevent.jetIdx[ij]];
      _minitree->dijet_eta[ij] = (*jetEta)[hggevent.jetIdx[ij]];
      _minitree->dijet_phi[ij] = (*jetPhi)[hggevent.jetIdx[ij]];
      _minitree->dijet_en[ij]  = (*jetEn)[hggevent.jetIdx[ij]];
      _minitree->dijet_rms[ij] = (*jetDR2Mean)[hggevent.jetIdx[ij]];
      _minitree->dijet_betaStar[ij]  = (*jetBetaStarClassicExt)[hggevent.jetIdx[ij]][hggevent.vertexIdx];
    
    }
    _minitree->dijet_ptj1 = _minitree->dijet_ptj[0];
    _minitree->dijet_ptj2 = _minitree->dijet_ptj[1];
    
  }

  
  //JM tth jets stuff
  
  _minitree->jets_pt.clear();
  _minitree->jets_en.clear();
  _minitree->jets_et.clear();
  _minitree->jets_eta.clear();
  _minitree->jets_phi.clear();
  _minitree->jets_z.clear();
  _minitree->jets_rms.clear();
  _minitree->jets_betaStar.clear();
  _minitree->jets_bjet_csv.clear();
  _minitree->jets_index.clear();
  
  int njetsTTH=hggevent.njetsTTH;    
  vector<int> goodJetsIndex=hggevent.jetIdxTTH;

  for( unsigned ij = 0; ij < goodJetsIndex.size(); ij++ ) {
    
    _minitree->jets_pt.push_back((*jetPt)[goodJetsIndex[ij]]);
    _minitree->jets_en.push_back((*jetEn)[goodJetsIndex[ij]]);
    _minitree->jets_et.push_back((*jetEt)[goodJetsIndex[ij]]);
    _minitree->jets_eta.push_back((*jetEta)[goodJetsIndex[ij]]);
    _minitree->jets_phi.push_back((*jetPhi)[goodJetsIndex[ij]]);
    _minitree->jets_z.push_back((*jetDZ)[goodJetsIndex[ij]]);
    _minitree->jets_rms.push_back((*jetDR2Mean)[goodJetsIndex[ij]]);
    _minitree->jets_betaStar.push_back((*jetBetaStarClassicExt)[goodJetsIndex[ij]][hggevent.vertexIdx]);
    _minitree->jets_bjet_csv.push_back((*jetCombinedSecondaryVtxBJetTags)[goodJetsIndex[ij]]);
  }
  
  _minitree->jets_njets = njetsTTH;
  _minitree->jets_nbtag=hggevent.nBjets;
  _minitree->jets_index=goodJetsIndex;

}


