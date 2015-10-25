#include "interface/xAna.hh"
#include "interface/MinitreeOut.hh"

using namespace std;


int xAna::photonID(int ipho, int ivtx, float corEt_ipho, bool invEleVeto ) {
  //  cout <<"ipho  " << ipho   << "    invEleVeto  "<<  invEleVeto <<endl;
  if( _config->analysisType() == "MVA"      ) return mvaPhotonID(ipho,ivtx,corEt_ipho,invEleVeto);
  if( _config->analysisType() == "baseline" ) return cicPhotonID(ipho,ivtx,corEt_ipho);
  if( _config->analysisType() == "baselineCiC4PF" ) {
    if( !mvaPhotonID(ipho,ivtx,corEt_ipho,invEleVeto     ) ) return 0;
    return cic4PFPhotonID(ipho,ivtx,corEt_ipho,invEleVeto);  
  }
  
  return 0;
}


int xAna::cic4PFPhotonID(int ipho, int ivtx, float corEt_ipho, bool invEleVeto){
    
  int passID = 1;
  float Et = (*phoEt)[ipho]; 
  if(corEt_ipho != -1) {
    Et = corEt_ipho;
  }
 
   
  Et = (*phoE)[ipho] / cosh((*phoEtaVtx)[ipho][ivtx]);
  float effectiveArea         = 0.09;
  float effectiveAreaWorstVtx = 0.23;
  float pfiso_charged_chosen_vertex = (*phoCiCPF4chgpfIso03)[ipho][ivtx];
  float pfiso_photon03 = (*phoCiCPF4phopfIso03)[ipho];
  float pfiso_photon04 = (*phoCiCPF4phopfIso04)[ipho];

  int badVtx = 0;
  float  pfiso_charged_worst_vertex = -99;
  for( int iv=0; iv<nVtxBS ; iv++) {
    if((*phoCiCPF4chgpfIso04)[ipho][iv] > pfiso_charged_worst_vertex) {
      pfiso_charged_worst_vertex = (*phoCiCPF4chgpfIso04)[ipho][iv];
      badVtx = iv;
    }
  }
  
  /// replica of the bug in gglobe
  //  float pho_tkiso_badvtx = worstSumTrackPtInCone(ipho,badVtx,0,0.4,0.02,0.0,1.0,0.1);
  
  float pfiso_charged_over_et        = (pfiso_charged_chosen_vertex)*50./Et;  
  float isosum_chosen_vertex_over_et = (pfiso_charged_chosen_vertex + pfiso_photon03 + 2.5 - rho2012*effectiveArea)*50./Et;  
  float isosum_worst_vertex_over_et = 
    (pfiso_charged_worst_vertex + pfiso_photon04 + 2.5 - rho2012*effectiveAreaWorstVtx)*50./( (*phoE)[ipho]/cosh((*phoEtaVtx)[ipho][badVtx]) );


  int cat = phocat(ipho);
  if( isosum_chosen_vertex_over_et > vcicST[0][cat] ) return 0;
  if( isosum_worst_vertex_over_et  > vcicST[1][cat] ) return 0;
  if( pfiso_charged_over_et  > vcicST[2][cat] ) return 0;
  if( (*phoSigmaIEtaIEta)[ipho] > vcicST[3][cat] ) return 0;
  if( (*phoHoverE)[ipho]        > vcicST[4][cat] ) return 0;
  if( (*phoR9)[ipho]            < vcicST[5][cat] ) return 0;

  if ((*phoCiCdRtoTrk)[ipho] == 0) (*phoCiCdRtoTrk)[ipho] = 99;
  //  hDRtoTrk->Fill((*phoCiCdRtoTrk)[ipho],wei);
  
  if( _config->noElectronVeto() ) return passID;

  int isElec = (*phoEleVeto)[ipho];
  if(  invEleVeto && isElec == 0 ) return 0;
  if( !invEleVeto && isElec == 1 ) return 0;

  // if(_config->isOnlyLepTag()){
  //   if( !_config->invertElectronVeto() && !invEleVeto_){
  //     //  if( (*phoEleVeto)[ipho]==1 || (*phoCiCdRtoTrk)[ipho]<1  ||  phohasPixelSeed[ipho]==1 )return 0;
  //    if( (*phoEleVeto)[ipho]==1)return 0;
  //   }
  //   if(invEleVeto_ ||  _config->invertElectronVeto()){
  //     // if( (*phoEleVeto)[ipho]==0  &&  (*phoCiCdRtoTrk)[ipho]>1  && phohasPixelSeed[ipho]==0) return 0;
  //     if( (*phoEleVeto)[ipho]==0  ) return 0;
  //   }  
  //}
  return passID;  
}

void xAna::fillPhotonVariablesToMiniTree(int ipho, int ivtx, int index ){ 

  float Et = (*phoE)[ipho]/cosh((*phoEtaVtx)[ipho][ivtx]);


  _minitree->mtree_enScale[index]    = (*phoEnergyScale)[ipho]; // CF
  _minitree->mtree_phoSCRawE[index]  = (*phoSCRawE)[ipho];      // CF
  _minitree->mtree_phoRegrE[index]   = (*phoRegrE)[ipho];      // CF
  _minitree->mtree_r9[index]    = (*phoR9)[ipho];
  _minitree->mtree_sceta[index] = (*phoSCEta)[ipho];
  _minitree->mtree_eg[index]    = (*phoE)[ipho];
  _minitree->mtree_ptg[index]   = (*phoE)[ipho]/cosh((*phoEtaVtx)[ipho][ivtx]);

  _minitree->mtree_eta[index]   = (*phoEtaVtx)[ipho][ivtx];
  _minitree->mtree_phi[index]   = (*phoPhiVtx)[ipho][ivtx];
  /* /// debugging */
  /* _minitree->mtree_eta[index]   = phoSeedDetId1[ipho]; */
  /* _minitree->mtree_phi[index]   = phoSeedDetId2[ipho];; */
  
  _minitree->mtree_relResOverE[index] = (*phoRegrErr)[ipho]/(*phoE)[ipho]; 
  _minitree->mtree_relSmearing[index] = (*phoRegrSmear)[ipho]/(*phoE)[ipho]; 

  //cout << " from photonID [" << ipho << "]  phoE = " << (*phoE)[ipho] << " phoRegrS = " << (*phoRegrSmear)[ipho] << " relSmearing = " << _minitree->mtree_relSmearing[index] << endl; 


  _minitree->mtree_coviEiE[index] = (*phoSigmaIEtaIEta)[ipho];
  _minitree->mtree_HoE[index]     = (*phoHoverE)[ipho];
  _minitree->mtree_drTrk[index]   = (*phoCiCdRtoTrk)[ipho];
  
  _minitree->mtree_cat[index]     = phocatCIC(ipho); // (r9 < 0.90 / r9>0.90 )

  /// preselection photon id cuts
  _minitree->mtree_isElec[index]            = (*phoEleVeto)[ipho];
  _minitree->mtree_etCorrEcalIso[index]     = (*phoEcalIsoDR03)[ipho] - 0.012*Et;
  _minitree->mtree_etCorrHcalIso[index]     = (*phoHcalIsoDR03)[ipho] - 0.005*Et;
  _minitree->mtree_etCorrTrkIso[index]      = (*phoTrkIsoHollowDR03)[ipho] - 0.002*Et;
  _minitree->mtree_pfChgIso02_selVtx[index] = (*phoCiCPF4chgpfIso02)[ipho][ivtx];

  /// mva photon id specific 2012
  _minitree->mtree_coviEiP[index]      = (*phoSigmaIEtaIPhi)[ipho];
  (*phoS4ratio)[ipho] = (*phoE2x2)[ipho]/(*phoE5x5)[ipho];
  _minitree->mtree_s4Ratio[index]      = (*phoS4ratio)[ipho];
  _minitree->mtree_esEffSigmaRR[index] = (*phoESEffSigmaRR_x)[ipho];

  //// PF isolation dedicated to CiCpf id (also used in photon Id MVA )
  /// compute isolation
  float effectiveArea         = 0.09;
  float effectiveAreaWorstVtx = 0.23;
  float pfiso_charged_chosen_vertex = (*phoCiCPF4chgpfIso03)[ipho][ivtx];
  float pfiso_photon03 = (*phoCiCPF4phopfIso03)[ipho];
  float pfiso_photon04 = (*phoCiCPF4phopfIso04)[ipho]; 
  
  int badVtx03   = 0;
  int badVtx04   = 0;
  float  pfiso_charged03_worst_vertex = -99;
  float  pfiso_charged04_worst_vertex = -99;
  for( int iv=0; iv<nVtxBS ; iv++) {
    // worst vtx for chIso03  
      
  if((*phoCiCPF4chgpfIso03)[ipho][iv] > pfiso_charged03_worst_vertex) {
      pfiso_charged03_worst_vertex = (*phoCiCPF4chgpfIso03)[ipho][iv];
      badVtx03 = iv;
    }
   
    // worst vtx for chIso04
    if((*phoCiCPF4chgpfIso04)[ipho][iv] > pfiso_charged04_worst_vertex) {
      pfiso_charged04_worst_vertex = (*phoCiCPF4chgpfIso04)[ipho][iv];
      badVtx04 = iv;
    }
  } 

  float pfiso_charged_over_et        = (pfiso_charged_chosen_vertex)*50./(*phoE)[ipho]/cosh((*phoEtaVtx)[ipho][ivtx]);
  float isosum_chosen_vertex_over_et = (pfiso_charged_chosen_vertex  + pfiso_photon03 + 2.5 
					- rho2012*effectiveArea)*50./Et;
  float isosum_worst_vertex_over_et  = (pfiso_charged04_worst_vertex + pfiso_photon04 + 2.5 
					- rho2012*effectiveAreaWorstVtx)*50./(*phoE)[ipho]/cosh((*phoEtaVtx)[ipho][badVtx04]);
  
  _minitree->mtree_phoIso1[index] = isosum_chosen_vertex_over_et; // CombinedPFIso1
  _minitree->mtree_phoIso2[index] = isosum_worst_vertex_over_et;  // CombinedPFIso2
  _minitree->mtree_phoIso3[index] = pfiso_charged_over_et;        // CombinedPFIso3

  /// particle flow photon isolation
  _minitree->mtree_pfPhoIso03[index] = (*phoCiCPF4phopfIso03)[ipho]; /// mvaPhotonID
  _minitree->mtree_pfPhoIso04[index] = (*phoCiCPF4phopfIso04)[ipho]; 
  _minitree->mtree_pfChgIso03_selVtx[index] = (*phoCiCPF4chgpfIso03)[ipho][ivtx];     /// mvaPhotonID
  _minitree->mtree_pfChgIso03_badVtx[index] = (*phoCiCPF4chgpfIso03)[ipho][badVtx03]; /// mvaPhotonID
  _minitree->mtree_pfChgIso03_defVtx[index] = (*phoCiCPF4chgpfIso03)[ipho][0];

  _minitree->mtree_pfChgIso04_selVtx[index] = (*phoCiCPF4chgpfIso04)[ipho][ivtx];
  _minitree->mtree_pfChgIso04_badVtx[index] = (*phoCiCPF4chgpfIso04)[ipho][badVtx04];
  _minitree->mtree_pfChgIso04_defVtx[index] = (*phoCiCPF4chgpfIso04)[ipho][0];

  PhoID_MVA(ipho,ivtx);
  _minitree->mtree_mvaid[index] = -1;
  if( fabs((*phoSCEta)[ipho]) < 1.45 ) _minitree->mtree_mvaid[index] = phoID_mva[0]->EvaluateMVA("BDT");
  else                              _minitree->mtree_mvaid[index] = phoID_mva[1]->EvaluateMVA("BDT");

  _minitree->mtree_mvaid1 = _minitree->mtree_mvaid[0];
  _minitree->mtree_mvaid2 = _minitree->mtree_mvaid[1];
  _minitree->mtree_eta1   = _minitree->mtree_sceta[0];
  _minitree->mtree_eta2   = _minitree->mtree_sceta[1];
  _minitree->mtree_r91    = _minitree->mtree_r9[0];
  _minitree->mtree_r92    = _minitree->mtree_r9[1];

  ///------ synch only variables.
  if( _minitree->addSyncVariables ) {
    int iconv = matchPhotonToConversion(ipho);
    _minitree->convInd[index] = iconv;
    _minitree->convNtk[index] = iconv>=0 ? (*convNTracks)[iconv] : -1;
    _minitree->convValidVtx[index] = iconv>=0 && (*convNTracks)[iconv] == 2 ? (*convValidVtx)[iconv] : -1;
    _minitree->convChi2Prob[index] = iconv>=0 && (*convNTracks)[iconv] == 2 ? (*convChi2Probability)[iconv] : -1;
    if( iconv >= 0 ) {
      TVector3 qconv;
      if     ( (*convNTracks)[iconv] == 2 ) qconv.SetXYZ( (*convRefittedMomentum_x)[iconv],
						       (*convRefittedMomentum_y)[iconv],
						       (*convRefittedMomentum_z)[iconv] );
      else if( (*convNTracks)[iconv] == 1 ) qconv.SetXYZ( (*convPairMomentum_x)[iconv],
						       (*convPairMomentum_y)[iconv],
						       (*convPairMomentum_z)[iconv] ); 
      _minitree->convPt[index] = qconv.Pt();	
    }
  }
  _minitree->phoID_covIEtaIPhi[index]    = (*phoSigmaIEtaIPhi)[ipho];   
  _minitree->phoID_S4[index]             = (*phoS4ratio)[ipho];        // CF
  _minitree->phoID_ESEffSigmaRR[index]   = (*phoESEffSigmaRR_x)[ipho];
  _minitree->phoID_pfPhoIso03[index]     = (*phoCiCPF4phopfIso03)[ipho];
  _minitree->phoID_pfChIso03[index]      = (*phoCiCPF4chgpfIso03)[ipho][ivtx];
  _minitree->phoID_pfChIso03worst[index] = pfiso_charged03_worst_vertex;
  _minitree->phoID_E2x2[index]           = (*phoE2x2)[ipho];
  _minitree->phoID_E5x5[index]           = (*phoE5x5)[ipho];

  
}


int xAna::cicPhotonID(int ipho, int ivtx, float corEt_ipho ) {
  int passID = 1;
  float Et = (*phoEt)[ipho]; 
  if(corEt_ipho != -1) {
    Et = corEt_ipho;
  }
  
  float rhoAbsFactor = 1;
  //  if( DoSetRho0_ ) rhoAbsFactor = 0;

  float rhofracbad = 0.52*rhoAbsFactor, rhofrac = 0.17*rhoAbsFactor, combIso = 0*rhoAbsFactor, trkIso = 0*rhoAbsFactor;
  int cat = phocat(ipho);
  
  /*
  // find out the tracker isolation with "worst vertex"
  float maxTrkIso = -99;
  for (int i=0; i<nVtxBS; ++i) {
  if ((*phoCiCPF4chgpfIso04)[ipho][i] >  maxTrkIso) {
  maxTrkIso = (*phoCiCPF4chgpfIso04)[ipho][i];
  }
  }
  
  // relative combIso (bad vertex)
  combIso = (maxTrkIso + phoEcalIsoDR04[ipho] + phoHcalIsoDR04[ipho] - rho25*rhofracbad) * 50./Et;
  */
  
  // relative combIso (good vertex)
  combIso = ((*phoCiCPF4chgpfIso03)[ipho][ivtx] + (*phoEcalIsoDR03)[ipho] + (*phoHcalIsoDR04)[ipho] - rho25*rhofrac) * 50./Et; 
  if (combIso > vcicSTmva[0][cat]) return 0;
  
  
  // find out the tracker isolation with "worst vertex"
  float maxTrkIso = -99;
  int badVtx = 0;
  for (int i=0; i<nVtxBS; ++i) {
    if ((*phoCiCPF4chgpfIso03)[ipho][i] >  maxTrkIso) {
      maxTrkIso = (*phoCiCPF4chgpfIso04)[ipho][i];
      badVtx = i;
    }
  }
  
  // relative combIso (bad vertex)
  combIso = (maxTrkIso + (*phoEcalIsoDR04)[ipho] + (*phoHcalIsoDR04)[ipho] - rho25*rhofracbad) * 50./((*phoE)[ipho]/cosh((*phoEtaVtx)[ipho][badVtx]));
  //combIso = (maxTrkIso + phoEcalIsoDR04[ipho] + phoHcalIsoDR04[ipho] - rho25*rhofracbad) * 50./phoEtVtx[ipho][badVtx];
  //combIso = (maxTrkIso + phoEcalIsoDR04[ipho] + phoHcalIsoDR04[ipho] - rho25*rhofracbad) * 50./Et;
  if (combIso > vcicSTmva[1][cat]) return 0;
  
  // relative trkIso (good vertex)
  trkIso = (*phoCiCPF4chgpfIso03)[ipho][ivtx] * 50./Et;
  if (trkIso > vcicSTmva[2][cat]) return 0;
  
  // Sigma_IetaIeta
  if ((*phoSigmaIEtaIEta)[ipho] > vcicSTmva[14][cat]) return 0;
  // H/E
  if ((*phoHoverE)[ipho] > vcicSTmva[15][cat]) return 0;
  // R9
  if ((*phoR9)[ipho] < vcicSTmva[5][cat])  return 0;
  // dR to electron track
  if ((*phoCiCdRtoTrk)[ipho] == 0) (*phoCiCdRtoTrk)[ipho] = 99;
  
  if(  _config->invertElectronVeto() && (*phoCiCdRtoTrk)[ipho] > vcicSTmva[6][cat] ) return 0;
  if( !_config->invertElectronVeto() && (*phoCiCdRtoTrk)[ipho] < vcicSTmva[6][cat] ) return 0;
  
  return passID;
}

int xAna::mvaPhotonID(int ipho, int ivtx, float corEt_ipho, bool invEleVeto ) {


  
  // cout <<" preselection     ipho  "<< ipho <<"  invEleVeto  "<<   invEleVeto <<endl;
  int passID = 1;
  int cat = phocatCIC(ipho);

  float Et = (*phoEt)[ipho]; 
  if(corEt_ipho != -1) {
    Et = corEt_ipho;  
  }

  Et = (*phoE)[ipho] / cosh((*phoEtaVtx)[ipho][ivtx]);
  //  float EtCorrEcalIso = (*phoEcalIsoDR03)[ipho]      - 0.012*Et;
  float EtCorrHcalIso = (*phoHcalIsoDR03)[ipho]      - 0.005*Et;
  float EtCorrTrkIso  = (*phoTrkIsoHollowDR03)[ipho] - 0.002*Et;
  float  trkPfIso02   = (*phoCiCPF4chgpfIso02)[ipho][ivtx];
  int isElec =  (*phoEleVeto)[ipho];//isMatchedToAPromptElec(ipho, dr );
  //Notify();

  if( (*phoSigmaIEtaIEta)[ipho]  >= vcicSTmva[ 3][cat] ) return 0;
  if( (*phoHoverE)[ipho]         >= vcicSTmva[ 4][cat] ) return 0;
  //if( EtCorrEcalIso           >= vcicSTmva[ 8][cat] ) return 0;
  if( EtCorrHcalIso           >= vcicSTmva[ 9][cat] ) return 0;
  if( EtCorrTrkIso            >= vcicSTmva[10][cat] ) return 0;


  if(run==191247  && event== 550253379){
    cout << " invEleVeto:  " <<  invEleVeto << " (*phoEleVeto)[ipho]  "   <<  (*phoEleVeto)[ipho] 
	 << " isElec "       << isElec      << " (*phoCiCdRtoTrk)[ipho] " <<  (*phoCiCdRtoTrk)[ipho] <<endl;   
  }


  /// 2012 preselection only
  if( trkPfIso02            >=  vcicSTmva[16][cat] ) return 0;

  /// 2011 preselection only
  //  float PUCorrHcalEcal       = EtCorrEcalIso + EtCorrHcalIso - 0.17*rho25;
  //  float HollowConeTrkIsoDr03 = phoTrkIsoHollowDR03[ipho];
  /// the 2 are not exactly identical though very very close
  /// the second gives slightly closer results to gglobe
  // float AbsTrkIsoCIC         = sumTrackPtInCone(ipho,ivtx,0,0.3,0.02,0.0,1.0,0.1);
  //  float AbsTrkIsoCIC         = phoCiCTrkIsoDR03[ipho][ivtx];
  //  if( PUCorrHcalEcal         >= vcicSTmva[11][cat] ) return 0;
  //  if( AbsTrkIsoCIC           >= vcicSTmva[12][cat] ) return 0;
  //  if( HollowConeTrkIsoDr03   >= vcicSTmva[13][cat] ) return 0;

  if ((*phoCiCdRtoTrk)[ipho] == 0) (*phoCiCdRtoTrk)[ipho] = 99;

  // if( (*phoCiCdRtoTrk)[ipho] > 0.15 ) (*phoCiCdRtoTrk)[ipho] = 0.15;

  // float dr = -1;
  //  hDRtoTrk->Fill((*phoCiCdRtoTrk)[ipho],wei);

  if( _config->noElectronVeto() ) return passID;
    
  if(  invEleVeto && isElec == 0 ) return 0;
  if( !invEleVeto && isElec == 1 ) return 0;

  return passID;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void xAna::ApplyEnergyCorrections( bool DoOverSmearing, int i ) {
  /// this is per photon energy correction

  float phoEnScaletmp = -999; // filled only if data

  //// ********************* define S4 variable **********************////
  if( _config->setup() == "Legacy2013_7TeV" ) (*phoS4ratio)[i] = 1;
  else                                        (*phoS4ratio)[i] = (*phoE2x2)[i] / (*phoE5x5)[i];
  //// ************************************************************* ////
  
  if( fabs((*phoSCEta)[i]) <= 2.5 ) {
    
    // specific corrections for MC: **********
    if( !isData) {
      
      float smearing = 0;
      overSmear_.setSmearSeed(event, run, lumis, (*phoSCPhi)[i], (*phoSCEta)[i]);
      if (DoOverSmearing) smearing = overSmear_.randOverSmearing((*phoSCEta)[i],(*phoR9)[i], (*phoRegrE)[i]/cosh((*phoSCEta)[i]),isInGAP_EB(i),overSmearSyst_);
      //note for Synchronization purposes the above is centered at 1.0
      float ecorr, esigma;
      CorrLikelihood(i, ecorr, esigma,true);
      (*phoRegrE)[i]   = ecorr;
      (*phoRegrErr)[i] = esigma;
      (*phoRegrE)[i]  *= (smearing); 
      (*phoE)[i]      *= (smearing); 
      /// from MassFactorized in gglobe:   energyCorrectedError[ipho] *=(l.pho_isEB[ipho]) ? 1.07 : 1.045 ;
      float smearFactor = 1;
      if( _config->setup() == "Legacy2013_7TeV" ) smearFactor = fabs((*phoSCEta)[i]) < 1.45 ? 1.07: 1.045;
      (*phoRegrErr)[i] *= smearFactor;
      
      rescalePhoIdVar(i);      
    } // end of if is MC  ***********
    
    // specific corrections for data: **********
    
    if( isData ){
      
      float enCorrSkim = 1;//enScaleSkimEOS.energyScale( phoR9[i], phoSCEta[i], run);
      float phoEnScale = enScale.energyScale( (*phoR9)[i], (*phoSCEta)[i],(*phoE)[i]/cosh((*phoEta)[i]), run, 0)/enCorrSkim;
      
      float ecorr, esigma;
      CorrLikelihood(i, ecorr, esigma, true);
      (*phoRegrE)[i]   = ecorr;
      (*phoRegrErr)[i] = esigma*phoEnScale;  // CF multiply the regression error by enscale for consistency with globe
      (*phoRegrE)[i]  *= phoEnScale;
      (*phoE)[i]      *= phoEnScale;
      phoEnScaletmp    = phoEnScale;  // will make sense only for data

    } // end of if is data **********
    
    (*phoStdE)[i]        = (*phoE)[i];
    (*phoE)[i]           = (*phoRegrE)[i];	 
    (*phoEnergyScale)[i] = phoEnScaletmp;  // will make sense only for data
    
    /// transform calo position abd etaVtx, phiVtx with SC position
    (*phoCaloPos_x)[i] = (*phoSCPos_x)[i];
    (*phoCaloPos_y)[i] = (*phoSCPos_y)[i];
    (*phoCaloPos_z)[i] = (*phoSCPos_z)[i];
    for( int ivtx = 0 ; ivtx < nVtxBS; ivtx++ ) {
      TVector3 xxi = getCorPhotonTVector3(i,ivtx);
      (*phoEtaVtx)[i][ivtx] = xxi.Eta();
      (*phoPhiVtx)[i][ivtx] = xxi.Phi();
    }
    
    /// additionnal smearing to go to data energy resolution
    (*phoEt)[i]        = (*phoE)[i] / cosh((*phoEta)[i]);  
    //(*phoRegrSmear)[i] = (*phoE)[i]*overSmear_.meanOverSmearing((*phoSCEta)[i],(*phoR9)[i],(*phoEt)[i],isInGAP_EB(i),0);
    (*phoRegrSmear)[i] = (*phoE)[i]*overSmear_.meanOverSmearing((*phoSCEta)[i],(*phoR9)[i],isInGAP_EB(i),0);

  } // end of eta cut
}

int xAna::phocatCIC(int i){
  int r9 = ((*phoR9)[i] > 0.9) ? 0 : 1;
  int eta = (fabs((*phoSCEta)[i]) < 1.479) ? 0 : 1;

  return r9 + 2*eta;
}

int xAna::phocat(int i) {

  int r9 = ((*phoR9)[i] > 0.94) ? 0 : 1;
  int eta = (fabs((*phoSCEta)[i]) < 1.479) ? 0 : 1;

  return r9 + 2*eta;
}



bool xAna::isInGAP_EB( int ipho ) {
  if( abs((*phoSCEta)[ipho]) > 1.45 ) return false;
  //  if( phoPos[ipho] == 2 || phoPos[ipho] == 4 ) return true;
  int detid1 = abs((*phoSeedDetId1)[ipho]);
  int detid2 = abs((*phoSeedDetId2)[ipho]);
  if( fabs(detid1-  0.5) < 5 || fabs(detid1- 25.5) < 5 || fabs(detid1- 45.5) < 5 ||
      fabs(detid1- 65.5) < 5 || fabs(detid1- 85.5) < 5 ) return true;

  if( fabs( detid2%20-10.5 ) > 5 ) return true;
  return false;
}

bool xAna::preselectPhoton(int i, float Et){

  bool rc = true;

  if (fabs((*phoSCEta)[i]) > 1.4442 && fabs((*phoSCEta)[i]) < 1.566) return false; 
  if (fabs((*phoSCEta)[i]) > 2.5) return false; 

  if (Et < 20) return false; 

  return rc;

}

TVector3 xAna::getCorPhotonTVector3(int i, int selVtx) {
  float vtx_X, vtx_Y, vtx_Z;
  vtx_X = (*vtxbs_x)[selVtx];
  vtx_Y = (*vtxbs_y)[selVtx];
  vtx_Z = (*vtxbs_z)[selVtx];
  TVector3 vtxVec(vtx_X,vtx_Y,vtx_Z);

  // float pho1_X, pho1_Y, pho1_Z;
  // pho1_X = phoCaloPos[i][0];
  // pho1_Y = phoCaloPos[i][1];
  // pho1_Z = phoCaloPos[i][2];
  // TVector3 phoVec(pho1_X,pho1_Y,pho1_Z);
  // // Temporary SC pos for old ggNtuple
  // phoVec.SetPtEtaPhi(sqrt(pow(pho1_X,2)+pow(pho1_Y,2)), phoSCEta[i], phoSCPhi[i]);

  // For new version ggNtuple with SC pos
  float pho1_X = (*phoSCPos_x)[i];
  float pho1_Y = (*phoSCPos_y)[i];
  float pho1_Z = (*phoSCPos_z)[i];
  TVector3 phoVec(pho1_X,pho1_Y,pho1_Z);


  return phoVec - vtxVec;
}





//--------------------------------------------------------------------------------------------------------------
///--------------------------- FC: some old code not use anymore (does not even know if this is really working
//--------------------------------------------------------------------------------------------------------------

float xAna::etaTransformation(  float EtaParticle , float vertexZ)  {


  //---Definitions for ECAL
  const float R_ECAL           = 136.5;
  const float Z_Endcap         = 328.0;
  const float etaBarrelEndcap  = 1.479; 
   
  //---ETA correction

  float Theta = 0.0  ; 
  float ZEcal = R_ECAL*sinh(EtaParticle)+vertexZ;

  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
  if(Theta<0.0) Theta = Theta+ TMath::Pi() ;
  double ETA = - log(tan(0.5*Theta));
         
  if( fabs(ETA) > etaBarrelEndcap )
    {
      float Zend = Z_Endcap ;
      if(EtaParticle<0.0 )  Zend = -Zend ;
      float Zlen = Zend - vertexZ ;
      float RR = Zlen/sinh(EtaParticle); 
      Theta = atan(RR/Zend);
      if(Theta<0.0) Theta = Theta+ TMath::Pi() ;
      ETA = - log(tan(0.5*Theta));		      
    } 
  //---Return the result
  return ETA;
  //---end
}


bool xAna::isAGoodConv( int iconv  ) {
  /// take this part of the code from CMSSW ConvertionTools

  float probMin = 1e-6;
  float lxymin  = 2.0;

  /// if it is matched to a conversion, check its quality
  /// vertex validity is there by default
  //  if( !phoConvValidVtx[ipho]          ) return false;

  if( (*convChi2Probability)[iconv] < probMin ) return false; 

  //compute transverse decay length
  //  XYZVector mom(convRefittedMomentum_x[iconv],convRefittedMomentum_y[iconv],convRefittedMomentum_z[iconv]);
  double xx = (*convRefittedMomentum_x)[iconv];
  double yy = (*convRefittedMomentum_y)[iconv];
  double dbsx = (*convVtx_x)[iconv] - bspotPos[0];
  double dbsy = (*convVtx_y)[iconv] - bspotPos[1];
  double lxy = (xx*dbsx + yy*dbsy)/sqrt(xx*xx+yy*yy);
  
  if( lxy < lxymin ) return false;

  if( (*convNHitsBeforeVtx_0)[iconv] > 0 ) return false;
  if( (*convNHitsBeforeVtx_1)[iconv] > 0 ) return false;
    
  return true;
}

bool xAna::isMatchedToAPromptElec( int ipho, float &dRclosestElec ) {
  /// take this part of the code from CMSSW ConvertionTools
  TLorentzVector pele;  
  bool match = false;
  dRclosestElec = 999;
  for( int iele = 0 ; iele < nEle; iele++ ) {
    if( (*eleMissHits)[iele] > 0 ) continue;

    float deta = fabs( (*eleSCEta)[iele] - (*phoSCEta)[ipho] );
    float dphi = acos( cos( (*eleSCPhi)[iele] - (*phoSCPhi)[ipho] ) );
    float dR   = sqrt(deta*deta + dphi*dphi );
    if( dR < dRclosestElec ) dRclosestElec = dR;

    /// matched the photon and the electron from SC position
    if( deta > 0.01 ) continue;
    if( dphi > 0.01 ) continue;
    
    /// if this is matched to an electron then check if it is conversion
    /// this is slightly different from the original code
    /// but I have only these information here (should bre really close though)
    // if( !(*phoIsConv)[ipho] ) {
    //   /// if the photon has no conversion track, then match
    //   match = true;
    //   continue;
    // } else {
    //   /// else try to find a conversion that would match the electron and check it
    //   pele.SetPtEtaPhiM( elePt[iele], eleEta[iele], elePhi[iele], 0 );
    //   bool isConv = false;
    //   for( int iconv = 0; iconv < nConv; iconv++ ) {
    // 	TLorentzVector pconv;
    // 	pconv.SetPxPyPzE( convP4[iconv][0],convP4[iconv][1],convP4[iconv][2],convP4[iconv][3]);
    // 	if( pconv.Pt() < 1 ) continue;
    // 	if( pele.DeltaR( pconv ) < 0.1 && isAGoodConv( iconv ) ) {isConv = true; break;}
    //   }
    //   if( isConv ) continue;
    //   match = true;
    // }
    match = true;
  }

  if( match ) {
    /// now check this is not a conversion
    /// try something different: see conversion tools
    TVector3 calpos((*phoCaloPos_x)[ipho], (*phoCaloPos_y)[ipho],(*phoCaloPos_z)[ipho]);
    bool isConv = false;
    for( int iconv = 0; iconv < nConv; iconv++ ) { 
      TVector3 convP3( (*convRefittedMomentum_x)[iconv], (*convRefittedMomentum_y)[iconv], (*convRefittedMomentum_z)[iconv]) ;
      TVector3 convVTX((*convVtx_x)[iconv], (*convVtx_y)[iconv],(*convVtx_z)[iconv]);
      TVector3 phoFromConvVtx = calpos - convVTX;
      float drConv = phoFromConvVtx.DeltaR( convP3 ); 
      if( drConv < 0.1 && isAGoodConv( iconv ) ) { isConv = true; break;}
    }
    if( isConv ) match = false;
  }

  return match ;
}
 

