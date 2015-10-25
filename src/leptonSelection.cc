#include "interface/xAna.hh"
#include "interface/MinitreeOut.hh"

#include <map>
#include <vector>
using namespace std;


//// 2012 muon selection

vector<int> xAna::selectMuonsHCP2012( int &vtxLep , double ptcut ) {
  // JM default ptcut is 20, ptcut added for VH dilep
  
  map<float,int,greater<float> > goodMuonIndex;
  int vtxInd = 0;
  
  for (int i = 0; i < nMu; ++i) {
    
    
    if( fabs((*muEta)[i]) > 2.4 ) continue;
    if( (*muPt)[i] < ptcut        ) continue;
    
    //if( !muID[i][1]        ) continue; 
    if( (*muType)[i]!=46      ) continue; 
    if( ((*muChi2NDF))[i] >= 10 ) continue; 
    if((*muNumberOfValidMuonHits)[i]<1)continue; 
    //    if(muChambers[i]<2 )continue;
    if( (*muStations)[i] < 2 ) continue; 
    /// find vertex:
    vtxInd = 0;
    double dzmin = 9999;
    for( int ivtx = 0 ; ivtx < nVtxBS; ivtx++ ) {
      if(  fabs( (*muVtx_z)[i]- (*vtxbs_z)[ivtx] ) < dzmin ) {
	dzmin = fabs( (*muVtx_z)[i] - (*vtxbs_z)[ivtx] );
	vtxInd = ivtx;
      }
    }
    if(fabs((*muD0Vtx)[i][vtxInd])>0.2)continue;    
    if(fabs((*muDzVtx)[i][vtxInd])>0.5)continue;    
    if((*muNumberOfValidPixelHits)[i] <=0 )continue;
    if((*muNumberOfValidTrkLayers)[i] <=5 )continue;  
    double muPFIsoSumRel =  ( (*muPFIsoR04_CH)[i]+   max(0., (*muPFIsoR04_NH)[i]+(*muPFIsoR04_Pho)[i]
							 - 0.5 *  ((*muPFIsoR04_PU)[i])))/(*muPt)[i];   // deltaBeta correction   
    if( muPFIsoSumRel>0.2  )continue;
							 
    /// order muon in decreasing pT
    goodMuonIndex.insert( make_pair( (*muPt)[i], i ) );
  }
  vector<int> out; 
  map<float,int>::iterator p;
  for( p = goodMuonIndex.begin(); p != goodMuonIndex.end(); ++p ) out.push_back( p->second );

  /// find the lepton vertex
  vtxLep = -1;
  if( goodMuonIndex.size() > 0 ) {
    int iBestMu = goodMuonIndex.begin()->second;
    double dzmin = 9999;
    for( int ivtx = 0 ; ivtx < nVtxBS; ivtx++ )
     if(  fabs( (*muVtx_z)[iBestMu] - (*vtxbs_z)[ivtx] ) < dzmin ) {
	dzmin = fabs( (*muVtx_z)[iBestMu] - (*vtxbs_z)[ivtx] );
	vtxLep = ivtx;
      }
    
    if( vtxLep == -1 ) vtxLep = 0; 

   // CF for synch: ------------------------------
    _minitree->mu_pt[0]  = (*muPt)[0];
    _minitree->mu_eta[0] = (*muEta)[0];
    _minitree->mu_phi[0] = (*muPhi)[0];
    if( goodMuonIndex.size() > 1 ) {
      _minitree->mu_pt[1]  = (*muPt)[1];
      _minitree->mu_eta[1] = (*muEta)[1];
      _minitree->mu_phi[1] = (*muPhi)[1];
    }
    //---------------------------------------------

  }

  return out;
}


vector<int> xAna::selectElectronsHCP2012( int &vtxLep , double ptcut ) { 
  // JM default ptcut is 20, ptcut added for VH dilep
 
  map<float,int,greater<float> > goodElectronIndex;
  int vtxInd = 0;

  vector<int> vtxLepIndex; vtxLepIndex.resize(nEle,-1);
  for( int i = 0; i< nEle; ++i) {
    float fSCEta = fabs((*eleSCEta)[i]);
    
    if( (*elePt)[i] < ptcut ) continue;
    if( fSCEta > 2.5  ) continue;
    if( fSCEta > 1.442 && fSCEta < 1.566 ) continue;

    float EffectiveArea =0;
    if      (fSCEta >= 0.000 && fSCEta < 1.000 ) EffectiveArea = 0.013 + 0.122;
    else if (fSCEta >= 1.000 && fSCEta < 1.479 ) EffectiveArea = 0.021 + 0.147;
    else if (fSCEta >= 1.479 && fSCEta < 2.000 ) EffectiveArea = 0.013 + 0.055;
    else if (fSCEta >= 2.000 && fSCEta < 2.200 ) EffectiveArea = 0.010 + 0.106;
    else if (fSCEta >= 2.200 && fSCEta < 2.300 ) EffectiveArea = 0.024 + 0.138;
    else if (fSCEta >= 2.300 && fSCEta < 2.400 ) EffectiveArea = 0.020 + 0.221;
    else                                         EffectiveArea = 0.019 + 0.211;

    float elePFIsoSumRel = ((*elePFChIso03)[i] + TMath::Max(float(0), (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i]
							 - EffectiveArea*rho2012))/(*elePt)[i];

    /// fill vertex lepton index for each lepton
    float dzmin = 9999;
    for( int ivtx = 0 ; ivtx < nVtxBS; ivtx++ )
      if(  fabs( (*eleVtx_z)[i] - (*vtxbs_z)[ivtx] ) < dzmin ) {
	dzmin = fabs( (*eleVtx_z)[i] - (*vtxbs_z)[ivtx] );
	vtxLepIndex[i] = ivtx;
      }
    if( vtxLepIndex[i] == -1 ) vtxLepIndex[i] = 0;
	
    vtxInd = vtxLepIndex[i];
    if( elePFIsoSumRel > 0.15 ) continue;
    if( (*eleMissHits)[i] >    1 ) continue;
    if( fabs((*eleD0Vtx)[i][vtxInd] ) > 0.02 ) continue;
    if( fabs((*eleDzVtx)[i][vtxInd] ) > 0.20 ) continue;
    if( (*eleIDMVANonTrig)[i] < 0.90 ) continue;

    if( (*eleConvVtxFit)[i] == 1 ) continue;

    goodElectronIndex.insert( make_pair( (*eleIDMVANonTrig)[i], i ) );
  }
  
  vector<int> out; 
  map<float,int>::iterator p;
  if( goodElectronIndex.size() > 0 ) {
    for( p = goodElectronIndex.begin(); p != goodElectronIndex.end(); ++p )  out.push_back( p->second );
    vtxLep = vtxLepIndex[p->second];
  }

  // CF for synch: ------------------------------
  if( out.size() > 0 ) {

    _minitree->ele_pt[0]  = (*elePt)[out[0]];
    _minitree->ele_eta[0] = (*eleEta)[out[0]];
    _minitree->ele_phi[0] = (*elePhi)[out[0]];
    if( out.size() > 1 ) {
      _minitree->ele_pt[1]  = (*elePt)[out[1]];
      _minitree->ele_eta[1] = (*eleEta)[out[1]];
      _minitree->ele_phi[1] = (*elePhi)[out[1]];
    }
  }   
  //--------------------------------------------- 

  return out;
}

/////////////////////////////////////// ICHEP 2012 selection /////////

vector<int> xAna::selectElectrons(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd) {

  vector<int> goodElectronIndex;
  //             Veto   Loose   Medium  Tight 
  //  dEtaIn    0.007   0.007   0.004   0.004     
  //  dPhiIn    0.8     0.15    0.06    0.03    
  //  sigmaIEtaIEta     0.01    0.01    0.01    0.01 
  //  H/E       0.15    0.12    0.12    0.12   
  // d0 (vtx)   0.04    0.02    0.02    0.02   
  // dZ (vtx)   0.2     0.2     0.1     0.1   
  // fabs(1/E - 1/p)    N/A     0.05    0.05    0.05  
  // PF isolation / pT  0.15    0.15    0.15    0.10  
  // Conversion rejection: vertex fit probability       N/A     1e-6    1e-6    1e-6 
  // Conversion rejection: missing hits         N/A     1       1       0                  

  Float_t dEtaInCut_B[4]   = {0.007, 0.007, 0.004, 0.004};
  Float_t dPhiInCut_B[4]   = {0.8,   0.15,  0.06,  0.03};
  Float_t sIeIeCut_B[4]    = {0.01,  0.01,  0.01,  0.01};
  Float_t HoECut_B[4]      = {0.15,  0.12,  0.12,  0.12};
  Float_t d0Cut_B[4]       = {0.04,  0.02,  0.02,  0.02};
  Float_t dzCut_B[4]       = {0.2,   0.2,   0.1,   0.1};
  Float_t diffEpCut_B[4]   = {999.,  0.05,  0.05,  0.05};
  Float_t relPFisoCut_B[4] = {0.15,  0.15,  0.15,  0.10};
  Int_t VtxFitProbCut_B[4] = {999,   1,     1,     1};
  Int_t missHitsCut_B[4]   = {99,    1,     1,     0};



  /// pT > 20 (pT < 20) 	Veto 	Loose 	Medium 	Tight
  // dEtaIn 	0.01 	0.009 	0.007 	0.005
  // dPhiIn 	0.7 	0.10 	0.03 	0.02
  // sigmaIEtaIEta 	0.03 	0.03 	0.03 	0.03
  // H/E 	N/A 	0.10 	0.10 	0.10
  // d0 (vtx) 	0.04 	0.02 	0.02 	0.02
  // dZ (vtx) 	0.2 	0.2 	0.1 	0.1
  // fabs(1/E - 1/p) 	N/A 	0.05 	0.05 	0.05
  // PF isolation / pT 	0.15 	0.15(0.10) 	0.15(0.10) 	0.10(0.07)
  // Conversion rejection: vertex fit probability 	N/A 	1e-6 	1e-6 	1e-6
  // Conversion rejection: missing hits 	N/A 	1 	1 	0 

  Float_t dEtaInCut_EC[4]   = {0.01 , 0.009,  0.007,  0.005};
  Float_t dPhiInCut_EC[4]   = {0.7,   0.10,   0.03,   0.02};
  Float_t sIeIeCut_EC[4]    = {0.03,  0.03,   0.03,   0.03};
  Float_t HoECut_EC[4]      = {999,   0.10,   0.10,   0.10};
  Float_t d0Cut_EC[4]       = {0.04,  0.02,   0.02,   0.02};
  Float_t dzCut_EC[4]       = {0.2,   0.2,    0.1,    0.1};
  Float_t diffEpCut_EC[4]   = {999,   0.05,   0.05,   0.05};
  Float_t relPFisoCut_EC[4] = {0.15,  0.15,   0.15,   0.10};
  Int_t VtxFitProbCut_EC[4] = {999,   1,      1,      1};
  Int_t missHitsCut_EC[4]   = {99,    1,      1,      0};

  Int_t nGoodEle = 0;


  for (int i=0; i<nEle; ++i) {
    double elePFIsoSumRel = 0.;

    double eleIsoRho = max(double(0),double(rho2012));

    float EffectiveArea =0;
    float fabsSCEta = fabs((*eleSCEta)[i]);

    // deltaR=0.3, 2011
   //  if (fabsSCEta >= 0.0 && fabsSCEta < 1.0) EffectiveArea = 0.100;
//     if (fabsSCEta >= 1.0 && fabsSCEta < 1.479) EffectiveArea = 0.120;
//     if (fabsSCEta >= 1.479 && fabsSCEta < 2.0) EffectiveArea = 0.085;
//     if (fabsSCEta >= 2.0 && fabsSCEta < 2.2) EffectiveArea = 0.110;
//     if (fabsSCEta >= 2.2 && fabsSCEta < 2.3) EffectiveArea = 0.120;
//     if (fabsSCEta >= 2.3 && fabsSCEta < 2.4) EffectiveArea = 0.120;
//     if (fabsSCEta >= 2.4) EffectiveArea = 0.130;

    if (fabsSCEta >= 0.0   && fabsSCEta < 1.000 ) EffectiveArea = 0.013 + 0.122;
    if (fabsSCEta >= 1.0   && fabsSCEta < 1.479 ) EffectiveArea = 0.021 + 0.147;
    if (fabsSCEta >= 1.479 && fabsSCEta < 2.000 ) EffectiveArea = 0.013 + 0.055;
    if (fabsSCEta >= 2.0   && fabsSCEta < 2.200 ) EffectiveArea = 0.010 + 0.106;
    if (fabsSCEta >= 2.2   && fabsSCEta < 2.300 ) EffectiveArea = 0.024 + 0.138;
    if (fabsSCEta >= 2.3   && fabsSCEta < 2.400 ) EffectiveArea = 0.020 + 0.221;
    if (fabsSCEta >= 2.4 ) EffectiveArea = 0.019 + 0.211;

    elePFIsoSumRel = ((*elePFChIso03)[i] + max(0., (*elePFPhoIso03)[i] + (*elePFNeuIso03)[i] - EffectiveArea*eleIsoRho))/(*elePt)[i];

    if ((*elePt)[i] < 20) continue; 
    if (fabsSCEta > 2.5) continue;
    if (fabsSCEta > 1.4442 && fabs((*eleSCEta)[i]) < 1.566) continue;

    int idWP = 1; // Loose eID

    if (fabs((*eleSCEta)[i]) < 1.4442) {
      helePFRelIso_BID_B->Fill(elePFIsoSumRel,wei);
      if (fabs((*eledEtaAtVtx)[i]) > dEtaInCut_B[idWP]) continue;
      if (fabs((*eledPhiAtVtx)[i]) > dPhiInCut_B[idWP]) continue;
      if ((*eleSigmaIEtaIEta)[i] > sIeIeCut_B[idWP]) continue;
      if ((*eleHoverE)[i] > HoECut_B[idWP]) continue;
      if (fabs((*eleD0Vtx)[i][vtxInd]) > d0Cut_B[idWP]) continue;
      if (fabs((*eleDzVtx)[i][vtxInd]) > dzCut_B[idWP]) continue;
      if (fabs(1./(*eleEcalEn)[i] - 1./(*elePin)[i]) > diffEpCut_B[idWP]) continue;
      if ((*eleConvVtxFit)[i] == VtxFitProbCut_B[idWP]) continue;
      if ((*eleMissHits)[i] > missHitsCut_B[idWP]) continue;

      helePFRelIso_AID_B->Fill(elePFIsoSumRel,wei);
      if (elePFIsoSumRel > relPFisoCut_B[idWP]) continue;
    } else {
      helePFRelIso_BID_EC->Fill(elePFIsoSumRel,wei);
      if (fabs((*eledEtaAtVtx)[i]) > dEtaInCut_EC[idWP]) continue;
      if (fabs((*eledPhiAtVtx)[i]) > dPhiInCut_EC[idWP]) continue;
      if ((*eleSigmaIEtaIEta)[i] > sIeIeCut_EC[idWP]) continue;
      if ((*eleHoverE)[i] > HoECut_EC[idWP]) continue;
      if (fabs((*eleD0Vtx)[i][vtxInd]) >d0Cut_EC[idWP]) continue;
      if (fabs((*eleDzVtx)[i][vtxInd]) >dzCut_EC[idWP]) continue;
      if (fabs(1./(*eleEcalEn)[i] - 1./(*elePin)[i]) > diffEpCut_EC[idWP]) continue;
      if ((*eleConvVtxFit)[i] == VtxFitProbCut_EC[idWP]) continue;
      if ((*eleMissHits)[i] > missHitsCut_EC[idWP]) continue;
      helePFRelIso_AID_EC->Fill(elePFIsoSumRel,wei);
      if (elePFIsoSumRel > relPFisoCut_EC[idWP]) continue;
    }

    nGoodEle++;

    double dR1_FSR  = deltaR((*phoEtaVtx)[iPho1][vtxInd], (*phoPhiVtx)[iPho1][vtxInd], (*eleEtaVtx)[i][vtxInd], (*elePhiVtx)[i][vtxInd]);
    double dR2_FSR  = deltaR((*phoEtaVtx)[iPho2][vtxInd], (*phoPhiVtx)[iPho2][vtxInd], (*eleEtaVtx)[i][vtxInd], (*elePhiVtx)[i][vtxInd]);
    double dRminFSR = min(dR1_FSR, dR2_FSR);

    double dR1_SC = deltaR((*phoSCEta)[iPho1], (*phoSCPhi)[iPho1], (*eleSCEta)[i], (*eleSCPhi)[i]);
    double dR2_SC = deltaR((*phoSCEta)[iPho2], (*phoSCPhi)[iPho2], (*eleSCEta)[i], (*eleSCPhi)[i]);
    double dRminSC = min(dR1_SC,dR2_SC);


    hDR_Pho1Ele->Fill(dR1_FSR,wei);     //  
    hDR_Pho2Ele->Fill(dR2_FSR,wei);

    hDR_SCPho1Ele->Fill(dR1_SC,wei);     //  
    hDR_SCPho2Ele->Fill(dR2_SC,wei);

    hDRmin_PhoEle->Fill(dRminFSR,wei);
    hDRmin_SCPhoEle->Fill(dRminSC,wei);

    if( dRminFSR > 1.0 )  goodElectronIndex.push_back( i );;
  }

  return goodElectronIndex;
}


vector<int> xAna::selectMuons(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd) {
  vector<int> goodMuonIndex;

  for (int i=0; i<nMu; ++i) {
    hAllMuPt->Fill((*muPt)[i],wei);
    hAllMuEta->Fill((*muEta)[i],wei);
    if( fabs((*muEta)[i]) > 2.4 ) continue;
    if( (*muPt)[i] < 20         ) continue;

    double muPFIsoSumRel =  ( (*muPFIsoR04_CH)[i]+   max(0., (*muPFIsoR04_NH)[i]+(*muPFIsoR04_Pho)[i]- 0.5 *  (*muPFIsoR04_PU)[i]))/(*muPt)[i];   // deltaBeta correction
    hmuPFRelIso_BID->Fill(muPFIsoSumRel,wei);
    //if(!muID[i][1])continue;          
    if((*muType)[i]!=46)continue;        
    if( (*muChi2NDF)[i] >= 10)continue;  
    if((*muNumberOfValidMuonHits)[i]<1)continue;
    //    if(muChambers[i]<2 )continue;
    if( (*muStations)[i] < 2 ) continue; 
    if(fabs((*muD0Vtx)[i][vtxInd])>0.2)continue;
    if(fabs((*muDzVtx)[i][vtxInd])>0.5)continue;
    if((*muNumberOfValidPixelHits)[i]<1)continue;
    if((*muNumberOfValidTrkLayers)[i]<6)continue;
    hmuPFRelIso_AID->Fill(muPFIsoSumRel,wei);
    if(muPFIsoSumRel>0.2)continue;
  
    double dR1_FSR = deltaR((*phoEta)[iPho1], (*phoPhi)[iPho1], (*muEta)[i], (*muPhi)[i]);
    double dR2_FSR = deltaR((*phoEta)[iPho2], (*phoPhi)[iPho2], (*muEta)[i], (*muPhi)[i]);
    double dRminFSR = min(dR1_FSR,dR2_FSR);
    hDRmin_PhoMu->Fill(dRminFSR,wei);
    //  hIdMuPt->Fill((*muPt)[i],wei);
    if( dRminFSR > 1.0 ) { goodMuonIndex.push_back( i ); }
  }
  return goodMuonIndex;
}




// JM ***************************************************************************************
 // WARNING: Globe does it that way! This might be stupid however

vector<int> xAna::selectElectronsAwayFromTwoPhotons(int i, int j, vector<int> elecIndex ){
  
  vector<int> selElecIndex;
  int ilead = i, itrail = j;

  for(unsigned ie=0;ie<elecIndex.size();ie++){
    
    int iel=elecIndex[ie];
    
    double dR1 = deltaR((*phoEta)[ilead], (*phoPhi)[ilead], (*eleEta)[iel], (*elePhi)[iel]);
    double dR2 = deltaR((*phoEta)[itrail], (*phoPhi)[itrail], (*eleEta)[iel], (*elePhi)[iel]);
    double dRmin = min(dR1,dR2);

    double dRTrk1=(*phoCiCdRtoTrk)[ilead];
    double dRTrk2=(*phoCiCdRtoTrk)[itrail];
    double dRminTrk = min(dRTrk1,dRTrk2);

    
    //cout<<" JM DEBUG ========= dR1="<<dR1<<" dR2="<< dR2<< endl;

    if( dRmin < 1.0 ) continue;
    if( dRminTrk < 1.0 ) continue;

    selElecIndex.push_back(iel);
  }
  return selElecIndex;  
}





Int_t* xAna::selectMuons2011(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd) {


  Int_t *results = new Int_t[5];
  for (int i=0; i<4; ++i) results[i+1] = -99;
  Int_t nGoodMu = 0;
  // cout <<" puwei  "<<  puwei << endl;

  for (int i=0; i<nMu; ++i) {
    hAllMuPt->Fill((*muPt)[i],wei);
    hAllMuEta->Fill((*muEta)[i],wei);
    if ((*muPt)[i] < 20) continue;
    
    if (fabs((*muEta)[i]) > 2.4 ) continue;
    double combRelIsoMu = ((*muIsoTrk)[i]+(*muIsoEcal)[i]+(*muIsoHcal)[i]- rho25 * 3.14 * 0.09)/(*muPt)[i] ;

    hMuCombRelIso->Fill(combRelIsoMu,wei);
    
    if ((*muType)[i] >= 6 && // GlobalMuon && TrackerMuon 
	(*muChi2NDF)[i] < 10 &&
	(*muNumberOfValidPixelHits)[i] > 0 &&
	(*muNumberOfValidTrkHits)[i] > 10  &&
	(*muNumberOfValidMuonHits)[i] > 0  &&
	(*muStations)[i] > 1 && 
	//	fabs(muD0[i]) < 0.02 &&
	//	fabs(muDz[i]) < 0.1 &&
	fabs((*muD0Vtx)[i][vtxInd]) < 0.02 &&
	fabs((*muDzVtx)[i][vtxInd]) < 0.1 )
	//	&&muID[i][1] == 1) 
	{
      
      if (((*muIsoTrk)[i]+(*muIsoEcal)[i]+(*muIsoHcal)[i]- rho25 * 3.14 * 0.09 ) < 0.1 * (*muPt)[i]) {
		
	//	double dR1_FSR = deltaR((*(*phoEta)Vtx)[iPho1][vtxInd], (*phoPhi)[iPho1], (*eleEtaVtx)[i][vtxInd], (*elePhiVtx)[i][vtxInd]);
	//	double dR2_FSR = deltaR((*(*phoEta)Vtx)[iPho2][vtxInd], (*phoPhi)[iPho2], (*eleEtaVtx)[i][vtxInd], (*elePhiVtx)[i][vtxInd]);
	double dR1_FSR = deltaR((*phoEta)[iPho1], (*phoPhi)[iPho1], (*muEta)[i], (*muPhi)[i]);
	double dR2_FSR = deltaR((*phoEta)[iPho2], (*phoPhi)[iPho2], (*muEta)[i], (*muPhi)[i]);
	double dRminFSR = min(dR1_FSR,dR2_FSR);
	hDRmin_PhoMu->Fill(dRminFSR,wei);
	
	if(dRminFSR>1){  	
	  results[nGoodMu+1] = i;
	  nGoodMu++;
	}
      }
    }
  }
  results[0] = nGoodMu;
  //  cout <<" results[0]  "<< results[0]   <<endl;
  return results;
 }




Int_t* xAna::selectElectrons2011(Int_t iPho1, Int_t iPho2, Double_t wei, Int_t vtxInd) {

  Int_t *results = new Int_t[5];
  for (int i=0; i<4; ++i) results[i+1] = -99;
  Int_t nFinalEle = 0;
  Int_t nGoodEle = 0;

  // cout <<" nele "<<  nEle <<endl;
  for (int i=0; i<nEle; ++i) {
    // cout << i  <<"  ele pt  "<< (*elePt)[i] <<endl;
    
    hAllElePt->Fill((*elePt)[i],wei);
    hAllSCEleEta->Fill((*eleSCEta)[i],wei);
    // hDEtaAllEle->Fill((*eledEtaAtVtx)[i],wei);
    // hDPhiAllEle->Fill((*eledPhiAtVtx)[i],wei);

    if ((*elePt)[i] < 20) continue;   // pt>35 ?
    if (fabs((*eleSCEta)[i]) > 2.5) continue;
    if (fabs((*eleSCEta)[i]) > 1.4442 && fabs((*eleSCEta)[i]) < 1.566) continue;
    
    //spike cleaning in EB
    if (fabs((*eleSCEta)[i]) < 1.4442) {    // is it still needed ?????????????
      if ((*eleSigmaIEtaIEta)[i]<0.002) continue;
    } 
   
    if(fabs((*eleSCEta)[i]) < 1.4442){
      float eleCombIso = ( (*eleIsoTrkDR03)[i] + max(0., (*eleIsoEcalDR03)[i] - 1.) +  (*eleIsoHcalDR03)[i]  - rho25 * TMath::Pi() *0.3*0.3) / (*elePt)[i]; 
      hEleCombRelIso_EB->Fill(eleCombIso,wei);
      //  cout <<" (*eledPhiAtVtx)[i] "<<(*eledPhiAtVtx)[i]   <<" (*eledEtaAtVtx)[i] "<<(*eledEtaAtVtx)[i]   <<     endl;
      if(eleCombIso> 0.053) continue;
      if((*eleSigmaIEtaIEta)[i]>  0.01) continue;
      if(fabs((*eledPhiAtVtx)[i])> 0.039 )continue;
      if(fabs((*eledEtaAtVtx)[i])> 0.005 )continue;
      //if((*eleHoverE)[i]> 0.04  )continue;
    }else{
      float eleCombIso = ( (*eleIsoTrkDR03)[i] + (*eleIsoEcalDR03)[i]  + (*eleIsoHcalDR03)[i] - rho25 * TMath::Pi() *0.3*0.3) / (*elePt)[i];
      hEleCombRelIso_EE->Fill(eleCombIso,wei);

      // cout <<" (*eledPhiAtVtx)[i] "<<(*eledPhiAtVtx)[i]   <<" (*eledEtaAtVtx)[i] "<<(*eledEtaAtVtx)[i]   <<     endl;
      if(eleCombIso> 0.042) continue;
      if((*eleSigmaIEtaIEta)[i]>  0.03) continue;
      if(fabs((*eledPhiAtVtx)[i])> 0.028 )continue;
      if(fabs((*eledEtaAtVtx)[i])> 0.007 )continue;
      //   if((*eleHoverE)[i]> 0.15  )continue;
    }
    if ((*eleMissHits)[i] != 0) continue;
    if (fabs((*eleConvDist)[i]) < 0.02 && fabs((*eleConvDcot)[i]) < 0.02) continue;
    //if(fabs(eleD0[i]) >= 0.02)continue;
    //if(fabs(eleDz[i]) >= 0.1)continue; 

    if(fabs((*eleD0Vtx)[i][vtxInd]) >= 0.02)continue;
    if(fabs((*eleDzVtx)[i][vtxInd]) >= 0.1)continue; 
    nGoodEle++;
    
    double dR1_SC = deltaR((*phoSCEta)[iPho1], (*phoSCPhi)[iPho1], (*eleSCEta)[i], (*eleSCPhi)[i]);
    double dR2_SC = deltaR((*phoSCEta)[iPho2], (*phoSCPhi)[iPho2], (*eleSCEta)[i], (*eleSCPhi)[i]);

    double dR1_FSR = deltaR((*phoEtaVtx)[iPho1][vtxInd], (*phoPhiVtx)[iPho1][vtxInd], (*eleEtaVtx)[i][vtxInd], (*elePhiVtx)[i][vtxInd]);
    double dR2_FSR = deltaR((*phoEtaVtx)[iPho2][vtxInd], (*phoPhiVtx)[iPho2][vtxInd], (*eleEtaVtx)[i][vtxInd], (*elePhiVtx)[i][vtxInd]);

    // cout <<" ele far from the photons? "<<  i  <<" dR1 "   << dR1 << " dR2 "   << dR2 <<  endl;
    hDR_Pho1Ele->Fill(dR1_FSR,wei);     //  
    hDR_Pho2Ele->Fill(dR2_FSR,wei);

    hDR_SCPho1Ele->Fill(dR1_SC,wei);     //  
    hDR_SCPho2Ele->Fill(dR2_SC,wei);

    double dRminFSR = min(dR1_FSR,dR2_FSR);
    double dRminSC = min(dR1_SC,dR2_SC);

    // cout  << "  dR1 = "  << dR1   << "  dR2 = "  << dR2   <<"   --> dRmin = "<<  dRmin  <<endl;
    
    hDRmin_PhoEle->Fill(dRminFSR,wei);
    hDRmin_SCPhoEle->Fill(dRminSC,wei);
    //   if (dR1 > 0.3 && dR2 > 0.3){

    //  if(dRmin>0.5){
    if(dRminFSR>1){    //changed the 7th October
      
      results[nFinalEle+1] = i;      
      nFinalEle++;
    }   
  }
  results[0] = nFinalEle;
  //  if(results[0]>0)cout <<" ele results[1] "<<    results[1]   <<endl;
  return results;  
}

