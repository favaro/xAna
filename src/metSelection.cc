#include "interface/xAna.hh"
#include "interface/MinitreeOut.hh"

#include "JetMETObjects/FactorizedJetCorrector.h"
#include "JetMETObjects/JetCorrectionUncertainty.h"
#include "JetMETObjects/JetCorrectorParameters.h"
#include "JetMETObjects/JetResolution.h"
#include "JetMETObjects/SimpleJetCorrectionUncertainty.h"
#include "JetMETObjects/SimpleJetCorrector.h"


TLorentzVector xAna::corMET( int ilead, int itrail ) {
  
   _minitree->mtree_rawMet    = recoPfMET;
   _minitree->mtree_rawMetPhi = recoPfMETPhi;

  //3.1 Soft Jet correction (FC?)
  applyJEC4SoftJets();

   //3.2 smearing
  if( !isData ) METSmearCorrection(ilead, itrail);

  //3.3 shifting (even in data? but different in data and MC)
  METPhiShiftCorrection();

  //3.4 scaling
  if( isData) METScaleCorrection(ilead, itrail);

  cout << "   met 4 " << endl;
  
  _minitree->mtree_corMet    = recoPfMET; // FIXME: move this to fillMETToMinitree?
  _minitree->mtree_corMetPhi = recoPfMETPhi;

  TLorentzVector met; 
  met.SetPtEtaPhiM( _minitree->mtree_corMet, 0, _minitree->mtree_corMetPhi = recoPfMETPhi, 0);  

  return met;
}

bool xAna::MetTagSelection( const TLorentzVector &g1, const TLorentzVector &g2, const vector<int>& jind ) {
  TLorentzVector gg = g1+g2;
  
  if( g1.Pt()/gg.M() < 45./120. ) return false;
  if( g2.Pt() < 25. ) return false;   
  
  TLorentzVector met; 
  met.SetPtEtaPhiM( _minitree->mtree_corMet, 0, _minitree->mtree_corMetPhi = recoPfMETPhi, 0);  

  //-- met back-to-bakc with gg
  if( fabs( gg.DeltaPhi( met ) ) <= 2.1 ) return false;

  //-- no jet bak-to-back with the met
  TLorentzVector jet;
  float jetThreshold = 50;
  float jetHighestPt = -1;
  int lJet = -1;
  for( unsigned ijet = 0; ijet < jind.size(); ++ijet ) {
    int j_ind = jind[ijet];
    if( (*jetPt)[j_ind] > jetThreshold ) {
     jet.SetPtEtaPhiE( (*jetPt)[j_ind], (*jetEta)[j_ind], (*jetPhi)[j_ind], (*jetEn)[j_ind] );
     if( jet.DeltaR( g1 ) < 0.5 || jet.DeltaR( g2 ) < 0.5 ) continue;
     if( jet.Pt() > jetHighestPt ) { jetHighestPt = jet.Pt(); lJet = j_ind; }
    }
  }

  if( jetHighestPt > jetThreshold ) {
    jet.SetPtEtaPhiE( (*jetPt)[lJet], (*jetEta)[lJet], (*jetPhi)[lJet], (*jetEn)[lJet] );
    if( fabs(jet.DeltaPhi(met)) > 2.7 ) return false;
  }

  return true;
}

void xAna::applyJEC4SoftJets() {
  setupJEC();

  for( int i = 0 ; i < nLowPtJet; ++i ) {
    /// step 4
    JetCorrector->setJetEta((*jetLowPtEta)[i]);
    JetCorrector->setJetPt((*jetLowPtRawPt)[i]);
    JetCorrector->setJetA((*jetLowPtArea)[i]);
    JetCorrector->setRho(rho2012); 
    /// step 5
    double correction = JetCorrector->getCorrection();
    // vector<double> factors = JetCorrector->getSubCorrections();

    /// step 6
    (*jetLowPtPt)[i] = (*jetLowPtRawPt)[i] * correction;
    (*jetLowPtEn)[i] = (*jetLowPtRawEn)[i] * correction;
    /// --------------------------------------------------------------------------------------//
  }
}


// MET Correction Step 1: MC smearing
void xAna::METSmearCorrection(Int_t iPho1, Int_t iPho2) {
  TRandom3 genRan((UInt_t)event);

  TLorentzVector jetSumSmeared; 
  jetSumSmeared.SetXYZT(0.,0.,0.,0.);
  TLorentzVector jetSumRaw;
  jetSumRaw.SetXYZT(0.,0.,0.,0.);
    
  // jets of 5<rawPt<15
  for (int j_ind=0; j_ind<nLowPtJet; j_ind++) {
    // jet basic cuts
    if ((*jetLowPtRawPt)[j_ind]<10. || fabs((*jetLowPtEta)[j_ind])>4.7) continue;
    
    // remove selected photons
    float dR_jetg1 = deltaR((*phoEta)[iPho1], (*phoPhi)[iPho1],(*jetLowPtEta)[j_ind],(*jetLowPtPhi)[j_ind]);
    float dR_jetg2 = deltaR((*phoEta)[iPho2], (*phoPhi)[iPho2],(*jetLowPtEta)[j_ind],(*jetLowPtPhi)[j_ind]);
    if( dR_jetg1 < 0.5 || dR_jetg2 < 0.5 ) continue;
 
    double SmearC = getMETJetSmearFactor((*jetLowPtEta)[j_ind]);
    Double_t minDRgenJet(999.);
    int      genMatched(-999);
    int      typeMatched(0); // 0: no matched, 1: LowPt jets, 2: normal jets


    cout << "    submet nlowptjet " << nLowPtJet << endl;
 
    // Search for gen jet, soft jets
    for (int iGenJet=0; iGenJet<nLowPtJet; iGenJet++) {

      cout << "    submet 2A " << jetLowPtGenJetPt << endl;

      if ((*jetLowPtGenJetPt)[iGenJet]<0.) continue;

      cout << "    submet 2B " << endl;
      Double_t dR_jg  = deltaR((*jetLowPtGenJetEta)[iGenJet],(*jetLowPtGenJetPhi)[iGenJet],(*jetLowPtEta)[j_ind],(*jetLowPtPhi)[j_ind]);

      cout << "    submet 2C " << endl;
      double   expres = ErrEt((*jetLowPtPt)[j_ind],(*jetLowPtEta)[j_ind]);

      cout << "    submet 3" << endl;

      if(dR_jg < minDRgenJet && ((*jetLowPtPt)[j_ind]-(*jetLowPtGenJetPt)[iGenJet])/(*jetLowPtPt)[j_ind] < 5. * expres) 
{
        typeMatched = 1;
        genMatched  = iGenJet;
        minDRgenJet = dR_jg;
      }
    }
    // Search for gen jet
    for (int iGenJet=0; iGenJet<nJet; iGenJet++) {
      if ((*jetGenJetPt)[iGenJet]<0.) continue;

      Double_t dR_jg  = deltaR((*jetGenJetEta)[iGenJet],(*jetGenJetPhi)[iGenJet],(*jetLowPtEta)[j_ind],(*jetLowPtPhi)[j_ind]);
      double   expres = ErrEt((*jetLowPtPt)[j_ind],(*jetLowPtEta)[j_ind]);

      if(dR_jg < minDRgenJet && ((*jetLowPtPt)[j_ind]-(*jetGenJetPt)[iGenJet])/(*jetLowPtPt)[j_ind] < 5. * expres) {
        typeMatched = 2;
        genMatched  = iGenJet;
        minDRgenJet = dR_jg;
      }
    }

    if ((genMatched>-1) && (typeMatched==1) && (minDRgenJet > 0.1 + 0.3 * exp(-0.05*((*jetLowPtGenJetPt)[genMatched]-10)))) genMatched = -999;
    if ((genMatched>-1) && (typeMatched==2) && (minDRgenJet > 0.1 + 0.3 * exp(-0.05*((*jetGenJetPt)[genMatched]-10))))      
genMatched = -999;
        
    double SmearRes = 0.;
    if(genMatched>-1 && typeMatched>0) {
      if (     typeMatched==1) SmearRes = (SmearC-1.) * 
((*jetLowPtPt)[j_ind]-(*jetLowPtGenJetPt)[genMatched])/(*jetLowPtPt)[j_ind];
      else if (typeMatched==2) SmearRes = (SmearC-1.) * ((*jetLowPtPt)[j_ind]-(*jetGenJetPt)[genMatched])/(*jetLowPtPt)[j_ind];
    } else {
       double expres = ErrEt((*jetLowPtRawPt)[j_ind], (*jetLowPtEta)[j_ind]);
       double relsmear = expres * sqrt(SmearC*SmearC-1);
       genRan.SetSeed(event+(Int_t)(*jetLowPtEta)[j_ind]*1000);
       SmearRes = genRan.Gaus(0.,relsmear);
    }
    
    float rawPtSmeared = (*jetLowPtRawPt)[j_ind];
    float rawEnSmeared = (*jetLowPtRawEn)[j_ind];

    if(SmearRes>-1 && SmearRes < 2) {
      rawPtSmeared *= (1 + SmearRes);
      rawEnSmeared *= (1 + SmearRes);
    }
 
    TLorentzVector smearedJet;
    smearedJet.SetPtEtaPhiE(rawPtSmeared,(*jetLowPtEta)[j_ind],(*jetLowPtPhi)[j_ind],rawEnSmeared);
        
    TLorentzVector rawJet;
    rawJet.SetPtEtaPhiE((*jetLowPtRawPt)[j_ind],(*jetLowPtEta)[j_ind],(*jetLowPtPhi)[j_ind],(*jetLowPtRawEn)[j_ind]);
        
    jetSumSmeared += smearedJet;
    jetSumRaw     += rawJet;

  }

  // jets of rawPt>=15
  for (int j_ind=0; j_ind<nJet; j_ind++) {
    (*jetRawPtSmeared)[j_ind] = (*jetRawPt)[j_ind];
    (*jetRawEnSmeared)[j_ind] = (*jetRawEn)[j_ind];

    // jet basic cuts
    if ((*jetRawPt)[j_ind]<10. || fabs((*jetEta)[j_ind])>4.7) continue;

    // remove selected photons
    float dR_jetg1 = deltaR((*phoEta)[iPho1], (*phoPhi)[iPho1],(*jetEta)[j_ind],(*jetPhi)[j_ind]);
    float dR_jetg2 = deltaR((*phoEta)[iPho2], (*phoPhi)[iPho2],(*jetEta)[j_ind],(*jetPhi)[j_ind]);
    if( dR_jetg1 < 0.5 || dR_jetg2 < 0.5 ) continue;

    double SmearC = getMETJetSmearFactor((*jetEta)[j_ind]);
    Double_t minDRgenJet(999.);
    int      genMatched(-999);
    int      typeMatched(0);

    // Search for gen jet, soft jets
    for (int iGenJet=0; iGenJet<nLowPtJet; iGenJet++) {
      if ((*jetLowPtGenJetPt)[iGenJet]<0.) continue;

      Double_t dR_jg  = 
deltaR((*jetLowPtGenJetEta)[iGenJet],(*jetLowPtGenJetPhi)[iGenJet],(*jetEta)[j_ind],(*jetPhi)[j_ind]);
      double   expres = ErrEt((*jetPt)[j_ind],(*jetEta)[j_ind]);

      if(dR_jg < minDRgenJet && ((*jetPt)[j_ind]-(*jetLowPtGenJetPt)[iGenJet])/(*jetPt)[j_ind] < 5. * expres) {
        typeMatched = 1;
        genMatched  = iGenJet;
        minDRgenJet = dR_jg;
      }
    }
    // Search for gen jet
    for (int iGenJet=0; iGenJet<nJet; iGenJet++) {
      if ((*jetGenJetPt)[iGenJet]<0.) continue;
      
      Double_t dR_jg  = deltaR((*jetGenJetEta)[iGenJet],(*jetGenJetPhi)[iGenJet],(*jetEta)[j_ind],(*jetPhi)[j_ind]);
      double   expres = ErrEt((*jetPt)[j_ind],(*jetEta)[j_ind]);

      if(dR_jg < minDRgenJet && ((*jetPt)[j_ind]-(*jetGenJetPt)[iGenJet])/(*jetPt)[j_ind] < 5. * expres) {
        typeMatched = 2;
        genMatched  = iGenJet;
        minDRgenJet = dR_jg;
      }
    }

    if ((genMatched>-1) && (typeMatched==1) && (minDRgenJet > 0.1 + 0.3 * 
						exp(-0.05*((*jetLowPtGenJetPt)[genMatched]-10)))) genMatched = -999;
    if ((genMatched>-1) && (typeMatched==2) && (minDRgenJet > 0.1 + 0.3 * exp(-0.05*((*jetGenJetPt)[genMatched]-10))))      
      genMatched = -999;

    double SmearRes = 0.;
    if(genMatched>-1 && typeMatched>0) {
      if (     typeMatched==1) SmearRes = (SmearC-1.) * 
((*jetPt)[j_ind]-(*jetLowPtGenJetPt)[genMatched])/(*jetPt)[j_ind];
      else if (typeMatched==2) SmearRes = (SmearC-1.) * ((*jetPt)[j_ind]-(*jetGenJetPt)[genMatched])/(*jetPt)[j_ind];
    } else {
       double expres = ErrEt((*jetRawPt)[j_ind], (*jetEta)[j_ind]);
       double relsmear = expres * sqrt(SmearC*SmearC-1);
       genRan.SetSeed(event+(Int_t)(*jetEta)[j_ind]*1000);
       SmearRes = genRan.Gaus(0.,relsmear);
    }

    float rawPtSmeared = (*jetRawPt)[j_ind];
    float rawEnSmeared = (*jetRawEn)[j_ind];

    if(SmearRes>-1 && SmearRes < 2) {
      rawPtSmeared *= (1 + SmearRes);
      rawEnSmeared *= (1 + SmearRes);
    }
    (*jetRawPtSmeared)[j_ind] = rawPtSmeared;
    (*jetRawEnSmeared)[j_ind] = rawEnSmeared;

    TLorentzVector smearedJet;
    smearedJet.SetPtEtaPhiE(rawPtSmeared,(*jetEta)[j_ind],(*jetPhi)[j_ind],rawEnSmeared);

    TLorentzVector rawJet;
    rawJet.SetPtEtaPhiE((*jetRawPt)[j_ind],(*jetEta)[j_ind],(*jetPhi)[j_ind],(*jetRawEn)[j_ind]);

    jetSumSmeared += smearedJet;
    jetSumRaw     += rawJet;

  }

  double corrMetPx = recoPfMET*(TMath::Cos(recoPfMETPhi));
  double corrMetPy = recoPfMET*(TMath::Sin(recoPfMETPhi));
  TLorentzVector smearCorrectedMET;
  smearCorrectedMET.SetPxPyPzE(corrMetPx, corrMetPy, 0., sqrt(corrMetPx*corrMetPx+corrMetPy*corrMetPy));
  smearCorrectedMET = smearCorrectedMET + jetSumRaw - jetSumSmeared;
  recoPfMET    = smearCorrectedMET.Pt();
  recoPfMETPhi = smearCorrectedMET.Phi();

}

// MET Correction Step 1b: get MC smear factor
double xAna::getMETJetSmearFactor(Float_t metJetEta) {
  double jetSmearFactor = 0.;

  if      (fabs(metJetEta)<=1.1) jetSmearFactor = 1.06177;
  else if (fabs(metJetEta)<=1.7) jetSmearFactor = 1.08352;
  else if (fabs(metJetEta)<=2.3) jetSmearFactor = 1.02911;
  else                           jetSmearFactor = 1.15288;

  //std::cout << "--> metCorr: smear factor = " << jetSmearFactor << std::endl;
  return jetSmearFactor;
}

// MET Correction Step 1c: get MC expected resolution
double xAna::ErrEt( double Et, double Eta) {
  double InvPerr2;

  double N(0), S(0), C(0), m(0);
  if(fabs(Eta) < 0.5 ) {
    N = 3.96859;
    S = 0.18348;
    C = 0.;
    m = 0.62627;
  } else if( fabs(Eta) < 1. ) {
    N = 3.55226;
    S = 0.24026;
    C = 0.;
    m = 0.52571;
  } else if( fabs(Eta) < 1.5 ) {
    N = 4.54826;
    S = 0.22652;
    C = 0.;
    m = 0.58963;
  } else if( fabs(Eta) < 2. ) {
    N = 4.62622;
    S = 0.23664;
    C = 0.;
    m = 0.48738;
  } else if( fabs(Eta) < 3. ) {
    N = 2.53324;
    S = 0.34306;
    C = 0.;
    m = 0.28662;
  } else if( fabs(Eta) < 5. ) {
    N = 2.95397;
    S = 0.11619;
    C = 0.;
    m = 0.96086;
  }

  // this is the absolute resolution (squared), not sigma(pt)/pt
  // so have to multiply by pt^2, thats why m+1 instead of m-1
  InvPerr2 =  (N * fabs(N) ) + (S * S) * pow(Et, m+1) + (C * C) * Et * Et ;

  return sqrt(InvPerr2)/Et;
} 


// MET Correction Step 2: phi shifting
void xAna::METPhiShiftCorrection() {
  double corr2METx = 0.;
  double corr2METy = 0.;

 
  // OLD (Moriond)
  /*
  if(_weight_manager->getCrossSection()->mode() < 0){ // <----- DATA
    corr2METx = recoPfMET*(TMath::Cos(recoPfMETPhi))-0.006239*recoPfMETsumEt+0.662;
    corr2METy = recoPfMET*(TMath::Sin(recoPfMETPhi))+0.004613*recoPfMETsumEt-0.673;    
  } else { // <----------- MC
    corr2METx = recoPfMET*(TMath::Cos(recoPfMETPhi))+0.00135*recoPfMETsumEt-0.021;
    corr2METy = recoPfMET*(TMath::Sin(recoPfMETPhi))+0.00371*recoPfMETsumEt-0.826;
  }
  */


  // NEW (Legacy)
  if(_weight_manager->getCrossSection()->mode() < 0){ // <----- DATA
    corr2METx = recoPfMET*(TMath::Cos(recoPfMETPhi))- (-0.83015 + 0.0051172*recoPfMETsumEt);
    corr2METy = recoPfMET*(TMath::Sin(recoPfMETPhi))- (0.384595 - 0.00281327*recoPfMETsumEt);    
  } else { // <----------- MC
    corr2METx = recoPfMET*(TMath::Cos(recoPfMETPhi))- ( -0.0844029 + 0.000685386 * recoPfMETsumEt);
    corr2METy = recoPfMET*(TMath::Sin(recoPfMETPhi))- (  0.415907  -0.00395011 *recoPfMETsumEt);
  }


  TVector3 shiftMET(corr2METx,corr2METy,0.);

  recoPfMET    = shiftMET.Pt();
  recoPfMETPhi = shiftMET.Phi();

  //std::cout << "--> metCorr: Px = " << corr2METx << ", Py = " << corr2METy << "  " << shiftMET.Pt() << "  " << sqrt(corr2METx*corr2METx + corr2METy*corr2METy) << std::endl;
}


// MET Correction Step 3: Data scaling
void xAna::METScaleCorrection(Int_t iPho1, Int_t iPho2) {
  TLorentzVector jetSumScaled;  
  jetSumScaled.SetXYZT(0.,0.,0.,0.);
  TLorentzVector jetSumRaw;
  jetSumRaw.SetXYZT(0.,0.,0.,0.);

  // jets with 10<pt<20
  for (int j_ind=0; j_ind<nLowPtJet; j_ind++) {
    // jet basic cuts
    if ((*jetLowPtRawPt)[j_ind]<10. || fabs((*jetLowPtEta)[j_ind])>4.7) continue;

    // remove selected photons
    float dR_jetg1 = deltaR((*phoEta)[iPho1], (*phoPhi)[iPho1],(*jetLowPtEta)[j_ind],(*jetLowPtPhi)[j_ind]);
    float dR_jetg2 = deltaR((*phoEta)[iPho2], (*phoPhi)[iPho2],(*jetLowPtEta)[j_ind],(*jetLowPtPhi)[j_ind]);
    if( dR_jetg1 < 0.5 || dR_jetg2 < 0.5 ) continue;

    double ScaleFactor = getMETJetScaleFactor((*jetLowPtEta)[j_ind]);
    double rawPtScaled = (*jetLowPtRawPt)[j_ind] * ScaleFactor;
    double rawEnScaled = (*jetLowPtRawEn)[j_ind] * ScaleFactor;

    TLorentzVector scaledJet;
    scaledJet.SetPtEtaPhiE(rawPtScaled,(*jetLowPtEta)[j_ind],(*jetLowPtPhi)[j_ind],rawEnScaled);
        
    TLorentzVector rawJet;
    rawJet.SetPtEtaPhiE((*jetLowPtRawPt)[j_ind],(*jetLowPtEta)[j_ind],(*jetLowPtPhi)[j_ind],(*jetLowPtRawEn)[j_ind]);
        
    jetSumScaled += scaledJet;
    jetSumRaw    += rawJet;

  }

  // jets with pt>=20
  for (int j_ind=0; j_ind<nJet; j_ind++) {
    // jet basic cuts
    if ((*jetRawPt)[j_ind]<10. || fabs((*jetEta)[j_ind])>4.7) continue;

    // remove selected photons
    float dR_jetg1 = deltaR((*phoEta)[iPho1], (*phoPhi)[iPho1],(*jetEta)[j_ind],(*jetPhi)[j_ind]);
    float dR_jetg2 = deltaR((*phoEta)[iPho2], (*phoPhi)[iPho2],(*jetEta)[j_ind],(*jetPhi)[j_ind]);
    if( dR_jetg1 < 0.5 || dR_jetg2 < 0.5 ) continue;

    double ScaleFactor = getMETJetScaleFactor((*jetEta)[j_ind]);
    double rawPtScaled = (*jetRawPt)[j_ind] * ScaleFactor;
    double rawEnScaled =(*jetRawEn)[j_ind] * ScaleFactor;

    TLorentzVector scaledJet;
    scaledJet.SetPtEtaPhiE(rawPtScaled,(*jetEta)[j_ind],(*jetPhi)[j_ind],rawEnScaled);

    TLorentzVector rawJet;
    rawJet.SetPtEtaPhiE((*jetRawPt)[j_ind],(*jetEta)[j_ind],(*jetPhi)[j_ind],(*jetRawEn)[j_ind]);

    jetSumScaled += scaledJet;
    jetSumRaw    += rawJet;
  }

  double corrMetPx = recoPfMET*(TMath::Cos(recoPfMETPhi));
  double corrMetPy = recoPfMET*(TMath::Sin(recoPfMETPhi));
  TLorentzVector scaleCorrectedMET;
  scaleCorrectedMET.SetPxPyPzE(corrMetPx, corrMetPy, 0., sqrt(corrMetPx*corrMetPx+corrMetPy*corrMetPy));
  scaleCorrectedMET = scaleCorrectedMET + jetSumRaw - jetSumScaled;

  recoPfMET    = scaleCorrectedMET.Pt();
  recoPfMETPhi = scaleCorrectedMET.Phi();
}

// MET Correction Step 3b: get Data scale factor
double xAna::getMETJetScaleFactor(Float_t metJetEta) {
  double jetScaleFactor = 1.;

  if(fabs(metJetEta)<1.5)    jetScaleFactor = 1.015;
  else if(fabs(metJetEta)<3) jetScaleFactor = 1.04;
  else                       jetScaleFactor = 1.15;

  //std::cout << "--> metCorr: scale factor = " << jetScaleFactor << std::endl;
  return jetScaleFactor;
}

