#include "interface/xAna.hh"
#include "interface/hggCandidate.hh"

#include <TString.h>

using namespace std;

float xAna::stringToFloat( const TString varname, HggEvtCandidate& hggevent ) {
  float varName = 0;

  TLorentzVector Hcand = hggevent.gammaK_lead + hggevent.gammaK_trail; 

  if( varname == "diphotonMva" )        varName = hggevent.diphotonMva;
  else if( varname == "nmuon" )         varName = (float)hggevent.nmu;
  else if( varname == "nelec" )         varName = (float)hggevent.nelec;
 
  // ttH   ********************************************
  
  else if( varname == "nlep" )         varName = (float)hggevent.nelec+(float)hggevent.nmu;
  else if( varname == "njetsTTH" )     varName = (float)hggevent.njetsTTH; 
  else if( varname == "njetsVH" )      varName = (float)hggevent.njetsVH; 
  else if( varname == "nBjets" )       varName = (float)hggevent.nBjets; 

  
  //  VH ********************************************
  
  //   JM:  now deltaR cut is also done in the selectTwoPhotonsAwayFromLeptons function in UtilsDipho

  else if( varname == "deltaR_leptonToLead" ) {
    TLorentzVector lepton;  
    float min=999;
    for( int nm=0;nm<hggevent.nmu;nm++){
      lepton.SetPtEtaPhiM( (*muPt)[hggevent.muIdx[nm]], (*muEta)[hggevent.muIdx[nm]], (*muPhi)[hggevent.muIdx[nm]], 0.105 );
      if(lepton.DeltaR( hggevent.gammaK_lead )<min) min=lepton.DeltaR( hggevent.gammaK_lead );
    }
     for( int ne=0;ne<hggevent.nelec;ne++){
       lepton.SetPtEtaPhiM( (*elePt)[hggevent.eleIdx[ne]], (*eleEta)[hggevent.eleIdx[ne]], (*elePhi)[hggevent.eleIdx[ne]], 0.0005 );
       if(lepton.DeltaR( hggevent.gammaK_lead )<min) min=lepton.DeltaR( hggevent.gammaK_lead );
     }
     varName = min; 
  } 
  else if( varname == "deltaR_leptonToTrail" ) {
    TLorentzVector lepton;  
    float min=999;
    for( int nm=0;nm<hggevent.nmu;nm++){
      lepton.SetPtEtaPhiM( (*muPt)[hggevent.muIdx[nm]], (*muEta)[hggevent.muIdx[nm]], (*muPhi)[hggevent.muIdx[nm]], 0.105 );
      if(lepton.DeltaR( hggevent.gammaK_trail )<min) min=lepton.DeltaR( hggevent.gammaK_trail );
    }
     for( int ne=0;ne<hggevent.nelec;ne++){
       lepton.SetPtEtaPhiM( (*elePt)[hggevent.eleIdx[ne]], (*eleEta)[hggevent.eleIdx[ne]], (*elePhi)[hggevent.eleIdx[ne]], 0.0005 );
       if(lepton.DeltaR( hggevent.gammaK_trail )<min) min=lepton.DeltaR( hggevent.gammaK_trail );
     }
     varName = min;
  }
  
  else if( varname == "nmuonVHDilep") varName = (float) hggevent.nmuVHDilep;
  else if( varname == "nelecVHDilep") varName = (float) hggevent.nelecVHDilep;
  else if( varname == "mllVHDilep"){ 
    //JM: fixme - not nice, can't find a better way to do this for now 
    TLorentzVector lep1;  TLorentzVector lep2; TLorentzVector leplep;    
    if(hggevent.nmuVHDilep>=2){
      lep1.SetPtEtaPhiM( (*muPt)[hggevent.muIdxVHDilep[0]], (*muEta)[hggevent.muIdxVHDilep[0]], (*muPhi)[hggevent.muIdxVHDilep[0]], 0.105 );
      lep2.SetPtEtaPhiM( (*muPt)[hggevent.muIdxVHDilep[1]], (*muEta)[hggevent.muIdxVHDilep[1]], (*muPhi)[hggevent.muIdxVHDilep[1]], 0.105 );
      leplep=lep1+lep2;
      varName=leplep.M();   
    } else if (hggevent.nelecVHDilep>=2){
      lep1.SetPtEtaPhiM( (*elePt)[hggevent.eleIdxVHDilep[0]], (*eleEta)[hggevent.eleIdxVHDilep[0]], (*elePhi)[hggevent.eleIdxVHDilep[0]], 0.0005 );
      lep2.SetPtEtaPhiM( (*elePt)[hggevent.eleIdxVHDilep[1]], (*eleEta)[hggevent.eleIdxVHDilep[1]], (*elePhi)[hggevent.eleIdxVHDilep[1]], 0.0005 );
      leplep=lep1+lep2;
      varName=leplep.M();
    }else{
      cout<<" WARNING: no same flavour dilepton "<< endl;
    }
  }
  else if( varname == "mEGMinusMZ"){
    float dMmin = 999;
    for( int ne=0;ne<hggevent.nelec;ne++){

      TLorentzVector lepton;  
      lepton.SetPtEtaPhiM( (*elePt)[hggevent.eleIdx[ne]], (*eleEta)[hggevent.eleIdx[ne]], (*elePhi)[hggevent.eleIdx[ne]], 0.0005 );
      TLorentzVector EG_lead=hggevent.gammaK_lead+lepton;
      TLorentzVector EG_trail=hggevent.gammaK_trail+lepton;

      float deltaM=min(fabs(EG_lead.M()-91.2),fabs(EG_trail.M()-91.2));      
      if(deltaM<dMmin) dMmin=deltaM;
    }
    varName = dMmin;
  }
  else if( varname == "corMET" ) { 
    TLorentzVector MET=hggevent.corMet;
    varName = MET.Pt();
  }
  else if( varname == "deltaPhiMETGG" ) { 
    TLorentzVector MET=hggevent.corMet; 
    varName =  fabs(MET.DeltaPhi(Hcand));
  }
  else if( varname == "deltaPhiLeadJetGG"){ 

    // NB: Globe cuts on deltaPhiLeadJetMET not on deltaPhiLeadJetGG as said in the note

    if( hggevent.njets == 0 ) varName = 0;
    else if ((*jetPt)[hggevent.jetIdx[0]]< 50 ) varName = 0;
    // JM: fixme hardcoded pt cut, not very nice
    else{
      TLorentzVector jet1; 
      jet1.SetPtEtaPhiE((*jetPt)[hggevent.jetIdx[0]], (*jetEta)[hggevent.jetIdx[0]], (*jetPhi)[hggevent.jetIdx[0]], (*jetEn)[hggevent.jetIdx[0]]);      
      varName =  fabs(jet1.DeltaPhi(Hcand)); 
    }
  }
  else if( varname == "deltaPhiLeadJetMET"){ 
    
    TLorentzVector MET=hggevent.corMet; 
    if( hggevent.njets == 0 ) varName = 0;
    else if ((*jetPt)[hggevent.jetIdx[0]]< 50 ) varName = 0; 
    // JM: fixme hardcoded pt cut, not very nice
    else{
      TLorentzVector jet1; 
      jet1.SetPtEtaPhiE((*jetPt)[hggevent.jetIdx[0]], (*jetEta)[hggevent.jetIdx[0]], (*jetPhi)[hggevent.jetIdx[0]], (*jetEn)[hggevent.jetIdx[0]]);      
      varName =  fabs(jet1.DeltaPhi(MET));       
    }
  }else if( varname == "njetsVHhad"){     
    varName=hggevent.njetsVHhad;
  }else if( varname == "MJJVHhad" ){
    if(hggevent.njetsVHhad<2) varName = -1.0;
    else{
      TLorentzVector jet1, jet2; 
      jet1.SetPtEtaPhiE((*jetPt)[hggevent.jetIdxVHhad[0]], (*jetEta)[hggevent.jetIdxVHhad[0]], (*jetPhi)[hggevent.jetIdxVHhad[0]], (*jetEn)[hggevent.jetIdxVHhad[0]]);
      jet2.SetPtEtaPhiE((*jetPt)[hggevent.jetIdxVHhad[1]], (*jetEta)[hggevent.jetIdxVHhad[1]], (*jetPhi)[hggevent.jetIdxVHhad[1]], (*jetEn)[hggevent.jetIdxVHhad[1]]);
      TLorentzVector jetjet = jet1 + jet2;
      varName = jetjet.M();
    } 
  }else if( varname == "cosThetaStarVHhad"){
     if(hggevent.njetsVHhad<2) varName = -999.0;
     else{
       TLorentzVector jet1, jet2; 
       jet1.SetPtEtaPhiE((*jetPt)[hggevent.jetIdxVHhad[0]], (*jetEta)[hggevent.jetIdxVHhad[0]], (*jetPhi)[hggevent.jetIdxVHhad[0]], (*jetEn)[hggevent.jetIdxVHhad[0]]);
       jet2.SetPtEtaPhiE((*jetPt)[hggevent.jetIdxVHhad[1]], (*jetEta)[hggevent.jetIdxVHhad[1]], (*jetPhi)[hggevent.jetIdxVHhad[1]], (*jetEn)[hggevent.jetIdxVHhad[1]]);
       TLorentzVector jetjet = jet1 + jet2;
       TLorentzVector VStar = jetjet + Hcand;      

       TLorentzVector H_VStar(Hcand);
       H_VStar.Boost(-VStar.BoostVector());              
       varName = fabs(H_VStar.CosTheta());   
     }
  }

  // CiC untag ********************************************
 
  else if( varname == "minphotonR9" )     varName = TMath::Min(hggevent.gammaR9_lead, hggevent.gammaR9_trail);
  else if( varname == "maxphotonEta" )    varName = TMath::Max( fabs(hggevent.gammaK_lead.Eta()), fabs(hggevent.gammaK_trail.Eta()));
  else if( varname == "maxSCEta") {
    float scetalead  = (*phoSCEta)[hggevent.gammaIdx_lead];
    float scetatrail = (*phoSCEta)[hggevent.gammaIdx_trail];
    varName = TMath::Max( fabs(scetalead), fabs(scetatrail) );
  }
  else if( varname == "ptHoM" )           varName = Hcand.Pt()/Hcand.M();


  //  VBF ***********************************************

  else if( varname == "njets" )           varName = (float)hggevent.njets;
  else if( varname == "jetPt_lead" )      varName = (*jetPt)[hggevent.jetIdx[0]];    // assume jets are pt-ordered (true when using the selectJetJEC function)
  else if( varname == "jetPt_trail" )     varName = (*jetPt)[hggevent.jetIdx[1]];    // assume jets are pt-ordered 
  else if( varname == "Zeppenfeld" )      varName = fabs(Hcand.Eta() - ((*jetEta)[hggevent.jetIdx[0]] + (*jetEta)[hggevent.jetIdx[1]])/2.); 
  else if( varname == "DeltaEtaJJ" )      varName = fabs( (*jetEta)[hggevent.jetIdx[0]] - (*jetEta)[hggevent.jetIdx[1]] );
  else if( varname == "MJJ" ) {
    TLorentzVector jet1, jet2; 
    jet1.SetPtEtaPhiE((*jetPt)[hggevent.jetIdx[0]], (*jetEta)[hggevent.jetIdx[0]], (*jetPhi)[hggevent.jetIdx[0]], (*jetEn)[hggevent.jetIdx[0]]);
    jet2.SetPtEtaPhiE((*jetPt)[hggevent.jetIdx[1]], (*jetEta)[hggevent.jetIdx[1]], (*jetPhi)[hggevent.jetIdx[1]], (*jetEn)[hggevent.jetIdx[1]]);
    TLorentzVector jetjet = jet1 + jet2;
    varName = jetjet.M();
  }
  else if( varname == "DeltaPhiJJGG" ) {
    TLorentzVector jet1, jet2; 
    jet1.SetPtEtaPhiE((*jetPt)[hggevent.jetIdx[0]], (*jetEta)[hggevent.jetIdx[0]], (*jetPhi)[hggevent.jetIdx[0]], (*jetEn)[hggevent.jetIdx[0]]);
    jet2.SetPtEtaPhiE((*jetPt)[hggevent.jetIdx[1]], (*jetEta)[hggevent.jetIdx[1]], (*jetPhi)[hggevent.jetIdx[1]], (*jetEn)[hggevent.jetIdx[1]]);
    TLorentzVector jetjet = jet1 + jet2; 
    varName =  fabs(jetjet.DeltaPhi(Hcand));
  }
  else if( varname == "diphodijetMva" ) {
    varName = combineMVAforVBF( hggevent );
  }

 // ******************************************************

  return varName; 
}

