#include "interface/xAna.hh"
#include "interface/MinitreeOut.hh"

#include <iostream>
using namespace std;



bool xAna::isZgamma(void) {

  for ( int i=0; i<nMC; ++i) { 
    if( ((*mcPID)[i] == 22 && (*mcMomPID)[i] == 22 &&  (*mcGMomPID)[i] == 22) || 
	((*mcPID)[i] == 22 && (*mcMomPID)[i] == 22 && ( (*mcGMomPID)[i] == 21 ||  fabs((*mcGMomPID)[i])<6) ) ||
	((*mcPID)[i] == 22 && (*mcMomPID)[i] == 22 &&  ( (*mcGMomPID)[i] == 23 || fabs((*mcGMomPID)[i])==24)) ||
	((*mcPID)[i] == 22 &&( (*mcMomPID)[i] == 23 || fabs((*mcMomPID)[i])==24)) ||
	((*mcPID)[i] == 22 && (fabs((*mcMomPID)[i]) == 11 || fabs((*mcMomPID)[i]) == 13 ) && (fabs((*mcGMomPID)[i]) == 24 || fabs((*mcGMomPID)[i]) == 23 ) ) 
	){
      cout << "----- ientry = " << -1 <<" removing  Zgamma  -------------" <<endl;
      cout <<i<<"  (*mcPID)[i] "<<   (*mcPID)[i]  <<   "  (*mcMomPID)[i] "<<   (*mcMomPID)[i]   <<"  (*mcGMomPID)[i]  "<< (*mcGMomPID)[i] << endl;
      return true;
    } // if DC
  } // loop over mc
  // mode Zjets
  return false;
}


void xAna::fillMCtruthInfo_HiggsSignal( void ) {
  
  if( mcHardPt == 0 ) {
    
    int hind=-1;
    vector<int> gind;
    for ( int i=0; i<nMC; ++i) { 
      if(       (*mcPID)[i]    == 25 || (*mcPID)[i]  == 39 ) hind = i;         //higgs 
      else if(  (*mcMomPID)[i] == 25 && (*mcPID)[i]  == 22 ) gind.push_back(i); //gammas
      else if(  (*mcMomPID)[i] == 39 && (*mcPID)[i]  == 22 ) gind.push_back(i); //gammas
    }
    
    if( gind.size() != unsigned(2) ) {
      cerr << " fillMCtruthInfo_HiggsSignal[ERROR]: # of photons != 2, nphotons = "
	   << gind.size() << endl;
      if( gind.size() < 2 ) return;    
    }

    TLorentzVector g[2],h;
    
    for( unsigned ip = 0 ; ip < 2; ip++ ) {
      g[ip].SetPtEtaPhiE( (*mcPt)[ gind[ip]], (*mcEta)[gind[ip]], 
			  (*mcPhi)[gind[ip]], (*mcE)[  gind[ip]] );
    }

    if( hind >= 0 )  h.SetPtEtaPhiE( (*mcPt)[hind], (*mcEta)[hind], (*mcPhi)[hind], (*mcE)[hind] );
    else             h = g[0]+g[1];
    if( hind >=0 && fabs((g[0]+g[1]).M() - h.M() ) > 0.01 )  {    
      cerr << " fillMCtruthInfo_HiggsSignal[ERROR]: the diphoton mass is not the same mass as the Higgs mass " << endl
	   << "  - higgs mass   : " << h.M() << endl
	   << "  - diphoton mass: " << (g[0]+g[1]).M() << endl;
      //return;
    }
    _minitree->mc_cThetaStar_CS =     2*(g[0].E()*g[1].Pz() - g[1].E()*g[0].Pz())/(h.M()*sqrt(h.M2()+h.Pt()*h.Pt()));
    _minitree->mc_mH = h.M();
  
  } else {

    vector<TLorentzVector> pin;
    vector<TLorentzVector> pou;
    TLorentzVector h, inPartons;



    cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>       MC truth info    mcHardPt.size() = " << mcHardPt->size() << endl; 

    for( unsigned ihard = 0 ; ihard < mcHardPt->size(); ihard++ ) {
      TLorentzVector ptmp;
      ptmp.SetPtEtaPhiM( (*mcHardPt)[ihard], (*mcHardEta)[ihard],(*mcHardPhi)[ihard], (*mcHardM)[ihard] );
      if     ( (*mcHardFun)[ihard] == -1 ) pin.push_back( ptmp );
      else if( (*mcHardFun)[ihard] == +1 ) pou.push_back( ptmp );
      else if( (*mcHardFun)[ihard] ==  0 ) h = ptmp;  

      cout << "  ihard = " << ihard << " mcHardPID = " << (*mcHardPID)[ihard] << endl;
    }
    
    if( pou.size() != unsigned(2) ) {
      cerr << " fillMCtruthInfo_HiggsSignal[ERROR]: # of photons != 2, nphotons = "
	   << pou.size() << endl;
      if( pou.size() < 2 ) return;    
    }

    // some of the variables used to unboost are present only on JHU ntuples produced by CF and are missing in FC's files
    if( mcHardOutPt != 0 ) { 

      //vector<TLorentzVector> partonout;

      cout << "  hardout.size() = " << mcHardOutPt->size() << endl;


      for( size_t ipart = 0; ipart < mcHardOutPt->size(); ipart++ ) {
	TLorentzVector partmp;
	partmp.SetPtEtaPhiE( (*mcHardOutPt)[ipart], (*mcHardOutEta)[ipart], (*mcHardOutPhi)[ipart], (*mcHardOutE)[ipart] );
	pou.push_back( partmp );

	cout << "   ihardout = " << ipart << " id = " << (*mcHardOutPdgId)[ipart] << endl;
      }
      
      inPartons = pin[0] + pin[1];
      // FILL UNBOOSTED VARIABLES HERE..code for boosting in the macro?
      _minitree->mc_H_Pt  = h.Pt();
      _minitree->mc_H_Eta = h.Eta();
      _minitree->mc_H_Phi = h.Phi();
      _minitree->mc_H_E   = h.E();   
      _minitree->mc_partonOut_Pt1  = pou[2].Pt();
      _minitree->mc_partonOut_Pt2  = pou[3].Pt(); 
      _minitree->mc_partonOut_Phi1 = pou[2].Phi();
      _minitree->mc_partonOut_Phi2 = pou[3].Phi();  
      _minitree->mc_partonOut_Eta1 = pou[2].Eta();
      _minitree->mc_partonOut_Eta2 = pou[3].Eta(); 
      _minitree->mc_partonOut_E1   = pou[2].E();
      _minitree->mc_partonOut_E2   = pou[3].E();            
      
      _minitree->mc_partonOut_deltaPhi = fabs(pou[2].Phi()-pou[3].Phi())<TMath::Pi()? fabs(pou[2].Phi()-pou[3].Phi()):2*TMath::Pi()-fabs(pou[2].Phi()-pou[3].Phi());
      //cout << "       deltaPhi    = " <<  _minitree->mc_partonOut_deltaPhi << "    p = " << h.P() << endl;
      pou[2].Boost( -inPartons.BoostVector() );
      pou[3].Boost( -inPartons.BoostVector() );
      _minitree->mc_partonOut_deltaPhiRF = fabs(pou[2].Phi()-pou[3].Phi())<TMath::Pi()? fabs(pou[2].Phi()-pou[3].Phi()):2*TMath::Pi()-fabs(pou[2].Phi()-pou[3].Phi());
      //cout << "       deltaPhi RF = " <<  _minitree->mc_partonOut_deltaPhi << endl;
      
      //double cosStar = -99;
      if( pin.size() > 0 ) {
	/// boost to Higgs rest frame
	//      cout << "========================" << endl;
	//      cout << " pz[IN] : pz[0] + pz[1]= " << (pin[0] + pin[1]).Pz() << endl;
	//      cout << " pz[OUT]: pz[0] + pz[1]= " << (pou[0] + pou[1]).Pz() << endl;
	//      cout << "                  pz[H]= " << h.Pz() << endl;
	pou[0].Boost( -h.BoostVector() );
	pin[0].Boost( -h.BoostVector() );
	//      pin[1].Boost( -h.BoostVector() );
	
	//TVector3 zBeam = pin[0].Vect();
	//cosStar = fabs(cos( pou[0].Angle(zBeam)));
      }
      _minitree->mc_cThetaStar_CS = 2*(pou[0].E()*pou[1].Pz() - pou[1].E()*pou[0].Pz())/(h.M()*sqrt(h.M2()+h.Pt()*h.Pt()));
      //_minitree->mc_cThetaStar    = cosStar;
      _minitree->mc_mH = h.M();
    }
    
  }
  
}


// ***************************************************************************************
// ********************** Reweighting functions
// ***************************************************************************************

void xAna::rescalePhoIdVar( int i ) {    
  
  if       ( _config->setup() == "Legacy2013_7TeV" && _config->analysisType() == "baseline" ) {
    if(fabs((*phoSCEta)[i])<1.479)  (*phoR9)[i] = (*phoR9)[i]*1.0048;
    else                         (*phoR9)[i] = (*phoR9)[i]*1.00492;
  } else if( _config->setup() == "Legacy2013_7TeV" && _config->analysisType() == "MVA"      ) {
    (*phoR9)[i]=(*phoR9)[i]*1.0035;
    if( fabs((*phoSCEta)[i])<1.479  ) (*phoSigmaIEtaIEta)[i] = (*phoSigmaIEtaIEta)[i]*0.87+0.0011;
    else                           (*phoSigmaIEtaIEta)[i] = (*phoSigmaIEtaIEta)[i]*0.99;
    (*phoSCEtaWidth)[i] *= 0.98;
    (*phoSCPhiWidth)[i] *= 0.99;
  } else if( (_config->setup() == "Prompt2012_ichep" || _config->setup() == "Legacy2013_8TeV") && 
	     _config->analysisType() == "baselineCiC4PF" ) {
    /// need to be clarify: gglobe ue 
    if(fabs((*phoSCEta)[i])<1.479)   {
      (*phoR9)[i]=(*phoR9)[i]*1.0053;
      (*phoSigmaIEtaIEta)[i]=(*phoSigmaIEtaIEta)[i]*0.998;
    } else{
      (*phoR9)[i]=(*phoR9)[i]*1.0075;
      (*phoSigmaIEtaIEta)[i]=(*phoSigmaIEtaIEta)[i]*1.00;
    }
  } else if( (_config->setup() == "Prompt2012_ichep" || _config->setup() == "Legacy2013_8TeV")&& _config->analysisType() == "MVA" ) {
    if(fabs((*phoSCEta)[i])<1.479) {
      (*phoR9)[i]            = 1.0045*(*phoR9)[i] + 0.0010;
      (*phoSigmaIEtaIEta)[i] = 0.891832 * (*phoSigmaIEtaIEta)[i] + 0.0009133;
      (*phoSCPhiWidth)[i]    = 1.00002 * (*phoSCPhiWidth)[i] - 0.000371;
      (*phoSCEtaWidth)[i]    = 1.04302 * (*phoSCEtaWidth)[i] - 0.000618;
     // (*phoS4ratio)[i]       = 1.01894 * (*phoS4ratio)[i]    - 0.01034;
    } else {
      (*phoR9)[i]            = 1.0086*(*phoR9)[i]- 0.0007 ;
      (*phoSigmaIEtaIEta)[i] = 0.99470 * (*phoSigmaIEtaIEta)[i]+ 0.00003;
      (*phoSCPhiWidth)[i]    = 0.99992 * (*phoSCPhiWidth)[i] - 0.00000048;
      (*phoSCEtaWidth)[i]    = 0.903254 * (*phoSCEtaWidth)[i] + 0.001346;
    //  (*phoS4ratio)[i]       = 1.04969 * (*phoS4ratio)[i]    - 0.03642;
    }
  } else {
    throw std::runtime_error(string("Supported configs (setup,analysisType): \n") +
			     " (Legacy2013_7TeV,baseline) ; (Legacy2013_7TeV,MVA) ; (Legacy2013_8TeV,MVA);(Prompt2012_ichep,baselineCiC4PF) ; (Prompt2012_ichep,MVA) ") ;
  }
}

void xAna::rescalePhoIdVarBaseline( int i ) {
  if(fabs((*phoSCEta)[i])<1.479)
    (*phoR9)[i]=(*phoR9)[i]*1.0048;
  else
    (*phoR9)[i]=(*phoR9)[i]*1.00492;

}



float xAna::computeMCweights(const HggEvtCandidate &hcand, const bool isprompt) {
  //------------ MC weighting factors and corrections   
  /// needs to be recomputed for each event (because reinit var for each event)
  _minitree->mc_wNgen = 100000./_weight_manager->getNevts();
  if( ( mode_ == 3 ||  mode_ == 19 ) && isprompt == 1 ) _minitree->mc_wXsec  *= 1.3;
  _minitree->mc_wBSz   = _weight_manager->bszW((*vtxbs_z)[hcand.vertexIdx]-(*mcVtx_z)[0]);
  _minitree->mc_wVtxId = _weight_manager->vtxPtCorrW( (hcand.gammaK_lead+hcand.gammaK_trail) .Pt() );
  
  /// photon identification
  float wPhoId[] = { 1., 1.};
  for( int i = 0 ; i < 2; i++ ) {
    wPhoId[i] = 1;
    int index = -1;
    if( i == 0 ) index = hcand.gammaIdx_lead;
    if( i == 1 ) index = hcand.gammaIdx_trail;
    wPhoId[i] *= _weight_manager->phoIdPresel((*phoR9)[index],(*phoSCEta)[index]);
    if( _config->analysisType() == "baselineCiC4PF" ) wPhoId[i] *= _weight_manager->phoIdCiC((*phoR9)[index],(*phoSCEta)[index]);
    if( _config->analysisType() == "MVA"            ) wPhoId[i] *= _weight_manager->phoIdMVA((*phoR9)[index],(*phoSCEta)[index]);
    if( vetoElec[i] ) wPhoId[i] *= _weight_manager->phoIdEleVeto((*phoR9)[index],(*phoSCEta)[index]);	      
  }
  _minitree->mc_wPhoEffi   = wPhoId[0]*wPhoId[1];
  
  /// trigger efficiency
  _minitree->mc_wTrigEffi = 0.9968; /// FIX ME
  
  /// cross section volontary not included in total weight
  _minitree->mc_wei = 
    _minitree->mc_wPU       *
    _minitree->mc_wBSz      *
    _minitree->mc_wHQT      *  /// = 1 but in 2011
    _minitree->mc_wVtxId    *
    _minitree->mc_wPhoEffi  *
    _minitree->mc_wTrigEffi *
    _minitree->mc_wNgen;
  
  return _minitree->mc_wei;
  
}
// ***********************************************************************************************************************
// ***********************************************************************************************************************
  

void xAna::recoToTruthMatching( HggEvtCandidate& event )
{

  vector<int>            genPhoFromHIdx;
  vector<TLorentzVector> genPhoFromH;

  float deltaRlead  = 999;
  float deltaRtrail = 999;
  int closestGenL = -999; int closestGenT = -999;
  //int closestGenIdxL = -999; int closestGenIdxT = -999;
  
  //cout << " before loop over mcE" << endl;
  for( size_t ipart = 0; ipart < mcE->size(); ipart++ ) {

    //cout << "   in loop " << ipart << endl;
    
    if( (*mcPID)[ipart] == 22 &&  (*mcMomPID)[ipart] == 25 ) {
      TLorentzVector gpho;  gpho.SetPtEtaPhiE( (*mcPt)[ipart], (*mcEta)[ipart], (*mcPhi)[ipart], (*mcE)[ipart] );
      genPhoFromH.push_back(gpho);  
      genPhoFromHIdx.push_back((int)ipart);
    }
  } 

  // **********************************

  if( genPhoFromH.size() == 2 ) {
    
  for( unsigned int igp = 0; igp < genPhoFromH.size(); igp++ ) {
    
    float deltaRleadtmp = genPhoFromH[igp].DeltaR( event.gammaK_lead );
    if( deltaRleadtmp < deltaRlead ) {
      deltaRlead = deltaRleadtmp;
      closestGenL    = igp;
      //closestGenIdxL = genPhoFromHIdx[igp];
    }
  }  // end loop on gen photons from H

  if( deltaRlead < 0.1 ) {
    _minitree->mc_genPhoEMatched[0] = genPhoFromH[closestGenL].E();
    genPhoFromH.erase( genPhoFromH.begin() + closestGenL );
    genPhoFromHIdx.erase( genPhoFromHIdx.begin() + closestGenL );
  }
  
  for( unsigned int jgp = 0; jgp < genPhoFromH.size(); jgp++ ) {
    
    float deltaRtrailtmp = genPhoFromH[jgp].DeltaR( event.gammaK_trail );
    if( deltaRtrailtmp < deltaRtrail ) {
      deltaRtrail = deltaRtrailtmp;
      closestGenT    = jgp;
      //closestGenIdxT = genPhoFromHIdx[jgp];
    }
  }  // end loop on gen photons from H
  
  if( deltaRtrail < 0.1 ) {
    _minitree->mc_genPhoEMatched[1] = genPhoFromH[closestGenT].E();
    genPhoFromH.erase( genPhoFromH.begin() + closestGenT );
  }
} // end if two gen photons from H

}


// ***********************************************************************************************************************
