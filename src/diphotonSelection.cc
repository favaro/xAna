#include "interface/xAna.hh"
#include "interface/MinitreeOut.hh"

#include <stdexcept>
#include <iostream>
using namespace std;
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


double xAna::evaluateDiPhotonMVA( int ilead, int itrail, 
				  const TLorentzVector & gKlead, const TLorentzVector & gKtrail, MassResolution & massResoCalc  )
{

  vector<int> sortedVertex = getSelectedVertex( ilead, itrail, true);
 
  if( sortedVertex.size() < 2 ) sortedVertex.push_back(-1);
  if( sortedVertex.size() < 3 ) sortedVertex.push_back(-1);
  int selVtx = sortedVertex[0];

  // MASS RESO CALC? can it be defined once for all events?

  TLorentzVector Hcandidate = gKlead + gKtrail;

  /*
  _minitree->mtree_rho   = rho2012;
  _minitree->mtree_rho25 = rho25;    
  _minitree->mtree_zVtx  = (*vtxbs_z)[selVtx];
  _minitree->mtree_nVtx  = nVtxBS;
  _minitree->mtree_nVtxNoBS = nVtx;
  
  
  _minitree->mtree_ivtx1    = selVtx;
  _minitree->mtree_ivtx2    = sortedVertex[1];
  _minitree->mtree_ivtx3    = sortedVertex[2];
  _minitree->mtree_vtxProb  = vertProb;
  _minitree->mtree_vtxMva   = vertMVA;
  
  _minitree->mtree_mass  = hcand.M();
  _minitree->mtree_pt    = hcand.Pt();
  _minitree->mtree_piT   = hcand.Pt()/hcand.M();
  _minitree->mtree_y     = hcand.Rapidity();
  */
  
  _minitree->mtree_minR9      = +999;
  _minitree->mtree_minPhoIdEB = +999;
  _minitree->mtree_minPhoIdEE = +999;
  _minitree->mtree_maxSCEta   =   -1;
  _minitree->mtree_minSCEta   = +999;
  
  
  for( int i = 0; i < 2; i++ ) {

    int igamma = -1;
    if( i == 0 ) igamma = ilead;
    if( i == 1 ) igamma = itrail;
    
    fillPhotonVariablesToMiniTree( igamma, selVtx, i ); 
    if( _minitree->mtree_r9[i] < _minitree->mtree_minR9 ) _minitree->mtree_minR9 = _minitree->mtree_r9[i];
    if( fabs( _minitree->mtree_sceta[i] ) <  1.5 &&  _minitree->mtree_mvaid[i] <  _minitree->mtree_minPhoIdEB ) _minitree->mtree_minPhoIdEB = _minitree->mtree_mvaid[i];
    if( fabs( _minitree->mtree_sceta[i] ) >= 1.5 &&  _minitree->mtree_mvaid[i] <  _minitree->mtree_minPhoIdEE ) _minitree->mtree_minPhoIdEE = _minitree->mtree_mvaid[i];
    if( fabs( _minitree->mtree_sceta[i] ) > _minitree->mtree_maxSCEta ) _minitree->mtree_maxSCEta =  fabs(_minitree->mtree_sceta[i]);
    if( fabs( _minitree->mtree_sceta[i] ) < _minitree->mtree_minSCEta ) _minitree->mtree_minSCEta =  fabs(_minitree->mtree_sceta[i]); 
  }


 //------------ compute diphoton masResEng for by varying smearing up and down ( and store output to minitree)----------------//
  float relResOverEPosErr[2];
  float relResOverENegErr[2];
  
  for( int i=0; i<2; ++i) {
    relResOverEPosErr[i] = _minitree->mtree_relResOverE[i]*1.1;
    relResOverENegErr[i] = _minitree->mtree_relResOverE[i]*0.9;
  } 
  
  // set Pos
  massResoCalc.setP4CalPosVtxResoSmear( gKlead,gKtrail, 
					TVector3((*phoSCPos_x)[ilead ], (*phoSCPos_y)[ilead ],(*phoSCPos_z)[ilead ]),
					TVector3((*phoSCPos_x)[itrail ], (*phoSCPos_y)[itrail ],(*phoSCPos_z)[itrail ]),
					TVector3((*vtxbs_x)[selVtx], (*vtxbs_y)[selVtx], (*vtxbs_z)[selVtx]),
					relResOverEPosErr, _minitree->mtree_relSmearing ) ;
  _minitree->mtree_massResoEngPosErr = massResoCalc.relativeMassResolutionFab_energy( );    
  
  // set Neg
  massResoCalc.setP4CalPosVtxResoSmear( gKlead,gKtrail, 
					TVector3((*phoSCPos_x)[ilead ], (*phoSCPos_y)[ilead ],(*phoSCPos_z)[ilead ]),
					TVector3((*phoSCPos_x)[itrail ], (*phoSCPos_y)[itrail ],(*phoSCPos_z)[itrail ]),
					TVector3((*vtxbs_x)[selVtx], (*vtxbs_y)[selVtx], (*vtxbs_z)[selVtx]),
					relResOverENegErr, _minitree->mtree_relSmearing ) ;
  _minitree->mtree_massResoEngNegErr = massResoCalc.relativeMassResolutionFab_energy( );    
  
  
  massResoCalc.setP4CalPosVtxResoSmear( gKlead, gKtrail, 
					TVector3((*phoCaloPos_x)[ilead ], (*phoCaloPos_y)[ilead ],(*phoCaloPos_z)[ilead ]),
					TVector3((*phoCaloPos_x)[itrail], (*phoCaloPos_y)[itrail],(*phoCaloPos_z)[itrail]),
					TVector3((*vtxbs_x)[selVtx], (*vtxbs_y)[selVtx], (*vtxbs_z)[selVtx]),
					_minitree->mtree_relResOverE, _minitree->mtree_relSmearing );
  
  
  _minitree->mtree_massResoTot = massResoCalc.relativeMassResolutionFab_total( vertProb );    
  _minitree->mtree_massResoEng = massResoCalc.relativeMassResolutionFab_energy( );    
  _minitree->mtree_massResoAng = massResoCalc.relativeMassResolutionFab_angular();    
  
  double diphotonmva = DiPhoID_MVA( gKlead, gKtrail, Hcandidate, massResoCalc, vertProb, 
				    _minitree->mtree_mvaid[0],_minitree->mtree_mvaid[1] );
 
  return diphotonmva;

}

// ****************************************************************************

bool xAna::selectTwoPhotons( int &i, int &j, int selVtx,
					   TLorentzVector &g1, TLorentzVector &g2) {

  if( selVtx < 0 || selVtx > nVtx ) 
    cout << " selected vertex no defined; ivtx : " << selVtx << " / " << nVtx << endl;

  /// do photonid preselection
  float etai = (*phoEtaVtx)[i][selVtx];
  float phii = (*phoPhiVtx)[i][selVtx];
  float etaj = (*phoEtaVtx)[j][selVtx];
  float phij = (*phoPhiVtx)[j][selVtx];
  g1.SetPtEtaPhiM( (*phoE)[i]/cosh(etai),etai,phii, 0 ); 
  g2.SetPtEtaPhiM( (*phoE)[j]/cosh(etaj),etaj,phij, 0 ); 
  int ilead = i, itrail = j;	  
  if( g1.Pt() < g2.Pt() ) {
    ilead = j; itrail = i;
    g1.SetPtEtaPhiM( (*phoE)[j]/cosh(etaj),etaj,phij, 0 ); 
    g2.SetPtEtaPhiM( (*phoE)[i]/cosh(etai),etai,phii, 0 ); 
  }
  
  if( photonID(ilead , selVtx, g1.Pt(), !vetoElec[0] ) == 0 ) return false; 
  if( photonID(itrail, selVtx, g2.Pt(), !vetoElec[1] ) == 0 ) return false; 

  /// apply loose mva photon id cut inside the loop for mva selection
  if ( _config->analysisType() == "MVA" ) 
  for( int ii = 0 ; ii < 2; ii++ ) {
      int index = -1;
      if( ii == 0 ) index = ilead ;
      if( ii == 1 ) index = itrail;
      
      PhoID_MVA(index,selVtx);
      float phoid = -1;
      if( fabs((*phoSCEta)[index]) < 1.45 ) phoid = phoID_mva[0]->EvaluateMVA("BDT");
      else                               phoid = phoID_mva[1]->EvaluateMVA("BDT");
      
      if( phoid < _config->photonIdMvaCut() ) return false;
  }


  i = ilead;
  j = itrail;  

  return true;
}


// JM ***************************************************************************************
// WARNING: Globe does not do it that way!

bool xAna::selectTwoPhotonsAwayFromLeptons( int i, int j,   vector<int> muonIndex, vector<int> elecIndex ) {


  if( elecIndex.size()== 0 && muonIndex.size() == 0 ) return true;

  int ilead = i, itrail = j;

  for(unsigned im=0;im<muonIndex.size();im++){
    
    int imu=muonIndex[im];
       
    double dR1 = deltaR((*phoEta)[ilead], (*phoPhi)[ilead], (*muEta)[imu], (*muPhi)[imu]);
    double dR2 = deltaR((*phoEta)[itrail], (*phoPhi)[itrail], (*muEta)[imu], (*muPhi)[imu]);
    double dRmin = min(dR1,dR2);
    
    if( dRmin < 0.5 ) return false;  

    // WARNING Selection might be looser here than in some VH categories... will be redone inside categories anyway    
  }
  
  for(unsigned ie=0;ie<elecIndex.size();ie++){
    

    int iel=elecIndex[ie];


    double dR1 = deltaR((*phoEta)[ilead], (*phoPhi)[ilead], (*eleEta)[iel], (*elePhi)[iel]);
    double dR2 = deltaR((*phoEta)[itrail], (*phoPhi)[itrail], (*eleEta)[iel], (*elePhi)[iel]);
    double dRmin = min(dR1,dR2);

    if( dRmin < 1.0 ) return false; 
    
  }
  

  return true;
}










