#include "interface/massResolutionCalculator.hh"


MassResolution::MassResolution(void) {
  _addSmearingTerm = false;
  _vtxSigZ[0] = 0.10;        /// correct vertex
  _vtxSigZ[1] = sqrt(2)*5.8; /// random  vertex
  _vtxSigZ[1] = sqrt(2)*5.0; /// random  vertex (update 15Oct2012)
  
}
void  MassResolution::addSmearingTerm( bool addSmearing  ) {
  _addSmearingTerm = addSmearing;
}

void MassResolution::setP4CalPosVtxResoSmear(const TLorentzVector & p1, const TLorentzVector &p2,
					     const TVector3   &calpos1,	const TVector3  &calpos2,
					     const TVector3 &vtx,
					     float relReso[2], float relSmear[2] ) {

  _vtx = vtx;
  for( int ie = 0 ; ie < 2 ;ie++ ) {
    _relResE[ie] = relReso[ie];
    _relSmearE[ie] = relSmear[ie];

    if( ie == 0 ) {
      _cal[0] = calpos1;
      _p4[0]  = p1;
    } else {
      _cal[1] = calpos2;
      _p4[1]  = p2;
    }       
  }
}

double MassResolution::angularTerm(int vtxType)  const  {
  if( vtxType < 0 || vtxType > 1 ) return 0;
  double sigma_dz = _vtxSigZ[vtxType];

  TVector3 epos[2];
  double r[2];
  for( int ie = 0 ; ie < 2; ie++ ) {
    epos[ie] = _cal[ie] - _vtx;
    r[ie]    = epos[ie].Mag(); 
  } 
  
  /// take this from gglobe
  /// need to cross check the actual formula
  double cos_term = TMath::Cos(epos[0].Phi()-epos[1].Phi());
  double sech1 = SecH(epos[0].Eta());
  double sech2 = SecH(epos[1].Eta());
  double tanh1 = TanH(epos[0].Eta());
  double tanh2 = TanH(epos[1].Eta());

  double numerator1 = sech1*(sech1*tanh2-tanh1*sech2*cos_term);
  double numerator2 = sech2*(sech2*tanh1-tanh2*sech1*cos_term);
  double denominator = 1. - tanh1*tanh2 - sech1*sech2*cos_term;

  double ResTerm = (-1.*sigma_dz/denominator)*(numerator1/r[0] + numerator2/r[1]);
  return ResTerm;
}


// utility functions
double MassResolution::SecH(double x)  const {
  return 1.0/TMath::CosH(x);
}

double MassResolution::TanH(double x)  const {
  return TMath::TanH(x);
}
 

double MassResolution::relativeMassResolution( int vtxType )  const {

  double massResolution2 = 0;

  //// energy term
  for( int ie = 0 ; ie < 2; ie++ ) {
    massResolution2 += _relResE[ie]*_relResE[ie];
    if( _addSmearingTerm ) 
          massResolution2 += _relSmearE[ie]*_relSmearE[ie];
  }
  /// angular term
  double alpha = _p4[0].Angle( _p4[1].Vect() );
  double angular_derivative = TMath::Sin(alpha)/(1.-TMath::Cos(alpha));
  /// FC: I don't understand why the derivative is set to 1 for wrong vertex but follow gglobe...
  if( vtxType == 1 ) angular_derivative = 1;
  double sigma_dalpha = angularTerm( vtxType ) * angular_derivative;

  massResolution2 += (sigma_dalpha*sigma_dalpha);

  //if( vtxType == 0 ) cout << "     from Roberta's function instead = " << 0.5*sqrt( massResolution2 ) << endl;

  return 0.5*sqrt( massResolution2 );
  
}

double MassResolution::relativeMassResolutionCorrVertex(void)  const {
  return relativeMassResolution(0);
}

double MassResolution::relativeMassResolutionWrongVertex(void)  const {
  return relativeMassResolution(1);
}


double MassResolution::relativeMassResolutionFab_total( double pCorrectVtx ) const {
  

  double sigma2_PVz = pCorrectVtx*_vtxSigZ[0]*_vtxSigZ[0] + (1-pCorrectVtx)*_vtxSigZ[1]*_vtxSigZ[1];
  
  /// M = sqrt ( 2 * E1/chEta1 * E2/chEta2 * (cosh Deta - cos Dphi ) )
  /// (sigmaM / M)^2 = 1/4 * { (sigmaE1 / E1)^2 + (sigmaE2 / E2)^2
  ///                          + ( sinhDeta / (cosh Deta - cos Dphi) * (dEta1/dPVz - dEta2/dPVz )
  ///                              -tanhEta1 *dEta1/dPVz -  tanhEta2 *dEta2/dPVz )^2 * sigma_PVz^2  }
  ///  dEta / dPVz = -1/ sqrt( Rcalo^2 + (zSC- PVz)^2 )

  //// energy term
  double energyTerm2 = 0;
  for( int ie = 0 ; ie < 2; ie++ ) {
    energyTerm2 += _relResE[ie]*_relResE[ie];

    if( _addSmearingTerm ) 
          energyTerm2 += _relSmearE[ie]*_relSmearE[ie];
  }

  /// add angular term
  TVector3 epos[2];
  for( int ie = 0 ; ie < 2; ie++ ) 
    epos[ie] = _cal[ie] - _vtx;
  
  float Deta = epos[0].Eta() - epos[1].Eta();
  float cos_dphi = cos( epos[0].Phi() - epos[1].Phi() ); 
  float deriv = 
    + sinh( Deta ) /(cosh(Deta)-cos_dphi) * (-1./epos[0].Mag() + 1./epos[1].Mag()) 
    - tanh( epos[0].Eta() )* (-1./epos[0].Mag() )
    - tanh( epos[1].Eta() )* (-1./epos[1].Mag() );
  
  double angularTerm2 = deriv*deriv*sigma2_PVz;

  // cout << "                       angular Term = " << angularTerm2  << endl;  

  //cout << "                       combination  = " << 0.5*sqrt( energyTerm2 +  angularTerm2) << endl;

  return 0.5 *sqrt( energyTerm2 +  angularTerm2);

}

double MassResolution::relativeMassResolutionFab_angular( void ) const {
  
  double sigma_PVz = _vtxSigZ[1]; /// only wrong vertex

  /// add angular term
  TVector3 epos[2];
  for( int ie = 0 ; ie < 2; ie++ ) 
    epos[ie] = _cal[ie] - _vtx;
  
  float Deta = epos[0].Eta() - epos[1].Eta();
  float cos_dphi = cos( epos[0].Phi() - epos[1].Phi() ); 
  float deriv = 
    + sinh( Deta ) /(cosh(Deta)-cos_dphi) * (-1./epos[0].Mag() + 1./epos[1].Mag()) 
    - tanh( epos[0].Eta() )* (-1./epos[0].Mag() )
    - tanh( epos[1].Eta() )* (-1./epos[1].Mag() );
  
  return 0.5 *fabs(deriv)*sigma_PVz;

}

double MassResolution::relativeMassResolutionFab_energy( void ) const {
  //// energy term
  double energyTerm2 = 0;
  for( int ie = 0 ; ie < 2; ie++ ) {
    energyTerm2 += _relResE[ie]*_relResE[ie];
    if( _addSmearingTerm ) 
          energyTerm2 += _relSmearE[ie]*_relSmearE[ie];
  }

  return 0.5 *sqrt( energyTerm2 );

}

