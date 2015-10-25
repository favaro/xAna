#ifndef massResoCalc_hh__
#define massResoCalc_hh__

#include <TLorentzVector.h>

class MassResolution {
 public:
  MassResolution(void);
  void addSmearingTerm( bool addSmearing = true );
  void setP4CalPosVtxResoSmear(const TLorentzVector & p1, const TLorentzVector &p2,
                               const TVector3   &calpos1, const TVector3  &calpos2,
                               const TVector3 &vtx,
                               float relReso[2], float relSmear[2] );


  /// vtxType:
  /// - -1: no angular term
  /// - 0: use correct vertex dZ RMS
  /// - 1: use random vertex dZ RMS
  double angularTerm(int vtxType = 0) const ;
  double relativeMassResolution( int vtxType ) const ;
  double relativeMassResolutionCorrVertex( ) const ;
  double relativeMassResolutionWrongVertex( ) const ;

  /// FC calc.
  double relativeMassResolutionFab_total( double pCorrectVtx )  const ;
  double relativeMassResolutionFab_angular( void ) const;
  double relativeMassResolutionFab_energy(  void ) const;

 private:
  TLorentzVector _p4[2];
  TVector3 _cal[2];
  TVector3 _vtx;
  double _relResE[2];
  double _relSmearE[2];
  double _vtxSigZ[2];

  bool  _addSmearingTerm;

  double SecH(double x) const ;
  double TanH(double x) const ;

};

#endif
