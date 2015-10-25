#include "interface/photonScaleAndSmearing.hh"
#include <cmath>
#include <iostream>
using namespace std; 

photonEnergySmearing::photonEnergySmearing( ) {}


void photonEnergySmearing::Initialize( string setupType ) {
  _oversmearing.resize(10,-1);
  _oversmearing_err.resize(10,-1);
  _oversmearingstoch_rho.resize(6,-1);
  _oversmearingstoch_rhoerr.resize(6,-1);
  _oversmearingstoch_phi.resize(6,-1);
  _oversmearingstoch_phierr.resize(6,-1);
  _oversmearingstoch_EMean.resize(6,-1);
  _oversmearingstoch_EMeanerr.resize(6,-1);
  setupType_= setupType;
  /// given in percents

  /// Jan16 smearing numbers
  /*
  _oversmearing[photonEnergySmearing::EBlowEtaGold_CM ] = 0.67;
  _oversmearing[photonEnergySmearing::EBlowEtaGold_GAP] = 0.77;
  _oversmearing[photonEnergySmearing::EBlowEtaBad_CM  ] = 1.00;
  _oversmearing[photonEnergySmearing::EBlowEtaBad_GAP ] = 0.89;
  _oversmearing[photonEnergySmearing::EBhighEtaGold   ] = 1.37;
  _oversmearing[photonEnergySmearing::EBhighEtaBad    ] = 1.88;
  _oversmearing[photonEnergySmearing::EElowEtaGold    ] = 2.68;
  _oversmearing[photonEnergySmearing::EElowEtaBad     ] = 2.79;
  _oversmearing[photonEnergySmearing::EEhighEtaGold   ] = 2.93;
  _oversmearing[photonEnergySmearing::EEhighEtaBad    ] = 3.01;
*/
  _oversmearing_err[photonEnergySmearing::EBlowEtaGold_CM ] = 0.36;
  _oversmearing_err[photonEnergySmearing::EBlowEtaGold_GAP] = 0.26;
  _oversmearing_err[photonEnergySmearing::EBlowEtaBad_CM  ] = 0.25;
  _oversmearing_err[photonEnergySmearing::EBlowEtaBad_GAP ] = 0.25;
  _oversmearing_err[photonEnergySmearing::EBhighEtaGold   ] = 0.73;
  _oversmearing_err[photonEnergySmearing::EBhighEtaBad    ] = 0.61;
  _oversmearing_err[photonEnergySmearing::EElowEtaGold    ] = 0.95;
  _oversmearing_err[photonEnergySmearing::EElowEtaBad     ] = 0.34;
  _oversmearing_err[photonEnergySmearing::EEhighEtaGold   ] = 0.37;
  _oversmearing_err[photonEnergySmearing::EEhighEtaBad    ] = 0.55;

  if( setupType == "Prompt2012_ichep" ) {
    /// numbers from sync http: https://twiki.cern.ch/twiki/bin/viewauth/CMS/UnblindingChecklist#Photon_energy_error 
    /* _oversmearing[photonEnergySmearing::EBlowEtaGold_GAP] = 1.03; */
    /* _oversmearing[photonEnergySmearing::EBlowEtaBad_GAP ] = 1.13; */
    /* _oversmearing[photonEnergySmearing::EBlowEtaGold_CM ] = 0.81; */
    /* _oversmearing[photonEnergySmearing::EBlowEtaBad_CM  ] = 1.26; */
    /* _oversmearing[photonEnergySmearing::EBhighEtaGold   ] = 1.95; */
    /* _oversmearing[photonEnergySmearing::EBhighEtaBad    ] = 2.22; */

    /* _oversmearing[photonEnergySmearing::EElowEtaGold    ] = 3.66; */
    /* _oversmearing[photonEnergySmearing::EElowEtaBad     ] = 3.34; */
    /* _oversmearing[photonEnergySmearing::EEhighEtaGold   ] = 5.28; */
    /* _oversmearing[photonEnergySmearing::EEhighEtaBad    ] = 5.58; */

    /// new numbers from HN https://hypernews.cern.ch/HyperNews/CMS/get/higgs2g/829.html
    /* _oversmearing[photonEnergySmearing::EBlowEtaGold_GAP] = 1.08; */
    /* _oversmearing[photonEnergySmearing::EBlowEtaGold_CM ] = 1.06; */
    /* _oversmearing[photonEnergySmearing::EBlowEtaBad_GAP ] = 1.19; */
    /* _oversmearing[photonEnergySmearing::EBlowEtaBad_CM  ] = 1.19; */

    /* _oversmearing[photonEnergySmearing::EBhighEtaGold   ] = 1.49; */
    /* _oversmearing[photonEnergySmearing::EBhighEtaBad    ] = 2.40; */

    /* _oversmearing[photonEnergySmearing::EElowEtaGold    ] = 3.75; */
    /* _oversmearing[photonEnergySmearing::EEhighEtaGold   ] = 6.07; */
    /* _oversmearing[photonEnergySmearing::EElowEtaBad     ] = 3.30; */
    /* _oversmearing[photonEnergySmearing::EEhighEtaBad    ] = 6.02; */

    /// new numbers for 53X: see https://twiki.cern.ch/twiki/pub/CMS/EcalEnergyScaleWithZeeForHgg/Hgg-scale-RUN2012ABC-53-1.pdf
    /* _oversmearing[photonEnergySmearing::EBlowEtaGold_GAP] = 0.98; */
    /* _oversmearing[photonEnergySmearing::EBlowEtaGold_CM ] = 0.99; */
    /* _oversmearing[photonEnergySmearing::EBlowEtaBad_GAP ] = 0.94; */
    /* _oversmearing[photonEnergySmearing::EBlowEtaBad_CM  ] = 1.13; */

    /* _oversmearing[photonEnergySmearing::EBhighEtaGold   ] = 1.95; */
    /* _oversmearing[photonEnergySmearing::EBhighEtaBad    ] = 1.98; */

    /* _oversmearing[photonEnergySmearing::EElowEtaGold    ] = 2.91; */
    /* _oversmearing[photonEnergySmearing::EElowEtaBad     ] = 2.75; */
    /* _oversmearing[photonEnergySmearing::EEhighEtaGold   ] = 3.25;     */
    /* _oversmearing[photonEnergySmearing::EEhighEtaBad    ] = 3.56; */

    //FOR MORIOND RUN ABCD SMEARING 
    _oversmearing[photonEnergySmearing::EBlowEtaGold_GAP] = 1.11;
    _oversmearing[photonEnergySmearing::EBlowEtaGold_CM ] = 1.11;
    _oversmearing[photonEnergySmearing::EBlowEtaBad_GAP ] = 1.07;
    _oversmearing[photonEnergySmearing::EBlowEtaBad_CM  ] = 1.07;
    _oversmearing[photonEnergySmearing::EBhighEtaGold   ] = 1.55;
    _oversmearing[photonEnergySmearing::EBhighEtaBad    ] = 1.94;
    
    _oversmearing[photonEnergySmearing::EElowEtaGold    ] = 2.95;
    _oversmearing[photonEnergySmearing::EElowEtaBad     ] = 2.76;
    _oversmearing[photonEnergySmearing::EEhighEtaGold   ] = 3.70;    
    _oversmearing[photonEnergySmearing::EEhighEtaBad    ] = 3.71;
    _oversmearing_err[photonEnergySmearing::EBlowEtaGold_GAP] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonEnergySmearing::EBlowEtaGold_CM ] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonEnergySmearing::EBlowEtaBad_GAP ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonEnergySmearing::EBlowEtaBad_CM  ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonEnergySmearing::EBhighEtaGold   ] = sqrt(0.40*0.40+0.60*0.60);
    _oversmearing_err[photonEnergySmearing::EBhighEtaBad    ] = sqrt(0.11*0.11+0.59*0.59);
    _oversmearing_err[photonEnergySmearing::EElowEtaGold    ] = sqrt(0.25*0.25+0.90*0.90);
    _oversmearing_err[photonEnergySmearing::EElowEtaBad     ] = sqrt(0.13*0.13+0.30*0.30);
    _oversmearing_err[photonEnergySmearing::EEhighEtaGold   ] = sqrt(0.11*0.11+0.34*0.34);
    _oversmearing_err[photonEnergySmearing::EEhighEtaBad    ] = sqrt(0.16*0.16+0.52*0.52);

  }

  if( setupType == "oversmear_hcp2012" ) {
    /// HCP2012 numbers
    _oversmearing[photonEnergySmearing::EBlowEtaGold_GAP] = 1.10;
    _oversmearing[photonEnergySmearing::EBlowEtaGold_CM ] = 1.00;
    _oversmearing[photonEnergySmearing::EBlowEtaBad_GAP ] = 1.06;
    _oversmearing[photonEnergySmearing::EBlowEtaBad_CM  ] = 1.06;
    _oversmearing[photonEnergySmearing::EBhighEtaGold   ] = 1.86;
    _oversmearing[photonEnergySmearing::EBhighEtaBad    ] = 1.96;

    _oversmearing[photonEnergySmearing::EElowEtaGold    ] = 2.83;
    _oversmearing[photonEnergySmearing::EElowEtaBad     ] = 2.67;
    _oversmearing[photonEnergySmearing::EEhighEtaGold   ] = 3.43;    
    _oversmearing[photonEnergySmearing::EEhighEtaBad    ] = 3.45;
    _oversmearing_err[photonEnergySmearing::EBlowEtaGold_GAP] = sqrt(0.14*0.14+0.22*0.22);
    _oversmearing_err[photonEnergySmearing::EBlowEtaGold_CM ] = sqrt(0.28*0.28+0.22*0.22);
    _oversmearing_err[photonEnergySmearing::EBlowEtaBad_GAP ] = sqrt(0.08*0.08+0.24*0.24);
    _oversmearing_err[photonEnergySmearing::EBlowEtaBad_CM  ] = sqrt(0.08*0.08+0.24*0.24);
    _oversmearing_err[photonEnergySmearing::EBhighEtaGold   ] = sqrt(0.42*0.42+0.60*0.60);
    _oversmearing_err[photonEnergySmearing::EBhighEtaBad    ] = sqrt(0.14*0.14+0.59*0.59);
    _oversmearing_err[photonEnergySmearing::EElowEtaGold    ] = sqrt(0.31*0.31+0.90*0.90);
    _oversmearing_err[photonEnergySmearing::EElowEtaBad     ] = sqrt(0.17*0.17+0.30*0.30);
    _oversmearing_err[photonEnergySmearing::EEhighEtaGold   ] = sqrt(0.15*0.15+0.34*0.34);
    _oversmearing_err[photonEnergySmearing::EEhighEtaBad    ] = sqrt(0.18*0.18+0.52*0.52);

  }



  if( setupType == "oversmear_moriond2013" ) {
    _oversmearing[photonEnergySmearing::EBlowEtaGold_GAP] = 1.11;
    _oversmearing[photonEnergySmearing::EBlowEtaGold_CM ] = 1.11;
    _oversmearing[photonEnergySmearing::EBlowEtaBad_GAP ] = 1.07;
    _oversmearing[photonEnergySmearing::EBlowEtaBad_CM  ] = 1.07;
    _oversmearing[photonEnergySmearing::EBhighEtaGold   ] = 1.55;
    _oversmearing[photonEnergySmearing::EBhighEtaBad    ] = 1.94;
    
    _oversmearing[photonEnergySmearing::EElowEtaGold    ] = 2.95;
    _oversmearing[photonEnergySmearing::EElowEtaBad     ] = 2.76;
    _oversmearing[photonEnergySmearing::EEhighEtaGold   ] = 3.70;    
    _oversmearing[photonEnergySmearing::EEhighEtaBad    ] = 3.71;
    _oversmearing_err[photonEnergySmearing::EBlowEtaGold_GAP] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonEnergySmearing::EBlowEtaGold_CM ] = sqrt(0.07*0.07+0.22*0.22);
    _oversmearing_err[photonEnergySmearing::EBlowEtaBad_GAP ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonEnergySmearing::EBlowEtaBad_CM  ] = sqrt(0.06*0.06+0.24*0.24);
    _oversmearing_err[photonEnergySmearing::EBhighEtaGold   ] = sqrt(0.40*0.40+0.60*0.60);
    _oversmearing_err[photonEnergySmearing::EBhighEtaBad    ] = sqrt(0.11*0.11+0.59*0.59);
    _oversmearing_err[photonEnergySmearing::EElowEtaGold    ] = sqrt(0.25*0.25+0.90*0.90);
    _oversmearing_err[photonEnergySmearing::EElowEtaBad     ] = sqrt(0.13*0.13+0.30*0.30);
    _oversmearing_err[photonEnergySmearing::EEhighEtaGold   ] = sqrt(0.11*0.11+0.34*0.34);
    _oversmearing_err[photonEnergySmearing::EEhighEtaBad    ] = sqrt(0.16*0.16+0.52*0.52);

  }

  if( setupType == "Legacy2013_8TeV" ) {
    
    _oversmearing[photonEnergySmearing::EBlowEtaGold_GAP] = 0.75;
    _oversmearing[photonEnergySmearing::EBlowEtaGold_CM ] = 0.75;
    _oversmearing[photonEnergySmearing::EBlowEtaBad_GAP ] = 0.86;
    _oversmearing[photonEnergySmearing::EBlowEtaBad_CM  ] = 0.86;
    _oversmearing[photonEnergySmearing::EBhighEtaGold   ] = 1.22;
    _oversmearing[photonEnergySmearing::EBhighEtaBad    ] = 1.88;
    _oversmearing[photonEnergySmearing::EElowEtaGold    ] = 1.63;
    _oversmearing[photonEnergySmearing::EElowEtaBad     ] = 1.98;
    _oversmearing[photonEnergySmearing::EEhighEtaGold   ] = 1.86;
    _oversmearing[photonEnergySmearing::EEhighEtaBad    ] = 1.92;
    _oversmearing_err[photonEnergySmearing::EBlowEtaGold_GAP] = sqrt(0.03*0.03+0.23*0.23);
    _oversmearing_err[photonEnergySmearing::EBlowEtaGold_CM ] = sqrt(0.03*0.03+0.23*0.23);
    _oversmearing_err[photonEnergySmearing::EBlowEtaBad_GAP ] = sqrt(0.02*0.02+0.25*0.25);
    _oversmearing_err[photonEnergySmearing::EBlowEtaBad_CM  ] = sqrt(0.02*0.02+0.25*0.25);
    _oversmearing_err[photonEnergySmearing::EBhighEtaGold   ] = sqrt(0.09*0.09+0.72*0.72);
    _oversmearing_err[photonEnergySmearing::EBhighEtaBad    ] = sqrt(0.02*0.02+0.60*0.60);
    _oversmearing_err[photonEnergySmearing::EElowEtaGold    ] = 0.903;
    _oversmearing_err[photonEnergySmearing::EElowEtaBad     ] = 0.301;
    _oversmearing_err[photonEnergySmearing::EEhighEtaGold   ] = 0.341;
    _oversmearing_err[photonEnergySmearing::EEhighEtaBad    ] = 0.522;

   //above are previously derived numbers for old smearing method 


   //here compute stochastic smearing term  rho*meansqrtEt*cosphi
   _oversmearingstoch_rho[photonEnergySmearing::EBlowEtaGold_GAP] =0.74;
   _oversmearingstoch_rhoerr[photonEnergySmearing::EBlowEtaGold_GAP]=0.019;
   _oversmearingstoch_rho[photonEnergySmearing::EBlowEtaGold_CM]=0.74;
   _oversmearingstoch_rhoerr[photonEnergySmearing::EBlowEtaGold_CM]=0.019;
   _oversmearingstoch_rho[photonEnergySmearing::EBlowEtaBad_GAP]=0.77;
   _oversmearingstoch_rhoerr[photonEnergySmearing::EBlowEtaBad_GAP]=0.0145;
   _oversmearingstoch_rho[photonEnergySmearing::EBlowEtaBad_CM]=0.77;
   _oversmearingstoch_rhoerr[photonEnergySmearing::EBlowEtaBad_CM]=0.0145;
   
   _oversmearingstoch_rho[photonEnergySmearing::EBhighEtaGold   ]=1.12;
   _oversmearingstoch_rhoerr[photonEnergySmearing::EBhighEtaGold   ]=0.0704;
   _oversmearingstoch_rho[photonEnergySmearing::EBhighEtaBad   ]=1.26;
   _oversmearingstoch_rhoerr[photonEnergySmearing::EBhighEtaBad ]=0.0204;
      
   _oversmearingstoch_phi[photonEnergySmearing::EBlowEtaGold_GAP]=0.0;
   _oversmearingstoch_phierr[photonEnergySmearing::EBlowEtaGold_GAP]=15.6;
   _oversmearingstoch_phi[photonEnergySmearing::EBlowEtaGold_CM]=0.0;
   _oversmearingstoch_phierr[photonEnergySmearing::EBlowEtaGold_CM]=15.6;
   _oversmearingstoch_phi[photonEnergySmearing::EBlowEtaBad_GAP]=0.0;
   _oversmearingstoch_phierr[photonEnergySmearing::EBlowEtaBad_GAP]=16.2;
   _oversmearingstoch_phi[photonEnergySmearing::EBlowEtaBad_CM]=0.0;
   _oversmearingstoch_phierr[photonEnergySmearing::EBlowEtaBad_CM]=16.2;   

   _oversmearingstoch_phi[photonEnergySmearing::EBhighEtaGold   ]=0.0;
   _oversmearingstoch_phierr[photonEnergySmearing::EBhighEtaGold   ]=22.3;
   _oversmearingstoch_phi[photonEnergySmearing::EBhighEtaBad   ]=0.0;
   _oversmearingstoch_phierr[photonEnergySmearing::EBhighEtaBad ]=7.17;

//these are not in percent absolute energy
   _oversmearingstoch_EMean[photonEnergySmearing::EBlowEtaGold_GAP]=6.60119;
   _oversmearingstoch_EMeanerr[photonEnergySmearing::EBlowEtaGold_GAP]=0.31923;
   _oversmearingstoch_EMean[photonEnergySmearing::EBlowEtaGold_CM]=6.60119;
   _oversmearingstoch_EMeanerr[photonEnergySmearing::EBlowEtaGold_CM]=0.31923;
   _oversmearingstoch_EMean[photonEnergySmearing::EBlowEtaBad_GAP]=6.7276;
   _oversmearingstoch_EMeanerr[photonEnergySmearing::EBlowEtaBad_GAP]=0.347156;
   _oversmearingstoch_EMean[photonEnergySmearing::EBlowEtaBad_CM]=6.7276;
   _oversmearingstoch_EMeanerr[photonEnergySmearing::EBlowEtaBad_CM]=0.347156;

   _oversmearingstoch_EMean[photonEnergySmearing::EBhighEtaGold   ]=6.52397;
   _oversmearingstoch_EMeanerr[photonEnergySmearing::EBhighEtaGold   ]=0.584685;
   _oversmearingstoch_EMean[photonEnergySmearing::EBhighEtaBad   ]=6.73261;
   _oversmearingstoch_EMeanerr[photonEnergySmearing::EBhighEtaBad ]=0.29007;  
   
  }

  if( setupType == "Legacy2013_7TeV" ) {
    _oversmearing[photonEnergySmearing::EBlowEtaGold_GAP] = 0.68;
    _oversmearing[photonEnergySmearing::EBlowEtaGold_CM ] = 0.68;
    _oversmearing[photonEnergySmearing::EBlowEtaBad_GAP ] = 0.96;
    _oversmearing[photonEnergySmearing::EBlowEtaBad_CM  ] = 0.96;
    _oversmearing[photonEnergySmearing::EBhighEtaGold   ] = 1.01;
    _oversmearing[photonEnergySmearing::EBhighEtaBad    ] = 1.85;

    _oversmearing[photonEnergySmearing::EElowEtaGold    ] = 1.58;
    _oversmearing[photonEnergySmearing::EElowEtaBad     ] = 1.85;
    _oversmearing[photonEnergySmearing::EEhighEtaGold   ] = 1.83;
    _oversmearing[photonEnergySmearing::EEhighEtaBad    ] = 2.01;
    //store error (absolute proofed against globe file)
    _oversmearing_err[photonEnergySmearing::EBlowEtaGold_GAP] = 0.221;
    _oversmearing_err[photonEnergySmearing::EBlowEtaGold_CM ] = 0.221;
    _oversmearing_err[photonEnergySmearing::EBlowEtaBad_GAP ] = 0.241;
    _oversmearing_err[photonEnergySmearing::EBlowEtaBad_CM  ] = 0.241;
    _oversmearing_err[photonEnergySmearing::EBhighEtaGold   ] = 0.616;
    _oversmearing_err[photonEnergySmearing::EBhighEtaBad    ] = 0.591;
    _oversmearing_err[photonEnergySmearing::EElowEtaGold    ] = 0.903;
    _oversmearing_err[photonEnergySmearing::EElowEtaBad     ] = 0.301;
    _oversmearing_err[photonEnergySmearing::EEhighEtaGold   ] = 0.341;
    _oversmearing_err[photonEnergySmearing::EEhighEtaBad    ] = 0.522;
  }


  

  /// absolute:
  for( unsigned ios = 0 ; ios < _oversmearing.size()    ; ios++ ) _oversmearing[ios]     /= 100.;
  for( unsigned ios = 0 ; ios < _oversmearing_err.size(); ios++ ) _oversmearing_err[ios] /= 100.;
  for( unsigned ios = 0 ; ios <_oversmearingstoch_rho.size(); ios++ ) {
  _oversmearingstoch_rho[ios]/=100.;
  _oversmearingstoch_rhoerr[ios]/=100.;
  _oversmearingstoch_phi[ios]/=100.;
  _oversmearingstoch_phierr[ios]/=100.;
}
}

void  photonEnergySmearing::setSmearSeed(int evtNum, int runNum, int lumisec, float SCPhi, float SCEta){
        UInt_t seedBase = (UInt_t) evtNum + (UInt_t) runNum + (UInt_t) lumisec;
        UInt_t seed1    = seedBase + 100000*(UInt_t) (TMath::Abs(100.*SCPhi)) +1000*(UInt_t) (TMath::Abs(100.*SCEta));
       overSmearingGen.SetSeed(seed1);       
}

photonEnergySmearing::category photonEnergySmearing::photonEnergySmearingCategory( float scEta, float r9, bool isEBGap  ) {
  scEta = fabs(scEta);
  if     ( scEta < 1.0 ) {
    if( r9 > 0.94 ) {
      if( isEBGap ) return photonEnergySmearing::EBlowEtaGold_GAP;
      else          return photonEnergySmearing::EBlowEtaGold_CM;
    } else {
      if( isEBGap ) return photonEnergySmearing::EBlowEtaBad_GAP;
      else          return photonEnergySmearing::EBlowEtaBad_CM;
    }
  } else if( scEta < 1.45 )  {
    if( r9 > 0.94 ) return photonEnergySmearing::EBhighEtaGold;
    else            return photonEnergySmearing::EBhighEtaBad;
  } else if( scEta < 2.0 ) {
    if( r9 > 0.94 ) return photonEnergySmearing::EElowEtaGold;
    else            return photonEnergySmearing::EElowEtaBad;
  } else if( scEta < 2.5 ) {
    if( r9 > 0.94 ) return photonEnergySmearing::EEhighEtaGold;
    else            return photonEnergySmearing::EEhighEtaBad;
  }

  return  photonEnergySmearing::UNKNOWN;
}


float photonEnergySmearing::meanOverSmearing( float scEta, float r9, bool isEBGap, int syst) {
  photonEnergySmearing::category cat = photonEnergySmearingCategory(scEta, r9, isEBGap);
  if( cat == photonEnergySmearing::UNKNOWN ) return 0;  
  
  return _oversmearing[cat] + syst * _oversmearing_err[cat];
}

float photonEnergySmearing::meanOverSmearing( float scEta, float r9, float et, bool isEBGap, int syst) {
  photonEnergySmearing::category cat = photonEnergySmearingCategory(scEta, r9, isEBGap);
  if( cat == photonEnergySmearing::UNKNOWN ) return 0;
  float phi   = _oversmearingstoch_phi[cat]+ syst*_oversmearingstoch_phierr[cat];
  float rho   = _oversmearingstoch_rho[cat]+ syst*_oversmearingstoch_rhoerr[cat];
  float EMean = _oversmearingstoch_EMean[cat]+ syst*_oversmearingstoch_EMeanerr[cat];

  float constSmear = rho*sin(phi);
  float stochSmear = rho*cos(phi)*EMean;   
  if(cat == photonEnergySmearing::EElowEtaGold  || cat == photonEnergySmearing::EElowEtaBad ||
     cat == photonEnergySmearing::EEhighEtaGold || cat == photonEnergySmearing::EEhighEtaBad ) 
    return _oversmearing[cat] + syst * _oversmearing_err[cat];
  return sqrt( constSmear*constSmear + (stochSmear*stochSmear)/et );
}


float photonEnergySmearing::randOverSmearing( float scEta, float r9,  float et, bool isEBGap, int syst) {
  float oversmearing = meanOverSmearing( scEta, r9, isEBGap, syst );
  if(setupType_== "Legacy2013_8TeV"){
  oversmearing=meanOverSmearing( scEta, r9, et,isEBGap, syst );
}

  Double_t scale=-1.0;
  while(scale<0.0){
   scale=overSmearingGen.Gaus( 1.0, oversmearing);
  }

  return scale;
}



