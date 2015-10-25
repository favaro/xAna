#include "interface/photonScaleAndSmearing.hh"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;

int photonEnergyScale::EcalPart( string ecal ) {
  
  /// internal code convention
  if( ecal == "EBlowEtaBad"   ) return 0;
  if( ecal == "EBlowEtaGold"  ) return 1;
  if( ecal == "EBhighEtaBad"  ) return 2;
  if( ecal == "EBhighEtaGold" ) return 3;
  if( ecal == "EElowEtaBad"   ) return 4;
  if( ecal == "EElowEtaGold"  ) return 5;
  if( ecal == "EEhighEtaBad"  ) return 6;
  if( ecal == "EEhighEtaGold" ) return 7;

  /// read other shervin's table format
  if( ecal == "EB-absEta_0_1-bad"      ) return 0;
  if( ecal == "EB-absEta_0_1-gold"     ) return 1;
  if( ecal == "EB-absEta_1_1.4442-bad" ) return 2;
  if( ecal == "EB-absEta_1_1.4442-gold") return 3;
  if( ecal == "EE-absEta_1.566_2-bad"  ) return 4;
  if( ecal == "EE-absEta_1.566_2-gold" ) return 5;
  if( ecal == "EE-absEta_2_2.5-bad"    ) return 6;
  if( ecal == "EE-absEta_2_2.5-gold"   ) return 7;
  
  return -1;
}

string photonEnergyScale::EcalPartString( float r9, float scEta ) {
  string ecal = "Unknown";
  if     ( fabs(scEta) <= 1.00 ) ecal = "EBlowEta";
  else if( fabs(scEta) <= 1.45 ) ecal = "EBhighEta";
  else if( fabs(scEta) <= 2.00 ) ecal = "EElowEta";
  else if( fabs(scEta) <= 2.50 ) ecal = "EEhighEta";
  
  if( r9 > 0.94 ) ecal += "Gold";
  else            ecal += "Bad";
  
  return ecal;
}

int photonEnergyScale::EcalPart(float r9, float scEta ) {
  return EcalPart( EcalPartString( r9, scEta ) );
}


bool photonEnergyScale::setup( string energyscaleName ) {

   if( energyscaleName.find( "EtDep" ) != string::npos ) {
     _isEtDep = true;
     return setup_EtDep(energyscaleName);
   }
   _isEtDep = false;

   ifstream input(energyscaleName.c_str());

   //// default energyscale values


   cout << " Opening energy scale file: " << energyscaleName<< endl;
   while ( input.good() && !input.eof() ) {
     string line;
    getline(input,line,'\n');
    
    /// comment
    if( line.find("#") != string::npos ) continue;

    istringstream isstream1(line);
    istringstream isstream2(line);
    string ecal;
    int runMin,runMax;
    float eScale, err_eScale;
    isstream1 >> ecal >> runMin >> runMax >> eScale >> err_eScale ;

    if( runMin <= 0 ) {
      string dummy;
      /// most likely we are reading the new Format
      isstream2 >> ecal >> dummy >> runMin >> runMax >> eScale >> err_eScale ;
      
    }
    int iEcal = EcalPart(ecal);
    //cout <<" iEcal "<<iEcal <<endl;

    if( iEcal < 0 ) continue;
    _runStart[  iEcal].push_back( runMin );
    _runStop[   iEcal].push_back( runMax );
    _eScale[    iEcal].push_back( eScale );
    _err_eScale[iEcal].push_back( err_eScale );
   }
   _isSetup = true;
   return true;
}

photonEnergyScale::photonEnergyScale(void) {_isSetup = false; _nWarnings = 0;}
photonEnergyScale::~photonEnergyScale(void) {}




// *************************************  combine ******************************************

float photonEnergyScale::energyScale( float r9, float scEta, float Et, int run, int systShift ) {
  
  if( _isEtDep ) 
    return energyScale_EtDep( r9, scEta, Et, run, systShift );
  return energyScale_noEtDep( r9, scEta, run, systShift );
}


// *********************************energy dependent correction******************************************

float photonEnergyScale::energyScale_EtDep( float r9, float scEta, float Et, int run, int systShift ) {
  
  if( !_isSetup ) {
    cout << " photonEnergyScale is not setup! use function setup with the ad hoc setup file" << endl;
    return 1;
  }

  if( Et < 20 || fabs( scEta) > 2.5 ) return 1;
  
  int iEcal = EcalPart_EtDep( r9, scEta, Et );
  
  //cout <<"EcalPart_EtDep( r9, scEta, Et ) "<< EcalPart_EtDep( r9, scEta, Et ) <<"  r9 "<< r9  << "  scEta "<<  scEta <<"  Et "<<Et   << endl;
  if( iEcal < 0 ) {
    //cout << " photonEnergyScale unknown Ecal part for: r9 = " << r9 << " eta = " << scEta << " et = " << Et << endl;
    return 1;
  }
  
  int irun = runRange( run, iEcal );

  if( irun < 0 ) {
    cout << " photonEnergyScale: can not find run: " << run << " for ECAL: " << EcalPartString( r9, scEta ) << endl;
    return 1;
  }
  
  //cout <<"inside energy scale --->   r9, scEta "<<  r9 <<", "<<  scEta  <<"  iEcal  "<< iEcal << "   _eScale[iEcal][irun] + systShift * _err_eScale[iEcal][irun] "<<  _eScale[iEcal][irun] + systShift * _err_eScale[iEcal][irun] <<endl;

  return _eScale[iEcal][irun] + systShift * _err_eScale[iEcal][irun];
}

// *********************************energy independent correction************************************

float photonEnergyScale::energyScale_noEtDep( float r9, float scEta, int run, int systShift ) {
  
  if( !_isSetup ) {
    cout << " photonEnergyScale is not setup! use function setup with the ad hoc setup file" << endl;
    return 1;
  }

  int iEcal = EcalPart( r9, scEta );

  if( iEcal < 0 ) {
    // cout << " photonEnergyScale unknown Ecal part for: r9 = " << r9 << " eta = " << scEta << endl;
    return 1;
  }
    

  int irun = runRange( run, iEcal );
  if( irun < 0 ) {
    cout << " photonEnergyScale: can not find run: " << run << " for ECAL: " << EcalPartString( r9, scEta ) << endl;
    return 1;
  }
  
  //  cout <<"inside energy scale --->   r9, scEta "<<  r9 <<", "<<  scEta  <<"  iEcal  "<< iEcal << "   _eScale[iEcal][irun] + systShift * _err_eScale[iEcal][irun] "<<  _eScale[iEcal][irun] + systShift * _err_eScale[iEcal][irun] <<endl;

  return _eScale[iEcal][irun] + systShift * _err_eScale[iEcal][irun];
}

// *****************************************************************************************

int photonEnergyScale::runRange( int run, int iEcal ) {
  for( unsigned irun = 0; irun < _runStart[iEcal].size(); irun++ )
    if( run >= _runStart[iEcal][irun] && run <= _runStop[iEcal][irun] ) return int(irun);


  if( run > _runStop[iEcal][_runStart[iEcal].size()-1] ) {
    if( _nWarnings < 20 ) {
      cout << "run " << run << " > last run known: " << _runStop[iEcal][_runStart[iEcal].size()-1] << " for iEcal " << iEcal << endl;
      _nWarnings++;
    }
    return _runStart[iEcal].size()-1;
  }

  return -1;
}


void testEnergyScale(int run,  string ecalCalibInput = "ecalCalibFiles/EnergyScale2012_prompt.txt") {

  
  photonEnergyScale myEnergyScale;
  myEnergyScale.setup( ecalCalibInput );

  float r9[]  = { 0.8, 0.8, 0.83, 0.85, 0.98, 0.96, 0.99, 1.05, 0.7 };
  float eta[] = { 0.5, 1.4, 1.8, 2.4, 0.6, 1.2, 1.9, 2.23, 3.0 };
  for( int itest = 0 ; itest < 9; itest++ )
    cout << " eScale( r9 = " << r9[itest] << ", eta = "<< eta[itest] << ") = " << myEnergyScale.energyScale(r9[itest],eta[itest],-1.,run) << endl;
  
}
  

//=====================================================================
//================  Et dependent corrections ================
///=====================================================================




int photonEnergyScale::EcalPart_EtDep( string ecal ) {


  if( ecal == "absEta_0_1-gold-Et_20_35"      ) return 0;
  if( ecal == "absEta_0_1-gold-Et_35_43"      ) return 1;
  if( ecal == "absEta_0_1-gold-Et_43_50"      ) return 2;
  if( ecal == "absEta_0_1-gold-Et_50_55"      ) return 3;
  if( ecal == "absEta_0_1-gold-Et_55_100"    ) return 4;

  if( ecal == "absEta_0_1-bad-Et_20_33"       ) return 5;
  if( ecal == "absEta_0_1-bad-Et_33_39"       ) return 6;
  if( ecal == "absEta_0_1-bad-Et_39_45"       ) return 7;
  if( ecal == "absEta_0_1-bad-Et_45_50"       ) return 8;
  if( ecal == "absEta_0_1-bad-Et_50_58"       ) return 9;
  if( ecal == "absEta_0_1-bad-Et_58_100"     ) return 10;


  if( ecal == "absEta_1_1.4442-gold-Et_20_40"      ) return 11;
  if( ecal == "absEta_1_1.4442-gold-Et_40_50"      ) return 12;
  if( ecal == "absEta_1_1.4442-gold-Et_50_100"    ) return 13;


  if( ecal == "absEta_1_1.4442-bad-Et_20_33"       ) return 14;
  if( ecal == "absEta_1_1.4442-bad-Et_33_39"       ) return 15;
  if( ecal == "absEta_1_1.4442-bad-Et_39_45"       ) return 16;
  if( ecal == "absEta_1_1.4442-bad-Et_45_50"       ) return 17;
  if( ecal == "absEta_1_1.4442-bad-Et_50_58"       ) return 18;
  if( ecal == "absEta_1_1.4442-bad-Et_58_100"     ) return 19;

  if( ecal == "absEta_1.566_2-gold" ) return 20;
  if( ecal == "absEta_1.566_2-bad"  ) return 21;
  if( ecal == "absEta_2_2.5-gold"   ) return 22;
  if( ecal == "absEta_2_2.5-bad"    ) return 23;
  
  return -1;
}


string photonEnergyScale::EcalPartString_EtDep( float r9, float scEta, float Et ) {
  string ecal = "Unknown";
  
  if     ( fabs(scEta) <= 1.00 ) ecal = "absEta_0_1";
  else if( fabs(scEta) <= 1.45 ) ecal = "absEta_1_1.4442";
  else if( fabs(scEta) <= 2.00 ) ecal = "absEta_1.566_2";
  else if( fabs(scEta) <= 2.50 ) ecal = "absEta_2_2.5";
  
  if( r9 > 0.94 ) ecal += "-gold";
  else            ecal += "-bad";
  
  if(ecal == "absEta_0_1-gold"){
    if(Et>=20 && Et<35)ecal +="-Et_20_35";
    else if (Et>=35 && Et<43)ecal +="-Et_35_43";
    else if (Et>=43 && Et<50)ecal +="-Et_43_50";
    else if (Et>=50 && Et<55)ecal +="-Et_50_55";
    else if (Et>=55)ecal +="-Et_55_100";
  }else if(ecal == "absEta_0_1-bad" || ecal == "absEta_1_1.4442-bad" ){
    if(Et>=20 && Et<33)ecal +="-Et_20_33";
    else if (Et>=33 && Et<39)ecal +="-Et_33_39";
    else if (Et>=39 && Et<45)ecal +="-Et_39_45";
    else if (Et>=45 && Et<50)ecal +="-Et_45_50";
    else if (Et>=50 && Et<58)ecal +="-Et_50_58";
    else if (Et>=58)ecal +="-Et_58_100";
  }else if(ecal == "absEta_1_1.4442-gold" ){
    if(Et>=20 && Et<40)ecal +="-Et_20_40";
    else if (Et>=40 && Et<50)ecal +="-Et_40_50";
    else if (Et>=50)ecal +="-Et_50_100";
  }

  return ecal;
}

int photonEnergyScale::EcalPart_EtDep(float r9, float scEta, float Et ) {
  return EcalPart_EtDep( EcalPartString_EtDep( r9, scEta, Et ) );
}

bool photonEnergyScale::setup_EtDep( string energyscaleName ) {
   ifstream input(energyscaleName.c_str());

   //// default energyscale values
   
   
   if(!input.good()){
     std::cerr << "[ERROR] file " << energyscaleName.c_str() << " not readable" << std::endl;
     return false;
   }


   // //  int runMin, runMax;
   // TString category, region2;
   // double deltaP, err_deltaP;
   
   // string ecal;
   // int runMin,runMax;
   // float eScale, err_eScale;
   
   // string dummy;
   // for(input >> category; input.good(); input >> category){
   //   //  input >> region2
   //   //	   >> runMin >> runMax
   //   //	   >> deltaP >> err_deltaP;
     

   //   input >> dummy >> runMin >> runMax >> eScale >> err_eScale ;
   //   cout <<"inside setup ---->ecal: "<< category   <<"  "<< dummy <<" "<< runMin <<" "<< runMax <<" "<<  eScale <<"  "<< err_eScale  <<endl;

   //   //cout << region2 <<"  "
   //   //	  <<       runMin <<"  "<< runMax <<"  "
   //   //	  << deltaP <<" "<< err_deltaP<<endl;
   //   // Add(category, runMin, runMax, deltaP, err_deltaP);
     
   //   int iEcal = EcalPart_EtDep(std::string(category));
   //   cout <<" iEcal "<<iEcal <<endl;
     
   //   if( iEcal < 0 ) continue;
   //   _runStart[  iEcal].push_back( runMin );
   //   _runStop[   iEcal].push_back( runMax );
   //   _eScale[    iEcal].push_back( eScale );
   //   _err_eScale[iEcal].push_back( err_eScale );

   // }
  
   // input.close();


   cout << " Opening energy scale file: " << energyscaleName << endl;
   while ( input.good() && !input.eof() ) {
     string line;
     getline(input,line,'\n');
    
     
     /// comment
     if( line.find("#") != string::npos ) continue;
     
     //  istringstream isstream1(line);
     istringstream isstream2(line);
     string ecal;
     int runMin,runMax;
     float eScale, err_eScale;
     
     string dummy;
     /// most likely we are reading the new Format
     isstream2 >> ecal >> dummy >> runMin >> runMax >> eScale >> err_eScale ;
     
     
     int iEcal = EcalPart_EtDep(ecal);
     
     
     if( iEcal < 0 ) continue;
     _runStart[  iEcal].push_back( runMin );
     _runStop[   iEcal].push_back( runMax );
     _eScale[    iEcal].push_back( eScale );
     _err_eScale[iEcal].push_back( err_eScale );

     //   cout <<"  _runStart[  iEcal]  "<<   _runStart[  iEcal] <<"  _runStop[   iEcal] "<< _runStop[   iEcal] << "  _eScale[ iEcal] "<< _eScale[iEcal] <<endl;

   }
   _isSetup = true;
   return true;
}

