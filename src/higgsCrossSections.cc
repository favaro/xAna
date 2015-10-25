#include "interface/higgsCrossSections.hh"
#include <iostream>
using namespace std;

string xSecDir = "etc/inputs/sigCrossSections/";

CrossSection::CrossSection( void ) {
  _xsec1 = _xsec2 = -1;
  _isMixedWHZH = false;
  _mode = -1;
  _higgsMass = -1;
  _br =1;
  _mcType = "unknown";
}

CrossSection::CrossSection( float xsec1, float xsec2, float br, int thismode, float higgsMass, string mcType ) {
  _xsec1 = xsec1; _xsec2 = xsec2;
  if( xsec2 > 0 ) _isMixedWHZH = true;
  _mode = thismode;
  _higgsMass = higgsMass;
  _br = br;
  _mcType = mcType;
}

CrossSection::CrossSection( const CrossSection &xsec ) {
  _xsec1 = xsec._xsec1;
  _xsec2 = xsec._xsec2;
  _br    = xsec._br;
  _mode  = xsec._mode;
  _higgsMass = xsec._higgsMass;
  _isMixedWHZH = xsec._isMixedWHZH;
  _mcType = xsec._mcType;
}

CrossSection& CrossSection::operator=(const CrossSection &xsec ) {
  _xsec1 = xsec._xsec1;
  _xsec2 = xsec._xsec2;
  _br    = xsec._br;
  _mode  = xsec._mode;
  _higgsMass = xsec._higgsMass;
  _isMixedWHZH = xsec._isMixedWHZH;
  _mcType = xsec._mcType;
  return *this;
}


SMHiggsCrossSection::SMHiggsCrossSection(void) {
  _file_7TeV_ggh  = xSecDir + "higgsCrossSection_7TeV_ggH.txt";
  _file_8TeV_ggh  = xSecDir + "higgsCrossSection_8TeV_ggH.txt";
  _file_14TeV_ggh = xSecDir + "higgsCrossSection_14TeV_ggH.txt";
  _file_7TeV_vbf  = xSecDir + "higgsCrossSection_7TeV_VBF.txt";
  _file_8TeV_vbf  = xSecDir + "higgsCrossSection_8TeV_VBF.txt";
  _file_14TeV_vbf = xSecDir + "higgsCrossSection_14TeV_VBF.txt";
  _file_7TeV_wh   = xSecDir + "higgsCrossSection_7TeV_WH.txt";
  _file_8TeV_wh   = xSecDir + "higgsCrossSection_8TeV_WH.txt";
  _file_14TeV_wh  = xSecDir + "higgsCrossSection_14TeV_WH.txt";
  _file_7TeV_zh   = xSecDir + "higgsCrossSection_7TeV_ZH.txt";
  _file_8TeV_zh   = xSecDir + "higgsCrossSection_8TeV_ZH.txt";
  _file_14TeV_zh  = xSecDir + "higgsCrossSection_14TeV_ZH.txt";
  _file_7TeV_tth  = xSecDir + "higgsCrossSection_7TeV_TTH.txt";
  _file_8TeV_tth  = xSecDir + "higgsCrossSection_8TeV_TTH.txt";
  _file_14TeV_tth = xSecDir + "higgsCrossSection_14TeV_TTH.txt";
  
  _file_br        = xSecDir + "higgsBranchingRatio.txt";

  _xsec_at_nTeV = 7;
}
SMHiggsCrossSection::~SMHiggsCrossSection(void) {};

float SMHiggsCrossSection::HiggsBR( float mass ) {
  if( _br.size() == 0 ) {
    string filein = _file_br;
    _br = ReadHiggsBranchingRatioFile( filein );
  }

  /// now the map should be correctly filled
  if( _br.size() == 0 ) {
    cout << " ========== Higgs SM branching ratios are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _br.find( mass );
  if( massit == _br.end() ) {
    cout << "  --- did not find mass: " << mass << " for branching ratio (return +1)" << endl;
    return +1;
  }

  return massit->second;
}



float SMHiggsCrossSection::HiggsSMxsec_ggh( float mass ) {
  if( _xsec_ggh.size() == 0 ) {
    string filein = "unknown";
    if( _xsec_at_nTeV ==  7 ) filein = _file_7TeV_ggh;
    if( _xsec_at_nTeV ==  8 ) filein = _file_8TeV_ggh;
    if( _xsec_at_nTeV == 14 ) filein = _file_14TeV_ggh;
    _xsec_ggh = ReadHiggsCrossSecionFile( filein );
  }

  /// now the map should be correctly filled
  if( _xsec_ggh.size() == 0 ) {
    cout << " ========== Higgs SM xsec are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _xsec_ggh.find( mass );
  if( massit == _xsec_ggh.end() ) {
    cout << "  --- did not find mass: " << mass << " for process ggh (return -1)" << endl;
    return -1;
  }

  return massit->second;
}



float SMHiggsCrossSection::HiggsSMxsec_vbf( float mass ) {
  if( _xsec_vbf.size() == 0 ) {
    string filein = "unknown";
    if( _xsec_at_nTeV ==  7 ) filein = _file_7TeV_vbf;
    if( _xsec_at_nTeV ==  8 ) filein = _file_8TeV_vbf;
    if( _xsec_at_nTeV == 14 ) filein = _file_14TeV_vbf;
    _xsec_vbf = ReadHiggsCrossSecionFile( filein );
  }

  /// now the map should be correctly filled
  if( _xsec_vbf.size() == 0 ) {
    cout << " ========== Higgs SM xsec are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _xsec_vbf.find( mass );
  if( massit == _xsec_vbf.end() ) {
    cout << "  --- did not find mass: " << mass << " for process vbf (return -1)" << endl;
    return -1;
  }

  return massit->second;
}


float SMHiggsCrossSection::HiggsSMxsec_wh( float mass ) {
  if( _xsec_wh.size() == 0 ) {
    string filein = "unknown";
    if( _xsec_at_nTeV ==  7 ) filein = _file_7TeV_wh;
    if( _xsec_at_nTeV ==  8 ) filein = _file_8TeV_wh;
    if( _xsec_at_nTeV == 14 ) filein = _file_14TeV_wh;
    _xsec_wh = ReadHiggsCrossSecionFile( filein );
  }

  /// now the map should be correctly filled
  if( _xsec_wh.size() == 0 ) {
    cout << " ========== Higgs SM xsec are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _xsec_wh.find( mass );
  if( massit == _xsec_wh.end() ) {
    cout << "  --- did not find mass: " << mass << " for process wh (return -1)" << endl;
    return -1;
  }

  return massit->second;
}


float SMHiggsCrossSection::HiggsSMxsec_zh( float mass ) {
  if( _xsec_zh.size() == 0 ) {
    string filein = "unknown";
    if( _xsec_at_nTeV ==  7 ) filein = _file_7TeV_zh;
    if( _xsec_at_nTeV ==  8 ) filein = _file_8TeV_zh;
    if( _xsec_at_nTeV == 14 ) filein = _file_14TeV_zh;
    _xsec_zh = ReadHiggsCrossSecionFile( filein );
  }

  /// now the map should be correctly filled
  if( _xsec_zh.size() == 0 ) {
    cout << " ========== Higgs SM xsec are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _xsec_zh.find( mass );
  if( massit == _xsec_zh.end() ) {
    cout << "  --- did not find mass: " << mass << " for process zh (return -1)" << endl;
    return -1;
  }

  return massit->second;
}



float SMHiggsCrossSection::HiggsSMxsec_tth( float mass ) {
  if( _xsec_tth.size() == 0 ) {
    string filein = "unknown";
    if( _xsec_at_nTeV ==  7 ) filein = _file_7TeV_tth;
    if( _xsec_at_nTeV ==  8 ) filein = _file_8TeV_tth;
    if( _xsec_at_nTeV == 14 ) filein = _file_14TeV_tth;
    _xsec_tth = ReadHiggsCrossSecionFile( filein );
  }

  /// now the map should be correctly filled
  if( _xsec_tth.size() == 0 ) {
    cout << " ========== Higgs SM xsec are not filled properly... =========" << endl;
    return -1;
  }

  map<float,float>::iterator massit = _xsec_tth.find( mass );
  if( massit == _xsec_tth.end() ) {
    cout << "  --- did not find mass: " << mass << " for process tth (return -1)" << endl;
    return -1;
  }

  return massit->second;
}


map<float,float> ReadHiggsCrossSecionFile(string file) {
  ifstream input(file.c_str());
  map<float,float> higgsCrossSection;

   cout << " --- read higgs cross section file: " << file << endl;
   if( !input.good()) {
     cout << "    ===> can not open file: " << file << endl;
     return higgsCrossSection;
   }

   while ( input.good() && !input.eof() ) {
     string line;
     getline(input,line,'\n');
    
     /// comment
     if( line.find("#") != string::npos ) continue;
     
     /// read line
     float mass, xsec;
     istringstream isstream(line);
     isstream >> mass >> xsec;
     higgsCrossSection.insert( make_pair( mass, xsec ) );
   }

   return higgsCrossSection;
}




map<float,float> ReadHiggsBranchingRatioFile(string file) {
  ifstream input(file.c_str());
  map<float,float> higgsBranchingRatio;

   cout << " --- read higgs cross section file: " << file << endl;
   if( !input.good()) {
     cout << "    ===> can not open file: " << file << endl;
     return higgsBranchingRatio;
   }

   while ( input.good() && !input.eof() ) {
     string line;
     getline(input,line,'\n');
    
     /// comment
     if( line.find("#") != string::npos ) continue;
     
     /// read line
     float mass, br_gluglu, br_gamgam;
     istringstream isstream(line);
     isstream >> mass >> br_gluglu >> br_gamgam;
     higgsBranchingRatio.insert( make_pair( mass, br_gamgam ) );
   }

   return higgsBranchingRatio;
}


void testXsecReader( int eLHC = 8) {
  SMHiggsCrossSection signalXsecReader;
  if( eLHC ==  7 )  signalXsecReader.is7TeV();
  if( eLHC ==  8 )  signalXsecReader.is8TeV();
  if( eLHC == 14 )  signalXsecReader.is14TeV();

  float mass[5] = {100 ,114.5, 120, 125, 140 };
  for( int im = 0; im < 5; im++ ) 
   cout << " - mass = " << mass[im] 
	<< "; ggh : " << signalXsecReader.HiggsSMxsec_ggh( mass[im] )
	<< "; vbf : " << signalXsecReader.HiggsSMxsec_vbf( mass[im] )
	<< "; wh  : " << signalXsecReader.HiggsSMxsec_wh( mass[im] )
	<< "; zh  : " << signalXsecReader.HiggsSMxsec_zh( mass[im] )
	<< "; tth : " << signalXsecReader.HiggsSMxsec_tth( mass[im] )
	<< "  ---> BR: " << signalXsecReader.HiggsBR( mass[im] )
	<< endl;

}

