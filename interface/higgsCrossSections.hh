#ifndef higgsXsection_hh_
#define higgsXsection_hh_

#include <map>
#include <string>
#include <fstream>
#include <sstream>

std::map<float,float> ReadHiggsCrossSecionFile(    std::string file);
std::map<float,float> ReadHiggsBranchingRatioFile( std::string file);

class CrossSection {
public:
  /// xsec2 is usually not relevant
  /// nevertheless for WH/ZH mixed signal: xsec1=WH ; xsec2=ZH
  
  CrossSection( void );
  CrossSection( float xsec1, float xsec2, float br,
		int mode, float higgMass, std::string mcType );
  CrossSection( const CrossSection &xsec );
  
  float getXsec(void)      const { return _xsec1; }
  float getBR(void)        const { return _br; }
  float getHiggsMass(void) const { return _higgsMass;}
  float getWHxsec(void)    const { return _isMixedWHZH ? _xsec1: -1;}
  float getZHxsec(void)    const { return _isMixedWHZH ? _xsec2: -1;}
  int   mode(void)         const { return _mode; }
  std::string getMCType(void)   const { return _mcType; }
  CrossSection& operator=(const CrossSection &xsec);
  
private:
  float _xsec1, _xsec2;
  float _br;
  int   _mode;
  float _higgsMass;
  bool _isMixedWHZH;
  std::string _mcType;
  
  float _vtxID_sf;
};


class SMHiggsCrossSection {
public:
  SMHiggsCrossSection( void );
  ~SMHiggsCrossSection(void);
  
  void is7TeV( void ) { _xsec_at_nTeV =  7; }
  void is8TeV( void ) { _xsec_at_nTeV =  8; }
  void is14TeV(void ) { _xsec_at_nTeV = 14; }
  
  float HiggsSMxsec_ggh( float mass );
  float HiggsSMxsec_vbf( float mass );
  float HiggsSMxsec_wh(  float mass );
  float HiggsSMxsec_zh(  float mass );
  float HiggsSMxsec_tth( float mass );
  float HiggsBR(  float mass );
  
private:
  int     _xsec_at_nTeV;
  std::string _file_7TeV_ggh;
  std::string _file_8TeV_ggh;
  std::string _file_14TeV_ggh;
  std::string _file_7TeV_vbf;
  std::string _file_8TeV_vbf;
  std::string _file_14TeV_vbf;
  std::string _file_7TeV_wh;
  std::string _file_8TeV_wh;
  std::string _file_14TeV_wh;
  std::string _file_7TeV_zh;
  std::string _file_8TeV_zh;
  std::string _file_14TeV_zh;
  std::string _file_7TeV_tth;
  std::string _file_8TeV_tth;
  std::string _file_14TeV_tth;
  
  std::string _file_br;
  
  std::map<float,float> _xsec_ggh;
  std::map<float,float> _xsec_vbf;
  std::map<float,float> _xsec_wh;
  std::map<float,float> _xsec_zh;
  std::map<float,float> _xsec_tth;
  std::map<float,float> _br;
};



#endif
