#ifndef vertexUtils_hh__
#define vertexUtils_hh__

#include <vector>
#include <string>

class ggVertexInfo : public VertexInfoAdapter {
private:
  int _vtx_n;
  std::vector<float> _vtx_x,_vtx_y,_vtx_z;
  std::vector<int>   _vtx_ntks;
  int _tk_n;
  std::vector<std::vector<unsigned short> > _vtx_tkind;
  std::vector<std::vector<float> > _vtx_tkweight;
  std::vector<float> _tk_px, _tk_py, _tk_pz, _tk_pterr;
  std::vector<int>   _tk_quality;

public:
  inline ggVertexInfo( const std::vector<float> &vtx_x, const std::vector<float> &vtx_y, const std::vector<float> &vtx_z,
		       const std::vector<int> &vtx_ntks,
		       const std::vector<float> &tk_px   , const std::vector<float> &tk_py, const std::vector<float> &tk_pz,
		       const std::vector<float> &tk_pterr, const std::vector<int> &tk_quality,
		       std::vector<std::vector<unsigned short> > & vtx_tkind,
		       std::vector<std::vector<float> >    & vtx_tkweight
		       )  {
    _vtx_x = vtx_x;
    _vtx_y = vtx_y;
    _vtx_z = vtx_z;
    _vtx_ntks = vtx_ntks;
    
    _tk_px = tk_px;
    _tk_py = tk_py;
   _tk_pz = tk_pz;
   _tk_pterr   = tk_pterr;
   _tk_quality = tk_quality;
   
   _vtx_tkind    = vtx_tkind;
   _vtx_tkweight = vtx_tkweight;
   
   
   _tk_n  = int(tk_px.size());
   _vtx_n = int(vtx_x.size());
  }
  
  
  virtual int nvtx() const    { return _vtx_n; };
  virtual int ntracks() const { return _tk_n; };
  
  virtual bool hasVtxTracks() const { return true; }
  virtual const unsigned short * vtxTracks(int ii) const { return &(_vtx_tkind)[ii][0]; };
  virtual int vtxNTracks(int ii) const { return _vtx_ntks[ii]; };
  virtual const float * vtxTkWeights(int ii) const { return &(_vtx_tkweight)[ii][0]; };
  
  virtual float tkpx(int ii) const { return _tk_px[ii]; };
  virtual float tkpy(int ii) const { return _tk_py[ii]; };
  virtual float tkpz(int ii) const { return _tk_pz[ii]; };
  
  virtual float tkPtErr(int ii) const { return _tk_pterr[ii]; };
  virtual int   tkVtxId(int ii) const { return -1; };
  
  virtual float tkWeight(int ii, int jj) const { return _vtx_tkweight[jj][ii]; };
  
  virtual float vtxx(int ii) const { return _vtx_x[ii]; };
  virtual float vtxy(int ii) const { return _vtx_y[ii]; };
  virtual float vtxz(int ii) const { return _vtx_z[ii]; };
  
  virtual float tkd0(int ii, int jj) const { return 0.; }; // FIXME
  virtual float tkd0Err(int ii, int jj) const { return 1.; };  // FIXME
  
  virtual float tkdz(int ii, int jj) const { return 0.; };  // FIXME
  virtual float tkdzErr(int ii, int jj) const { return 1.; };  // FIXME
  
  virtual bool tkIsHighPurity(int ii) const { return ( _tk_quality[ii] & (1<<2) ) >> 2; };
  
  // virtual ~ggVertexInfo();
  
};



//---- function completely defined inline for commodity
#include "interface/setupReader.hh"
#include "VertexAnalysis/VertexAlgoParameters.h"
inline VertexAlgoParameters VertexAlgoParametersReader( std::string file ){
  VertexAlgoParameters params;
  SetupReader reader( file );
  bool btmp(false);
  int   itmp;
  float ftmp;
  std::string stmp;
  
  itmp = reader.getInt( "fixTkIndex"        , btmp); if( btmp ) params.fixTkIndex         = bool(itmp);
  itmp = reader.getInt( "rescaleTkPtByError", btmp); if( btmp ) params.rescaleTkPtByError = bool(itmp);
  itmp = reader.getInt( "highPurityOnly"    , btmp); if( btmp ) params.highPurityOnly     = bool(itmp);
  itmp = reader.getInt( "removeTracksInCone", btmp); if( btmp ) params.removeTracksInCone = bool(itmp);
  itmp = reader.getInt( "useAllConversions" , btmp); if( btmp ) params.useAllConversions  = itmp;
  ftmp = reader.getFloat( "trackCountThr", btmp); if( btmp ) params.trackCountThr = ftmp;
  ftmp = reader.getFloat( "maxD0Signif"  , btmp); if( btmp ) params.maxD0Signif = ftmp;
  ftmp = reader.getFloat( "maxDzSignif"  , btmp); if( btmp ) params.maxDzSignif = ftmp;
  ftmp = reader.getFloat( "coneSize"     , btmp); if( btmp ) params.coneSize = ftmp;
  ftmp = reader.getFloat( "sigma1Pix"    , btmp); if( btmp ) params.sigma1Pix = ftmp;
  ftmp = reader.getFloat( "sigma1Tib"    , btmp); if( btmp ) params.sigma1Tib = ftmp;
  ftmp = reader.getFloat( "sigma1Tob"    , btmp); if( btmp ) params.sigma1Tob = ftmp;
  ftmp = reader.getFloat( "sigma1PixFwd" , btmp); if( btmp ) params.sigma1PixFwd = ftmp;
  ftmp = reader.getFloat( "sigma1Tid"    , btmp); if( btmp ) params.sigma1Tid = ftmp;
  ftmp = reader.getFloat( "sigma1Tec"    , btmp); if( btmp ) params.sigma1Tec = ftmp;
  ftmp = reader.getFloat( "sigma2Pix"    , btmp); if( btmp ) params.sigma2Pix = ftmp;
  ftmp = reader.getFloat( "sigma2Tib"    , btmp); if( btmp ) params.sigma2Tib = ftmp;
  ftmp = reader.getFloat( "sigma2Tob"    , btmp); if( btmp ) params.sigma2Tob = ftmp;
  ftmp = reader.getFloat( "sigma2PixFwd" , btmp); if( btmp ) params.sigma2PixFwd = ftmp;
  ftmp = reader.getFloat( "sigma2Tid"    , btmp); if( btmp ) params.sigma2Tid = ftmp;
  ftmp = reader.getFloat( "sigma2Tec"    , btmp); if( btmp ) params.sigma2Tec = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1Pix"   , btmp); if( btmp ) params.singlelegsigma1Pix = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1Tib"   , btmp); if( btmp ) params.singlelegsigma1Tib = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1Tob"   , btmp); if( btmp ) params.singlelegsigma1Tob = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1PixFwd", btmp); if( btmp ) params.singlelegsigma1PixFwd = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1Tid"   , btmp); if( btmp ) params.singlelegsigma1Tid = ftmp;
  ftmp = reader.getFloat( "singlelegsigma1Tec"   , btmp); if( btmp ) params.singlelegsigma1Tec = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2Pix"   , btmp); if( btmp ) params.singlelegsigma2Pix = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2Tib"   , btmp); if( btmp ) params.singlelegsigma2Tib = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2Tob"   , btmp); if( btmp ) params.singlelegsigma2Tob = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2PixFwd", btmp); if( btmp ) params.singlelegsigma2PixFwd = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2Tid"   , btmp); if( btmp ) params.singlelegsigma2Tid = ftmp;
  ftmp = reader.getFloat( "singlelegsigma2Tec"   , btmp); if( btmp ) params.singlelegsigma2Tec = ftmp;

  stmp = reader.getString("vtxProbFormula", btmp); if( btmp ) params.vtxProbFormula = stmp;

  return params;
}


#endif


