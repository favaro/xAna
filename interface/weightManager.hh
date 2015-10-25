#ifndef weightmanager_hh__
#define weightmanager_hh__

#include <string>
#include "interface/higgsCrossSections.hh"

class TH1F;

class WeightManager {
 public:
   WeightManager(std::string filename, std::string filelist, int LHC_sTeV = 8 );
  
  /// xsec
  /// xsecType == 0 : normal mode
  /// xsecType == 1 : WH in case of mixed WHZH MC
  /// xsecType == 2 : ZH in case of mixed WHZH MC
  /// this can be adapted to any type of mixed MC signal
   float xSecW(int xsecType = 0);
   CrossSection* getCrossSection(void) {return &_xsec;}
   void setCrossSection(const CrossSection &xsec);
   void setLumi(float);
   void setPUTarget(std::string);
   void setRDPUMC(TH1F*);
  /// Primary Vertex ReWeighting
   float puW(int npv); 
   float pvzW(float pvz); 
   float bszW(float diff_pvz);
  
  // vertex correction
   float vtxPtCorrW( float hPt );

   float getNevts(void) const {return _nEvtTot;}

  /// Z Pt ReWeighting
   float pTzW(float pTz, bool isZmc ); 

  /// photon identification
   float phoIdPresel( float r9, float sceta );
   float phoIdMVA(    float r9, float sceta );
   float phoIdCiC(    float r9, float sceta );
   float phoIdEleVeto(float r9, float sceta );

  /// electron identification
   float elecId( float sceta );

 private:
  std::string _weight_file;
  std::string _pu_data_file;
  std::string _vtxCorr_weight_file;
  void SetPUandNevt(  std::string filename, std::string filelist );

  float _nEvtTot;
  float _lumi;
  CrossSection _xsec;
  TH1F *hnPV, *hPVz;
  TH1F *hVtxId;
  TH1F *hPU;
  TH1F *hpTz;

  int iLHC_sTeV;
  float photonIdCorrMVA[4][3];
  float photonIdCorrPresel[4][3];
  float photonIdCorrCiC[4][3];
  float photonIdCorrEleVeto[4][3];

  float electronIdCorr[4][3];
  
};

#endif
