#ifndef configReader_hh_
#define configReader_hh_

#include <TSystem.h>
#include <TLorentzVector.h>


#include <string>


void removeSpaces(std::string& str);
class ConfigReader {
 public:
  ConfigReader( int nargc, char **argv );
  ConfigReader( std::string configName );
  ~ConfigReader(void);

  /// accessors
  inline bool go(void) const { return _go;}
  
  inline std::string inputfile(void) const {return _inputFile;}
  inline int iJob(void) const {return _iJob;}

  inline std::string setup(void)          const {return _setup;}
  inline std::string outDir(void)         const {return _outDir;}
  inline std::string categoryDefDir(void) const {return _categoryDefDir; }
  inline std::string analysisType(void)   const {return _analysisType;}
  inline std::string configName(void)     const {return _configName;}
  inline std::string energyScaleFile(void) const {return _energyScaleCorrectionFile; }
  inline bool   invertElectronVeto(void) const {return _invertElectronVeto;}
  inline bool   noElectronVeto(void) const {return _noElectronVeto;}
  inline bool   doSkimming(void) const {return _doSkimming;}
  inline bool   doVBFmvaCat(void) const {return _vbfMvaCat; }
  inline bool   doVBFCombmvaCat(void) const {return _vbfCombMvaCat; }
  inline std::string storageElementDir(void) const {return _se;}
  inline float  pt1Cut(void) const {return _pt1Cut;}
  inline float  pt2Cut(void) const {return _pt2Cut;}
  inline float  loosePt1MCut() const { return _loosePt1MCut; }
  inline float  loosePt2MCut() const { return _loosePt2MCut; }
  inline float  VHPt1MCut() const { return _VHPt1MCut; }   //JM
  inline float  VHhadPt1MCut() const { return _VHhadPt1MCut; }   //JM
  inline float  TTHPt1MCut() const { return _TTHPt1MCut; } //JM

  inline float  mggCut(void) const {return _mggCut;}
  inline float  photonIdMvaCut(void) const { return _phoIdMvaCut; }

  inline float  luminosity(void) const {return _luminosity;}

  inline int    nEvtsPerJob(void) const {return _nEvtsPerJob;}
  inline int    getEnergyOverSmearingSyst(void) const {return _systEnergyResolution;}
  inline int    doPUJetID(void) const {return _doPUJetID;}

  //  std::string setupType(void) const {return _setupType;}
  /// dijetTree cuts
  inline float  diJetPtg1M(void) const { return _diJetPtg1M; }
  inline float  diJetPtg2M(void) const { return _diJetPtg2M; }
  inline float  diJetMjj(void)   const { return _diJetMjj;   }

  
  inline bool   isEleGamma(void) const {return _runEleGamma;}
  inline bool   isGammaEle(void) const {return _runGammaEle;}

  void printConfig(void) const;

 private:
  std::string _configName;

  bool _go;
  bool _runLepTag;
  bool _runOnlyLepTag;
  bool _runEleGamma;
  bool _runGammaEle;

  /// VBF category, use Mva?
  bool _vbfMvaCat;
  bool _vbfCombMvaCat;

  /// analysis type: baseline, MVA
  /// this mostly changes the photonID selection
  std::string _analysisType;
  std::string _setup;
  std::string _inputFile;
  std::string _categoryDefDir;
  int _iJob;

  /// EnergyScale file
  std::string _energyScaleCorrectionFile;
  
  /// output dir
  std::string _outDir;
  std::string _se; /// for skimming store on eos

  /// electron veto: Z analysis: invert the electron veto
  bool _invertElectronVeto;

  /// do not cut on electron veto, usefull if on needs to do skimming
  bool _noElectronVeto;
  bool _doSkimming;

  /// in case one want to cut on pt.
  float _pt1Cut;
  float _pt2Cut;
  float _mggCut;
  float _phoIdMvaCut;
  float _loosePt1MCut; 
  float _loosePt2MCut; 
  float _VHPt1MCut; 
  float _VHhadPt1MCut; 
  float _TTHPt1MCut;
  

  /// # of evts per job (allows to split one file in several jobs)
  int _nEvtsPerJob;

  /// luminosity
  float _luminosity;

  /// energy resolution over smearing systematic
  int _systEnergyResolution;

  /// Jet ID for 2012 high PU
  bool  _doPUJetID;

  /// dijet cuts to fill dijet tree
  float _diJetPtg1M,_diJetPtg2M,_diJetMjj;
  float _diJetPtg1MCiC;//JM
  
};

#endif
