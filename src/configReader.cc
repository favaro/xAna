#include "interface/configReader.hh"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
using namespace std;

#include <boost/program_options.hpp>

ConfigReader::ConfigReader( int nargc, char **argv ) {
  //// default config values
  /// analysis: baseline / MVA / ...
  _analysisType = "baseline";
  /// setup: Legacy2013_7TeV / prompt2012_ichep /...
  _setup = "Legacy2013_7TeV";
  
  _diJetPtg1M = 1./3.;
  _diJetPtg2M = 1./4.;
  _diJetMjj   = 250;
  _diJetPtg1MCiC = 1./2.; //JM
  
  _outDir = "./output_test_baseline/";
  _se = "unknown";
  _runEleGamma = 0;
  _runGammaEle = 0;
  _invertElectronVeto = 0;
  _pt1Cut = 25;
  _pt2Cut = 25;
  _loosePt1MCut = 1./3.;
  _loosePt2MCut = 1./4.;
  _VHPt1MCut = 45./120.; 
  _VHhadPt1MCut = 60./120.; 
  _TTHPt1MCut = 1./2.;   
  
  _mggCut = 90;
  
  _nEvtsPerJob = 500000;
  _systEnergyResolution = 0;
  _doPUJetID = 0;
  _noElectronVeto = 0;
  
  _phoIdMvaCut = -0.2; /// 2012 setup (-0.3 in 2011).
  _doSkimming = false;
  _luminosity = 5089;
  _vbfMvaCat = false;
  _vbfCombMvaCat = false;
  _energyScaleCorrectionFile = "unknown";
  
  namespace po = boost::program_options;
  po::options_description config("Configuration");    
    

  config.add_options()
    ("energyScaleCorrection", po::value<string>(&_energyScaleCorrectionFile),"ecal scale correction file")
    ("VBFmvaCategory"       , po::value<bool>(&_vbfMvaCat), "VBF cat" )
    ("VBFCombmvaCategory"   , po::value<bool>(&_vbfCombMvaCat), "VBF comb cat" )
    ("doPUJetID"            , po::value<bool>(&_doPUJetID), "PU resilient jet ID" )
    ("pt1_cut"  , po::value<float>(&_pt1Cut), "pt cut leading" )
    ("pt2_cut"  , po::value<float>(&_pt2Cut), "pt cut trailing")
    ("mgg_cut"  , po::value<float>(&_mggCut), "mgg cut")
    ("phoIdMva_cut"  , po::value<float>(&_phoIdMvaCut ),  "loose mva pho ID cut in mva Ana" )
    ("loosePt1M_cut" , po::value<float>(&_loosePt1MCut), "loose pt/m cut leading" )
    ("loosePt2M_cut" , po::value<float>(&_loosePt2MCut), "loose pt/m cut trailing")
    ("VHPt1M_cut"    , po::value<float>(&_VHPt1MCut   ), "pt/m leading for VH"    )
    ("VHhadPt1M_cut" , po::value<float>(&_VHhadPt1MCut), "pt/m leading for VHhad" )
    ("TTHPt1M_cut"   , po::value<float>(&_TTHPt1MCut  ), "pt/m leading for tth"   )
    ("diJetPtg1M_cut", po::value<float>(&_diJetPtg1M  ), "pt/m leading  gamma for VBF")
    ("diJetPtg2M_cut", po::value<float>(&_diJetPtg2M  ), "pt/m trailing gamma for VBF")
    ("diJetMjj_cut"  , po::value<float>(&_diJetMjj    ), "Mjj cut")
    ("doSkimming"    , po::value<bool>(  &_doSkimming)->default_value(false),"skimming")
    ("storageElement", po::value<string>(&_se), "Storage element (skimming only)")
    ("invEleVeto"    , po::value<bool>(&_invertElectronVeto),"invert electron veto")
    ("noElectronVeto", po::value<bool>(&_noElectronVeto    ),"remove electron veto")
    ("runEleGamma"   , po::value<bool>(&_runEleGamma       ),"e-gamma bkg study" )
    ("runGammaEle"   , po::value<bool>(&_runGammaEle       ),"gamma-e bkg study" )
    ("systEnergyResolution", po::value<int>(&_systEnergyResolution)->default_value(0),"syst on energy resolution -1/0/+1")
    ;

  po::options_description generic("Generic options");
  generic.add_options()
    ("help,h", "produce help message")
    ("inputFile,i"   ,po::value<string>(&_inputFile    )->default_value("unknown"), "input file")
    ("iJob"          ,po::value<int>   (&_iJob         )->default_value(0), "job number (when splitting a single input file)")
    ("analysisType,a",po::value<string>(&_analysisType )->default_value("MVA"),"MVA/CiC" )
    ("setup,s"       ,po::value<string>(&_setup        )->default_value("Legacy2013_8TeV"),"setup" )
    ("outputDir,o"   ,po::value<string>(&_outDir       )->default_value("tmp/"),"output directory")
    ("nEvtsPerJob,n" ,po::value<int>   (&_nEvtsPerJob  )->default_value(-1),"#evts / job")
    ("luminosity"    ,po::value<float> (&_luminosity   )->default_value(1),"luminosity")
    ("categoryDefDir,c"   ,po::value<string>(&_categoryDefDir   )->default_value("notdefined"),"category definition")
    ;
  
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("config-file", po::value<string>(&_configName), "config file")
    ;
  
  po::positional_options_description p;
  p.add("config-file", -1);
  
  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config).add(hidden);
  
  po::options_description config_file_options;
  config_file_options.add(config).add(generic);
  
  po::options_description visible("Allowed options");
  visible.add(generic).add(config);
  
  po::variables_map vm;
  po::store(po::command_line_parser(nargc, argv).
	    options(cmdline_options).positional(p).run(), vm);
  po::notify(vm);
  
  _go = true;
  if( vm.count("help") ) {
    cout << visible << endl;
    cout << "Usage: " << endl
	 << "xAnaExe [options] configFile" << endl;
    _go = false;
    return;
  }
  
  
  ifstream ifs(_configName.c_str());
  bool setupOK = true;
  if (!ifs) {
    cout << "can not open config file: " << _configName << endl;
    setupOK = false;
  } else {
    po::store(parse_config_file(ifs, config_file_options), vm);
    notify(vm);
  }
  cout << " Opening config file: " << _configName << endl;
  
  
  if( !setupOK ) _go = false;
  
  
  if( ! ( _setup == "Legacy2013_7TeV" || _setup == "Prompt2012_ichep" || _setup == "Legacy2013_8TeV")  ) 
    throw std::runtime_error("setup can only be: \n" +
			     string("    - Legacy2013_7TeV \n"  ) +
			     string("    - Legacy2013_8TeV \n"  ) );  

  
  string mkdir = "mkdir -p " + _outDir;
  gSystem->Exec( mkdir.c_str() );
  printConfig();
}

void ConfigReader::printConfig(void) const {
  ///// summary of the config file content...
  cout << " ============= Config file Content ================" << endl;
  cout << " -------- config file   : " << _configName << endl;
  cout << "    - inputFile         : " << _inputFile << endl; 
  cout << "    - iJob              : " << _iJob << endl; 
  cout << "    - analysis type     : " << _analysisType << endl; 
  cout << "    - category def      : " << _categoryDefDir << endl; 
  cout << "    - setup             : " << _setup << endl; 
  cout << "    - output directory  : " << _outDir << endl; 
  cout << "    - invert elec veto  : " << _invertElectronVeto << endl; 
  cout << "    - no elec veto      : " << _noElectronVeto << " (do not use for real analysis)" <<endl; 
  cout << "    - pt1 cut           : " << _pt1Cut << endl;
  cout << "    - pt2 cut           : " << _pt2Cut << endl;
  cout << "    - loose pt1/m cut   : " << _loosePt1MCut << endl;
  cout << "    - loose pt2/m cut   : " << _loosePt2MCut << endl;
  cout << "    - mgg cut           : " << _mggCut << endl;
  cout << "    - lumiosity         : " << _luminosity << endl;
  cout << "    - #evts per job     : " << _nEvtsPerJob << endl;
  cout << "    - oversmearing syst : " << _systEnergyResolution << endl;
  cout << "    - Jet PU ID cut     : " << _doPUJetID  << endl;   
  cout << "    - DiJet Tree pt1/M  : " << _diJetPtg1M << endl;
  cout << "    - DiJet Tree pt2/M  : " << _diJetPtg2M << endl;
  cout << "    - VH pt1/m cut      : " << _VHPt1MCut << endl; 
  cout << "    - VH had pt1/m cut  : " << _VHhadPt1MCut << endl; 
  cout << "    - TTH pt1/m cut     : " << _TTHPt1MCut << endl;
  cout << "    - DiJet Tree Mjj    : " << _diJetMjj   << endl;
  cout << "    - skimming          : " << _doSkimming << endl;
  cout << "    - storage dir (eos) : " << _se << "  (used for skimming only)" << endl;
  cout << " ==================================================" << endl;
}



ConfigReader::ConfigReader( string myConfigName ) {
  _configName = myConfigName;
   ifstream input(_configName.c_str());

   //// default config values
   /// analysis: baseline / MVA / ...
   _analysisType = "baseline";
   /// setup: Legacy2013_7TeV / prompt2012_ichep /...
   _setup = "Legacy2013_7TeV";
   
   _diJetPtg1M = 1./3.;
   _diJetPtg2M = 1./4.;
   _diJetMjj   = 250;
   _diJetPtg1MCiC = 1./2.; //JM

   _outDir = "./output_test_baseline/";
   _se = "unknown";
   _runEleGamma = 0;
   _runGammaEle = 0;
   _invertElectronVeto = 0;
   _pt1Cut = 25;
   _pt2Cut = 25;
   _loosePt1MCut = 1./3.;
   _loosePt2MCut = 1./4.;
   _VHPt1MCut = 45./120.; 
   _VHhadPt1MCut = 60./120.; 
   _TTHPt1MCut = 1./2.;   

   _mggCut = 90;

   _nEvtsPerJob = 500000;
   _systEnergyResolution = 0;
   _doPUJetID = 0;
   _noElectronVeto = 0;

   _phoIdMvaCut = -0.2; /// 2012 setup (-0.3 in 2011).
   _doSkimming = false;
   _luminosity = 5089;
   _vbfMvaCat = false;
   _vbfCombMvaCat = false;
   _energyScaleCorrectionFile = "unknown";
   cout << " Opening config file: " << _configName<< endl;
   while ( input.good() && !input.eof() ) {
     string line;
     getline(input,line,'\n');
     
     /// comment
     if( line.find("#") != string::npos ) continue;
     
     int posSeparator = line.find("=");
     
     if( posSeparator < 0 ) continue;
     
     /// split the line in key : val
     string key = line.substr(0,posSeparator);
     string val = line.substr(posSeparator+1,line.size());
     removeSpaces(key);
     removeSpaces(val);
     
     if( key.find( "analysisType" ) != string::npos ) _analysisType = val; 
     
     if( key.find("pt1_cut" ) != string::npos ) _pt1Cut = atof(val.c_str());
     if( key.find("pt2_cut" ) != string::npos ) _pt2Cut = atof(val.c_str());
     if( key.find("loosePt1M_cut") != string::npos ) _loosePt1MCut = atof(val.c_str());
     if( key.find("loosePt2M_cut") != string::npos ) _loosePt2MCut = atof(val.c_str());
     if( key.find("VHPt1M_cut") != string::npos ) _VHPt1MCut = atof(val.c_str()); 
     if( key.find("VHhadPt1M_cut") != string::npos ) _VHhadPt1MCut = atof(val.c_str()); 
     if( key.find("TTHPt1M_cut") != string::npos ) _TTHPt1MCut = atof(val.c_str()); 
     if( key.find("mgg_cut" ) != string::npos ) _mggCut = atof(val.c_str());
     if( key.find("phoIdMva_cut" ) != string::npos ) _phoIdMvaCut = atof(val.c_str());
     
     if( key.find("luminosity" ) != string::npos ) _luminosity = atof(val.c_str());
     
     if( key.find("outputDir") != string::npos ||
	 key.find("outDir") != string::npos ||
	 key.find("outdir") != string::npos ) _outDir = val;
     
     if( key.find("invertElectronVeto" ) != string::npos ||
	 key.find("invertEleVeto")       != string::npos ||
	 key.find("invEleVeto")          != string::npos ) _invertElectronVeto = atoi(val.c_str());
     
     if( key.find("noElectronVeto" ) != string::npos ) _noElectronVeto = atoi(val.c_str());
     if( key.find("doSkimming"     ) != string::npos ) _doSkimming = atoi(val.c_str());
     if( key.find("storageElement" ) != string::npos ) _se = val;
     
     if( key.find("VBFmvaCategory" ) != string::npos ) _vbfMvaCat = atoi(val.c_str());
     if( key.find("VBFCombmvaCategory" ) != string::npos ) _vbfCombMvaCat = atoi(val.c_str());
     
     if( key.find("nEvtsPerJob") != string::npos ) _nEvtsPerJob = atoi(val.c_str());
     if( key.find("systEnergyResolution") != string::npos ) _systEnergyResolution = atoi(val.c_str());
     if( key.find("doPUJetID") != string::npos ) _doPUJetID = atoi(val.c_str());
     
     if( key.find("runEleGamma") != string::npos ) _runEleGamma = atoi(val.c_str());
     if( key.find("runGammaEle") != string::npos ) _runGammaEle = atoi(val.c_str());
     
     if( key.find("energyScaleCorrection") != string::npos ) _energyScaleCorrectionFile = val;
       
     if( key.find("diJetPtg1M_cut") != string::npos ){
       _diJetPtg1M    = atof(val.c_str());
       _diJetPtg1MCiC    = atof(val.c_str()); 

     }
     if( key.find("diJetPtg2M_cut") != string::npos ) _diJetPtg2M = atof(val.c_str());
     if( key.find("diJetMjj_cut"  ) != string::npos ) _diJetMjj   = atof(val.c_str());
     
     
     if( key.find("setup") != string::npos || key.find("Setup") != string::npos ) {
       /// first verify the setup is allowed
       if( ! ( val == "Legacy2013_7TeV" || val == "Prompt2012_ichep" || val == "Legacy2013_8TeV")  ) throw std::runtime_error("setup can only be: \n" +
															      string("    - Legacy2013_7TeV \n"      ) +
															      string("    - Prompt2012_ichep \n") +
															      string("    - Legacy2013_8TeV \n") );     
       _setup = val;
     }  
   }
   
   // JM: commented because seems dangerous
   //if( _diJetPtg1M < _loosePt1MCut ) _loosePt1MCut = _diJetPtg1M;
   //if( _diJetPtg2M < _loosePt2MCut ) _loosePt2MCut = _diJetPtg2M;
   
   input.close();
   if(_analysisType=="baselineCiC4PF") _diJetPtg1M=_diJetPtg1MCiC; //JM: can't find a better way to do this for now 


   string mkdir = "mkdir -p " + _outDir;
   gSystem->Exec( mkdir.c_str() );
   



   ///// summary of the config file content...
   printConfig();
   
}



ConfigReader::~ConfigReader(void) {}


void removeSpaces(string &str) {
  unsigned long ipos = -1;
  while( (ipos = str.find(" ")) != string::npos ) str.replace(ipos,1,"");
}

void testConfigFileReader(string configFile ) {
  ConfigReader config( configFile );

}
