#include "interface/xAna.hh"
#include "interface/higgsCrossSections.hh"
#include "interface/weightManager.hh"
#include "interface/setupReader.hh"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>

#include <iostream>
#include <stdio.h>
#include <string>
using namespace std;


CrossSection FindXsec( string inputfile, int enLHC);

int main( int nargc, char **argv ) {
  int ijob = 0;
  string fileListForNorm = "unknown";

  ConfigReader config(nargc,argv);
  if( !config.go() ) return 1;

  TString jobName       = config.inputfile().c_str(); 
  jobName.Remove(0,jobName.Last('/')+1);  /// this is the fileName
  TString inputFileName = config.inputfile().c_str(); 

  TString outDir        = config.outDir().c_str();
  TString outFileBase   = jobName + "_" + itostr(config.iJob()) + ".root";
  TString outFileName   = outDir + outFileBase;

  cout << " ----- overwrite output directory ---------  " << endl;
  cout << " ----- new dir (from cmd line): " << outDir << endl;
  gSystem->Exec( "mkdir -p " + outDir );
    
  /// FC: add xsec finder for MC (Data xsec = -1 by default)
  int enLHC = 2; /// ==> 14TeV!!!
  if( config.setup().find("2011") != string::npos ) enLHC = 0;
  if( config.setup().find("2012") != string::npos ) enLHC = 1;
  if( config.setup().find("2013") != string::npos ) enLHC = 1;
  CrossSection xsec = FindXsec( config.inputfile(), enLHC );
  Double_t lumi = config.luminosity();

  /// FC: change the logic to manage weights
  WeightManager weight_manager( config.inputfile(), fileListForNorm);
  weight_manager.setCrossSection(xsec);
  weight_manager.setLumi(lumi);
  Int_t   nEvtsGen = weight_manager.getNevts();
 
  //JM: to disable MC branches when running on data 
  bool isItData=false;
  int mode = weight_manager.getCrossSection()->mode();
  if(mode==-1) isItData = true;     

  TFile* f = TFile::Open(inputFileName);   
  cout << "Opened file : " << inputFileName << endl;
  
  f->cd("ggNtuplizer");    
  TTree* tree = (TTree*) gDirectory->Get("EventTree");
  cout <<"ggNtuplizer opened  "<<endl;

  xAna* ggAna = new xAna(tree,isItData); 
  ggAna->SetWeightManager(&weight_manager);

  ggAna->SetConfigReader(config);

  cout << " --- Cross Section (pb): " << xsec.getXsec() << endl;
  cout << " --- Lumi        (pb-1): " << lumi << endl;
  cout << " --- Nevts             : " << nEvtsGen << endl;
  cout <<"     ====== MODE = "<< xsec.mode()  <<endl;

  //cout << " *******************************" <<  weight_manager.getCrossSection()->getMCType() << endl;
  
  int nevtpass = ggAna->HggTreeWriteLoop(outFileName, config.iJob());

  cout << "============= job ending ===========" << endl;
  cout <<"outFileName:  "<< outFileName <<endl; 
  if( nevtpass > 0 && config.doSkimming() ) {
    TH1F* h[3]= { 0, 0, 0 };
    h[0] = (TH1F*) TFile::Open(inputFileName)->Get("ggNtuplizer/hEvents");   
    h[1] = (TH1F*) TFile::Open(inputFileName)->Get("ggNtuplizer/hPU");   
    h[2] = (TH1F*) TFile::Open(inputFileName)->Get("ggNtuplizer/hPUTrue");   
    
    TFile *fout = TFile::Open(outFileName,"update");
    fout->cd();
    for( int ih = 0 ; ih < 3; ih++ ) 
      if( h[ih] != 0 ) h[ih]->Write();
    fout->Close();
    
    //// store skimmed file of SE
    string dirSE = config.storageElementDir();
    if( dirSE != "unknown" ) {
      cout << " --- trying to store on storage element (must be on eos!!!): " << dirSE << endl;
      TString mkdir  = "cmsMkdir " + TString( dirSE.c_str() );
      TString rmfile = "cmsRm "    + TString( dirSE.c_str() ) + outFileBase;
      TString mvfile = "cmsStage " + outFileName + " " + TString( dirSE.c_str() );
      TString rmloc  = "rm -f " + outFileName;
      cout << mkdir  << endl; gSystem->Exec( mkdir  );
      cout << rmfile << endl; gSystem->Exec( rmfile );
      cout << mvfile << endl; int out = gSystem->Exec( mvfile );
      cout << "  output com = " << out << endl;
      if( out == 0 ) {
	/// rm file only is stored command exit properly
	cout << rmloc  << endl; gSystem->Exec( rmloc  );
      }
    }
  }

  delete ggAna;
  return 0;
}



CrossSection FindXsec( string inputfile, int enLHC  ) {
  /// findout the sigmaxkFactor from the file name

  Double_t xsec    = -1;
  Double_t xsec_zh = -1;
  int mode = -1;
  float higgsMass = -1;
  float br = 1;
  string sigType = "background";
  if( enLHC == 0 ) {
    if( inputfile.find("gjet_pt20_doubleEM"   ) != string::npos ) { xsec = 501.15 *1.3; mode = 2;}
    if( inputfile.find("diphoton_box_10to25"  ) != string::npos ) { xsec = 358.2  *1.3; mode = 5;}
    if( inputfile.find("diphoton_box_25to250" ) != string::npos ) { xsec = 12.37  *1.3; mode = 5;}
    if( inputfile.find("diphoton_box_250toInf") != string::npos ) { xsec = 2.08e-4*1.3; mode = 5;}
    if( inputfile.find("DYJetsToLL")            != string::npos ) { xsec = 2475  *1.15; mode = 20;}
    if( inputfile.find("diphoton_jets")         != string::npos ) { xsec = 154.7 *1.15; mode = 4;} 
    
    /// some subtleties for kFactor depending if there are prompt photons
    /// for now kFactor is 1 in both cases
    if( inputfile.find("qcd_pt30to40")          != string::npos ) { xsec = 10868; mode = 19;}
    if( inputfile.find("qcd_pt40")              != string::npos ) { xsec = 43571; mode = 19;}
    
    /// Higgs signal file (mode is not yet setup properly)
    if( inputfile.find("madgraph_ggF_012jet")   != string::npos ) { xsec = 1; mode = 10004;} 
  } else if( enLHC == 1 )  {
    if( inputfile.find("gjet_pt20to40_doubleEM") != string::npos ) { xsec = 8.19*1.84*10 *1.3; mode = 2;} ///10 = 10e4*10e-3
    if( inputfile.find("gjet_pt20to40_doubleEM") != string::npos ) { xsec = 8.84*5.39*10 *1.3; mode = 2;} ///1- = 10e3*10e-2
    if( inputfile.find("diphoton_box_10to25"   ) != string::npos ) { xsec = 424.8  *1.3; mode = 5;}
    if( inputfile.find("diphoton_box_25to250"  ) != string::npos ) { xsec = 15.54  *1.3; mode = 5;}
    if( inputfile.find("diphoton_box_250toInf" ) != string::npos ) { xsec = 1.18e-3*1.3; mode = 5;}
    if( inputfile.find("DYJetsToLL")             != string::npos ) { xsec = 3530; mode = 20;} /// already NNLO
    if( inputfile.find("diphoton_jets")          != string::npos ) { xsec = 75.4 *1.15; mode = 4;} 
    if( inputfile.find("diphotonjets")           != string::npos ) { xsec = 75.4 *1.15; mode = 4;} 
    
    if( inputfile.find("DiPhotonJetsBox_M60_8TeV-sherpa" ) != string::npos ) { xsec = 122.1 *1.20; mode = 4;} 
    
    /// some subtleties for kFactor depending if there are prompt photons
    /// for now kFactor is 1 in both cases
    if( inputfile.find("qcd_pt30to40")          != string::npos ) { xsec = 5.20*2.35*1000 ; mode = 19;}
    if( inputfile.find("qcd_pt40")              != string::npos ) { xsec = 2.37*2.18*10000; mode = 19;}
    
    /// Higgs signal file (mode is not yet setup properly)
    if( inputfile.find("madgraph_ggF_012jet")   != string::npos ) { xsec = 75.4; mode = 10004;} 
  }

  if( inputfile.find("hgg")   != string::npos ) sigType = "ggh";
  if( inputfile.find("ggH")   != string::npos ) sigType = "ggh";
  if( inputfile.find("VBFH")  != string::npos ) sigType = "vbf";
  if( inputfile.find("VBF_")  != string::npos ) sigType = "vbf_";
  if( inputfile.find("WH_ZH") != string::npos ) sigType = "vh";
  if( inputfile.find("TTH")   != string::npos ) sigType = "tth";

  if( sigType != "background" ) {
    SMHiggsCrossSection SMhiggs_xsec; 
    if( enLHC == 0 ) SMhiggs_xsec.is7TeV();
    if( enLHC == 1 ) SMhiggs_xsec.is8TeV();
    if( enLHC == 2 ) SMhiggs_xsec.is14TeV();

    int p1 = inputfile.find("H_")   + 2;
    int p2 = inputfile.find("root") - 1;
    if( sigType == "vh" ) p1 = inputfile.find("ZH_") + 3;
    if( sigType == "vbf_" ) {p1 = inputfile.find("VBF_") + 4; sigType = "vbf";}
    string mass_str = inputfile.substr(p1,p2-p1);
    higgsMass = atof(mass_str.c_str());    
    br = SMhiggs_xsec.HiggsBR( higgsMass );
    if( sigType == "ggh" ) xsec = SMhiggs_xsec.HiggsSMxsec_ggh( higgsMass );
    if( sigType == "vbf" ) xsec = SMhiggs_xsec.HiggsSMxsec_vbf( higgsMass );
    if( sigType == "tth" ) xsec = SMhiggs_xsec.HiggsSMxsec_tth( higgsMass );
    if( sigType == "vh"  ) {
      xsec    = SMhiggs_xsec.HiggsSMxsec_wh( higgsMass );
      xsec_zh = SMhiggs_xsec.HiggsSMxsec_zh( higgsMass );
    }
    if( sigType == "ggh" ) mode = 11000;
    if( sigType == "vbf" ) mode = 555;
    if( sigType == "vh"  ) mode = 1300;
    if( sigType == "tth" ) mode = 20000; // FC: I don't know the mode???

  }
  
  CrossSection crossSectionMC( xsec, xsec_zh, br, mode, higgsMass, sigType );
  return crossSectionMC;
}

