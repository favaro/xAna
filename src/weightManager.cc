#include "interface/weightManager.hh"

#include <TH1.h>
#include <TFile.h>

#include <iostream>
#include <cmath>

using namespace std;

WeightManager::WeightManager( string filename, string filelist, int LHC_sTeV){

  if( LHC_sTeV ==  7 ) iLHC_sTeV = 0;
  if( LHC_sTeV ==  8 ) iLHC_sTeV = 1;
  if( LHC_sTeV == 14 ) iLHC_sTeV = 2;

  _weight_file = "PU/zPtReweigthing.root";
  if(LHC_sTeV==7) _pu_data_file = "pileup_nov08/nov08_rereco.json.68000.pileup.root"; 
  _pu_data_file = "etc/inputs/pileup/pu2012/Moriond_pileup_190456-208686_minBiasXsec69400_corr_observed.root";

  string vtxId_corrFactor_7TeV  = "???";
  string vtxId_corrFactor_8TeV  = "etc/inputs/pileup/pu2012/vtxIdScaleFactorFromZmumu_PUweights_minBiasXsec69400_observed_Run2012ABC_BSreweight_new.root";
  string vtxId_corrFactor_14TeV = "???";
  if     ( LHC_sTeV == 7 ) _vtxCorr_weight_file = vtxId_corrFactor_8TeV;
  else if( LHC_sTeV == 8 ) _vtxCorr_weight_file = vtxId_corrFactor_8TeV;
  else if( LHC_sTeV == 8 ) _vtxCorr_weight_file = vtxId_corrFactor_14TeV;

  /// vertex identification correction (MVA)
  hVtxId = 0; 
  hVtxId = (TH1F*) TFile::Open( _vtxCorr_weight_file.c_str(),"read")->Get("hscaleFactor");
  if( hVtxId == 0 ) {
    cout << " Vertex identification file: " << _vtxCorr_weight_file << " does not contain hist: hscalefactor" << endl;
    cout << "       ===> most likely I will crash." << endl;
  }
  
  TH1F *h_data = 0;
  TH1F *h_mc   = 0;
  
  hnPV = 0;
  _nEvtTot = -1;
  _lumi = -1;
  SetPUandNevt(  filename,  filelist );
  
  hPVz = 0;
  hpTz = 0;
  h_data = 0; h_mc = 0;
  TFile *ftmp = 0; ftmp = TFile::Open(_weight_file.c_str(),"read");
  if( ftmp != 0 && ftmp->IsOpen() ) {
    h_data = (TH1F* )ftmp->Get("PVz_data");
    h_mc   = (TH1F*) ftmp->Get("PVz_mc");
    
    if( h_data ) {
      hPVz = (TH1F*) h_data->Clone();
      h_mc->Scale(h_data->Integral()/h_mc->Integral());
      hPVz->Divide(h_mc);
    }
    
    h_data = 0; h_mc = 0;
    h_data = (TH1F*) ftmp->Get("pTz_data");
    h_mc   = (TH1F*) ftmp->Get("pTz_mc");
    if( h_data ) {
      hpTz = (TH1F*) h_data->Clone();
      h_mc->Scale(h_data->Integral()/h_mc->Integral());
      hpTz->Divide(h_mc);
    }
  }
  //---------------------- photon identification corrections ----------------------//
  /// 7TeV photon id
  photonIdCorrPresel[0][0] = -999;
  photonIdCorrPresel[1][0] = -999;
  photonIdCorrPresel[2][0] = -999;
  photonIdCorrPresel[3][0] = -999;
  photonIdCorrMVA[0][0] = 1.;
  photonIdCorrMVA[1][0] = 1.;
  photonIdCorrMVA[2][0] = 1.;
  photonIdCorrMVA[3][0] = 1.;
  photonIdCorrCiC[0][0] = -999;
  photonIdCorrCiC[1][0] = -999;
  photonIdCorrCiC[2][0] = -999;
  photonIdCorrCiC[3][0] = -999;
  photonIdCorrEleVeto[0][0] = 1.;
  photonIdCorrEleVeto[1][0] = 1.;
  photonIdCorrEleVeto[2][0] = 1.;
  photonIdCorrEleVeto[3][0] = 1.;
  /// 8 TeV photon id: https://hypernews.cern.ch/HyperNews/CMS/get/higgs2g/995/1.html
  photonIdCorrPresel[0][1] = 0.999;
  photonIdCorrPresel[1][1] = 0.995;
  photonIdCorrPresel[2][1] = 1.005;
  photonIdCorrPresel[3][1] = 1.011;
  photonIdCorrMVA[0][1] = 1.0001;
  photonIdCorrMVA[1][1] = 0.9995;
  photonIdCorrMVA[2][1] = 0.9998;
  photonIdCorrMVA[3][1] = 0.9968;
  photonIdCorrCiC[0][1] = 0.985;
  photonIdCorrCiC[1][1] = 0.968;
  photonIdCorrCiC[2][1] = 1.028;
  photonIdCorrCiC[3][1] = 1.001;
  photonIdCorrEleVeto[0][1] = 0.9984;
  photonIdCorrEleVeto[1][1] = 0.9867;
  photonIdCorrEleVeto[2][1] = 0.9893;
  photonIdCorrEleVeto[3][1] = 0.9639;
  /// 14 TeV photon id
  photonIdCorrPresel[0][2] = -999;
  photonIdCorrPresel[1][2] = -999;
  photonIdCorrPresel[2][2] = -999;
  photonIdCorrPresel[3][2] = -999;
  photonIdCorrMVA[0][2] = 1.;
  photonIdCorrMVA[1][2] = 1.;
  photonIdCorrMVA[2][2] = 1.;
  photonIdCorrMVA[3][2] = 1.;
  photonIdCorrCiC[0][2] = -999;
  photonIdCorrCiC[1][2] = -999;
  photonIdCorrCiC[2][2] = -999;
  photonIdCorrCiC[3][2] = -999;
  photonIdCorrEleVeto[0][2] = 1.;
  photonIdCorrEleVeto[1][2] = 1.;
  photonIdCorrEleVeto[2][2] = 1.;
  photonIdCorrEleVeto[3][2] = 1.;


  //---------------------- electron identification corrections ----------------------//
  /// 7TeV
  electronIdCorr[0][0] = 1;
  electronIdCorr[1][0] = 1;
  electronIdCorr[2][0] = 1;
  electronIdCorr[3][0] = 1;
  /// 8TeV:  https://hypernews.cern.ch/HyperNews/CMS/get/higgs2g/995/1.html
  electronIdCorr[0][1] = 0.984;
  electronIdCorr[1][1] = 0.990;
  electronIdCorr[2][1] = 1.013;
  electronIdCorr[3][1] = 0.993;
  /// 14TeV
  electronIdCorr[0][2] = 1;
  electronIdCorr[1][2] = 1;
  electronIdCorr[2][2] = 1;
  electronIdCorr[3][2] = 1;
}

void  WeightManager::setCrossSection(const CrossSection & x) { _xsec = x; }
void  WeightManager::setLumi(float l) { _lumi = l; }
void  WeightManager::setPUTarget(string target){target=_pu_data_file;}
float WeightManager::puW(int npv) {
  if( hnPV == 0 ) return 1;
  int ib = hnPV->GetXaxis()->FindBin(npv);
  return hnPV->GetBinContent(ib);
}
  
float WeightManager::pvzW(float pvz) {
  if( hPVz == 0 ) return 1;
  int ib = hPVz->GetXaxis()->FindBin(pvz);
  return hPVz->GetBinContent(ib);
}

float WeightManager::bszW(float diff_pvz) {
  /// FC: do not really understand why we reweight vs diff MC gen - MC reco
  ///     but ok...
  if( fabs(diff_pvz) < 0.1 ) return 1;

   float newBSmean1  = 9.9391e-02;
   float newBSmean2  = 1.8902e-01;
   float newBSnorm1  = 5.3210e+00;
   float newBSnorm2  = 4.1813e+01;
   float newBSsigma1 = 9.7530e-01;
   float newBSsigma2 = 7.0811e+00;

   float oldBSmean1  = 7.2055e-02;
   float oldBSmean2  = 4.9986e-01;
   float oldBSnorm1  = 3.5411e+00;
   float oldBSnorm2  = 4.0258e+01;
   float oldBSsigma1 = 7.9678e-01;
   float oldBSsigma2 = 8.5356e+00;


   float newBSgaus1 = newBSnorm1*exp(-0.5*pow((diff_pvz-newBSmean1)/newBSsigma1,2));
   float newBSgaus2 = newBSnorm2*exp(-0.5*pow((diff_pvz-newBSmean2)/newBSsigma2,2));
   float oldBSgaus1 = oldBSnorm1*exp(-0.5*pow((diff_pvz-oldBSmean1)/oldBSsigma1,2));
   float oldBSgaus2 = oldBSnorm2*exp(-0.5*pow((diff_pvz-oldBSmean2)/oldBSsigma2,2));
   return 1.1235* (newBSgaus1+newBSgaus2)/(oldBSgaus1+oldBSgaus2);
}


float WeightManager::pTzW(float pTz, bool isZmc ) {
  if( !isZmc    ) return 1;
  if( hpTz == 0 ) return 1;
  int ib = hpTz->GetXaxis()->FindBin(pTz);
  return hpTz->GetBinContent(ib);
}

float WeightManager::xSecW( int xsecType ) {
  if( _lumi < 0 ) return 1;

  float xsec = _xsec.getXsec();
  if( xsecType == 1 )  xsec = _xsec.getWHxsec();
  if( xsecType == 2 )  xsec = _xsec.getZHxsec();
  
  /// remove total number of events from the weight to change easily signal cross section
  return xsec*_xsec.getBR() *_lumi;
}


/// vtx correction.
float WeightManager::vtxPtCorrW( float hPt ) {
  return hVtxId->GetBinContent( hVtxId->GetXaxis()->FindBin( hPt ) );

}

void WeightManager::setRDPUMC(TH1F*hPU_){
      hPU=0;
      hPU=hPU_;
}
void WeightManager::SetPUandNevt(  string filename, string filelist ) {
  TDirectory *mydir = gDirectory;
  TFile *f = TFile::Open(filename.c_str(),"read");

  TH1F *hEvents = 0;
  f->cd("ggNtuplizer"); /// in case ggNtuplizer does not exit: cd 
  hEvents = (TH1F*) gDirectory->Get("hskim");
  if( hEvents == 0 )
    hEvents = (TH1F*)  gDirectory->Get("hEvents");
  hPU = 0;
//  TH1F *hPU = new TH1F("hPU",     "number of pileup",  100, 0, 100);
//  TTree*EventTree=(TTree*)gDirectory->Get("EventTree");
//  EventTree->Draw("nPU>>hPU");
  //hPU = (TH1F*) gDirectory->Get("hPU");
  TFile *froot = 0; 
  if( filelist != "unknown" ) {
    string flistroot = "PU/" + filelist + ".root";
    froot = TFile::Open(flistroot.c_str(),"read");
    if( froot && !froot->IsZombie() && froot->IsOpen() ) {
      /// one file containing all required info has already been setup
      hEvents = (TH1F*) froot->Get("hEvents");
      hPU     = (TH1F*) froot->Get("hPU");
    } else {
      /// setup a root file with all info for this filelist
      froot = TFile::Open(flistroot.c_str(),"recreate");
      
      ifstream file(filelist.c_str());
      cout << " Opening file: " << filelist << endl;
      while ( file.good() && !file.eof() ) {
	string fname;
	getline(file,fname,'\n');
	cout << " --> PUandNevt: Adding data from file: " << fname << endl;
	TH1F *hEvtLoc = 0;
	TH1F *hPULoc  = 0;
	if( fname.find(filename ) != string::npos ) continue;
	TFile *fLoc = TFile::Open(fname.c_str(),"read");
	if( fLoc && !fLoc->IsZombie() && fLoc->IsOpen()  ) {
	  gDirectory->cd("ggNtuplizer");  
	  hEvtLoc = (TH1F*)  gDirectory->Get("hskim");
	  if( hEvtLoc == 0 )
	    hEvtLoc = (TH1F*)  gDirectory->Get("hEvents");
	  
	  hPULoc = (TH1F*) gDirectory->Get("hPU");
	  if( hEvtLoc ) cout << " adding nevts = " << hEvtLoc->GetBinContent(1) << endl;
	  if( hEvtLoc ) hEvents->Add( hEvtLoc );
	  if( hPULoc  ) hPU->Add(hPULoc );
	  fLoc->Close();
	} 
      }
      
      froot->cd();
      if( hPU     ) hPU->Write("hPU");
      if( hEvents ) hEvents->Write("hEvents");
    }
  }
    
  TFile *fpudata = TFile::Open(_pu_data_file.c_str(), "read");
  TH1F *hPUdata  = 0; hPUdata = (TH1F*) fpudata->Get("pileup"); 
  
  /// normalise both MC and data distributions
  if( hPUdata ) hPUdata->Scale( 1./hPUdata->Integral() );
  if( hPU     ) hPU->Scale( 1./ hPU->Integral() );
    
  mydir->cd();
  TH1F *hResult = 0;
  if( hPU != 0 ) {
    hResult = (TH1F*) hPU->Clone("puWeights");
    hResult->Reset();
    for( int ib = 1; ib < hResult->GetNbinsX(); ib++ ) {
    /// MC bining no so appropriate becase beginning and end of bins are integer
    /// this might induce confusion via machine error limits (add 0.0001 to be sure which pu).
      float npu = hResult->GetXaxis()->GetBinLowEdge(ib) + 0.0001;
      int ib_data = hPUdata->GetXaxis()->FindBin( npu  );
      float pu_data = -1;
      if( ib_data >= hPUdata->GetNbinsX() ) pu_data = 0;
      else if( ib_data == 0 ) pu_data = 0; /// this one should not happen
      else  pu_data = hPUdata->GetBinContent(ib_data);
      
      float pu_mc = -1;
      if( hPU ) pu_mc = hPU->GetBinContent(ib);
      if( pu_mc > 0.000000001 ) 
	hResult->SetBinContent(ib, pu_data/pu_mc);
    }
  }

  if( hEvents ) _nEvtTot = hEvents->GetBinContent(1);

  if( hPU == 0 ) {
    hResult = 0;
    cout << " SetPUandNevts: no PU reweighting available... (no problem for data of course) " << endl;
  }
  fpudata->Close();
  if( froot && froot->IsOpen() ) froot->Close();
  delete hPU;
  hnPV = hResult;
}



//--------------------- object identification ----------------------------//
/// photon identification
float WeightManager::phoIdPresel( float r9, float sceta ) {
  int cat = -1;
  if( r9 >= 0.90 && fabs(sceta) <  1.45 ) cat = 0;
  if( r9 <  0.90 && fabs(sceta) <  1.45 ) cat = 1;
  if( r9 >= 0.90 && fabs(sceta) >= 1.45 ) cat = 2;
  if( r9 <  0.90 && fabs(sceta) >= 1.45 ) cat = 3;	      
  return photonIdCorrPresel[cat][iLHC_sTeV];
}
float WeightManager::phoIdMVA(    float r9, float sceta ) {

  /// FC: for now return 1, this is just a check
  return 1;

  int cat = -1;
  if( r9 >= 0.94 && fabs(sceta) <  1.45 ) cat = 0;
  if( r9 <  0.94 && fabs(sceta) <  1.45 ) cat = 1;
  if( r9 >= 0.94 && fabs(sceta) >= 1.45 ) cat = 2;
  if( r9 <  0.94 && fabs(sceta) >= 1.45 ) cat = 3;		      
  return photonIdCorrMVA[cat][iLHC_sTeV];

}
float WeightManager::phoIdCiC(    float r9, float sceta ){
  int cat = -1;
  if( r9 >= 0.94 && fabs(sceta) <  1.45 ) cat = 0;
  if( r9 <  0.94 && fabs(sceta) <  1.45 ) cat = 1;
  if( r9 >= 0.94 && fabs(sceta) >= 1.45 ) cat = 2;
  if( r9 <  0.94 && fabs(sceta) >= 1.45 ) cat = 3;		      
  return photonIdCorrCiC[cat][iLHC_sTeV];
}

float WeightManager::phoIdEleVeto(    float r9, float sceta ){
  int cat = -1;
  if( r9 >= 0.94 && fabs(sceta) <  1.45 ) cat = 0;
  if( r9 <  0.94 && fabs(sceta) <  1.45 ) cat = 1;
  if( r9 >= 0.94 && fabs(sceta) >= 1.45 ) cat = 2;
  if( r9 <  0.94 && fabs(sceta) >= 1.45 ) cat = 3;		      
  return photonIdCorrEleVeto[cat][iLHC_sTeV];
}

/// electron identification
float WeightManager::elecId( float sceta ) {
  int cat = -1;
  if     ( fabs(sceta) < 1.0000 ) cat = 0;
  else if( fabs(sceta) < 1.4442 ) cat = 1;
  else if( fabs(sceta) < 2.1000 ) cat = 2;
  else                            cat = 3;
 
  return electronIdCorr[cat][iLHC_sTeV];
}
