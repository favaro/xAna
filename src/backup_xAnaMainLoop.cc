#include "interface/xAna.hh"
#include "interface/MinitreeOut.hh"

#include <TH2.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <Math/VectorUtil.h>
#include <TGraphAsymmErrors.h>

#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace std;

#include "TSystem.h" // For memory check
ProcInfo_t pi; // For memory check


string pileupDir = "etc/inputs/pileup/pu2012/";

int xAna::HggTreeWriteLoop(const char* filename, int ijob  ) {

  //bool invDRtoTrk = false;


  TString fileN = "debug.txt";
  debugfile.open(fileN);


  if( _config == 0 ) {
    cout << " config file was not set properly... bail out" << endl;
    return -1;
  }


  cout << " ============= Define Categories ================" << endl;
  vector<category> tagCategories;
  vector<category> untagCategories;
  createCategories( tagCategories, untagCategories );
  cout << " number of tag = " << tagCategories.size() << " number of untag = " << untagCategories.size() << endl;
  cout << " ============= Categories defined ================" << endl;
  
  


  
  //  enScaleSkimEOS.setup( "ecalCalibFiles/EnergyScale2012_Lisbon_9fb.txt" );
  enScale.setup( _config->energyScaleFile() );

  Float_t HiggsMCMass =  _weight_manager->getCrossSection()->getHiggsMass();
  Float_t HiggsMCPt   = -1;
  bool isHiggsSignal = false;
  if( HiggsMCMass > 0 ) isHiggsSignal = true;

  mode_ = _weight_manager->getCrossSection()->mode();
  isData = false;
  if(mode_==-1) isData = true;        

  Setup_MVA();
  DiscriVBF_UseDiPhoPt = true;
  DiscriVBF_UsePhoPt   = true;
  DiscriVBF_useMvaSel      = _config->doVBFmvaCat();
  DiscriVBF_useCombMvaSel  = _config->doVBFCombmvaCat();

  bool DoOverSmearing = true;  
  MassResolution massResoCalc;
  if(DoOverSmearing) massResoCalc.addSmearingTerm();



  /// depending on the selection veto or not on electrons (can do muele, elemu,eleele)
  vetoElec[0] = true; vetoElec[1] = true;   
  if( _config->invertElectronVeto() ) { vetoElec[0] = false; vetoElec[1] = false; }
  if( _config->isEleGamma()         ) { vetoElec[0] = false; vetoElec[1] = true ; }
  if( _config->isGammaEle()         ) { vetoElec[0] = true ; vetoElec[1] = false; }
  cout << " --------- veto electron config -----------" << endl;
  cout << " Leading  Pho VetoElec: " << vetoElec[0] << endl;
  cout << " Trailing Pho VetoElec: " << vetoElec[1] << endl;
  

  DoDebugEvent = false;  
  bool DoPreselection = true;
  
  
  TRandom3 *rnd = new TRandom3();
  rnd->SetSeed(0);
  
  _xcheckTextFile.open(TString(filename)+".xcheck.txt");

  if (fChain == 0) return -1;
  if( _config->setup() == "Legacy2013_8TeV" && !isData ){
    hPURunAB= new TH1F("hPURunAB", "number of pileup", 100, 0, 100);
    hPURunC=new TH1F("hPURunC", "number of pileup", 100, 0, 100);
    hPURunD=new TH1F("hPURunD", "number of pileup", 100, 0, 100);
    fChain->Draw("nPU>>hPURunAB", "run <= 197495","goff");
    fChain->Draw("nPU>>hPURunC" , "run > 197495 && run <= 203767","goff");
    fChain->Draw("nPU>>hPURunD" , "run > 203767","goff");
  }

  Long64_t nentries = fChain->GetEntriesFast();
  //  nentries = 10000;
  cout << "nentries: " << nentries << endl;   
  Long64_t entry_start = ijob    *_config->nEvtsPerJob();
  Long64_t entry_stop  = (ijob+1)*_config->nEvtsPerJob();
  if( _config->nEvtsPerJob() < 0 ) {
    entry_stop = nentries;
  }
  if( entry_stop  > nentries ) entry_stop = nentries;
  cout << "   *** doing entries from: " << entry_start << " -> " << entry_stop << endl;
  if( entry_start > entry_stop ) return -1;

  _minitree = new MiniTree( filename );
  TTree * tSkim = 0;
  if( _config->doSkimming() ) tSkim = (TTree*) fChain->CloneTree(0);

  InitHists();
  _minitree->mc_wXsec = _weight_manager->xSecW();
  _minitree->mc_wNgen = 100000./_weight_manager->getNevts();
  if( isData ) {
    _minitree->mc_wXsec = 1;
    _minitree->mc_wNgen = 1;
  }
    
  Int_t isprompt0 = -1;
  Int_t isprompt1 = -1;
  
  set<Long64_t> syncEvt;

  cout <<" ================ mode  "  << mode_   <<" ===============================  "<<endl;
  /// setupType has to be passed via config file
  //  photonOverSmearing overSmearICHEP("Test52_ichep");
  //photonOverSmearing overSmearHCP( "oversmear_hcp2012" );

  //photonOverSmearing overSmear(    _config->setup()    );
  //int overSmearSyst = _config->getEnergyOverSmearingSyst();
  overSmear_.Initialize( _config->setup() );
  overSmearSyst_ = _config->getEnergyOverSmearingSyst();

  Long64_t nbytes = 0, nb = 0;

  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////// Start the loop ////////////////////////////////***************************************************
  ////////////////////////////////////////////////////////////////////////////

  vector<int>    nEvts;
  vector<string> nCutName;
  nCutName.push_back("TOTAL                      : "); nEvts.push_back(0);
  nCutName.push_back("2 gammas                   : "); nEvts.push_back(0);
  nCutName.push_back("triggers                   : "); nEvts.push_back(0);
  nCutName.push_back("nan weight                 : "); nEvts.push_back(0);
  nCutName.push_back("presel kin cuts            : "); nEvts.push_back(0);
  nCutName.push_back("pass all                   : "); nEvts.push_back(0);

  // JM add category counters
  int lmax=7;
  for(unsigned ict=0;ict< tagCategories.size();ict++){
    TString catstring=" --> tag category ";
    catstring+=tagCategories[ict].name; catstring+=" ";
    catstring+=tagCategories[ict].subcat;
    int lenght=tagCategories[ict].name.Length();     
    if (lenght>lmax) lenght=lmax;
    for(int i=0;i<lmax-lenght;i++)  catstring+=" ";
    catstring+=": ";
    nCutName.push_back(catstring.Data()); nEvts.push_back(0);
  } 
  lmax-=2;
  for(unsigned icu=0;icu< untagCategories.size();icu++){
    TString catstring=" --> untag category ";
    //catstring+=untagCategories[icu].name; catstring+=" ";
    catstring+=untagCategories[icu].subcat;
    catstring+=" ";
    int lenght=0;
    if (lenght>lmax) lenght=lmax;
    for(int i=0;i<lmax-lenght;i++)  catstring+=" ";
    catstring+=": ";
    nCutName.push_back(catstring.Data()); nEvts.push_back(0);   
  }
  


  for (Long64_t jentry=entry_start; jentry< entry_stop ; ++jentry) {

    cout << "debug 1 " << endl;
    // **************************
    TString filen = fChain->GetCurrentFile()->GetName();
    if (nb < 0){cerr << "corrupted events jentry in file: " << filen << endl;  continue; }

    Long64_t ientry = LoadTree(jentry);
    if (ientry < -9999) break;
    if( jentry % 100 == 0) cout<<"processing event "<<jentry<<endl;
    nb = GetEntry(jentry);   nbytes += nb;

    /// reset minitree variables
    _minitree->initEvent();
    /// study mc truth block
    if( !isData ) {
      fillMCtruthInfo_HiggsSignal();

      //cout << " done" << endl;
      _minitree->fillMCtrueOnly();
    }

    cout << "debug 1 A" << endl;
    /// reco analysis
    unsigned icutlevel = 0;
    nEvts[icutlevel++]++;
    if( nPho < 2 ) continue; 
    nEvts[icutlevel++]++;
    
    /// set synchronisation flag
    if( syncEvt.find(event) != syncEvt.end() ) DoDebugEvent = true;
    else                                       DoDebugEvent = false;
   
          /// PU & BSz reweightings
    if( !isData ) {
      int itpu = 1;  /// 0 without OOT PU - 1 with OOT PU

      if(_config->setup() == "Legacy2013_8TeV"){
	if (run <= 197495){
	  _weight_manager->setRDPUMC(hPURunAB);
	  _weight_manager->setPUTarget( pileupDir + "AB.json.69400.observed.pileup.root");
	}
	else if (run > 197495 && run <= 203767){
	  _weight_manager->setRDPUMC(hPURunC);
	  _weight_manager->setPUTarget( pileupDir + "C.json.69400.observed.pileup.root");
     
	}
	else if (run > 203767){
	  _weight_manager->setRDPUMC(hPURunD);
	  _weight_manager->setPUTarget( pileupDir + "D.json.69400.observed.pileup.root");
	}
      }

      cout << "debug 1 D " <<  endl;

      if( _config->setup() == "Legacy2013_7TeV") _weight_manager->setPUTarget( pileupDir +  "nov08_rereco.json.68000.pileup.root"); 
      cout << " pile up " << nPU << endl;
      //_minitree->mc_wPU  = _weight_manager->puW( (*nPU)[itpu] );
      // PUwei  = _weight_manager->puWTrue( puTrue[itpu] ); 
    }
    hTotEvents->Fill(1,_minitree->mc_wPU);

    cout << "debug 2" << endl;    

    bool sigWH = false ;
    bool sigZH = false ;    
    int  mc_whzh_type = 0;
    _minitree->mc_wHQT = 1;
    if( isHiggsSignal ) {
      
	 if ( _weight_manager->getCrossSection()->getMCType() == "vh" ) 
	 for( Int_t i=0; i < nMC && i <= 1; ++i ) { 
	 if( abs((*mcPID)[i]) == 24 ) sigWH=true;
	 if( abs((*mcPID)[i]) == 23 ) sigZH=true;
	 }
	 if( sigWH ) mc_whzh_type = 1;
	 if( sigZH ) mc_whzh_type = 2;
	 
      for( Int_t i=0; i<nMC; ++i)
	if ( abs((*mcPID)[i]) == 25 ) HiggsMCPt   = (*mcPt)[i];
      
      if( _weight_manager->getCrossSection()->getMCType() == "ggh" && 
	  _config->setup().find( "Legacy2013_7TeV" ) != string::npos )
	_minitree->mc_wHQT = getHqTWeight(HiggsMCMass, HiggsMCPt);      
    }
    _minitree->mc_wXsec  = _weight_manager->xSecW(mc_whzh_type);      

    if ((mode_ == 2 || mode_ ==1 ||  mode_ == 18 || mode_ == 19) && processID==18) continue;      
    
    // Remove double counting in gamma+jets and QCDjets
    
    Int_t mcIFSR_pho = 0;
    Int_t mcPartonic_pho = 0;
    if (mode_ == 1 || mode_ == 2 || mode_ == 3  ||  mode_ == 18 ||  mode_ == 19 ) {
      for (Int_t i=0; i<nMC; ++i) {    
	if ((*mcPID)[i] == 22 && (fabs((*mcMomPID)[i]) < 6 || (*mcMomPID)[i] == 21)) mcIFSR_pho++;
	if ((*mcPID)[i] == 22 && (*mcMomPID)[i] == 22) mcPartonic_pho++;
      }
    }
    
    // if pythia is used for diphoton.. no IFSR removing from QCD and Gjets!!!!!!      
    if ((mode_==1 || mode_ == 2  ||   mode_ == 18 ) && mcIFSR_pho >= 1 && mcPartonic_pho >=1) continue;
    if ((mode_ == 3 ||  mode_ == 19 )&& mcIFSR_pho == 2) continue;   

    
    //    bool prompt2= false;
    bool prompt1= false;
    bool prompt0= false;
    vertProb = -1;
    
    if (mode_ == 2  || mode_ == 1    ||   mode_ == 18     ){
      if ( mcPartonic_pho >= 1 &&  mcIFSR_pho == 0 ) prompt1 = true;
      //else if      ( mcPartonic_pho >= 1 &&  mcIFSR_pho >= 1 ) prompt2 = true;
      //else if ( mcPartonic_pho == 0 &&  mcIFSR_pho == 0 ) prompt0 = true;
    } else if(mode_ == 3 ||  mode_ == 19   ){
      if ( mcIFSR_pho == 1 ) prompt1 = true;
      //else if  ( mcIFSR_pho >= 2 ) prompt2 = true;
      //else if ( mcIFSR_pho == 0 ) prompt0 = true;
      if( prompt1 ) _minitree->mc_wXsec = 1.3*_weight_manager->xSecW();
    }
    
    
      if(mode_==1 || mode_==2 || mode_==3 || mode_==18 || mode_==19){
      if(prompt0)isprompt0=1;
      else isprompt0=0;
      if(prompt1)isprompt1=1;
      else isprompt1=0;
    }
      
    
    if( mode_ == 20 && isZgamma() ) continue;

    /// wei weight is just temporary and may not contain all info.
    float wei = _minitree->mc_wXsec * _minitree->mc_wPU  
      * _minitree->mc_wNgen * _minitree->mc_wHQT;

    if( isData && !PassTriggerSelection2012() ) continue;  nEvts[icutlevel++]++;
    if( std::isinf( wei ) || std::isnan( wei ) )continue;  nEvts[icutlevel++]++;
    
    //// ********** Apply energy corrections to data and MC ************* ////
    for( int iphoton = 0; iphoton < nPho; iphoton++ )  ApplyEnergyCorrections( DoOverSmearing, iphoton );
    
    cout << " debug 3 "<< endl;

    //// ************************************************************************************************** ////
    ////                                 EVENT SELECTION AND CATEGORIES
    //// ************************************************************************************************** ////

    vector<HggEvtCandidate> selectedCand =  buildHggCandidates( DoPreselection );

    cout << "debug 3 A AA" << endl;
    if( selectedCand.size() == 0 ) continue;
    nEvts[icutlevel++]++;

    cout << "debug 3 A" << endl;

    HggEvtCandidate selectedEvt = sortHggCandidates( selectedCand, tagCategories, untagCategories, massResoCalc );
    cout << "debug 3 C" << endl;
    if( selectedEvt.gammaIdx_lead == -1 && selectedEvt.gammaIdx_trail == -1 ) continue;
    nEvts[icutlevel++]++;

    cout << "debug 3 D" << endl;
        
    // reapply all functions just to fill the minitree:
    evaluateDiPhotonMVA( selectedEvt.gammaIdx_lead, selectedEvt.gammaIdx_trail, selectedEvt.gammaK_lead,  selectedEvt.gammaK_trail, massResoCalc  );
    
    cout << "debug 3 E" << endl;

    vector<int> selVtxAgain;      
    selVtxAgain = getSelectedVertex( selectedEvt.gammaIdx_lead, selectedEvt.gammaIdx_trail, true );
    
    cout << "debug 4" << endl;
    // fill minitree: **************************************************
    
    _minitree->mtree_runNum  = run;
    _minitree->mtree_evtNum  = event;
    _minitree->mtree_lumiSec = lumis;
    
    _minitree->mtree_rho   = rho2012;
    _minitree->mtree_rho25 = rho25;    
    _minitree->mtree_zVtx     = (*vtxbs_z)[selectedEvt.vertexIdx];
    _minitree->mtree_nVtx     = nVtxBS;
    _minitree->mtree_nVtxNoBS = nVtx;
    
    _minitree->mtree_ivtx1    = selVtxAgain[0];
    _minitree->mtree_ivtx2    = selVtxAgain[1];
    _minitree->mtree_ivtx3    = selVtxAgain[2];
    _minitree->mtree_vtxProb  = vertProb;
    _minitree->mtree_vtxMva   = vertMVA;

   
    TLorentzVector Hcand;
    Hcand = selectedEvt.gammaK_lead + selectedEvt.gammaK_trail;
    _minitree->mtree_mass     = Hcand.M();
    _minitree->mtree_pt       = Hcand.Pt();
    _minitree->mtree_piT      = Hcand.Pt()/Hcand.M();
    _minitree->mtree_y        = Hcand.Rapidity();
  
    /// spin variables

    double eta2[2];
    eta2[0] = +log( Hcand.Pt() ) - log(Hcand.E()-Hcand.Pz());
    eta2[1] = -log( Hcand.Pt() ) + log(Hcand.E()+Hcand.Pz());
    TLorentzVector gtmp1, gtmp2;
    gtmp1.SetPtEtaPhiM( Hcand.Pt(), eta2[0], Hcand.Phi(),0); gtmp1.Boost( -Hcand.BoostVector() );
    gtmp2.SetPtEtaPhiM( Hcand.Pt(), eta2[1], Hcand.Phi(),0); gtmp2.Boost( -Hcand.BoostVector() );
    _minitree->mtree_cThetaLead_heli  = cos( gtmp1.Angle(Hcand.BoostVector()) );
    _minitree->mtree_cThetaTrail_heli = cos( gtmp2.Angle(Hcand.BoostVector()) );
    _minitree->mtree_cThetaStar_CS    =
      2*(selectedEvt.gammaK_lead.E()*selectedEvt.gammaK_trail.Pz() - selectedEvt.gammaK_trail.E()*selectedEvt.gammaK_lead.Pz())   / 
      (Hcand.M()*sqrt(Hcand.M2()+Hcand.Pt()*Hcand.Pt()));


    fillJetVariablesToMinitree( selectedEvt );
   
    int finalCat = convertToOfficialCat( selectedEvt.category, selectedEvt.subcategory );
    _minitree->mtree_tagCat   = finalCat;
    _minitree->mtree_untagCat = selectedEvt.untagsubcat;
    
    for(unsigned ict=0;ict< tagCategories.size();ict++){          
      if(selectedEvt.category==tagCategories[ict].name && selectedEvt.subcategory==tagCategories[ict].subcat )  nEvts[icutlevel]++;
       icutlevel++;
    }
    for(unsigned icu=0;icu< untagCategories.size();icu++){
      if(finalCat==-1 && selectedEvt.untagsubcat==untagCategories[icu].subcat )  nEvts[icutlevel]++;
      icutlevel++;
    }

    cout << "debug 5" << endl;
    //----- MC weights
    if( !isData ) {

      wei = computeMCweights(selectedEvt,isprompt1);
      cout << "debug 5bis" << endl;
      recoToTruthMatching( selectedEvt );
    }
    
    cout << "debug 6" << endl;

    _minitree->fill();    
    if( _config->doSkimming() && tSkim ) {
      fChain->GetEntry(jentry);
      tSkim->Fill();
    }
    //---------------      
    
  }// end for entry
  
  
  
  
  if( _config->doSkimming() ) {
    _minitree->mtree_file->cd();
    TH1F *hEvents = new TH1F("hEvents","hEvents",2,0,2);
    hEvents->SetBinContent(1,_weight_manager->getNevts());
    hEvents->SetBinContent(2,nEvts[nEvts.size()-1]);
    hEvents->Write();
    tSkim->Write();
    
    _minitree->mtree_file->Close();
  } else  _minitree->end();
  

  cout<<"============= Events summary ==========="<<endl;

  for( unsigned icut = 0 ; icut < nEvts.size(); icut++ )
    cout << nCutName[icut] << nEvts[icut] << endl;
  
  delete rnd;  
  
  return nEvts[nEvts.size()-1];
  
}



