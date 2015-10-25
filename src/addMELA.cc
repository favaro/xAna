#include "JHUGenMELA/TEvtProb.hh"
#include "JHUGenMELA/TVar.hh"
#include "JHUGenMELA/TEvtProb.hh"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>

#include <iostream>
#include <stdio.h>
#include <string>
using namespace std;

int main() {

  //TString infileName  = "minitrees_data/data_8TeV_skimMVA_runABCD_withMELA.root";
  //TString infileName = "/afs/cern.ch/work/f/favaro/private/amina_ggH_0p.root";
  //  TString infileName  = "minitree_jhu_8TeV_0p_VBF_125_noduplicates.root";
  //TString infileName  = "fortest_runB.root";
  TString infileName = "/afs/cern.ch/work/f/favaro/public/minitreesMELAStretch3/CIC/data/data_8TeV_skimMVA_runABCD_withMELA.root";
  

  TFile filein( infileName, "update" );


  // matrix elements:
  float mela_dXsec_HJJ;
  float mela_dXsec_HJJVBF;
  float mela_dXsec_HJJ_PS;
  float mela_dXsec_HJJVBF_PS;
  // discriminants
  float mela_VBFvsgg;
  float mela_SMvsPS_VBF;
  float mela_SMvsPS_gg;  

  
  TTree *treein = (TTree*)filein.Get("HToGG");
  
  // branches in **************************************************

  TBranch *bHJJ      = treein->Branch( "dXsec_HJJ",       &mela_dXsec_HJJ,       "dXsec_HJJ/F");
  TBranch *bHJJVBF   = treein->Branch( "dXsec_HJJVBF",    &mela_dXsec_HJJVBF,    "dXsec_HJJVBF/F");
  TBranch *bHJJPS    = treein->Branch( "dXsec_HJJ_PS",    &mela_dXsec_HJJ_PS,     "dXsec_HJJ_PS/F");
  TBranch *bHJJVBFPS = treein->Branch( "dXsec_HJJVBF_PS", &mela_dXsec_HJJVBF_PS, "dXsec_HJJVBF_PS/F");

  TBranch *bVBFvsGG  = treein->Branch( "mela_VBFvsgg",    &mela_VBFvsgg,    "mela_VBFvsgg/F");
  TBranch *bSMvsPSVBF = treein->Branch( "mela_SMvsPS_VBF", &mela_SMvsPS_VBF, "mela_SMvsPS_VBF/F");
  TBranch *bSMvsPSGG = treein->Branch( "mela_SMvsPS_gg",  &mela_SMvsPS_gg,  "mela_SMvsPS_gg/F");

  // branches out **************************************************

  Float_t         nJet;
  Float_t         mass;
  Int_t           tagCat;
  Float_t         ptj[2];   //[nent]
  Float_t         phij[2];   //[nent]
  Float_t         etaj[2];   //[nent]
  Float_t         enj[2];   //[nent]
  Float_t         ptg[2];   //[nent]
  Float_t         phi[2];   //[nent]
  Float_t         eta[2];   //[nent]
    
  // List of branches

  TBranch        *b_mass;
  TBranch        *b_tagCat;
  TBranch        *b_nJet;  
  TBranch        *b_ptj;   
  TBranch        *b_phij;   
  TBranch        *b_etaj;   
  TBranch        *b_enj;   
  TBranch        *b_ptg;   
  TBranch        *b_phi;   
  TBranch        *b_eta;  
  
  treein->SetBranchAddress("nJet", &nJet, &b_nJet);
  treein->SetBranchAddress("ptj", ptj, &b_ptj);
  treein->SetBranchAddress("phij", phij, &b_phij);
  treein->SetBranchAddress("etaj", etaj, &b_etaj);
  treein->SetBranchAddress("enj", enj, &b_enj);
  treein->SetBranchAddress("ptg", ptg, &b_ptg);
  treein->SetBranchAddress("phi", phi, &b_phi);
  treein->SetBranchAddress("eta", eta, &b_eta);
  treein->SetBranchAddress("mass", &mass, &b_mass);
  treein->SetBranchAddress("tagCat", &tagCat, &b_tagCat);
  
  TVar::VerbosityLevel verbosity = TVar::ERROR;    //INFO
  TEvtProb Xcal2;

  unsigned int nentries = treein->GetEntries();

  for( unsigned int ie = 0; ie < nentries; ie++ ) {
    
    treein->GetEntry(ie);
    
    mela_dXsec_HJJ = -9;
    mela_dXsec_HJJVBF = -9; 
    mela_dXsec_HJJ_PS = -9;
    mela_dXsec_HJJVBF_PS= -9;
    mela_VBFvsgg = -9;
    mela_SMvsPS_VBF = -9;
    mela_SMvsPS_gg = -9; 
        
    //if( ptj[0]<0 || ptj[1]<0 || mass < 120 || mass > 140 ) continue;
    if( tagCat ==8 ) {

      //cout << " pt " <<  ptj[0] << " " <<  ptj[1] << endl;
      
      unsigned int jL, jR;
      if( etaj[0] > etaj[1] ) { jL = 0; jR = 1; }
      else{ jL = 1; jR = 0; }
      
      TLorentzVector pho1, pho2;
      pho1.SetPtEtaPhiM( ptg[0], eta[0], phi[0], 0);  pho2.SetPtEtaPhiM( ptg[1], eta[1], phi[1], 0); 
      TLorentzVector kinematics[3];
      kinematics[0].SetPtEtaPhiE( ptj[jL], etaj[jL], phij[jL], enj[jL] );
      kinematics[1].SetPtEtaPhiE( ptj[jR], etaj[jR], phij[jR], enj[jR] );
      kinematics[2] = pho1 + pho2;

      mela_dXsec_HJJ       = Xcal2.XsecCalcXJJ(TVar::HJJNONVBF, kinematics, verbosity);
      mela_dXsec_HJJVBF    = Xcal2.XsecCalcXJJ(TVar::HJJVBF, kinematics, verbosity);
      
      mela_dXsec_HJJ_PS    = Xcal2.XsecCalcXJJ(TVar::HJJNONVBF_PS, kinematics, verbosity);
      mela_dXsec_HJJVBF_PS = Xcal2.XsecCalcXJJ(TVar::HJJVBF_PS, kinematics, verbosity);
      
      mela_VBFvsgg    = ( mela_dXsec_HJJVBF/(mela_dXsec_HJJVBF + 3*mela_dXsec_HJJ) );   
      mela_SMvsPS_VBF = ( mela_dXsec_HJJVBF/(mela_dXsec_HJJVBF+0.1*mela_dXsec_HJJVBF_PS) );
      mela_SMvsPS_gg  = ( mela_dXsec_HJJ/(mela_dXsec_HJJ+mela_dXsec_HJJ_PS) );

      
      /* cout << " mela_SMvsPS_gg " << mela_SMvsPS_gg  
	  << " 1: pt " << kinematics[0].Pt() 
	   << " eta " << kinematics[0].Eta() 
	   << " phi " << kinematics[0].Phi()
	   << " 2: pt " << kinematics[1].Pt() 
	   << " eta " << kinematics[1].Eta() 
	   << " phi " << kinematics[1].Phi()
	   << " gamma: pt " << kinematics[2].Pt() 
	   << " eta " << kinematics[2].Eta() 
	   << " phi " << kinematics[2].Phi()
	   << endl;
      */
    }
    
    bHJJ->Fill();
    bHJJVBF->Fill();
    bHJJPS->Fill();
    bHJJVBFPS->Fill();
    
    bVBFvsGG->Fill();
    bSMvsPSVBF->Fill();
    bSMvsPSGG->Fill();
    
  }
  
  treein->Write("",TObject::kOverwrite);

}
