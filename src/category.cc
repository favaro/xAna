#include "interface/category.hh"
#include "interface/configReader.hh"

#include <iostream>
#include <fstream>

using namespace std;

category::category( TString Name, int subCat, bool isTag ):
  name(Name),
  subcat(subCat),
  istag(isTag),
  ncuts(-1) { 
 
}

category::~category() {}

//*********************************************************

void category::Init( const ConfigReader *config ) {
  setCommonCuts( config );
  setSpecificCuts( config );
}
  



void category::setCommonCuts( const ConfigReader *config ) {

  cut_pt_lead    = config->pt1Cut();
  cut_pt_trail   = config->pt2Cut();
  cut_ptOm_lead  = config->loosePt1MCut();
  cut_ptOm_trail = config->loosePt2MCut(); 
 
  if( name == "electag" || name == "muontag" )  cut_ptOm_lead  = config->VHPt1MCut();
  if( name == "vh" )  cut_ptOm_lead  = config->VHPt1MCut();// JM
  if( name == "vhhad" )  cut_ptOm_lead  = config->VHhadPt1MCut();// JM
  if( name == "tth")  cut_ptOm_lead  = config->TTHPt1MCut();  
  if( name == "vbf"){
    cut_ptOm_lead   = config->diJetPtg1M();  // JM
    cut_ptOm_trail  = config->diJetPtg2M();  // JM
  }
  if( name == "electag" || name == "muontag" )  cut_ptOm_lead  = config->VHPt1MCut();
  cut_mgg        = config->mggCut(); 

}

void category::setSpecificCuts( const ConfigReader *config ){

  ifstream listOfCuts;
  TString filename = config->categoryDefDir() + "/" + name + Form("_%d",subcat) + ".txt";
 
  listOfCuts.open(filename);
  listOfCuts >> ncuts;
  cout << "- Category: " << name << endl;

  if(cut_pt_lead!= config->pt1Cut()) cout << "      * pt1 > " <<cut_pt_lead<< endl;
  if(cut_pt_trail!= config->pt2Cut()) cout << "      * pt2 > " <<cut_pt_trail<< endl;
  if(cut_ptOm_lead != config->loosePt1MCut())  cout << "      * pt1/M > " <<cut_ptOm_lead<< endl;
  if(cut_ptOm_trail != config->loosePt2MCut()) cout << "      * pt2/M > " <<cut_ptOm_trail<< endl;

  for( size_t ic = 1; ic < ncuts+1; ic++ ) {
    cut onecut;
    TString side;
    TString sign;
    listOfCuts >> onecut.first >> side >> onecut.second;
    specificCuts.push_back(onecut);
    specificCuts_side.push_back(side);
    if(side=="max")sign=" < ";
    else if(side=="min")sign=" > ";
    
    cout << "      * " << onecut.first << sign << onecut.second << endl; 
  }
}

//*********************************************************

bool category::isGoodPhotonPair( HggEvtCandidate hggevent ){


  double hMass = ( hggevent.gammaK_lead + hggevent.gammaK_trail ).M();

  if( hggevent.gammaK_lead.Pt() / hMass < cut_ptOm_lead || hggevent.gammaK_trail.Pt() / hMass < cut_ptOm_trail ) return false;
  if( hggevent.gammaK_lead.Pt() < cut_pt_lead || hggevent.gammaK_trail.Pt() < cut_pt_trail ) return false;
  if( hMass < cut_mgg ) return false;

  // add lepton cuts here 
 
  return true;  
}
//*********************************************************

