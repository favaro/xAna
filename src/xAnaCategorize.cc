#include "interface/xAna.hh"
#include "interface/hggCandidate.hh"
#include "interface/MinitreeOut.hh"

using namespace std;


//******************************************************************************************************

int xAna::convertToOfficialCat( const TString catname, const int subcatidx  ) {

  int idxtmp = -1;

  if( _config->analysisType() == "MVA" ) {
    if( catname == "muontag" )      idxtmp = 8;   // THIS IS NOT A CORRECT CATEGORY 
    else if( catname == "electag" ) idxtmp = 9;   // THIS IS NOT A CORRECT CATEGORY 
    else if( catname == "vbf" )     idxtmp = subcatidx + 5;
    else if( catname == "vh"  ){
      idxtmp = subcatidx + 8;
      if( subcatidx==3 || subcatidx==4 ) idxtmp = 8; // JM dilepton categories 
    }
    else if( catname == "tth" )     idxtmp = subcatidx + 11;  // JM
    else if( catname == "vhhad"  )     idxtmp = subcatidx + 13;
  } 
  else if( _config->analysisType() == "baselineCiC4PF" ||  _config->analysisType() == "" ) {
    if( catname == "muontag" )      idxtmp = 10;  // THIS IS NOT A CORRECT CATEGORY 
    else if( catname == "electag" ) idxtmp = 11;   // THIS IS NOT A CORRECT CATEGORY
    else if( catname == "vbf" )     idxtmp = subcatidx + 8;
    else if( catname == "vh"  ){
      idxtmp = subcatidx + 10;
      if( subcatidx==3 || subcatidx==4 ) idxtmp = 10; // JM dilepton categories 
    }
    else if( catname == "tth" )     idxtmp = subcatidx + 13;  // JM
    else if( catname == "vhhad"  )  idxtmp = subcatidx + 15;
  }

  return idxtmp;
}

//******************************************************************************************************


void xAna::createCategories( std::vector<category>& tag, std::vector<category>& untag ) {
  
  ifstream listOfCategories;
  TString listname; listname =  _config->categoryDefDir() + "/listOfCategories.txt";
  if( _config->doSkimming() )
    listname = "etc/categoryDef/" + _config->analysisType() + "/Skimming/listOfCategories.txt";
  cout << "file path " << listname << endl;
  listOfCategories.open(listname);

  size_t ncat;
  listOfCategories >> ncat;

  cout << "number of categories " << ncat << endl;
  for( size_t il = 1; il < ncat+1; il++ ){
   
    TString Name; int subCat; int isTag;
    listOfCategories >> Name >> subCat >> isTag;
    category oneCat( Name, subCat, isTag );
    
    if( oneCat.istag != 0 ) tag.push_back(oneCat);
    else if(  oneCat.istag == 0 ){
      untag.push_back(oneCat);
    }    
  }

  for( size_t ic = 0; ic < untag.size(); ic++ ){
    untag[ic].Init(_config);
  }
  for( size_t ic = 0; ic < tag.size(); ic++ ){
    tag[ic].Init(_config);
  }
   
}

//******************************************************************************************************


bool xAna::isEvtInCategory( HggEvtCandidate hggevent, category cat) {
  for( size_t icut = 0; icut < cat.ncuts; icut++ ){

    float thisvar = stringToFloat( cat.specificCuts[icut].first, hggevent ); 

    if( cat.specificCuts_side[icut] == "max" ){
      if( thisvar > cat.specificCuts[icut].second ) return false;     
    }
    else{
      if( thisvar < cat.specificCuts[icut].second ) return false;
    }  
  } 
  

  return true;
}


// *******************************************************************************************
// *******************************************************************************************

std::vector<HggEvtCandidate> xAna::buildHggCandidates(  bool DoPreselection ) {

  std::vector<HggEvtCandidate> hggcandidates;

  int iElecVtx(-1), iMuonVtx(-1);
  int iElecVtxDilep(-1), iMuonVtxDilep(-1);

  vector<int> elecIndex = selectElectronsHCP2012( iElecVtx );  
  vector<int> muonIndex = selectMuonsHCP2012(    iMuonVtx ); 
  
  // NB: cut at 20 GeV already implemented inside functions

  vector<int> elecIndexVHDilep = selectElectronsHCP2012( iElecVtxDilep , 10 );  
  vector<int> muonIndexVHDilep = selectMuonsHCP2012(    iMuonVtxDilep , 10 );

  for( int ii = 0     ; ii < nPho ; ++ii ) {
    for( int jj = (ii+1); jj < nPho ; ++jj ) {

      if (DoPreselection && !preselectPhoton(ii,(*phoEt)[ii])) continue;
      if (DoPreselection && !preselectPhoton(jj,(*phoEt)[jj])) continue;

      int i = ii; 
      int j = jj;
      if((*phoEt)[j] > (*phoEt)[i]){ i = jj; j = ii; }	 
      TLorentzVector gi,gj;  
      
      // Select vertex  for each photon pair ---------------

      vector<int> selVtxIncl; 
      int selVtx  = 0;
      
      selVtxIncl = getSelectedVertex( i, j, true );

      if( selVtxIncl.size() < 2 ) selVtxIncl.push_back(-1);
      if( selVtxIncl.size() < 3 ) selVtxIncl.push_back(-1); 
      selVtx = selVtxIncl[0];
      if( selVtx < 0 ) continue;

      // ---------------------------------------------------
      
      if( !selectTwoPhotons(i, j, selVtx, gi, gj ) ) continue;    
      
      // Get MET ---------------------------------------------

      //TLorentzVector corMet = corMET(i, j);
     
      HggEvtCandidate hggtmp; 

      
      hggtmp.vertexIdx = selVtx;
      hggtmp.gammaIdx_lead  = i;
      hggtmp.gammaIdx_trail = j;
      hggtmp.gammaK_lead    = gi;
      hggtmp.gammaK_trail   = gj;
      // hggtmp.corMet         = corMet; //JM
      hggtmp.corMet         = -10.; //JM
      hggtmp.gammaR9_lead   = (*phoR9)[i];
      hggtmp.gammaR9_trail  = (*phoR9)[j];
      hggtmp.nmu   = muonIndex.size();
      for( size_t im = 0; im < muonIndex.size(); im++ ) hggtmp.muIdx.push_back( muonIndex[im] );
      hggtmp.category    = "none";
      hggtmp.subcategory = -1;
      hggtmp.untagsubcat = -1;


      // JM: Do as globe --------------------------------------
      // Remove electrons matching photons from electron list
      // This might be stupid but globe does it this way    
      // this selectElectronsAwayFromTwoPhotons should be removed and 
      // selElecIndex should be replaced elecIndex by when we are done with Synch studies

      vector<int> selElecIndex= selectElectronsAwayFromTwoPhotons(i, j, elecIndex );
      vector<int> selElecIndexVHDilep= selectElectronsAwayFromTwoPhotons(i, j, elecIndexVHDilep );

      hggtmp.nelec = selElecIndex.size(); //JM
      hggtmp.nmuVHDilep =  muonIndexVHDilep.size(); //JM
      hggtmp.nelecVHDilep = selElecIndexVHDilep.size();  //JM
      
      for( size_t ie = 0; ie < selElecIndex.size(); ie++ ) hggtmp.eleIdx.push_back( selElecIndex[ie] ); //JM
      for( size_t ie = 0; ie < selElecIndexVHDilep.size(); ie++ ) hggtmp.eleIdxVHDilep.push_back( selElecIndexVHDilep[ie]);//JM
      for( size_t im = 0; im < muonIndexVHDilep.size(); im++ ) hggtmp.muIdxVHDilep.push_back( muonIndexVHDilep[im] ); //JM

      
      // Select Hgg candidates away from leptons -------------------------
      
      // if( ! selectTwoPhotonsAwayFromLeptons(i, j, muonIndex , elecIndex ) ) continue; //WHAT SHOULD BETTER BE DONE? 
      if( ! selectTwoPhotonsAwayFromLeptons(i, j, muonIndex , selElecIndex ) ) continue; 
      
      if( ! selectTwoPhotonsAwayFromLeptons(i, j, muonIndexVHDilep , selElecIndexVHDilep ) ){ 
	//JM: fixme this could be done better 
	hggtmp.nmuVHDilep =0; 
	hggtmp.nelecVHDilep = 0;
      }

      // Select jets  --------------------------------------

      int nVtxJetID = -1;
      if( _config->doPUJetID() ) nVtxJetID = nVtxBS;
      vector<int> allGoodJetsIndex;
      vector<int> goodJetsIndex = selectJetsJEC(  i, j, nVtxJetID, selVtx , allGoodJetsIndex );
    
      hggtmp.njets = goodJetsIndex.size();
      hggtmp.jetIdx.clear();
      for( size_t iji = 0; iji < goodJetsIndex.size(); iji++ ) hggtmp.jetIdx.push_back(goodJetsIndex[iji]);
        
      // Select ttH jets and btag -----------------------  

      int nBjettmp      = 0;
      
      vector<int> goodTTHJetsIndex = selectJetsTTH( allGoodJetsIndex , selElecIndex, muonIndex, nBjettmp );
      hggtmp.nBjets     = nBjettmp;
      hggtmp.njetsTTH = goodTTHJetsIndex.size();
      hggtmp.jetIdxTTH.clear();
      for( size_t iji = 0; iji < goodTTHJetsIndex.size(); iji++ ) hggtmp.jetIdxTTH.push_back(goodTTHJetsIndex[iji]);
      
      
      // Select VH jets -----------------------------------  
      vector<int> goodVHJetsIndex = selectJetsVHLep( allGoodJetsIndex , selElecIndex, muonIndex ); // changed to selElecIndex 24 avril
      hggtmp.njetsVH = goodVHJetsIndex.size();
     
      // Select VH had jets -----------------------------------  
      vector<int> goodVHhadJetsIndex = selectJetsVHHad( allGoodJetsIndex );
      hggtmp.njetsVHhad = goodVHhadJetsIndex.size();
      hggtmp.jetIdxVHhad.clear();
      for( size_t iji = 0; iji < goodVHhadJetsIndex.size(); iji++ ) hggtmp.jetIdxVHhad.push_back(goodVHhadJetsIndex[iji]);
      
      hggcandidates.push_back(hggtmp);     

    }     
  }

  
  if( hggcandidates.size() > 0 ) {
    compareHggCandbyEt sortH;
    sort( hggcandidates.begin(), hggcandidates.end(), sortH );
  } 
  
  return hggcandidates;
}


//*******************************************************************************************************


void xAna::analyzeHggCandidates( HggEvtCandidate& hggevt, std::vector<HggEvtCandidate>& hggevents, std::vector<category> allCategories,  MassResolution massresoc )
{
   vector <HggEvtCandidate> hggevtPassing;
   HggEvtCandidate goodHggEvt;

   for( size_t icat = 0; icat < allCategories.size(); icat++ ) { 

     hggevtPassing.clear();
     
     for( size_t ievt = 0; ievt < hggevents.size(); ievt++ ) {
       
       if( !(allCategories[icat].isGoodPhotonPair(hggevents[ievt])) ) continue;
       hggevtPassing.push_back(hggevents[ievt]);

     } 
     
     if( hggevtPassing.size() == 0 ) continue;    // go to lower category if no candidate passes
     else goodHggEvt = hggevtPassing[0];       
     
     float diphotonMVA = evaluateDiPhotonMVA( goodHggEvt.gammaIdx_lead, goodHggEvt.gammaIdx_trail, goodHggEvt.gammaK_lead, goodHggEvt.gammaK_trail, massresoc );
     goodHggEvt.diphotonMva = diphotonMVA;

     if( !isEvtInCategory(goodHggEvt, allCategories[icat] ) ) continue;

     if( allCategories[icat].istag ) {
       goodHggEvt.category    = allCategories[icat].name; 
       goodHggEvt.subcategory = allCategories[icat].subcat; 
     }
     else{
       goodHggEvt.untagsubcat    = allCategories[icat].subcat; 
     }
     
     
     hggevt = goodHggEvt;
     return;
   }
   
}

//********************************************************************************************


HggEvtCandidate xAna::sortHggCandidates( std::vector<HggEvtCandidate>& hggevents, std::vector<category> tagCategories, std::vector<category> untagCategories,  MassResolution massresoc )
{  

  HggEvtCandidate goodHggEvt;
  analyzeHggCandidates( goodHggEvt, hggevents, tagCategories, massresoc );

  if( goodHggEvt.gammaIdx_lead != -1 ) {
  
    for( size_t ica = 0; ica < untagCategories.size(); ica++ ) {
   
      if( !isEvtInCategory( goodHggEvt, untagCategories[ica] ) ) continue;
      goodHggEvt.untagsubcat = untagCategories[ica].subcat; 
      break;     
    }
  }
  else{
    analyzeHggCandidates( goodHggEvt, hggevents, untagCategories, massresoc );
  }
  return goodHggEvt;
}













