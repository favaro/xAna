#include "interface/xAna.hh"
#include "interface/vertexUtils.hh"
#include "interface/MinitreeOut.hh"

#include "VertexAnalysis/PhotonInfo.h"


int  xAna::matchPhotonToConversion(int lpho) {
  
  if(DoDebugEvent)  _xcheckTextFile << " ------- xAna::matchPhotonToConversion (pho=" << lpho << " ) --------" << endl; 
  
  TVector3 pho_xyz; 
  pho_xyz.SetXYZ((*phoCaloPos_x)[lpho],
		 (*phoCaloPos_y)[lpho],
		 (*phoCaloPos_z)[lpho]);
  
  float dRMin = 999.;
  int iMatch=-1;     
  
  //// defining a good conversion
  float probMin = 1e-6;
  for(int iconv=0; iconv<nConv; iconv++) {
    TVector3 qconv; 
    if     ( (*convNTracks)[iconv] == 2 ) qconv.SetXYZ( (*convRefittedMomentum_x)[iconv],
							(*convRefittedMomentum_y)[iconv],
							(*convRefittedMomentum_z)[iconv] );
    else if( (*convNTracks)[iconv] == 1 ) qconv.SetXYZ( (*phoConvPairMomentum_x)[iconv],
							(*phoConvPairMomentum_y)[iconv],
							(*phoConvPairMomentum_z)[iconv] ); 
      


    TVector3 vtx_xyz( (*convVtx_x)[iconv], (*convVtx_y)[iconv], (*convVtx_z)[iconv] );
    TVector3 qpho = pho_xyz - vtx_xyz;        
    float dR = qpho.DeltaR( qconv ); 
    if(DoDebugEvent)    {
      _xcheckTextFile  << " check conv[" << iconv << "]: pt = " << qconv.Pt() << endl
		       << "              : nTrk = " <<  (*convNTracks)[iconv]<< endl
		       << "              : vtxValid? = " << (*convValidVtx)[iconv] << endl
		       << "              : chi2 = " << (*convChi2Probability)[iconv] << endl
		       << "              : z   = " << vtx_xyz.Z() << " vs PV: " << (*vtxbs_z)[0]<< endl
		       << "              : phoDR = " << dR << endl;	
    }

    if( qconv.Pt() < 10 ) continue;

    if( _vtxAlgoParams.useAllConversions == 1 && (*convNTracks)[iconv] != 1 ) continue;
    if( _vtxAlgoParams.useAllConversions == 2 && (*convNTracks)[iconv] != 2 ) continue;
    if( _vtxAlgoParams.useAllConversions == 3 && (*convNTracks)[iconv] != 1 && (*convNTracks)[iconv] != 2 ) continue;
    if( (*convNTracks)[iconv] == 2 &&
	! ( (*convValidVtx)[iconv] && (*convChi2Probability)[iconv] > probMin ) ) continue;
    if(DoDebugEvent)  _xcheckTextFile  << "              : pass quality cut. " << endl;

    if ( dR < dRMin )  {
      dRMin=dR;
      iMatch=iconv;
    }
  }
  

  if ( dRMin < 0.1 )  {
    if(DoDebugEvent)    _xcheckTextFile  << "       =====>  matched conversion index " << iMatch << endl;
    return iMatch;
  } 
 
  return -1;
}


PhotonInfo xAna::fillPhotonInfo( int iPho ) {

    bool pho_isEB_1 = false;
    if(fabs((*phoSCEta)[iPho]) < 1.4442) pho_isEB_1 = true;

    int iConv1 = -1;
    if( _vtxAlgoParams.useAllConversions > 0 ) 
      iConv1 = matchPhotonToConversion(iPho);

    if(DoDebugEvent)  _xcheckTextFile << "iPho = " << iPho << "  iConv1 = " << iConv1 << endl;
	
    TVector3 qpho; qpho.SetXYZ( (*phoCaloPos_x)[iPho], (*phoCaloPos_y)[iPho], (*phoCaloPos_z)[iPho] );
    TVector3 bs  ; bs  .SetXYZ( bspotPos[0], bspotPos[1], bspotPos[2] );


    if(iConv1>=0) {
      TVector3 vtxconv; vtxconv.SetXYZ((*convVtx_x)[iConv1],(*convVtx_y)[iConv1],(*convVtx_z)[iConv1]);
   
      if(DoDebugEvent)
	_xcheckTextFile << " fillPhotonInfo: conversion match, convZofPrimVtxFromTrks["<<iConv1<<"] = " 
			<<  (*convZofPrimVtxFromTrks)[iConv1] << " : defVtx = " << (*vtxbs_z)[0]<< endl;
      
      PhotonInfo pho(iPho, qpho, bs, vtxconv,
		    (*convNTracks)[iConv1] == 1 ? 
		     TVector3((*convPairMomentum_x)[iConv1], (*convPairMomentum_y)[iConv1], (*convPairMomentum_z)[iConv1]):
		     TVector3((*convRefittedMomentum_x)[iConv1], (*convRefittedMomentum_y)[iConv1], (*convRefittedMomentum_z)[iConv1]),
		    (*phoE)[iPho],  pho_isEB_1,
		    (*convNTracks)[iConv1], (*convValidVtx)[iConv1], (*convChi2Probability)[iConv1], (*convEoverP)[iConv1] ); 
      
      return pho;      
    } 

    float phoConvV_X = 0;//phoConvVtx[iPho][0];
    float phoConvV_Y = 0;//phoConvVtx[iPho][1];
    float phoConvV_Z = 0;//phoConvVtx[iPho][2];
    
    // PhotonInfo pho(iPho, qpho, bs, TVector3(phoConvV_X,phoConvV_Y,phoConvV_Z),
    // 		   TVector3(phoConvRefittedMomentum_x[iPho], phoConvRefittedMomentum_y[iPho], phoConvRefittedMomentum_z[iPho]),
    // 		  (*phoE)[iPho], pho_isEB_1,
    // 		   (*phoConvNTrks)[iPho]   , phoConvValidVtx[iPho],
    // 		   phoConvChi2Prob[iPho], (*phoConvEoverP)[iPho] );
    PhotonInfo pho(iPho, qpho, bs, TVector3(phoConvV_X,phoConvV_Y,phoConvV_Z),
		   TVector3((*phoConvRefittedMomentum_x)[iPho], (*phoConvRefittedMomentum_y)[iPho], (*phoConvRefittedMomentum_z)[iPho]),
		  (*phoE)[iPho], pho_isEB_1,
		   (*phoConvNTrks)[iPho]   , 0,
		   -1, (*phoConvEoverP)[iPho] ); 

    return pho;
}


vector<int>  xAna::getSelectedVertex(int iPho1, int iPho2, bool useMva) {

  PhotonInfo pho1 = fillPhotonInfo(iPho1);
  PhotonInfo pho2 = fillPhotonInfo(iPho2);

  // This is a temp fix to run on ggNtuples that did not fill trk info correctly
  // It is ok for the DijetTag only (Vladimir Rekovic 19Nov2013)
  bool PROBLEM_w_ggNuples = false;
  if(PROBLEM_w_ggNuples) {
    vector<int> temp_result;
    temp_result.push_back(0);
    return temp_result;
  }

  ///// FC: try to fill it differently since it does not give the same results as gglobe
  vector<float> vtx__x,vtx__y,vtx__z;
  vector<int>   vtx_ntks;
  vector<vector<unsigned short> > vtx_tkind;
  vector<vector<float> > vtx_tkweight;
  vector<float> tk_px, tk_py, tk_pz, tk_pterr;
  vector<int>   tk_quality;
  
  /// fill tracks
  for( int t = 0; t < nTrk; ++t) {
    tk_px.push_back( (*trkP_x)[t] );
    tk_py.push_back( (*trkP_y)[t] );
    tk_pz.push_back( (*trkP_z)[t] );
    tk_quality.push_back( (*trkQuality)[t] );
    tk_pterr  .push_back( (*trkPtErr)[t] );
  }
  
  /// fix in case the track branch does not exist
  bool trkBranchExist = fChain->GetBranchStatus("nTrk");

  /// fill vtx
  vtx_tkind.resize(nVtxBS);
  vtx_tkweight.resize(nVtxBS);

  for( int v = 0; v < nVtxBS; ++v) {
    vtx__x.push_back( (*vtxbs_x)[v] );
    vtx__y.push_back( (*vtxbs_y)[v] );
    vtx__z.push_back( (*vtxbs_z)[v] );

    unsigned tksize = 0;
    if( trkBranchExist ) tksize = (*vtxbsTkIndex)[v].size();


    vtx_ntks.push_back( tksize );

    /// thenindex should be already ok...    
    for( unsigned it = 0 ; it <  tksize; it++ ) {
      vtx_tkind[v].push_back((unsigned short)(*vtxbsTkIndex)[v][it]);
      vtx_tkweight[v].push_back( (*vtxbsTkWeight)[v][it]);
    }
  }
  
  ggVertexInfo vinfo(vtx__x, vtx__y,vtx__z,
		     vtx_ntks,
		     tk_px,tk_py,tk_pz,
		     tk_pterr, tk_quality,
		     vtx_tkind, vtx_tkweight
		     ) ;


  _vtxAna->clear();
  _vtxAna->analyze(vinfo, pho1, pho2);


  //int result = -1;
  vector<int> result;
  result.push_back(-1);
  if (useMva) { 
    vector<int> rankresults = _vtxAna->rank(*_perVtxReader, _perVtxMvaMethod);
    vertMVA  = _vtxAna->perEventMva(*_perEvtReader ,_perEvtMvaMethod,rankresults);
    vertProb = _vtxAna->vertexProbability(vertMVA,nVtxBS);

    if( ! trkBranchExist ) {
      vertProb = 0.7; 
      if( rankresults.size() > 0 ) rankresults[0];
      else  rankresults.push_back(0);
    }

    if( _minitree->addSyncVariables ) {
      for( unsigned iv = 0; iv < 3; iv++ ) {
	_minitree->vtxId[ iv] = iv < rankresults.size() ? rankresults[iv] : -1;
	_minitree->vtxMva[iv] = iv < rankresults.size() ?_vtxAna->mva(rankresults[iv]) : -2;
	_minitree->vtxDz[iv]  = iv < rankresults.size() ?_vtxAna->vertexz(rankresults[iv]) - _vtxAna->vertexz(rankresults[0]) : -999;
      }
      _minitree->vtxPtBal      = _vtxAna->ptbal(rankresults[0]);
      _minitree->vtxAsym       = _vtxAna->ptasym(rankresults[0]);
      _minitree->vtxLogSumPt2  = _vtxAna->logsumpt2(rankresults[0]);
      _minitree->vtxPullToConv = _vtxAna->pulltoconv(rankresults[0]);
      _minitree->vtxNConv      = _vtxAna->nconv(rankresults[0]);
 
      if( DoDebugEvent ) {
	_xcheckTextFile << "  vertex analysis   nConversion: " << _minitree->vtxNConv << endl;
	_xcheckTextFile << "  vertex analysis:  pullToConv: " << _minitree->vtxPullToConv << endl;
      }
   }

    return rankresults;
  }


  return result; 
}

