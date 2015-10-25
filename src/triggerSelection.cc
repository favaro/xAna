#include "interface/xAna.hh"

#include <iostream>
using namespace std;

///---------------------------------------------------------------------------------///
/// -------------------------- 8TeV Diphoton Trigger -------------------------------///
///---------------------------------------------------------------------------------///
bool xAna::PassTriggerSelection2012( void ) {
  int runNum = run;
  Int_t passHLT = 1;
  
  if (mode_ == -1) { // 2012 data
   
    if       ( runNum >= 190456 && runNum <= 194269 ) {

      //if(HLTIndex[0]==-1 || HLTIndex[8]==-1 )
      //	{
      //	  cout <<" run "<<  runNum <<"  event "<< event   <<" HLTIndex[0]  "<<  HLTIndex[0] <<"     HLTIndex[8]  "<<  HLTIndex[8] <<endl;
      //}
      //else{
      if( HLTIndex[0]!=-1 && HLTIndex[8]!=-1 ){
	if ( HLT[HLTIndex[0]] == 0 &&  HLT[HLTIndex[8]] == 0) passHLT=0;
      }

    } else if( runNum >= 194270 && runNum < 203894 ) {
      
      // JM: remove cout
      //if(HLTIndex[0]==-1 || HLTIndex[7]==-1 ||  HLTIndex[8]==-1)
      //	{
      //	  cout <<" run "<<  runNum <<"  event "<< event   <<" HLTIndex[0]  "<<  HLTIndex[0] <<
      //	    " HLTIndex[7]  "<<  HLTIndex[7] <<"  HLTIndex[8]  "<<  HLTIndex[8] <<endl;
      //	}
      //else {	
      
      if(HLTIndex[0]!=-1 && HLTIndex[7]!=-1 &&  HLTIndex[8]!=-1){
	if ( HLT[HLTIndex[0]] == 0 && HLT[HLTIndex[7]] == 0 &&  HLT[HLTIndex[8]] == 0 ) passHLT=0;  
      }
      
    } else if( runNum >= 203894 ) {
      
      // JM: remove cout
      // if( HLTIndex[7]==-1 ||  HLTIndex[8]== -1 ){ 
      //	cout <<" run "<<  runNum <<"  event "<< event   <<" HLTIndex[0]  "<<  HLTIndex[0] <<
      //	  " HLTIndex[7]  "<<  HLTIndex[7] <<"  HLTIndex[8]  "<<  HLTIndex[8] <<endl;
      //} else {	
      if( HLTIndex[7]!=-1 &&  HLTIndex[8]!= -1 ){
	if ( HLT[HLTIndex[7]] == 0 &&  HLT[HLTIndex[8]] == 0 ) passHLT=0;  
      }
    }
    
  } else {   // is MC
    passHLT=1;   //==================== NO HLT applied to MC ===============================
  }
  
  return passHLT;
  
}



///---------------------------------------------------------------------------------///
/// -------------------------- 7TeV Diphoton Trigger -------------------------------///
///---------------------------------------------------------------------------------///

bool xAna::PassTriggerSelection2011( void ) {
  int runNum = run;
  Int_t passHLT = 1;
  if (mode_ == -1) { // 2011 data
    if(runNum>=160404 && runNum <= 161176){
      if(HLTIndex[21]==-1 ) {cout <<" run "<<  runNum <<"  event "<< event   <<" HLTIndex[21]  "<<  HLTIndex[21] <<endl;}
      else {
	if ( HLT[HLTIndex[21]] == 0)passHLT=0;
      }
    } else if(runNum>=161216  && runNum <= 165208){
      if(HLTIndex[21]==-1 || HLTIndex[45]==-1 || HLTIndex[46] ==-1 || HLTIndex[23] ==-1 ){
	// cout <<" run "<<  runNum <<"  event "<< event   <<" HLTIndex[21]  "<<  HLTIndex[21]  
	// 	 <<"  HLTIndex[45] "<<  HLTIndex[45]  <<"  HLTIndex[46] "<<   HLTIndex[46] <<"  HLTIndex[23] "<<   HLTIndex[23] <<endl;
      }else{
	if ( HLT[HLTIndex[21]] == 0 && HLT[HLTIndex[23]] == 0  &&   HLT[HLTIndex[45]] == 0 && HLT[HLTIndex[46]] == 0 ) {
	  passHLT=0;
	}
      }
    } else if(runNum>=165970  && runNum <= 173198){
      if ( HLTIndex[21]==-1 || HLTIndex[45]==-1 || HLTIndex[46] ==-1 || HLTIndex[24] ==-1
	   || HLTIndex[54]==-1 || HLTIndex[51]==-1 || HLTIndex[56] ==-1 || HLTIndex[57] ==-1 )
	{
	  //cout <<" run "<<  run <<"  event "<< event   <<"  HLTIndex[21] "<<  HLTIndex[21] <<"  HLTIndex[45] "<<  HLTIndex[45] << "  HLTIndex[46] "<<  HLTIndex[46] <<"  HLTIndex[24] "<<  HLTIndex[24] <<
	  //                                               "  HLTIndex[54] "<<  HLTIndex[54] <<"  HLTIndex[55] "<<  HLTIndex[55] << "  HLTIndex[56] "<<  HLTIndex[56] <<"  HLTIndex[57] "<<  HLTIndex[57] << endl;
	  //********************** missing HLTIndex[55] ***************************************************************************//
	  // if ( HLTIndex[21]==-1 || HLTIndex[45]==-1 || HLTIndex[46] ==-1 || HLTIndex[24] ==-1
	  //   || HLTIndex[54]==-1 ||                     HLTIndex[56] ==-1 || HLTIndex[57] ==-1 )
	  // {
	  
	  if ( HLTIndex[51]==-1){
	    // cout << " run "<<  runNum <<"  event "<< event   <<"  HLTIndex[21] "<<  HLTIndex[21] 
	    // 	   << "  HLTIndex[45] "<<  HLTIndex[45] << "  HLTIndex[46] "<<  HLTIndex[46] <<"  HLTIndex[24] "<<  HLTIndex[24] 
	    // 	   << "  HLTIndex[54] "<<  HLTIndex[54] <<"  HLTIndex[55] "<<  HLTIndex[55] << "  HLTIndex[56] "<<  HLTIndex[56] 
	    //           << "  HLTIndex[57] "<<  HLTIndex[57] << endl;
              if ( HLT[HLTIndex[21]] == 0 && HLT[HLTIndex[45]] == 0  &&   HLT[HLTIndex[46]] == 0 && HLT[HLTIndex[24]] == 0
		   &&  HLT[HLTIndex[54]] == 0 && HLT[HLTIndex[55]] == 0  &&   HLT[HLTIndex[56]] == 0 && HLT[HLTIndex[57]] == 0 )
		{
		  passHLT=0;
		}
	  } else {
	    if ( HLT[HLTIndex[21]] == 0 && HLT[HLTIndex[45]] == 0  &&   HLT[HLTIndex[46]] == 0 && HLT[HLTIndex[24]] == 0
		 &&  HLT[HLTIndex[54]] == 0 && HLT[HLTIndex[51]] == 0  &&   HLT[HLTIndex[56]] == 0 && HLT[HLTIndex[57]] == 0 )
	      {
		passHLT=0;
	      }
	  }
	  //***********************************************************************************************************************//
	} else {
	if ( HLT[HLTIndex[21]] == 0 && HLT[HLTIndex[45]] == 0  &&   HLT[HLTIndex[46]] == 0 && HLT[HLTIndex[24]] == 0
             &&  HLT[HLTIndex[54]] == 0 && HLT[HLTIndex[55]] == 0  &&   HLT[HLTIndex[56]] == 0 && HLT[HLTIndex[57]] == 0 )
	  {
	    passHLT=0;
	  }
      }
      
    } else if(runNum>= 173236 && runNum <178420 ){
      if ( HLTIndex[48]==-1 || HLTIndex[49]==-1 || HLTIndex[50] ==-1 || HLTIndex[24] ==-1
	   || HLTIndex[54]==-1 || HLTIndex[55]==-1 || HLTIndex[56] ==-1 || HLTIndex[57] ==-1 )
	{
	  //cout <<" run "<<  run <<"  event "<< event   <<" HLTIndex[48]  "<<  HLTIndex[48]  <<"  HLTIndex[49] "<<  HLTIndex[49]  <<"  HLTIndex[50] "<<   HLTIndex[50] <<"  HLTIndex[24] "<<   HLTIndex[24] <<
	  //                                               " HLTIndex[54]  "<<  HLTIndex[54]  <<"  HLTIndex[55] "<<  HLTIndex[55]  <<"  HLTIndex[56] "<<   HLTIndex[56] <<"  HLTIndex[57] "<<   HLTIndex[57] << endl;
	  //********************** missing HLTIndex[55] ***************************************************************************//
	  if ( HLTIndex[48]==-1 || HLTIndex[49]==-1 || HLTIndex[50] ==-1 || HLTIndex[24] ==-1
	       || HLTIndex[54]==-1 ||                     HLTIndex[56] ==-1 || HLTIndex[57] ==-1 )
            {
              // cout <<" run "<<  runNum <<"  event "<< event   <<" HLTIndex[48]  "<<  HLTIndex[48]  <<"  HLTIndex[49] "<<  HLTIndex[49]  <<"  HLTIndex[50] "<<   HLTIndex[50] 
	      // 	   <<"  HLTIndex[24] "<<   HLTIndex[24] << " HLTIndex[54]  "<<  HLTIndex[54]  <<"  HLTIndex[55] "<<  HLTIndex[55]  <<"  HLTIndex[56] "
	      // 	   <<   HLTIndex[56] <<"  HLTIndex[57] "<<   HLTIndex[57] << endl;
            } else {
	    if ( HLT[HLTIndex[48]] == 0 && HLT[HLTIndex[49]] == 0  &&   HLT[HLTIndex[50]] == 0 && HLT[HLTIndex[24]] == 0
		 &&  HLT[HLTIndex[54]] == 0 &&                              HLT[HLTIndex[56]] == 0 && HLT[HLTIndex[57]] == 0 )
              {
                passHLT=0;
              }
	  }
	  //***********************************************************************************************************************//
	} else {
	if ( HLT[HLTIndex[48]] == 0 && HLT[HLTIndex[49]] == 0  &&   HLT[HLTIndex[50]] == 0 && HLT[HLTIndex[24]] == 0
             &&  HLT[HLTIndex[54]] == 0 && HLT[HLTIndex[55]] == 0  &&   HLT[HLTIndex[56]] == 0 && HLT[HLTIndex[57]] == 0 )
	  {
	    passHLT=0;
	  }
      }
    } else if(runNum>=178420){
      if ( HLTIndex[59]==-1 || HLTIndex[60]==-1 || HLTIndex[61] ==-1 || HLTIndex[62] ==-1
	   || HLTIndex[54]==-1 || HLTIndex[55]==-1 || HLTIndex[56] ==-1 || HLTIndex[57] ==-1 )
	{
	  //cout <<" run "<<  run <<"  event "<< event   <<" HLTIndex[59]  "<<  HLTIndex[59]  <<"  HLTIndex[60] "<<  HLTIndex[60]  <<"  HLTIndex[61] "<<   HLTIndex[61] <<"  HLTIndex[62] "<<   HLTIndex[62] <<
	  //                                               " HLTIndex[54]  "<<  HLTIndex[54]  <<"  HLTIndex[55] "<<  HLTIndex[55]  <<"  HLTIndex[56] "<<   HLTIndex[56] <<"  HLTIndex[57] "<<   HLTIndex[57] << endl;
	  //********************** missing HLTIndex[55] ***************************************************************************//
	  if ( HLTIndex[59]==-1 || HLTIndex[60]==-1 || HLTIndex[61] ==-1 || HLTIndex[62] ==-1
	       || HLTIndex[54]==-1 ||                     HLTIndex[56] ==-1 || HLTIndex[57] ==-1 )
            {
	      // cout <<" run "<<  runNum <<"  event "<< event   <<" HLTIndex[59]  "<<  HLTIndex[59]  <<"  HLTIndex[60] "<<  HLTIndex[60]  
	      // 	   <<"  HLTIndex[61] "<<   HLTIndex[61] <<"  HLTIndex[62] "<<   HLTIndex[62] << " HLTIndex[54]  "<<  HLTIndex[54]  <<"  HLTIndex[55] "
	      // 	   <<  HLTIndex[55]  <<"  HLTIndex[56] "<<   HLTIndex[56] <<"  HLTIndex[57] "<<   HLTIndex[57] << endl;
            } else {
	    if ( HLT[HLTIndex[59]] == 0 && HLT[HLTIndex[60]] == 0  &&   HLT[HLTIndex[61]] == 0 && HLT[HLTIndex[62]] == 0
		 &&  HLT[HLTIndex[54]] == 0 &&                              HLT[HLTIndex[56]] == 0 && HLT[HLTIndex[57]] == 0 )
              {
                passHLT=0;
              }
	  }
	  //***********************************************************************************************************************//
	}else{
	if ( HLT[HLTIndex[59]] == 0 && HLT[HLTIndex[60]] == 0  &&   HLT[HLTIndex[61]] == 0 && HLT[HLTIndex[62]] == 0
             &&  HLT[HLTIndex[54]] == 0 && HLT[HLTIndex[55]] == 0  &&   HLT[HLTIndex[56]] == 0 && HLT[HLTIndex[57]] == 0 )
	  {
	    passHLT=0;
	  }
      }
      //passHLT=1;
    }
  } else {   // is MC
    passHLT=1;   //==================== NO HLT applied to MC ===============================
  }
  
  return passHLT;

}
