#ifndef category_hh_
#define category_hh_

#include "interface/hggCandidate.hh"

#include <vector>

class ConfigReader; 


typedef std::pair<TString,float> cut;
class category {

 public:

  category( TString catname, int subcat, bool istag );
  ~category();
  
  void  Init( const ConfigReader *config  );
  void  setCommonCuts( const ConfigReader *config );
  void  setSpecificCuts( const ConfigReader *config );
  bool  isGoodPhotonPair( HggEvtCandidate hggevent );
  

  TString name;
  int     subcat;
  int     istag;
  size_t  ncuts;

  float cut_pt_lead, cut_pt_trail;
  float cut_ptOm_lead, cut_ptOm_trail;
  float cut_mgg;
  std::vector<cut> specificCuts;
  std::vector<TString> specificCuts_side;
};


#endif
