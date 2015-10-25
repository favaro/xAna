#ifndef hggCandidate_hh__
#define hggCandidate_hh__

#include <TLorentzVector.h>

class HggEvtCandidate {

 public:
  int            vertexIdx;
  int            gammaIdx_lead;
  int            gammaIdx_trail;
  TLorentzVector gammaK_lead;
  TLorentzVector gammaK_trail;
  TLorentzVector corMet;        
  float          gammaR9_lead;
  float          gammaR9_trail;
  double         diphotonMva;
  int            njets;
  int            njetsTTH;    
  int            njetsVH;     
  int            njetsVHhad;  
  int            nBjets;
  std::vector<int>    jetIdx;
  std::vector<int>    jetIdxTTH;      
  std::vector<int>    jetIdxVHhad;    
  int            nmu;
  int            nelec;
  std::vector<int>    muIdx;
  std::vector<int>    eleIdx;
  int            nmuVHDilep;   
  int            nelecVHDilep; 
  std::vector<int>    muIdxVHDilep;
  std::vector<int>    eleIdxVHDilep;

  TString        category;
  int            subcategory;
  int            untagsubcat;

 inline HggEvtCandidate():
   vertexIdx(-1),
   gammaIdx_lead(-1),
   gammaIdx_trail(-1),
   gammaR9_lead(-1),
   gammaR9_trail(-1),
   diphotonMva(-999),
   njets(-1),
   njetsTTH(-1),
   njetsVH(-1),
   njetsVHhad(-1),
   nBjets(-1),
   nmu(-1),
   nelec(-1),
   nmuVHDilep(-1),
   nelecVHDilep(-1),
   category("none"),
   subcategory(-1),
   untagsubcat(-1)
  {}
  
};

struct compareHggCandbyEt {
  inline bool operator()( const HggEvtCandidate hgg1, const HggEvtCandidate hgg2 ) const {
    
    float mass1 = hgg1.gammaK_lead.Et() + hgg1.gammaK_trail.Et();
    float mass2 = hgg2.gammaK_lead.Et() + hgg2.gammaK_trail.Et();
    
    return  mass1 > mass2;
  }
};

#endif
