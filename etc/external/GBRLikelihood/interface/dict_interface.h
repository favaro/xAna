/********************************************************************
* interface/dict_interface.h
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************************/
#ifdef __CINT__
#error interface/dict_interface.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#define G__PRIVATE_GVALUE
#include "G__ci.h"
#include "FastAllocString.h"
extern "C" {
extern void G__cpp_setup_tagtabledict_interface();
extern void G__cpp_setup_inheritancedict_interface();
extern void G__cpp_setup_typetabledict_interface();
extern void G__cpp_setup_memvardict_interface();
extern void G__cpp_setup_globaldict_interface();
extern void G__cpp_setup_memfuncdict_interface();
extern void G__cpp_setup_funcdict_interface();
extern void G__set_cpp_environmentdict_interface();
}


#include "TObject.h"
#include "TMemberInspector.h"
#include "interface/HybridGBRForest.h"
#include "interface/HybridGBRForestD.h"
#include "interface/HybridGBRTree.h"
#include "interface/HybridGBRTreeD.h"
#include "interface/RooCBExp.h"
#include "interface/RooCBFast.h"
#include "interface/RooDoubleCBFast.h"
#include "interface/RooGaussianFast.h"
#include "interface/RooHybridBDTAutoPdf.h"
#include "interface/RooRevCBFast.h"
#include <algorithm>
namespace std { }
using namespace std;

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__dict_interfaceLN_TClass;
extern G__linked_taginfo G__dict_interfaceLN_TBuffer;
extern G__linked_taginfo G__dict_interfaceLN_TMemberInspector;
extern G__linked_taginfo G__dict_interfaceLN_TObject;
extern G__linked_taginfo G__dict_interfaceLN_TNamed;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEunsignedsPcharcOallocatorlEunsignedsPchargRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEfloatcOallocatorlEfloatgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEdoublecOallocatorlEdoublegRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_string;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_HybridGBRTree;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEintcOallocatorlEintgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEintcOallocatorlEintgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEpairlEfloatcOfloatgRcOallocatorlEpairlEfloatcOfloatgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEpairlEfloatcOfloatgRcOallocatorlEpairlEfloatcOfloatgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEvectorlEpairlEfloatcOfloatgRcOallocatorlEpairlEfloatcOfloatgRsPgRsPgRcOallocatorlEvectorlEpairlEfloatcOfloatgRcOallocatorlEpairlEfloatcOfloatgRsPgRsPgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEvectorlEpairlEfloatcOfloatgRcOallocatorlEpairlEfloatcOfloatgRsPgRsPgRcOallocatorlEvectorlEpairlEfloatcOfloatgRcOallocatorlEpairlEfloatcOfloatgRsPgRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_HybridGBRForest;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEHybridGBRTreecOallocatorlEHybridGBRTreegRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEHybridGBRTreecOallocatorlEHybridGBRTreegRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEvectorlEHybridGBRTreecOallocatorlEHybridGBRTreegRsPgRcOallocatorlEvectorlEHybridGBRTreecOallocatorlEHybridGBRTreegRsPgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEvectorlEHybridGBRTreecOallocatorlEHybridGBRTreegRsPgRcOallocatorlEvectorlEHybridGBRTreecOallocatorlEHybridGBRTreegRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_HybridGBRTreeD;
extern G__linked_taginfo G__dict_interfaceLN_HybridGBRForestD;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEHybridGBRTreeDcOallocatorlEHybridGBRTreeDgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEHybridGBRTreeDcOallocatorlEHybridGBRTreeDgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEvectorlEHybridGBRTreeDcOallocatorlEHybridGBRTreeDgRsPgRcOallocatorlEvectorlEHybridGBRTreeDcOallocatorlEHybridGBRTreeDgRsPgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEvectorlEHybridGBRTreeDcOallocatorlEHybridGBRTreeDgRsPgRcOallocatorlEvectorlEHybridGBRTreeDcOallocatorlEHybridGBRTreeDgRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_iteratorlEbidirectional_iterator_tagcOTObjectmUcOlongcOconstsPTObjectmUmUcOconstsPTObjectmUaNgR;
extern G__linked_taginfo G__dict_interfaceLN_RooPrintable;
extern G__linked_taginfo G__dict_interfaceLN_listlEdoublecOallocatorlEdoublegRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_listlEstringcOallocatorlEstringgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_RooAbsArg;
extern G__linked_taginfo G__dict_interfaceLN_RooArgSet;
extern G__linked_taginfo G__dict_interfaceLN_RooArgList;
extern G__linked_taginfo G__dict_interfaceLN_setlEstringcOlesslEstringgRcOallocatorlEstringgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_maplEstringcOTObjArraymUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTObjArraymUgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_TTree;
extern G__linked_taginfo G__dict_interfaceLN_RooArgProxy;
extern G__linked_taginfo G__dict_interfaceLN_RooListProxy;
extern G__linked_taginfo G__dict_interfaceLN_RooRealProxy;
extern G__linked_taginfo G__dict_interfaceLN_maplEstringcOstringcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOstringgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_setlEpairlERooAbsArgmUcORooAbsArgmUgRcOlesslEpairlERooAbsArgmUcORooAbsArgmUgRsPgRcOallocatorlEpairlERooAbsArgmUcORooAbsArgmUgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_RooAbsReal;
extern G__linked_taginfo G__dict_interfaceLN_dequelERooAbsCachemUcOallocatorlERooAbsCachemUgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_RooDataSet;
extern G__linked_taginfo G__dict_interfaceLN_maplERooAbsArgmUcOTRefArraymUcOlesslERooAbsArgmUgRcOallocatorlEpairlERooAbsArgmUsPconstcOTRefArraymUgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_TVectorTlEfloatgR;
extern G__linked_taginfo G__dict_interfaceLN_TVectorTlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_RooRealVar;
extern G__linked_taginfo G__dict_interfaceLN_vectorlERooCurvemUcOallocatorlERooCurvemUgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlERooCurvemUcOallocatorlERooCurvemUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_RooAbsPdf;
extern G__linked_taginfo G__dict_interfaceLN_maplEstringcOTH1mUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOTH1mUgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_maplEstringcORooDataHistmUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcORooDataHistmUgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_maplEstringcORooDataSetmUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcORooDataSetmUgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_maplEstringcORooAbsDatamUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcORooAbsDatamUgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_listlERooAbsRealcLcLEvalErrorcOallocatorlERooAbsRealcLcLEvalErrorgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_pairlEstringcOlistlERooAbsRealcLcLEvalErrorcOallocatorlERooAbsRealcLcLEvalErrorgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_maplEconstsPRooAbsArgmUcOpairlEstringcOlistlERooAbsRealcLcLEvalErrorcOallocatorlERooAbsRealcLcLEvalErrorgRsPgRsPgRcOlesslEconstsPRooAbsArgmUgRcOallocatorlEpairlEconstsPRooAbsArgmUsPconstcOpairlEstringcOlistlERooAbsRealcLcLEvalErrorcOallocatorlERooAbsRealcLcLEvalErrorgRsPgRsPgRsPgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_maplEintcOstringcOlesslEintgRcOallocatorlEpairlEconstsPintcOstringgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlERooMsgServicecLcLStreamConfigcOallocatorlERooMsgServicecLcLStreamConfiggRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlERooMsgServicecLcLStreamConfigcOallocatorlERooMsgServicecLcLStreamConfiggRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_dequelEvectorlERooMsgServicecLcLStreamConfigcOallocatorlERooMsgServicecLcLStreamConfiggRsPgRcOallocatorlEvectorlERooMsgServicecLcLStreamConfigcOallocatorlERooMsgServicecLcLStreamConfiggRsPgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_stacklEvectorlERooMsgServicecLcLStreamConfigcOallocatorlERooMsgServicecLcLStreamConfiggRsPgRcOdequelEvectorlERooMsgServicecLcLStreamConfigcOallocatorlERooMsgServicecLcLStreamConfiggRsPgRcOallocatorlEvectorlERooMsgServicecLcLStreamConfigcOallocatorlERooMsgServicecLcLStreamConfiggRsPgRsPgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_maplEstringcObasic_ostreamlEcharcOchar_traitslEchargRsPgRmUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcObasic_ostreamlEcharcOchar_traitslEchargRsPgRmUgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlERooNormSetCachecOallocatorlERooNormSetCachegRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlERooNormSetCachecOallocatorlERooNormSetCachegRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlERooAbsCacheElementmUcOallocatorlERooAbsCacheElementmUgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlERooAbsCacheElementmUcOallocatorlERooAbsCacheElementmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_RooCBExp;
extern G__linked_taginfo G__dict_interfaceLN_RooCBFast;
extern G__linked_taginfo G__dict_interfaceLN_RooDoubleCBFast;
extern G__linked_taginfo G__dict_interfaceLN_RooGaussianFast;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTBaselEfloatgR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTBaselEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TElementActionTlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TElementPosActionTlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTSymlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTRow_constlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTRowlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTColumn_constlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTDiag_constlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTFlat_constlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTSub_constlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTSparseRow_constlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTSparselEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTSparseDiag_constlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTColumnlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTDiaglEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTFlatlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTSublEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTSparseRowlEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_TMatrixTSparseDiaglEdoublegR;
extern G__linked_taginfo G__dict_interfaceLN_RooTreeConvert;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEstringcOallocatorlEstringgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEstringcOallocatorlEstringgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_RooNormPdf;
extern G__linked_taginfo G__dict_interfaceLN_RooRealConstraint;
extern G__linked_taginfo G__dict_interfaceLN_RooPowerLaw;
extern G__linked_taginfo G__dict_interfaceLN_RooCondAddPdf;
extern G__linked_taginfo G__dict_interfaceLN_RooPdfAddReal;
extern G__linked_taginfo G__dict_interfaceLN_RooGBRFunction;
extern G__linked_taginfo G__dict_interfaceLN_RooGBRTarget;
extern G__linked_taginfo G__dict_interfaceLN_vectorlERooAbsDatamUcOallocatorlERooAbsDatamUgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlERooAbsDatamUcOallocatorlERooAbsDatamUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEHybridGBREventmUcOallocatorlEHybridGBREventmUgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEHybridGBREventmUcOallocatorlEHybridGBREventmUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlERooArgListcOallocatorlERooArgListgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlERooArgListcOallocatorlERooArgListgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlERooArgSetcOallocatorlERooArgSetgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlERooArgSetcOallocatorlERooArgSetgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRcOallocatorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRcOallocatorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEvectorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRcOallocatorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRsPgRsPgRcOallocatorlEvectorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRcOallocatorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRsPgRsPgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEvectorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRcOallocatorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRsPgRsPgRcOallocatorlEvectorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRcOallocatorlEvectorlERooAbsRealmUcOallocatorlERooAbsRealmUgRsPgRsPgRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlETMatrixTlEdoublegRcOallocatorlETMatrixTlEdoublegRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlETMatrixTlEdoublegRcOallocatorlETMatrixTlEdoublegRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlETVectorTlEdoublegRcOallocatorlETVectorTlEdoublegRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlETVectorTlEdoublegRcOallocatorlETVectorTlEdoublegRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEvectorlEintcOallocatorlEintgRsPgRcOallocatorlEvectorlEintcOallocatorlEintgRsPgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEvectorlEintcOallocatorlEintgRsPgRcOallocatorlEvectorlEintcOallocatorlEintgRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_setlEpairlEintcOintgRcOlesslEpairlEintcOintgRsPgRcOallocatorlEpairlEintcOintgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_vectorlEsetlEpairlEintcOintgRcOlesslEpairlEintcOintgRsPgRcOallocatorlEpairlEintcOintgRsPgRsPgRcOallocatorlEsetlEpairlEintcOintgRcOlesslEpairlEintcOintgRsPgRcOallocatorlEpairlEintcOintgRsPgRsPgRsPgRsPgR;
extern G__linked_taginfo G__dict_interfaceLN_reverse_iteratorlEvectorlEsetlEpairlEintcOintgRcOlesslEpairlEintcOintgRsPgRcOallocatorlEpairlEintcOintgRsPgRsPgRcOallocatorlEsetlEpairlEintcOintgRcOlesslEpairlEintcOintgRsPgRcOallocatorlEpairlEintcOintgRsPgRsPgRsPgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__dict_interfaceLN_RooRevCBFast;

/* STUB derived class for protected member access */
