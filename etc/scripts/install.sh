#!/bin/bash

source MTools/setupmtools.sh
source etc/scripts/setup_xAna.sh
chmod +x etc/scripts/setup_xAna.sh
chmod +x etc/scripts/xAnaBatch.sh
chmod +x etc/scripts/cleanupDist.sh
chmod +x etc/scripts/getNumberOfEvents.py
baseDir=`pwd`


if [ ! -e $xAnaExtLIBS/lib/  ]; then
    mkdir -p $xAnaExtLIBS/lib/
fi

if [ ! -e $xAnaExtLIBS/interface/  ]; then
    mkdir -p $xAnaExtLIBS/interface/
fi

rm -f $xAnaExtLIBS/VertexAnalysis
ln -s 611_VertexAnalysis $xAnaExtLIBS/VertexAnalysis

cd $baseDir/$xAnaExtLIBS/VertexAnalysis
make -j 4

cd 
cd $baseDir/$xAnaExtLIBS/JetMETObjects
make -j 4

cd $baseDir/$xAnaExtLIBS/GBRLikelihood
ln -s $baseDir/MTools .
source MTools/setupmtools.sh
make -j 4

cd $baseDir

rm -f $xAnaExtLIBS/lib/libHiggsAnalysisGBRLikelihood.so*
rm -f $xAnaExtLIBS/lib/libh2gglobeVertexAnalysis.so*
rm -f $xAnaExtLIBS/lib/libJetMETObjects.so*
rm -f $xAnaExtLIBS/lib/libmcfm_6p6.so*
rm -f $xAnaExtLIBS/lib/libME.so*
rm -f $xAnaExtLIBS/lib/libEG.so*
rm -f $xAnaExtLIBS/interface/VertexAnalysis
rm -f $xAnaExtLIBS/interface/GBRLikelihood
rm -f $xAnaExtLIBS/interface/JetMETObjects
rm -f $xAnaExtLIBS/interface/JHUGenMELA 

ln -s ../GBRLikelihood/lib/libHiggsAnalysisGBRLikelihood.so $xAnaExtLIBS/lib/libHiggsAnalysisGBRLikelihood.so
ln -s ../GBRLikelihood/lib/libHiggsAnalysisGBRLikelihood.so $xAnaExtLIBS/lib/libHiggsAnalysisGBRLikelihood.so.0.0
ln -s ../VertexAnalysis/lib/libh2gglobeVertexAnalysis.so    $xAnaExtLIBS/lib/libh2gglobeVertexAnalysis.so
ln -s ../JetMETObjects/lib/libJetMETObjects.so              $xAnaExtLIBS/lib/libJetMETObjects.so

ln -s ../MCFM-6.6/obj/libmcfm_6p6.so $xAnaExtLIBS/lib/libmcfm_6p6.so 
ln -s ../JHUGenMELA/ggZZ_MCFM/libME.so $xAnaExtLIBS/lib/libME.so
ln -s ../JHUGenMELA/ggZZ_MCFM/libEG.so $xAnaExtLIBS/lib/libEG.so

ln -s ../VertexAnalysis/interface $xAnaExtLIBS/interface/VertexAnalysis
ln -s ../GBRLikelihood/interface  $xAnaExtLIBS/interface/GBRLikelihood
ln -s ../JetMETObjects/interface  $xAnaExtLIBS/interface/JetMETObjects
ln -s ../JHUGenMELA/ggZZ_MCFM     $xAnaExtLIBS/interface/JHUGenMELA
make -j 4