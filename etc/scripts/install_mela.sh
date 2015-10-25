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

rm -f $xAnaExtLIBS/JHUGenMELA
ln -s JHUGenerator.v4.0.1/JHUGenMELA/ggZZ_MCFM $xAnaExtLIBS/JHUGenMELA

rm -f $xAnaExtLIBS/MCFM-6.6
ln -s JHUGenerator.v4.0.1/MCFM-6.6 $xAnaExtLIBS/MCFM-6.6


##### compile all external libraries
cd $baseDir/$xAnaExtLIBS/VertexAnalysis
make -j 4

cd $baseDir/$xAnaExtLIBS/JetMETObjects
make -j 4

cd $baseDir/$xAnaExtLIBS/GBRLikelihood
ln -s $baseDir/MTools .
source MTools/setupmtools.sh
make -j 4


cd $baseDir/$xAnaExtLIBS/JHUGenMELA
make clean
make -j 4

cd $baseDir/$xAnaExtLIBS/MCFM-6.6
chmod +x Install
./Install
cp $baseDir/$xAnaExtLIBS/JHUGenMELA/gg_ZZ_int.f  src/HZZ/
make -j 4 mcfmlib


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
ln -s ../JHUGenMELA/libME.so $xAnaExtLIBS/lib/libME.so
cd $xAnaExtLIBS/lib
for libmcfm in `find  ../MCFM-6.6/ -name "*.a"`; do
    ln -s $libmcfm `basename $libmcfm`
done
cd  $baseDir

ln -s ../VertexAnalysis/interface $xAnaExtLIBS/interface/VertexAnalysis
ln -s ../GBRLikelihood/interface  $xAnaExtLIBS/interface/GBRLikelihood
ln -s ../JetMETObjects/interface  $xAnaExtLIBS/interface/JetMETObjects
ln -s ../JHUGenMELA               $xAnaExtLIBS/interface/JHUGenMELA
make -j 4