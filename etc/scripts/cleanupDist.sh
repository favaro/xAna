#!/bin/bash


source etc/scripts/setup_xAna.sh

cd $xAnaExtLIBS/VertexAnalysis
echo "cleanup VertexAnalysis"
make clean

cd -
cd $xAnaExtLIBS/JetMETObjects
echo "cleanup JetMETObjects"
make clean

cd -
cd $xAnaExtLIBS/GBRLikelihood
echo "cleanup GBRLikelihood"
make clean
rm -f Makefile
rm -f MTools

cd -
cd $xAnaExtLIBS/MCFM-6.6
make clean
#rm -rf Bin/
for lib in `find  . -name "*.a"`; do
    rm -rf $lib
done

cd -
cd $xAnaExtLIBS/JHUGenMELA
make clean

cd -
make clean
rm -rf $xAnaExtLIBS/lib/
rm -rf $xAnaExtLIBS/interface/
rm -f  $xAnaExtLIBS/VertexAnalysis
rm -f  $xAnaExtLIBS/JHUGenMELA
rm -f  $xAnaExtLIBS/MCFM-6.6
rm -f Makefile
rm -f dictgen.mk gen.mk
rm -f etc/scripts/env.sh