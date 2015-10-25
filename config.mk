LIBRARYNAME := xAna

BOOST_LIBS := -lboost_program_options
EXTLIBS += -Letc/external/lib/ -Letc/external/JHUGenerator.v4.0.1/MCFM-6.6/obj/ -lh2gglobeVertexAnalysis -lHiggsAnalysisGBRLikelihood -lJetMETObjects  -lME -lgfortran -lPhysics -lEG -lm -lGenVector -lmcfm -lpv -lqcdloop -lov -lff -lsmallG -lsmallY -lsmallP -lsmallF
#-lRooFitCore -lRooFit 
EXTINC  += -Ietc/external/interface/
CXXSPECIFIC += -fopenmp
