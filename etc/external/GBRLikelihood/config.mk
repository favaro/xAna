LIBRARYNAME := HiggsAnalysisGBRLikelihood

BOOST_LIBS := -lboost_program_options 
EXTLIBS += -lRooFitCore -lRooFit -lRooStats 
EXTINC  += -I$(VDTDIR)/include/
CXXSPECIFIC += -fopenmp 