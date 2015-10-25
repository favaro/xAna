# Source files to compile
FILES := GBRArrayUtils GBRMath 
FILES += HybridGBREvent HybridGBRForest HybridGBRForestD HybridGBRTree HybridGBRTreeD
FILES += RooCBExp RooCBFast RooDoubleCBFast RooGaussianFast RooHybridBDTAutoPdf RooRevCBFast

# Header files to use for dictionary generation
# In this module: is the same as FILES
DICTFILES := 

# Executable files
PROGRAMS := 

NEEDS_ROOT  := yes
NEEDS_BOOST := yes
