# Source files to compile
FILES := setupReader configReader category higgsCrossSections weightManager massResolutionCalculator
FILES += photonEnergyScale photonEnergySmearing
FILES += triggerSelection vertexSelection
FILES += xAnaInitializer xAnaMVAsInit xAnaCuts xAnaCategorize xAnaMCutils xAnaMainLoop
FILES += photonSelection leptonSelection jetSelection metSelection
FILES += diphotonSelection

# Header files to use for dictionary generation
# In this module: is the same as FILES
DICTFILES :=

# Executable files
PROGRAMS := xAnaExe
PROGRAMS := addMELA

NEEDS_ROOT  := yes
NEEDS_BOOST := yes
