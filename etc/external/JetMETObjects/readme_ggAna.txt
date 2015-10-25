## Fabrice 22/01/2014
replace the JetMETObjects with the one from CMSSW_6_0_0 adding all the source so it can be compiled properly.
I had to change several things in the actual code, namely:
- cms::Exception --> hanleError
- remove some interface
- in Makefile:  standalonedir ../ --> ./

To compile this:
cd JetMETObjects
make




