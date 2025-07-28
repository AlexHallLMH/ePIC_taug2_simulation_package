#!/bin/csh -f

#Move into build filder
cd grape-dilepton_adapted.1k/build
setenv LD_LIBRARY_PATH ../../HepMC3-3.2.6/build/outputs/lib64

#Perform integration step
./integ

#Generate events
./spring

#Remove Unneeded outputs
rm bases.rz
rm grp.rz