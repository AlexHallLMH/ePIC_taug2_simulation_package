#!/bin/csh -f

cd grape-dilepton_adapted.1k/build

cp ../tauola_interface/clean.py .
cp ../tauola_interface/reader_tester.py .

#Duplicate the original output, to leave the initial output untouched
cp Grape.hepmc Grape_dup.hepmc

#Run Cleaning script (fixes corrupted file issues and )
python3 clean.py

#Remove redundant files
rm clean.py
rm reader_tester.py
rm Grape_dup.hepmc
rm Grape_temp.hepmc

#Copy, compile and run the Tauola processing script
cp ../tauola_interface/tauola_decay.cpp .
g++ tauola_decay.cpp -o tauola_decay -I../../HepMC3-3.2.6/include -I../../TAUOLA1.1.8-LHC/Tauola/include -L../../HepMC3-3.2.6/build/outputs/lib64 -Wl,-rpath,../../HepMC3-3.2.6/build/outputs/lib64 -lHepMC3 -L../../TAUOLA.1.1.8-LHC/TAUOLA/lib -lTauolaCxxInterface -lTauolaFortran -lTauolaHEPEVT -lTauolaHepMC3 -std=c++11
rm tauola_decay.cpp

setenv LD_LIBRARY_PATH ../../TAUOLA.1.1.8-LHC/TAUOLA/lib:$LD_LIBRARY_PATH
./tauola_decay

#Tidies final version ready for TauSpinner
cp ../tauola_interface/tidy.py .
python3 tidy.py
rm tidy.py