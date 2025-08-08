#!/bin/bash

rm -r build
rm -r HepMC3-3.2.6/build

cd HepMC3-3.2.6
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX=../hepmc3-install   \
        -DHEPMC3_ENABLE_ROOTIO:BOOL=OFF            \
        -DHEPMC3_ENABLE_PROTOBUFIO:BOOL=OFF        \
        -DHEPMC3_ENABLE_TEST:BOOL=OFF              \
        -DHEPMC3_INSTALL_INTERFACES:BOOL=ON        \
        -DHEPMC3_BUILD_STATIC_LIBS:BOOL=OFF        \
        -DHEPMC3_BUILD_DOCS:BOOL=OFF     \
        -DHEPMC3_ENABLE_PYTHON:BOOL=ON   \
        -DHEPMC3_PYTHON_VERSIONS=3.12     \
        -DHEPMC3_Python_SITEARCH312=../hepmc3-install/lib/python3.12/site-packages \
        ../ -DHEPMC3_BUILD_EXAMPLES=ON

make

cd ../..


#Compile Tauola
ls
cd TAUOLA.1.1.8-LHC/TAUOLA
make clean
./configure --with-hepmc3=../../HepMC3-3.2.6