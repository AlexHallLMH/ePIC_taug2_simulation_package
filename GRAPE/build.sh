#!/bin/csh -f

echo "Starting build process..."

# GRAPE BUILD
echo "Building grape-dilepton_adapted.1k..."
cd grape-dilepton_adapted.1k

rm -r build
mkdir build

cd basesv5.1
make clean
make

cd ../chanel
make clean
make

cd ../kinemlib
make clean
make

cd ../src
source set_grape_spring

make clean
cp ../../HepMC3-3.2.6/build/examples/Pythia6Example/CMakeFiles/pythia6_example.exe.dir/__/__/interfaces/pythia6/include/Pythia6/Pythia6ToHepMC3.cc.o .



make -f Makefile.spring integ
make -f Makefile.spring spring

cp integ ../build
cp spring ../build
cd ../build
cp ../cards/grape.cards .


echo "Build finished successfully!"