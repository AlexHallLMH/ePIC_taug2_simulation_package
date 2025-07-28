#!/bin/csh -f

echo "Starting build process..."

# GRAPE BUILD
echo "Building grape-dilepton_adapted.1k..."
cd grape-dilepton_adapted.1k

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
cp ../Pythia6toHepMC3.cc.o .

make -f Makefile.spring integ
make -f Makefile.spring spring

cp integ ../build
cp spring ../build
cp ../cards/grape.cards .


echo "Build finished successfully!"