#!/bin/bash

# INSTRUCTIONS:
# Run this using ./install.sh from the top directory
# The first time running, you may need to run chmod +x install.sh

# Make pythia
cd pythia8315
gmake clean
./configure
make

cd ../upcgen-master

# Set Pythia8 path
export PYTHIA8=../pythia8315

# Create build directory
rm -r build
mkdir -p build
cd build

# Configure with CMake (disable Pythia6)
cmake -DROOT_ROOT_DIR=$(/bin/root-config --prefix) -DBUILD_WITH_PYTHIA6=OFF ..

# Fix /lib/libEG.so path to just /libEG.so
sed -i 's|/lib/libEG.so|/usr/lib64/root/libEG.so|g' CMakeFiles/upcgen.dir/build.make

# Compile
make -j$(nproc)

cp ../parameters.in .
cp ../run_analysis.py .

echo "Build complete. Run in build directory with: ./upcgen -parfile parameters.in"
