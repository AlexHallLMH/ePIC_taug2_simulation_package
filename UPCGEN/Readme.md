# UPCGEN and Pythia8 installation Guide

1. Run: chmod +x install.sh
2. Run ./install.sh
3. The build will be found in upcgen-master/build. The parameters in paramters.in are set up to be similar to ePIC. To run a simulation, first run ./upcgen (this will take quite a while). Once this completes, the events can be found in events.root or events.hepmc
4. To create the standard histograms, run python3 analysis.py from the build directory. The histograms will be saved to upcgen-master/build/outputs

