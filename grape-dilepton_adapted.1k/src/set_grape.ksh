export F77=gfortran
export FFLAGS='-O -g -fno-automatic -ffixed-line-length-132 -std=legacy '
export USRLDIR="-L/network/software/el9/cernlib/2006a/new/lib"
export SYSLIB="-L/home/lady6801/Documents/new/p6152"
export USRLIB="${SYSLIB} -lkernlib -lpdflib804 -lmathlib -lpacklib -lpythia6152 -lkernlib -L/home/lady6801/Documents/new/TAUOLA.1.1.8-LHC/TAUOLA/lib/libTauolaCxxInterface.a -L/home/lady6801/Documents/new/TAUOLA.1.1.8-LHC/TAUOLA/lib/libTauolaFortran.a -L/home/lady6801/Documents/new/TAUOLA.1.1.8-LHC/TAUOLA/lib/libTauolaTauSpinner.a -L/home/lady6801/Documents/new/TAUOLA.1.1.8-LHC/TAUOLA/lib/libTauolaHEPEVT.a"

export GRACELDIR="../lib"
source ./set_grape.usrobj.ksh
