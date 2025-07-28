export F77=gfortran
export FFLAGS='-O -g -fno-automatic -ffixed-line-length-132 -std=legacy '
export USRLDIR="-L/network/software/el9/cernlib/2006a/new/lib"
export SYSLIB="-L/home/lady6801/Documents/new/p6152 -L/home/lady6801/Documents/new/TAUOLA.1.1.8-LHC/TAUOLA/lib/ -L/network/software/el9/cernlib/2006a/new/lib"
export USRLIB="${SYSLIB} -lkernlib -lpdflib804 -lmathlib -lpacklib -lpythia6152 -lkernlib "

export GRACELDIR="../lib"
source ./set_grape.usrobj.ksh
