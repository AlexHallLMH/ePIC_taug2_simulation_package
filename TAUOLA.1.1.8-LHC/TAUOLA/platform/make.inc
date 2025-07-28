exec_prefix   = /tmp/TAUOLA.1.1.8/TAUOLA

PREFIX        = /tmp/TAUOLA.1.1.8/TAUOLA
LIBDIR        = ${exec_prefix}/lib
HEPMCLOCATION = 
HEPMC3LOCATION = 
LHAPDFLOCATION   = 
PYTHIALOCATION   = 
MAKE_TAU_SPINNER = 

LIB_VER = 1.1.8

F77      = gfortran
F77FLAGS = -fPIC -fno-automatic -fno-backslash -ffixed-line-length-132
FFLAGS   = -O2 $(F77FLAGS)

LD       = gfortran
LDFLAGS  = 
SOFLAGS  = -shared

RANLIB   = ranlib
AR       = ar
DIFF    = diff -ib

CPP     = cpp
CPPFLAGS= -C -P -traditional-cpp

CC      = g++
CFLAGS  = -O2 -fPIC -Wall -Wno-write-strings -std=c++11

GCC     = gcc
GCCFLAGS= -O3 -fomit-frame-pointer -DHAVE_UNDERSCORE

#Uncomment to compile with debug mode turned on
#DEBUG = -D"_LOG_DEBUG_MODE_"
