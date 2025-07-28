include ../../../make.inc

OBJ=

HEPMC_LIBS     = -L$(HEPMCLOCATION)/lib -lHepMC
TAUOLA_LIBS    = -L$(PREFIX)/lib -L../lib -lTauolaTauSpinner -lTauolaFortran -lTauolaCxxInterface
LHAPDF_LIBS    = -L$(LHAPDFLOCATION)/lib -lLHAPDF
ROOT_LIBS      = $(shell root-config --glibs)
ROOT_INC       = $(shell root-config --cflags)

all:
	@make table-parsing-test.exe


%.o: %.cxx
	$(CC) $(CFLAGS) -c $< -o $@ -I. -I$(PREFIX)/include $(ROOT_INC)

%.o: %.f
	$(F77) $(F77FLAGS) -c $< -o $@

%.exe: %.o $(OBJ)
	$(CC) -o $@ $^ $(TAUOLA_LIBS) $(HEPMC_LIBS) $(LHAPDF_LIBS) $(ROOT_LIBS) \
	-Wl,-rpath,$(PREFIX)/lib \
	-Wl,-rpath,$(HEPMCLOCATION)/lib \
	-Wl,-rpath,$(LHAPDFLOCATION)/lib \
	-lgfortran

clean:
	rm -f *.exe *.o

