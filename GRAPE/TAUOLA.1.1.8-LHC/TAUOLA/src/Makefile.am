includedir=$(prefix)/include/Tauola
lib_LTLIBRARIES = libTauolaCxxInterface.la libTauolaHEPEVT.la

include_HEADERS =                            \
eventRecordInterfaces/TauolaHEPEVTParticle.h \
eventRecordInterfaces/TauolaHEPEVTEvent.h    \
tauolaCInterfaces/Version.h                  \
tauolaCInterfaces/Tauola.h                   \
tauolaCInterfaces/TauolaParticle.h           \
tauolaCInterfaces/TauolaParticlePair.h       \
tauolaCInterfaces/TauolaEvent.h              \
tauolaCInterfaces/DecayList.h                \
tauolaFortranInterfaces/f_FilHep.h           \
tauolaFortranInterfaces/f_Init.h             \
tauolaFortranInterfaces/f_Decay.h            \
tauolaFortranInterfaces/f_Variables.h        \
utilities/Log.h                              \
utilities/Plots.h

libTauolaCxxInterface_la_SOURCES =            \
tauolaCInterfaces/Tauola.cxx                  \
tauolaCInterfaces/TauolaParticle.cxx          \
tauolaCInterfaces/TauolaParticlePair.cxx      \
tauolaCInterfaces/TauolaEvent.cxx             \
tauolaCInterfaces/DecayList.cxx               \
tauolaFortranInterfaces/f_FilHep.cxx          \
tauolaFortranInterfaces/f_Init.cxx            \
tauolaFortranInterfaces/f_Decay.cxx           \
tauolaFortranInterfaces/tauola_extras.f       \
utilities/Log.cxx                             \
utilities/Plots.cxx

libTauolaCxxInterface_la_LIBADD = ../tauola-fortran/libTauolaFortran.la

libTauolaHEPEVT_la_SOURCES =                   \
eventRecordInterfaces/TauolaHEPEVTEvent.cxx    \
eventRecordInterfaces/TauolaHEPEVTParticle.cxx

if HAS_HEPMC
  lib_LTLIBRARIES += libTauolaHepMC.la

  libTauolaHepMC_la_SOURCES = \
  eventRecordInterfaces/TauolaHepMCEvent.cxx    \
  eventRecordInterfaces/TauolaHepMCParticle.cxx

  include_HEADERS += \
  eventRecordInterfaces/TauolaHepMCParticle.h  \
  eventRecordInterfaces/TauolaHepMCEvent.h
endif

if HAS_HEPMC3
  lib_LTLIBRARIES += libTauolaHepMC3.la

  libTauolaHepMC3_la_SOURCES = \
  eventRecordInterfaces/TauolaHepMC3Event.cxx    \
  eventRecordInterfaces/TauolaHepMC3Particle.cxx

  include_HEADERS += \
  eventRecordInterfaces/TauolaHepMC3Particle.h  \
  eventRecordInterfaces/TauolaHepMC3Event.h
endif


AM_CPPFLAGS = -I$(top_srcdir)/src/tauolaFortranInterfaces \
              -I$(top_srcdir)/src/tauolaCInterfaces \
              -I$(top_srcdir)/src/utilities \
              -I$(HEPMC_DIR)/include

libTauolaCxxInterface_la_FFLAGS = -O2 -fPIC -fno-automatic -fno-backslash -ffixed-line-length-132
