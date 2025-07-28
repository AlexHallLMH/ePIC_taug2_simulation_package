C...A simple skeleton program, illustrating a typical Pythia run:
C...Z0 production at LEP 1. 
C...Toy task: compare multiplicity distribution with matrix elements
C...and with parton showers (using same fragmentation parameters).
C...This code contains modifications for HepMC3 examples
C-----------------------------------------------------------------

C...Preamble: declarations.
 
C...All real arithmetic in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP

C...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

C...Commonblocks.
C...The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
C...Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
C...Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
C...Parameters. 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
C...Supersymmetry parameters.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
C...Generation and cross section statistics.
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
C...Random number generator information.
      COMMON/PYDATR/MRPY(6),RRPY(100)

C...HepMC3
      PARAMETER (NMXHEP=4000)
      COMMON /HEPEVT/  NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &                 JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
     &                 VHEP(4,NMXHEP)
      INTEGER          NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      DOUBLE PRECISION PHEP,VHEP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
include "Pythia6ToHepMC3.inc"
      INTEGER OUTID(2), HEPMC3STATUS

C-----------------------------------------------------------------

C...First section: initialization.
 

C...Create output writers
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      OUTID(1)=HepMC3_new_writer(0,1,'Grape.hepmc'//char(0))
c      HEPMC3STATUS=HepMC3_new_weight(OUTID(1),'Default'//char(0))      
c      HEPMC3STATUS=HepMC3_new_weight(OUTID(1),'weme1'//char(0))      
c      HEPMC3STATUS=HepMC3_new_weight(OUTID(1),'weme2'//char(0))

      NEVHEP=-123456
      HEPMC3STATUS=HepMC3_set_hepevt_address(NEVHEPL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C-----------------------------------------------------------------

C...Second section: event loop.

      ica = 1
 
C...Begin event loop.
      DO 200 IEV=1,4000
         read(77,*,END=99) NEVHEP,nhep
         NEVHEPL=IEV
         NHEPL=NHEP
         do 500 j = 1,nhep
	    read(77,10,END=99) isthep(j),idhep(j),jmohep(1,j),jmohep(2,j),
     +           jdahep(1,j),jdahep(2,j)
 10         FORMAT(6(1x,i6))
	    read(77,20) (phep(jz,j),jz = 1,4),(vhep(jz,j),jz = 1,4)
 20         FORMAT(4(1x,f10.5),4(1x,f10.5))
            ISTHEPL(J)=ISTHEP(J)
            IDHEPL(J)=IDHEP(J)
            JMOHEPL(1,J)=JMOHEP(1,J)
            JMOHEPL(2,J)=JMOHEP(2,J)
            JDAHEPL(1,J)=JDAHEP(1,J)
            JDAHEPL(2,J)=JDAHEP(2,J)
            PHEPL(1,J)=PHEP(1,J)
            PHEPL(2,J)=PHEP(2,J)
            PHEPL(3,J)=PHEP(3,J)
            PHEPL(4,J)=PHEP(4,J)
            PHEPL(5,J)=PHEP(5,J)
            VHEPL(1,J)=VHEP(1,J)
            VHEPL(2,J)=VHEP(2,J)
            VHEPL(3,J)=VHEP(3,J)
            VHEPL(4,J)=VHEP(4,J)
 500     CONTINUE                    
         

         HEPMC3STATUS=HepMC3_convert_event(OUTID(ICA))


C     Note there should be PDF ids
         HEPMC3STATUS=HepMC3_write_event(OUTID(ICA))
         HEPMC3STATUS=HepMC3_clear_event(OUTID(ICA))          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C...  End event loop.
 200  CONTINUE
 99   continue
C...  End outer loop.
 300  CONTINUE
C...  Delete output writers
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      HEPMC3STATUS=HepMC3_delete_writer(OUTID(1))
c     HEPMC3STATUS=HepMC3_delete_writer(OUTID(2))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C-----------------------------------------------------------------

C...Third section: produce output and end.

C...Cross section table.


      END
