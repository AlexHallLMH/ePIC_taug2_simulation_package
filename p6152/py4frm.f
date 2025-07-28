 
C*********************************************************************
 
C...PY4FRM
C...An interface from a four-fermion generator to include
C...parton showers and hadronization.
 
      SUBROUTINE PY4FRM(ATOTSQ,A1SQ,A2SQ,ISTRAT,IRAD,ITAU,ICOM)

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYJETS/,/PYDAT1/
C...Local arrays.
      DIMENSION IJOIN(2),INTAU(4)
 
C...Call PYHEPC to convert input from HEPEVT to PYJETS common.
      IF(ICOM.EQ.0) THEN 
        MSTU(28)=0
        CALL PYHEPC(2)
      ENDIF
 
C...Loop through entries and pick up all final fermions/antifermions.
      I1=0
      I2=0
      I3=0
      I4=0
      DO 100 I=1,N
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 100
      KFA=IABS(K(I,2))
      IF((KFA.GE.1.AND.KFA.LE.6).OR.(KFA.GE.11.AND.KFA.LE.16)) THEN
        IF(K(I,2).GT.0) THEN
          IF(I1.EQ.0) THEN
            I1=I
          ELSEIF(I3.EQ.0) THEN
            I3=I
          ELSE
            CALL PYERRM(16,'(PY4FRM:) more than two fermions')
          ENDIF
        ELSE
          IF(I2.EQ.0) THEN
            I2=I
          ELSEIF(I4.EQ.0) THEN
            I4=I
          ELSE
            CALL PYERRM(16,'(PY4FRM:) more than two antifermions')
          ENDIF
        ENDIF
      ENDIF
  100 CONTINUE
 
C...Check that event is arranged according to conventions.
      IF(I3.EQ.0.OR.I4.EQ.0) THEN
        CALL PYERRM(16,'(PY4FRM:) event contains too few fermions')
      ENDIF
      IF(I2.LT.I1.OR.I3.LT.I2.OR.I4.LT.I3) THEN
        CALL PYERRM(6,'(PY4FRM:) fermions arranged in wrong order')
      ENDIF
 
C...Check which fermion pairs are quarks and which leptons.
      IF(IABS(K(I1,2)).LT.10.AND.IABS(K(I2,2)).LT.10) THEN
        IQL12=1
      ELSEIF(IABS(K(I1,2)).GT.10.AND.IABS(K(I2,2)).GT.10) THEN
        IQL12=2
      ELSE
        CALL PYERRM(16,'(PY4FRM:) first fermion pair inconsistent')
      ENDIF
      IF(IABS(K(I3,2)).LT.10.AND.IABS(K(I4,2)).LT.10) THEN
        IQL34=1
      ELSEIF(IABS(K(I3,2)).GT.10.AND.IABS(K(I4,2)).GT.10) THEN
        IQL34=2
      ELSE
        CALL PYERRM(16,'(PY4FRM:) second fermion pair inconsistent')
      ENDIF
 
C...Decide whether to allow or not photon radiation in showers.
      MSTJ(41)=2
      IF(IRAD.EQ.0) MSTJ(41)=1
 
C...Decide on dipole pairing.
      IP1=I1
      IP2=I2
      IP3=I3
      IP4=I4
      IF(IQL12.EQ.IQL34) THEN
        R1SQ=A1SQ
        R2SQ=A2SQ
        DELTA=ATOTSQ-A1SQ-A2SQ
        IF(ISTRAT.EQ.1) THEN
          IF(DELTA.GT.0D0) R1SQ=R1SQ+DELTA
          IF(DELTA.LT.0D0) R2SQ=MAX(0D0,R2SQ+DELTA)
        ELSEIF(ISTRAT.EQ.2) THEN
          IF(DELTA.GT.0D0) R2SQ=R2SQ+DELTA
          IF(DELTA.LT.0D0) R1SQ=MAX(0D0,R1SQ+DELTA)
        ENDIF
        IF(R2SQ.GT.PYR(0)*(R1SQ+R2SQ)) THEN
          IP2=I4
          IP4=I2
        ENDIF
      ENDIF
 
C...Do colour joinings and parton showers.
      IF(IQL12.EQ.1) THEN 
        IJOIN(1)=IP1
        IJOIN(2)=IP2
        CALL PYJOIN(2,IJOIN)
      ENDIF
      IF(IQL12.EQ.1.OR.IRAD.EQ.1) THEN
        PM12S=(P(IP1,4)+P(IP2,4))**2-(P(IP1,1)+P(IP2,1))**2-
     &  (P(IP1,2)+P(IP2,2))**2-(P(IP1,3)+P(IP2,3))**2
        CALL PYSHOW(IP1,IP2,SQRT(MAX(0D0,PM12S)))
      ENDIF
      IF(IQL34.EQ.1) THEN
        IJOIN(1)=IP3
        IJOIN(2)=IP4
        CALL PYJOIN(2,IJOIN)
      ENDIF
      IF(IQL34.EQ.1.OR.IRAD.EQ.1) THEN        
        PM34S=(P(IP3,4)+P(IP4,4))**2-(P(IP3,1)+P(IP4,1))**2-
     &  (P(IP3,2)+P(IP4,2))**2-(P(IP3,3)+P(IP4,3))**2
        CALL PYSHOW(IP3,IP4,SQRT(MAX(0D0,PM34S)))
      ENDIF
 
C...Do fragmentation and decays. Possibly except tau decay.
      IF(ITAU.EQ.0) THEN
        NTAU=0
        DO 110 I=1,N
        IF(IABS(K(I,2)).EQ.15.AND.K(I,1).EQ.1) THEN
          NTAU=NTAU+1
          INTAU(NTAU)=I
          K(I,1)=11
        ENDIF
  110   CONTINUE
      ENDIF       
      CALL PYEXEC
      IF(ITAU.EQ.0) THEN
        DO 120 I=1,NTAU
        K(INTAU(I),1)=1
  120   CONTINUE
      ENDIF       
 
C...Call PYHEPC to convert output from PYJETS to HEPEVT common.
      IF(ICOM.EQ.0) THEN 
        MSTU(28)=0
        CALL PYHEPC(1)
      ENDIF
 
      END
