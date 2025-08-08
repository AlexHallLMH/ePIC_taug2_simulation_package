 
C*********************************************************************
 
C...PYSAVE
C...Saves and restores parameter and cross section values for the
C...3 gamma-p and 6 (or 4, or 9, or 13) gamma-gamma alnternatives.
C...Also makes random choice between alternatives.
 
      SUBROUTINE PYSAVE(ISAVE,IGA)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      SAVE /PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT5/,/PYINT7/
C...Local arrays and saved variables.
      DIMENSION NCP(15),NSUBCP(15,20),MSUBCP(15,20),COEFCP(15,20,20),
     &NGENCP(15,0:20,3),XSECCP(15,0:20,3),SIGTCP(15,0:6,0:6,0:5),
     &INTCP(15,20),RECP(15,20)
      SAVE NCP,NSUBCP,MSUBCP,COEFCP,NGENCP,XSECCP,SIGTCP,INTCP,RECP
 
C...Save list of subprocesses and cross-section information.
      IF(ISAVE.EQ.1) THEN
        ICP=0
        DO 120 I=1,500
          IF(MSUB(I).EQ.0.AND.I.NE.96.AND.I.NE.97) GOTO 120
          ICP=ICP+1
          NSUBCP(IGA,ICP)=I
          MSUBCP(IGA,ICP)=MSUB(I)
          DO 100 J=1,20
            COEFCP(IGA,ICP,J)=COEF(I,J)
  100     CONTINUE
          DO 110 J=1,3
            NGENCP(IGA,ICP,J)=NGEN(I,J)
            XSECCP(IGA,ICP,J)=XSEC(I,J)
  110     CONTINUE
  120   CONTINUE
        NCP(IGA)=ICP
        DO 130 J=1,3
          NGENCP(IGA,0,J)=NGEN(0,J)
          XSECCP(IGA,0,J)=XSEC(0,J)
  130   CONTINUE
        DO 136 I1=0,6
          DO 134 I2=0,6
            DO 132 J=0,5
              SIGTCP(IGA,I1,I2,J)=SIGT(I1,I2,J)
  132       CONTINUE
  134     CONTINUE
  136   CONTINUE

C...Save various common process variables.
        DO 140 J=1,10
          INTCP(IGA,J)=MINT(40+J)
  140   CONTINUE
        INTCP(IGA,11)=MINT(101)
        INTCP(IGA,12)=MINT(102)
        INTCP(IGA,13)=MINT(107)
        INTCP(IGA,14)=MINT(108)
        INTCP(IGA,15)=MINT(123)
        RECP(IGA,1)=CKIN(3)
        RECP(IGA,2)=VINT(318)
 
C...Save cross-section information only.
      ELSEIF(ISAVE.EQ.2) THEN
        DO 160 ICP=1,NCP(IGA)
          I=NSUBCP(IGA,ICP)
          DO 150 J=1,3
            NGENCP(IGA,ICP,J)=NGEN(I,J)
            XSECCP(IGA,ICP,J)=XSEC(I,J)
  150     CONTINUE
  160   CONTINUE
        DO 170 J=1,3
          NGENCP(IGA,0,J)=NGEN(0,J)
          XSECCP(IGA,0,J)=XSEC(0,J)
  170   CONTINUE
 
C...Choose between allowed alternatives.
      ELSEIF(ISAVE.EQ.3.OR.ISAVE.EQ.4) THEN
        IF(ISAVE.EQ.4) THEN
          XSUMCP=0D0
          DO 180 IG=1,MINT(121)
            XSUMCP=XSUMCP+XSECCP(IG,0,1)
  180     CONTINUE
          XSUMCP=XSUMCP*PYR(0)
          DO 190 IG=1,MINT(121)
            IGA=IG
            XSUMCP=XSUMCP-XSECCP(IG,0,1)
            IF(XSUMCP.LE.0D0) GOTO 200
  190     CONTINUE
  200     CONTINUE
        ENDIF
 
C...Restore cross-section information.
        DO 210 I=1,500
          MSUB(I)=0
  210   CONTINUE
        DO 240 ICP=1,NCP(IGA)
          I=NSUBCP(IGA,ICP)
          MSUB(I)=MSUBCP(IGA,ICP)
          DO 220 J=1,20
            COEF(I,J)=COEFCP(IGA,ICP,J)
  220     CONTINUE
          DO 230 J=1,3
            NGEN(I,J)=NGENCP(IGA,ICP,J)
            XSEC(I,J)=XSECCP(IGA,ICP,J)
  230     CONTINUE
  240   CONTINUE
        DO 250 J=1,3
          NGEN(0,J)=NGENCP(IGA,0,J)
          XSEC(0,J)=XSECCP(IGA,0,J)
  250   CONTINUE
        DO 256 I1=0,6
          DO 254 I2=0,6
            DO 252 J=0,5
              SIGT(I1,I2,J)=SIGTCP(IGA,I1,I2,J)
  252       CONTINUE
  254     CONTINUE
  256   CONTINUE
 
C...Restore various common process variables.
        DO 260 J=1,10
          MINT(40+J)=INTCP(IGA,J)
  260   CONTINUE
        MINT(101)=INTCP(IGA,11)
        MINT(102)=INTCP(IGA,12)
        MINT(107)=INTCP(IGA,13)
        MINT(108)=INTCP(IGA,14)
        MINT(123)=INTCP(IGA,15)
        CKIN(3)=RECP(IGA,1)
        CKIN(1)=2D0*CKIN(3)
        VINT(318)=RECP(IGA,2)
 
C...Sum up cross-section info (for PYSTAT).
      ELSEIF(ISAVE.EQ.5) THEN
        DO 270 I=1,500
          MSUB(I)=0
          NGEN(I,1)=0
          NGEN(I,3)=0
          XSEC(I,3)=0D0
  270   CONTINUE
        NGEN(0,1)=0
        NGEN(0,2)=0
        NGEN(0,3)=0
        XSEC(0,3)=0
        DO 290 IG=1,MINT(121)
          DO 280 ICP=1,NCP(IG)
            I=NSUBCP(IG,ICP)
            IF(MSUBCP(IG,ICP).EQ.1) MSUB(I)=1
            NGEN(I,1)=NGEN(I,1)+NGENCP(IG,ICP,1)
            NGEN(I,3)=NGEN(I,3)+NGENCP(IG,ICP,3)
            XSEC(I,3)=XSEC(I,3)+XSECCP(IG,ICP,3)
  280     CONTINUE
          NGEN(0,1)=NGEN(0,1)+NGENCP(IG,0,1)
          NGEN(0,2)=NGEN(0,2)+NGENCP(IG,0,2)
          NGEN(0,3)=NGEN(0,3)+NGENCP(IG,0,3)
          XSEC(0,3)=XSEC(0,3)+XSECCP(IG,0,3)
  290   CONTINUE
      ENDIF
 
      RETURN
      END
