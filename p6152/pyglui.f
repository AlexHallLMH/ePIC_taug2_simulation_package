 
C*********************************************************************
 
C...PYGLUI
C...Calculates gluino decay modes.
 
      SUBROUTINE PYGLUI(KFIN,XLAM,IDLAM,IKNT)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KEXCIT=4000000)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4)
      COMMON/PYINTS/XXM(20)
      SAVE /PYDAT1/,/PYDAT2/,/PYMSSM/,/PYSSMT/,/PYINTS/
 
C...Local variables.
      INTEGER KFIN,KCIN,KF
      DOUBLE PRECISION XMI,XMJ,XMF,XMSF1,XMSF2,XMW,XMW2,
     &XMZ,XMZ2,AXMJ,AXMI
      DOUBLE PRECISION XMI2,XMI3,XMJ2,XMA2,XMB2,XMFP
      DOUBLE PRECISION C1L,C1R,D1L,D1R
      DOUBLE PRECISION C2L,C2R,D2L,D2R
      DOUBLE PRECISION PYLAMF,XL
      DOUBLE PRECISION TANW,XW,AEM,C1,AS,S12MAX,S12MIN
      DOUBLE PRECISION CA,CB,AL,AR,BL,BR
      DOUBLE PRECISION ALFA,BETA
      DOUBLE PRECISION SW,CW,SINB,COSB,QT,T3
      DOUBLE PRECISION XLAM(0:200)
      INTEGER IDLAM(200,3)
      INTEGER LKNT,IX,IC,ILR,IDU,J,IJ,I,IKNT,IFL
      DOUBLE PRECISION SR2
      DOUBLE PRECISION GAM
      DOUBLE PRECISION PYALEM,PI,PYALPS,EI
      EXTERNAL PYGAUS,PYXXZ5,PYXXW5,PYXXZ2
      DOUBLE PRECISION PYGAUS,PYXXZ5,PYXXW5,PYXXZ2
      DOUBLE PRECISION PREC
      INTEGER KFNCHI(4),KFCCHI(2)
      DATA PI/3.141592654D0/
      DATA SR2/1.4142136D0/
      DATA PREC/1D-2/
      DATA KFNCHI/1000022,1000023,1000025,1000035/
      DATA KFCCHI/1000024,1000037/
 
C...COUNT THE NUMBER OF DECAY MODES
      LKNT=0
      IF(KFIN.NE.KSUSY1+21) RETURN
      KCIN=PYCOMP(KFIN)
 
      XMW=PMAS(24,1)
      XMW2=XMW**2
      XMZ=PMAS(23,1)
      XMZ2=XMZ**2
      XW=PARU(102)
      TANW = SQRT(XW/(1D0-XW))
 
      XMI=PMAS(KCIN,1)
      AXMI=ABS(XMI)
      XMI2=XMI**2
      AEM=PYALEM(XMI2)
      AS =PYALPS(XMI2)
      C1=AEM/XW
      XMI3=XMI**3
      BETA=ATAN(RMSS(5))
 
C...2-BODY DECAYS OF GLUINO -> GRAVITINO GLUON
 
      IF(IMSS(11).EQ.1) THEN
        XMP=RMSS(29)
        IDG=39+KSUSY1
        XMGR=PMAS(PYCOMP(IDG),1)
        XFAC=(XMI2/(XMP*XMGR))**2*XMI/48D0/PI
        IF(AXMI.GT.XMGR) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=21
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC
        ENDIF
      ENDIF
 
C...2-BODY DECAYS OF GLUINO -> QUARK SQUARK
 
      DO 110 IFL=1,6
        DO 100 ILR=1,2
          XMJ=PMAS(PYCOMP(ILR*KSUSY1+IFL),1)
          AXMJ=ABS(XMJ)
          XMF=PMAS(IFL,1)
          IDU=3-(1+MOD(IFL,2))
          IF(XMI.GE.AXMJ+XMF) THEN
C...Minus sign difference from gluino-quark-squark feynman rules
            AL=SFMIX(IFL,1)
            BL=-SFMIX(IFL,3)
            AR=SFMIX(IFL,2)
            BR=-SFMIX(IFL,4)
C...F1 -> F CHI
            IF(ILR.EQ.1) THEN
              CA=AL
              CB=BL
C...F2 -> F CHI
            ELSE
              CA=AR
              CB=BR
            ENDIF
            LKNT=LKNT+1
            XMA2=XMJ**2
            XMB2=XMF**2
            XL=PYLAMF(XMI2,XMA2,XMB2)
            XLAM(LKNT)=4D0/8D0*AS/4D0/XMI3*SQRT(XL)*((XMI2+XMB2-XMA2)*
     &      (CA**2+CB**2)-4D0*CA*CB*XMI*XMF)
            IDLAM(LKNT,1)=ILR*KSUSY1+IFL
            IDLAM(LKNT,2)=-IFL
            IDLAM(LKNT,3)=0
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=0
          ENDIF
  100   CONTINUE
  110 CONTINUE
 
C...3-BODY DECAYS TO GAUGINO FERMION-FERMION
C...GLUINO -> NI Q QBAR
      DO 160 IX=1,4
        XMJ=SMZ(IX)
        AXMJ=ABS(XMJ)
        IF(XMI.GE.AXMJ) THEN
          XXM(1)=0D0
          XXM(2)=XMJ
          XXM(3)=0D0
          XXM(4)=XMI
          XXM(5)=PMAS(PYCOMP(KSUSY1+1),1)
          XXM(6)=PMAS(PYCOMP(KSUSY2+1),1)
          XXM(7)=1D6
          XXM(8)=0D0
          XXM(9)=0D0
          XXM(10)=0D0
          S12MIN=0D0
          S12MAX=(XMI-AXMJ)**2
C...D-TYPE QUARKS
          XXM(11)=0D0
          XXM(12)=0D0
          XXM(13)=1D0
          XXM(14)=-SR2*(-0.5D0*ZMIX(IX,2)+TANW*ZMIX(IX,1)/6D0)
          XXM(15)=1D0
          XXM(16)=SR2*(-TANW*ZMIX(IX,1)/3D0)
          IF( XXM(5).LT.AXMI .OR. XXM(6).LT.AXMI ) GOTO 120
          IF(XMI.GE.AXMJ+2D0*PMAS(1,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1*AS/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-2)
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=1
            IDLAM(LKNT,3)=-1
          ENDIF
          IF(XMI.GE.AXMJ+2D0*PMAS(3,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=3
            IDLAM(LKNT,3)=-3
          ENDIF
  120     CONTINUE
          IF( XXM(5).LT.AXMI .OR. XXM(6).LT.AXMI ) GOTO 130
          IF(XMI.GE.AXMJ+2D0*PMAS(5,1)) THEN
            CALL PYTBBN(IX,80,-1D0/3D0,AXMI,GAM)
            LKNT=LKNT+1
            XLAM(LKNT)=GAM
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=5
            IDLAM(LKNT,3)=-5
          ENDIF
C...U-TYPE QUARKS
  130     CONTINUE
          XXM(5)=PMAS(PYCOMP(KSUSY1+2),1)
          XXM(6)=PMAS(PYCOMP(KSUSY2+2),1)
          XXM(13)=1D0
          XXM(14)=-SR2*(0.5D0*ZMIX(IX,2)+TANW*ZMIX(IX,1)/6D0)
          XXM(15)=1D0
          XXM(16)=SR2*(2D0*TANW*ZMIX(IX,1)/3D0)
          IF( XXM(5).LT.AXMI .OR. XXM(6).LT.AXMI ) GOTO 140
          IF(XMI.GE.AXMJ+2D0*PMAS(2,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1*AS/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-2)
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=2
            IDLAM(LKNT,3)=-2
          ENDIF
          IF(XMI.GE.AXMJ+2D0*PMAS(4,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=4
            IDLAM(LKNT,3)=-4
          ENDIF
  140     CONTINUE
C...INCLUDE THE DECAY GLUINO -> NJ + T + T~
C...IF THE DECAY GLUINO -> ST + T CANNOT OCCUR
          IF(XMI.GE.PMAS(PYCOMP(KSUSY1+6),1)+PMAS(6,1)) GOTO 150
          XMF=PMAS(6,1)
          IF(XMI.GE.AXMJ+2D0*XMF) THEN
            CALL PYTBBN(IX,80,2D0/3D0,AXMI,GAM)
            LKNT=LKNT+1
            XLAM(LKNT)=GAM
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=6
            IDLAM(LKNT,3)=-6
          ENDIF
  150     CONTINUE
        ENDIF
  160 CONTINUE
 
C...GLUINO -> CI Q QBAR'
      DO 190 IX=1,2
        XMJ=SMW(IX)
        AXMJ=ABS(XMJ)
        IF(XMI.GE.AXMJ) THEN
          S12MIN=0D0
          S12MAX=(AXMI-AXMJ)**2
          XXM(1)=0D0
          XXM(2)=XMJ
          XXM(3)=0D0
          XXM(4)=XMI
          XXM(5)=0D0
          XXM(6)=0D0
          XXM(9)=1D6
          XXM(10)=0D0
          XXM(7)=UMIX(IX,1)*SR2
          XXM(8)=VMIX(IX,1)*SR2
          XXM(11)=PMAS(PYCOMP(KSUSY1+1),1)
          XXM(12)=PMAS(PYCOMP(KSUSY1+2),1)
          IF( XXM(11).LT.AXMI .OR. XXM(12).LT.AXMI ) GOTO 170
          IF(XMI.GE.AXMJ+PMAS(1,1)+PMAS(2,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=0.5D0*C1*AS/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXW5,S12MIN,S12MAX,PREC)
            IDLAM(LKNT,1)=KFCCHI(IX)
            IDLAM(LKNT,2)=1
            IDLAM(LKNT,3)=-2
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
          ENDIF
          IF(XMI.GE.AXMJ+PMAS(3,1)+PMAS(4,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFCCHI(IX)
            IDLAM(LKNT,2)=3
            IDLAM(LKNT,3)=-4
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
          ENDIF
  170     CONTINUE
 
          IF(XMI.GE.PMAS(PYCOMP(KSUSY1+5),1)+PMAS(5,1)) GOTO 180
          IF(XMI.GE.PMAS(PYCOMP(KSUSY1+6),1)+PMAS(6,1)) GOTO 180
          XMF=PMAS(6,1)
          XMFP=PMAS(5,1)
          IF(XMI.GE.AXMJ+XMF+XMFP) THEN
            CALL PYTBBC(IX,80,AXMI,GAM)
            LKNT=LKNT+1
            XLAM(LKNT)=GAM
            IDLAM(LKNT,1)=KFCCHI(IX)
            IDLAM(LKNT,2)=5
            IDLAM(LKNT,3)=-6
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
          ENDIF
  180     CONTINUE
        ENDIF
  190 CONTINUE
 
      IKNT=LKNT
      XLAM(0)=0D0
      DO 200 I=1,IKNT
        IF(XLAM(I).LT.0D0) XLAM(I)=0D0
        XLAM(0)=XLAM(0)+XLAM(I)
  200 CONTINUE
      IF(XLAM(0).EQ.0D0) XLAM(0)=1D-6
 
      RETURN
      END
