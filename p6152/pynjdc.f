 
C*********************************************************************
 
C...PYNJDC
C...Calculates decay widths for the neutralinos (admixtures of
C...Bino, W3-ino, Higgs1-ino, Higgs2-ino)
 
C...Input:  KCIN = KF code for particle
C...Output: XLAM = widths
C...        IDLAM = KF codes for decay particles
C...        IKNT = number of decay channels defined
C...AUTHOR: STEPHEN MRENNA
C...Last change:
C...10-15-95:  force decay chi^0_2 -> chi^0_1 + gamma
C...when CHIGAMMA .NE. 0
C...10 FEB 96:  Calculate this decay for small tan(beta)
 
      SUBROUTINE PYNJDC(KFIN,XLAM,IDLAM,IKNT)
 
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
      INTEGER KFIN,KCIN
      DOUBLE PRECISION XMI,XMJ,XMF,XMSF1,XMSF2,XMW,XMW2,
     &XMZ,XMZ2,AXMJ,AXMI
      DOUBLE PRECISION XMFP,XMF1,XMF2,XMSL,XMG,XMK
      DOUBLE PRECISION S12MIN,S12MAX
      DOUBLE PRECISION XMI2,XMI3,XMJ2,XMH,XMH2,XMHP,XMHP2,XMA2,XMB2
      DOUBLE PRECISION PYLAMF,XL,QIJ,RIJ
      DOUBLE PRECISION TANW,XW,AEM,C1,AS,EI,T3
      DOUBLE PRECISION PYX2XH,PYX2XG
      DOUBLE PRECISION XLAM(0:200)
      INTEGER IDLAM(200,3)
      INTEGER LKNT,IX,IH,J,IJ,I,IKNT,FID
      INTEGER ITH(3),KF1,KF2
      INTEGER ITHC
      DOUBLE PRECISION ETAH(3),CH(3),DH(3),EH(3)
      DOUBLE PRECISION SR2
      DOUBLE PRECISION CBETA,SBETA,GR,GL,F12K,F21K
      DOUBLE PRECISION GAMCON,XMT1,XMT2
      DOUBLE PRECISION PYALEM,PI,PYALPS
      DOUBLE PRECISION AL,BL,AR,BR,ALP,ARP,BLP,BRP
      DOUBLE PRECISION RAT1,RAT2
      DOUBLE PRECISION T3T,CA,CB,FCOL
      DOUBLE PRECISION ALFA,BETA,TANB
      DOUBLE PRECISION PYXXGA
      EXTERNAL PYXXW5,PYGAUS,PYXXZ5
      DOUBLE PRECISION PYXXW5,PYGAUS,PYXXZ5
      DOUBLE PRECISION PREC
      INTEGER KFNCHI(4),KFCCHI(2)
      DATA ETAH/1D0,1D0,-1D0/
      DATA ITH/25,35,36/
      DATA ITHC/37/
      DATA PREC/1D-2/
      DATA PI/3.141592654D0/
      DATA SR2/1.4142136D0/
      DATA KFNCHI/1000022,1000023,1000025,1000035/
      DATA KFCCHI/1000024,1000037/
 
C...COUNT THE NUMBER OF DECAY MODES
      LKNT=0
 
      XMW=PMAS(24,1)
      XMW2=XMW**2
      XMZ=PMAS(23,1)
      XMZ2=XMZ**2
      XW=1D0-XMW2/XMZ2
      TANW = SQRT(XW/(1D0-XW))
 
C...IX IS 1 - 4 DEPENDING ON SEQUENCE NUMBER
      KCIN=PYCOMP(KFIN)
      IX=1
      IF(KFIN.EQ.KFNCHI(2)) IX=2
      IF(KFIN.EQ.KFNCHI(3)) IX=3
      IF(KFIN.EQ.KFNCHI(4)) IX=4
 
      XMI=SMZ(IX)
      XMI2=XMI**2
      AXMI=ABS(XMI)
      AEM=PYALEM(XMI2)
      AS =PYALPS(XMI2)
      C1=AEM/XW
      XMI3=ABS(XMI**3)
 
      TANB=RMSS(5)
      BETA=ATAN(TANB)
      ALFA=RMSS(18)
      CBETA=COS(BETA)
      SBETA=TANB*CBETA
      CALFA=COS(ALFA)
      SALFA=SIN(ALFA)
 
C...CHECK ALL 2-BODY DECAYS TO GAUGE AND HIGGS BOSONS
      IF(IX.EQ.1.AND.IMSS(11).EQ.0) GOTO 260
 
C...FORCE CHI0_2 -> CHI0_1 + GAMMA
      IF(IX.EQ.2 .AND. IMSS(10).NE.0 ) THEN
        XMJ=SMZ(1)
        AXMJ=ABS(XMJ)
        LKNT=LKNT+1
        GAMCON=AEM**3/8D0/PI/XMW2/XW
        XMT1=(PMAS(PYCOMP(KSUSY1+6),1)/PMAS(6,1))**2
        XMT2=(PMAS(PYCOMP(KSUSY2+6),1)/PMAS(6,1))**2
        XLAM(LKNT)=PYXXGA(GAMCON,AXMI,AXMJ,XMT1,XMT2)
        IDLAM(LKNT,1)=KSUSY1+22
        IDLAM(LKNT,2)=22
        IDLAM(LKNT,3)=0
        WRITE(MSTU(11),*) 'FORCED N2 -> N1 + GAMMA ',XLAM(LKNT)
        GOTO 300
      ENDIF
 
C...GRAVITINO DECAY MODES
 
      IF(IMSS(11).EQ.1) THEN
        XMP=RMSS(29)
        IDG=39+KSUSY1
        XMGR=PMAS(PYCOMP(IDG),1)
        SINW=SQRT(XW)
        COSW=SQRT(1D0-XW)
        XFAC=(XMI2/(XMP*XMGR))**2*AXMI/48D0/PI
        IF(AXMI.GT.XMGR+PMAS(22,1)) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=22
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*(ZMIX(IX,1)*COSW+ZMIX(IX,2)*SINW)**2
        ENDIF
        IF(AXMI.GT.XMGR+XMZ) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=23
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*((ZMIX(IX,1)*SINW-ZMIX(IX,2)*COSW)**2 +
     $  .5D0*(ZMIX(IX,3)*CBETA-ZMIX(IX,4)*SBETA)**2)*(1D0-XMZ2/XMI2)**4
        ENDIF
        IF(AXMI.GT.XMGR+PMAS(25,1)) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=25
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*((ZMIX(IX,3)*SALFA-ZMIX(IX,4)*CALFA)**2)*
     $  .5D0*(1D0-PMAS(25,1)**2/XMI2)**4
        ENDIF
        IF(AXMI.GT.XMGR+PMAS(35,1)) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=35
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*((ZMIX(IX,3)*CALFA+ZMIX(IX,4)*SALFA)**2)*
     $  .5D0*(1D0-PMAS(35,1)**2/XMI2)**4
        ENDIF
        IF(AXMI.GT.XMGR+PMAS(36,1)) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=36
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*((ZMIX(IX,3)*SBETA+ZMIX(IX,4)*CBETA)**2)*
     $  .5D0*(1D0-PMAS(36,1)**2/XMI2)**4
        ENDIF
        IF(IX.EQ.1) GOTO 260
      ENDIF
 
      DO 180 IJ=1,IX-1
        XMJ=SMZ(IJ)
        AXMJ=ABS(XMJ)
        XMJ2=XMJ**2
 
C...CHI0_I -> CHI0_J + GAMMA
        IF(AXMI.GE.AXMJ.AND.SBETA/CBETA.LE.2D0) THEN
          RAT1=ZMIX(IJ,1)**2+ZMIX(IJ,2)**2
          RAT1=RAT1/( 1D-6+ZMIX(IX,3)**2+ZMIX(IX,4)**2 )
          RAT2=ZMIX(IX,1)**2+ZMIX(IX,2)**2
          RAT2=RAT2/( 1D-6+ZMIX(IJ,3)**2+ZMIX(IJ,4)**2 )
          IF((RAT1.GT. 0.90D0 .AND. RAT1.LT. 1.10D0) .OR.
     &    (RAT2.GT. 0.90D0 .AND. RAT2.LT. 1.10D0)) THEN
            LKNT=LKNT+1
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=22
            IDLAM(LKNT,3)=0
            GAMCON=AEM**3/8D0/PI/XMW2/XW
            XMT1=(PMAS(PYCOMP(KSUSY1+6),1)/PMAS(6,1))**2
            XMT2=(PMAS(PYCOMP(KSUSY2+6),1)/PMAS(6,1))**2
            XLAM(LKNT)=PYXXGA(GAMCON,AXMI,AXMJ,XMT1,XMT2)
          ENDIF
        ENDIF
 
C...CHI0_I -> CHI0_J + Z0
        IF(AXMI.GE.AXMJ+XMZ) THEN
          LKNT=LKNT+1
          GL=-0.5D0*(ZMIX(IX,3)*ZMIX(IJ,3)-ZMIX(IX,4)*ZMIX(IJ,4))
          GR=-GL
          XLAM(LKNT)=PYX2XG(C1/XMW2,XMI,XMJ,XMZ,GL,GR)
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=23
          IDLAM(LKNT,3)=0
        ELSEIF(AXMI.GE.AXMJ) THEN
          FID=11
          EI=KCHG(FID,1)/3D0
          T3=-0.5D0
          XXM(1)=0D0
          XXM(2)=XMJ
          XXM(3)=0D0
          XXM(4)=XMI
          XXM(5)=PMAS(PYCOMP(KSUSY1+11),1)
          XXM(6)=PMAS(PYCOMP(KSUSY2+11),1)
          XXM(7)=XMZ
          XXM(8)=PMAS(23,2)
          XXM(9)=-0.5D0*(ZMIX(IX,3)*ZMIX(IJ,3)-ZMIX(IX,4)*ZMIX(IJ,4))
          XXM(10)=-XXM(9)
          XXM(11)=(T3-EI*XW)/(1D0-XW)
          XXM(12)=-EI*XW/(1D0-XW)
          XXM(13)=-SR2*(T3*ZMIX(IX,2)-TANW*(T3-EI)*ZMIX(IX,1))
          XXM(14)=-SR2*(T3*ZMIX(IJ,2)-TANW*(T3-EI)*ZMIX(IJ,1))
          XXM(15)=SR2*TANW*(EI*ZMIX(IX,1))
          XXM(16)=SR2*TANW*(EI*ZMIX(IJ,1))
          S12MIN=0D0
          S12MAX=(AXMI-AXMJ)**2
 
C...CHARGED LEPTONS
          IF( XXM(5).LT.AXMI ) THEN
            XXM(5)=1D6
          ENDIF
          IF(XXM(6).LT.AXMI ) THEN
            XXM(6)=1D6
          ENDIF
          IF(AXMI.GE.AXMJ+2D0*PMAS(11,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-3)
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=11
            IDLAM(LKNT,3)=-11
            IF(AXMI.GE.AXMJ+2D0*PMAS(13,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFNCHI(IJ)
              IDLAM(LKNT,2)=13
              IDLAM(LKNT,3)=-13
            ENDIF
          ENDIF
  100     CONTINUE
          IF(ABS(SFMIX(15,1)).GT.ABS(SFMIX(15,2))) THEN
            XXM(5)=PMAS(PYCOMP(KSUSY1+15),1)
            XXM(6)=PMAS(PYCOMP(KSUSY2+15),1)
          ELSE
            XXM(6)=PMAS(PYCOMP(KSUSY1+15),1)
            XXM(5)=PMAS(PYCOMP(KSUSY2+15),1)
          ENDIF
          IF( XXM(5).LT.AXMI ) THEN
            XXM(5)=1D6
          ENDIF
          IF(XXM(6).LT.AXMI ) THEN
            XXM(6)=1D6
          ENDIF
 
          IF(AXMI.GE.AXMJ+2D0*PMAS(15,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-3)
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=15
            IDLAM(LKNT,3)=-15
          ENDIF
 
C...NEUTRINOS
  110     CONTINUE
          FID=12
          EI=KCHG(FID,1)/3D0
          T3=0.5D0
          XXM(5)=PMAS(PYCOMP(KSUSY1+12),1)
          XXM(6)=1D6
          XXM(11)=(T3-EI*XW)/(1D0-XW)
          XXM(12)=-EI*XW/(1D0-XW)
          XXM(13)=-SR2*(T3*ZMIX(IX,2)-TANW*(T3-EI)*ZMIX(IX,1))
          XXM(14)=-SR2*(T3*ZMIX(IJ,2)-TANW*(T3-EI)*ZMIX(IJ,1))
          XXM(15)=SR2*TANW*(EI*ZMIX(IX,1))
          XXM(16)=SR2*TANW*(EI*ZMIX(IJ,1))
 
          IF( XXM(5).LT.AXMI ) THEN
            XXM(5)=1D6
          ENDIF
 
          LKNT=LKNT+1
          XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-3)
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=12
          IDLAM(LKNT,3)=-12
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=14
          IDLAM(LKNT,3)=-14
  120     CONTINUE
          XXM(5)=PMAS(PYCOMP(KSUSY1+16),1)
          IF( XXM(5).LT.AXMI ) THEN
            XXM(5)=1D6
          ENDIF
          LKNT=LKNT+1
          XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-3)
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=16
          IDLAM(LKNT,3)=-16
 
C...D-TYPE QUARKS
  130     CONTINUE
          XXM(5)=PMAS(PYCOMP(KSUSY1+1),1)
          XXM(6)=PMAS(PYCOMP(KSUSY2+1),1)
          FID=1
          EI=KCHG(FID,1)/3D0
          T3=-0.5D0
 
          XXM(11)=(T3-EI*XW)/(1D0-XW)
          XXM(12)=-EI*XW/(1D0-XW)
          XXM(13)=-SR2*(T3*ZMIX(IX,2)-TANW*(T3-EI)*ZMIX(IX,1))
          XXM(14)=-SR2*(T3*ZMIX(IJ,2)-TANW*(T3-EI)*ZMIX(IJ,1))
          XXM(15)=SR2*TANW*(EI*ZMIX(IX,1))
          XXM(16)=SR2*TANW*(EI*ZMIX(IJ,1))
 
          IF( XXM(5).LT.AXMI .AND. XXM(6).LT.AXMI ) GOTO 140
          IF( XXM(5).LT.AXMI ) THEN
            XXM(5)=1D6
          ELSEIF( XXM(6).LT.AXMI ) THEN
            XXM(6)=1D6
          ENDIF
          IF(AXMI.GE.AXMJ+2D0*PMAS(1,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-3)*3D0
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=1
            IDLAM(LKNT,3)=-1
            IF(AXMI.GE.AXMJ+2D0*PMAS(3,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFNCHI(IJ)
              IDLAM(LKNT,2)=3
              IDLAM(LKNT,3)=-3
            ENDIF
          ENDIF
  140     CONTINUE
          IF(ABS(SFMIX(5,1)).GT.ABS(SFMIX(5,2))) THEN
            XXM(5)=PMAS(PYCOMP(KSUSY1+5),1)
            XXM(6)=PMAS(PYCOMP(KSUSY2+5),1)
          ELSE
            XXM(6)=PMAS(PYCOMP(KSUSY1+5),1)
            XXM(5)=PMAS(PYCOMP(KSUSY2+5),1)
          ENDIF
          IF( XXM(5).LT.AXMI .AND. XXM(6).LT.AXMI ) GOTO 150
          IF(XXM(5).LT.AXMI) THEN
            XXM(5)=1D6
          ELSEIF(XXM(6).LT.AXMI) THEN
            XXM(6)=1D6
          ENDIF
          IF(AXMI.GE.AXMJ+2D0*PMAS(5,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-3)*3D0
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=5
            IDLAM(LKNT,3)=-5
          ENDIF
 
C...U-TYPE QUARKS
  150     CONTINUE
          XXM(5)=PMAS(PYCOMP(KSUSY1+2),1)
          XXM(6)=PMAS(PYCOMP(KSUSY2+2),1)
          FID=2
          EI=KCHG(FID,1)/3D0
          T3=0.5D0
 
          XXM(11)=(T3-EI*XW)/(1D0-XW)
          XXM(12)=-EI*XW/(1D0-XW)
          XXM(13)=-SR2*(T3*ZMIX(IX,2)-TANW*(T3-EI)*ZMIX(IX,1))
          XXM(14)=-SR2*(T3*ZMIX(IJ,2)-TANW*(T3-EI)*ZMIX(IJ,1))
          XXM(15)=SR2*TANW*(EI*ZMIX(IX,1))
          XXM(16)=SR2*TANW*(EI*ZMIX(IJ,1))
 
          IF( XXM(5).LT.AXMI .AND. XXM(6).LT.AXMI ) GOTO 160
          IF(XXM(5).LT.AXMI) THEN
            XXM(5)=1D6
          ELSEIF(XXM(6).LT.AXMI) THEN
            XXM(6)=1D6
          ENDIF
          IF(AXMI.GE.AXMJ+2D0*PMAS(2,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-3)*3D0
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=2
            IDLAM(LKNT,3)=-2
            IF(AXMI.GE.AXMJ+2D0*PMAS(4,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFNCHI(IJ)
              IDLAM(LKNT,2)=4
              IDLAM(LKNT,3)=-4
            ENDIF
          ENDIF
  160     CONTINUE
        ENDIF
 
C...CHI0_I -> CHI0_J + H0_K
        EH(1)=SIN(ALFA)
        EH(2)=COS(ALFA)
        EH(3)=-SIN(BETA)
        DH(1)=COS(ALFA)
        DH(2)=-SIN(ALFA)
        DH(3)=COS(BETA)
 
        QIJ=ZMIX(IX,3)*ZMIX(IJ,2)+ZMIX(IJ,3)*ZMIX(IX,2)-
     &  TANW*(ZMIX(IX,3)*ZMIX(IJ,1)+ZMIX(IJ,3)*ZMIX(IX,1))
        RIJ=ZMIX(IX,4)*ZMIX(IJ,2)+ZMIX(IJ,4)*ZMIX(IX,2)-
     &  TANW*(ZMIX(IX,4)*ZMIX(IJ,1)+ZMIX(IJ,4)*ZMIX(IX,1))
 
        DO 170 IH=1,3
          XMH=PMAS(ITH(IH),1)
          XMH2=XMH**2
          IF(AXMI.GE.AXMJ+XMH) THEN
            LKNT=LKNT+1
            XL=PYLAMF(XMI2,XMJ2,XMH2)
            F21K=0.5D0*(QIJ*EH(IH)+RIJ*DH(IH))
            F12K=F21K
C...SIGN OF MASSES I,J
            XMK=XMJ
            IF(IH.EQ.3) XMK=-XMK
            XLAM(LKNT)=PYX2XH(C1,XMI,XMK,XMH,F12K,F21K)
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=ITH(IH)
            IDLAM(LKNT,3)=0
          ENDIF
  170   CONTINUE
  180 CONTINUE
 
C...CHI0_I -> CHI+_J + W-
      DO 220 IJ=1,2
        XMJ=SMW(IJ)
        AXMJ=ABS(XMJ)
        XMJ2=XMJ**2
        IF(AXMI.GE.AXMJ+XMW) THEN
          LKNT=LKNT+1
          GL=ZMIX(IX,2)*VMIX(IJ,1)-ZMIX(IX,4)*VMIX(IJ,2)/SR2
          GR=ZMIX(IX,2)*UMIX(IJ,1)+ZMIX(IX,3)*UMIX(IJ,2)/SR2
          XLAM(LKNT)=PYX2XG(C1/XMW2,XMI,XMJ,XMW,GL,GR)
          IDLAM(LKNT,1)=KFCCHI(IJ)
          IDLAM(LKNT,2)=-24
          IDLAM(LKNT,3)=0
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=-KFCCHI(IJ)
          IDLAM(LKNT,2)=24
          IDLAM(LKNT,3)=0
        ELSEIF(AXMI.GE.AXMJ) THEN
          S12MIN=0D0
          S12MAX=(AXMI-AXMJ)**2
          XXM(5)=ZMIX(IX,2)*VMIX(IJ,1)-ZMIX(IX,4)*VMIX(IJ,2)/SR2
          XXM(6)=ZMIX(IX,2)*UMIX(IJ,1)+ZMIX(IX,3)*UMIX(IJ,2)/SR2
 
C...LEPTONS
          FID=11
          EI=KCHG(FID,1)/3D0
          T3=-0.5D0
          XXM(7)=-SR2*(T3*ZMIX(IX,2)-TANW*(T3-EI)*ZMIX(IX,1))*UMIX(IJ,1)
          FID=12
          EI=KCHG(FID,1)/3D0
          T3=0.5D0
          XXM(8)=-SR2*(T3*ZMIX(IX,2)-TANW*(T3-EI)*ZMIX(IX,1))*VMIX(IJ,1)
 
          XXM(1)=0D0
          XXM(2)=XMJ
          XXM(3)=0D0
          XXM(4)=XMI
          XXM(9)=PMAS(24,1)
          XXM(10)=PMAS(24,2)
          XXM(11)=PMAS(PYCOMP(KSUSY1+11),1)
          XXM(12)=PMAS(PYCOMP(KSUSY1+12),1)
          IF( XXM(11).LT.AXMI .AND. XXM(12).LT.AXMI ) GOTO 190
          IF(XXM(11).LT.AXMI) THEN
            XXM(11)=1D6
          ELSEIF(XXM(12).LT.AXMI) THEN
            XXM(12)=1D6
          ENDIF
          IF(AXMI.GE.AXMJ+PMAS(11,1)+PMAS(12,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXW5,S12MIN,S12MAX,PREC)
            IDLAM(LKNT,1)=KFCCHI(IJ)
            IDLAM(LKNT,2)=11
            IDLAM(LKNT,3)=-12
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
            IF(AXMI.GE.AXMJ+PMAS(13,1)+PMAS(14,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFCCHI(IJ)
              IDLAM(LKNT,2)=13
              IDLAM(LKNT,3)=-14
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
              IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
              IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
            ENDIF
          ENDIF
  190     CONTINUE
          IF(ABS(SFMIX(15,1)).GT.ABS(SFMIX(15,2))) THEN
            XXM(11)=PMAS(PYCOMP(KSUSY1+15),1)
            XXM(12)=PMAS(PYCOMP(KSUSY1+16),1)
          ELSE
            XXM(11)=PMAS(PYCOMP(KSUSY2+15),1)
            XXM(12)=PMAS(PYCOMP(KSUSY1+16),1)
          ENDIF
 
          IF(XXM(11).LT.AXMI) THEN
            XXM(11)=1D6
          ENDIF
          IF(XXM(12).LT.AXMI) THEN
            XXM(12)=1D6
          ENDIF
          IF(AXMI.GE.AXMJ+PMAS(15,1)+PMAS(16,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXW5,S12MIN,S12MAX,PREC)
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFCCHI(IJ)
            IDLAM(LKNT,2)=15
            IDLAM(LKNT,3)=-16
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
          ENDIF
 
C...NOW, DO THE QUARKS
  200     CONTINUE
          FID=1
          EI=KCHG(FID,1)/3D0
          T3=-0.5D0
          XXM(7)=-SR2*(T3*ZMIX(IX,2)-TANW*(T3-EI)*ZMIX(IX,1))*UMIX(IJ,1)
          FID=2
          EI=KCHG(FID,1)/3D0
          T3=0.5D0
          XXM(8)=-SR2*(T3*ZMIX(IX,2)-TANW*(T3-EI)*ZMIX(IX,1))*VMIX(IJ,1)
 
          XXM(11)=PMAS(PYCOMP(KSUSY1+1),1)
          XXM(12)=PMAS(PYCOMP(KSUSY1+2),1)
          IF( XXM(11).LT.AXMI .AND. XXM(12).LT.AXMI ) GOTO 210
          IF(XXM(11).LT.AXMI) THEN
            XXM(11)=1D6
          ELSEIF(XXM(12).LT.AXMI) THEN
            XXM(12)=1D6
          ENDIF
          IF(AXMI.GE.AXMJ+PMAS(2,1)+PMAS(1,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=3D0*C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXW5,S12MIN,S12MAX,PREC)
            IDLAM(LKNT,1)=KFCCHI(IJ)
            IDLAM(LKNT,2)=1
            IDLAM(LKNT,3)=-2
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
            IF(AXMI.GE.AXMJ+PMAS(3,1)+PMAS(4,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFCCHI(IJ)
              IDLAM(LKNT,2)=3
              IDLAM(LKNT,3)=-4
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
              IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
              IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
            ENDIF
          ENDIF
  210     CONTINUE
        ENDIF
  220 CONTINUE
  230 CONTINUE
 
C...CHI0_I -> CHI+_I + H-
      DO 240 IJ=1,2
        XMJ=SMW(IJ)
        AXMJ=ABS(XMJ)
        XMJ2=XMJ**2
        XMHP=PMAS(ITHC,1)
        XMHP2=XMHP**2
        IF(AXMI.GE.AXMJ+XMHP) THEN
          LKNT=LKNT+1
          GL=CBETA*(ZMIX(IX,4)*VMIX(IJ,1)+(ZMIX(IX,2)+
     &    ZMIX(IX,1)*TANW)*VMIX(IJ,2)/SR2)
          GR=SBETA*(ZMIX(IX,3)*UMIX(IJ,1)-(ZMIX(IX,2)+
     &    ZMIX(IX,1)*TANW)*UMIX(IJ,2)/SR2)
          XLAM(LKNT)=PYX2XH(C1,XMI,XMJ,XMHP,GL,GR)
          IDLAM(LKNT,1)=KFCCHI(IJ)
          IDLAM(LKNT,2)=-ITHC
          IDLAM(LKNT,3)=0
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
          IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
          IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
        ELSE
 
        ENDIF
  240 CONTINUE
 
C...2-BODY DECAYS TO FERMION SFERMION
      DO 250 J=1,16
        IF(J.GE.7.AND.J.LE.10) GOTO 250
        KF1=KSUSY1+J
        KF2=KSUSY2+J
        XMSF1=PMAS(PYCOMP(KF1),1)
        XMSF2=PMAS(PYCOMP(KF2),1)
        XMF=PMAS(J,1)
        IF(J.LE.6) THEN
          FCOL=3D0
        ELSE
          FCOL=1D0
        ENDIF
 
        EI=KCHG(J,1)/3D0
        T3T=SIGN(1D0,EI)
        IF(J.EQ.12.OR.J.EQ.14.OR.J.EQ.16) T3T=1D0
        IF(MOD(J,2).EQ.0) THEN
          BL=T3T*ZMIX(IX,2)+TANW*ZMIX(IX,1)*(2D0*EI-T3T)
          AL=XMF*ZMIX(IX,4)/XMW/SBETA
          AR=-2D0*EI*TANW*ZMIX(IX,1)
          BR=AL
        ELSE
          BL=T3T*ZMIX(IX,2)+TANW*ZMIX(IX,1)*(2D0*EI-T3T)
          AL=XMF*ZMIX(IX,3)/XMW/CBETA
          AR=-2D0*EI*TANW*ZMIX(IX,1)
          BR=AL
        ENDIF
 
C...D~ D_L
        IF(AXMI.GE.XMF+XMSF1) THEN
          LKNT=LKNT+1
          XMA2=XMSF1**2
          XMB2=XMF**2
          XL=PYLAMF(XMI2,XMA2,XMB2)
          CA=AL*SFMIX(J,1)+AR*SFMIX(J,2)
          CB=BL*SFMIX(J,1)+BR*SFMIX(J,2)
          XLAM(LKNT)=0.5D0*FCOL*C1/8D0/XMI3*SQRT(XL)*( (XMI2+XMB2-XMA2)*
     &    (CA**2+CB**2)+4D0*CA*CB*XMF*XMI)
          IDLAM(LKNT,1)=KF1
          IDLAM(LKNT,2)=-J
          IDLAM(LKNT,3)=0
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
          IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
          IDLAM(LKNT,3)=0
        ENDIF
 
C...D~ D_R
        IF(AXMI.GE.XMF+XMSF2) THEN
          LKNT=LKNT+1
          XMA2=XMSF2**2
          XMB2=XMF**2
          CA=AL*SFMIX(J,3)+AR*SFMIX(J,4)
          CB=BL*SFMIX(J,3)+BR*SFMIX(J,4)
          XL=PYLAMF(XMI2,XMA2,XMB2)
          XLAM(LKNT)=0.5D0*FCOL*C1/8D0/XMI3*SQRT(XL)*( (XMI2+XMB2-XMA2)*
     &    (CA**2+CB**2)+4D0*CA*CB*XMF*XMI)
          IDLAM(LKNT,1)=KF2
          IDLAM(LKNT,2)=-J
          IDLAM(LKNT,3)=0
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
          IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
          IDLAM(LKNT,3)=0
        ENDIF
  250 CONTINUE
  260 CONTINUE
C...3-BODY DECAY TO Q Q~ GLUINO
      XMJ=PMAS(PYCOMP(KSUSY1+21),1)
      IF(AXMI.GE.XMJ) THEN
        AXMJ=ABS(XMJ)
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
        S12MAX=(AXMI-AXMJ)**2
C...ALL QUARKS BUT T
        XXM(11)=0D0
        XXM(12)=0D0
        XXM(13)=1D0
        XXM(14)=-SR2*(-0.5D0*ZMIX(IX,2)+TANW*ZMIX(IX,1)/6D0)
        XXM(15)=1D0
        XXM(16)=SR2*(-TANW*ZMIX(IX,1)/3D0)
        IF( XXM(5).LT.AXMI .OR. XXM(6).LT.AXMI ) GOTO 270
        IF(AXMI.GE.AXMJ+2D0*PMAS(1,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=4D0*C1*AS/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-3)
          IDLAM(LKNT,1)=KSUSY1+21
          IDLAM(LKNT,2)=1
          IDLAM(LKNT,3)=-1
          IF(AXMI.GE.AXMJ+2D0*PMAS(3,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KSUSY1+21
            IDLAM(LKNT,2)=3
            IDLAM(LKNT,3)=-3
          ENDIF
        ENDIF
  270   CONTINUE
        IF(ABS(SFMIX(5,1)).GT.ABS(SFMIX(5,2))) THEN
          XXM(5)=PMAS(PYCOMP(KSUSY1+5),1)
          XXM(6)=PMAS(PYCOMP(KSUSY2+5),1)
        ELSE
          XXM(6)=PMAS(PYCOMP(KSUSY1+5),1)
          XXM(5)=PMAS(PYCOMP(KSUSY2+5),1)
        ENDIF
        IF( XXM(5).LT.AXMI .OR. XXM(6).LT.AXMI ) GOTO 280
        IF(AXMI.GE.AXMJ+2D0*PMAS(5,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=0.5D0*C1*AS/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-3)
          IDLAM(LKNT,1)=KSUSY1+21
          IDLAM(LKNT,2)=5
          IDLAM(LKNT,3)=-5
        ENDIF
C...U-TYPE QUARKS
  280   CONTINUE
        XXM(5)=PMAS(PYCOMP(KSUSY1+2),1)
        XXM(6)=PMAS(PYCOMP(KSUSY2+2),1)
        XXM(13)=1D0
        XXM(14)=-SR2*(0.5D0*ZMIX(IX,2)+TANW*ZMIX(IX,1)/6D0)
        XXM(15)=1D0
        XXM(16)=SR2*(2D0*TANW*ZMIX(IX,1)/3D0)
        IF( XXM(5).LT.AXMI .OR. XXM(6).LT.AXMI ) GOTO 290
        IF(AXMI.GE.AXMJ+2D0*PMAS(2,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=0.5D0*C1*AS/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ5,S12MIN,S12MAX,1D-3)
          IDLAM(LKNT,1)=KSUSY1+21
          IDLAM(LKNT,2)=2
          IDLAM(LKNT,3)=-2
          IF(AXMI.GE.AXMJ+2D0*PMAS(4,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KSUSY1+21
            IDLAM(LKNT,2)=4
            IDLAM(LKNT,3)=-4
          ENDIF
        ENDIF
  290   CONTINUE
      ENDIF
 
  300 IKNT=LKNT
      XLAM(0)=0D0
      DO 310 I=1,IKNT
        IF(XLAM(I).LT.0D0) XLAM(I)=0D0
        XLAM(0)=XLAM(0)+XLAM(I)
  310 CONTINUE
      IF(XLAM(0).EQ.0D0) XLAM(0)=1D-6
 
      RETURN
      END
