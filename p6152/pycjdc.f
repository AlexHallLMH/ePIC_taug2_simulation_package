 
C*********************************************************************
 
C...PYCJDC
C...Calculate decay widths for the charginos (admixtures of
C...charged Wino and charged Higgsino.
 
C...Input:  KCIN = KF code for particle
C...Output: XLAM = widths
C...        IDLAM = KF codes for decay particles
C...        IKNT = number of decay channels defined
C...AUTHOR: STEPHEN MRENNA
C...Last change:
C...10-16-95:  force decay chi^+_1 -> chi^0_1 e+ nu_e
C...when CHIENU .NE. 0
 
      SUBROUTINE PYCJDC(KFIN,XLAM,IDLAM,IKNT)
 
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
      DOUBLE PRECISION XMFP,XMF1,XMF2,XMSL,XMG
      DOUBLE PRECISION S12MIN,S12MAX
      DOUBLE PRECISION XMI2,XMI3,XMJ2,XMH,XMH2,XMHP,XMHP2,XMA2,XMB2,XMK
      DOUBLE PRECISION PYLAMF,XL
      DOUBLE PRECISION TANW,XW,AEM,C1,AS,EI,T3,BETA,ALFA
      DOUBLE PRECISION PYX2XH,PYX2XG
      DOUBLE PRECISION XLAM(0:200)
      INTEGER IDLAM(200,3)
      INTEGER LKNT,IX,IH,J,IJ,I,IKNT,FID
      INTEGER ITH(3)
      INTEGER ITHC
      DOUBLE PRECISION ETAH(3),CH(3),DH(3),EH(3)
      DOUBLE PRECISION SR2
      DOUBLE PRECISION CBETA,SBETA,GR,GL,F12K,F21K,TANB
 
      DOUBLE PRECISION PYALEM,PI,PYALPS
      DOUBLE PRECISION AL,BL,AR,BR,ALP,BLP,ARP,BRP
      DOUBLE PRECISION CA,CB,FCOL
      INTEGER KF1,KF2,ISF
      INTEGER KFNCHI(4),KFCCHI(2)
 
      DOUBLE PRECISION TEMP
      EXTERNAL PYGAUS,PYXXZ5,PYXXW5,PYXXZ2
      DOUBLE PRECISION PYGAUS,PYXXZ5,PYXXW5,PYXXZ2
      DOUBLE PRECISION PREC
      DATA ITH/25,35,36/
      DATA ITHC/37/
      DATA ETAH/1D0,1D0,-1D0/
      DATA SR2/1.4142136D0/
      DATA PI/3.141592654D0/
      DATA PREC/1D-2/
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
 
C...1 OR 2 DEPENDING ON CHARGINO TYPE
      IX=1
      IF(KFIN.EQ.KFCCHI(2)) IX=2
      KCIN=PYCOMP(KFIN)
 
      XMI=SMW(IX)
      XMI2=XMI**2
      AXMI=ABS(XMI)
      AEM=PYALEM(XMI2)
      AS =PYALPS(XMI2)
      C1=AEM/XW
      XMI3=ABS(XMI**3)
      TANB=RMSS(5)
      BETA=ATAN(TANB)
      CBETA=COS(BETA)
      SBETA=TANB*CBETA
      ALFA=RMSS(18)
 
C...GRAVITINO DECAY MODES
 
      IF(IMSS(11).EQ.1) THEN
        XMP=RMSS(29)
        IDG=39+KSUSY1
        XMGR=PMAS(PYCOMP(IDG),1)
        SINW=SQRT(XW)
        COSW=SQRT(1D0-XW)
        XFAC=(XMI2/(XMP*XMGR))**2*AXMI/48D0/PI
        IF(AXMI.GT.XMGR+XMW) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=24
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*(.5D0*(VMIX(IX,1)**2+UMIX(IX,1)**2)+
     &  .5D0*((VMIX(IX,2)*SBETA)**2+(UMIX(IX,2)*CBETA)**2))*
     &  (1D0-XMW2/XMI2)**4
        ENDIF
        IF(AXMI.GT.XMGR+PMAS(37,1)) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=37
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*(.5D0*((VMIX(IX,2)*CBETA)**2+
     &   (UMIX(IX,2)*SBETA)**2))
     &   *(1D0-PMAS(37,1)**2/XMI2)**4
       ENDIF
      ENDIF
 
C...CHECK ALL 2-BODY DECAYS TO GAUGE AND HIGGS BOSONS
      IF(IX.EQ.1) GOTO 150
      XMJ=SMW(1)
      AXMJ=ABS(XMJ)
      XMJ2=XMJ**2
 
C...CHI_2+ -> CHI_1+ + Z0
      IF(AXMI.GE.AXMJ+XMZ) THEN
        LKNT=LKNT+1
        GL=VMIX(2,1)*VMIX(1,1)+0.5D0*VMIX(2,2)*VMIX(1,2)
        GR=UMIX(2,1)*UMIX(1,1)+0.5D0*UMIX(2,2)*UMIX(1,2)
        XLAM(LKNT)=PYX2XG(C1/XMW2,XMI,XMJ,XMZ,GL,GR)
        IDLAM(LKNT,1)=KFCCHI(1)
        IDLAM(LKNT,2)=23
        IDLAM(LKNT,3)=0
 
C...CHARGED LEPTONS
      ELSEIF(AXMI.GE.AXMJ) THEN
        XXM(5)=-(VMIX(2,1)*VMIX(1,1)+0.5D0*VMIX(2,2)*VMIX(1,2))
        XXM(6)=-(UMIX(2,1)*UMIX(1,1)+0.5D0*UMIX(2,2)*UMIX(1,2))
        XXM(9)=XMZ
        XXM(10)=PMAS(23,2)
        XXM(1)=0D0
        XXM(2)=XMJ
        XXM(3)=0D0
        XXM(4)=XMI
        S12MIN=0D0
        S12MAX=(AXMJ-AXMI)**2
        XXM(7)= (-0.5D0+XW)/(1D0-XW)
        XXM(8)= XW/(1D0-XW)
        XXM(11)=PMAS(PYCOMP(KSUSY1+12),1)
        XXM(12)=VMIX(2,1)*VMIX(1,1)
        IF( XXM(11).LT.AXMI ) THEN
          XXM(11)=1D6
        ENDIF
        IF(AXMI.GE.AXMJ+2D0*PMAS(11,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ2,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=11
          IDLAM(LKNT,3)=-11
          IF(AXMI.GE.AXMJ+2D0*PMAS(13,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFCCHI(1)
            IDLAM(LKNT,2)=13
            IDLAM(LKNT,3)=-13
            IF(AXMI.GE.AXMJ+2D0*PMAS(15,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFCCHI(1)
              IDLAM(LKNT,2)=15
              IDLAM(LKNT,3)=-15
            ENDIF
          ENDIF
        ENDIF
 
C...NEUTRINOS
  100   CONTINUE
        XXM(7)= (0.5D0)/(1D0-XW)
        XXM(8)= 0D0
        XXM(11)=PMAS(PYCOMP(KSUSY1+11),1)
        XXM(12)=UMIX(2,1)*UMIX(1,1)
        IF( XXM(11).LT.AXMI ) THEN
          XXM(11)=1D6
        ENDIF
        IF(AXMI.GE.AXMJ+2D0*PMAS(12,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ2,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=12
          IDLAM(LKNT,3)=-12
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=14
          IDLAM(LKNT,3)=-14
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=16
          IDLAM(LKNT,3)=-16
        ENDIF
 
C...D-TYPE QUARKS
  110   CONTINUE
        XXM(7)= (-0.5D0+XW/3D0)/(1D0-XW)
        XXM(8)= XW/3D0/(1D0-XW)
        XXM(11)=PMAS(PYCOMP(KSUSY1+2),1)
        XXM(12)=VMIX(2,1)*VMIX(1,1)
        IF( XXM(11).LT.AXMI ) GOTO 120
        IF(AXMI.GE.AXMJ+2D0*PMAS(1,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=3D0*C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ2,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=1
          IDLAM(LKNT,3)=-1
          IF(AXMI.GE.AXMJ+2D0*PMAS(3,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFCCHI(1)
            IDLAM(LKNT,2)=3
            IDLAM(LKNT,3)=-3
            IF(AXMI.GE.AXMJ+2D0*PMAS(5,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFCCHI(1)
              IDLAM(LKNT,2)=5
              IDLAM(LKNT,3)=-5
            ENDIF
          ENDIF
        ENDIF
 
C...U-TYPE QUARKS
  120   CONTINUE
        XXM(7)= (0.5D0-2D0*XW/3D0)/(1D0-XW)
        XXM(8)= -2D0*XW/3D0/(1D0-XW)
        XXM(11)=PMAS(PYCOMP(KSUSY1+1),1)
        XXM(12)=UMIX(2,1)*UMIX(1,1)
        IF( XXM(11).LT.AXMI ) GOTO 130
        IF(AXMI.GE.AXMJ+2D0*PMAS(2,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=3D0*C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ2,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=2
          IDLAM(LKNT,3)=-2
          IF(AXMI.GE.AXMJ+2D0*PMAS(4,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFCCHI(1)
            IDLAM(LKNT,2)=4
            IDLAM(LKNT,3)=-4
          ENDIF
        ENDIF
  130   CONTINUE
      ENDIF
 
C...CHI_2+ -> CHI_1+ + H0_K
      EH(2)=COS(ALFA)
      EH(1)=SIN(ALFA)
      EH(3)=-SBETA
      DH(2)=-SIN(ALFA)
      DH(1)=COS(ALFA)
      DH(3)=COS(BETA)
      DO 140 IH=1,3
        XMH=PMAS(ITH(IH),1)
        XMH2=XMH**2
C...NO 3-BODY OPTION
        IF(AXMI.GE.AXMJ+XMH) THEN
          LKNT=LKNT+1
          XL=PYLAMF(XMI2,XMJ2,XMH2)
          F21K=(VMIX(2,1)*UMIX(1,2)*EH(IH) -
     &    VMIX(2,2)*UMIX(1,1)*DH(IH))/SR2
          F12K=(VMIX(1,1)*UMIX(2,2)*EH(IH) -
     &    VMIX(1,2)*UMIX(2,1)*DH(IH))/SR2
          XMK=XMJ*ETAH(IH)
          XLAM(LKNT)=PYX2XH(C1,XMI,XMK,XMH,F12K,F21K)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=ITH(IH)
          IDLAM(LKNT,3)=0
        ENDIF
  140 CONTINUE
 
C...CHI1 JUMPS TO HERE
  150 CONTINUE
 
C...CHI+_I -> CHI0_J + W+
      DO 180 IJ=1,4
        XMJ=SMZ(IJ)
        AXMJ=ABS(XMJ)
        XMJ2=XMJ**2
        IF(AXMI.GE.AXMJ+XMW) THEN
          LKNT=LKNT+1
          GL=ZMIX(IJ,2)*VMIX(IX,1)-ZMIX(IJ,4)*VMIX(IX,2)/SR2
          GR=ZMIX(IJ,2)*UMIX(IX,1)+ZMIX(IJ,3)*UMIX(IX,2)/SR2
          XLAM(LKNT)=PYX2XG(C1/XMW2,XMI,XMJ,XMW,GL,GR)
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=24
          IDLAM(LKNT,3)=0
 
C...LEPTONS
        ELSEIF(AXMI.GE.AXMJ) THEN
          XMF1=0D0
          XMF2=0D0
          S12MIN=(XMF1+XMF2)**2
          S12MAX=(AXMJ-AXMI)**2
          XXM(5)=-1D0/SR2*ZMIX(IJ,4)*VMIX(IX,2)+ZMIX(IJ,2)*VMIX(IX,1)
          XXM(6)= 1D0/SR2*ZMIX(IJ,3)*UMIX(IX,2)+ZMIX(IJ,2)*UMIX(IX,1)
          FID=11
          EI=KCHG(FID,1)/3D0
          T3=-0.5D0
          XXM(7)=-SR2*(T3*ZMIX(IJ,2)-TANW*(T3-EI)*ZMIX(IJ,1))*UMIX(IX,1)
          FID=12
          EI=KCHG(FID,1)/3D0
          T3=0.5D0
          XXM(8)=-SR2*(T3*ZMIX(IJ,2)-TANW*(T3-EI)*ZMIX(IJ,1))*VMIX(IX,1)
 
          XXM(4)=XMI
          XXM(1)=XMF1
          XXM(2)=XMJ
          XXM(3)=XMF2
          XXM(9)=PMAS(24,1)
          XXM(10)=PMAS(24,2)
          XXM(11)=PMAS(PYCOMP(KSUSY1+11),1)
          XXM(12)=PMAS(PYCOMP(KSUSY1+12),1)
 
C...1/(2PI)**3*/(32*M**3)*G^4, G^2/(4*PI)= AEM/XW,
C...--> 1/(16PI)/M**3*(AEM/XW)**2
 
          IF(XXM(11).LT.AXMI) THEN
            XXM(11)=1D6
          ENDIF
          IF(XXM(12).LT.AXMI) THEN
            XXM(12)=1D6
          ENDIF
          IF(AXMI.GE.AXMJ+PMAS(11,1)+PMAS(12,1)) THEN
            LKNT=LKNT+1
            TEMP=PYGAUS(PYXXW5,S12MIN,S12MAX,PREC)
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*TEMP
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=-11
            IDLAM(LKNT,3)=12
 
C...ONLY DECAY CHI+1 -> E+ NU_E
            IF( IMSS(12).NE. 0 ) GOTO 220
            IF(AXMI.GE.AXMJ+PMAS(13,1)+PMAS(14,1)) THEN
              LKNT=LKNT+1
              XXM(11)=PMAS(PYCOMP(KSUSY1+13),1)
              XXM(12)=PMAS(PYCOMP(KSUSY1+14),1)
              IF(XXM(11).LT.AXMI) THEN
                XXM(11)=1D6
              ELSEIF(XXM(12).LT.AXMI) THEN
                XXM(12)=1D6
              ENDIF
              TEMP=PYGAUS(PYXXW5,S12MIN,S12MAX,PREC)
              XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*TEMP
              IDLAM(LKNT,1)=KFNCHI(IJ)
              IDLAM(LKNT,2)=-13
              IDLAM(LKNT,3)=14
              IF(AXMI.GE.AXMJ+PMAS(15,1)+PMAS(16,1)) THEN
                LKNT=LKNT+1
                IF(ABS(SFMIX(15,1)).GT.ABS(SFMIX(15,2))) THEN
                  XXM(11)=PMAS(PYCOMP(KSUSY1+15),1)
                ELSE
                  XXM(11)=PMAS(PYCOMP(KSUSY2+15),1)
                ENDIF
                XXM(12)=PMAS(PYCOMP(KSUSY1+16),1)
                IF(XXM(11).LT.AXMI) THEN
                  XXM(11)=1D6
                ENDIF
                IF(XXM(12).LT.AXMI) THEN
                  XXM(12)=1D6
                ENDIF
                TEMP=PYGAUS(PYXXW5,S12MIN,S12MAX,PREC)
                XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*TEMP
                IDLAM(LKNT,1)=KFNCHI(IJ)
                IDLAM(LKNT,2)=-15
                IDLAM(LKNT,3)=16
              ENDIF
            ENDIF
          ENDIF
 
C...NOW, DO THE QUARKS
  160     CONTINUE
          FID=1
          EI=KCHG(FID,1)/3D0
          T3=-0.5D0
          XXM(7)=-SR2*(T3*ZMIX(IJ,2)-TANW*(T3-EI)*ZMIX(IJ,1))*UMIX(IX,1)
          FID=1
          EI=KCHG(FID,1)/3D0
          T3=0.5D0
          XXM(8)=-SR2*(T3*ZMIX(IJ,2)-TANW*(T3-EI)*ZMIX(IJ,1))*VMIX(IX,1)
 
          XXM(11)=PMAS(PYCOMP(KSUSY1+1),1)
          XXM(12)=PMAS(PYCOMP(KSUSY1+2),1)
          IF( XXM(11).LT.AXMI .AND. XXM(12).LT.AXMI ) GOTO 170
          IF(XXM(11).LT.AXMI) THEN
            XXM(11)=1D6
          ELSEIF(XXM(12).LT.AXMI) THEN
            XXM(12)=1D6
          ENDIF
          IF(AXMI.GE.AXMJ+PMAS(1,1)+PMAS(2,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=3D0*C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXW5,S12MIN,S12MAX,PREC)
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=-1
            IDLAM(LKNT,3)=2
            IF(AXMI.GE.AXMJ+PMAS(3,1)+PMAS(4,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFNCHI(IJ)
              IDLAM(LKNT,2)=-3
              IDLAM(LKNT,3)=4
            ENDIF
          ENDIF
  170     CONTINUE
        ENDIF
  180 CONTINUE
 
C...CHI+_I -> CHI0_J + H+
      DO 190 IJ=1,4
        XMJ=SMZ(IJ)
        AXMJ=ABS(XMJ)
        XMJ2=XMJ**2
        XMHP=PMAS(ITHC,1)
        XMHP2=XMHP**2
        IF(AXMI.GE.AXMJ+XMHP) THEN
          LKNT=LKNT+1
          GL=CBETA*(ZMIX(IJ,4)*VMIX(IX,1)+(ZMIX(IJ,2)+
     &    ZMIX(IJ,1)*TANW)*VMIX(IX,2)/SR2)
          GR=SBETA*(ZMIX(IJ,3)*UMIX(IX,1)-(ZMIX(IJ,2)+
     &    ZMIX(IJ,1)*TANW)*UMIX(IX,2)/SR2)
          XLAM(LKNT)=PYX2XH(C1,XMI,XMJ,XMHP,GL,GR)
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=ITHC
          IDLAM(LKNT,3)=0
        ELSE
 
        ENDIF
  190 CONTINUE
 
C...2-BODY DECAYS TO FERMION SFERMION
      DO 200 J=1,16
        IF(J.GE.7.AND.J.LE.10) GOTO 200
        IF(MOD(J,2).EQ.0) THEN
          KF1=KSUSY1+J-1
        ELSE
          KF1=KSUSY1+J+1
        ENDIF
        KF2=KF1+KSUSY1
        XMSF1=PMAS(PYCOMP(KF1),1)
        XMSF2=PMAS(PYCOMP(KF2),1)
        XMF=PMAS(J,1)
        IF(J.LE.6) THEN
          FCOL=3D0
        ELSE
          FCOL=1D0
        ENDIF
 
C...U~ D_L
        IF(MOD(J,2).EQ.0) THEN
          XMFP=PMAS(J-1,1)
          AL=UMIX(IX,1)
          BL=-XMF*VMIX(IX,2)/XMW/SBETA/SR2
          AR=-XMFP*UMIX(IX,2)/XMW/CBETA/SR2
          BR=0D0
          ISF=J-1
        ELSE
          XMFP=PMAS(J+1,1)
          AL=VMIX(IX,1)
          BL=-XMF*UMIX(IX,2)/XMW/CBETA/SR2
          BR=0D0
          AR=-XMFP*VMIX(IX,2)/XMW/SBETA/SR2
          ISF=J+1
        ENDIF
 
C...~U_L D
        IF(AXMI.GE.XMF+XMSF1) THEN
          LKNT=LKNT+1
          XMA2=XMSF1**2
          XMB2=XMF**2
          XL=PYLAMF(XMI2,XMA2,XMB2)
          CA=AL*SFMIX(ISF,1)+AR*SFMIX(ISF,2)
          CB=BL*SFMIX(ISF,1)+BR*SFMIX(ISF,2)
          XLAM(LKNT)=FCOL*C1/8D0/XMI3*SQRT(XL)*( (XMI2+XMB2-XMA2)*
     &    (CA**2+CB**2)+4D0*CA*CB*XMF*XMI)
          IDLAM(LKNT,3)=0
          IF(MOD(J,2).EQ.0) THEN
            IDLAM(LKNT,1)=-KF1
            IDLAM(LKNT,2)=J
          ELSE
            IDLAM(LKNT,1)=KF1
            IDLAM(LKNT,2)=-J
          ENDIF
        ENDIF
 
C...U~ D_R
        IF(AXMI.GE.XMF+XMSF2) THEN
          LKNT=LKNT+1
          XMA2=XMSF2**2
          XMB2=XMF**2
          CA=AL*SFMIX(ISF,3)+AR*SFMIX(ISF,4)
          CB=BL*SFMIX(ISF,3)+BR*SFMIX(ISF,4)
          XL=PYLAMF(XMI2,XMA2,XMB2)
          XLAM(LKNT)=FCOL*C1/8D0/XMI3*SQRT(XL)*( (XMI2+XMB2-XMA2)*
     &    (CA**2+CB**2)+4D0*CA*CB*XMF*XMI)
          IDLAM(LKNT,3)=0
          IF(MOD(J,2).EQ.0) THEN
            IDLAM(LKNT,1)=-KF2
            IDLAM(LKNT,2)=J
          ELSE
            IDLAM(LKNT,1)=KF2
            IDLAM(LKNT,2)=-J
          ENDIF
        ENDIF
  200 CONTINUE
 
C...3-BODY DECAY TO Q Q~' GLUINO, ONLY IF IT CANNOT PROCEED THROUGH
C...A 2-BODY -- 2-BODY CHAIN
      XMJ=PMAS(PYCOMP(KSUSY1+21),1)
      IF(AXMI.GE.XMJ) THEN
        AXMJ=ABS(XMJ)
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
        IF( XXM(11).LT.AXMI .OR. XXM(12).LT.AXMI ) GOTO 210
        IF(AXMI.GE.AXMJ+PMAS(1,1)+PMAS(2,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=4D0*C1*AS/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXW5,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KSUSY1+21
          IDLAM(LKNT,2)=-1
          IDLAM(LKNT,3)=2
          IF(AXMI.GE.AXMJ+PMAS(3,1)+PMAS(4,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KSUSY1+21
            IDLAM(LKNT,2)=-3
            IDLAM(LKNT,3)=4
          ENDIF
        ENDIF
  210   CONTINUE
      ENDIF
 
  220 IKNT=LKNT
      XLAM(0)=0D0
      DO 230 I=1,IKNT
        XLAM(0)=XLAM(0)+XLAM(I)
        IF(XLAM(I).LT.0D0) THEN
          WRITE(MSTU(11),*) ' XLAM(I) = ',XLAM(I),KCIN,
     &    (IDLAM(I,J),J=1,3)
          XLAM(I)=0D0
        ENDIF
  230 CONTINUE
      IF(XLAM(0).EQ.0D0) THEN
        XLAM(0)=1D-6
        WRITE(MSTU(11),*) ' XLAM(0) = ',XLAM(0)
        WRITE(MSTU(11),*) LKNT
        WRITE(MSTU(11),*) (XLAM(J),J=1,LKNT)
      ENDIF
 
      RETURN
      END
