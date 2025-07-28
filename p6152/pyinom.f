 
C*********************************************************************
 
C...PYINOM
C...Finds the mass eigenstates and mixing matrices for neutralinos
C...and charginos.
 
      SUBROUTINE PYINOM
 
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
      SAVE /PYDAT1/,/PYDAT2/,/PYMSSM/,/PYSSMT/
 
C...Local variables.
      DOUBLE PRECISION XMW,XMZ
      DOUBLE PRECISION AR(4,4),WR(4),ZR(4,4)
      DOUBLE PRECISION ZP(4,4)
      DOUBLE PRECISION DETX,XI(2,2)
      DOUBLE PRECISION XXX,YYY,XMH,XML
      DOUBLE PRECISION COSW,SINW
      DOUBLE PRECISION XMU
      DOUBLE PRECISION TERMB,TERMC,DISCR,XMH2,XML2
      DOUBLE PRECISION TANB,AL,BE,COSA,COSB,SINA,SINB,XW
      DOUBLE PRECISION XM1,XM2,XM3,BETA
      DOUBLE PRECISION Q2,AEM,A1,A2,A3,AQ,RM1,RM2
      DOUBLE PRECISION ARG,X0,X1,AX0,AX1,AT,BT
      DOUBLE PRECISION Y0,Y1,AMGX0,AM1X0,AMGX1,AM1X1
      DOUBLE PRECISION ARGX0,AR1X0,ARGX1,AR1X1
      DOUBLE PRECISION PYALPS,PYALEM
      DOUBLE PRECISION PYRNM3
      INTEGER IERR,INDEX(4),I,J,K,L,IOPT,ILR,KFNCHI(4)
      DATA KFNCHI/1000022,1000023,1000025,1000035/
 
      IOPT=IMSS(2)
      IF(IMSS(1).EQ.2) THEN
        IOPT=1
      ENDIF
C...M1, M2, AND M3 ARE INDEPENDENT
      IF(IOPT.EQ.0) THEN
        XM1=RMSS(1)
        XM2=RMSS(2)
        XM3=RMSS(3)
      ELSEIF(IOPT.GE.1) THEN
        Q2=PMAS(23,1)**2
        AEM=PYALEM(Q2)
        A2=AEM/PARU(102)
        A1=AEM/(1D0-PARU(102))
        XM1=RMSS(1)
        XM2=RMSS(2)
        IF(IMSS(1).EQ.2) XM1=RMSS(1)/RMSS(20)*A1*5D0/3D0
        IF(IOPT.EQ.1) THEN
          XM2=XM1*A2/A1*3D0/5D0
          RMSS(2)=XM2
        ELSEIF(IOPT.EQ.3) THEN
          XM1=XM2*5D0/3D0*A1/A2
          RMSS(1)=XM1
        ENDIF
        XM3=PYRNM3(XM2/A2)
        RMSS(3)=XM3
        IF(XM3.LE.0D0) THEN
          WRITE(MSTU(11),*) ' ERROR WITH M3 = ',XM3
          STOP
        ENDIF
      ENDIF
 
C...GLUINO MASS
      IF(IMSS(3).EQ.1) THEN
        PMAS(PYCOMP(KSUSY1+21),1)=XM3
      ELSE
        AQ=0D0
        DO 110 I=1,4
          DO 100 ILR=1,2
            RM1=PMAS(PYCOMP(ILR*KSUSY1+I),1)**2/XM3**2
            AQ=AQ+0.5D0*((2D0-RM1)*(RM1*LOG(RM1)-1D0)
     &      +(1D0-RM1)**2*LOG(ABS(1D0-RM1)))
  100     CONTINUE
  110   CONTINUE
 
        DO 130 I=5,6
          DO 120 ILR=1,2
            RM1=PMAS(PYCOMP(ILR*KSUSY1+I),1)**2/XM3**2
            RM2=PMAS(I,1)**2/XM3**2
            ARG=(RM1-RM2-1D0)**2-4D0*RM2**2
            IF(ARG.GE.0D0) THEN
              X0=0.5D0*(1D0+RM2-RM1-SQRT(ARG))
              AX0=ABS(X0)
              X1=0.5D0*(1D0+RM2-RM1+SQRT(ARG))
              AX1=ABS(X1)
              IF(X0.EQ.1D0) THEN
                AT=-1D0
                BT=0.25D0
              ELSEIF(X0.EQ.0D0) THEN
                AT=0D0
                BT=-0.25D0
              ELSE
                AT=0.5D0*LOG(ABS(1D0-X0))*(1D0-X0**2)+
     &          0.5D0*X0**2*LOG(AX0)
                BT=(-1D0-2D0*X0)/4D0
              ENDIF
              IF(X1.EQ.1D0) THEN
                AT=-1D0+AT
                BT=0.25D0+BT
              ELSEIF(X1.EQ.0D0) THEN
                AT=0D0+AT
                BT=-0.25D0+BT
              ELSE
                AT=0.5D0*LOG(ABS(1D0-X1))*(1D0-X1**2)+0.5D0*
     &          X1**2*LOG(AX1)+AT
                BT=(-1D0-2D0*X1)/4D0+BT
              ENDIF
              AQ=AQ+AT+BT
            ELSE
              X0=0.5D0*(1D0+RM2-RM1)
              Y0=-0.5D0*SQRT(-ARG)
              AMGX0=SQRT(X0**2+Y0**2)
              AM1X0=SQRT((1D0-X0)**2+Y0**2)
              ARGX0=ATAN2(-X0,-Y0)
              AR1X0=ATAN2(1D0-X0,Y0)
              X1=X0
              Y1=-Y0
              AMGX1=AMGX0
              AM1X1=AM1X0
              ARGX1=ATAN2(-X1,-Y1)
              AR1X1=ATAN2(1D0-X1,Y1)
              AT=0.5D0*LOG(AM1X0)*(1D0-X0**2+3D0*Y0**2)
     &        +0.5D0*(X0**2-Y0**2)*LOG(AMGX0)
              BT=(-1D0-2D0*X0)/4D0+X0*Y0*( AR1X0-ARGX0 )
              AT=AT+0.5D0*LOG(AM1X1)*(1D0-X1**2+3D0*Y1**2)
     &        +0.5D0*(X1**2-Y1**2)*LOG(AMGX1)
              BT=BT+(-1D0-2D0*X1)/4D0+X1*Y1*( AR1X1-ARGX1 )
              AQ=AQ+AT+BT
            ENDIF
  120     CONTINUE
  130   CONTINUE
        PMAS(PYCOMP(KSUSY1+21),1)=XM3*(1D0+PYALPS(XM3**2)/(2D0*PARU(2))*
     &  (15D0+AQ))
      ENDIF
 
C...NEUTRALINO MASSES
      XMZ=PMAS(23,1)
      XMW=PMAS(24,1)
      XMU=RMSS(4)
      SINW=SQRT(PARU(102))
      COSW=SQRT(1D0-PARU(102))
      TANB=RMSS(5)
      BETA=ATAN(TANB)
      COSB=COS(BETA)
      SINB=TANB*COSB
      AR(1,1) = XM1
      AR(2,2) = XM2
      AR(3,3) = 0D0
      AR(4,4) = 0D0
      AR(1,2) = 0D0
      AR(2,1) = 0D0
      AR(1,3) = -XMZ*SINW*COSB
      AR(3,1) = AR(1,3)
      AR(1,4) = XMZ*SINW*SINB
      AR(4,1) = AR(1,4)
      AR(2,3) = XMZ*COSW*COSB
      AR(3,2) = AR(2,3)
      AR(2,4) = -XMZ*COSW*SINB
      AR(4,2) = AR(2,4)
      AR(3,4) = -XMU
      AR(4,3) = -XMU
      CALL PYEIG4(AR,WR,ZR)
      DO 150 I=1,4
        SMZ(I)=WR(I)
        PMAS(PYCOMP(KFNCHI(I)),1)=ABS(SMZ(I))
        DO 140 J=1,4
          ZMIX(I,J)=ZR(I,J)
          IF(ABS(ZMIX(I,J)).LT.1D-6) ZMIX(I,J)=0D0
  140   CONTINUE
  150 CONTINUE
 
C...CHARGINO MASSES
      AR(1,1) = XM2
      AR(2,2) = XMU
      AR(1,2) = SQRT(2D0)*XMW*SINB
      AR(2,1) = SQRT(2D0)*XMW*COSB
      TERMB=AR(1,1)**2+AR(2,2)**2+AR(1,2)**2+AR(2,1)**2
      TERMC=(AR(1,1)**2-AR(2,2)**2)**2+(AR(1,2)**2-AR(2,1)**2)**2
      TERMC=TERMC+2D0*(AR(1,1)**2+AR(2,2)**2)*
     &(AR(1,2)**2+AR(2,1)**2)+
     &8D0*AR(1,1)*AR(2,2)*AR(1,2)*AR(2,1)
      DISCR=TERMC
      IF(DISCR.LT.0D0) THEN
        WRITE(MSTU(11),*) ' PROBLEM WITH DISCR '
      ELSE
        DISCR=SQRT(DISCR)
      ENDIF
      XML2=0.5D0*(TERMB-DISCR)
      XMH2=0.5D0*(TERMB+DISCR)
      XML=SQRT(XML2)
      XMH=SQRT(XMH2)
      PMAS(PYCOMP(KSUSY1+24),1)=XML
      PMAS(PYCOMP(KSUSY1+37),1)=XMH
      SMW(1)=XML
      SMW(2)=XMH
      XXX=AR(1,1)**2+AR(2,1)**2
      YYY=AR(1,1)*AR(1,2)+AR(2,2)*AR(2,1)
      VMIX(2,2) = YYY/SQRT(YYY**2+(XML2-XXX)**2)
      VMIX(1,1) = SIGN(VMIX(2,2),AR(1,1)*AR(2,2)-0.5D0*AR(1,2)**2)
      VMIX(2,1) = -(XML2-XXX)/SQRT(YYY**2+(XML2-XXX)**2)
      VMIX(1,2) = -SIGN(VMIX(2,1),AR(1,1)*AR(2,2)-0.5D0*AR(1,2)**2)
      ZR(1,1) = XML
      ZR(1,2) = 0D0
      ZR(2,1) = 0D0
      ZR(2,2) = XMH
      DETX = AR(1,1)*AR(2,2)-AR(1,2)*AR(2,1)
      XI(1,1) = AR(2,2)/DETX
      XI(2,2) = AR(1,1)/DETX
      XI(1,2) = -AR(1,2)/DETX
      XI(2,1) = -AR(2,1)/DETX
      DO 190 I=1,2
        DO 180 J=1,2
          UMIX(I,J)=0D0
          DO 170 K=1,2
            DO 160 L=1,2
              UMIX(I,J)=UMIX(I,J)+ZR(I,K)*VMIX(K,L)*XI(L,J)
  160       CONTINUE
  170     CONTINUE
  180   CONTINUE
  190 CONTINUE
 
      RETURN
      END
