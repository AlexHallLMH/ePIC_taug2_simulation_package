 
C*********************************************************************
 
C...PYXXW5
C...Calculates chi0(+) -> chi+(0) + f + ~f'.
 
      FUNCTION PYXXW5(X)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KEXCIT=4000000)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYINTS/XXM(20)
      SAVE /PYDAT1/,/PYINTS/
 
C...Local variables.
      DOUBLE PRECISION PYXXW5,X
      DOUBLE PRECISION XM12,XM22,XM32,S,S23,S13,S12,WPROP2
      DOUBLE PRECISION WW,WU,WD,WWU,WWD,WUD
      DOUBLE PRECISION SR2,OL,OR,FLD,FLU,XMV,XMG,XMSD,XMSU
      DOUBLE PRECISION SIJ
      DOUBLE PRECISION S23MIN,S23MAX,S23AVE,S23DEL
      INTEGER IK
      SAVE IK
      DATA IK/0/
      DATA SR2/1.4142136D0/
 
C...Statement functions.
C...Integral from x to y of (t-a)(b-t) dt.
      TINT(X,Y,A,B)=(X-Y)*(-(X**2+X*Y+Y**2)/3D0+(B+A)*(X+Y)/2D0-A*B)
C...Integral from x to y of (t-a)(b-t)/(t-c) dt.
      TINT2(X,Y,A,B,C)=(X-Y)*(-0.5D0*(X+Y)+(B+A-C))-
     &LOG(ABS((X-C)/(Y-C)))*(C-B)*(C-A)
C...Integral from x to y of (t-a)(b-t)/(t-c)**2 dt.
      TINT3(X,Y,A,B,C)=-(X-Y)+(C-A)*(C-B)*(Y-X)/(X-C)/(Y-C)+
     &(B+A-2D0*C)*LOG(ABS((X-C)/(Y-C)))
C...Integral from x to y of (t-a)/(b-t) dt.
      UTINT(X,Y,A,B)=LOG(ABS((X-A)/(B-X)*(B-Y)/(Y-A)))/(B-A)
C...Integral from x to y of 1/(t-a) dt.
      TPROP(X,Y,A)=LOG(ABS((X-A)/(Y-A)))
 
      XM12=XXM(1)**2
      XM22=XXM(2)**2
      XM32=XXM(3)**2
      S=XXM(4)**2
      S13=X
      IF(XXM(1).EQ.0.AND.XXM(3).EQ.0D0) THEN
        S23AVE=0.5D0*(XM22+S-S13)
        S23DEL=0.5D0*SQRT( (X-XM22-S)**2-4D0*XM22*S )
      ELSE
        S23AVE=XM22+XM32-0.5D0/X*(X+XM32-XM12)*(X+XM22-S)
        S23DEL=0.5D0/X*SQRT( ( (X-XM12-XM32)**2-4D0*XM12*XM32)*
     &  ( (X-XM22-S)**2  -4D0*XM22*S  ) )
      ENDIF
      S23MIN=(S23AVE-S23DEL)
      S23MAX=(S23AVE+S23DEL)
      IF(S23DEL.LT.1D-3) THEN
        PYXXW5=0D0
        RETURN
      ENDIF
      XMV=XXM(9)
      XMG=XXM(10)
      XMSD=XXM(11)**2
      XMSU=XXM(12)**2
      OL=XXM(5)
      OR=XXM(6)
      FLD=XXM(7)
      FLU=XXM(8)
 
      WPROP2=((S13-XMV**2)**2+(XMV*XMG)**2)
      SIJ=S13*XXM(2)*XXM(4)
      IF(XMV.LE.1000D0) THEN
        WW=(OR**2+OL**2)*TINT(S23MAX,S23MIN,XM22,S)
     &  -2D0*OL*OR*SIJ*(S23MAX-S23MIN)
        WW=WW/WPROP2
        IF(XXM(11).LE.10000D0) THEN
          WWD=OL*SIJ*TPROP(S23MAX,S23MIN,XMSD)
     &    -OR*TINT2(S23MAX,S23MIN,XM22,S,XMSD)
          WWD=-WWD*SR2*FLD
          WWD=WWD*(S13-XMV**2)/WPROP2
        ELSE
          WWD=0D0
        ENDIF
        IF(XXM(12).LE.10000D0) THEN
          WWU=OR*SIJ*TPROP(S23MAX,S23MIN,XMSU)
     &    -OL*TINT2(S23MAX,S23MIN,XM22,S,XMSU)
          WWU=WWU*SR2*FLU
          WWU=WWU*(S13-XMV**2)/WPROP2
        ELSE
          WWU=0D0
        ENDIF
      ELSE
        WW=0D0
        WWD=0D0
        WWU=0D0
      ENDIF
      IF(XXM(12).LE.10000D0) THEN
        WU=0.5D0*FLU**2*TINT3(S23MAX,S23MIN,XM22,S,XMSU)
      ELSE
        WU=0D0
      ENDIF
      IF(XXM(11).LE.10000D0) THEN
        WD=0.5D0*FLD**2*TINT3(S23MAX,S23MIN,XM22,S,XMSD)
      ELSE
        WD=0D0
      ENDIF
      IF(XXM(11).LE.10000D0.AND.XXM(12).LE.10000D0) THEN
        WUD=FLU*FLD*SIJ*UTINT(S23MAX,S23MIN,XMSD,XM22+S-S13-XMSU)
      ELSE
        WUD=0D0
      ENDIF
 
      PYXXW5=WW+WU+WD+WWU+WWD+WUD
 
      IF(PYXXW5.LT.0D0) THEN
        IF(IK.EQ.0) THEN
          WRITE(MSTU(11),*) ' NEGATIVE WT IN PYXXW5 '
          WRITE(MSTU(11),*) WW,WU,WD
          WRITE(MSTU(11),*) WWD,WWU,WUD
          WRITE(MSTU(11),*) SQRT(S13)
          WRITE(MSTU(11),*) TINT(S23MAX,S23MIN,XM22,S)
          IK=1
        ENDIF
        PYXXW5=0D0
      ENDIF
 
      RETURN
      END
