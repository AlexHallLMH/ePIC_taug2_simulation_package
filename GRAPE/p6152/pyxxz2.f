 
C*********************************************************************
 
C...PYXXZ2
C...Calculates chi+ -> chi+ + f + ~f.
 
      FUNCTION PYXXZ2(X)
 
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
      DOUBLE PRECISION PYXXZ2,X
      DOUBLE PRECISION XM12,XM22,XM32,S,S23,S13,S12,WPROP2
      DOUBLE PRECISION WW,WU,WD,WWU,WWD,WUD
      DOUBLE PRECISION SR2,OL,OR,FLD,FLU,XMV,XMG,XMSL
      DOUBLE PRECISION SIJ
      DOUBLE PRECISION LE,RE,LE2,RE2,OL2,OR2,CT
      DOUBLE PRECISION S23MIN,S23MAX,S23AVE,S23DEL
      INTEGER I
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
        PYXXZ2=0D0
        RETURN
      ENDIF
 
      XMV=XXM(9)
      XMG=XXM(10)
      XMSL=XXM(11)**2
      OL=XXM(5)
      OR=XXM(6)
      OL2=OL**2
      OR2=OR**2
      LE=XXM(7)
      RE=XXM(8)
      LE2=LE**2
      RE2=RE**2
      CT=XXM(12)
 
      WPROP2=(S13-XMV**2)**2+(XMV*XMG)**2
      SIJ=XXM(2)*XXM(4)*S13
      WW=(LE2+RE2)*(OR2+OL2)*2D0*TINT(S23MAX,S23MIN,XM22,S)
     &- 4D0*(LE2+RE2)*OL*OR*SIJ*(S23MAX-S23MIN)
      WW=WW/WPROP2
      IF(XMSL.GT.1D4*S) THEN
        WD=0D0
        WWD=0D0
      ELSE
        WD=0.5D0*CT**2*TINT3(S23MAX,S23MIN,XM22,S,XMSL)
        WWD=OL*TINT2(S23MAX,S23MIN,XM22,S,XMSL)-
     &  OR*SIJ*TPROP(S23MAX,S23MIN,XMSL)
        WWD=2D0*WWD*LE*CT*(S13-XMV**2)/WPROP2
      ENDIF
 
      PYXXZ2=(WW+WD+WWD)
      IF(PYXXZ2.LT.0D0) THEN
        WRITE(MSTU(11),*) ' NEGATIVE WT IN PYXXZ2 '
        WRITE(MSTU(11),*) WW,WD,WWD
        WRITE(MSTU(11),*) S23MIN,S23MAX
        WRITE(MSTU(11),*) (XXM(I),I=1,4)
        WRITE(MSTU(11),*) (XXM(I),I=5,8)
        WRITE(MSTU(11),*) (XXM(I),I=9,12)
        PYXXZ2=0D0
      ENDIF
 
      RETURN
      END
