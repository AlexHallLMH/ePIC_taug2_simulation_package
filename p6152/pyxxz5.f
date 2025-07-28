 
C*********************************************************************
 
C...PYXXZ5
C...Calculates chi0 -> chi0 + f + ~f.
 
      FUNCTION PYXXZ5(X)
 
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
      DOUBLE PRECISION PYXXZ5,X
      DOUBLE PRECISION XM12,XM22,XM32,S,S23,S13,WPROP2
      DOUBLE PRECISION WW,WF1,WF2,WFL1,WFL2
      DOUBLE PRECISION SIJ
      DOUBLE PRECISION SR2,OL,OR,FLD,FLU,XMV,XMG,XMSU,XMSD
      DOUBLE PRECISION LE,RE,LE2,RE2,OL2,OR2,FLI,FLJ,FRI,FRJ
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
C...Integral from x to y of (t-a)/(b-t) dt.
      UTINT(X,Y,A,B)=LOG(ABS((X-A)/(B-X)*(B-Y)/(Y-A)))/(B-A)
C...Integral from x to y of 1/(t-a) dt.
      TPROP(X,Y,A)=LOG(ABS((X-A)/(Y-A)))
 
      XM12=XXM(1)**2
      XM22=XXM(2)**2
      XM32=XXM(3)**2
      S=XXM(4)**2
      S13=X
 
      S23AVE=XM22+XM32-0.5D0/X*(X+XM32-XM12)*(X+XM22-S)
      S23DEL=0.5D0/X*SQRT( ( (X-XM12-XM32)**2-4D0*XM12*XM32)*
     &( (X-XM22-S)**2  -4D0*XM22*S  ) )
 
      S23MIN=(S23AVE-S23DEL)
      S23MAX=(S23AVE+S23DEL)
 
      XMV=XXM(7)
      XMG=XXM(8)
      XMSD=XXM(5)**2
      XMSU=XXM(6)**2
      OL=XXM(9)
      OR=XXM(10)
      OL2=OL**2
      OR2=OR**2
      LE=XXM(11)
      RE=XXM(12)
      LE2=LE**2
      RE2=RE**2
      FLI=XXM(13)
      FLJ=XXM(14)
      FRI=XXM(15)
      FRJ=XXM(16)
 
      WPROP2=(S13-XMV**2)**2+(XMV*XMG)**2
      SIJ=2D0*XXM(2)*XXM(4)*S13
 
      IF(XMV.LE.1000D0) THEN
        WW=2D0*(LE2+RE2)*(OL2)*( 2D0*TINT(S23MAX,S23MIN,XM22,S)
     &  +SIJ*(S23MAX-S23MIN) )/WPROP2
        IF(XXM(5).LE.10000D0) THEN
          WFL1=2D0*FLI*FLJ*OL*LE*( 2D0*TINT2(S23MAX,S23MIN,XM22,S,XMSD)
     &    + SIJ*TPROP(S23MAX,S23MIN,XMSD) )
          WFL1=WFL1*(S13-XMV**2)/WPROP2
        ELSE
          WFL1=0D0
        ENDIF
        IF(XXM(6).LE.10000D0) THEN
          WFL2=2D0*FRI*FRJ*OR*RE*( 2D0*TINT2(S23MAX,S23MIN,XM22,S,XMSU)
     &    + SIJ*TPROP(S23MAX,S23MIN,XMSU) )
          WFL2=WFL2*(S13-XMV**2)/WPROP2
        ELSE
          WFL2=0D0
        ENDIF
      ELSE
        WW=0D0
        WFL1=0D0
        WFL2=0D0
      ENDIF
      IF(XXM(5).LE.10000D0) THEN
        WF1=0.5D0*(FLI*FLJ)**2*( 2D0*TINT3(S23MAX,S23MIN,XM22,S,XMSD)
     &  + SIJ*UTINT(S23MAX,S23MIN,XMSD,XM22+S-S13-XMSD) )
      ELSE
        WF1=0D0
      ENDIF
      IF(XXM(6).LE.10000D0) THEN
        WF2=0.5D0*(FRI*FRJ)**2*( 2D0*TINT3(S23MAX,S23MIN,XM22,S,XMSU)
     &  + SIJ*UTINT(S23MAX,S23MIN,XMSU,XM22+S-S13-XMSU) )
      ELSE
        WF2=0D0
      ENDIF
 
C...WFL1=0.0
C...WFL2=0.0
      PYXXZ5=(WW+WF1+WF2+WFL1+WFL2)
      IF(PYXXZ5.LT.0D0) THEN
        WRITE(MSTU(11),*) ' NEGATIVE WT IN PYXXZ5 '
        WRITE(MSTU(11),*) XXM(1),XXM(2),XXM(3),XXM(4)
        WRITE(MSTU(11),*) (XXM(I),I=5,8)
        WRITE(MSTU(11),*) (XXM(I),I=9,12)
        WRITE(MSTU(11),*) (XXM(I),I=13,16)
        WRITE(MSTU(11),*) WW,WF1,WF2,WFL1,WFL2
        WRITE(MSTU(11),*) S23MIN,S23MAX
        PYXXZ5=0D0
      ENDIF
 
      RETURN
      END
