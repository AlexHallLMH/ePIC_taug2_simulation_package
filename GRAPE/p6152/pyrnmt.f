 
C*********************************************************************
 
C...PYRNMT
C...Determines the running mass of the top quark.
 
      FUNCTION PYRNMT(XMT)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblock.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      SAVE /PYMSSM/
 
C...Local variables.
      DOUBLE PRECISION XMT
      DOUBLE PRECISION PI,R
      DOUBLE PRECISION TOL
      EXTERNAL PYALPS
      DOUBLE PRECISION PYALPS
      DATA TOL/0.001D0/
      DATA PI,R/3.141592654D0,0.61803399D0/
 
      C=1D0-R
 
      BX=XMT
      AX=MIN(50D0,BX*0.5D0)
      CX=MAX(300D0,2D0*BX)
 
      X0=AX
      X3=CX
      IF(ABS(CX-BX).GT.ABS(BX-AX))THEN
        X1=BX
        X2=BX+C*(CX-BX)
      ELSE
        X2=BX
        X1=BX-C*(BX-AX)
      ENDIF
      AS1=PYALPS(X1**2)/PI
      F1=ABS(XMT/(1D0+4D0/3D0*AS1+11D0*AS1**2)-X1)
      AS2=PYALPS(X2**2)/PI
      F2=ABS(XMT/(1D0+4D0/3D0*AS2+11D0*AS2**2)-X2)
  100 IF(ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2))) THEN
        IF(F2.LT.F1) THEN
          X0=X1
          X1=X2
          X2=R*X1+C*X3
          F1=F2
          AS2=PYALPS(X2**2)/PI
          F2=ABS(XMT/(1D0+4D0/3D0*AS2+11D0*AS2**2)-X2)
        ELSE
          X3=X2
          X2=X1
          X1=R*X2+C*X0
          F2=F1
          AS1=PYALPS(X1**2)/PI
          F1=ABS(XMT/(1D0+4D0/3D0*AS1+11D0*AS1**2)-X1)
        ENDIF
        GOTO 100
      ENDIF
      IF(F1.LT.F2) THEN
        PYRNMT=X1
        XMIN=X1
      ELSE
        PYRNMT=X2
        XMIN=X2
      ENDIF
 
      RETURN
      END
