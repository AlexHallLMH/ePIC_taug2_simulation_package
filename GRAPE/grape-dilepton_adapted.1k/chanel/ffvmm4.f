      SUBROUTINE FFVMM4(LIND,AL,AR,Q,AALL,AM1,AM2)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 AALL
C----- Define the anomalous tau magnetic moment factor.
      DOUBLE PRECISION a_tau
      PARAMETER (a_tau = 0.001D0)
      DIMENSION AALL(4,2,2)
      DIMENSION Q(4)
      COMMON /CHWORK/
     .        ALL11R,ALL11I,ALL12R,ALL12I,ALL21R,ALL21I,ALL22R,ALL22I,
     .        R,RR,
     .        Q1RPP1,Q1R1Y,Q1R1Z,Q2RPP1,Q2R1Y,Q2R1Z,
     .        P1RPP2,P1R2Y,P1R2Z,P2RPP2,P2R2Y,P2R2Z,
     .        J1,J2,K1,K2,
     .        LQ1,LQ2,LP1,LP2
      LOGICAL    LQ1,LQ2,LP1,LP2
      COMMON /WRKFFX/P13RPP,P13R1Y,P13R1Z,P23RPP,P23R2Y,P23R2Z,LP13,LP23
      LOGICAL    LP13,LP23

*------------------------ Entry point ----------------------------------
         ALL11R = 0.0D0
         ALL11I = 0.0D0

      DO 90 I3= 1 , 2
      DO 90 I2= 1 , 2
         AALL(LIND,I2,I3)=(0.0D0,0.0D0)
   90 CONTINUE
C- Q1,P1
           IF( LP13.AND.LP23 ) THEN
               ALL11R = P13RPP*P23RPP*(Q(4)+Q(1))
     .                 -P13RPP*(Q(2)*P23R2Y+Q(3)*P23R2Z)
     .                 -P23RPP*(Q(2)*P13R1Y+Q(3)*P13R1Z)
     .                 +(Q(4)-Q(1))*(P13R1Y*P23R2Y+P13R1Z*P23R2Z)
               ALL11I = P13RPP*(Q(2)*P23R2Z-Q(3)*P23R2Y)
     .                 +P23RPP*(Q(3)*P13R1Y-Q(2)*P13R1Z)
     .                 +(Q(4)-Q(1))*(P13R1Z*P23R2Y-P13R1Y*P23R2Z)
           ENDIF

           AALL(LIND,1,1) = AL*DCMPLX(ALL11R,-ALL11I)
           AALL(LIND,2,2) = AR*DCMPLX(ALL11R, ALL11I)


C----- Add the anomalous magnetic moment contribution.
C     The extra term is ΔΓ^μ = a_tau*(i/2)[γ^μ,γ^ν]q_ν.
C     We add anomalous contributions to all elements of AALL.
C
C     Here the dummy assignments are as follows:
C       - The factor 0.5D0 comes from the (1/2) in (i/2) [γ^μ,γ^ν]q_ν.
C       - For the (J1,K1) element, we set:
C             DeltaALL11R = 0.5 * Q(1)
C             DeltaALL11I = 0.5 * Q(2)
C         because in this simplified example we assume the real part is
C         dominated by the time-component Q(1) and the imaginary part
C         by the x-component Q(2) of q_ν.
C       - Similar dummy assignments are made for the other matrix elements.
C         In a complete derivation each ΔALL term would be the projection
C         of (i/2)[γ^μ,γ^ν]q_ν onto the corresponding helicity basis.

      IF( AM1 .GT. 1.0D0 .AND. AM2 .GT. 1.0D0 ) THEN
          DeltaALL11R = Q(1)*0.5D0
          DeltaALL11I = Q(2)*0.5D0

          DeltaALL22R = Q(3)*0.5D0
          DeltaALL22I = Q(4)*0.5D0

          DeltaALL12R = (Q(1) + Q(3))*0.5D0
          DeltaALL12I = (Q(2) + Q(4))*0.5D0

          DeltaALL21R = (Q(1) - Q(3))*0.5D0
          DeltaALL21I = (Q(2) - Q(4))*0.5D0

          AALL(LIND,1,1) = AALL(LIND,1,1)
     &                 + a_tau * DCMPLX(DeltaALL11R, DeltaALL11I)
          AALL(LIND,2,2) = AALL(LIND,2,2) 
     &                 + a_tau * DCMPLX(DeltaALL22R, DeltaALL22I)
          AALL(LIND,1,1) = AALL(LIND,1,2) 
     &                 + a_tau * DCMPLX(DeltaALL21R, DeltaALL21I)
          AALL(LIND,2,1) = AALL(LIND,2,1) 
     &                 + a_tau * DCMPLX(DeltaALL12R, DeltaALL12I)

      ENDIF

      RETURN
      END


