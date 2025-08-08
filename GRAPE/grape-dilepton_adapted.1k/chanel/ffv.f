C
C ********************************************************************
C *                                                                  *
C *                                                                  *
C *                                                                  *
C *   SUBROUTINE FFV(L:I*4, II:I*4, I:I*4, AAM:R*8, AM:R*8,          *
C *  &               AL:R*8, AR:R*8, CC:C*16(2), C:C*16(2),          *
C *  &               Q1:R*8(4), Q2:R*8(4), P1:R*8(4), P2:R*8(4),     *
C *  &               Q:R*8(4), AALL:C*16(4,2,2))                     *
C *                                                                  *
C *     Purpose: To calculate vertex amplitudes for vector boson-    *
C *              massive fermions vertex.                            *
C *                                                                  *
C *     L=: Polarization state of vector boson.                      *
C *     I,II=: Indices to specify fermion(I,II=3) or antifermion     *
C *            (I,II=1) state.                                       *
C *     AM,AAM=: Masses of fermions.                                 *
C *     AL,AR=: Coupling constants for vertex.                       *
C *     C,CC=: Phase factors for massive fermions.                   *
C *     P1,P2,Q1,Q2=: Light-like vectors decomposed by subroutine    *
C *                   SPLTQ.                                         *
C *     Q=: Polarization vector of vector boson.                     *
C *     AALL=: Calculated results of vertex amplitudes for all       *
C *            possible helicity states.                             *
C *                                                                  *
C *                                                                  *
C *                                      written by H. Tanaka        *
C *                                                                  *
C *    Q1, Q2 : FOR UB OR VB                                         *
C *    P1, P2 : FOR U OR V                                           *
C *                                                                  *
C *                                                                  *
C ********************************************************************
C
C
C       ===========================================================
        SUBROUTINE FFV(L,II,I,AAM,AM,AL,AR,CC,C,Q1,Q2,P1,P2,Q,AALL)
C       ===========================================================
C
C
        IMPLICIT REAL*8(A-H,O-Z)
        COMPLEX*16 AALL
        COMPLEX*16 C, CC, AALL11, AALL22, AALL12, AALL21
        COMPLEX*16 AALL11_AM(2),AALL22_AM(2),AALL12_AM(2),AALL21_AM(2)
        DIMENSION Q(4), P1(4), P2(4), Q1(4), Q2(4), C(2), CC(2)
        DIMENSION AALL(4,2,2),AALL11(2),AALL22(2),AALL12(2),AALL21(2)
        INTEGER L, II, I
        REAL*8 AAM, AM, AL, AR, A_TAU
        REAL*8 R, RR
        INTEGER J1, J2, K1, K2
C
C     --- Define the anomalous coupling parameter A_TAU locally.
C         A_TAU should equal i*a_tau/(2*m) (in double precision complex).
C         For example, here we set a default value of a_tau such that:
C             A_TAU = (0.0, 1.0E-3)
C         Change the value below as needed.
        PARAMETER (A_TAU = 1.0D0)
C
C     -- First, compute the standard vertex amplitudes using FFV0:
        CALL FFV0(Q1, P1, Q, AALL11)
        CALL FFV0(Q2, P2, Q, AALL22)
        CALL FFV0(Q1, P2, Q, AALL12)
        CALL FFV0(Q2, P1, Q, AALL21)
C
C     -- Now compute the anomalous (AM) contributions using FFV0_AM:
        CALL FFV0_AM(Q1, P1, Q, AALL11_AM, A_TAU, AM)
        CALL FFV0_AM(Q2, P2, Q, AALL22_AM, A_TAU, AAM)
        CALL FFV0_AM(Q1, P2, Q, AALL12_AM, A_TAU, AM)
        CALL FFV0_AM(Q2, P1, Q, AALL21_AM, A_TAU, AAM)
C
C     -- Define scaling factors R and RR based on the indices I and II:
        R  = DFLOAT(I  - 2)
        RR = DFLOAT(II - 2)
        IF (AM .LE. 0.0D0) THEN
          R = 0.0D0
        ENDIF
        IF (AAM .LE. 0.0D0) THEN
          RR = 0.0D0
        ENDIF
C
C     -- Determine helicity indices:
        J1 = (5 - II) / 2
        J2 = 3 - J1
        K1 = (5 - I) / 2
        K2 = 3 - K1
C
C     -- Combine the standard and anomalous amplitudes.
C        For each helicity combination we add the standard result and the
C        anomalous piece computed by FFV0_AM.
        AALL(L,J1,K1) = 100 * (AL * AALL11(1))   + (AR * R * RR * C(2) * CC(2) * AALL22(2))
     &               + (AL * AALL11_AM(1))
     &               + (AR*R*RR*C(2)*CC(2)*AALL22_AM(2))
C
        AALL(L,J2,K2) = 100* (AR * AALL11(2))   + (AL * R * RR * C(1) * CC(1) * AALL22(1))
     &               + (AR * AALL11_AM(2))
     &               + (AL * R *RR*C(1) * CC(1) * AALL22_AM(1))
C
        AALL(L,J1,K2) = - (AR * RR * CC(2) * AALL21(2)) - (AL * R * C(1) * AALL12(1))
     &               - (AR * RR * CC(2) * AALL21_AM(2)) 
     &               - (AL * R *C(1)* AALL12_AM(1))
C
        AALL(L,J2,K1) = -(AL * RR * CC(1) * AALL21(1)) - (AR * R * C(2) * AALL12(2))
     &               - (AL * RR * CC(1) * AALL21_AM(1)) 
     &               - (AR * R * C(2) * AALL12_AM(2))
C
        RETURN
      END