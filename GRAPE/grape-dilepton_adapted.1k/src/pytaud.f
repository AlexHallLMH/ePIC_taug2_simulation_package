      SUBROUTINE PYTAUD
      IMPLICIT NONE

C      ------------------------------------------------------------
C      PYTAUD: Tau decay routine for PYTHIA6 using TAUOLA and TauSpinner.
C      It processes taus from the /PYJETS/ common block (HEPEVT style),
C      decaying them with TAUOLA and reweighting for spin correlations.
C      ------------------------------------------------------------

      INTEGER N, I
      INTEGER K(4000,5)
      REAL*8 P(4000,5), V(4000,5)
      COMMON /PYJETS/ N, K, P, V

C      Loop over the event record to process all taus (PDG code ±15, status = 1)
      DO I = 1, N
         IF (ABS(K(I,2)) .EQ. 15 .AND. K(I,1) .EQ. 1) THEN
C           Call TAUOLA’s decay routine.
            CALL DEKAY(P(I,1), P(I,2), P(I,3), P(I,4), P(I,5))
C           Apply spin reweighting via TauSpinner
         ENDIF
      ENDDO

      RETURN
      END
