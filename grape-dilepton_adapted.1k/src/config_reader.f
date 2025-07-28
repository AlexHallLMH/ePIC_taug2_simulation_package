      SUBROUTINE READ_CONFIG
      IMPLICIT REAL*8(A-H,O-Z)

      DOUBLE PRECISION a_tau
      COMMON /TAU_CONFIG/ a_tau

      a_tau = 0.00D0  ! Default value

      OPEN(UNIT=10, FILE='config.txt', STATUS='OLD', FORM='FORMATTED', ERR=100)
      READ(10,*,ERR=100) a_tau
 100  CLOSE(10)

      RETURN
      END
