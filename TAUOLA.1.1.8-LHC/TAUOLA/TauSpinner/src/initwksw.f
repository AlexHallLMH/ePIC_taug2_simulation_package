
      SUBROUTINE INITWKSWDELT(mode,IDEX,IDFX,SVAR,SWSQEFF, DELTSQ,  DeltV, GMU, ALPHAINV,  AMZi, GAMMZi, KEYGSW,
     &ReGSW1,CImGSW1,ReGSW2,CImGSW2,ReGSW3,CImGSW3,ReGSW4,CImGSW4,ReGSW6,CImGSW6 )
      

! initialization routine coupling masses etc., explicitly varying SWSQ
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / T_BEAMPM / ENE ,AMIN,AMFIN,IDE,IDF
      REAL*8              ENE ,AMIN,AMFIN
      COMMON / T_GAUSPM /SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI0  ,XUPZI0  ,XUPGF0  ,XUPZF0  
     &                  ,NDIAG0,NDIAGA,KEYA,KEYZ
     &                  ,ITCE,JTCE,ITCF,JTCF,KOLOR
      REAL*8            SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI0(2),XUPZI0(2),XUPGF0(2),XUPZF0(2)
      COMMON / T_GAUSPM1/VVcor, ZetVPi, GamVPi
     &                  ,XUPGI   ,XUPZI   ,XUPGF   ,XUPZF
      COMPLEX*16         VVcor, ZetVPi, GamVPi
      COMPLEX*16         XUPGI(2),XUPZI(2),XUPGF(2),XUPZF(2)
      
      COMMON / T_GSWPRMn /SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
      REAL*8             SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
      COMMON / T_EWn    / GMUn, ALPHAINVn
      REAL*8              GMUn, ALPHAINVn
      COMPLEX *16 GSW(10)
      REAL*8 PI     
      DATA PI /3.141592653589793238462643D0/
      GSW(1) = DCMPLX(ReGSW1,CImGSW1)
      GSW(2) = DCMPLX(ReGSW2,CImGSW2)
      GSW(3) = DCMPLX(ReGSW3,CImGSW3)
      GSW(4) = DCMPLX(ReGSW4,CImGSW4)
      ! GSW(5) out
      GSW(6) = DCMPLX(ReGSW6,CImGSW6)
   
C      PRINT *, ' initwksw GSW = ', SWSQEFF,  ReGSW1, CImGSW1, ReGSW2, CImGSW2, ReGSW6, CImGSW6
      
C     SWSQ        = sin2 (theta Weinberg)
C     AMW,AMZ     = W & Z boson masses respectively
C     AMH         = the Higgs mass
C     AMTOP       = the top mass
C     GAMMZ       = Z0 width
C
      ENE=SQRT(SVAR)/2
      AMIN=0.511D-3
      SWSQ=SWSQEFF
      AMZ=AMZi !91.1887
      GAMMZ=GAMMZi              !2.4952
      GMUn=GMU
      ALPHAINVn=ALPHAINV
       
      
C      Gfermi=1.16639d-5
      Gfermi=GMU

                ZetVPi =  Gfermi *AMZ**2 *ALPHAINV /(DSQRT(2.d0)*8.d0*PI)
     $          *(SWSQ*(1d0-SWSQ)) *16d0 
     $     * GSW(1)
C     updated following KK2f_defaults
C      IF( KEYGSW.NE.0) THEN
C         GAMMZ=2.50072032     
C      ENDIF



      
      GamVPi = 1d0   /(2d0-GSW(6))

C      PRINT *, ' initwksw ZetVPi, GamVPi = ', GSW(1), ZetVPi, GamVPi 
      
     
      IF  (IDFX.EQ. 11) then       
        IDF=2  ! denotes tau +2 tau-
        AMFIN=0.511D-3 !this mass is irrelevant if small, used in ME only
      ELSEIF (IDFX.EQ.-11) then
        IDF=-2  ! denotes tau -2 tau-
        AMFIN=0.511D-3 !this mass is irrelevant if small, used in ME only
      ELSEIF  (IDFX.EQ. 15) then       
        IDF=2  ! denotes tau +2 tau-
        AMFIN=1.77703 !this mass is irrelevant if small, used in ME only
      ELSEIF (IDFX.EQ.-15) then
        IDF=-2  ! denotes tau -2 tau-
        AMFIN=1.77703 !this mass is irrelevant if small, used in ME only
      ELSE
        WRITE(*,*) 'INITWKSW: WRONG IDFX'
        STOP
      ENDIF

      IF     (IDEX.EQ. 11) then      !electron
        IDE= 2
        AMIN=0.511D-3
      ELSEIF (IDEX.EQ.-11) then      !positron
        IDE=-2
        AMIN=0.511D-3
      ELSEIF (IDEX.EQ. 13) then      !mu+
        IDE= 2
        AMIN=0.105659
      ELSEIF (IDEX.EQ.-13) then      !mu-
        IDE=-2
        AMIN=0.105659
      ELSEIF (IDEX.EQ.  1) then      !d
        IDE= 4
        AMIN=0.05
      ELSEIF (IDEX.EQ.- 1) then      !d~
        IDE=-4
        AMIN=0.05
      ELSEIF (IDEX.EQ.  2) then      !u
        IDE= 3
        AMIN=0.02
      ELSEIF (IDEX.EQ.- 2) then      !u~
        IDE=-3
        AMIN=0.02
      ELSEIF (IDEX.EQ.  3) then      !s
        IDE= 4
        AMIN=0.3
      ELSEIF (IDEX.EQ.- 3) then      !s~
        IDE=-4
        AMIN=0.3
      ELSEIF (IDEX.EQ.  4) then      !c
        IDE= 3
        AMIN=1.3
      ELSEIF (IDEX.EQ.- 4) then      !c~
        IDE=-3
        AMIN=1.3
      ELSEIF (IDEX.EQ.  5) then      !b
        IDE= 4
        AMIN=4.5
      ELSEIF (IDEX.EQ.- 5) then      !b~
        IDE=-4
        AMIN=4.5
      ELSEIF (IDEX.EQ.  12) then     !nu_e
        IDE= 1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.- 12) then     !nu_e~
        IDE=-1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.  14) then     !nu_mu
        IDE= 1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.- 14) then     !nu_mu~
        IDE=-1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.  16) then     !nu_tau
        IDE= 1
        AMIN=0.1D-3
      ELSEIF (IDEX.EQ.- 16) then     !nu_tau~
        IDE=-1
        AMIN=0.1D-3

      ELSE
        WRITE(*,*) 'INITWKSW: WRONG IDEX'
        STOP
      ENDIF

C ----------------------------------------------------------------------
C
C     INITIALISATION OF COUPLING CONSTANTS AND FERMION-GAMMA / Z0 VERTEX
C
C     called by : KORALZ
C ----------------------------------------------------------------------
      ITCE=IDE/IABS(IDE)
      JTCE=(1-ITCE)/2
      ITCF=IDF/IABS(IDF)
      JTCF=(1-ITCF)/2
      CALL T_GIVIZO( IDE, 1,AIZOR,QE,KDUMM)
      CALL T_GIVIZO( IDE,-1,AIZOL,QE,KDUMM)
      XUPGI(1)=QE
      XUPGI(2)=QE
      T3E    = (AIZOL+AIZOR)/2.
      XUPZI(1)=(AIZOR-QE*(SWSQ+DeltSQ)*GSW(3)-QE*DeltV)/SQRT(SWSQ*(1-SWSQ))
      XUPZI(2)=(AIZOL-QE*(SWSQ+DeltSQ)*GSW(3)-QE*DeltV)/SQRT(SWSQ*(1-SWSQ))
      Ve      =(XUPZI(1)+XUPZI(2))/2.
      CALL T_GIVIZO( IDF, 1,AIZOR,QF,KOLOR)
      CALL T_GIVIZO( IDF,-1,AIZOL,QF,KOLOR)
      XUPGF(1)=QF
      XUPGF(2)=QF
      T3F    =  (AIZOL+AIZOR)/2.
      XUPZF(1)=(AIZOR-QF*(SWSQ+DeltSQ)*GSW(2)-QF*DeltV)/SQRT(SWSQ*(1-SWSQ))
      XUPZF(2)=(AIZOL-QF*(SWSQ+DeltSQ)*GSW(2)-QF*DeltV)/SQRT(SWSQ*(1-SWSQ))
      Vf      =(XUPZF(1)+XUPZF(2))/2.

* Coupling costants times EW form-factors
      Deno   = DSQRT(SWSQ*(1d0-SWSQ))
  !       Ve     = (2*T3e -4*Qe*m_Sw2*CorEle)/Deno
  !       Vf     = (2*T3f -4*Qf*m_Sw2*CorFin)/Deno
  !       Ae     =  2*T3e             /Deno
  !       Af     =  2*T3f             /Deno
* Angle dependent double-vector extra-correction
      VVCef  = ( (T3e) *(T3f) 
     $     -(QE*SWSQ+DeltSQ) *(T3f) *GSW(3) -QE*(T3f)*DeltV
     $     -(QF*SWSQ+DeltSQ) *(T3e) *GSW(2) -QF*(T3e)*DeltV
     $     + (QE*SWSQ) *(QF*SWSQ)  *GSW(4)
     $     + 2*QE*QF*DeltSQ*SWSQ + 2*QE*QF*DeltV*SWSQ )/Deno**2

      VVCor = 1D0
      IF(KEYGSW.NE.0.AND.KEYGSW.NE.4)  THEN
         VVCor  = VVCef/(Ve*Vf)
      ENDIF
C
C     PRINT *,' initwksw VVCor = ', VVCor
      NDIAG0=2
      NDIAGA=11
      KEYA  = 1
      KEYZ  = 1
C
C
      RETURN
      END
      FUNCTION T_BORNEW(MODE,KEYGSW,SVAR,COSTHE,TA,TB)
C ----------------------------------------------------------------------
C THIS ROUTINE PROVIDES BORN CROSS SECTION. IT HAS THE SAME         
C STRUCTURE AS FUNTIS AND FUNTIH, THUS CAN BE USED AS SIMPLER       
C EXAMPLE OF THE METHOD APPLIED THERE                               
C INPUT PARAMETERS ARE: SVAR    -- transfer
C                       COSTHE  -- cosine of angle between tau+ and 1st beam
C                       TA,TB   -- helicity states of tau+ tau-
C                       mode    -- parameter for mass terms; 1 means mass terms are on.
C                       keyGSW  -- keyGSW=0 gamma propagator is off 
C                                  keyGSW=10 running Z width 
C
C     called by : BORNY, BORAS, BORNV, WAGA, WEIGHT
C ----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / T_BEAMPM / ENE ,AMIN,AMFIN,IDE,IDF
      REAL*8              ENE ,AMIN,AMFIN
      COMMON / T_GAUSPM /SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI0   ,XUPZI0   ,XUPGF0   ,XUPZF0  
     &                  ,NDIAG0,NDIAGA,KEYA,KEYZ
     &                  ,ITCE,JTCE,ITCF,JTCF,KOLOR
      REAL*8            SS,POLN,T3E,QE,T3F,QF
     &                  ,XUPGI0(2),XUPZI0(2),XUPGF0(2),XUPZF0(2)
      COMMON / T_GAUSPM1/VVcor, ZetVPi, GamVPi
     &                  ,XUPGI   ,XUPZI   ,XUPGF   ,XUPZF
      COMPLEX*16         VVcor, ZetVPi, GamVPi
      COMPLEX*16         XUPGI(2),XUPZI(2),XUPGF(2),XUPZF(2)
      COMMON / T_EWn    / GMUn, ALPHAINVn
      REAL*8              GMUn, ALPHAINVn

      
      REAL*8            SEPS1,SEPS2
C=====================================================================
      COMMON / T_GSWPRMn /SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
      REAL*8             SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
C     SWSQ        = sin2 (theta Weinberg)
C     AMW,AMZ     = W & Z boson masses respectively
C     AMH         = the Higgs mass
C     AMTOP       = the top mass
C     GAMMZ       = Z0 width
      COMPLEX*16 ABORN(2,2),APHOT(2,2),AZETT(2,2)
      COMPLEX*16 XUPZFP(2),XUPZIP(2),XUPZIF(2,2)
      COMPLEX*16 ABORNM(2,2),APHOTM(2,2),AZETTM(2,2)
      COMPLEX*16 PROPA,PROPZ
      COMPLEX*16 XR,XI
      COMPLEX*16 XUPF,XUPI
      COMPLEX*16 XTHING
      DATA XI/(0.D0,1.D0)/,XR/(1.D0,0.D0)/
      DATA MODE0 /-5/
      DATA IDE0 /-55/
      DATA SVAR0,COST0 /-5.D0,-6.D0/
      DATA PI /3.141592653589793238462643D0/
      DATA SEPS1,SEPS2 /0D0,0D0/
 
C
C MEMORIZATION =========================================================
      IF ( MODE.NE.MODE0.OR.SVAR.NE.SVAR0.OR.COSTHE.NE.COST0
     $    .OR.IDE0.NE.IDE)THEN
C

   !      PRINT *,' T_BORN EW loop ( ',sqrt(svar),XUPGI(1),')= ', VVcor, ZetVPi!, GamVPi
   !      PRINT *,' T_BORN new( ',mode,')= ',SWSQ,AMW,AMZ,AMH,AMTOP,GAMMZ
C ** SWITCH OF MEMORISATION
C        IDE0=IDE
C        MODE0=MODE
C        SVAR0=SVAR
C        COST0=COSTHE
C ** PROPAGATORS
        SINTHE=SQRT(1.D0-COSTHE**2)
        BETA=SQRT(MAX(0D0,1D0-4D0*AMFIN**2/SVAR))
!        BETA=1.D0! Dec 10, 2019 mass term may need to be killed for EW tests 
C I MULTIPLY AXIAL COUPLING BY BETA FACTOR.
        XUPZFP(1)=0.5D0*(XUPZF(1)+XUPZF(2))+0.5D0*BETA*(XUPZF(1)-XUPZF(2))
        XUPZFP(2)=0.5D0*(XUPZF(1)+XUPZF(2))-0.5D0*BETA*(XUPZF(1)-XUPZF(2))
        XUPZIP(1)=0.5D0*(XUPZI(1)+XUPZI(2))+0.5D0*(XUPZI(1)-XUPZI(2))
        XUPZIP(2)=0.5D0*(XUPZI(1)+XUPZI(2))-0.5D0*(XUPZI(1)-XUPZI(2))
        XUPZIF(1,1)=(0.5D0*(XUPZI(1)+XUPZI(2))+0.5D0*(XUPZI(1)-XUPZI(2)))*(0.5D0*(XUPZF(1)+XUPZF(2))+0.5D0*BETA*(XUPZF(1)-XUPZF(2)))
     $             +(0.5D0*(XUPZI(1)+XUPZI(2)))*(0.5D0*(XUPZF(1)+XUPZF(2)))*(VVcor-1)
        XUPZIF(1,2)=(0.5D0*(XUPZI(1)+XUPZI(2))+0.5D0*(XUPZI(1)-XUPZI(2)))*(0.5D0*(XUPZF(1)+XUPZF(2))-0.5D0*BETA*(XUPZF(1)-XUPZF(2)))
     $             +(0.5D0*(XUPZI(1)+XUPZI(2)))*(0.5D0*(XUPZF(1)+XUPZF(2)))*(VVcor-1)
        XUPZIF(2,1)=(0.5D0*(XUPZI(1)+XUPZI(2))-0.5D0*(XUPZI(1)-XUPZI(2)))*(0.5D0*(XUPZF(1)+XUPZF(2))+0.5D0*BETA*(XUPZF(1)-XUPZF(2)))
     $             +(0.5D0*(XUPZI(1)+XUPZI(2)))*(0.5D0*(XUPZF(1)+XUPZF(2)))*(VVcor-1)
        XUPZIF(2,2)=(0.5D0*(XUPZI(1)+XUPZI(2))-0.5D0*(XUPZI(1)-XUPZI(2)))*(0.5D0*(XUPZF(1)+XUPZF(2))-0.5D0*BETA*(XUPZF(1)-XUPZF(2)))
     $             +(0.5D0*(XUPZI(1)+XUPZI(2)))*(0.5D0*(XUPZF(1)+XUPZF(2)))*(VVcor-1)
        
C FINAL STATE VECTOR COUPLING
        XUPF     =0.5D0*(XUPZF(1)+XUPZF(2))
        XUPI     =0.5D0*(XUPZI(1)+XUPZI(2))
        XTHING   =0D0
        

        PROPA =1D0/SVAR*GamVPi
C       use running width
        PROPZ =1D0/DCMPLX(SVAR-AMZ**2,SVAR/AMZ*GAMMZ)*ZetVPi


        IF( KEYGSW. EQ. 2) THEN
           Gfermi=GMUn
           ALPHAINV=ALPHAINVn
            ZetV =  Gfermi *AMZ**2 *ALPHAINV /(DSQRT(2.d0)*8.d0*PI)
     $          *(SWSQ*(1d0-SWSQ)) *16d0
            
!     variants of the Z propagators for the non-ew case
! ==1==            
 !        PROPZ =1D0/DCMPLX(SVAR-AMZ**2,AMZ*GAMMZ)*ZetV    !default 
! ==2==            
!         PROPZ =1D0/DCMPLX(SVAR-AMZ**2/(1+GAMMZ**2/AMZ**2),  ! alternative as
!     $                     AMZ*GAMMZ  /(1+GAMMZ**2/AMZ**2) ) ! running width
!     $                    *ZetV
!         PROPZ =PROPZ*DCMPLX(1,-GAMMZ/AMZ/(1+GAMMZ**2/AMZ**2))
! ==3==
         PROPZ =1D0/DCMPLX(SVAR-AMZ**2 ,                     ! running
     $                     GAMMZ*SVAR/AMZ   )*ZetV         
      ENDIF

C       use fixed width
      IF( KEYGSW. EQ. 10) THEN
!         PROPZ =1D0/DCMPLX(SVAR-AMZ**2,AMZ*GAMMZ)*ZetVPi ! this form need redefined M_Z and G_Z
!        below variant with this rescaling implemented
         PROPZ =1D0/DCMPLX(SVAR-AMZ**2/(1+GAMMZ**2/AMZ**2), ! alternative as
     $                     AMZ*GAMMZ  /(1+GAMMZ**2/AMZ**2) ) ! running width
     $                    *ZetV
         PROPZ =PROPZ*DCMPLX(1,-GAMMZ/AMZ/(1+GAMMZ**2/AMZ**2))

      ENDIF
        IF (KEYGSW.EQ.0) PROPA=0.D0
        DO 50 I=1,2
         DO 50 J=1,2
          REGULA= (3-2*I)*(3-2*J) + COSTHE
          REGULM=-(3-2*I)*(3-2*J) * SINTHE *2.D0*AMFIN/SQRT(SVAR)
          APHOT(I,J)=PROPA*(XUPGI(I)*XUPGF(J)*REGULA)
          AZETT(I,J)=PROPZ*(XUPZIP(I)*XUPZFP(J)+XTHING)*REGULA
          AZETT(I,J)=PROPZ*(XUPZIF(I,J)+XTHING)*REGULA         ! with electroweak effects in.          
          ABORN(I,J)=APHOT(I,J)+AZETT(I,J)
          APHOTM(I,J)=PROPA*DCMPLX(0D0,1D0)*XUPGI(I)*XUPGF(J)*REGULM
          AZETTM(I,J)=PROPZ*DCMPLX(0D0,1D0)*(XUPZIP(I)*XUPF+XTHING)*REGULM
          ABORNM(I,J)=APHOTM(I,J)+AZETTM(I,J)
   50   CONTINUE
      ENDIF
C
C******************
C* IN CALCULATING CROSS SECTION ONLY DIAGONAL ELEMENTS
C* OF THE SPIN DENSITY MATRICES ENTER (LONGITUD. POL. ONLY.)
C* HELICITY CONSERVATION EXPLICITLY OBEYED
      POLAR1=  (SEPS1)
      POLAR2= (-SEPS2)
      BORN=0D0
      DO 150 I=1,2
       HELIC= 3-2*I
       DO 150 J=1,2
        HELIT=3-2*J
        FACTOR=KOLOR*(1D0+HELIC*POLAR1)*(1D0-HELIC*POLAR2)/4D0
        FACTOM=FACTOR*(1+HELIT*TA)*(1-HELIT*TB)
        FACTOR=FACTOR*(1+HELIT*TA)*(1+HELIT*TB)

        BORN=BORN+CDABS(ABORN(I,J))**2*FACTOR
C      MASS TERM IN BORN
        IF (MODE.GE.1) THEN
         BORN=BORN+CDABS(ABORNM(I,J))**2*FACTOM
        ENDIF

  150 CONTINUE
C************
      FUNT=BORN
      IF(FUNT.LT.0.D0)  FUNT=BORN

C
      IF (SVAR.GT.4D0*AMFIN**2) THEN
C PHASE SPACE THRESHOLD FACTOR
        THRESH=SQRT(1-4D0*AMFIN**2/SVAR)
        T_BORNEW= FUNT*SVAR**2*THRESH
      ELSE
        THRESH=0.D0
        T_BORNEW=0.D0
      ENDIF
      END
