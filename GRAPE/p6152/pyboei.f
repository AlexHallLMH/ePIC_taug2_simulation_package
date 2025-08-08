 
C*********************************************************************
 
C...PYBOEI
C...Modifies an event so as to approximately take into account
C...Bose-Einstein effects according to a simple phenomenological
C...parametrization.
 
      SUBROUTINE PYBOEI(NSAV)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KEXCIT=4000000)
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
C...Local arrays and data.
      DIMENSION DPS(4),KFBE(9),NBE(0:10),BEI(100),BEI3(100),
     &BEIW(100),BEI3W(100)
      DATA KFBE/211,-211,111,321,-321,130,310,221,331/
C...Statement function: squared invariant mass.
      SDIP(I,J)=((P(I,4)+P(J,4))**2-(P(I,3)+P(J,3))**2-
     &(P(I,2)+P(J,2))**2-(P(I,1)+P(J,1))**2)
 
C...Boost event to overall CM frame. Calculate CM energy.
      IF((MSTJ(51).NE.1.AND.MSTJ(51).NE.2).OR.N-NSAV.LE.1) RETURN
      DO 100 J=1,4
        DPS(J)=0D0
  100 CONTINUE
      DO 120 I=1,N
        KFA=IABS(K(I,2))
        IF(K(I,1).LE.10.AND.((KFA.GT.10.AND.KFA.LE.20).OR.KFA.EQ.22)
     &  .AND.K(I,3).GT.0) THEN
          KFMA=IABS(K(K(I,3),2))
          IF(KFMA.GT.10.AND.KFMA.LE.80) K(I,1)=-K(I,1)
        ENDIF
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 120
        DO 110 J=1,4
          DPS(J)=DPS(J)+P(I,J)
  110   CONTINUE
  120 CONTINUE
      CALL PYROBO(0,0,0D0,0D0,-DPS(1)/DPS(4),-DPS(2)/DPS(4),
     &-DPS(3)/DPS(4))
      PECM=0D0
      DO 130 I=1,N
        IF(K(I,1).GE.1.AND.K(I,1).LE.10) PECM=PECM+P(I,4)
  130 CONTINUE
 
C...Reserve copy of particles by species at end of record.
      IWP=0
      IWN=0
      NBE(0)=N+MSTU(3)
      NMAX=NBE(0)
      SMMIN=PECM
      DO 180 IBE=1,MIN(10,MSTJ(52)+1)
        NBE(IBE)=NBE(IBE-1)
        DO 170 I=NSAV+1,N
          IF(IBE.EQ.MIN(10,MSTJ(52)+1)) THEN
            DO 140 IIBE=1,IBE-1
              IF(K(I,2).EQ.KFBE(IIBE)) GOTO 170
  140       CONTINUE
          ELSE
            IF(K(I,2).NE.KFBE(IBE)) GOTO 170
          ENDIF
          IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 170
          IF(NBE(IBE).GE.MSTU(4)-MSTU(32)-5) THEN
            CALL PYERRM(11,'(PYBOEI:) no more memory left in PYJETS')
            RETURN
          ENDIF
          NBE(IBE)=NBE(IBE)+1
          NMAX=NBE(IBE)
          K(NBE(IBE),1)=I
          K(NBE(IBE),5)=0
          SMMIN=MIN(SMMIN,P(I,5))
          IF(MSTJ(53).NE.0.OR.MSTJ(56).GT.0) THEN
            IM=I
  150       IF(K(IM,3).GT.0) THEN
              IM=K(IM,3)
              IF(ABS(K(IM,2)).NE.24) GOTO 150
              K(NBE(IBE),5)=K(IM,2)
              IF(IWP.EQ.0.AND.K(IM,2).EQ.24) IWP=IM
              IF(IWN.EQ.0.AND.K(IM,2).EQ.-24) IWN=IM
            ENDIF
          ENDIF
          DO 160 J=1,3
            P(NBE(IBE),J)=0D0
            V(NBE(IBE),J)=0D0
  160     CONTINUE
          P(NBE(IBE),5)=-1.0D0
  170   CONTINUE
  180 CONTINUE
      IF(NBE(MIN(9,MSTJ(52)))-NBE(0).LE.1) GOTO 500
 
C...Calculate separation between W+ and W-
      SIGW=PARJ(93)
      IF(IWP.GT.0.AND.IWN.GT.0.AND.MSTJ(56).GT.0) THEN
        DMW=PMAS(24,1)
        DGW=PMAS(24,2)
        DMP=P(IWP,5)
        DMN=P(IWN,5)
        TAUPD=DMP/SQRT((DMP**2-DMW**2)**2+(DGW*(DMP**2)/DMW)**2)
        TAUND=DMN/SQRT((DMN**2-DMW**2)**2+(DGW*(DMN**2)/DMW)**2)
        TAUP=-TAUPD*LOG(PYR(IDUM))
        TAUN=-TAUND*LOG(PYR(IDUM))
        DXP=TAUP*PYP(IWP,8)/DMP
        DXN=TAUN*PYP(IWN,8)/DMN
        DX=DXP+DXN
        SIGW=1.0D0/(1.0D0/PARJ(93)+REAL(MSTJ(56))*DX)
      ELSE
        SIGW=PARJ(93)
      ENDIF
 
      IF(MSTJ(57).EQ.1.AND.MSTJ(54).LT.0) THEN
        DO 210 IBE=1,MIN(9,MSTJ(52))
          DO 200 I1M=NBE(IBE-1)+1,NBE(IBE)-1
            Q2MIN=PECM**2
            I1=K(I1M,1)
            DO 190 I2M=NBE(IBE-1)+1,NBE(IBE)-1
              IF(I2M.EQ.I1M) GOTO 190
              I2=K(I2M,1)
              Q2=(P(I1,4)+P(I2,4))**2-(P(I1,1)+P(I2,1))**2-
     &        (P(I1,2)+P(I2,2))**2-(P(I1,3)+P(I2,3))**2-
     &        (P(I1,5)+P(I2,5))**2
              IF(Q2.GT.0.0D0.AND.Q2.LT.Q2MIN) THEN
                Q2MIN=Q2
              ENDIF
  190       CONTINUE
            P(I1M,5)=Q2MIN
  200     CONTINUE
  210   CONTINUE
      ENDIF
 
C...Tabulate integral for subsequent momentum shift.
      DO 390 IBE=1,MIN(9,MSTJ(52))
        IF(IBE.NE.1.AND.IBE.NE.4.AND.IBE.LE.7) GOTO 260
        IF(IBE.EQ.1.AND.MAX(NBE(1)-NBE(0),NBE(2)-NBE(1),NBE(3)-NBE(2))
     &  .LE.1) GOTO 260
        IF(IBE.EQ.4.AND.MAX(NBE(4)-NBE(3),NBE(5)-NBE(4),NBE(6)-NBE(5),
     &  NBE(7)-NBE(6)).LE.1) GOTO 260
        IF(IBE.GE.8.AND.NBE(IBE)-NBE(IBE-1).LE.1) GOTO 260
        IF(IBE.EQ.1) PMHQ=2D0*PYMASS(211)
        IF(IBE.EQ.4) PMHQ=2D0*PYMASS(321)
        IF(IBE.EQ.8) PMHQ=2D0*PYMASS(221)
        IF(IBE.EQ.9) PMHQ=2D0*PYMASS(331)
        QDEL=0.1D0*MIN(PMHQ,PARJ(93))
        QDEL3=0.1D0*MIN(PMHQ,PARJ(93)*3.0D0)
        QDELW=0.1D0*MIN(PMHQ,SIGW)
        QDEL3W=0.1D0*MIN(PMHQ,SIGW*3.0D0)
        IF(MSTJ(51).EQ.1) THEN
          NBIN=MIN(100,NINT(9D0*PARJ(93)/QDEL))
          NBIN3=MIN(100,NINT(27D0*PARJ(93)/QDEL3))
          NBINW=MIN(100,NINT(9D0*SIGW/QDELW))
          NBIN3W=MIN(100,NINT(27D0*SIGW/QDEL3W))
          BEEX=EXP(0.5D0*QDEL/PARJ(93))
          BEEX3=EXP(0.5D0*QDEL3/(3.0D0*PARJ(93)))
          BEEXW=EXP(0.5D0*QDELW/SIGW)
          BEEX3W=EXP(0.5D0*QDEL3W/(3.0D0*SIGW))
          BERT=EXP(-QDEL/PARJ(93))
          BERT3=EXP(-QDEL3/(3.0D0*PARJ(93)))
          BERTW=EXP(-QDELW/SIGW)
          BERT3W=EXP(-QDEL3W/(3.0D0*SIGW))
        ELSE
          NBIN=MIN(100,NINT(3D0*PARJ(93)/QDEL))
          NBIN3=MIN(100,NINT(9D0*PARJ(93)/QDEL3))
          NBINW=MIN(100,NINT(3D0*SIGW/QDELW))
          NBIN3W=MIN(100,NINT(9D0*SIGW/QDEL3W))
        ENDIF
        DO 220 IBIN=1,NBIN
          QBIN=QDEL*(IBIN-0.5D0)
          BEI(IBIN)=QDEL*(QBIN**2+QDEL**2/12D0)/SQRT(QBIN**2+PMHQ**2)
          IF(MSTJ(51).EQ.1) THEN
            BEEX=BEEX*BERT
            BEI(IBIN)=BEI(IBIN)*BEEX
          ELSE
            BEI(IBIN)=BEI(IBIN)*EXP(-(QBIN/PARJ(93))**2)
          ENDIF
          IF(IBIN.GE.2) BEI(IBIN)=BEI(IBIN)+BEI(IBIN-1)
  220   CONTINUE
        DO 230 IBIN=1,NBIN3
          QBIN=QDEL3*(IBIN-0.5D0)
          BEI3(IBIN)=QDEL3*(QBIN**2+QDEL3**2/12D0)/SQRT(QBIN**2+PMHQ**2)
          IF(MSTJ(51).EQ.1) THEN
            BEEX3=BEEX3*BERT3
            BEI3(IBIN)=BEI3(IBIN)*BEEX3
          ELSE
            BEI3(IBIN)=BEI3(IBIN)*EXP(-(QBIN/(3.0D0*PARJ(93)))**2)
          ENDIF
          IF(IBIN.GE.2) BEI3(IBIN)=BEI3(IBIN)+BEI3(IBIN-1)
  230   CONTINUE
        DO 240 IBIN=1,NBINW
          QBIN=QDELW*(IBIN-0.5D0)
          BEIW(IBIN)=QDELW*(QBIN**2+QDELW**2/12D0)/SQRT(QBIN**2+PMHQ**2)
          IF(MSTJ(51).EQ.1) THEN
            BEEXW=BEEXW*BERTW
            BEIW(IBIN)=BEIW(IBIN)*BEEXW
          ELSE
            BEIW(IBIN)=BEIW(IBIN)*EXP(-(QBIN/SIGW)**2)
          ENDIF
          IF(IBIN.GE.2) BEIW(IBIN)=BEIW(IBIN)+BEIW(IBIN-1)
  240   CONTINUE
        DO 250 IBIN=1,NBIN3W
          QBIN=QDEL3W*(IBIN-0.5D0)
          BEI3W(IBIN)=QDEL3W*(QBIN**2+QDEL3W**2/12D0)/
     &    SQRT(QBIN**2+PMHQ**2)
          IF(MSTJ(51).EQ.1) THEN
            BEEX3W=BEEX3W*BERT3W
            BEI3W(IBIN)=BEI3W(IBIN)*BEEX3W
          ELSE
            BEI3W(IBIN)=BEI3W(IBIN)*EXP(-(QBIN/(3.0D0*SIGW))**2)
          ENDIF
          IF(IBIN.GE.2) BEI3W(IBIN)=BEI3W(IBIN)+BEI3W(IBIN-1)
  250   CONTINUE
 
C...Loop through particle pairs and find old relative momentum.
  260   DO 380 I1M=NBE(IBE-1)+1,NBE(IBE)-1
          I1=K(I1M,1)
          DO 370 I2M=I1M+1,NBE(IBE)
            IF(MSTJ(53).EQ.1.AND.K(I1M,5).NE.K(I2M,5)) GOTO 370
            IF(MSTJ(53).EQ.2.AND.K(I1M,5).EQ.K(I2M,5)) GOTO 370
            I2=K(I2M,1)
            Q2OLD=(P(I1,4)+P(I2,4))**2-(P(I1,1)+P(I2,1))**2-(P(I1,2)+
     &      P(I2,2))**2-(P(I1,3)+P(I2,3))**2-(P(I1,5)+P(I2,5))**2
            IF(Q2OLD.LE.0.0D0) GOTO 370
            QOLD=SQRT(Q2OLD)
 
C...Calculate new relative momentum.
            QMOV=0.0D0
            QMOV3=0.0D0
            QMOVW=0.0D0
            QMOV3W=0.0D0
            IF(QOLD.LT.1D-3*QDEL) THEN
              GOTO 270
            ELSEIF(QOLD.LE.QDEL) THEN
              QMOV=QOLD/3D0
            ELSEIF(QOLD.LT.(NBIN-0.1D0)*QDEL) THEN
              RBIN=QOLD/QDEL
              IBIN=RBIN
              RINP=(RBIN**3-IBIN**3)/(3*IBIN*(IBIN+1)+1)
              QMOV=(BEI(IBIN)+RINP*(BEI(IBIN+1)-BEI(IBIN)))*
     &        SQRT(Q2OLD+PMHQ**2)/Q2OLD
            ELSE
              QMOV=BEI(NBIN)*SQRT(Q2OLD+PMHQ**2)/Q2OLD
            ENDIF
  270       Q2NEW=Q2OLD*(QOLD/(QOLD+3D0*PARJ(92)*QMOV))**(2D0/3D0)
            IF(QOLD.LT.1D-3*QDEL3) THEN
              GOTO 280
            ELSEIF(QOLD.LE.QDEL3) THEN
              QMOV3=QOLD/3D0
            ELSEIF(QOLD.LT.(NBIN3-0.1D0)*QDEL3) THEN
              RBIN3=QOLD/QDEL3
              IBIN3=RBIN3
              RINP3=(RBIN3**3-IBIN3**3)/(3*IBIN3*(IBIN3+1)+1)
              QMOV3=(BEI3(IBIN3)+RINP3*(BEI3(IBIN3+1)-BEI3(IBIN3)))*
     &        SQRT(Q2OLD+PMHQ**2)/Q2OLD
            ELSE
              QMOV3=BEI3(NBIN3)*SQRT(Q2OLD+PMHQ**2)/Q2OLD
            ENDIF
  280       Q2NEW3=Q2OLD*(QOLD/(QOLD+3D0*PARJ(92)*QMOV3))**(2D0/3D0)
            RSCALE=1.0D0
            IF(MSTJ(54).EQ.2)
     &      RSCALE=1.0D0-EXP(-(QOLD/(2D0*PARJ(93)))**2)
            IF(MSTJ(56).LE.0.OR.IWP.EQ.0.OR.IWN.EQ.0.OR.
     &      K(I1M,5).EQ.K(I2M,5)) GOTO 310
 
            IF(QOLD.LT.1D-3*QDELW) THEN
              GOTO 290
            ELSEIF(QOLD.LE.QDELW) THEN
              QMOVW=QOLD/3D0
            ELSEIF(QOLD.LT.(NBINW-0.1D0)*QDELW) THEN
              RBINW=QOLD/QDELW
              IBINW=RBINW
              RINPW=(RBINW**3-IBINW**3)/(3*IBINW*(IBINW+1)+1)
              QMOVW=(BEIW(IBINW)+RINPW*(BEIW(IBINW+1)-BEIW(IBINW)))*
     &        SQRT(Q2OLD+PMHQ**2)/Q2OLD
            ELSE
              QMOVW=BEIW(NBINW)*SQRT(Q2OLD+PMHQ**2)/Q2OLD
            ENDIF
  290       Q2NEW=Q2OLD*(QOLD/(QOLD+3D0*PARJ(92)*QMOVW))**(2D0/3D0)
            IF(QOLD.LT.1D-3*QDEL3W) THEN
              GOTO 300
            ELSEIF(QOLD.LE.QDEL3W) THEN
              QMOV3W=QOLD/3D0
            ELSEIF(QOLD.LT.(NBIN3W-0.1D0)*QDEL3W) THEN
              RBIN3W=QOLD/QDEL3W
              IBIN3W=RBIN3W
              RINP3W=(RBIN3W**3-IBIN3W**3)/(3*IBIN3W*(IBIN3W+1)+1)
              QMOV3W=(BEI3W(IBIN3W)+RINP3W*(BEI3W(IBIN3W+1)-
     &        BEI3W(IBIN3W)))*SQRT(Q2OLD+PMHQ**2)/Q2OLD
            ELSE
              QMOV3W=BEI3W(NBIN3W)*SQRT(Q2OLD+PMHQ**2)/Q2OLD
            ENDIF
  300       Q2NEW3=Q2OLD*(QOLD/(QOLD+3D0*PARJ(92)*QMOV3W))**(2D0/3D0)
            IF(MSTJ(54).EQ.2)
     &      RSCALE=1.0D0-EXP(-(QOLD/(2D0*SIGW))**2)
 
  310       CALL PYBESQ(I1,I2,NMAX,Q2OLD,Q2NEW)
            DO 320 J=1,3
              P(I1M,J)=P(I1M,J)+P(NMAX+1,J)
              P(I2M,J)=P(I2M,J)+P(NMAX+2,J)
  320       CONTINUE
            IF(MSTJ(54).GE.1) THEN
              CALL PYBESQ(I1,I2,NMAX,Q2OLD,Q2NEW3)
              DO 330 J=1,3
                V(I1M,J)=V(I1M,J)+P(NMAX+1,J)*RSCALE
                V(I2M,J)=V(I2M,J)+P(NMAX+2,J)*RSCALE
  330         CONTINUE
            ELSEIF(MSTJ(54).LE.-1) THEN
              EDEL=P(I1,4)+P(I2,4)-
     &        SQRT(MAX(Q2NEW-Q2OLD+(P(I1,4)+P(I2,4))**2,0.0D0))
              A2=(P(I1,1)-P(I2,1))**2+(P(I1,2)-P(I2,2))**2+
     &        (P(I1,3)-P(I2,3))**2
              WMAX=-1.0D20
              MI3=0
              MI4=0
              S12=SDIP(I1,I2)
              SM1=(P(I1,5)+SMMIN)**2
              DO 350 I3M=NBE(0)+1,NBE(MIN(10,MSTJ(52)+1))
                IF(I3M.EQ.I1M.OR.I3M.EQ.I2M) GOTO 350
                IF(MSTJ(53).EQ.1.AND.K(I3M,5).NE.K(I1M,5)) GOTO 350
                IF(MSTJ(53).EQ.-2.AND.K(I1M,5).EQ.K(I2M,5).AND.
     &          K(I3M,5).NE.K(I1M,5)) GOTO 350
                I3=K(I3M,1)
                IF(K(I3,2).EQ.K(I1,2)) GOTO 350
                S13=SDIP(I1,I3)
                S23=SDIP(I2,I3)
                SM3=(P(I3,5)+SMMIN)**2
                IF(MSTJ(54).EQ.-2) THEN
                  WI=(MIN(S12*SM3,S13*MIN(SM1,SM3),
     &            S23*MIN(SM1,SM3))*SM1)
                ELSE
                  WI=((P(I1,4)+P(I2,4)+P(I3,4))**2-
     &            (P(I1,3)+P(I2,3)+P(I3,3))**2-
     &            (P(I1,2)+P(I2,2)+P(I3,2))**2-
     &            (P(I1,1)+P(I2,1)+P(I3,1))**2)
                ENDIF
                IF(MSTJ(57).EQ.1.AND.P(I3M,5).GT.0) THEN
                  IF (WMAX*WI.GE.(1.0D0-EXP(-P(I3M,5)/(PARJ(93)**2))))
     &                 GOTO 350
                ELSE
                  IF(WMAX*WI.GE.1.0) GOTO 350
                ENDIF
                DO 340 I4M=I3M+1,NBE(MIN(10,MSTJ(52)+1))
                  IF(I4M.EQ.I1M.OR.I4M.EQ.I2M) GOTO 340
                  IF(MSTJ(53).EQ.1.AND.K(I4M,5).NE.K(I1M,5)) GOTO 340
                  IF(MSTJ(53).EQ.-2.AND.K(I1M,5).EQ.K(I2M,5).AND.
     &            K(I4M,5).NE.K(I1M,5)) GOTO 340
                  I4=K(I4M,1)
                  IF(K(I3,2).EQ.K(I4,2).OR.K(I4,2).EQ.K(I1,2))
     &            GOTO 340
                  IF((P(I3,4)+P(I4,4)+EDEL)**2.LT.
     &            (P(I3,1)+P(I4,1))**2+(P(I3,2)+P(I4,2))**2+
     &            (P(I3,3)+P(I4,3))**2+(P(I3,5)+P(I4,5))**2)
     &            GOTO 340
                  IF(MSTJ(54).EQ.-2) THEN
                    S14=SDIP(I1,I4)
                    S24=SDIP(I2,I4)
                    S34=SDIP(I3,I4)
                    W=S12*MIN(MIN(S23,S24),MIN(S13,S14))*S34
                    W=MIN(W,S13*MIN(MIN(S23,S34),S12)*S24)
                    W=MIN(W,S14*MIN(MIN(S24,S34),S12)*S23)
                    W=MIN(W,MIN(S23,S24)*S13*S14)
                    W=1.0D0/W
                  ELSE
C...weight=1-cos(theta)/mtot2
                    S1234=(P(I1,4)+P(I2,4)+P(I3,4)+P(I4,4))**2-
     &              (P(I1,3)+P(I2,3)+P(I3,3)+P(I4,3))**2-
     &              (P(I1,2)+P(I2,2)+P(I3,2)+P(I4,2))**2-
     &              (P(I1,1)+P(I2,1)+P(I3,1)+P(I4,1))**2
                    W=1.0D0/S1234
                    IF(W.LE.WMAX) GOTO 340
                  ENDIF
                  IF(MSTJ(57).EQ.1.AND.P(I3M,5).GT.0)
     &            W=W*(1.0D0-EXP(-P(I3M,5)/(PARJ(93)**2)))
                  IF(MSTJ(57).EQ.1.AND.P(I4M,5).GT.0)
     &            W=W*(1.0D0-EXP(-P(I4M,5)/(PARJ(93)**2)))
                  IF(W.LE.WMAX) GOTO 340
                  MI3=I3M
                  MI4=I4M
                  WMAX=W
  340           CONTINUE
  350         CONTINUE
              IF(MI4.EQ.0) GOTO 370
              I3=K(MI3,1)
              I4=K(MI4,1)
              EOLD=P(I3,4)+P(I4,4)
              ENEW=EOLD+EDEL
              P2=(P(I3,1)+P(I4,1))**2+(P(I3,2)+P(I4,2))**2+
     &        (P(I3,3)+P(I4,3))**2
              Q2NEWP=MAX(0.0D0,ENEW**2-P2-(P(I3,5)+P(I4,5))**2)
              Q2OLDP=MAX(0.0D0,EOLD**2-P2-(P(I3,5)+P(I4,5))**2)
              CALL PYBESQ(I3,I4,NMAX,Q2OLDP,Q2NEWP)
              DO 360 J=1,3
                V(MI3,J)=V(MI3,J)+P(NMAX+1,J)
                V(MI4,J)=V(MI4,J)+P(NMAX+2,J)
  360         CONTINUE
            ENDIF
  370     CONTINUE
  380   CONTINUE
  390 CONTINUE
 
C...Shift momenta and recalculate energies.
      ESUMP=0.0D0
      ESUM=0.0D0
      PROD=0.0D0
      DO 420 IM=NBE(0)+1,NBE(MIN(10,MSTJ(52)+1))
        I=K(IM,1)
        ESUMP=ESUMP+P(I,4)
        DO 400 J=1,3
          P(I,J)=P(I,J)+P(IM,J)
  400   CONTINUE
        P(I,4)=SQRT(P(I,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2)
        ESUM=ESUM+P(I,4)
        DO 410 J=1,3
          PROD=PROD+V(IM,J)*P(I,J)/P(I,4)
  410   CONTINUE
  420 CONTINUE
 
      PARJ(96)=0.0D0
      IF(MSTJ(54).NE.0.AND.PROD.NE.0.0D0) THEN
  430   ALPHA=(ESUMP-ESUM)/PROD
        PARJ(96)=PARJ(96)+ALPHA
        PROD=0.0D0
        ESUM=0.0D0
        DO 460 IM=NBE(0)+1,NBE(MIN(10,MSTJ(52)+1))
          I=K(IM,1)
          DO 440 J=1,3
            P(I,J)=P(I,J)+ALPHA*V(IM,J)
  440     CONTINUE
          P(I,4)=SQRT(P(I,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2)
          ESUM=ESUM+P(I,4)
          DO 450 J=1,3
            PROD=PROD+V(IM,J)*P(I,J)/P(I,4)
  450     CONTINUE
  460   CONTINUE
        IF(PROD.NE.0.0D0.AND.ABS(ESUMP-ESUM)/PECM.GT.0.00001D0)
     &  GOTO 430
      ENDIF
 
C...Rescale all momenta for energy conservation.
      PES=0D0
      PQS=0D0
      DO 470 I=1,N
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 470
        PES=PES+P(I,4)
        PQS=PQS+P(I,5)**2/P(I,4)
  470 CONTINUE
      PARJ(95)=PES-PECM
      FAC=(PECM-PQS)/(PES-PQS)
      DO 490 I=1,N
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 490
        DO 480 J=1,3
          P(I,J)=FAC*P(I,J)
  480   CONTINUE
        P(I,4)=SQRT(P(I,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2)
  490 CONTINUE
 
C...Boost back to correct reference frame.
  500 CALL PYROBO(0,0,0D0,0D0,DPS(1)/DPS(4),DPS(2)/DPS(4),DPS(3)/DPS(4))
      DO 510 I=1,N
        IF(K(I,1).LT.0) K(I,1)=-K(I,1)
  510 CONTINUE
 
      RETURN
      END
