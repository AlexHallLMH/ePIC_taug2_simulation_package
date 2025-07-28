 
C*********************************************************************
 
C...PYSHOW
C...Generates timelike parton showers from given partons.
 
      SUBROUTINE PYSHOW(IP1,IP2,QMAX)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
C...Local arrays.
      DIMENSION PMTH(5,50),PS(5),PMA(4),PMSD(4),IEP(4),IPA(4),
     &KFLA(4),KFLD(4),KFL(4),ITRY(4),ISI(4),ISL(4),DP(4),DPT(5,4),
     &KSH(0:40),KCII(2),NIIS(2),IIIS(2,2),THEIIS(2,2),PHIIIS(2,2),
     &ISII(2),ISSET(3)
 
C...Check that QMAX not too low.
      IF(MSTJ(41).LE.0) THEN
        RETURN
      ELSEIF(MSTJ(41).EQ.1) THEN
        IF(QMAX.LE.PARJ(82).AND.IP2.GT.-5) RETURN
      ELSE
        IF(QMAX.LE.MIN(PARJ(82),PARJ(83),PARJ(90)).AND.IP2.GT.-5) 
     &  RETURN
      ENDIF
 
C...Initialization of cutoff masses etc.
      DO 100 IFL=0,40
        KSH(IFL)=0
  100 CONTINUE
      KSH(21)=1
      PMTH(1,21)=PYMASS(21)
      PMTH(2,21)=SQRT(PMTH(1,21)**2+0.25D0*PARJ(82)**2)
      PMTH(3,21)=2D0*PMTH(2,21)
      PMTH(4,21)=PMTH(3,21)
      PMTH(5,21)=PMTH(3,21)
      PMTH(1,22)=PYMASS(22)
      PMTH(2,22)=SQRT(PMTH(1,22)**2+0.25D0*PARJ(83)**2)
      PMTH(3,22)=2D0*PMTH(2,22)
      PMTH(4,22)=PMTH(3,22)
      PMTH(5,22)=PMTH(3,22)
      PMQTH1=PARJ(82)
      IF(MSTJ(41).GE.2) PMQTH1=MIN(PARJ(82),PARJ(83))
      PMQT1E=MIN(PMQTH1,PARJ(90))
      PMQTH2=PMTH(2,21)
      IF(MSTJ(41).GE.2) PMQTH2=MIN(PMTH(2,21),PMTH(2,22))
      PMQT2E=MIN(PMQTH2,0.5D0*PARJ(90))
      DO 110 IFL=1,8
        KSH(IFL)=1
        PMTH(1,IFL)=PYMASS(IFL)
        PMTH(2,IFL)=SQRT(PMTH(1,IFL)**2+0.25D0*PMQTH1**2)
        PMTH(3,IFL)=PMTH(2,IFL)+PMQTH2
        PMTH(4,IFL)=SQRT(PMTH(1,IFL)**2+0.25D0*PARJ(82)**2)+PMTH(2,21)
        PMTH(5,IFL)=SQRT(PMTH(1,IFL)**2+0.25D0*PARJ(83)**2)+PMTH(2,22)
  110 CONTINUE
      DO 120 IFL=11,17,2
        IF(MSTJ(41).GE.2) KSH(IFL)=1
        PMTH(1,IFL)=PYMASS(IFL)
        PMTH(2,IFL)=SQRT(PMTH(1,IFL)**2+0.25D0*PARJ(90)**2)
        PMTH(3,IFL)=PMTH(2,IFL)+0.5D0*PARJ(90)
        PMTH(4,IFL)=PMTH(3,IFL)
        PMTH(5,IFL)=PMTH(3,IFL)
  120 CONTINUE
      PT2MIN=MAX(0.5D0*PARJ(82),1.1D0*PARJ(81))**2
      ALAMS=PARJ(81)**2
      ALFM=LOG(PT2MIN/ALAMS)
 
C...Store positions of shower initiating partons.
      MPSPD=0
      IF(IP1.GT.0.AND.IP1.LE.MIN(N,MSTU(4)-MSTU(32)).AND.IP2.EQ.0) THEN
        NPA=1
        IPA(1)=IP1
      ELSEIF(MIN(IP1,IP2).GT.0.AND.MAX(IP1,IP2).LE.MIN(N,MSTU(4)-
     &  MSTU(32))) THEN
        NPA=2
        IPA(1)=IP1
        IPA(2)=IP2
      ELSEIF(IP1.GT.0.AND.IP1.LE.MIN(N,MSTU(4)-MSTU(32)).AND.IP2.LT.0
     &  .AND.IP2.GE.-3) THEN
        NPA=IABS(IP2)
        DO 130 I=1,NPA
          IPA(I)=IP1+I-1
  130   CONTINUE
      ELSEIF(IP1.GT.0.AND.IP1.LE.MIN(N,MSTU(4)-MSTU(32)).AND.
     &IP2.EQ.-8) THEN
        MPSPD=1
        NPA=2
        IPA(1)=IP1+6
        IPA(2)=IP1+7  
      ELSE
        CALL PYERRM(12,
     &  '(PYSHOW:) failed to reconstruct showering system')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
 
C...Check on phase space available for emission.
      IREJ=0
      DO 140 J=1,5
        PS(J)=0D0
  140 CONTINUE
      PM=0D0
      DO 160 I=1,NPA
        KFLA(I)=IABS(K(IPA(I),2))
        PMA(I)=P(IPA(I),5)
C...Special cutoff masses for t, l, h with variable masses.
        IFLA=KFLA(I)
        IF(KFLA(I).GE.6.AND.KFLA(I).LE.8) THEN
          IFLA=37+KFLA(I)+ISIGN(2,K(IPA(I),2))
          PMTH(1,IFLA)=PMA(I)
          PMTH(2,IFLA)=SQRT(PMTH(1,IFLA)**2+0.25D0*PMQTH1**2)
          PMTH(3,IFLA)=PMTH(2,IFLA)+PMQTH2
          PMTH(4,IFLA)=SQRT(PMTH(1,IFLA)**2+0.25D0*PARJ(82)**2)+
     &    PMTH(2,21)
          PMTH(5,IFLA)=SQRT(PMTH(1,IFLA)**2+0.25D0*PARJ(83)**2)+
     &    PMTH(2,22)
        ENDIF
        IF(KFLA(I).LE.40) THEN
          IF(KSH(KFLA(I)).EQ.1) PMA(I)=PMTH(3,IFLA)
        ENDIF
        PM=PM+PMA(I)
        IF(KFLA(I).GT.40) THEN
          IREJ=IREJ+1
        ELSE
          IF(KSH(KFLA(I)).EQ.0.OR.PMA(I).GT.QMAX) IREJ=IREJ+1
        ENDIF
        DO 150 J=1,4
          PS(J)=PS(J)+P(IPA(I),J)
  150   CONTINUE
  160 CONTINUE
      IF(IREJ.EQ.NPA.AND.IP2.GT.-5) RETURN
      PS(5)=SQRT(MAX(0D0,PS(4)**2-PS(1)**2-PS(2)**2-PS(3)**2))
      IF(NPA.EQ.1) PS(5)=PS(4)
      IF(PS(5).LE.PM+PMQT1E) RETURN
 
C...Check if 3-jet matrix elements to be used.
      M3JC=0
      IF(NPA.EQ.2.AND.MSTJ(47).GE.1.AND.MPSPD.EQ.0) THEN
        IF(KFLA(1).GE.1.AND.KFLA(1).LE.8.AND.KFLA(2).GE.1.AND.
     &  KFLA(2).LE.8) M3JC=1
        IF((KFLA(1).EQ.11.OR.KFLA(1).EQ.13.OR.KFLA(1).EQ.15.OR.
     &  KFLA(1).EQ.17).AND.KFLA(2).EQ.KFLA(1)) M3JC=1
        IF((KFLA(1).EQ.11.OR.KFLA(1).EQ.13.OR.KFLA(1).EQ.15.OR.
     &  KFLA(1).EQ.17).AND.KFLA(2).EQ.KFLA(1)+1) M3JC=1
        IF((KFLA(1).EQ.12.OR.KFLA(1).EQ.14.OR.KFLA(1).EQ.16.OR.
     &  KFLA(1).EQ.18).AND.KFLA(2).EQ.KFLA(1)-1) M3JC=1
        IF(MSTJ(47).EQ.2.OR.MSTJ(47).EQ.4) M3JC=1
        M3JCM=0
        IF(M3JC.EQ.1.AND.MSTJ(47).GE.3.AND.KFLA(1).EQ.KFLA(2)) THEN
          M3JCM=1
          PQMES=PMTH(1,KFLA(1))**2
          QME=4D0*PQMES/PS(5)**2
          RESCZ=MIN(1D0,LOG(PMTH(2,KFLA(1))/PS(5))/
     &    LOG(PMTH(2,21)/PS(5)))
        ENDIF
      ENDIF
 
C...Find if interference with initial state partons.
      MIIS=0
      IF(MSTJ(50).GE.1.AND.MSTJ(50).LE.3.AND.NPA.EQ.2.AND.MPSPD.EQ.0) 
     &MIIS=MSTJ(50)
      IF(MIIS.NE.0) THEN
        DO 180 I=1,2
          KCII(I)=0
          KCA=PYCOMP(KFLA(I))
          IF(KCA.NE.0) KCII(I)=KCHG(KCA,2)*ISIGN(1,K(IPA(I),2))
          NIIS(I)=0
          IF(KCII(I).NE.0) THEN
            DO 170 J=1,2
              ICSI=MOD(K(IPA(I),3+J)/MSTU(5),MSTU(5))
              IF(ICSI.GT.0.AND.ICSI.NE.IPA(1).AND.ICSI.NE.IPA(2).AND.
     &        (KCII(I).EQ.(-1)**(J+1).OR.KCII(I).EQ.2)) THEN
                NIIS(I)=NIIS(I)+1
                IIIS(I,NIIS(I))=ICSI
              ENDIF
  170       CONTINUE
          ENDIF
  180   CONTINUE
        IF(NIIS(1)+NIIS(2).EQ.0) MIIS=0
      ENDIF
 
C...Boost interfering initial partons to rest frame
C...and reconstruct their polar and azimuthal angles.
      IF(MIIS.NE.0) THEN
        DO 200 I=1,2
          DO 190 J=1,5
            K(N+I,J)=K(IPA(I),J)
            P(N+I,J)=P(IPA(I),J)
            V(N+I,J)=0D0
  190     CONTINUE
  200   CONTINUE
        DO 220 I=3,2+NIIS(1)
          DO 210 J=1,5
            K(N+I,J)=K(IIIS(1,I-2),J)
            P(N+I,J)=P(IIIS(1,I-2),J)
            V(N+I,J)=0D0
  210     CONTINUE
  220   CONTINUE
        DO 240 I=3+NIIS(1),2+NIIS(1)+NIIS(2)
          DO 230 J=1,5
            K(N+I,J)=K(IIIS(2,I-2-NIIS(1)),J)
            P(N+I,J)=P(IIIS(2,I-2-NIIS(1)),J)
            V(N+I,J)=0D0
  230     CONTINUE
  240   CONTINUE
        CALL PYROBO(N+1,N+2+NIIS(1)+NIIS(2),0D0,0D0,-PS(1)/PS(4),
     &  -PS(2)/PS(4),-PS(3)/PS(4))
        PHI=PYANGL(P(N+1,1),P(N+1,2))
        CALL PYROBO(N+1,N+2+NIIS(1)+NIIS(2),0D0,-PHI,0D0,0D0,0D0)
        THE=PYANGL(P(N+1,3),P(N+1,1))
        CALL PYROBO(N+1,N+2+NIIS(1)+NIIS(2),-THE,0D0,0D0,0D0,0D0)
        DO 250 I=3,2+NIIS(1)
          THEIIS(1,I-2)=PYANGL(P(N+I,3),SQRT(P(N+I,1)**2+P(N+I,2)**2))
          PHIIIS(1,I-2)=PYANGL(P(N+I,1),P(N+I,2))
  250   CONTINUE
        DO 260 I=3+NIIS(1),2+NIIS(1)+NIIS(2)
          THEIIS(2,I-2-NIIS(1))=PARU(1)-PYANGL(P(N+I,3),
     &    SQRT(P(N+I,1)**2+P(N+I,2)**2))
          PHIIIS(2,I-2-NIIS(1))=PYANGL(P(N+I,1),P(N+I,2))
  260   CONTINUE
      ENDIF
 
C...Define imagined single initiator of shower for parton system.
      NS=N
      IF(N.GT.MSTU(4)-MSTU(32)-5) THEN
        CALL PYERRM(11,'(PYSHOW:) no more memory left in PYJETS')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
  265 N=NS
      IF(NPA.GE.2) THEN
        K(N+1,1)=11
        K(N+1,2)=21
        K(N+1,3)=0
        K(N+1,4)=0
        K(N+1,5)=0
        P(N+1,1)=0D0
        P(N+1,2)=0D0
        P(N+1,3)=0D0
        P(N+1,4)=PS(5)
        P(N+1,5)=PS(5)
        V(N+1,5)=PS(5)**2
        N=N+1
      ENDIF
 
C...Loop over partons that may branch.
      NEP=NPA
      IM=NS
      IF(NPA.EQ.1) IM=NS-1
  270 IM=IM+1
      IF(N.GT.NS) THEN
        IF(IM.GT.N) GOTO 510
        KFLM=IABS(K(IM,2))
        IF(KFLM.GT.40) GOTO 270
        IF(KSH(KFLM).EQ.0) GOTO 270
        IFLM=KFLM
        IF(KFLM.GE.6.AND.KFLM.LE.8) IFLM=37+KFLM+ISIGN(2,K(IM,2))
        IF(P(IM,5).LT.PMTH(2,IFLM)) GOTO 270
        IGM=K(IM,3)
      ELSE
        IGM=-1
      ENDIF
      IF(N+NEP.GT.MSTU(4)-MSTU(32)-5) THEN
        CALL PYERRM(11,'(PYSHOW:) no more memory left in PYJETS')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
 
C...Position of aunt (sister to branching parton).
C...Origin and flavour of daughters.
      IAU=0
      IF(IGM.GT.0) THEN
        IF(K(IM-1,3).EQ.IGM) IAU=IM-1
        IF(N.GE.IM+1.AND.K(IM+1,3).EQ.IGM) IAU=IM+1
      ENDIF
      IF(IGM.GE.0) THEN
        K(IM,4)=N+1
        DO 280 I=1,NEP
          K(N+I,3)=IM
  280   CONTINUE
      ELSE
        K(N+1,3)=IPA(1)
      ENDIF
      IF(IGM.LE.0) THEN
        DO 290 I=1,NEP
          K(N+I,2)=K(IPA(I),2)
  290   CONTINUE
      ELSEIF(KFLM.NE.21) THEN
        K(N+1,2)=K(IM,2)
        K(N+2,2)=K(IM,5)
      ELSEIF(K(IM,5).EQ.21) THEN
        K(N+1,2)=21
        K(N+2,2)=21
      ELSE
        K(N+1,2)=K(IM,5)
        K(N+2,2)=-K(IM,5)
      ENDIF
 
C...Reset flags on daughters and tries made.
      DO 300 IP=1,NEP
        K(N+IP,1)=3
        K(N+IP,4)=0
        K(N+IP,5)=0
        KFLD(IP)=IABS(K(N+IP,2))
        IF(KCHG(PYCOMP(KFLD(IP)),2).EQ.0) K(N+IP,1)=1
        ITRY(IP)=0
        ISL(IP)=0
        ISI(IP)=0
        IF(KFLD(IP).LE.40) THEN
          IF(KSH(KFLD(IP)).EQ.1) ISI(IP)=1
        ENDIF
  300 CONTINUE
      ISLM=0
 
C...Maximum virtuality of daughters.
      IF(IGM.LE.0) THEN
        DO 310 I=1,NPA
          IF(NPA.GE.3) P(N+I,4)=(PS(4)*P(IPA(I),4)-PS(1)*P(IPA(I),1)-
     &    PS(2)*P(IPA(I),2)-PS(3)*P(IPA(I),3))/PS(5)
          P(N+I,5)=MIN(QMAX,PS(5))
          IF(IP2.LE.-5) P(N+I,5)=MAX(P(N+I,5),
     &    2D0*PMTH(3,IABS(K(N+I,2))))
          IF(NPA.GE.3) P(N+I,5)=MIN(P(N+I,5),P(N+I,4))
          IF(ISI(I).EQ.0) P(N+I,5)=P(IPA(I),5)
  310   CONTINUE
      ELSE
        IF(MSTJ(43).LE.2) PEM=V(IM,2)
        IF(MSTJ(43).GE.3) PEM=P(IM,4)
        P(N+1,5)=MIN(P(IM,5),V(IM,1)*PEM)
        P(N+2,5)=MIN(P(IM,5),(1D0-V(IM,1))*PEM)
        IF(K(N+2,2).EQ.22) P(N+2,5)=PMTH(1,22)
      ENDIF
      DO 320 I=1,NEP
        PMSD(I)=P(N+I,5)
        IF(ISI(I).EQ.1) THEN
          IFLD=KFLD(I)
          IF(KFLD(I).GE.6.AND.KFLD(I).LE.8) IFLD=37+KFLD(I)+
     &    ISIGN(2,K(N+I,2))
          IF(P(N+I,5).LE.PMTH(3,IFLD)) P(N+I,5)=PMTH(1,IFLD)
        ENDIF
        V(N+I,5)=P(N+I,5)**2
  320 CONTINUE
 
C...Choose one of the daughters for evolution.
  330 INUM=0
      IF(NEP.EQ.1) INUM=1
      DO 340 I=1,NEP
        IF(INUM.EQ.0.AND.ISL(I).EQ.1) INUM=I
  340 CONTINUE
      DO 350 I=1,NEP
        IF(INUM.EQ.0.AND.ITRY(I).EQ.0.AND.ISI(I).EQ.1) THEN
          IFLD=KFLD(I)
          IF(KFLD(I).GE.6.AND.KFLD(I).LE.8) IFLD=37+KFLD(I)+
     &    ISIGN(2,K(N+I,2))
          IF(P(N+I,5).GE.PMTH(2,IFLD)) INUM=I
        ENDIF
  350 CONTINUE
      IF(INUM.EQ.0) THEN
        RMAX=0D0
        DO 360 I=1,NEP
          IF(ISI(I).EQ.1.AND.PMSD(I).GE.PMQT2E) THEN
            RPM=P(N+I,5)/PMSD(I)
            IFLD=KFLD(I)
            IF(KFLD(I).GE.6.AND.KFLD(I).LE.8) IFLD=37+KFLD(I)+
     &      ISIGN(2,K(N+I,2))
            IF(RPM.GT.RMAX.AND.P(N+I,5).GE.PMTH(2,IFLD)) THEN
              RMAX=RPM
              INUM=I
            ENDIF
          ENDIF
  360   CONTINUE
      ENDIF

C...Cancel choice of predetermined daughter already treated.
      INUM=MAX(1,INUM)
      INUMT=INUM 
      IF(MPSPD.EQ.1.AND.IGM.EQ.0.AND.ITRY(INUMT).GE.1) THEN
        IF(K(IP1-1+INUM,4).GT.0) INUM=3-INUM
      ELSEIF(MPSPD.EQ.1.AND.IM.EQ.NS+2.AND.ITRY(INUMT).GE.1) THEN
        IF(KFLD(INUMT).NE.21.AND.K(IP1+2,4).GT.0) INUM=3-INUM
        IF(KFLD(INUMT).EQ.21.AND.K(IP1+3,4).GT.0) INUM=3-INUM
      ENDIF
       
C...Store information on choice of evolving daughter.
      IEP(1)=N+INUM
      DO 370 I=2,NEP
        IEP(I)=IEP(I-1)+1
        IF(IEP(I).GT.N+NEP) IEP(I)=N+1
  370 CONTINUE
      DO 380 I=1,NEP
        KFL(I)=IABS(K(IEP(I),2))
  380 CONTINUE
      ITRY(INUM)=ITRY(INUM)+1
      IF(ITRY(INUM).GT.200) THEN
        CALL PYERRM(14,'(PYSHOW:) caught in infinite loop')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      Z=0.5D0
      IF(KFL(1).GT.40) GOTO 430
      IF(KSH(KFL(1)).EQ.0) GOTO 430
      IFL=KFL(1)
      IF(KFL(1).GE.6.AND.KFL(1).LE.8) IFL=37+KFL(1)+
     &ISIGN(2,K(IEP(1),2))
      IF(P(IEP(1),5).LT.PMTH(2,IFL)) GOTO 430

C...Check if evolution already predetermined for daughter.
      IPSPD=0
      IF(MPSPD.EQ.1.AND.IGM.EQ.0) THEN
        IF(K(IP1-1+INUM,4).GT.0) IPSPD=IP1-1+INUM
      ELSEIF(MPSPD.EQ.1.AND.IM.EQ.NS+2) THEN
        IF(KFL(1).NE.21.AND.K(IP1+2,4).GT.0) IPSPD=IP1+2
        IF(KFL(1).EQ.21.AND.K(IP1+3,4).GT.0) IPSPD=IP1+3
      ENDIF
      ISSET(INUM)=0
      IF(IPSPD.NE.0) ISSET(INUM)=1  
 
C...Select side for interference with initial state partons.
      IF(MIIS.GE.1.AND.IEP(1).LE.NS+3) THEN
        III=IEP(1)-NS-1
        ISII(III)=0
        IF(IABS(KCII(III)).EQ.1.AND.NIIS(III).EQ.1) THEN
          ISII(III)=1
        ELSEIF(KCII(III).EQ.2.AND.NIIS(III).EQ.1) THEN
          IF(PYR(0).GT.0.5D0) ISII(III)=1
        ELSEIF(KCII(III).EQ.2.AND.NIIS(III).EQ.2) THEN
          ISII(III)=1
          IF(PYR(0).GT.0.5D0) ISII(III)=2
        ENDIF
      ENDIF
 
C...Calculate allowed z range.
      IF(NEP.EQ.1) THEN
        PMED=PS(4)
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN
        PMED=P(IM,5)
      ELSE
        IF(INUM.EQ.1) PMED=V(IM,1)*PEM
        IF(INUM.EQ.2) PMED=(1D0-V(IM,1))*PEM
      ENDIF
      IF(MOD(MSTJ(43),2).EQ.1) THEN
        ZC=PMTH(2,21)/PMED
        ZCE=PMTH(2,22)/PMED
        IF(KFL(1).GE.11.AND.KFL(1).LE.18) ZCE=0.5D0*PARJ(90)/PMED
      ELSE
        ZC=0.5D0*(1D0-SQRT(MAX(0D0,1D0-(2D0*PMTH(2,21)/PMED)**2)))
        IF(ZC.LT.1D-6) ZC=(PMTH(2,21)/PMED)**2
        PMTMPE=PMTH(2,22)
        IF(KFL(1).GE.11.AND.KFL(1).LE.18) PMTMPE=0.5D0*PARJ(90)
        ZCE=0.5D0*(1D0-SQRT(MAX(0D0,1D0-(2D0*PMTMPE/PMED)**2)))
        IF(ZCE.LT.1D-6) ZCE=(PMTMPE/PMED)**2
      ENDIF
      ZC=MIN(ZC,0.491D0)
      ZCE=MIN(ZCE,0.49991D0)
      IF(((MSTJ(41).EQ.1.AND.ZC.GT.0.49D0).OR.(MSTJ(41).GE.2.AND.
     &MIN(ZC,ZCE).GT.0.4999D0)).AND.IPSPD.EQ.0) THEN
        P(IEP(1),5)=PMTH(1,IFL)
        V(IEP(1),5)=P(IEP(1),5)**2
        GOTO 430
      ENDIF
 
C...Integral of Altarelli-Parisi z kernel for QCD.
      IF(MSTJ(49).EQ.0.AND.KFL(1).EQ.21) THEN
        FBR=6D0*LOG((1D0-ZC)/ZC)+MSTJ(45)*0.5D0
      ELSEIF(MSTJ(49).EQ.0) THEN
        FBR=(8D0/3D0)*LOG((1D0-ZC)/ZC)
 
C...Integral of Altarelli-Parisi z kernel for scalar gluon.
      ELSEIF(MSTJ(49).EQ.1.AND.KFL(1).EQ.21) THEN
        FBR=(PARJ(87)+MSTJ(45)*PARJ(88))*(1D0-2D0*ZC)
      ELSEIF(MSTJ(49).EQ.1) THEN
        FBR=(1D0-2D0*ZC)/3D0
        IF(IGM.EQ.0.AND.M3JC.EQ.1) FBR=4D0*FBR
 
C...Integral of Altarelli-Parisi z kernel for Abelian vector gluon.
      ELSEIF(KFL(1).EQ.21) THEN
        FBR=6D0*MSTJ(45)*(0.5D0-ZC)
      ELSE
        FBR=2D0*LOG((1D0-ZC)/ZC)
      ENDIF
 
C...Reset QCD probability for lepton.
      IF(KFL(1).GE.11.AND.KFL(1).LE.18) FBR=0D0
 
C...Integral of Altarelli-Parisi kernel for photon emission.
      IF(MSTJ(41).GE.2.AND.KFL(1).GE.1.AND.KFL(1).LE.18) THEN
        FBRE=(KCHG(KFL(1),1)/3D0)**2*2D0*LOG((1D0-ZCE)/ZCE)
        IF(MSTJ(41).EQ.10) FBRE=PARJ(84)*FBRE
      ENDIF
 
C...Inner veto algorithm starts. Find maximum mass for evolution.
  390 PMS=V(IEP(1),5)
      IF(IGM.GE.0) THEN
        PM2=0D0
        DO 400 I=2,NEP
          PM=P(IEP(I),5)
          IF(KFL(I).LE.40) THEN
            IFLI=KFL(I)
            IF(KFL(I).GE.6.AND.KFL(I).LE.8) IFLI=37+KFL(I)+
     &      ISIGN(2,K(IEP(I),2))
            IF(KSH(KFL(I)).EQ.1) PM=PMTH(2,IFLI)
          ENDIF
          PM2=PM2+PM
  400   CONTINUE
        PMS=MIN(PMS,(P(IM,5)-PM2)**2)
      ENDIF
 
C...Select mass for daughter in QCD evolution.
      B0=27D0/6D0
      DO 410 IFF=4,MSTJ(45)
        IF(PMS.GT.4D0*PMTH(2,IFF)**2) B0=(33D0-2D0*IFF)/6D0
  410 CONTINUE
C...Already predetermined choice.
      IF(IPSPD.NE.0) THEN
        PMSQCD=P(IPSPD,5)**2
      ELSEIF(FBR.LT.1D-3) THEN
        PMSQCD=0D0
      ELSEIF(MSTJ(44).LE.0) THEN
        PMSQCD=PMS*EXP(MAX(-50D0,LOG(PYR(0))*PARU(2)/(PARU(111)*FBR)))
      ELSEIF(MSTJ(44).EQ.1) THEN
        PMSQCD=4D0*ALAMS*(0.25D0*PMS/ALAMS)**(PYR(0)**(B0/FBR))
      ELSE
        PMSQCD=PMS*EXP(MAX(-50D0,ALFM*B0*LOG(PYR(0))/FBR))
      ENDIF
      IF(ZC.GT.0.49D0.OR.PMSQCD.LE.PMTH(4,IFL)**2) PMSQCD=
     &  PMTH(2,IFL)**2
      V(IEP(1),5)=PMSQCD
      MCE=1
 
C...Select mass for daughter in QED evolution.
      IF(MSTJ(41).GE.2.AND.KFL(1).GE.1.AND.KFL(1).LE.18.AND.
     &IPSPD.EQ.0) THEN
        PMSQED=PMS*EXP(MAX(-50D0,LOG(PYR(0))*PARU(2)/(PARU(101)*FBRE)))
        IF(ZCE.GT.0.4999D0.OR.PMSQED.LE.PMTH(5,IFL)**2) PMSQED=
     &  PMTH(2,IFL)**2
        IF(PMSQED.GT.PMSQCD) THEN
          V(IEP(1),5)=PMSQED
          MCE=2
        ENDIF
      ENDIF
 
C...Check whether daughter mass below cutoff.
      P(IEP(1),5)=SQRT(V(IEP(1),5))
      IF(P(IEP(1),5).LE.PMTH(3,IFL)) THEN
        P(IEP(1),5)=PMTH(1,IFL)
        V(IEP(1),5)=P(IEP(1),5)**2
        GOTO 430
      ENDIF

C...Already predetermined choice of z, and flavour in g -> qqbar.
      IF(IPSPD.NE.0) THEN
        IPSGD1=K(IPSPD,4)
        IPSGD2=K(IPSPD,5)
        PMSGD1=P(IPSGD1,5)**2
        PMSGD2=P(IPSGD2,5)**2
        ALAMPS=SQRT(MAX(1D-10,(PMSQCD-PMSGD1-PMSGD2)**2-
     &  4D0*PMSGD1*PMSGD2))
        Z=0.5D0*(PMSQCD*(2D0*P(IPSGD1,4)/P(IPSPD,4)-1D0)+ALAMPS-
     &  PMSGD1+PMSGD2)/ALAMPS
        Z=MAX(0.00001D0,MIN(0.99999D0,Z))
        IF(KFL(1).NE.21) THEN
          K(IEP(1),5)=21
        ELSE
          K(IEP(1),5)=IABS(K(IPSGD1,2))
        ENDIF
 
C...Select z value of branching: q -> qgamma.
      ELSEIF(MCE.EQ.2) THEN
        Z=1D0-(1D0-ZCE)*(ZCE/(1D0-ZCE))**PYR(0)
        IF(1D0+Z**2.LT.2D0*PYR(0)) GOTO 390
        K(IEP(1),5)=22
 
C...Select z value of branching: q -> qg, g -> gg, g -> qqbar.
      ELSEIF(MSTJ(49).NE.1.AND.KFL(1).NE.21) THEN
        Z=1D0-(1D0-ZC)*(ZC/(1D0-ZC))**PYR(0)
        IF(IGM.EQ.0.AND.M3JCM.EQ.1) Z=1D0-(1D0-Z)**RESCZ
        IF(1D0+Z**2.LT.2D0*PYR(0)) GOTO 390
        K(IEP(1),5)=21
      ELSEIF(MSTJ(49).EQ.0.AND.MSTJ(45)*0.5D0.LT.PYR(0)*FBR) THEN
        Z=(1D0-ZC)*(ZC/(1D0-ZC))**PYR(0)
        IF(PYR(0).GT.0.5D0) Z=1D0-Z
        IF((1D0-Z*(1D0-Z))**2.LT.PYR(0)) GOTO 390
        K(IEP(1),5)=21
      ELSEIF(MSTJ(49).NE.1) THEN
        Z=PYR(0)
        IF(Z**2+(1D0-Z)**2.LT.PYR(0)) GOTO 390
        KFLB=1+INT(MSTJ(45)*PYR(0))
        PMQ=4D0*PMTH(2,KFLB)**2/V(IEP(1),5)
        IF(PMQ.GE.1D0) GOTO 390
        IF(MSTJ(44).LE.2) THEN
          IF(Z.LT.ZC.OR.Z.GT.1D0-ZC) GOTO 390
          PMQ0=4D0*PMTH(2,21)**2/V(IEP(1),5)
          IF(MOD(MSTJ(43),2).EQ.0.AND.(1D0+0.5D0*PMQ)*SQRT(1D0-PMQ)
     &    .LT.PYR(0)*(1D0+0.5D0*PMQ0)*SQRT(1D0-PMQ0)) GOTO 390
        ELSE
          IF((1D0+0.5D0*PMQ)*SQRT(1D0-PMQ).LT.PYR(0)) GOTO 390
        ENDIF
        K(IEP(1),5)=KFLB
 
C...Ditto for scalar gluon model.
      ELSEIF(KFL(1).NE.21) THEN
        Z=1D0-SQRT(ZC**2+PYR(0)*(1D0-2D0*ZC))
        K(IEP(1),5)=21
      ELSEIF(PYR(0)*(PARJ(87)+MSTJ(45)*PARJ(88)).LE.PARJ(87)) THEN
        Z=ZC+(1D0-2D0*ZC)*PYR(0)
        K(IEP(1),5)=21
      ELSE
        Z=ZC+(1D0-2D0*ZC)*PYR(0)
        KFLB=1+INT(MSTJ(45)*PYR(0))
        PMQ=4D0*PMTH(2,KFLB)**2/V(IEP(1),5)
        IF(PMQ.GE.1D0) GOTO 390
        K(IEP(1),5)=KFLB
      ENDIF

C...Correct to alpha_s(pT^2) (optionally m^2/4 for g -> q qbar).
      IF(MCE.EQ.1.AND.MSTJ(44).GE.2.AND.IPSPD.EQ.0) THEN
        IF(KFL(1).EQ.21.AND.K(IEP(1),5).LT.10.AND.MSTJ(44).EQ.3) THEN
          IF(ALFM/LOG(V(IEP(1),5)*0.25D0/ALAMS).LT.PYR(0)) GOTO 390
        ELSE
          IF(Z*(1D0-Z)*V(IEP(1),5).LT.PT2MIN) GOTO 390
          IF(ALFM/LOG(V(IEP(1),5)*Z*(1D0-Z)/ALAMS).LT.PYR(0)) GOTO 390
        ENDIF
      ENDIF
 
C...Check if z consistent with chosen m.
      IF(KFL(1).EQ.21) THEN
        KFLGD1=IABS(K(IEP(1),5))
        KFLGD2=KFLGD1
      ELSE
        KFLGD1=KFL(1)
        KFLGD2=IABS(K(IEP(1),5))
      ENDIF
      IF(NEP.EQ.1) THEN
        PED=PS(4)
      ELSEIF(NEP.GE.3) THEN
        PED=P(IEP(1),4)
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN
        PED=0.5D0*(V(IM,5)+V(IEP(1),5)-PM2**2)/P(IM,5)
      ELSE
        IF(IEP(1).EQ.N+1) PED=V(IM,1)*PEM
        IF(IEP(1).EQ.N+2) PED=(1D0-V(IM,1))*PEM
      ENDIF
      IF(MOD(MSTJ(43),2).EQ.1) THEN
        IFLGD1=KFLGD1
        IF(KFLGD1.GE.6.AND.KFLGD1.LE.8) IFLGD1=IFL
        PMQTH3=0.5D0*PARJ(82)
        IF(KFLGD2.EQ.22) PMQTH3=0.5D0*PARJ(83)
        IF(KFL(1).GE.11.AND.KFL(1).LE.18) PMQTH3=0.5D0*PARJ(90)
        PMQ1=(PMTH(1,IFLGD1)**2+PMQTH3**2)/V(IEP(1),5)
        PMQ2=(PMTH(1,KFLGD2)**2+PMQTH3**2)/V(IEP(1),5)
        ZD=SQRT(MAX(0D0,(1D0-V(IEP(1),5)/PED**2)*((1D0-PMQ1-PMQ2)**2-
     &  4D0*PMQ1*PMQ2)))
        ZH=1D0+PMQ1-PMQ2
      ELSE
        ZD=SQRT(MAX(0D0,1D0-V(IEP(1),5)/PED**2))
        ZH=1D0
      ENDIF
      IF(KFL(1).EQ.21.AND.K(IEP(1),5).LT.10.AND.MSTJ(44).EQ.3) THEN
      ELSEIF(IPSPD.NE.0) THEN
      ELSE 
        ZL=0.5D0*(ZH-ZD)
        ZU=0.5D0*(ZH+ZD)
        IF(Z.LT.ZL.OR.Z.GT.ZU) GOTO 390
      ENDIF
      IF(KFL(1).EQ.21) V(IEP(1),3)=LOG(ZU*(1D0-ZL)/MAX(1D-20,ZL*
     &(1D0-ZU)))
      IF(KFL(1).NE.21) V(IEP(1),3)=LOG((1D0-ZL)/MAX(1D-10,1D0-ZU))
 
C...Width suppression for q -> q + g.
      IF(MSTJ(40).NE.0.AND.KFL(1).NE.21.AND.IPSPD.EQ.0) THEN
        IF(IGM.EQ.0) THEN
          EGLU=0.5D0*PS(5)*(1D0-Z)*(1D0+V(IEP(1),5)/V(NS+1,5))
        ELSE
          EGLU=PMED*(1D0-Z)
        ENDIF
        CHI=PARJ(89)**2/(PARJ(89)**2+EGLU**2)
        IF(MSTJ(40).EQ.1) THEN
          IF(CHI.LT.PYR(0)) GOTO 390
        ELSEIF(MSTJ(40).EQ.2) THEN
          IF(1D0-CHI.LT.PYR(0)) GOTO 390
        ENDIF
      ENDIF
 
C...Three-jet matrix element correction (on both sides).
      IF(IGM.EQ.0.AND.M3JC.EQ.1) THEN
        X1=Z*(1D0+V(IEP(1),5)/V(NS+1,5))
        X2=1D0-V(IEP(1),5)/V(NS+1,5)
        X3=(1D0-X1)+(1D0-X2)
        IF(MCE.EQ.2) THEN
          KI1=K(IPA(INUM),2)
          KI2=K(IPA(3-INUM),2)
          QF1=KCHG(IABS(KI1),1)*ISIGN(1,KI1)/3D0
          QF2=KCHG(IABS(KI2),1)*ISIGN(1,KI2)/3D0
          WSHOW=QF1**2*(1D0-X1)/X3*(1D0+(X1/(2D0-X2))**2)+
     &    QF2**2*(1D0-X2)/X3*(1D0+(X2/(2D0-X1))**2)
          WME=(QF1*(1D0-X1)/X3-QF2*(1D0-X2)/X3)**2*(X1**2+X2**2)
        ELSEIF(MSTJ(49).NE.1.AND.M3JCM.NE.1) THEN
          WSHOW=1D0+(1D0-X1)/X3*(X1/(2D0-X2))**2+
     &    (1D0-X2)/X3*(X2/(2D0-X1))**2
          WME=X1**2+X2**2
        ELSEIF(MSTJ(49).NE.1) THEN
          X1=(1D0+(V(IEP(1),5)-PQMES)/V(NS+1,5))*
     &    (Z+(1D0-Z)*PQMES/V(IEP(1),5))
          X2=1D0-(V(IEP(1),5)-PQMES)/V(NS+1,5)
          X3=(1D0-X1)+(1D0-X2)
          Z1SH=(X1-(PQMES/V(NS+1,5))*(X3/MAX(1D-10,1D0-X2)))/(2D0-X2)
          Z2SH=(X2-(PQMES/V(NS+1,5))*(X3/MAX(1D-10,1D0-X1)))/(2D0-X1)
          WSHOW=(((1D0-X1)/(2D0-X2))*(1D0+Z1SH**2)/MAX(1D-10,1D0-Z1SH)+
     &    ((1D0-X2)/(2D0-X1))*(1D0+Z2SH**2)/MAX(1D-10,1D0-Z2SH))/RESCZ
          WME=X1**2+X2**2-QME*X3-0.5D0*QME**2-
     &    (0.5D0*QME+0.25D0*QME**2)*((1D0-X2)/MAX(1D-10,1D0-X1)+
     &    (1D0-X1)/MAX(1D-10,1D0-X2))
        ELSE
          WSHOW=4D0*X3*((1D0-X1)/(2D0-X2)**2+(1D0-X2)/(2D0-X1)**2)
          WME=X3**2
          IF(MSTJ(102).GE.2) WME=X3**2-2D0*(1D0+X3)*(1D0-X1)*(1D0-X2)*
     &    PARJ(171)
        ENDIF
        IF(WME.LT.PYR(0)*WSHOW) GOTO 390
 
C...Impose angular ordering by rejection of nonordered emission.
      ELSEIF(MCE.EQ.1.AND.IGM.GT.0.AND.MSTJ(42).GE.2.AND.IPSPD.EQ.0) 
     &THEN
        PEMAO=V(IM,1)*P(IM,4)
        IF(IEP(1).EQ.N+2) PEMAO=(1D0-V(IM,1))*P(IM,4)
        IF(KFL(1).EQ.21.AND.K(IEP(1),5).LE.10.AND.MSTJ(42).EQ.4) THEN
          MAOD=0
        ELSEIF(KFL(1).EQ.21.AND.K(IEP(1),5).LE.10.AND.MSTJ(42).EQ.3) 
     &  THEN
          MAOD=1
          PMDAO=PMTH(2,K(IEP(1),5))
          THE2ID=Z*(1D0-Z)*PEMAO**2/(V(IEP(1),5)-4D0*PMDAO**2)
        ELSE
          MAOD=1
          THE2ID=Z*(1D0-Z)*PEMAO**2/V(IEP(1),5)
        ENDIF
        MAOM=1
        IAOM=IM
  420   IF(K(IAOM,5).EQ.22) THEN
          IAOM=K(IAOM,3)
          IF(K(IAOM,3).LE.NS) MAOM=0
          IF(MAOM.EQ.1) GOTO 420
        ENDIF
        IF(MAOM.EQ.1.AND.MAOD.EQ.1) THEN
          THE2IM=V(IAOM,1)*(1D0-V(IAOM,1))*P(IAOM,4)**2/V(IAOM,5)
          IF(THE2ID.LT.THE2IM) GOTO 390
        ENDIF
      ENDIF
 
C...Impose user-defined maximum angle at first branching.
      IF(MSTJ(48).EQ.1.AND.IPSPD.EQ.0) THEN
        IF(NEP.EQ.1.AND.IM.EQ.NS) THEN
          THE2ID=Z*(1D0-Z)*PS(4)**2/V(IEP(1),5)
          IF(PARJ(85)**2*THE2ID.LT.1D0) GOTO 390
        ELSEIF(NEP.EQ.2.AND.IEP(1).EQ.NS+2) THEN
          THE2ID=Z*(1D0-Z)*(0.5D0*P(IM,4))**2/V(IEP(1),5)
          IF(PARJ(85)**2*THE2ID.LT.1D0) GOTO 390
        ELSEIF(NEP.EQ.2.AND.IEP(1).EQ.NS+3) THEN
          THE2ID=Z*(1D0-Z)*(0.5D0*P(IM,4))**2/V(IEP(1),5)
          IF(PARJ(86)**2*THE2ID.LT.1D0) GOTO 390
        ENDIF
      ENDIF
 
C...Impose angular constraint in first branching from interference
C...with initial state partons.
      IF(MIIS.GE.2.AND.IEP(1).LE.NS+3) THEN
        THE2D=MAX((1D0-Z)/Z,Z/(1D0-Z))*V(IEP(1),5)/(0.5D0*P(IM,4))**2
        IF(IEP(1).EQ.NS+2.AND.ISII(1).GE.1) THEN
          IF(THE2D.GT.THEIIS(1,ISII(1))**2) GOTO 390
        ELSEIF(IEP(1).EQ.NS+3.AND.ISII(2).GE.1) THEN
          IF(THE2D.GT.THEIIS(2,ISII(2))**2) GOTO 390
        ENDIF
      ENDIF
 
C...End of inner veto algorithm. Check if only one leg evolved so far.
  430 V(IEP(1),1)=Z
      ISL(1)=0
      ISL(2)=0
      IF(NEP.EQ.1) GOTO 460
      IF(NEP.EQ.2.AND.P(IEP(1),5)+P(IEP(2),5).GE.P(IM,5)) GOTO 330
      DO 440 I=1,NEP
        IF(ITRY(I).EQ.0.AND.KFLD(I).LE.40) THEN
          IF(KSH(KFLD(I)).EQ.1) THEN
            IFLD=KFLD(I)
            IF(KFLD(I).GE.6.AND.KFLD(I).LE.8) IFLD=37+KFLD(I)+
     &      ISIGN(2,K(N+I,2))
            IF(P(N+I,5).GE.PMTH(2,IFLD)) GOTO 330
          ENDIF
        ENDIF
  440 CONTINUE
 
C...Check if chosen multiplet m1,m2,z1,z2 is physical.
      IF(NEP.EQ.3) THEN
        PA1S=(P(N+1,4)+P(N+1,5))*(P(N+1,4)-P(N+1,5))
        PA2S=(P(N+2,4)+P(N+2,5))*(P(N+2,4)-P(N+2,5))
        PA3S=(P(N+3,4)+P(N+3,5))*(P(N+3,4)-P(N+3,5))
        PTS=0.25D0*(2D0*PA1S*PA2S+2D0*PA1S*PA3S+2D0*PA2S*PA3S-
     &  PA1S**2-PA2S**2-PA3S**2)/PA1S
        IF(PTS.LE.0D0) GOTO 330
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2.OR.MOD(MSTJ(43),2).EQ.0) THEN
        DO 450 I1=N+1,N+2
          KFLDA=IABS(K(I1,2))
          IF(KFLDA.GT.40) GOTO 450
          IF(KSH(KFLDA).EQ.0) GOTO 450
          IFLDA=KFLDA
          IF(KFLDA.GE.6.AND.KFLDA.LE.8) IFLDA=37+KFLDA+
     &    ISIGN(2,K(I1,2))
          IF(P(I1,5).LT.PMTH(2,IFLDA)) GOTO 450
          IF(KFLDA.EQ.21) THEN
            KFLGD1=IABS(K(I1,5))
            KFLGD2=KFLGD1
          ELSE
            KFLGD1=KFLDA
            KFLGD2=IABS(K(I1,5))
          ENDIF
          I2=2*N+3-I1
          IF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN
            PED=0.5D0*(V(IM,5)+V(I1,5)-V(I2,5))/P(IM,5)
          ELSE
            IF(I1.EQ.N+1) ZM=V(IM,1)
            IF(I1.EQ.N+2) ZM=1D0-V(IM,1)
            PML=SQRT((V(IM,5)-V(N+1,5)-V(N+2,5))**2-
     &      4D0*V(N+1,5)*V(N+2,5))
            PED=PEM*(0.5D0*(V(IM,5)-PML+V(I1,5)-V(I2,5))+PML*ZM)/
     &      V(IM,5)
          ENDIF
          IF(MOD(MSTJ(43),2).EQ.1) THEN
            PMQTH3=0.5D0*PARJ(82)
            IF(KFLGD2.EQ.22) PMQTH3=0.5D0*PARJ(83)
            IF(KFLDA.GE.11.AND.KFLDA.LE.18) PMQTH3=0.5D0*PARJ(90)
            IFLGD1=KFLGD1
            IF(KFLGD1.GE.6.AND.KFLGD1.LE.8) IFLGD1=IFLDA
            PMQ1=(PMTH(1,IFLGD1)**2+PMQTH3**2)/V(I1,5)
            PMQ2=(PMTH(1,KFLGD2)**2+PMQTH3**2)/V(I1,5)
            ZD=SQRT(MAX(0D0,(1D0-V(I1,5)/PED**2)*((1D0-PMQ1-PMQ2)**2-
     &      4D0*PMQ1*PMQ2)))
            ZH=1D0+PMQ1-PMQ2
          ELSE
            ZD=SQRT(MAX(0D0,1D0-V(I1,5)/PED**2))
            ZH=1D0
          ENDIF
          IF(KFLDA.EQ.21.AND.KFLGD1.LT.10.AND.MSTJ(44).EQ.3) THEN
          ELSE 
            ZL=0.5D0*(ZH-ZD)
            ZU=0.5D0*(ZH+ZD)
            IF(I1.EQ.N+1.AND.(V(I1,1).LT.ZL.OR.V(I1,1).GT.ZU).AND.
     &      ISSET(1).EQ.0) THEN
              ISL(1)=1
            ELSEIF(I1.EQ.N+2.AND.(V(I1,1).LT.ZL.OR.V(I1,1).GT.ZU).AND.
     &      ISSET(2).EQ.0) THEN 
              ISL(2)=1
            ENDIF
          ENDIF
          IF(KFLDA.EQ.21) V(I1,4)=LOG(ZU*(1D0-ZL)/MAX(1D-20,
     &    ZL*(1D0-ZU)))
          IF(KFLDA.NE.21) V(I1,4)=LOG((1D0-ZL)/MAX(1D-10,1D0-ZU))
  450   CONTINUE
        IF(ISL(1).EQ.1.AND.ISL(2).EQ.1.AND.ISLM.NE.0) THEN
          ISL(3-ISLM)=0
          ISLM=3-ISLM
        ELSEIF(ISL(1).EQ.1.AND.ISL(2).EQ.1) THEN
          ZDR1=MAX(0D0,V(N+1,3)/MAX(1D-6,V(N+1,4))-1D0)
          ZDR2=MAX(0D0,V(N+2,3)/MAX(1D-6,V(N+2,4))-1D0)
          IF(ZDR2.GT.PYR(0)*(ZDR1+ZDR2)) ISL(1)=0
          IF(ISL(1).EQ.1) ISL(2)=0
          IF(ISL(1).EQ.0) ISLM=1
          IF(ISL(2).EQ.0) ISLM=2
        ENDIF
        IF(ISL(1).EQ.1.OR.ISL(2).EQ.1) GOTO 330
      ENDIF
      IFLD1=KFLD(1)
      IF(KFLD(1).GE.6.AND.KFLD(1).LE.8) IFLD1=37+KFLD(1)+
     &ISIGN(2,K(N+1,2))
      IFLD2=KFLD(2)
      IF(KFLD(2).GE.6.AND.KFLD(2).LE.8) IFLD2=37+KFLD(2)+
     &ISIGN(2,K(N+2,2))
      IF(IGM.GT.0) THEN
        IF(MOD(MSTJ(43),2).EQ.1.AND.(P(N+1,5).GE.
     &  PMTH(2,IFLD1).OR.P(N+2,5).GE.PMTH(2,IFLD2))) THEN
          PMQ1=V(N+1,5)/V(IM,5)
          PMQ2=V(N+2,5)/V(IM,5)
          ZD=SQRT(MAX(0D0,(1D0-V(IM,5)/PEM**2)*((1D0-PMQ1-PMQ2)**2-
     &    4D0*PMQ1*PMQ2)))
          ZH=1D0+PMQ1-PMQ2
          ZL=0.5D0*(ZH-ZD)
          ZU=0.5D0*(ZH+ZD)
          IF(V(IM,1).LT.ZL.OR.V(IM,1).GT.ZU) GOTO 330
        ENDIF
      ENDIF
 
C...Accepted branch. Construct four-momentum for initial partons.
  460 MAZIP=0
      MAZIC=0
      IF(NEP.EQ.1) THEN
        P(N+1,1)=0D0
        P(N+1,2)=0D0
        P(N+1,3)=SQRT(MAX(0D0,(P(IPA(1),4)+P(N+1,5))*(P(IPA(1),4)-
     &  P(N+1,5))))
        P(N+1,4)=P(IPA(1),4)
        V(N+1,2)=P(N+1,4)
      ELSEIF(IGM.EQ.0.AND.NEP.EQ.2) THEN
        PED1=0.5D0*(V(IM,5)+V(N+1,5)-V(N+2,5))/P(IM,5)
        P(N+1,1)=0D0
        P(N+1,2)=0D0
        P(N+1,3)=SQRT(MAX(0D0,(PED1+P(N+1,5))*(PED1-P(N+1,5))))
        P(N+1,4)=PED1
        P(N+2,1)=0D0
        P(N+2,2)=0D0
        P(N+2,3)=-P(N+1,3)
        P(N+2,4)=P(IM,5)-PED1
        V(N+1,2)=P(N+1,4)
        V(N+2,2)=P(N+2,4)
      ELSEIF(NEP.EQ.3) THEN
        P(N+1,1)=0D0
        P(N+1,2)=0D0
        P(N+1,3)=SQRT(MAX(0D0,PA1S))
        P(N+2,1)=SQRT(PTS)
        P(N+2,2)=0D0
        P(N+2,3)=0.5D0*(PA3S-PA2S-PA1S)/P(N+1,3)
        P(N+3,1)=-P(N+2,1)
        P(N+3,2)=0D0
        P(N+3,3)=-(P(N+1,3)+P(N+2,3))
        V(N+1,2)=P(N+1,4)
        V(N+2,2)=P(N+2,4)
        V(N+3,2)=P(N+3,4)
 
C...Construct transverse momentum for ordinary branching in shower.
      ELSE
        ZM=V(IM,1)
        LOOPPT=0
  465   LOOPPT=LOOPPT+1
        PZM=SQRT(MAX(0D0,(PEM+P(IM,5))*(PEM-P(IM,5))))
        PMLS=(V(IM,5)-V(N+1,5)-V(N+2,5))**2-4D0*V(N+1,5)*V(N+2,5)
        IF(PZM.LE.0D0) THEN
          PTS=0D0
        ELSEIF(K(IM,2).EQ.21.AND.IABS(K(N+1,2)).LE.10.AND.
     &  MSTJ(44).EQ.3) THEN
          PTS=PMLS*ZM*(1D0-ZM)/V(IM,5) 
        ELSEIF(MOD(MSTJ(43),2).EQ.1) THEN
          PTS=(PEM**2*(ZM*(1D0-ZM)*V(IM,5)-(1D0-ZM)*V(N+1,5)-
     &    ZM*V(N+2,5))-0.25D0*PMLS)/PZM**2
        ELSE
          PTS=PMLS*(ZM*(1D0-ZM)*PEM**2/V(IM,5)-0.25D0)/PZM**2
        ENDIF
        IF(PTS.LT.0D0.AND.LOOPPT.LT.10) THEN
          ZM=0.05D0+0.9D0*ZM
          GOTO 465
        ELSEIF(PTS.LT.0D0) THEN 
          GOTO 265
        ENDIF 
        PT=SQRT(MAX(0D0,PTS))
 
C...Find coefficient of azimuthal asymmetry due to gluon polarization.
        HAZIP=0D0
        IF(MSTJ(49).NE.1.AND.MOD(MSTJ(46),2).EQ.1.AND.K(IM,2).EQ.21
     &  .AND.IAU.NE.0) THEN
          IF(K(IGM,3).NE.0) MAZIP=1
          ZAU=V(IGM,1)
          IF(IAU.EQ.IM+1) ZAU=1D0-V(IGM,1)
          IF(MAZIP.EQ.0) ZAU=0D0
          IF(K(IGM,2).NE.21) THEN
            HAZIP=2D0*ZAU/(1D0+ZAU**2)
          ELSE
            HAZIP=(ZAU/(1D0-ZAU*(1D0-ZAU)))**2
          ENDIF
          IF(K(N+1,2).NE.21) THEN
            HAZIP=HAZIP*(-2D0*ZM*(1D0-ZM))/(1D0-2D0*ZM*(1D0-ZM))
          ELSE
            HAZIP=HAZIP*(ZM*(1D0-ZM)/(1D0-ZM*(1D0-ZM)))**2
          ENDIF
        ENDIF
 
C...Find coefficient of azimuthal asymmetry due to soft gluon
C...interference.
        HAZIC=0D0
        IF(MSTJ(49).NE.2.AND.MSTJ(46).GE.2.AND.(K(N+1,2).EQ.21.OR.
     &  K(N+2,2).EQ.21).AND.IAU.NE.0) THEN
          IF(K(IGM,3).NE.0) MAZIC=N+1
          IF(K(IGM,3).NE.0.AND.K(N+1,2).NE.21) MAZIC=N+2
          IF(K(IGM,3).NE.0.AND.K(N+1,2).EQ.21.AND.K(N+2,2).EQ.21.AND.
     &    ZM.GT.0.5D0) MAZIC=N+2
          IF(K(IAU,2).EQ.22) MAZIC=0
          ZS=ZM
          IF(MAZIC.EQ.N+2) ZS=1D0-ZM
          ZGM=V(IGM,1)
          IF(IAU.EQ.IM-1) ZGM=1D0-V(IGM,1)
          IF(MAZIC.EQ.0) ZGM=1D0
          IF(MAZIC.NE.0) HAZIC=(P(IM,5)/P(IGM,5))*
     &    SQRT((1D0-ZS)*(1D0-ZGM)/(ZS*ZGM))
          HAZIC=MIN(0.95D0,HAZIC)
        ENDIF
      ENDIF
 
C...Construct energies for ordinary branching in shower.
  470 IF(NEP.EQ.2.AND.IGM.GT.0) THEN
        IF(K(IM,2).EQ.21.AND.IABS(K(N+1,2)).LE.10.AND.
     &  MSTJ(44).EQ.3) THEN
          P(N+1,4)=0.5D0*(PEM*(V(IM,5)+V(N+1,5)-V(N+2,5))+
     &    PZM*SQRT(MAX(0D0,PMLS))*(2D0*ZM-1D0))/V(IM,5)
        ELSEIF(MOD(MSTJ(43),2).EQ.1) THEN
          P(N+1,4)=PEM*V(IM,1)
        ELSE
          P(N+1,4)=PEM*(0.5D0*(V(IM,5)-SQRT(PMLS)+V(N+1,5)-V(N+2,5))+
     &    SQRT(PMLS)*ZM)/V(IM,5)
        ENDIF

C...Already predetermined choice of phi angle or not
        PHI=PARU(2)*PYR(0)
        IF(MPSPD.EQ.1.AND.IGM.EQ.NS+1) THEN
          IPSPD=IP1+IM-NS-2
          IF(K(IPSPD,4).GT.0) THEN
            IPSGD1=K(IPSPD,4)
            IF(IM.EQ.NS+2) THEN
              PHI=PYANGL(P(IPSGD1,1),P(IPSGD1,2))
            ELSE
              PHI=PYANGL(-P(IPSGD1,1),P(IPSGD1,2))
            ENDIF
          ENDIF
        ELSEIF(MPSPD.EQ.1.AND.IGM.EQ.NS+2) THEN
          IPSPD=IP1+IM-NS-2
          IF(K(IPSPD,4).GT.0) THEN
            IPSGD1=K(IPSPD,4)
            PHIPSM=PYANGL(P(IPSPD,1),P(IPSPD,2))
            THEPSM=PYANGL(P(IPSPD,3),SQRT(P(IPSPD,1)**2+P(IPSPD,2)**2))
            CALL PYROBO(IPSGD1,IPSGD1,0D0,-PHIPSM,0D0,0D0,0D0)
            CALL PYROBO(IPSGD1,IPSGD1,-THEPSM,0D0,0D0,0D0,0D0)
            PHI=PYANGL(P(IPSGD1,1),P(IPSGD1,2))
            CALL PYROBO(IPSGD1,IPSGD1,THEPSM,PHIPSM,0D0,0D0,0D0)   
          ENDIF
        ENDIF

C...Construct momenta for ordinary branching in shower.
        P(N+1,1)=PT*COS(PHI)
        P(N+1,2)=PT*SIN(PHI)
        IF(K(IM,2).EQ.21.AND.IABS(K(N+1,2)).LE.10.AND.
     &  MSTJ(44).EQ.3) THEN
          P(N+1,3)=0.5D0*(PZM*(V(IM,5)+V(N+1,5)-V(N+2,5))+
     &    PEM*SQRT(MAX(0D0,PMLS))*(2D0*ZM-1D0))/V(IM,5)
        ELSEIF(PZM.GT.0D0) THEN
          P(N+1,3)=0.5D0*(V(N+2,5)-V(N+1,5)-V(IM,5)+
     &    2D0*PEM*P(N+1,4))/PZM
        ELSE
          P(N+1,3)=0D0
        ENDIF
        P(N+2,1)=-P(N+1,1)
        P(N+2,2)=-P(N+1,2)
        P(N+2,3)=PZM-P(N+1,3)
        P(N+2,4)=PEM-P(N+1,4)
        IF(MSTJ(43).LE.2) THEN
          V(N+1,2)=(PEM*P(N+1,4)-PZM*P(N+1,3))/P(IM,5)
          V(N+2,2)=(PEM*P(N+2,4)-PZM*P(N+2,3))/P(IM,5)
        ENDIF
      ENDIF
 
C...Rotate and boost daughters.
      IF(IGM.GT.0) THEN
        IF(MSTJ(43).LE.2) THEN
          BEX=P(IGM,1)/P(IGM,4)
          BEY=P(IGM,2)/P(IGM,4)
          BEZ=P(IGM,3)/P(IGM,4)
          GA=P(IGM,4)/P(IGM,5)
          GABEP=GA*(GA*(BEX*P(IM,1)+BEY*P(IM,2)+BEZ*P(IM,3))/(1D0+GA)-
     &    P(IM,4))
        ELSE
          BEX=0D0
          BEY=0D0
          BEZ=0D0
          GA=1D0
          GABEP=0D0
        ENDIF
        PTIMB=SQRT((P(IM,1)+GABEP*BEX)**2+(P(IM,2)+GABEP*BEY)**2)
        THE=PYANGL(P(IM,3)+GABEP*BEZ,PTIMB)
        IF(PTIMB.GT.1D-4) THEN
          PHI=PYANGL(P(IM,1)+GABEP*BEX,P(IM,2)+GABEP*BEY)
        ELSE
          PHI=0D0
        ENDIF 
        DO 480 I=N+1,N+2
          DP(1)=COS(THE)*COS(PHI)*P(I,1)-SIN(PHI)*P(I,2)+
     &    SIN(THE)*COS(PHI)*P(I,3)
          DP(2)=COS(THE)*SIN(PHI)*P(I,1)+COS(PHI)*P(I,2)+
     &    SIN(THE)*SIN(PHI)*P(I,3)
          DP(3)=-SIN(THE)*P(I,1)+COS(THE)*P(I,3)
          DP(4)=P(I,4)
          DBP=BEX*DP(1)+BEY*DP(2)+BEZ*DP(3)
          DGABP=GA*(GA*DBP/(1D0+GA)+DP(4))
          P(I,1)=DP(1)+DGABP*BEX
          P(I,2)=DP(2)+DGABP*BEY
          P(I,3)=DP(3)+DGABP*BEZ
          P(I,4)=GA*(DP(4)+DBP)
  480   CONTINUE
      ENDIF
 
C...Weight with azimuthal distribution, if required.
      IF(MAZIP.NE.0.OR.MAZIC.NE.0) THEN
        DO 490 J=1,3
          DPT(1,J)=P(IM,J)
          DPT(2,J)=P(IAU,J)
          DPT(3,J)=P(N+1,J)
  490   CONTINUE
        DPMA=DPT(1,1)*DPT(2,1)+DPT(1,2)*DPT(2,2)+DPT(1,3)*DPT(2,3)
        DPMD=DPT(1,1)*DPT(3,1)+DPT(1,2)*DPT(3,2)+DPT(1,3)*DPT(3,3)
        DPMM=DPT(1,1)**2+DPT(1,2)**2+DPT(1,3)**2
        DO 500 J=1,3
          DPT(4,J)=DPT(2,J)-DPMA*DPT(1,J)/MAX(1D-10,DPMM)
          DPT(5,J)=DPT(3,J)-DPMD*DPT(1,J)/MAX(1D-10,DPMM)
  500   CONTINUE
        DPT(4,4)=SQRT(DPT(4,1)**2+DPT(4,2)**2+DPT(4,3)**2)
        DPT(5,4)=SQRT(DPT(5,1)**2+DPT(5,2)**2+DPT(5,3)**2)
        IF(MIN(DPT(4,4),DPT(5,4)).GT.0.1D0*PARJ(82)) THEN
          CAD=(DPT(4,1)*DPT(5,1)+DPT(4,2)*DPT(5,2)+
     &    DPT(4,3)*DPT(5,3))/(DPT(4,4)*DPT(5,4))
          IF(MAZIP.NE.0) THEN
            IF(1D0+HAZIP*(2D0*CAD**2-1D0).LT.PYR(0)*(1D0+ABS(HAZIP)))
     &      GOTO 470
          ENDIF
          IF(MAZIC.NE.0) THEN
            IF(MAZIC.EQ.N+2) CAD=-CAD
            IF((1D0-HAZIC)*(1D0-HAZIC*CAD)/(1D0+HAZIC**2-2D0*HAZIC*CAD)
     &      .LT.PYR(0)) GOTO 470
          ENDIF
        ENDIF
      ENDIF
 
C...Azimuthal anisotropy due to interference with initial state partons.
      IF(MOD(MIIS,2).EQ.1.AND.IGM.EQ.NS+1.AND.(K(N+1,2).EQ.21.OR.
     &K(N+2,2).EQ.21)) THEN
        III=IM-NS-1
        IF(ISII(III).GE.1) THEN
          IAZIID=N+1
          IF(K(N+1,2).NE.21) IAZIID=N+2
          IF(K(N+1,2).EQ.21.AND.K(N+2,2).EQ.21.AND.
     &    P(N+1,4).GT.P(N+2,4)) IAZIID=N+2
          THEIID=PYANGL(P(IAZIID,3),SQRT(P(IAZIID,1)**2+P(IAZIID,2)**2))
          IF(III.EQ.2) THEIID=PARU(1)-THEIID
          PHIIID=PYANGL(P(IAZIID,1),P(IAZIID,2))
          HAZII=MIN(0.95D0,THEIID/THEIIS(III,ISII(III)))
          CAD=COS(PHIIID-PHIIIS(III,ISII(III)))
          PHIREL=ABS(PHIIID-PHIIIS(III,ISII(III)))
          IF(PHIREL.GT.PARU(1)) PHIREL=PARU(2)-PHIREL
          IF((1D0-HAZII)*(1D0-HAZII*CAD)/(1D0+HAZII**2-2D0*HAZII*CAD)
     &    .LT.PYR(0)) GOTO 470
        ENDIF
      ENDIF
 
C...Continue loop over partons that may branch, until none left.
      IF(IGM.GE.0) K(IM,1)=14
      N=N+NEP
      NEP=2
      IF(N.GT.MSTU(4)-MSTU(32)-5) THEN
        CALL PYERRM(11,'(PYSHOW:) no more memory left in PYJETS')
        IF(MSTU(21).GE.1) N=NS
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      GOTO 270
 
C...Set information on imagined shower initiator.
  510 IF(NPA.GE.2) THEN
        K(NS+1,1)=11
        K(NS+1,2)=94
        K(NS+1,3)=IP1
        IF(IP2.GT.0.AND.IP2.LT.IP1) K(NS+1,3)=IP2
        K(NS+1,4)=NS+2
        K(NS+1,5)=NS+1+NPA
        IIM=1
      ELSE
        IIM=0
      ENDIF
 
C...Reconstruct string drawing information.
      DO 520 I=NS+1+IIM,N
        IF(K(I,1).LE.10.AND.K(I,2).EQ.22) THEN
          K(I,1)=1
        ELSEIF(K(I,1).LE.10.AND.IABS(K(I,2)).GE.11.AND.
     &    IABS(K(I,2)).LE.18) THEN
          K(I,1)=1
        ELSEIF(K(I,1).LE.10) THEN
          K(I,4)=MSTU(5)*(K(I,4)/MSTU(5))
          K(I,5)=MSTU(5)*(K(I,5)/MSTU(5))
        ELSEIF(K(MOD(K(I,4),MSTU(5))+1,2).NE.22) THEN
          ID1=MOD(K(I,4),MSTU(5))
          IF(K(I,2).GE.1.AND.K(I,2).LE.8) ID1=MOD(K(I,4),MSTU(5))+1
          ID2=2*MOD(K(I,4),MSTU(5))+1-ID1
          K(I,4)=MSTU(5)*(K(I,4)/MSTU(5))+ID1
          K(I,5)=MSTU(5)*(K(I,5)/MSTU(5))+ID2
          K(ID1,4)=K(ID1,4)+MSTU(5)*I
          K(ID1,5)=K(ID1,5)+MSTU(5)*ID2
          K(ID2,4)=K(ID2,4)+MSTU(5)*ID1
          K(ID2,5)=K(ID2,5)+MSTU(5)*I
        ELSE
          ID1=MOD(K(I,4),MSTU(5))
          ID2=ID1+1
          K(I,4)=MSTU(5)*(K(I,4)/MSTU(5))+ID1
          K(I,5)=MSTU(5)*(K(I,5)/MSTU(5))+ID1
          IF(IABS(K(I,2)).LE.10.OR.K(ID1,1).GE.11) THEN
            K(ID1,4)=K(ID1,4)+MSTU(5)*I
            K(ID1,5)=K(ID1,5)+MSTU(5)*I
          ELSE
            K(ID1,4)=0
            K(ID1,5)=0
          ENDIF
          K(ID2,4)=0
          K(ID2,5)=0
        ENDIF
  520 CONTINUE
 
C...Transformation from CM frame.
      IF(NPA.GE.2) THEN
        BEX=PS(1)/PS(4)
        BEY=PS(2)/PS(4)
        BEZ=PS(3)/PS(4)
        GA=PS(4)/PS(5)
        GABEP=GA*(GA*(BEX*P(IPA(1),1)+BEY*P(IPA(1),2)+BEZ*P(IPA(1),3))
     &  /(1D0+GA)-P(IPA(1),4))
      ELSE
        BEX=0D0
        BEY=0D0
        BEZ=0D0
        GABEP=0D0
      ENDIF
      THE=PYANGL(P(IPA(1),3)+GABEP*BEZ,SQRT((P(IPA(1),1)
     &+GABEP*BEX)**2+(P(IPA(1),2)+GABEP*BEY)**2))
      PHI=PYANGL(P(IPA(1),1)+GABEP*BEX,P(IPA(1),2)+GABEP*BEY)
      IF(NPA.EQ.3) THEN
        CHI=PYANGL(COS(THE)*COS(PHI)*(P(IPA(2),1)+GABEP*BEX)+COS(THE)*
     &  SIN(PHI)*(P(IPA(2),2)+GABEP*BEY)-SIN(THE)*(P(IPA(2),3)+GABEP*
     &  BEZ),-SIN(PHI)*(P(IPA(2),1)+GABEP*BEX)+COS(PHI)*(P(IPA(2),2)+
     &  GABEP*BEY))
        MSTU(33)=1
        CALL PYROBO(NS+1,N,0D0,CHI,0D0,0D0,0D0)
      ENDIF
      MSTU(33)=1
      CALL PYROBO(NS+1,N,THE,PHI,BEX,BEY,BEZ)
 
C...Decay vertex of shower.
      DO 540 I=NS+1,N
        DO 530 J=1,5
          V(I,J)=V(IP1,J)
  530   CONTINUE
  540 CONTINUE
 
C...Delete trivial shower, else connect initiators.
      IF(N.LE.NS+NPA+IIM) THEN
        N=NS
      ELSE
        DO 550 IP=1,NPA
          K(IPA(IP),1)=14
          K(IPA(IP),4)=K(IPA(IP),4)+NS+IIM+IP
          K(IPA(IP),5)=K(IPA(IP),5)+NS+IIM+IP
          K(NS+IIM+IP,3)=IPA(IP)
          IF(IIM.EQ.1.AND.MSTU(16).NE.2) K(NS+IIM+IP,3)=NS+1
          IF(K(NS+IIM+IP,1).NE.1) THEN
            K(NS+IIM+IP,4)=MSTU(5)*IPA(IP)+K(NS+IIM+IP,4)
            K(NS+IIM+IP,5)=MSTU(5)*IPA(IP)+K(NS+IIM+IP,5)
          ENDIF
  550   CONTINUE
      ENDIF
 
      RETURN
      END
