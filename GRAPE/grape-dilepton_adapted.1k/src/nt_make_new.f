      subroutine HB_init(LUN, NTPYT_flag, NTVEC_flag)
include "Pythia6ToHepMC3.inc"
      INTEGER OUTID, HEPMC3STATUS
      common/moj/ OUTID(2), HEPMC3STATUS
c      implicit NONE
* ---------- Argument ----------
      integer  LUN      
      logical  NTPYT_flag, NTVEC_flag
* ------------------------------
* ---------- HBOOK stuff ----------
      integer    NWRD_HLIMIT
       parameter(NWRD_HLIMIT = 500000)
      real*4        P
       common /PAWC/P(NWRD_HLIMIT)
      integer*4      IQUEST(100)
       common /QUEST/IQUEST
* ---------------------------------
* ---------- Ntuple varialbes ----------
      include 'nt_inc.h'
* --------------------------------------
* ------ Local variables ------
      integer  ierr
* -----------------------------
C      include './inc/graepia.h'
* --------- Initialization of HBOOK ----------
      call hlimit(NWRD_HLIMIT)
      IQUEST(10) = IQ10           
      call hbset('BSIZE', LREC, ierr)   
      if (ierr .NE. 0) then
         write(6,*) '!!!ERROR from HBSET in HB_init!!!'
         write(6,*) '  ---> Ierr =', ierr
         write(6,*) '  ---> Good-bye!'
         STOP
      endif
* ------------ Opening a Ntuple_file ------------
      call hropen(LUN , 'grp' , NT_NAME , 'NQ' , LREC , ierr)
      write(6,*) ' HROPEN called for file: ',NT_NAME
      write(6,*) ' LREC: ',LREC,' Error code: ', ierr
      if (ierr .NE. 0) then
         write(6,*) '!!!ERROR from HROPEN in HB_init!!!'
         write(6,*) '  ---> Ierr =', ierr
         write(6,*) '  ---> Good-bye!'
         STOP
      endif
      if (NTPYT_flag) then
        call hbnt(NTID, NT_TITLE, ' ')         
        call hbname(NTID, 'PYTHIA',  Npy,  PYTHIA)
      endif
      call hbnt(NTID+10, 'PROCESS parameters', ' ')
      call hbname(NTID+10, 'INTEGER4', nt_jproc,  NT_PRC_i4)
      call hbname(NTID+10, 'REAL8',    nt_P1_lab, NT_PRC_r8)
      if (NTVEC_flag) then
        call hbnt(NTID+20, 'EVENT variables', ' ')
        call hbname(NTID+20, 'INTEGER4',  nt_Nisr, NT_EVT_i4)
        call hbnt(NTID+30, '4-vectors', ' ')
        call hbname(NTID+30, 'VECTOR',  vec_px, NT_VEC)
        call hbnt(NTID+40, 'KF code', ' ')
        call hbname(NTID+40, 'KF code',  vec_kf, NT_KF)
        call hbnt(NTID+50, 'ISR photons', ' ')
        call hbname(NTID+50, 'ISR', nt_vec_isr(1), NT_ISR)
      endif
      OUTID(1)=hepmc3_new_writer(0,1,'Grape.hepmc'//char(0))
      NEVHEP=-123456
      HEPMC3STATUS=hepmc3_set_hepevt_address(NEVHEP)
      return
      end
* ==========================================================================
      subroutine  FILL_nt(LUN, Ngen, Ievt, merge, Nextn, Nisr
     &                                       ,NTPYT_flag,NTVEC_flag)
include "Pythia6ToHepMC3.inc"
      INTEGER OUTID, HEPMC3STATUS
      common/moj/ OUTID(2), HEPMC3STATUS
c      implicit NONE
* ---------- Argument ----------
      integer  LUN, Ngen, Ievt, merge, Nextn, Nisr
      logical  NTPYT_flag, NTVEC_flag
* ------------------------------
* -------- Ntuple variables --------
      include 'nt_inc.h'
* ----------------------------------
* ---------- GRACE stuff ----------
      integer  mxextn
       parameter(mxextn=10)
      double precision  vec(4,mxextn)
       common /sp4vec/  vec
      double precision  amass1(mxextn), amass2(mxextn)
       common /kmmass/  amass1,         amass2
      integer          kcharg(mxextn), kfcode(mxextn)
       common /kminfo/ kcharg,         kfcode
      integer          jproc
       common /amjprc/ jproc
* ---------------------------------
* ------------ BASES common on its result ------------
      integer           ITG,ITF
      real*4            STIME
      double precision  AVGI,SD,CHI2A
      COMMON /BSRSLT/AVGI,SD,CHI2A,STIME,ITG,ITF
* ----------------------------------------------------
* ---------- PYTHIA stuff ----------
      include './inc/py_common.h'
      integer nmxhep
      PARAMETER (NMXHEP=4000)
      integer NEVHEP, NHEP, ISTHEP, IDHEP, JMOHEP, JDAHEP
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP, VHEP
* ----------------------------------
* ---------- Kinematical variables ----------
      double precision  P1_lab,P2_lab, E1_lab,E2_lab  
     &                 ,Pcms_lab(3),Ecms_lab(3)
     &                 ,GAMMAcms_lab(3),BETGAMcms_lab(3)
     &                 ,vec_isr(4)
       common /GEP_LAB/ P1_lab,P2_lab, E1_lab,E2_lab
     &                 ,Pcms_lab,Ecms_lab
     &                 ,GAMMAcms_lab,   BETGAMcms_lab
     &                 ,vec_isr


       double precision fse(4),theta_e,beam_e(4),x_e
       double precision beam_p(4),x_p,fsp(4),ptp
       double precision fslp(4),eta_lp,theta_lp,fslm(4),eta_lm,theta_lm
* -------------------------------------------
* -------- Local variables --------
      integer  i,jc,jz, ientry
      data ientry /0/
* ---------------------------------
C      include './inc/graepia.h'
      ica = 1
      if(ientry.lt.1) then
        NEVHEP=-123456
        HEPMC3STATUS=hepmc3_set_hepevt_address(NEVHEP)
        ientry = 10
      endif
      call hcdir('//grp', ' ')   
******** Initialization of HBOOK *******
      if (Ievt .LE. 1) then
C         call HB_init(LUN, NTPYT_flag, NTVEC_flag)
         nt_jproc = jproc
         nt_merge = merge
         nt_Nextn = Nextn
         nt_Ngen  = Ngen
         nt_P1_lab = P1_lab
         nt_P2_lab = P2_lab
         nt_xsec(1) = AVGI   !!! x-sec
         nt_xsec(2) = SD     !!! error
         call hfnt(NTID+10)
      endif
****************************************
      if ((Ievt .GE. 1).and.(Ievt .LE. Ngen)) then !!!!!!!!!!!!!!
         if (NTPYT_flag) then
            Npy = min(N,N_max)
            do 1000 i=1, Npy          
              px(i) = p(i,1)
              py(i) = p(i,2)
              pz(i) = p(i,3)
              pe(i) = p(i,4)
              pm(i) = p(i,5)
              kf(i)  = K(i,2)
              sta(i) = K(i,1)
              mot(i) = K(i,3)
 1000       continue

c
c fs lepton cuts  0.5<E'/E<0.9
C
c  Theta_e>pi-0.010
c     
            beam_e(1) = px(2)
            beam_e(2) = py(2)
            beam_e(3) = pz(2)
            beam_e(4) = pe(2)
            beam_p(1) = px(1)
            beam_p(2) = py(1)
            beam_p(3) = pz(1)
            beam_p(4) = pe(1)
            do 1001 i = 1,Npy
               if(sta(i).eq.1) then
                  if(abs(kf(i)).eq.2212) then
                     fsp(1) = px(i)
                     fsp(2) = py(i)
                     fsp(3) = pz(i)
                     fsp(4) = pe(i)
                     ptp = sqrt(fsp(1)*fsp(1)+fsp(2)*fsp(2))
                     x_p = fsp(3)/beam_p(3)
                  endif
                  if(abs(kf(i)).eq.11) then
                     fse(1) = px(i)
                     fse(2) = py(i)
                     fse(3) = pz(i)
                     fse(4) = pe(i)
                  endif
                  if(kf(i).eq.15) then
                     fslm(1) = px(i)
                     fslm(2) = py(i)
                     fslm(3) = pz(i)
                     fslm(4) = pe(i)
                  endif
                  if(kf(i).eq.-15) then
                     fslp(1) = px(i)
                     fslp(2) = py(i)
                     fslp(3) = pz(i)
                     fslp(4) = pe(i)
                  endif            
               endif
 1001       continue
            theta_e = dacos(fse(3)
     +      /dsqrt(fse(1)*fse(1)+fse(2)*fse(2)+fse(3)*fse(3)))
            theta_lp = dacos(fslp(3)
     +      /dsqrt(fslp(1)*fslp(1)+fslp(2)*fslp(2)+fslp(3)*fslp(3)))
            eta_lp = -dlog(dtan(theta_lp/2.0))
            theta_lm = dacos(fslm(3)
     +      /dsqrt(fslm(1)*fslm(1)+fslm(2)*fslm(2)+fslm(3)*fslm(3)))
            eta_lm = -dlog(dtan(theta_lm/2.0))
            x_e = fse(4)/beam_e(4)
            if(0.47.lt.x_e.and.x_e.lt.0.93) then
               if(theta_e.gt.dacos(-1.d0)-0.01) then
                  if(abs(eta_lm).lt.4.1.and.abs(eta_lp).lt.4.1) then
                     if(x_p<0.98.or.ptp>0.09) then
                        call hfnt(NTID)
            call PYHEPC(1)
c            write(77,*) NEVHEP,nhep
c            do jc = 1,nhep
c	    write(77,10) isthep(jc),idhep(jc),jmohep(1,jc),jmohep(2,jc),
c     +      jdahep(1,jc),jdahep(2,jc)
c  10        FORMAT(6(1x,i6))
c	    write(77,20) (phep(jz,jc),jz = 1,4),(vhep(jz,jc),jz = 1,4)
c  20        FORMAT(4(1x,f10.5),4(1x,f10.5))
c            enddo
c          NEVHEP=IEV
C...One can copy to some predefined block size
          NEVHEPL=NEVHEP
          NHEPL=NHEP
           DO 500 J=1,NHEP
          ISTHEPL(J)=ISTHEP(J)
          IDHEPL(J)=IDHEP(J)
          JMOHEPL(1,J)=JMOHEP(1,J)
          JMOHEPL(2,J)=JMOHEP(2,J)
          JDAHEPL(1,J)=JDAHEP(1,J)
          JDAHEPL(2,J)=JDAHEP(2,J)
          PHEPL(1,J)=PHEP(1,J)
          PHEPL(2,J)=PHEP(2,J)
          PHEPL(3,J)=PHEP(3,J)
          PHEPL(4,J)=PHEP(4,J)
          PHEPL(5,J)=PHEP(5,J)
          VHEPL(1,J)=VHEP(1,J)
          VHEPL(2,J)=VHEP(2,J)
          VHEPL(3,J)=VHEP(3,J)
          VHEPL(4,J)=VHEP(4,J)
  500     CONTINUE
          HEPMC3STATUS=hepmc3_convert_event(OUTID(ICA))
C Note there should be PDF ids
          HEPMC3STATUS=hepmc3_write_event(OUTID(ICA))
          HEPMC3STATUS=hepmc3_clear_event(OUTID(ICA))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                     endif
                  endif
               endif
            endif
         endif

         if (NTVEC_flag) then
            nt_Nisr = Nisr
            call hfnt(NTID+20)
            do 1100 i=1,Nextn      
               vec_px = vec(1,i)
               vec_py = vec(2,i)
               vec_pz = vec(3,i)
               vec_e  = vec(4,i)
               vec_m  = amass1(i) 
               vec_kf = kfcode(i)
               call hfnt(NTID+30)
               call hfnt(NTID+40)
 1100       continue
            do 1200 i=1,Nisr
               nt_vec_isr(1) = vec_isr(1)
               nt_vec_isr(2) = vec_isr(2)
               nt_vec_isr(3) = vec_isr(3)
               nt_vec_isr(4) = vec_isr(4)
               call hfnt(NTID+50)
 1200       continue
         endif 
      endif   
********* Termination of HBOOK *********
      if (Ievt .GE. Ngen) then
        call HB_term(NTPYT_flag, NTVEC_flag)
      endif
****************************************
      return
      end
* ==========================================================================
      subroutine  HB_term(NTPYT_flag, NTVEC_flag)
c      implicit NONE
include "Pythia6ToHepMC3.inc"
      INTEGER OUTID, HEPMC3STATUS
      common/moj/ OUTID(2), HEPMC3STATUS

* ---------- Argument ----------
      logical  NTPYT_flag, NTVEC_flag
* ------------------------------
* ---------- Ntuple varialbes ----------
      include 'nt_inc.h'
* --------------------------------------
* ------ Local variables ------
      integer  icycle
* -----------------------------
C      include './inc/graepia.h'
      call hldir('//grp', 'T')   
                                 
      call hcdir('//grp', ' ')   
      if (NTPYT_flag) then
        call hrout(NTID, icycle, ' ')
      endif
      call hrout(NTID+10, icycle, ' ')
      if (NTVEC_flag) then
        call hrout(NTID+20, icycle, ' ')
        call hrout(NTID+30, icycle, ' ')
        call hrout(NTID+40, icycle, ' ')
        call hrout(NTID+50, icycle, ' ')
      endif
      call hrend('grp')          
C...Delete output writers
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      HEPMC3STATUS=hepmc3_delete_writer(OUTID(1))
c      HEPMC3STATUS=hepmc3_delete_writer(OUTID(2))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-----------------------------------------------------------------
    
      return
      end
* ==========================================================================
