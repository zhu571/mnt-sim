      SUBROUTINE   ROT (NSTEP,NROW,AMASS,ZEE,SHELK,ITEST,NOF)
C
C  Liest bzw. berechnet Yrastlinien (auch am Sattelpunkt),
C  liest Yrastlinien von logischer Einheit 'IYR' wenn IRAT=1;
C  modifiziert fuer Sierk-Barrieren bei schweren Elementen (8.2.1989)
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/YR/EJAY(1500),RATIOS(6),RATIOF(6),FINERT,FINERF,
     &          JDIM,JUPYR,IYR,IRAST,INERF,NOROTF
      COMMON/BARR/BFLDM(25,25),BFS(25,25),BAR0,BARFAC,IBF,IFISRT
      COMMON/DNS/YRSMO(1500),CS,CE,UCRIT,RPAGAP,ENHANC(400),
     1           ENHANS(400),IASYM,IASYMS
      COMMON/FISRT/ROT0,A2MS,RKAPPA,X0,Y0,RKAPA4
      DIMENSION ARRAY(500),SHELK(5)
      CHARACTER*1 AW(72)
      EPS=0.01
      AMPAR1=0.
C----------------------------------------------------------------------
C             'IRAST' NOT ZERO: YRAST FROM DISC LOG UNIT IYR
      IF(IRAST   .EQ.0) GOTO 150
      ITEST=0
      J0=JDIM
      JDIM1=JDIM-1
      IRC=0
      CALL DECODN(ARRAY,AW,IRC,5,72)
      NJ =0
      IF(IRC.GT.0) NJ =ARRAY(1)+EPS
      NJ1=0
      IF(IRC.GT.1) NJ1=ARRAY(2)+EPS
      NJ2=0
      IF(IRC.GT.2) NJ2=ARRAY(3)+EPS
      NJ3=0
      IF(IRC.GT.3) NJ3=ARRAY(4)+EPS
      NJ4=0
      IF(IRC.GT.4) NJ4=ARRAY(5)+EPS
C
      IF(NJ .GT.0) CALL MYRD(EJAY   ,NJ,-1,IYR)
      DO 40 K=NJ ,JDIM1
 40   EJAY(K+1)=9999.
C
      IF(NJ1.GT.0) CALL MYRD(EJAY(J0+1),NJ1,-1,IYR)
      DO 41 K=NJ1,JDIM1
      KK=K+1+J0
 41   EJAY(  KK )=9999.
      J0=J0+JDIM
C
      IF(NJ2.GT.0) CALL MYRD(EJAY(J0+1),NJ2,-1,IYR)
      DO 42 K=NJ2,JDIM1
      KK=K+1+J0
 42   EJAY(  KK )=9999.
      J0=J0+JDIM
C
      IF(NJ3.GT.0) CALL MYRD(EJAY(J0+1),NJ3,-1,IYR)
      DO 43 K=NJ3,JDIM1
      KK=K+1+J0
 43   EJAY(  KK )=9999.
      IF(NOF.EQ.1) GOTO 150
      J0=J0+JDIM
C
      IF(NJ4.GT.0) CALL MYRD(EJAY(J0+1),NJ4,-1,IYR)
      DO 44 K=NJ4,JDIM1
      KK=K+1+J0
 44   EJAY(  KK )=9999.
C-----------------------------------------------------------------------
C                               LIQUID DROP ROTOR (CPS 1974)
 150  IF(FINERT.LE.0.) FINERT=1.
      IF(FINERF.LE.0.) FINERF=1.
      DUM=0.
      IF(JUPYR.LE.0) JUPYR=JDIM-1
      AN=AMASS-ZEE
      EJAY(1)=0.
      YRSMO(1)=0.
      J0F=4*JDIM
C----------------------------------------------------
C     FISSION BARRIER FOR ZERO AM
      BF0=BFS(NSTEP,NROW)
      IF(IRAST.EQ.0) EJAY(J0F+1)=BF0
      SPT0=BFLDM(NSTEP,NROW)
      DM0=SPT0*BARFAC
      DBF0=SPT0 * (BARFAC-1.)
      YRSMO(J0F+1)=DM0
C
C------------------------------------------------------
C                          YRAST OF GAMMA DAUGHTER AND SADDLE POINT
      ALIM=RATIOS(6)+16.
      IRAT =RATIOS(6)+1.01
      DIFF=0.
      SHELK0=SHELK(1)
      JUPYR1=MAX0(JUPYR,31)
      CALL FISROT(AMASS,ZEE,AN,0.D0,DELR,BFCPS,EROT,0.D0)
!----------The following is added by N.Wang 2007-07-02 -------------------
	DM0=BFCPS
!-------------------------------------------------------------------------
      DO 160 J=2,JUPYR1
      AA=J-1
      CALL FISROT(AMASS,ZEE,AN,AA,DELR,SPT,EROT,DUM)
      IF(IFISRT.EQ.2) THEN
         IZEE=ZEE+0.01
         MASS=AMASS+0.01
         IF(IZEE.LE.102) THEN
             CALL BARFIT(IZEE,MASS,J-1,SPT,DELR,ELMAX)
             SPT=SPT+DELR
         ELSE
             SPT=SPT-BFCPS +SPT0
         ENDIF
      ENDIF
C
      DELR=DELR*(AA+1.)/(AA*FINERT)
      SPT=(SPT-SPT0)*(AA+1.)/(AA*FINERF) + SPT0
      IF(NOROTF.EQ.2) THEN
         IF(J.EQ.2) DUMG=DELR*0.5
         DELR=AA*(AA+1.)*DUMG
         IF(J.EQ.2)DUMF=(SPT-BFLDM(NSTEP,NROW))*0.5
         SPT=AA*(AA+1.)*DUMF + BFLDM(NSTEP,NROW)
      ENDIF
      DM=(SPT-DELR)*BARFAC
      DM=DMAX1(DM,AMPAR1)
      YRSMO(J)=DELR
      IF(AA.GT.RATIOS(6).OR. IRAST.NE.0) GOTO 154
      CALL IP2(RATIOS(2),RATIOS(4),RATIOS(6),AA,RATIOS(1),RATIOS(3),
     1 RATIOS(5),EJAY(J))
      IF( J.EQ.IRAT     ) SHELK0=(EJAY(J)-YRSMO(J)+SHELK(1))*DM0/DM
      GOTO 155
 154  IF(IRAST.EQ.0) EJAY(J)=YRSMO(J)-SHELK(1)+SHELK0*DM /DM0
 155  YRSMO(J0F+J)=DELR+DM
      IF(NOROTF.EQ.1) YRSMO(J0F+J)=SPT+DBF0
      IF(IRAST.EQ.0)EJAY(J0F+J)=YRSMO(J0F+J)-SHELK(1)+SHELK(5)*DM/DM0
 160  CONTINUE
C
C-----------------------------------------------------------------
C                                    YRAST OF N,P,A DAUGHTERS
      DUM=0.
      K=1
 161  K=K+1
      GOTO(163,163,166,169,190),K
 163  AF=AMASS-1
      ZF=ZEE
      GOTO 175
 166  ZF=ZEE-1
      GOTO 175
 169  AF=AMASS-4.
      ZF=ZEE-2.
 175  J0=(K-1)*JDIM
      SHELK0=SHELK(K)
      AN=AF-ZF
      MASS=AF+0.01
      IZEE=ZF+0.01
      AA=0.
      CALL FISROT(AF,ZF,AN,AA,DELR,DM0,EROT,0.D0)
      BFCPS=DM0
      IF(IFISRT.EQ.2)
     1CALL BARFIT(IZEE,MASS,0,DM0,DELR,ELMAX)
      EJAY(J0+1)=0.
      YRSMO(J0+1)=0.
      DIFF=0.
C---------------------------------
      DO 185 J=2,JUPYR1
      AA=J-1
      CALL FISROT(AF,ZF,AN,AA,DELR,SPT,EROT,DUM)
      IF(IFISRT.EQ.2) THEN
         IF(IZEE.LE.102) THEN
             CALL BARFIT(IZEE,MASS,J-1,SPT,DELR,ELMAX)
             SPT=SPT+DELR
         ELSE
             SPT=SPT-BFCPS +DM0
         ENDIF
      ENDIF
C
      DELR=DELR*(AA+1.)/(AA*FINERT)
      SPT=(SPT-DM0)*(AA+1.)/(AA*FINERF) + DM0
      IF(NOROTF.EQ.2) THEN
         IF(J.EQ.2) DUMG=DELR*0.5
         DELR=AA*(AA+1.)*DUMG
         IF(J.EQ.2)DUMF=(SPT-DM0)*0.5
         SPT=AA*(AA+1.)*DUMF + DM0
      ENDIF
      DM =SPT-DELR
      DM=DMAX1(AMPAR1,DM)
      YRSMO(J0+J)=DELR
      IF(AA.GT.RATIOS(6).OR. IRAST.NE.0) GOTO 178
      CALL IP2(RATIOS(2),RATIOS(4),RATIOS(6),AA,RATIOS(1),RATIOS(3),
     1 RATIOS(5),EJAY(J0+J))
      IF( J.EQ.IRAT     ) SHELK0=(EJAY(J0+J)-YRSMO(J0+J) +SHELK(K))*
     1                            DM0/DM
      GOTO 185
 178  IF(IRAST.EQ.0)EJAY(J0+J)=YRSMO(J0+J)-SHELK(K)+SHELK0*DM/DM0
 185  CONTINUE
C----------------------------
      GOTO 161
 190  IF(IRAST.NE.0) GOTO 999
      JUP1=JUPYR1+1
      DO 205 K=1,5
      J0=(K-1)*JDIM
      DO 200 J=JUP1,JDIM
      JJ=J0+J
 200  EJAY(JJ)=9999.
 205  CONTINUE
C----------------------------------------------------------------
C                                      TEST OUTPUT
      IF(ITEST.EQ.0) GOTO 999
      ITEST=0
      J1=1
      J11=0
      J21=1
      WRITE(16,210) J1,ITEST,ITEST
 210  FORMAT(' ROT OUTPUT  FORMAT:',3I4/' EJAY,EJAYF')
 215  FORMAT(15F8.2)
 220  FORMAT(6X)
      DO 240 K=1,5,4
      J2=J1+JUPYR1-1
      WRITE(16,230) (EJAY(J),J=J1,J2)
 230  FORMAT(5F10.4)
      WRITE(16,220)
      WRITE(16,221)J11,J21
 221  FORMAT(2I6/)
      J1=J1+4*JDIM
 240  CONTINUE
      J0=4*JDIM
      DO 235 J=1,JUPYR1
 235  ARRAY(J)=EJAY  (J0+J)-EJAY  (J)
      WRITE(16,236)
 236  FORMAT(' BARRIER')
      WRITE(16,230) (ARRAY(J),J=1,JUPYR1)
      WRITE(16,220)
      WRITE(16,221)J11,J21
      WRITE(16,241)
 241  FORMAT(' ROT OUTPUT  YRSMO,YRSMOF')
      J1=1
      DO 250 K=1,5,4
      J2=J1+JUPYR1-1
      WRITE(16,230) (YRSMO(J),J=J1,J2)
      WRITE(16,220)
      WRITE(16,221)J11,J21
 232  FORMAT(1X,10F10.4)
      J1=J1+4*JDIM
 250  CONTINUE
      J0=4*JDIM
      DO 255 J=1,JUPYR1
 255  ARRAY(J)=YRSMO (J0+J)-YRSMO (J)
      WRITE(16,236)
      WRITE(16,230) (ARRAY(J),J=1,JUPYR1)
      WRITE(16,220)
      WRITE(16,221)J11,J21
      WRITE(16,220)
 999  RETURN
      END