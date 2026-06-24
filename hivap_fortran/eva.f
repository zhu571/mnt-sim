      SUBROUTINE EVA(NROW,NSTEP,CUTPOP,NOLEP,NOLJI,NOLJF,LEVPRN,
     1               LDBM,KHS)
C
C  Berechnet ED-JD - Bevoelkerung der Tochter aus Mutter - ED-JD -
C  Bevoelkerung, ruft POPUL auf
C  Option zur Behinderung der Spaltung
C  Stand: 2.5.1994
C
      INCLUDE 'common.f'
      COMMON/DELFIS/EDELFIS,ADELFIS,INOF
      CHARACTER*1 DOT/'.'/, OU/'O'/
      INTEGER*2 NULL
      INTEGER*2 TABL(4096)
C-----------------------------------------------------------------------
      BELOW=0.
      NULL=0
      DO 120 I=1,IPOPS
 120  TABLE(I)=BLANK
      HCST=41.35
C     PLANCKS CONSTANT (MEVSEC*10E22)
      DO 150 K=1,KEG
      IF(EFIVE(KEG-K).GT.0.) GOTO 151
 150  CONTINUE
      KE2MX=0
      GOTO 152
 151  KE2MX=K
C--------------------------------
 152  IF(LPRINT.LT.3)
     1WRITE(6,100) M1,N1,NOP,NOA,NOF,NOG,KE2MX,NOLEP,NOLJI,NOLJF,JFJI
 100  FORMAT(' EVA M1,N1=',2I6/' NOP,NOA,NOF,NOG,KE2MX=',5I2,
     1       '  NOLEP,NOLJI,NOLJF,JFJI=',4I2)
C--------------------------------
      EWG=FLOAT(MWG)+EXMING
      EWN=FLOAT(MWN)+EXMINN
      EWP=FLOAT(MWP)+EXMINP
      EWA=FLOAT(MWA)+EXMINA
      EWF=FLOAT(MWF)+EXMINF
      MG=EXMAX1+2.-DEL(1)
      MN=EXMAX1+Q(2)+2.-DEL(2)
      MP=EXMAX1+Q(3)+2 -DEL(3)
      MA=EXMAX1+Q(4)+2.-DEL(4)
      MF=EXMAX1+2.-DEL(5)
      ICU=0
      IOFF=0
      CUTS=0.
C--------------------
      IF(JOSC1.NE.1)JOSC1=0
      IF(JOSC .NE.1)JOSC =0
      IF(IADD1.NE.1)IADD1=0
      IF(IADD2.NE.1)IADD2=0
      IF(IADD3.NE.1)IADD3=0
C---------------------------------LEP LJI
      LEP=1
      IF(NOLEP.EQ.0 .AND. IAMX.LT.2000) LEP=2
      LJI=1
      IF(NOLJI.EQ.0 .AND. AVGJ(NSTEP).GT.10.)LJI=2
C
      IF(LEP.EQ.1 .AND. LJI.EQ.1) GOTO 195
      IOSCE=MOD(NSTEP,2 )
      IF(LEP.EQ.1) IOSCE=1
      IF(LJI.EQ.1) JOSC=1
C
      DO 180 IEP=1,M1
      IF(LEP.EQ.2) IOSCE=1-IOSCE
      E=IEP-1
      EXCURR=EXMAX1-E
      JPEP=IEP
      IF(IOSCE.EQ.1)JOSC1=1-JOSC1
      JOSC=JOSC1
      DO 170 JP1=1,N1
      IF(LJI.EQ.2)JOSC=1-JOSC
      IEF =EXCURR-EJAY(JP1)
      TABL(JPEP)=1
      IF(((IOSCE.EQ.0) .OR. (JOSC.EQ.0)).AND. (IEF .GT.20))
     1   TABL(JPEP)=NULL
      IF(POP(JPEP).LE.0.) TABL(JPEP)=NULL
      JPEP=JPEP+M1
 170  CONTINUE
 172  FORMAT(1X,120I1)
 180  CONTINUE
C---------------------REPOP
      M1N1=M1*N1
      DO 192 IEP=1,M1
      JPEP=IEP-M1
      DO 190 JP1=1,N1
      JPEP=JPEP+M1
      IF(TABL(JPEP).GT.NULL) GOTO 190
      IF(POP(JPEP).LE.0.)GOTO 190
      ISU2=0
      L=JPEP-M1
      IF(L.GT.0)ISU2=ISU2+TABL(L)
      L=JPEP+M1
      IF(L.LE.M1N1) ISU2=ISU2+TABL(L)
      IF(IEP.EQ.M1) GOTO 1810
      L=JPEP+1
      IF(L.LE.M1N1) ISU2=ISU2+TABL(L)
      L=JPEP-M1+1
      IF(L.GT.0)ISU2=ISU2+TABL(L)
      L=JPEP+M1+1
      IF(L.LE.M1N1)ISU2=ISU2+TABL(L)
 1810 IF(IEP.EQ.1) GOTO 183
      L=JPEP-1
      IF(L.GT.0)ISU2=ISU2+TABL(L)
      L=JPEP-M1-1
      IF(L.GT.0)ISU2=ISU2+TABL(L)
      L=JPEP+M1-1
      IF(L.LE.M1N1)ISU2=ISU2+TABL(L)
C
 183  ISUM=ISU2
      IF(ISUM.EQ.0)GOTO  190
      DUM=POP(JPEP)/FLOAT(ISUM)
      POP(JPEP)=0.
 1830 FORMAT(5I6,E12.4)
C
      L=JPEP-M1
      IF(L.LE.0)GOTO 1840
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM
 1840 L=JPEP+M1
      IF(L.GT.M1N1) GOTO 185
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM
 185  IF(IEP.EQ.M1) GOTO 1841
      L=JPEP+1
      IF(L.GT.M1N1) GOTO 1842
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM
 1842 L=JPEP-M1+1
      IF(L.LE.0) GOTO 1844
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM
 1844 L=JPEP+M1+1
      IF(L.GT.M1N1) GOTO 1841
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM
 1841 IF(IEP.EQ.1) GOTO 190
      L=JPEP-1
      IF(L.LE.0) GOTO 1843
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM
 1843 L=JPEP-M1-1
      IF(L.LE.0)GOTO 1845
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM
 1845 L=JPEP+M1-1
      IF(L.GT.M1N1) GOTO 190
      IF(TABL(L).GT.NULL) POP(L)=POP(L)+DUM
 190  CONTINUE
 192  CONTINUE
 195  NOGG=NOG
C------------------------------
C
C******************************  EP LOOP  ******************************
C
      DO 2000 IEP=1    ,M1
C------------------------------
      E=FLOAT(IEP-1)+0.0001
      EXCURR=EXMAX1-E
      ILAND=NROW+NSTEP+IEP
C------------------------------------------------------------
      NOG=NOGG
      MAXG=MIN0(MG-IEP,KEG)
      IF(MAXG.LT.0)  MAXG=0
      IF(NOG    .EQ.1) MAXG=0
      IGL=MAXG+MOD(MAXG,2)
      IGL=MAX0(IGL,4)
      ITG=5*IGL
C
      MAXN=MIN0(MN-IEP,KEN)
      IF(MAXN.LT.0)  MAXN=0
      IF(NON.EQ.1) MAXN=0
      INL=MAXN+MOD(MAXN+1,2)
      INL=MAX0(INL,5)
      ITN=(2*LN -1)*INL
      IF(JFJI.EQ.1)ITN=INL
C------------------------------------------------------------
      MAXP=MIN0(MP-IEP,KEP)
      IF(MAXP.LT.0)  MAXP=0
      IF(NOP.EQ.1) MAXP=0
      IPL=MAXP+MOD(MAXP,2)
      IPL=MAX0(IPL,4)
      ITP=(2*LP -1)*IPL
      IF(JFJI.EQ.1)ITP=IPL
C------------------------------------------------------------
      MAXA=MIN0(MA-IEP,KEA)
      IF(MAXA.LT.0)  MAXA=0
      IF(NOA     .EQ.1) MAXA=0
      IAL=MAXA+MOD(MAXA,2)
      IAL=MAX0(IAL,4)
      ITA=(2*LA -1)*IAL
      IF(JFJI.EQ.1)ITA=IAL
C-------------------------------------------------------------
      MAXF=MIN0(MF-IEP,KEF)
      IF(EXCURR.LT.BFS(NSTEP,NROW))  MAXF=0
      IF(NOF.EQ.1) MAXF=0
C--------------------------------------------------------------
C
C************ READJUST LEVEL DENSITY ***********************************
C
  200 IF(IEP.EQ.1) GOTO 230
  201 IF(MWG.LT.IDEN/NWG.AND.MWA.LT.IDNA/NWA.AND.MWN.LT.IDEN/NWN.AND.
     1    MWP.LT.IDEN/NWP.AND.MWF.LT.IDEN/NWF) GOTO 230
      IDUM=0
      IF(EXMAXG.LT.EWG)  GOTO 205
      LADJ(1)=1
      IDUM=1
 205  IF(EXMAXN.LT.EWN)  GOTO 210
      LADJ(2)=1
      IDUM=1
 210  IF(EXMAXP.LT.EWP)  GOTO 215
      LADJ(3)=1
      IDUM=1
 215  IF(EXMAXA.LT.EWA)  GOTO 220
      LADJ(4)=1
      IDUM=1
 220  IF(EXMAXF.LT.EWF)  GOTO 225
      LADJ(5)=1
      IDUM=1
 225  IF(IDUM .EQ.0) GO TO 230
      CALL DENSTY(NROW,NSTEP,LDBM,1,LEVPRN,KHS)
C     GO TO 200
C***********************************************************************
  230 JPEPMX=(N1-1)*M1+IEP
C                                    OFFSETS FOR LEVDNS
      KONSTG=EXMAXG-EXMAX1+E+DEL(1)
      KONSTN=EXMAXN-EXMAX1+E-Q(2)+DEL(2)
      KONSTP=EXMAXP-EXMAX1+E-Q(3)+DEL(3)
      KONSTA=EXMAXA-EXMAX1+E-Q(4)+DEL(4)
      KONSTF=EXMAXF-EXMAX1+E+DEL(5)
C
C******************* JP (JI,J) LOOP ************************************
C
      DO 1000 JPEP= IEP,JPEPMX, M1
C--------------------------------------FIRST CHECK CUTOFF-------
      IF(POP(JPEP).GT.0.)GOTO 350
      ICU=ICU+1
      TABLE(JPEP)=BLANK
      GOTO 1000
 350  JI=(JPEP-IEP)/M1
      IF(POP(JPEP).GT.CUTPOP )GOTO 360
 355  IOFF=IOFF+1
      CUTS=CUTS+POP(JPEP)
      PENTRY(JPEP)=0.
      TABLE(JPEP)=DOT
C     TOTAL(NSTEP)=TOTAL(NSTEP)-POP(JPEP)
C     POP(JPEP)=0.
      GOTO 1000
C-----------------------------------------------------------------------
 360  JI1=JI+1
      ETH=EXCURR-EJAY(JI1)
      IF(ETH.GT.0.) GOTO 364
      TABLE(JPEP)=OU
      BELOW=BELOW+POP(JPEP)
      POP(JPEP)=0.
      GOTO 1000
 364  IEF=ETH + 0.5
C     IF(TABL(JPEP).EQ.NULL .AND. IEF.GT.20) GOTO 355
C     ------------------------------ISOMER
      IF(ISO.EQ.0) GOTO 380
      IF(JI.EQ.ISOJ(ISO) .AND. DUM.LT.1.1) SIGISO(ISO)=SIGISO(ISO)
     1                                                 +POP(JPEP)
C     ------------------------------
 380  IF(IOPT.EQ.0) GOTO 400
      IF(DUM                .GT.STRIPE) GOTO 400
      P(1)=-1.
      GOTO 900
C-----------------------------------------------NEUTRONS-------
  400 IF(MAXN.EQ.0) GOTO 490
      IF(EXMAX2.LE.0.) GOTO 490
      DO 404 K=1,ITN
  404 TMTBN (K)=0.
      DO 408 K=1,INL
  408 TSPECN(K)=0.
      IF(JFJI.EQ.1) GOTO 407
      JFMAX1=MIN0(NWN,JI+LN )
      JFMIN1=MAX0(1,JI+2-LN )
      GOTO 409
 407  RJI=JI
      DJJ2=DJ(2)
 401  IF(DJJ2.LE.0.) GOTO 402
      JDUM=10.*DJJ2+0.5
      IF(JI.GT.JDUM) GOTO 403
      DJJ2=DJJ2-1.
      GOTO 401
 402  DJJ2=0.
 403  JFMIN1=RJI-DJJ2  +1.5
      JFMIN1=MAX0(1,JFMIN1)
      JFMAX1=JFMIN1
      LMAX2=LN+1
      LMIN1=LMAX2
 409  LJF=1
      JFMINN=JFMIN1
      MMWN=MWN
      IF(NOLJF.EQ.1 .OR. IEF.LT.20) GOTO 410
      IDUM=JFMAX1-JFMIN1
      IF(IDUM.LT.6) GOTO 410
      LJF=2
      MMWN=MWN+MWN
      IADD1=1-IADD1
      IDUM=JI1-JFMIN1+IADD1
      JFMIN1=JFMIN1+MOD(IDUM,2)
      JFMINN=-JFMIN1
 410  KNDEX=0
      JFMIN = JFMIN1-1
      JROWS=JFMIN*MWN
C
C************************************** JD LOOP  N    ******************
C
      DO 440 JF1=JFMIN1,JFMAX1,LJF
      IF(JFJI.EQ.1)GOTO 420
      LMAX2=JI1+JF1
      LMIN1=IABS(JI1-JF1)+1
C---------------------------------------------
 420  DO 430 IE=1,MAXN
      IW=KONSTN+IE
      IF(IW.GT.MWN) GO TO 435
      IW=IW+JROWS
      IF(OMEGN(IW).LE.0.)   GO TO 435
      B=TCOFN(IE,LMIN1)
      IF(LMAX2.LE.LN ) B=B-TCOFN(IE,LMAX2)
      B=B*OMEGN(IW)
      TSPECN(IE)      =TSPECN(IE)      +B
      INDEX=KNDEX+IE
      TMTBN(INDEX)=B
  430 CONTINUE
C---------------------------------------------
  435 KNDEX=KNDEX+MAXN
      JROWS=JROWS+MMWN
  440 CONTINUE
C
C     IF(NSTEP.GT.1 .OR. JI.GT.7) GOTO 459
C     WRITE(6,450) JFMIN1,JFMAX1,LJF,JFJI,MAXN,KONSTN,MWN,NWN,M2,N2
C450  FORMAT(' EVA',10I6)
C     IT1=1
C     DO 456 JF1=JFMIN1,JFMAX1,LJF
C     IT2=IT1+MAXN-1
C     WRITE(6,455) (TMTBN(I8),I8=IT1,IT2)
C455  FORMAT(10E12.4)
C     IT1=IT1+MAXN
C456  CONTINUE
C459  CONTINUE
C**************************************** END JD LOOP  N   *************
C
C------------------------------------NEUTRON INTEGRATIONS------
      P(2)=0.
      SUM=0.
C     NEUTRON INTEGRATIONS HAVE NO IMPLIED ZEROS
      DO 470 K= 2 ,INL,2
  470 P(2)=P(2)+TSPECN(K)
      HIST(2)=P(2)
      P(2)=P(2)*2.
      DO 480 K= 3 ,INL,2
  480 SUM =SUM +TSPECN(K)
      HIST(2)=HIST(2)+SUM+TSPECN(1)
      P(2)=P(2)+SUM
      P(2)=(2.*P(2)+TSPECN(1))/3.
      P(2)=P(2)*LJF
      HIST(2)=HIST(2)*LJF
      GOTO 500
  490 P(2)=0.
      HIST(2)=0.
C-----------------------------------------------PROTONS---------
  500 IF(MAXP.EQ.0) GOTO 590
      IF(EXMAX3.LE.0.) GOTO 590
      DO 504 K=1,IPL
  504 TSPECP(K)=0.
      DO 508 K=1,ITP
  508 TMTBP (K)=0.
      IF(JFJI.EQ.1) GOTO 507
      JFMAX1=MIN0(NWP,JI+LP )
      JFMIN1=MAX0(1,JI+2-LP )
      GOTO 509
 507  RJI=JI
      DJJ3=DJ(3)
 501  IF(DJJ3.LE.0.) GOTO 502
      JDUM=10.*DJJ3+0.5
      IF(JI.GT.JDUM) GOTO 503
      DJJ3=DJJ3-1.
      GOTO 501
 502  DJJ3=0.
 503  JFMIN1=RJI-DJJ3  +1.5
      JFMIN1=MAX0(1,JFMIN1)
      JFMAX1=JFMIN1
      LMAX2=LP+1
      LMIN1=LMAX2
 509  LJF=1
      JFMINP=JFMIN1
      MMWP=MWP
      IF(NOLJF.EQ.1 .OR. IEF.LT.20) GOTO 510
      IDUM=JFMAX1-JFMIN1
      IF(IDUM.LT.6) GOTO 510
      LJF=2
      MMWP=MWP+MWP
      IADD2=1-IADD2
      IDUM=JI1-JFMIN1+IADD2
      JFMIN1=JFMIN1+MOD(IDUM,2)
      JFMINP=-JFMIN1
 510  KNDEX=0
      JFMIN = JFMIN1-1
      JROWS=JFMIN*MWP
C
C************************************** JD LOOP  P    ******************
C
      DO 540 JF1=JFMIN1,JFMAX1,LJF
      IF(JFJI.EQ.1)GOTO 520
      LMAX2=JI1+JF1
      LMIN1=IABS(JI1-JF1)+1
C-----------------------------------------
 520  DO 530 IE=1,MAXP
      IW=KONSTP+IE
      IF(IW.GT.MWP) GO TO 535
      IW=IW+JROWS
      IF(OMEGP(IW).LE.0.)   GO TO 535
      B=TCOFP(IE,LMIN1)
      IF(LMAX2.LE.LP ) B=B-TCOFP(IE,LMAX2)
      B=B*OMEGP(IW)
      TSPECP(IE)      =TSPECP(IE)      +B
      INDEX=KNDEX+IE
      TMTBP(INDEX)=B
  530 CONTINUE
C-----------------------------------------
  535 KNDEX=KNDEX+MAXP
      JROWS=JROWS+MMWP
  540 CONTINUE
C
C**************************************** END JD LOOP  P   *************
C
C------------------------------------PROTON INTEGRATIONS---------
      P(3)=0.
      SUM=0.
C     ZEROS IMPLIED AT BOTH ENDS OF CHARGED PARTICLE SPECTRA
      DO 570 K= 1 ,IPL,2
      P(3)=P(3)+TSPECP(K)
  570 CONTINUE
      HIST(3)=P(3)
      P(3)=2.*P(3)
      DO 580 K= 2 ,IPL,2
      SUM =SUM +TSPECP(K)
  580 CONTINUE
      HIST(3)=HIST(3)+SUM
      P(3)=(P(3)+SUM)*0.6666667
      P(3)=P(3)*LJF
      HIST(3)=HIST(3)*LJF
C     END PROTON EVAPORATION
      GOTO 600
  590 P(3)=0.
      HIST(3)=0.
C---------------------------------------------ALPHAS---------------
 600  IF(MAXA.EQ.0) GOTO 690
      IF(EXMAX4.LE.0.) GOTO 690
      DO 604 K=1,IAL
 604  TSPECA(K)=0.
      DO 608  K=1,ITA
 608  TMTBA(K)=0.
      IF(JFJI.EQ.1) GOTO 607
      JFMIN1=MAX0(1,JI+2-LA )
      JFMAX1=MIN0(JI+LA ,NWA)
      GOTO 609
 607  RJI=JI
      DJJ4=DJ(4)
 601  IF(DJJ4.LE.0.) GOTO 602
      JDUM=5.*DJJ4+0.5
      IF(JI.GT.JDUM) GOTO 603
      DJJ4=DJJ4-1.
      GOTO 601
 602  DJJ4=0.
 603  JFMIN1=RJI-DJJ4  +1.5
      JFMIN1=MAX0(1,JFMIN1)
      JFMAX1=JFMIN1
      LMAX2=LA+1
      LMIN1=LMAX2
 609  LJF=1
      JFMINA=JFMIN1
      MMWA=MWA
      IF(NOLJF.EQ.1 .OR. IEF.LT.20) GOTO 610
      IDUM=JFMAX1-JFMIN1
      IF(IDUM.LT.6) GOTO 610
      LJF=2
      MMWA=MWA+MWA
      IADD3=1-IADD3
      IDUM=JI1-JFMIN1+IADD3
      JFMIN1=JFMIN1+MOD(IDUM,2)
      JFMINA=-JFMIN1
 610  KNDEX=0
      JFMIN = JFMIN1-1
      JROWS=JFMIN*MWA
C
C************************************* JD LOOP ALPHAS ******************
C
      DO 640 JF1=JFMIN1,JFMAX1,LJF
      IF(JFJI.EQ.1)GOTO 620
      LMAX2=JI1+JF1
      LMIN1=IABS(JI1-JF1)+1
C------------------------------------------------
 620  DO 630 KE=1,MAXA
      IW=KONSTA+KE
      IF(IW.GT.MWA) GO TO 635
      IW=IW+JROWS
      IF(OMEGA(IW).LE.0.)   GO TO 635
      B=TCOFA(KE,LMIN1)
      IF(LMAX2.LE.LA )B=B-TCOFA(KE,LMAX2)
      B=B*OMEGA(IW)
      TSPECA(KE)=TSPECA(KE)+B
      INDEX=KNDEX+KE
      TMTBA(INDEX)=B
  630 CONTINUE
C------------------------------------------------
  635 KNDEX=KNDEX+MAXA
      JROWS=JROWS+MMWA
  640 CONTINUE
C
C**************************************** END JD LOOP ALPHAS ***********
C----------------------------------ALPHA INTEGRATIONS-------
      P(4)=0.
      SUM=0.
      DO 670 KE= 1 ,IAL,2
      P(4)=P(4)+TSPECA(KE)
  670 CONTINUE
      HIST(4)=P(4)
      P(4)=2.*P(4)
      DO 680 KE=2  ,IAL,2
      SUM =SUM +TSPECA(KE)
  680 CONTINUE
      HIST(4)=HIST(4)+SUM
      P(4)=(P(4)+SUM)*0.6666667
      P(4)=P(4)*LJF
      HIST(4)=HIST(4)*LJF
      GOTO 700
  690 P(4)=0.
      HIST(4)=0.
C------------------------------------------ GAMMAS --------------------
C
  700 IF(NOG.EQ.1) GOTO 790
      IF(MAXG.LE.0) GO TO 790
      DO 704 KE=1,IGL
      TSPECQ(KE)=0.
  704 TSPECG(KE)=0.
      DO 708 K=1,ITG
  708 TMTBG(K)=0.
      DUMQ=0.
      LGD=2
      KE2M=KE2MX
C     IF(IEF.GT.20) KE2M=0
      IF(KE2M.EQ.0) LGD=1
      JFMIN=MAX0(JI-LGD,0)
      JFMIN1=JFMIN+1
      JFMAX1=MIN0(JI+1+LGD,NWG)
      GAMJ=1.
      IF(JFACTR.EQ.1)GAMJ=JI+JI+1
      IF(JFMAX1.LT.JFMIN1) GOTO 790
      JROWS=JFMIN*MWG
      KNDEX=0
C************************************** JD LOOP GAMMAS **************
C
      DO 740 JF1=JFMIN1,JFMAX1
      JDG=IABS(JI1-JF1)
C----------------------------------------------
      DO 730 KE=1,MAXG
      IW=KONSTG+KE
      IF(IW.GT.MWG)GOTO 735
      IW=IW+JROWS
      IF(OMEGG(IW).LE.0.)GOTO 735
      DUMD=0.
      IF(JDG.LE.1) DUMD=ECUBE(KE)*OMEGG(IW)*GAMJ
      TSPECG(KE)=TSPECG(KE)+DUMD
      IF(KE.GT.KE2M) GOTO 720
      DUMQ=EFIVE(KE)*OMEGG(IW)
      TSPECQ(KE)=TSPECQ(KE)+DUMQ
      DUMD=DUMD+DUMQ
 720  INDEX=KNDEX+KE
      TMTBG(INDEX)=DUMD
 730  CONTINUE
C----------------------------------------------
 735  KNDEX=KNDEX+MAXG
      JROWS=JROWS+MWG
 740  CONTINUE
C*********************************** END JD LOOP GAMMAS ************
C
C--------------------------- GAMMA INTEGRATIONS-----------
      P(1)=0.
      SUM=0.
      DO 770 KE=1,IGL,2
 770  P(1)=P(1)+TSPECQ(KE)+TSPECG(KE)
      HIST(1)=P(1)
      P(1)=2.*P(1)
      DO 780 KE=2,IGL,2
 780  SUM =SUM +TSPECG(KE)+TSPECQ(KE)
      HIST(1)=HIST(1)+SUM
      P(1)=(P(1)+SUM)*0.6666667
      GO TO 800
 790  P(1)=0.
      HIST(1)=0.
C----------------------------------------------- FISSION-------------
C
 800  IF(MAXF.EQ.0)GOTO 890
      IF(JI.GE.NWF) GOTO 890
      DO 804 KE=1,KEF
 804  TSPECF(KE)=0.
      JROWS=JI*MWF
      B=0.
      DO 830 KE=1,MAXF
      IW=KONSTF+KE
      IF(IW.GT.MWF) GOTO 835
      IW=IW+JROWS
      IF(OMEGF(IW).LE.0.)GOTO 835
      TSPECF(KE)=OMEGF(IW)
 830  CONTINUE
C------------------------------------
 835  SUM=0.
      SUM1=0.
      IF(MAXF.EQ.1) GOTO 885
      DO 870 KE=2,MAXF,2
 870  SUM=SUM+TSPECF(KE)
      IF(MAXF.EQ.2) GOTO 885
      DO 880 KE=3,MAXF,2
 880  SUM1=SUM1+TSPECF(KE)
 885  P(5)=0.333333*TSPECF(1)+1.333333*SUM+0.6666667*SUM1
C
C Behinderung der Spaltung, fuer E* > EDELFIS Spaltung behindert
C
      IF(INOF.NE.1) GOTO 889
      IF(EXMAX1.LT.EDELFIS) GOTO 889
      BEHIF1=0.6931
      EXPFAK=BEHIF1*(EXMAX1-EDELFIS)/ADELFIS
      BEHIFAK=1./EXP(EXPFAK)
      P(5)=P(5)*BEHIFAK
C
 889  HIST(5)=SUM+SUM1+TSPECF(1)
      GOTO 900
 890  P(5)=0.
      HIST(5)=0.
C**********************************************************************
C
 900  IF(M1.LE.0) WRITE (6,125) JPEP,IEP,NSTEP,JI
      CALL POPUL(JPEP,IEP,NSTEP,JI)
      IF(M1.LE.0) WRITE (6,125) JPEP,IEP,NSTEP,JI
 125  FORMAT(' POPUL ',10I10)
C
      IF(ILAND.NE.3 .OR.IOPT.NE.0) GOTO 1000
      PALL=P(1)+P(2)+P(3)+P(4)+P(5)
      IF(PALL.LE.0.) GOTO 1000
      IW=IEP+JI*MWG
      IDUM=IW+DEL(1)+0.5
      IW1=IW
      IF(MWG.GE.IDUM)IW1=IDUM
      IF(OMEGG(IW1).LE.0.)GOTO 1000
      DUM=OMEGG(IW)/OMEGG(IW1)
      FLAND(JI1)=HCST*DUM*OMEGG(IW)/ PALL
 1000 CONTINUE
C
C********************************* END JP LOOP *************************
C
 2000 CONTINUE
C
C********************************* END EP LOOP *************************
C
      IF(RLOSTF(NSTEP).GT.0.) AVGJF(NSTEP)=AVGJF(NSTEP)/RLOSTF(NSTEP)
      IF(RLOSTF(NSTEP).GT.0.) AVGEF(NSTEP)=AVGEF(NSTEP)/RLOSTF(NSTEP)
      IF(SSP(2)       .GT.0.) AVESP(2    )=AVESP(2    )/SSP(2)
      IF(SSP(3)       .GT.0.) AVESP(3    )=AVESP(3    )/SSP(3)
      IF(SSP(4)       .GT.0.) AVESP(4    )=AVESP(4    )/SSP(4)
      EAVN(NSTEP)=AVESP(2)
      EAVP(NSTEP)=AVESP(3)
      EAVA(NSTEP)=AVESP(4)
      MN=M1*N1
      TRIM(NSTEP)=TRIM(NSTEP)+CUTS
      IF(LPRINT.LT.3)
     1WRITE(6,5000) MN,ICU,IOFF,CUTPOP,CUTS,BELOW
 5000 FORMAT(/' EVA  M1*N1',I6,4X,'ZEROS',I6,4X,'CUTS',I6,4X,'CUTPOP',
     1        E12.4,4X,'SUM CUTS',E12.4 /' LOSS BELOW YRAST LINE',E12.4)
      NOG=NOGG
      RETURN
      END