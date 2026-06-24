C        KC02.GRO.FORT(POPUL)
      SUBROUTINE POPUL(JPEP,IEP,NSTEP,JI)
C
C  Wird von 'EVA' aufgerufen ( innnerhalb EP,JP - Schleife )
C  Stand: 3.5.1994
C
       INCLUDE 'common.f'
C-----------------------------------------------------------------------
C
      RJI=JI
      JP1=JI+1
      E=IEP-1
      EXCURR=EXMAX1-E
C
      IF(P(1).GE.-0.5) GOTO 40
      FLAND(JP1)=FLAND(JP1)+POP(JPEP)
      PALL=0.
      GOTO 45
 40   PALL=P(2)+P(3)+P(4)+P(5)
      PER0=P(1)+P(2)+P(3)+P(4)
      IF(P(1).GT.0.) GOTO 50
C
 45   IF(NOG.EQ.0) GOTO 47
      IF(PALL.GT.0.)GOTO 50
 47   TABLE(JPEP)=SYMB(6)
C                           SYMB(6) *   SYMB(7) X
      IF((EXCURR+.1).LT.EJAY(JI+1)) GOTO 48
C     IF(P(5).GT.0.) GOTO 48
      SIG(NSTEP)=SIG(NSTEP)+POP(JPEP)
      IF(JP1.LT.3) GOTO 46
      JP11=JP1+2
      DO 43 JPP=3,JP1,2
      JP11=JP11-2
      IEG=EJAY(JP11)-EJAY(JP11-2)+1.
      IF(IEG.LT.11) SPECQ(IEG)=SPECQ(IEG)+POP(JPEP)
 43   CONTINUE
 46   IF(JI.GT.ISOJ(ISO)) SIGISO(ISO)=SIGISO(ISO)+POP(JPEP)
      IF(PALL.EQ.0.)RETURN
 48   RP(NSTEP)=RP(NSTEP)+POP(JPEP)
      TABLE(JPEP)=SYMB(7)
      PENTRY(JPEP)=0.
      IF(P(5).GT.0.) GOTO 67
      RETURN
C---------------------------------------------------------------
 50   PMAX=0.
      DO 60 K=1,5
 60   PMAX=DMAX1(PMAX,P(K))
      DO 65 K=1,5
      IF(PMAX.GT.P(K)) GOTO 65
      TABLE(JPEP)=SYMB(K)
 65   CONTINUE
      IF(PALL.GT.0.) PENTRY(JPEP) = PENTRY(JPEP) * P(1)/(PALL+P(1))
 67   PALL=PALL+P(1)
      FACT=POP(JPEP)/PALL
C--------------------------------------------------------------------
C                                         GAMMAS
 100  IF(HIST(1).LE.0.) GOTO 200
      FACTOR=FACT*P(1)/HIST(1)
      JFMAX1=MIN0(NWG,JI+LGD+1)
      JF1UP=JFMAX1
      JFMIN1=MAX0(1,JI+1-LGD)
      IF(JFMIN1.GT.JF1UP) GOTO 130
      MUP=M1-IEP
      IF(MUP.EQ.0  .AND.  TABLE(JPEP).EQ.SYMB(1)) GOTO 47
      KEUP=MIN0(MAXG,MUP)
      IF(KEUP.LT.1) GOTO 130
      IPOP=(JFMIN1-1)*M1+IEP
      ITM=0
C
      DO 120 JF1=JFMIN1,JF1UP
      DO 110 KE=1,KEUP
      IPO =IPOP+KE
      IT =ITM+KE
      POP (IPO )=POP (IPO )+TMTBG(IT )*FACTOR
      TMTBG(IT )=0.
 110  CONTINUE
      ITM=ITM+MAXG
      IF(JF1.LT.N1) IPOP=IPOP+M1
 120  CONTINUE
C
 130  SUM=0.
      DO 140 ITM=1,ITG
 140  SUM=SUM+TMTBG(ITM)
      RLOSTG(NSTEP)=RLOSTG(NSTEP)+SUM*FACTOR
      IE=DEL(1)
C
      DO 150 KE=1,MAXG
      IE=IE+1
      SPECG(IE)=SPECG(IE)+ TSPECG(KE)*FACT
      SPECQ(IE)=SPECQ(IE)+ TSPECQ(KE)*FACT
 150  CONTINUE
C--------------------------------------------------------------------
C                                         NEUTRONS
 200  IF(HIST(2).LE.0.) GOTO 300
      FACTOR=FACT*P(2)/HIST(2)
      LJF=1
      ITWO=0
      IF(JFJI.EQ.1) GOTO 204
      JFMAX1=MIN0(NWN,JI+LN )
      JF1UP=MIN0(N2,JFMAX1)
      JFMIN1=MAX0(1,JI+2-LN)
      IF(JFMINN.LT.0)LJF=2
      IF(LJF.EQ.2) ITWO=1
      IDUM=IABS(JFMINN)
      IF(IDUM.GT.JFMIN1) ITWO=-1
      JFMIN1=IDUM
      IF(JFMIN1.GT.JF1UP) GOTO 230
      GOTO 206
 204  IDUM=DJJ2
      FDUM=IDUM
      FDUM=DJJ2-FDUM
C     IF(FDUM.GT.0.25 .AND. FDUM.LT.0.75) ITWO=1
      JF1=RJI-DJJ2+1.5
      JFMIN1=MAX0(1,JF1)
      JF1UP=JFMIN1
 206  MUP=M2-IEP+1
      KEUP=MIN0(MAXN,MUP)
      IF(KEUP.LT.1) GOTO 230
      IPOP=(JFMIN1-1)*M2 +IEP-1
      ITM=0
      IDUM=JI+1-LN
      MM2=M2
      IF(LJF.EQ.2)MM2=M2+M2
C
      DO 220 JF1=JFMIN1,JF1UP,LJF
      JID=JF1-IDUM
      IF(JF1.EQ.N2)ITWO=0
      DO 210 KE=1,KEUP
      IPO =IPOP+KE
      IT =ITM+KE
      DUM=TMTBN(IT)*FACTOR
      POPN(IPO )= POPN(IPO )+DUM
      DJN(JID)=DJN(JID)+DUM
      TMTBN(IT )=0.
      IF(ITWO) 205,210,207
 207  IPO1=IPO+M2
      JID1=JID+1
      GOTO 208
 205  IPO1=IPO-M2
      JID1=MAX0(1,JID-1)
 208  POPN(IPO1)=POPN(IPO1)+DUM
      DJN(JID1)=DJN(JID1)+DUM
 210  CONTINUE
      ITM=ITM+MAXN
      IPOP=IPOP+MM2
 220  CONTINUE
C
 230  SUM=0.
      IDUM=ITN/LJF
      DO 240 ITM=1,IDUM
 240  SUM=SUM+TMTBN(ITM)
      RLOSTN(NSTEP)=RLOSTN(NSTEP)+SUM*FACTOR*FLOAT(LJF)
C
      FACT1=FACT*FLOAT(LJF)
      DUM=0.
      DO 250 KE=1,MAXN
      FIE= DEL(2)+KE -1
      IE=FIE+1
      DUM=DUM+TSPECN(KE)*FIE
 250  SPECN(IE)=SPECN(IE)+TSPECN(KE)*FACT1
      AVESP(2)=AVESP(2)+DUM*FACT1
      SSP(2)=SSP(2)+HIST(2)*FACT
C--------------------------------------------------------------------
C                                         PROTONS
 300  IF(HIST(3).LE.0.) GOTO 400
      FACTOR=FACT*P(3)/HIST(3)
      LJF=1
      ITWO=0
      IF(JFJI.EQ.1) GOTO 304
      JFMAX1=MIN0(NWP,JI+LP )
      JF1UP=MIN0(N3,JFMAX1)
      JFMIN1=MAX0(1,JI+2-LP)
      IF(JFMINP.LT.0)LJF=2
      IF(LJF.EQ.2)ITWO=1
      IDUM=IABS(JFMINP)
      IF(IDUM.GT.JFMIN1) ITWO=-1
      JFMIN1=IDUM
      IF(JFMIN1.GT.JF1UP) GOTO 330
      GOTO 306
 304  IDUM=DJJ3
      FDUM=IDUM
      FDUM=DJJ3-FDUM
C     IF(FDUM.GT.0.25 .AND. FDUM.LT.0.75) ITWO=1
      JF1=RJI-DJJ3+1.5
      JFMIN1=MAX0(1,JF1)
      JF1UP=JFMIN1
 306  MUP=M3-IEP+1
      KEUP=MIN0(MAXP,MUP)
      IF(KEUP.LT.1) GOTO 330
      IPOP=(JFMIN1-1)*M3+IEP-1
      ITM=0
      IDUM=JI+1-LP
      MM3=M3
      IF(LJF.EQ.2)MM3=M3+M3
C
      DO 320 JF1=JFMIN1,JF1UP,LJF
      JID=JF1-IDUM
      IF(JF1.EQ.N3) ITWO=0
      DO 310 KE=1,KEUP
      IPO =IPOP+KE
      IT =ITM+KE
      DUM=TMTBP(IT)*FACTOR
      POPP(IPO )= POPP(IPO )+DUM
      DJP(JID)=DJP(JID)+DUM
      TMTBP(IT )=0.
      IF(ITWO) 305,310,307
 307  IPO1=IPO+M3
      JID1=JID+1
      GOTO 308
 305  IPO1=IPO-M3
      JID1=MAX0(1,JID-1)
 308  POPP(IPO1)=POPP(IPO1)+DUM
      DJP(JID1)=DJP(JID1)+DUM
 310  CONTINUE
      ITM=ITM+MAXP
      IPOP=IPOP+MM3
 320  CONTINUE
C
 330  SUM=0.
      IDUM=ITP/LJF
      DO 340 ITM=1,IDUM
 340  SUM=SUM+TMTBP(ITM)
      RLOSTP(NSTEP)=RLOSTP(NSTEP)+SUM*FACTOR*FLOAT(LJF)
C
      FACT1=FACT*FLOAT(LJF)
      DUM=0.
      DO 350 KE=1,MAXP
      FIE= DEL(3)+KE-1
      IE=FIE+1
      DUM=DUM+TSPECP(KE)*FIE
 350  SPECP(IE)=SPECP(IE)+TSPECP(KE)*FACT1
      AVESP(3)=AVESP(3)+DUM*FACT1
      SSP(3)=SSP(3)+HIST(3)*FACT
C--------------------------------------------------------------------
C                                         ALPHAS
 400  IF(HIST(4).LE.0.) GOTO 500
      FACTOR=FACT*P(4)/HIST(4)
      LJF=1
      ITWO=0
      IF(JFJI.EQ.1) GOTO 404
      JFMAX1=MIN0(NWA,JI+LA )
      JF1UP=MIN0(N4,JFMAX1)
      JFMIN1=MAX0(1,JI+2-LA)
      IF(JFMINA.LT.0)LJF=2
      IF(LJF.EQ.2)ITWO=1
      IDUM=IABS(JFMINA)
      IF(IDUM.GT.JFMIN1) ITWO=-1
      JFMIN1=IDUM
      IF(JFMIN1.GT.JF1UP) GOTO 430
      GOTO 406
 404  IDUM=DJJ4
      FDUM=IDUM
      FDUM=DJJ4-FDUM
C     IF(FDUM.GT.0.25 .AND. FDUM.LT.0.75) ITWO=1
      JF1=RJI-DJJ4+1.5
      JFMIN1=MAX0(1,JF1)
      JF1UP=JFMIN1
 406  MUP=M4-IEP+1
      KEUP=MIN0(MAXA,MUP)
      IF(KEUP.LT.1) GOTO 430
      IPOP=(JFMIN1-1)*M4+IEP-1
      ITM=0
      IDUM=JI+1-LA
      MM4=M4
      IF(LJF.EQ.2)MM4=M4+M4
C
      DO 420 JF1=JFMIN1,JF1UP,LJF
      JID=JF1-IDUM
      IF(JF1.EQ.N4) ITWO=0
      DO 410 KE=1,KEUP
      IPO =IPOP+KE
      IT =ITM+KE
      DUM=TMTBA(IT)*FACTOR
      POPA(IPO )= POPA(IPO )+DUM
      DJA(JID)=DJA(JID)+DUM
      TMTBA(IT )=0.
      IF(ITWO) 405,410,407
 407  IPO1=IPO+M4
      JID1=JID+1
      GOTO 408
 405  IPO1=IPO-M4
      JID1=MAX0(1,JID-1)
 408  POPA(IPO1)=POPA(IPO1)+DUM
      DJA(JID1)=DJA(JID1)+DUM
 410  CONTINUE
      ITM=ITM+MAXA
      IPOP=IPOP+MM4
 420  CONTINUE
C
 430  SUM=0.
      IDUM=ITA/LJF
      DO 440 ITM=1,IDUM
 440  SUM=SUM+TMTBA(ITM)
      RLOSTA(NSTEP)=RLOSTA(NSTEP)+SUM*FACTOR *FLOAT(LJF)
C
      FACT1=FACT*FLOAT(LJF)
      DUM=0.
      DO 450 KE=1,MAXA
      FIE= DEL(4)+KE -1
      IE=FIE+1
      DUM=DUM+TSPECA(KE)*FIE
 450  SPECA(IE)=SPECA(IE)+TSPECA(KE)*FACT1
      AVESP(4)=AVESP(4)+DUM*FACT1
      SSP(4)=SSP(4)+HIST(4)*FACT
C----------------------------------------------------------------------
C                                      FISSION
 500  IF(      P(5).LE.0.) GOTO 600
      DUM=POP   (JPEP)*   P(5)/PALL
      RLOSTF(NSTEP)=RLOSTF(NSTEP)+DUM
      JI1=JI+1
      YFISJ(JI1)=YFISJ(JI1)+DUM
      YPOPJ(JI1)=YPOPJ(JI1)+POP(JPEP)
      YERJ(JI1)=YERJ(JI1)+     POP(JPEP) * PER0/PALL
      AVGJF(NSTEP)=AVGJF(NSTEP)+DUM*FLOAT(JI)
      AVGEF(NSTEP)=AVGEF(NSTEP)+DUM*(EXMAX1+1-IEP)
C
      DO 550 KE=1,MAXF
      IE=KE+DEL(5)
 550  SPECF(IE)=SPECF(IE)+TSPECF(KE)*FACT
C----------------------------------------------------------------------
 600  CONTINUE
      RETURN
      END