      SUBROUTINE CUTOFF(N1,M1,POP,SPECJ,TOTAL,NSTEP,AMX,CUT,PRCNT,
     1                  ICUT,TRIM,IPOPS,EXMAX1,ALEVEL,IPRINT)
C
C  Eliminiert vernachlaessigbar kleine Elemente in der E-J-Bevoelkerung
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION POP(1),SPECJ(1),TOTAL(1)
C
C------------------------------------------------------------------
      TEST=TOTAL(NSTEP)*CUT
      TRIMJ=0.
      TRIME=0.
      TRIMEL=0.
C------------------------------------------------------------------
C                                        TRIM OFF J STRIPE
  200 SUM=0.
      DO 210 J1=1,N1
      INDEX=N1+1-J1
      SUM=SUM+SPECJ(INDEX)
      IF(SUM-TEST)210,220,220
  210 CONTINUE
  220 N11=N1-J1+1
      IF(N1.EQ.N11)GOTO 235
      K1=N11*M1+1
      KLAST=N1*M1
      DO 230  K=K1,KLAST
  230 POP(K)=0.
      N1=N11
      TRIMJ=SUM-SPECJ(N1)
C-------------------------------------------------------------------
C                                        TRIM OFF E STRIPE
C
  235 SUM=0.
      DO 260 KE=1,M1
      SUM1=0.
      INDEX=KE-M1
      DO 250  J1=1,N1
      INDEX=INDEX+M1
  250 SUM1=SUM1+POP(INDEX)
      SUM=SUM+SUM1
      IF(SUM.GE.TEST) GO TO 270
  260 CONTINUE
  270 M11=M1-KE+1
      M11=MAX0(M11,1)
      IF(M11.EQ.M1) GO TO 300
      TRIME=    SUM-SUM1
C----------------------  REARRANGE POP FOR NEW VALUE OF M1 ----------
      I=KE
      IROW=-M11
      ITROW=-M1
      DO 290 J1=1,N1
      IROW= IROW+M11
      ITROW= ITROW+M1
      DO 285 KE=1,M11
      INDEX=IROW+KE
      ITM=ITROW+I-1+KE
      POP(INDEX)=POP(ITM)
  285 CONTINUE
  290 CONTINUE
      EXMAX1=EXMAX1-FLOAT(I-1)
      M1=M11
C---------------------------------------------------------------------
C                                TRIM OFF INDIVIDUAL SMALL ELEMENTS
C
 300  ALEVEL=0.
      IF(AMX.LE.0.) GOTO 350
      IAGAIN=0
      TRIMEL=0.
      IF(PRCNT.GE.1.) GOTO 350
      IF(PRCNT.LE..5) GOTO 350
      ALEVEL=AMX
      FRACTN=PRCNT*(TOTAL(NSTEP)-TRIME-TRIMJ)
C
 310  ALEVEL=ALEVEL*0.5
      IAGAIN=IAGAIN+1
      INDEX=-M1
      SUM=0.
C
      DO 330 J1=1,N1
      INDEX=INDEX+M1
      DO 320 KE=1,M1
      DUM=POP(INDEX+KE)
      IF(DUM.GT.ALEVEL)SUM=SUM+DUM
 320  CONTINUE
 330  CONTINUE
C
      IF(IAGAIN.GT.20) GOTO 350
      IF(SUM.LT.FRACTN) GOTO 310
C
      ICUT=0
      IF(M1.LT.2) GOTO 350
      INDEX=-M1
      DO 340 J1=1,N1
      INDEX=INDEX+M1
      DO 335 KE=1,M1
      I=INDEX+KE
      IF(POP(I).GT.ALEVEL) GOTO 335
      IF(POP(I).GT.0) ICUT=ICUT+1
      TRIMEL=TRIMEL+POP(I)
      POP(I)=0.
 335  CONTINUE
 340  CONTINUE
 350  TRIM=TRIM+TRIMJ+TRIME+TRIMEL
      IF(IPRINT.LT.3)
     1WRITE(6,351) NSTEP,TOTAL(NSTEP),TRIM,TRIMJ,TRIME,TRIMEL,ICUT,AMX,
     2             ALEVEL,IAGAIN ,CUT,PRCNT,M1,N1
 351  FORMAT(' CUTOFF   NSTEP,TOTAL,TRIM,TRIMJ,TRIME,TRIMLEV',
     1 I4,5E12.4/10X,'ICUT,AMX,ALEVEL,IAGAIN',I4,2E12.4,I4,
     2 4X,'CUT,PRCNT',2E11.3,'  M1,N1=',2I4)
C--------------------------------------------------------------------
C                     ADJUST M1 SO IT EXTENDS AS FAR DOWN AS POSSIBLE
C                           REARRANGE POP FOR NEW VALUES OF M1 AND N1
C
  400 IA=EXMAX1+1.
      M11=MIN0(IA,IPOPS/N1)
      M11=MAX0(M11,1)
      IF(M1.EQ.M11)GOTO 500
      IF(M1.GT.M11)GOTO 440
C                                               M11>M1
      IDIFF=M11-M1
      MN=M11*N1
      MN1=M1*N1
 410  J=MN+1
      K=MN1+1
      DO 420 I=1,IDIFF
      J=J-1
 420  POP(J)=0.
      IF(J.LE.(M1+1)) GOTO 430
      DO 425 I=1,M1
      J=J-1
      K=K-1
 425  POP(J)=POP(K)
      MN=MN-M11
      MN1=MN1-M1
      GOTO 410
 430  IF(IPRINT.LT.3)
     1WRITE(6,435) M1,M11
 435  FORMAT(' CUTOFF OLD,NEW M1',2I4)
      M1=M11
      RETURN
C                                                M1>M11
 440  IDIFF=M1-M11
      SUM=0.
      J=0
      K=-M11
      MN=M11*N1
      MN1=M1*N1
 450  J=J+M11
      K=K+M11
      DO 460 I=1,IDIFF
      J=J+1
      IF(J.GT.IPOPS) GOTO 460
      SUM=SUM+POP(J)
 460  CONTINUE
      IF(J.GE.MN1) GOTO 480
      DO 470 I=1,M11
      J=J+1
      K=K+1
      IF(J.GT.IPOPS) GOTO 470
      POP(K)=POP(J)
 470  CONTINUE
      GOTO 450
 480  M1=M11
      TRIM=TRIM+SUM
 500  CONTINUE
      RETURN
      END