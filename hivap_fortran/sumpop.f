      SUBROUTINE SUMPOP(N,M,EXMAX,POP,TOT,ISIZE,N1,M1,EXMAX1,POP1,TOT1,
     1                  POP0,FLOST,IPRINT)
C  Addiert Bevoelkerungswahrscheinlichkeiten aus verschiedenen
C  Mutterkernen, die ueber xn, pxn, alphaxn in den gleichen Restkern
C  (Tochterkern) zerfallen
C  Stand: 3.5.1994
C
C     ADDS POP AND POP1 INTO POP0,NEW DIMENSIONS INTO N1,M1,EXMAX1
C                   LOST PARTS INTO FLOST
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION POP(1),POP1(1),POP0(1),SUME(120),SUMJ(400)
C
      GLOST=0.
      IF(TOT .LE.0. .OR. M .EQ.0) GOTO 300
      IF(TOT1.LE.0. .OR. M1.EQ.0) GOTO 400
C
C                            GET ENCOMPASSING DIMENSIONS M0,N0,EXMAX0
      N0=MAX0(N,N1)
      EXMAX0=DMAX1(EXMAX,EXMAX1)
      EXMIN=EXMAX-M+1
      EXMIN1=EXMAX1-M1+1
      EXMIN0=DMIN1(EXMIN,EXMIN1)
      M0=EXMAX0-EXMIN0+1.
      MN0=M0*N0
      MN=M*N
      MN1=M1*N1
      GLOST=0.
      IF(MN0.LE.ISIZE) GOTO 150
C                                  GET SUME(J) OF (POP+POP1)
      JM=-M
      JM1=-M1
      DO 40 J=1,N0
      SUM=0.
      IF(J.GT.N) GOTO 25
      JM=JM+M
      DO 20 K=1,M
 20   SUM=SUM+POP(JM+K)
 25   IF(J.GT.N1) GOTO 35
      JM1=JM1+M1
      DO 30 K=1,M1
 30   SUM=SUM+POP1(JM1+K)
 35   SUME(J)=SUM
 40   CONTINUE
C                                GET SUMJ(KE) OF (POP+POP1)
      IK=EXMAX0-EXMAX   +0.5
      IK1=EXMAX0-EXMAX1 +0.5
      DO 70 K=1,M0
      SUM=0.
      JM=K-IK
      IF(JM.LE.0) GOTO 55
      IF(JM.GT.M) GOTO 55
      DO 50 J=1,N
      SUM=SUM+POP(JM)
 50   JM=JM+M
 55   JM1=K-IK1
      IF(JM1.LE.0) GOTO 65
      IF(JM1.GT.M1) GOTO 65
      DO 60 J=1,N1
      SUM=SUM+POP1(JM1)
 60   JM1=JM1+M1
 65   SUMJ(K)=SUM
 70   CONTINUE
C
      I=1
 80   SUM=DMIN1(SUME(N0),SUMJ(I),SUMJ(M0))
 85   IF(SUME(N0).GT.SUM) GOTO 90
      N0=N0-1
      MN0=M0*N0
      IF(MN0.LE.ISIZE) GOTO 150
      GOTO 80
 90   IF(SUMJ(M0).GT.SUM) GOTO 100
      M0=M0-1
      MN0=M0*N0
      IF(MN0.LE.ISIZE) GOTO 150
      GOTO 80
 100  IF(SUMJ(I).GT.SUM) GOTO 85
      DO 105 K=2,M0
 105  SUMJ(K-1)=SUMJ(K)
      M0=M0-1
      EXMAX0=EXMAX0-1
      MN0=M0*N0
      IF(MN0.LE.ISIZE) GOTO 150
      GOTO 80
C---------------------------------------------------------------------
C                                     SUM POP+POP1 INTO POP0
 150  M0=MAX0(M0,1)
      N0=MAX0(N0,1)
      DUM=EXMAX0-EXMAX+0.5
      IK=DUM
      IK=-IK
      DUM=EXMAX0-EXMAX1+0.5
      IK1=DUM
      IK1=-IK1
      SUM=0.
      JM0=0
      JM1=0
      JM=0
      DO 200 J0=1,N0
      J=J0-1
      DO 180 K0=1,M0
      DUM1=0.
      K1=K0+IK1
      IF(K1.LE.0) GOTO 170
      IF(K1.GT.M1)GOTO 170
      JK1=JM1+K1
      IF(JK1.GT.MN1  ) GOTO 170
      DUM1=POP1(JK1)
 170  DUM=0.
      K=K0+IK
      IF(K.LE.0) GOTO 175
      IF(K.GT.M) GOTO 175
      JK=JM+K
      IF(JK.GT.MN   ) GOTO 175
      DUM=POP(JK)
 175  JK0=JM0+K0
      IF(JK0.GT.ISIZE)GOTO 180
      POP0(JK0)=DUM+DUM1
      SUM=SUM+POP0(JK0)
 180  CONTINUE
      JM0=JM0+M0
      JM1=JM1+M1
      JM=JM+M
 200  CONTINUE
C-------------------------------------------------------------------
      GLOST=TOT1+TOT-SUM
      IF(IPRINT.LT.3)
     1WRITE(6,210) GLOST,SUM,EXMAX0,M0,N0,TOT,EXMAX,M,N,TOT1,EXMAX1,M1,
     2             N1,IK,IK1
 210  FORMAT(' SUMPOP FLOST,SUM,EXMAX,M,N',E12.4,3(E12.4,F8.2,2I4),2I3)
      TOT1=SUM
      N1=N0
      M1=M0
      EXMAX1=EXMAX0
      TOT1=SUM
 215  I1=M1*N1+1
      DO 220 I=I1,ISIZE
 220  POP0(I)=0.
      FLOST=FLOST+GLOST
      RETURN
C
 300  MN1=M1*N1
      DO 310 I=1,MN1
 310  POP0(I)=POP1(I)
      IF(IPRINT.LT.3)
     1WRITE(6,210) GLOST,TOT1,EXMAX1,M1,N1
      GOTO 215
C
 400  MN=M*N
      DO 410  I=1,MN
 410  POP0(I)=POP(I)
      N1=N
      M1=M
      EXMAX1=EXMAX
      TOT1=TOT
      IF(IPRINT.LT.3)
     1WRITE(6,210) GLOST,TOT1,EXMAX1,M1,N1
      GOTO 215
      END