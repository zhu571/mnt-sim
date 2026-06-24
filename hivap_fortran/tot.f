      SUBROUTINE TOT(N,M,POP,SPECJ,AMX,AVGJ,TOTAL,SPECE,AVGE,EXMAX)
C
C  Bestimmt Summe und erstes Moment der E-J-Bevoelkerung
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION POP(1),SPECJ(1),SPECE(1)
      JDIM=100
      IEDIM=400
C----------------------------------------------------------
      SUM=0.
      AMX=0.
      INDEX=0
      AVG=0.
      AVGE=0.
      DO 5 I=1,JDIM
 5    SPECJ(I)=0.
      DO 6 I=1,IEDIM
 6    SPECE(I)=0.
C----------------------------------------------------------------
      NN=MIN0(JDIM,N)
      DO 50 JF1=1,NN
      SUMJ=0.
      FJ=JF1-1
      DO 45 KE=1,M
      JFKE=INDEX+KE
      DUM=POP(JFKE)
      SUMJ=SUMJ+DUM
      AMX=DMAX1(DUM,AMX)
 45   CONTINUE
      SPECJ(JF1)=SUMJ
      SUM=SUM+SUMJ
      AVG=AVG+SUMJ*FJ
      INDEX=INDEX+M
 50   CONTINUE
C
      MM=MIN0(IEDIM,M)
      DO 80 KE=1,MM
      FKE=1-KE
      SUME=0.
      JFKE=KE-M
      DO 70 JF1=1,N
      JFKE=JFKE+M
      SUME=SUME+POP(JFKE)
 70   CONTINUE
      SPECE(MM+1-KE)=SUME
      AVGE=AVGE+FKE*SUME
 80   CONTINUE
C
      TOTAL=SUM
      IF(SUM.NE.0.)AVGJ=AVG/SUM
      IF(SUM.NE.0.)AVGE=AVGE/SUM +EXMAX
      RETURN
      END