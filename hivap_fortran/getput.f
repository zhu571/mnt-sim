      SUBROUTINE GETPUT(POP,IPOPS,IUNIT,IREAD,M,N,EXMAX,MZ,MA,ACN,ZCN)
C
C  Speichert oder liest die Puffer der E-J - Bevoelkerung
C  (logische Einheiten 10-13)
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION POP(4096)
      L=IPOPS/256
      I1=1
      IF(IREAD.EQ.1) GOTO 80
      WRITE(IUNIT,101) M,N,EXMAX,MZ,MA,ACN,ZCN,IUNIT
      IF(M.EQ.0) RETURN
      DO 50 NREC=1,L
      I2=I1+255
      DO 79 I=I1,I2
      WRITE(IUNIT,102) POP(I)
 79   CONTINUE
      I1=I1+256
 50   CONTINUE
      RETURN
 80   READ (IUNIT,101) M,N,EXMAX,MZ,MA,ACN,ZCN,IDUM
      IF(M.EQ.0) RETURN
      DO 90 NREC=1,L
      I2=I1+255
      DO 81 I=I1,I2
      READ(IUNIT,102) POP(I)
 81   CONTINUE
      I1=I1+256
 90   CONTINUE
      RETURN
 101  FORMAT(2I4,F8.2,2I4,2F7.3,I4)
 102  FORMAT(E10.5)
      END