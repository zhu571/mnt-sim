      SUBROUTINE IP2(X0,X1,X2,X,F0,F1,F2,F)
C
C  Berechnet F durch quadratische Interpolation
C  X0,X1,X2 in aufsteigender Folge
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
C     SPECIAL CASES
      IF(X0.EQ.X1 .AND. X1.EQ.X2) GOTO 20
      IF(X0.EQ.X1 .AND. X1.NE.X2) GOTO 25
      IF(X0.EQ.X2 .AND. X1.NE.X2) GOTO 25
      IF(X1.EQ.X2 .AND. X1.NE.X0) GOTO 35
      GOTO 50
 20   F=(F0+F1+F2)/3.
      RETURN
 25   CALL IP1(X1,X2,X,F1,F2,F)
      RETURN
 35   CALL IP1(X0,X1,X,F0,F1,F)
      RETURN
C
C     NORMAL CASE  FROM APPL.NUM.METH.(1969)11
C
 50   F10=(F1-F0)/(X1-X0)
      F21=(F2-F1)/(X2-X1)
      F210=(F21-F10)/(X2-X0)
      F=F0+(X-X0)*F10+(X-X0)*(X-X1)*F210
      RETURN
      END