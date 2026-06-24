      SUBROUTINE IP1(X0,X1,X,Y0,Y1,Y)
C
C
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(X0.EQ.X1)GOTO 10
      F10=(Y1-Y0)/(X1-X0)
      Y=Y0+F10*(X-X0)
      RETURN
 10   Y=(Y0+Y1)/2.
      RETURN
      END