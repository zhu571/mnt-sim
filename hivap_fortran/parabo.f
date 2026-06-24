      SUBROUTINE PARABO(X,Y,A,X0,Y0)
C
C  Sucht Ort und Wert des Maximums einer Parabel
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(3),Y(3),A(3)
      X0=0.
      Y0=0.
      DUM=X(2)-X(1)
      IF(DUM.EQ.0.)RETURN
      F21=(Y(2)-Y(1))/DUM
      DUM=X(3)-X(2)
      IF(DUM.EQ.0.)RETURN
      F32=(Y(3)-Y(2))/DUM
      DUM=X(3)-X(1)
      IF(DUM.EQ.0.)RETURN
      F321=(F32-F21)/DUM
      A(2)=F21-F321*(X(1)+X(2))
      A(3)=F321
      A(1)=Y(1)-X(1)*(F21-X(2)*F321)
      IF(A(3).EQ.0.)RETURN
      X0=-0.5*A(2)/A(3)
      Y0=-0.25*A(2)*A(2)/A(3) + A(1)
      RETURN
      END