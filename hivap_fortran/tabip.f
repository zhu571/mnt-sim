      FUNCTION TABIP(TABLX,TABLY,NTAB,IDEG,ILOG,X,IPRINT)
C
C
C  Stand: 3.5.1994
C
C     INTERPOLATE FROM A TABLE
C     TABLX( ) X-VALUES OF TABLE,MUST BE IN ASCENDING ORDER
C     TABLY( ) Y-VALUES OF TABLE,LIMITED FOR THE MOMENT TO POSITIVE Y
C     NTAB     SIZE OF TABLE
C     X        ARGUMENT
C     IDEG     0,1   LINEAR IP
C              2,>2  QUADRATIC IP
C     ILOG     >0    LOGARITHMIC (IN Y) IP
C-------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TABLX(NTAB),TABLY(NTAB)
      IF(NTAB.LT.2) GOTO 999
      IDG=IDEG
      IF(NTAB.LT.3) IDG=1
C-------------------------------------------------
      IF(X.LT.TABLX(2)) GOTO 40
      IF(X.GT.TABLX(NTAB-1)) GOTO 50
      GOTO 60
 40   I1=1
      GOTO 120
 50   I1=NTAB-2
      IF(IDG.LE.1) I1=NTAB-1
      GOTO 120
 60   DO 80 K=1,NTAB
      IF(X.LT.TABLX(K)) GOTO 90
 80   CONTINUE
      GOTO 999
 90   I1=K-1
      I1=MIN0(I1,NTAB-1)
      IF(IDG.LE.1) GOTO 120
      IF(X-TABLX(K) .LT. TABLX(K+1)-X) I1=K-1
      I1=MAX0(1,I1)
      I1=MIN0(I1,NTAB-2)
C-------------------------------------------------
 120  I2=I1+1
      I3=I1+2
      X1=TABLX(I1)
      Y1=TABLY(I1)
      X2=TABLX(I2)
      Y2=TABLY(I2)
      IF(IDG.GE.2) X3=TABLX(I3)
      IF(IDG.GE.2) Y3=TABLY(I3)
      ILOG1=0
      IF(ILOG.LE.0) GOTO 125
      IF(Y1.EQ.0.) GOTO 125
      Y1=DLOG(ABS(Y1))
      Y2=DLOG(ABS(Y2))
      IF(IDG.GE.2) Y3=DLOG(ABS(Y3))
      ILOG1=1
 125  IF(IDG.LE.1) CALL IP1(X1,X2,X,Y1,Y2,Y)
      IF(IDG.GE.2) CALL IP2(X1,X2,X3,X,Y1,Y2,Y3,Y)
      IF(IPRINT.GT.0) WRITE(6,130) X1,X2,X3,X,Y1,Y2,Y3,Y
 130  FORMAT(' X1',4F8.2,4E11.3)
      TABIP=Y
      IF(ILOG1.GT.0) TABIP=EXP(Y)
      RETURN
 999  TABIP=0.
      RETURN
      END