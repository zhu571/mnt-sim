        SUBROUTINE BETA(IZ,IA,DU0,DP0,EI,EDAMP,BETAMX,FROT,UEFF,
     1			AMPAR1,AMPAR2)
C
C  Berechnet wahrscheinlichste Deformation als Funktion der Kerntemperat
C  frei nach Vermeulen 1983
C  Stand: 3.5.1994
C       CALCULATES MOST PROBABLE DEFO BETAMX AS FUNCTION OF EXCIT EI
C       BETA0 GROUND STATE DEFORMATION BETA
C       DU0   GROUND STATE SHELL CORRECTION
C       DP0   PAIRING SHIFT(NOT USED FOR THE MOMENT)
C       EDAMP IGNATYUK DAMPING ENERGY
C       FORTRAN VERSION OF KC52.THEO.PLI(BETA)  SEE COMMENTS THERE
C       PROG. D.VERMEULEN,TECHNICALLY MODIFIED BY WR
C
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 N,IMAX,IMIN,M0
        COMMON/DNS/YRSMO(1500),CS,CE,UCRIT,BETA0,ENHANC(400),
     1             ENHANS(400),IENH,IENHP
C                ONLY BETA0 IS NEEDED IN THIS ROUTINE
        DATA AR0,M0,ECRIT,C3
     1    /.444,938.,5.,0.7714/
C------------------------------------------------------------------
      Z=IZ
      A=IA
      A3=A**.333333
      A6=A3*A3
      N = A - Z
      DP0 = 0
      ALPHAR = 2.23606 * AR0 / A3
      ALPHA0 = 0.63078 * BETA0
      THETA0 = ALPHA0 / ALPHAR
      THETA2=THETA0*THETA0
      R0 = 1.12E0 * A3
      C2 = 18.36E0*(1.E0-1.78E0*((N-Z)/A)**2.)
      ELD = (.4E0 *C2*A6 - .2*C3*Z*Z/A3) * ALPHAR*ALPHAR
      BPHASE=0.
      IF(THETA0.GT.0.1E-3)
     1BPHASE = THETA0/.47E0 - ATAN((2.E0*ELD*THETA0/(DU0-ELD*
     2     THETA2)-0.90*THETA0)*.47E0)
      DUK = (DU0-ELD*THETA2) / (COS(THETA0/.47-BPHASE)*
     1     EXF(-THETA2/2.2))
      RHOMX = -1.E10
      DBETA=0.04
      IBETA1=1
      IBETA2=15
      BET=-0.04
      IREP=0
C-------------------------------------------------------------
 40   DO 200 IBETA=IBETA1,IBETA2
        BET=BET+DBETA
        ALPHA = 0.63078 * BET
        DELTA = 1.5E0 * ALPHA
        THETA = ALPHA / ALPHAR
        THETA2= THETA*THETA
        AELL = R0 * (1.E0 - .27951 * BET)
        CELL = R0 * (1.E0 + .55902 * BET)
        IMAX = .2E0 * A * M0 * (AELL*AELL + CELL*CELL)
        IMIN = .4E0 * A * M0 *  AELL*AELL
        BS = 1.E0 + .4E0 * ALPHA*ALPHA
        BK=BS
        ATIL = A/14.61E0 * (1.E0 + 3.114E0 * BS / A3
     1                            + 5.626E0 * BK / A6)
        GAMIGN = ATIL / (0.4E0 * A*A3)
        IF(EDAMP.GT.0.) GAMIGN=1./EDAMP
        EDAMP=1./GAMIGN
        IF(DELTA.LE.0.) GOTO 50
          TCROT = (40.E0 * DELTA) / A3
          GAMROT = 1.E0 / (ATIL * TCROT*TCROT)
 50     DU = DUK * COS(THETA/.47-BPHASE)*EXF(-THETA2/2.2)
        EDROP = ELD * THETA2
        E = EI + DU0 - DU - EDROP
        IF(E.GE.0.) GOTO 70
        RHO = 0.
        GOTO 90
C
 70       F = 1.E0 - EXF(-GAMIGN * E)
          UEFF=E+DU*F
          UEFF=DMAX1(AMPAR1,UEFF)
          T = SQRT(UEFF/ATIL)
          SIGMX2 = 2.577E-5*T*IMAX
          SIGMIN = 5.076E-3*SQRT(T * IMIN)
          IF(DELTA.LE.0.) GOTO 75
          FROT = 1.E0 + SIGMX2 * EXF(-GAMROT * UEFF)
          GOTO 77
 75       FROT = 1.E0
 77       IF(E.LT.ECRIT)
     1      H = 1.E0 - (1.E0 - E/ECRIT)**2.
          IF(E.GE.ECRIT) H = 1.E0
          S = 2.E0 * SQRT(ATIL * (UEFF + DP0 * H) )
          DET =  SQRT(ATIL) * UEFF**2.5
          RHO = EXP(S) / (  SQRT(DET) ) * FROT /
     1          (SIGMX2 * SIGMIN)
C
 90     IF(RHO.LE.RHOMX) GOTO 200
          RHOMX = DMAX1(RHO,RHOMX)
          BETAMX = BET
 200  CONTINUE
C        BETA LOOP END-----------------------------------------
C         TRY FINER LOOP
         IF(IREP.EQ.1) GOTO 999
         IREP=1
         DBETA=0.01
         IBETA2=9
         BET=DMAX1(AMPAR2,BETAMX-0.05)
         GOTO 40
 999  RETURN
      END
C------------------------
      FUNCTION EXF(X)
      REAL*8 X,EXF
      EXF=0.
      IF(X.LT.-150.) RETURN
      EXF=EXP(X)
      RETURN
      END