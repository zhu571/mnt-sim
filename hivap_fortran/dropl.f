      SUBROUTINE DROPL(Z,A,ICODE,DROP,DEL,EPS,LDM,IPRNT)
C
C
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BEES/BS,BC,BK,BR,BV,BW
      IF(LDM.EQ.1) GOTO 15
C----------------------------------------
             A1=               15.96
             A2=               20.69
             A3=               0.
             ZJ=               36.8
             Q=                17.0
             RZ=               1.18
             ZK=               240.
             ZL=               100.
             ZM=               0.
C-----------------------------------------------------
          RNULL = RZ
         R3 = 1.E0/3.E0
         C1 = 0.6E0*1.4399784E0/RZ
         C2 = (C1**2.E0/168.E0) * (0.5E0/ZJ + 9.E0/ZK)
         C3 = 2.5E0 * C1 * (0.99E0/RZ)**2
         C4 = 1.25E0 * (1.5E0 / 3.14159265359E0)**(2.E0/3.E0)*C1
         C4P = C4/1.25992105E0
         C5 = (C1**2/Q)/64.E0
         T3 = 0.1875E0 * C1 / Q
         T4 = 2.25E0 * ZJ/Q
         T5 = 2.E0 * A2 / ZK
         T6 = ZL / ZK
         T7 = C1 / ZK
         T8 = ZJ * T4
         T9 = 0.5E0 * ZK
         T10 = 0.5E0 * ZM
C------------------------------------
 15      IF(ICODE.EQ.2) GOTO 20
         BS=1.
         BC=1.
         BK=1.
         BR=1.
         BV=1.
         BW=1.
C-------------------------------------
 20      X = A**.333333
         X2 = X*X
         AS = 1.E0 - 2.E0 * Z/A
C
         IF(LDM.EQ.1) GOTO 40
         DEL = (AS + T3 * Z / X2 * BV)/(1.E0 + T4 / X * BS )
         D2 = DEL*DEL
         EPS = (- T5 / X * BS + T6 * D2 + T7 * Z * Z / X2**2 * BC)
         E2 = EPS*EPS
C
         V1=-A1*A
         V2=ZJ*D2*A
         V3=-  E2*A *T9
         V4=T10*D2*D2*A
         S1=A2*X2*BS
         S2=T8*D2*X2*BS
         CURV=A3*X*BK
         COUL1=Z*Z*C1*BC/X
         COUL2=- Z*C2*X*BR *Z
         COUL3=-Z*Z*C3/A
         COUL4=-C4P*Z
         COUL5=- Z*C5*BW*Z
C--------------------------------------------------------
      IF(LDM.EQ.0) GOTO 50
 40   A1=15.4941
      RKAPPA=1.7826
C     RKAPPA=2.815
      ZJJ=A1*RKAPPA
      A2=17.9439
C     A2=18.7081
      C1=0.7053
C     C1=0.6765
      T88=-A2*RKAPPA
      AS2=AS*AS
      V1=   -A1*A
      V2=   ZJJ*AS2*A
      V3=0.
      V4=0.
      S1= A2*X2*BS
      S2=T88*AS2*X2*BS
      CURV=0.
      COUL1= Z*Z*C1*BC/X
      COUL2=0.
      COUL3=-Z*Z*1.15303/A
      COUL4=0.
      COUL5=0.
C----------------------------------------------------------
  50  DROP=COUL1+COUL2+COUL3+COUL4+COUL5+V1+V2+V3+V4+S1+S2+CURV
      IF(IPRNT.NE.1)RETURN
      WRITE(6,101) Z,A,LDM,DROP,EPS,DEL
 101  FORMAT(/' Z,A',2F6.0,4X,'LDM',I4,4X,'DROP',F10.2/' EPS,DEL',2E12.4
     1)
      WRITE(6,102) V1
 102  FORMAT(' BULK VOLUME', 6X,F10.2)
      WRITE(6,103) V2
 103  FORMAT(' ASYM VOLUME', 6X,F10.2)
      WRITE(6,104) V3
 104  FORMAT(' COMP VOLUME', 6X,F10.2)
      WRITE(6,105) V4
 105  FORMAT(' ASYM VOLUME 2', 4X,F10.2)
      VV=V1+V2+V3+V4
      WRITE(6,115) VV
 115  FORMAT(' TOTAL VOLUME ',14X,F10.2)
      WRITE(6,106) S1
 106  FORMAT(' BULK SURFACE ', 4X,F10.2)
      WRITE(6,107) S2
 107  FORMAT(' ASYM SURFACE ', 4X,F10.2)
      SS=S1+S2
      WRITE(6,116) SS
 116  FORMAT(' TOTL SURFACE ',14X,F10.2)
      WRITE(6,108) COUL1
 108  FORMAT(' COULOMB      ', 4X,F10.2)
      WRITE(6,109) COUL2
 109  FORMAT(' COUL VOL RED ', 4X,F10.2)
      WRITE(6,112) COUL5
 112  FORMAT(' COUL SURF RED', 4X,F10.2)
      COUL=COUL1+COUL2+COUL5
      WRITE(6,117) COUL
 117  FORMAT(' TOTAL COUL   ',14X,F10.2)
      WRITE(6,113) CURV
 113  FORMAT(' CURVATURE    ', 4X,F10.2)
      WRITE(6,110) COUL3
 110  FORMAT(' COUL DIFFUS  ', 4X,F10.2)
      WRITE(6,111) COUL4
 111  FORMAT(' COUL EXCHANGE', 4X,F10.2)
         RETURN
         END