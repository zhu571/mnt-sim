      SUBROUTINE  GAMMAS(ACN,IOUT,DEL,IGAM,NOE2,STRIPE,IOPT)
C
C  Berechnet Gamma - Staerkefunktionen
C  modifiziert am 1.6.1982 und am 25.4.1983
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/GAM/ECUBE(26),EFIVE(26),CGIANT,EGIANT,GAMMA,GAMQ,KEG,JFACTR
      DIMENSION ARRAY(10)
C     GAMMA NORMALISATION
C--------------------------------------------------------------
C     DEFAULTS  E2:TEN WEISSKOPF UNITS, 1 TO 4 MEV
C               E1:GIANT DIPOLE RESONANCE 1 TO (EGIANT+WGIANT) MEV
      DEL=1.
      IOPT=0
      STRIPE=0.
      JFACTR=0
      EG1WU= 0.
      EG2WU=10.
      EG2MIN=1.
      EG2MAX=4.
      CGIANT=  11.5/3.
      EGIANT=80./(ACN**0.33333)
      WGIANT=5.
      DUM=EGIANT+WGIANT
      EG1MIN= 1.
      EG1MAX=-1.
      DO 10 KE=1,KEG
      ECUBE(KE)=0.
 10   EFIVE(KE)=0.
      IF(IGAM.EQ.0)GOTO 60
C---------------------------------------------------------------
      CALL MYRD(ARRAY,7,23,5)
      EG1WU =ARRAY(1)
      EG2WU =ARRAY(2)
      EG1MIN=ARRAY(3)
      EG1MAX=ARRAY(4)
      EG2MIN=ARRAY(5)
      EG2MAX=ARRAY(6)
      JFACTR=ARRAY(7) + 0.01
      CALL MYRD(ARRAY,5,23,5)
      CGIANT=ARRAY(1)
      EGIANT=ARRAY(2)
      WGIANT=ARRAY(3)
      STRIPE=ARRAY(4)
      IOPT  =ARRAY(5)+0.01
      IF(EG2WU.LE.0.)NOE2=1
      IF(EG2MAX.LE.0.)EG2MAX=4
      IF(EGIANT.LE.0.)EGIANT=80/(ACN**0.33333)
      IF(WGIANT.LE.0.)WGIANT=5.
C---------------------------------------------------------------
 60   IDEL=DEL+0.01
      DEL=IDEL
      IF(DEL.LE.1.) DEL=1.
      IF(KEG.EQ.0)KEG=26
C---------------------------------------------------
      IF(CGIANT.LE.0.) GOTO 180
C
C**************************** LORENTZIAN ******************************
 160  CONTINUE
      EG2=EGIANT*EGIANT
      W2=WGIANT*WGIANT
      CONST=(1.63E-6)*CGIANT*ACN*WGIANT
      FK=DEL
      DO 163 K=1,KEG
      FK2=FK*FK
      EFIVE(K)=0.
      DUM=EG2-FK2
      ECUBE(K)=CONST*FK2*FK2/(DUM*DUM+FK2*W2)
      FK=FK+1.
 163  CONTINUE
C-------------------------------------------------------------------
 180  GNORM=(6.75E-8)*(ACN**0.6667)*6.283*EG1WU
      QNOR= (4.77E-14)*(ACN**1.3333)*6.283*EG2WU
 95   CONTINUE
      FK=DEL
      DO 80 K=1,KEG
      IF(FK.LT.EG1MIN .OR. FK.GT.EG1MAX)GOTO 70
      ECUBE(K)=GNORM*FK**3
 70   IF(FK.LT.EG2MIN .OR. FK.GT.EG2MAX)GOTO 80
      EFIVE(K)=QNOR*FK**5
 80   FK=FK+1.
      FKSI1=GNORM/6.283
      FKSI2=QNOR/6.283
C--------------------------------------------------------------
      IF(IOUT.EQ.0) THEN
         WRITE(6,230) (ECUBE(I),I=1,20)
 230     FORMAT(/' GAMMAS:   ECUBE'/4(5E12.3/))
      ENDIF
      RETURN
      END