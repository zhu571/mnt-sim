      SUBROUTINE  MSBEN(ZCN,ACN,NA,NZ,CMASS,SHLL)
C
C  Berechnet Separationsenergien unter Benutzung der Lysekyl -
C  Parametrisierung (1967) von Myers-Swiatecki
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BIND/BE(25,25,5),EXC(27,27),DELT,MC,MP,NUMB,MASSES
      DIMENSION EM(10),XK(10),Y(2),F(2),EMP(10),XQ(20)
C     NEEDED EXTERNAL ACN,ZCN,BE,MC,MP,NZ,NA
C     MP=1 NO PAIRING FLAG
C     MC=1 NO SHELL CORRECTION FLAG
C     DELIVERS SEPARATION ENERGIES OF N,P,A,D,T IN A SQUARE NA*NZ OF THE
C     ISOTOPE TABLE;(ACN,ZCN) IS UPPER RIGHT CORNER;USES MYERS-SWIATECKI
C     1967 MASS FORMULA;BE ARRAY FOR SEPARATION ENERGIES,EXC=ARRAY
C     CONTAINING MASS EXCESSES
C-----------------------------------------------------------------------
      DATA EM/0.0,2.,8.,14.,28.,50.,82.,126.,184.,258./,
     1     CAY1,CAY2,CAY3,CAY4,CAY5,CAY6/1.15303,0.0,200.,11.,8.07144,
     2     7.28899/,GAMMA/1.7826/,A1,A2,A3/15.4941,17.9439,0.7053/,
     3     D,C,SMALC/0.444,5.8,0.325/
C-----------------------------------------------------------------------
      IF(DELT.LE.0.) DELT=CAY4
C-----------------------------------------------------------------------
      IFLAG=0
      IF( (NA+NZ).LE.2 ) IFLAG=1
      DO 3 I=1,10
      EMP(I)=EM(I)**(5.0/3.0)
 3    CONTINUE
      DO 4 I=1,9
      XK(I)=0.6*(EMP(I+1)-EMP(I))/(EM(I+1)-EM(I))
 4    CONTINUE
      RZ=.863987/A3
      L=0
      Z=1.0
      KZ=ZCN+0.01
      KA=ACN+0.01
      NNZ=1
      NNA=1
      IF(IFLAG.EQ.1) GOTO 8
      NNZ=NZ+2
      NNA=NA+2
C
C********************************** LOOPS NZ,NA ***********************
 8    DO 32 JZ=1,NNZ
      DO 32 JA=1,NNA
 13   IA=KA+2-JA-JZ
      IZ=KZ+1-JZ
      N=IA-IZ
 14   Z=IZ
      UN=N
      A=IA
 15   A3RT=A**(1.0/3.0)
      A2RT=SQRT(A)
      A3RT2=A3RT**2.0
      ZSQ=Z**2.0
      SYM=((UN-Z)/A)**2
      ACOR=1.0-GAMMA*SYM
      PARMAS=CAY5*UN+CAY6*Z
      VOLNUC=-1.0*A1*ACOR*A
      SUFNUC=A2*ACOR*A3RT2
      COULMB=A3*ZSQ/A3RT
      FUZSUR=-1.0*CAY1*ZSQ/A
      ODDEV=0.
      IF(MP.EQ.0)
     1ODDEV=-1.0*(1.0+2.0*(N/2)-UN+2.*(IZ/2)-Z)/SQRT(A)*DELT
      WTERM=-1.*CAY2*A3RT2*EXP(-1.*CAY3*SYM)
      WOTNUC=PARMAS+COULMB+FUZSUR+ODDEV+WTERM
      SMASS=WOTNUC+VOLNUC+SUFNUC
      SHLL=0.
      XQ(JZ)=SMASS
C------------- The following added by N.Wang 2007-07-02-----------------
	CMASS=SMASS+SHLL
C-----------------------------------------------------------------------
C                                         MS   SHELL CORRECTION
      IF(MC.NE.0) GOTO 31
 16   CONTINUE
      C2=(SUFNUC+WTERM)/(A**(2.0/3.0))
      X=COULMB/(2.0*(SUFNUC+WTERM))
 17   BARR=0.0
 18   Y(1)=UN
      Y(2)=Z
      DO 22 J=1,2
      DO 19 I=1,9
      IF (EM(I+1).GT.(Y(J)-0.1)) GOTO 21
 19   CONTINUE
 20   STOP
 21   F(J)=XK(I)*(Y(J)-EM(I))-.6*(Y(J)**(5./3.)-EMP(I))
 22   CONTINUE
      S=(2.0/A)**(2.0/3.0)*(F(1)+F(2))-SMALC*A**(1./3.)
      EE=2.*C2*D**2*(1.0-X)
      FF=.42591771*C2*D**3*(1.+2.*X)/A3RT
      SSHELL=C*S
      V=SSHELL/EE
      EPS=1.5*FF/EE
      IF(EE*(1.-3.*V).LE.0.0) GO TO 23
      QCALC=0.0
      THETA=0.0
      SHLL=SSHELL
      GO TO 31
 23   TO=1.0
 24   DO 25 IPQ=1,10
      EXPDUM=0.
      DUM=TO**2
      IF(DUM.LT.25.) EXPDUM=EXP(-DUM)
      T=TO-(1.-EPS*TO-V*(3.-2.*TO**2)*EXPDUM     )/(-EPS+V*(10.*TO-4.
     1 *TO**3)*EXPDUM     )
      IF (T.LE.0.0) GO TO 27
      IF (ABS(T-TO) .LT.0.0001) GO TO 26
      TO=T
 25   CONTINUE
      GO TO 29
 26   EXPDUM=0.
      DUM=T **2
      IF(DUM.LT.25.) EXPDUM=EXP(-DUM)
      IF (2.*EE*(1.-2.*EPS*T-V*(3.-12.*T**2+4.*T**4)*EXPDUM    )
     1 .GT.0.0) GO TO 30
 27   DO 28 I=1,20
      TO=FLOAT(I)/10.
      GL=EE*(1.-EPS*TO-V*(3.-2.*TO**2)*EXP(-TO**2))
      IF (GL.GE.0.0) GO TO 24
 28   CONTINUE
 29   CONTINUE
      GO TO 32
 30   THETA=T
      ALPHA0=D*SQRT(5.)/A**(1./3.)
      ALPHA=ALPHA0*THETA
      SIGMA=ALPHA*(1.+ALPHA/14.)
      QCALC=.004*Z*(RZ*A3RT)**2*(EXP(2.*SIGMA)-EXP(-SIGMA))
      SHLL=EE*T**2-FF*T**3+SSHELL*(1.-2.*T**2)*EXP(-T**2)
 31   CMASS=SMASS+SHLL
      IF(IFLAG.EQ.0)      EXC(JA,JZ)=CMASS
      XQ(JZ)=CMASS
 32   CONTINUE
C
C******************************* END LOOPS NA,NZ **********************
 36   IF(IFLAG.EQ.1)RETURN
      DO 37 JZ=1,NZ
      DO 37 JA=1,NA
      BE(JZ,JA,1)=8.07+EXC(JA+1,JZ)-EXC(JA,JZ)
      BE(JZ,JA,2)=7.29+EXC(JA,JZ+1)-EXC(JA,JZ)
      BE(JZ,JA,4)=13.14+EXC(JA+1,JZ+1)-EXC(JA,JZ)
      BE(JZ,JA,5)=14.95+EXC(JA+2,JZ+1)-EXC(JA,JZ)
 37   BE(JZ,JA,3)=2.42+EXC(JA+2,JZ+2)-EXC(JA,JZ)
 38   RETURN
      END