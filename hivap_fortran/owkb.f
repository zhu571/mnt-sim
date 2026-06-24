      SUBROUTINE OWKB(K,KEMAX,LMAX,TLE,APAR,ZPAR,TLIM,SIGLOW,TLOW,
     1                KDIM,Q2,IPRNT)
C
C  Berechnet Transmissionskoeefizienten fuer n,p,alphas fuer Abdampfung
C  nach der WKB - Methode, ruft FUSION auf
C  Stand: 3.5.1994
C
C     CALCULATES TRANSMISSION COEFF. N,P,A FOR EVAP.THEORY WITH WKB
C     JE=1 CORRESPONDS TO TLIM(LOWEST ENERGY)
C     K=1 NEUTRONS,=2 PROTONS,=3 ALPHAS
C     APAR,ZPAR PARENT NUCLEUS
C     KEMAX=MAXIMUM NUMBER OF ENERGIES FOR WHICH TRANSM.IS CALCULATED
C     LMAX=MAX.ORBITAL L FOR EMITTED N,P,A
C---------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TLE(1),SIG(100),SIGML(60)
      COMMON/POTS/ VRN,R0RN,ADIFRN,VIN,RIN,ADIFIN,
     1             VRP,R0RP,ADIFRP,VIP,RIP,ADIFIP,RCLMBP,CBFACP,
     2             VRA,R0RA,ADIFRA,VIA,RIA,ADIFIA,RCLMBA,CBFACA,QQ2
      COMMON/FUS/VBFUS,RBFUS,HWFUS
      K=IABS(K)
      HBAR=.65819D-21
      AMU=1.03641D-44
      R0FUS=0.
      CRED=1.
      NORUTH=1
      NOPROX=1
      NOCURV=1
      ION=1
      ITEST=0
      IOPT=-1
C---------------------------
      GOTO (10,20,30),K
C                                   NEUTRONS PEREY
 10   AP=1.
      ZP=0.
      AT=APAR-AP
      V0=VRN
      R0=R0RN
      D=ADIFRN
      IF(V0.LE.0.) V0=47.01
      IF(R0.LE.0.) R0=1.322-7.6E-4* AT+4.E-6* AT* AT-8.E-9* AT**3.
      R0T=R0
      IF(D.LE.0.)  D=0.660
      IF(IPRNT.NE.0) WRITE(6,15) V0,R0,D
 15   FORMAT(/' OWKB NEUTRONS',4X,'V0=',F7.2,4X,'R0=',F7.3,4X,'D=',F7.3)
      GOTO 35
C-------------------PROTONS OWN FIT
 20   AP=1.
      ZP=1.
      AT=APAR-AP
      V0=VRP
      R0=R0RP
      D=ADIFRP
C     PEREY 1963 POT
C     IF(V0.LE.0.) V0=53.3+27.*( AT-2.*ZT)/AT +0.4*ZT/AT**0.33333
C     IF(R0.LE.0.) R0=1.25
C     FOLLOWING PARAMETERS FOR TEST PURPOSES: FITTED TO DATA HOFMANN
C      V0=60
C      R0=1.40
C     WKB POT FITTED TO KURCEWICZ DATA TH232+P(1981)
      IF(V0.LE.0.) V0=62
      IF(R0.LE.0.) R0=1.254
      R0T=R0
      IF(D.LE.0.)  D=0.750
      IF(IPRNT.NE.0) WRITE(6,25)
 25   FORMAT(/' OWKB  PROTONS')
      GOTO 35
C--ALPH-MODIFIED  SATCHLER (IGO USE V0=1100.;R0=1.17-1.6/AT3;D=0.574)
C             FITTED TO ALPHA DATA SEP 81
 30   AP=4.
      ZP=2.
      AT =APAR-AP
      AT3=AT**0.333333
      V0=VRA
      R0=R0RA
      D=ADIFRA
      IF(V0.LE.0.) V0=50.2
      IF(R0.LE.0.) R0=1.2067
      R0T=R0+1.6/AT3
      IF(D.LE.0.)  D=0.564
      IF(IPRNT.NE.0) WRITE(6,34)
 34   FORMAT(/' OWKB  ALPHAS ')
C---------------------------
 35   Q2=QQ2
C---------------------------
 45   ELAB=TLIM -1.
      AT=APAR-AP
      ZT=ZPAR-ZP
      ARED=AMU*AP*AT/(AP+AT)
      HW0=HBAR*HBAR/ARED
      KE=0
C-----------------------------
 50   ELAB=ELAB+1.
      IF(ELAB.GT.150) GOTO 1999
      ECM=ELAB*AT/(AP+AT)
      SIG0=15.708*HW0/ECM
      LMX =-IABS(LMAX)
      CALL FUSION(ELAB,AP,ZP,AT,ZT,Q2,V0,R0,D,0.D0,0.D0,0.D0,
     1            CRED,NORUTH,
     2            NOPROX,NOCURV,ION,SIGML,FLGRAZ,LMX ,SIGF,ITEST,IOPT)
      IF(ITEST.EQ.99) GOTO 1999
      IF(SIGF.GT.SIGLOW) GOTO 55
      TLIM=TLIM+1.
      GOTO 50
 55   KE=KE+1
      DO 60 L=1,LMX
      IDUM=KE + KDIM *(L-1)
      FL=L+L-1
      TLE(IDUM)  =SIGML(L)/(SIG0*FL)
      IF(TLE(IDUM).LT.TLOW) TLE(IDUM)=0.
      SIG(KE)=SIGF
 60   CONTINUE
      IF(KE.LT.KEMAX) GOTO 50
C------------------------------
      IF(IPRNT.NE.0) WRITE(6,315) K,LMAX,APAR,ZPAR,TLIM,V0,R0,D,Q2,R0T
 315  FORMAT(/' OWKB K,LMAX,A,Z,EMIN',2I3,3F8.1/  ' V0,R0,D,Q2=', F8.2,
     1 2F9.4,F8.1,4X,'R0T=',F9.4)
 999  RETURN
1999  WRITE(6,2000) ELAB,ITEST
2000  FORMAT(' CALCULATION STOPS BECAUSE PROBLEM WITH WKB TRANSMISSION'/
     1       ' ELAB,ITEST=',F10.1,I4)
      STOP
      END