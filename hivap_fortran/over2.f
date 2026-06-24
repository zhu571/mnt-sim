      SUBROUTINE OVER2(K,KEMAX,LMAX,TLE,AMASS,ZEE,TLIM,SIGLOW,TLOW,
     1                 KDIM,IPRNT)
C
C  Berechnet Transmissionskoeffizienten fuer n,p,alphas fuer Abdampfung
C  nach dem optischen Modell, ruft TDL auf
C  Stand: 3.5,1994
C
C     MOD 2 NOV 81 TO FEED IN POTS FROM MAIN
C     MOD 20 JUN 90 SATCHLER POT FOR ALPHAS IS DEFAULT
C     CALCULATES TRANSMISSION COEFFICIENTS FOR HIVAP
C     FOR DEUTERON AND TRITON EVAPORATION THIS MUST BE COMPLETED BY INSE
C     RELEVANT OPTICAL MODEL PARAMETERS
C     JE=1 CORRESPONDS TO E=TLIM
C     K=1 NEUTRONS,=2 PROTONS,=3 ALPHAS,=4 DEUTERONS,=5 TRITONS
C     KEMAX=MAXIMUM NUMBER OF ENERGIES FOR WHICH TRANSMISSION IS CALCUL.
C     IF 'KEMAX' IS NEGATIVE,WILL INTERPRETE 'SIGLOW' AS ENERGY STEP(MEV
C     ELSE STEP IS 1MEV AND SIGLOW IS CUTOFF ABSORB. CROSS SECTION(MB)
C     'TLIM' IS LOWEST ENERGY,'LMAX' IS MAX ORBITAL ANG.MOMENTUM
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TLE (1) ,ARRAY(10)
      COMMON/TLJ/ T(2,31),V(15),V1(3)
      COMMON/POTS/ VRN,R0RN,ADIFRN,VIN,R0IN,ADIFIN,
     1             VRP,R0RP,ADIFRP,VIP,R0IP,ADIFIP,RCLMBP,CBFACP,
     2             VRA,R0RA,ADIFRA,VIA,R0IA,ADIFIA,RCLMBA,CBFACA,QQ2
C      DOUBLE PRECISION T,V,V1,H4,H1,Z,W1,STPLTH
C-------------------------------------------------------------------
 1001 FORMAT(6E10.4)
 1002 FORMAT(1I1)
 1003 FORMAT(1H ,8F10.3  )
 1004 FORMAT(6F10.0)
 1008 FORMAT(9F8.3)
      DDE=1.
      IF(KEMAX.GT.0.) GOTO 21
      KEMAX=IABS(KEMAX)
      DDE=SIGLOW
      SIGLOW=0.
 21   IDUM=KEMAX*LMAX
      DO 20 I=1,IDUM
 20   TLE(I)=0.
C--------------------------------------------------------------------
C                                      PARTICLE LOOP
C
      K=IABS(K)
      E=TLIM
      JE1=TLIM+1.
      E1=0.
      E2=0.
      E3=0.
      TSUM1=0.
      TSUM2=0.
      TSUM3=0.
      TSUM = 0.
      XLM=LMAX
      GOTO (1,2,3,320,340),K
C-----------------------------------   NEUTRON  --------------
C     DEFAULT IS    WILMORE-HODGSON NP55(64)673,ALSO PEREY NDTA10(72)539
 1    XJP=0.50
      XMP=1.00
      XMT=AMASS-1.0
      ZP=0.0
      ZT=ZEE
      VR   =VRN
      IDEFN=0        ! USER POTENTIALS NOT ENERGY DEPENDENT
      RV   =R0RN
      AV   =ADIFRN
      W    =VIN
      IDEFW=0        ! USER POTENTIALS NOT ENERGY DEPENDENT
      RW   =R0IN
      AW   =ADIFIN
      CBFAC=1.
C
 200  IF(VR.LE.0.) IDEFN=1     ! ENERGY DEPENDENCE IN ENERGY LOOP
      IF(VR.LE.0.) VR=47.01
      IF(RV.LE.0.) RV=1.322  -7.6E-4*XMT  +4.E-6*XMT*XMT  -8.E-9*XMT**3.
      IF(AV.LE.0.) AV=.660
      IF(W.LE.0.) IDEFW=1       ! ENERGY DEPENDENCE IN ENERGY LOOP
      IF(W.LE.0.) W=9.520
      IF(RW.LE.0.) RW=1.266-3.7E-4*XMT  +2.E-6*XMT*XMT  -4.E-9*XMT **3.
      IF(AW.LE.0.) AW=0.48
      S=1.00         ! SURFACE ABSORPTION FLAG
C
 205  VSO=7.00       ! SPIN ORBIT
      RSO=RW
      ASO=AW
      RCLMB=.001
      GO TO 4
C  ---------------------------------   PROTON ---------------------
C      DEFAULT:    PROTON POT BECCHETTI-GREENLEES PR182(1969)1190
 2    XJP=0.50
      XMP=1.00
      XMT=AMASS-1.0
      ZP=1.0
      ZT=ZEE-1.0
      VR   =VRP
      IDEFP=0
      RV   =R0RP
      AV   =ADIFRP
      W    =VIP
      IDEFPW=0
      RW   =R0IP
      AW   =ADIFIP
      RCLMB=RCLMBP
      CBFAC =CBFACP
 210  IF(VR.LE.0.) IDEFP=1       ! ENERGY DEPENDENCE IN ENERGY LOOP
      IF(VR.LE.0.)  VR=54.  +24.*(XMT -2.*ZT)/XMT  +0.4*ZT/XMT**0.33333
      IF(RV.LE.0.)  RV=1.17
      IF(AV.LE.0.)  AV=.750
      IF(W.LE.0.)   IDEFPW=1      ! ENERGY DEPENDENCE IN ENERGY LOOP
      IF(W.LE.0.)   W=11.8 + 12.*(XMT-2.*ZT)/XMT
      IF(RW.LE.0.)  RW=1.32
      IF(AW.LE.0.)  AW=.51 + 0.7*(XMT-2.*ZT)/XMT
      IF(RCLMB.LE.0.) RCLMB=1.16
      IF(CBFAC.LE.0.) CBFAC=1.
      S=1.00                        ! SURFACE ABSORPTION FLAG
C
 215  VSO=6.2                       ! SPIN ORBIT
      RSO=1.01
      ASO=0.75
      GO TO 4
C------------------------------------  ALPHA ---------------
C     DEFAULT                          SATCHLER,NP70(1965)177
 3    XJP=0.0
      XMP=4.0
      XMT=AMASS-4.0
      ZP=2.00
      ZT=ZEE-2.0
      VR   =VRA
      RV   =R0RA
      AV   =ADIFRA
      W    =VIA
      RW   =R0IA
      AW   =ADIFIA
      RCLMB=RCLMBA
      CBFAC=CBFACA
C
      IF(VR.LE.0.)    VR=50.2
      IF(RV.LE.0.)    RV=1.2  +1.5/XMT**0.33333
      IF(AV.LE.0.)    AV=.564
      IF(W.LE.0.)     W=12.30
      IF(RW.LE.0.)    RW=1.2  +1.5/XMT**0.33333
      IF(AW.LE.0.)    AW=AV
      IF(RCLMB.LE.0.) RCLMB=1.30
      IF(CBFAC.LE.0.) CBFAC=1.
      S=0.0           ! VOLUME ABSORPTION FLAG
      RSO=1.
      VSO=.001
      ASO=1.0
      GOTO 4
C------------------------------------  DEUTERON ----------------
 320  CONTINUE
      GOTO 4
C------------------------------------  TRITON ---------------
 340  CONTINUE
C-----------------------------------------------------------------------
 4    IF(IPRNT.NE.1) GOTO 350
      WRITE (6,1005)
 1005 FORMAT( /' PROJ. SPIN  A(PROJ.) A(TARGET)  Z(PROJ.) Z(TARGET)')
      WRITE(6,1003)XJP,XMP,XMT,ZP,ZT
      WRITE (6,1006)
 1006 FORMAT( /  '    V(REAL)   R(REAL)   A(REAL)  W(IMAG.)  R(IMAG.)  A
     1(IMAG.) ')
      WRITE(6,1003)VR,RV,AV,W,RW,AW
      WRITE (6,1007)
 1007 FORMAT (/  '      V(SO)     R(SO)     A(SO)   R(COUL)        S   L
     1(LIMIT) ')
      WRITE(6,1003)VSO,RSO,ASO,RCLMB,S,XLM
 350  CONTINUE
      XF=XMT**.3333
      XM=XMT/(XMT + XMP)
      XJT=0.0
      IC=0
C-----------------------------------------------------------------------
C                                      ENERGY LOOP
C
      KE=0
      DO 130 JE=JE1,80
 370  V(1)=RV*XF
      V(2)=AV
      V(3)=RW*XF
      V(4)=AW
      V(5)=S
      V(6)=VSO
      V(9)=RCLMB*XF
      V(10)=RSO*XF
      V(11)=ASO
      V1(1)=XJP
      V1(2)=+0.0
      W1=0.04783258*XMP*XM
      STPLTH=0.1
      V1(3)=XLM
      V(7)=VR
      IF(K.EQ.1 .AND. IDEFN.EQ.1) V(7)=VR-0.267*E   ! ENERGY DEP
      IF(K.EQ.2 .AND. IDEFP.EQ.1) V(7)=VR-0.32 *E
      V(8)=W
      IF(K.EQ.1 .AND. IDEFW.EQ.1)V(8)=W-0.053*E
      IF(K.EQ.2 .AND. IDEFPW.EQ.1)V(8)=W-0.25*E
      H4=0.04783*XMP*(XM**2)*E
      H1=DSQRT(H4)
      Z=CBFAC*(0.03478*ZP*ZT*XMP*XM)/H1
      C6=31.42/H4
      DO 5 J=1,2
      DO 5 I=1,31
 5    T(J,I)=0.0
C--------------------------------------------------------------------
      CALL TLD(H4,N1,J1,H1,Z,W1,STPLTH)
C--------------------------------------------------------------------
      TSUM=0.0
      AVSUM=0.
      INTRPO=2.*XJP
C
      DO 7 I=1,LMAX
      AL=I+I-1
      IF(INTRPO.EQ.0) GOTO 6
      TK=I-1
      FL=1./(2.*TK+1)
      T(1,I)=FL*(T(1,I)*TK + T(2,I)*(TK+1.))
 6    TSUM=TSUM + T(1,I)*AL*C6
 7    CONTINUE
C
      IF(TSUM.GE.SIGLOW  ) GOTO 70
      TLIM=TLIM+1.
      GOTO 10
 70   KE=KE+1
      DO 8 I=1,LMAX
      IDUM=KE + KDIM *(I-1)
      TLE(IDUM)=T(1,I)
      IF(TLE(IDUM).LT.TLOW) TLE(IDUM)=0.
 8    CONTINUE
      IF(KE.GE.KEMAX) GOTO 14
 10   CONTINUE
      E=E+DDE
 130  CONTINUE
C
C                                      STOP ENERGY LOOP
C
 14   CONTINUE
C
C                                      STOP PARTICLE LOOP
C
 15   RETURN
      END