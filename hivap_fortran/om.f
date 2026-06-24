      SUBROUTINE OM (K,ELAB,LMAX,ATARG,ZTARG,SIGTOT,SIGML,IPOT,IPRNT)
C
C  Berechnet Absorptionsquerschnitte fuer n,p,a nach dem optischen
C  Modell
C  Stand: 3.5.1994
C
C     CALCULATES ABSORPTION CROSS SECTIONS FOR NEUTRONS,PROTONS,ALPHAS
C     FOR DEUTERON AND TRITON  THIS MUST BE COMPLETED BY INSERTING
C     RELEVANT OPTICAL MODEL PARAMETERS
C     K=1 NEUTRONS,=2 PROTONS,=3 ALPHAS,=4 DEUTERONS,=5 TRITONS
C     ELAB    =PROJECTILES ENERGY IN MEV
C     LMAX    MAXIMUM+1 ORBITAL ANGULAR MOMENTUM
C     ATARG,ZTARG  MASS AND Z OF TARGET NUCLEUS
C     SIGTOT  TOTAL ABSORPTION CROSS SECTION IN MILLIBARN(OUTPUT)
C     SIGML   PARTIAL CROSS SECTIONS(AVERAGED FOR SPIN ORBIT) (OUTPUT)
C     IPOT    =0  TAKE DEFAULT POTENTIAL
C             =1  READ IN NEW POTENTIAL
C             =2  KEEP OLD POTENTIAL
C     IPRNT   PRINT FLAG
C     T       TRANSMISSION COEFFICIENTS (OUTPUT)
C----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIGML (1) ,ARRAY(10)
      COMMON/TLJ/ T(2,31),V(15),V1(3)
      DOUBLE PRECISION T,V,V1,H4,H1,Z,W1,STPLTH
C-------------------------------------------------------------------
 1001 FORMAT(6E10.4)
 1002 FORMAT(1I1)
 1003 FORMAT(1H ,8F10.3  )
 1009 FORMAT(1H ,8F10.3//)
 1004 FORMAT(6F10.0)
 1008 FORMAT(9F8.3)
      IF(LMAX.GT.31) LMAX=31
      DO 20 I=1,LMAX
 20   SIGML(I)=0.
C--------------------------------------------------------------------
C ------------------------------------ PARTICLE SWITCH
C
      K=IABS(K)
      ISAME=0
      IF(K.EQ.KOLD) ISAME=1
      KOLD=K
      IPOT1=IPOT
      IF(ISAME.EQ.0 .AND.IPOT.EQ.2) IPOT1=0
      E=ELAB
      TSUM = 0.
      XLM=LMAX
      GOTO (1,2,3,320,340),K
C------------------------------------- NEUTRON
 1    XJP=0.50
      XMP=1.00
      XMT=ATARG
      ZP=0.0
      ZT=ZTARG
      IF(IPOT1.EQ.0) GOTO 200
      IF(IPOT1.EQ.2) GOTO 205
C     READ(5,1008) P,RV,AV,W,RW,AW
      IF(IPRNT.NE.13) CALL MYRD(ARRAY,7,6,5)
      IF(IPRNT.EQ.13) CALL MYRD(ARRAY,7,13,5)
      P    =ARRAY(1)
      RV   =ARRAY(2)
      AV   =ARRAY(3)
      W    =ARRAY(4)
      RW   =ARRAY(5)
      AW   =ARRAY(6)
C     RCLMB=ARRAY(7)
      GOTO 205
 200  P=47.01 - 0.267*E
      RV=1.322-7.6E-4*XMT  +4.E-6*XMT*XMT      -8.E-9* XMT **3.
      AV=.660
      W=9.520 -0.053*E
      RW=1.266-3.7E-4*XMT  +2.E-6*XMT  *XMT  -4.E-9*XMT **3.
      AW=0.48
 205  VSO=7.00
      RSO=RW
      ASO=AW
      RCLMB=.001
      S=1.00                 ! SURFACE IMAG
      GO TO 4
C ------------------------------------ PROTON
 2    XJP=0.50
      XMP=1.00
      XMT=ATARG
      ZP=1.0
      ZT=ZTARG
      IF(IPOT1.EQ.0) GO TO 210
      IF(IPOT1.EQ.2) GO TO 215
C     READ(5,1008) P,RV,AV,W,RW,AW,RCLMB
      IF(IPRNT.NE.13) CALL MYRD(ARRAY,7,6,5)
      IF(IPRNT.EQ.13) CALL MYRD(ARRAY,7,13,5)
      P    =ARRAY(1)
      RV   =ARRAY(2)
      AV   =ARRAY(3)
      W    =ARRAY(4)
      RW   =ARRAY(5)
      AW   =ARRAY(6)
      RCLMB=ARRAY(7)
      GOTO 215
C     BECCHETTI-GREENLEES PR182(1969)1190
 210  P=54.        +24.*(XMT  -2.*ZT)/XMT  +0.4*ZT/XMT  **0.33333
     1             -0.32*E
      RV=1.17
      AV=.750
      W=11.8 +12.*(XMT - 2.*ZT)/XMT-0.25*E
      RW=1.32
      AW=.51 + 0.7*(XMT - 2.*ZT)/XMT
      RCLMB=1.16
C     SPIN ORBIT PEREY 1963
 215  RSO=1.01
      VSO=6.2
      ASO=.75
      S=1.00             ! SURFACE IMAG
      GO TO 4
C--------------------------------------ALPHA
 3    XJP=0.0
      XMP=4.0
      XMT=ATARG
      ZP=2.00
      ZT=ZTARG
      IF(IPOT1.EQ.0) GO TO 220
      IF(IPOT1.EQ.2) GO TO 225
C     READ(5,1008) P,RV,AV,W,RW,AW,RCLMB
      IF(IPRNT.NE.13) CALL MYRD(ARRAY,7,6,5)
      IF(IPRNT.EQ.13) CALL MYRD(ARRAY,7,13,5)
      P    =ARRAY(1)
      RV   =ARRAY(2) +1.5/XMT**0.33333
      AV   =ARRAY(3)
      W    =ARRAY(4)
      RW   =ARRAY(5)
      AW   =ARRAY(6)
      RCLMB=ARRAY(7)
      GOTO 225
 220  P=50.2
C     RV=1.2+1.5/XMT**0.33333          SATCHLER
      RV=1.275+1.5/XMT**0.33333
      AV=.564
      W=12.30
C     RW=RV                            SATCHLER
      RW=1.25
      AW=.564
      RCLMB=1.30
 225  RSO=1.
      VSO=.001
      ASO=1.0
      S=0.0              ! VOLUME IMAG
      GOTO 4
C                                      DEUTERON
 320  CONTINUE
      GOTO 4
C                                      TRITON
 340  CONTINUE
C------------------------------------------------------------------
 4    IF(IPRNT.EQ.0) GOTO 350
      WRITE (6,1005)
 1005 FORMAT(/' PROJ. SPIN  A(PROJ.) A(TARGET)  Z(PROJ.) Z(TARGET)')
      WRITE(6,1003)XJP,XMP,XMT,ZP,ZT
      WRITE (6,1006)
 1006 FORMAT( '    V(REAL)   R(REAL)   A(REAL)  W(IMAG.)  R(IMAG.)  A
     1(IMAG.) ')
      WRITE(6,1003)P,RV,AV,W,RW,AW
      WRITE (6,1007)
 1007 FORMAT ('      V(SO)     R(SO)     A(SO)   R(COUL)        S   L
     1(LIMIT) ')
      WRITE(6,1009)VSO,RSO,ASO,RCLMB,S,XLM
C---------------------------------------------------------------------
      IF(IPRNT.NE.13) GOTO 350
      WRITE (13,1005)
      WRITE(13,1003)XJP,XMP,XMT,ZP,ZT
      WRITE (13,1006)
      WRITE(13,1003)P,RV,AV,W,RW,AW
      WRITE(13,1007)
      WRITE(13,1009)VSO,RSO,ASO,RCLMB,S,XLM
C---------------------------------------------------------------------
 350  CONTINUE
      XF=XMT**.3333
      XM=XMT/(XMT + XMP)
      XJT=0.0
      IC=0
C
C--------------------------------------------------------------------
C
      V(1)=RV*XF
      V(2)=AV
      V(3)=RW*XF
      V(4)=AW
      V(5)=S       ! IS ONE FOR SURFACE IMAG
      V(6)=VSO
      V(9)=RCLMB*XF
      V(10)=RSO*XF
      V(11)=ASO
      V1(1)=XJP
      V1(2)=+0.0
      W1=0.04783258*XMP*XM
      STPLTH=0.1
      V1(3)=XLM
      V(7)=P
      V(8)=W
      H4=0.04783*XMP*(XM**2)*E
      H1=DSQRT(H4)
      Z=(0.03478*ZP*ZT*XMP*XM)/H1
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
 6    SIGML(I)=T(1,I)*AL*C6
      TSUM=TSUM+SIGML(I)
 7    CONTINUE
      SIGTOT=TSUM
C
C
C------------------------------------------------------------------
C
      RETURN
      END