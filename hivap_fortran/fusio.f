      SUBROUTINE FUSIO (ELAB,AP,ZP,AT,ZT,Q2,V0,R0,D,XTH,APUSH,FPUSH,
     2            SIGR0,CRED,NORUTH,NOPROX,NOCURV,
     1            ION,SIGML,FLGRAZ,LMAX,SIGF,ITEST,IOPT,MPT,CUTOFF)
C
C
C   modifiziert am 8.7.1987 (new fusion) und am 8.10.1987 (limbar)
C   Stand: 3.5.1994
C
C                      PARAMETERS ARE AS VLRO, AND FUSION
C     SIGR0 IS THE PERCENT CHANGE OF R0 (STAND.DEV.)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NPOINT(7),KEY(8),Z(24),WEIGHT(24)
      DATA NPOINT/2,3,4,5,6,10,15/,KEY/1,2,4,6,9,12,17,25/,
     1     Z/0.5773502,0.0,      0.7745966,
     2       0.3399810,0.8611363,0.0,      0.5384693,
     3       0.9061798,0.2386191,0.6612093,0.9324695,
     4       0.1488743,0.4333953,0.6794095,0.8650633,
     5       0.9739065,0.0,      0.2011940,0.3941513,
     6       0.5709721,0.7244177,0.8482065,0.9372733,
     7       0.9879925/,
     8WEIGHT/1.0000000,0.8888888,0.5555555,
     9       0.6521451,0.3478548,0.5688888,0.4786286,
     A       0.2369268,0.4679139,0.3607615,0.1713244,
     B       0.2955242,0.2692667,0.2190863,0.1494513,
     C       0.0666713,0.2025782,0.1984314,0.1861610,
     D       0.1662692,0.1395706,0.1071592,0.0703660,
     E       0.0307532/
C
C--------------------------------------------------
      DIMENSION SIGML(1),SIGL(200),SIGL1(200)
      COMMON/FUS/VBFUS,RBFUS,HWFUS
C------------------------------------------------
C     M-POINT GAUSS QUADRATURE,M=2,3,4,5,6,10,15;DEFAULT IS 5
      M=MPT
      MTIMES=1
      DO 7  L=1,7
      IF(M.EQ.NPOINT(L)) GOTO 9
 7    CONTINUE
      IF(MPT.LT.16) GOTO 8
      MTIMES=M/5
      L=4
      GOTO 9
 8    M=10
      L=6
C     DEFAULT
 9    JFIRST=KEY(L)
      JLAST=KEY(L+1)-1
      IF(CUTOFF.LE.0.) CUTOFF=2.5
C----------------------------------------------
      LIMBAR=0
      IF (LMAX.LT.0) LIMBAR=1
      R0FUS=0.
      DO 40 L=1,200
 40   SIGML(L)=0.
      LMAX2=0
      NOQ2=0
      IF(ABS(Q2).LT.0.001) NOQ2=1
 50   ION1=ION
      IF(NOQ2.EQ.1)ION1=4
      IF(SIGR0.GT.0.) GOTO 80
C----------------------------------------------
      LMAX=0
      IF(LIMBAR.EQ.1) LMAX=-1
      CALL FUSION(ELAB,AP,ZP,AT,ZT,Q2,V0,R0,D,XTH,APUSH,FPUSH,
     1            CRED,NORUTH,NOPROX,
     2            NOCURV,ION1,SIGML,FLGRAZ,LMAX,SIGF,ITEST,IOPT)
      RETURN
C-----------------------------------------------
 80   AP3=AP**0.333333
      AT3=AT**0.333333
      FLUCT=0.01*SIGR0
      XLIM=CUTOFF*FLUCT
      LMAX2=0
C---------------------------------------------------------------------
C     RANGE OF INTEGRATION IS XIN TO XFI
C--------------------------------------------------

 100  XIN=R0*(1.-XLIM)
      XFI=R0*(1.+XLIM)
      XXDIFF=(XFI-XIN)/2.
      DXIN=(XFI-XIN)/FLOAT(MTIMES)
      WIDTH=R0*FLUCT
      XXAV=(XFI+XIN)/2.
      SUM=0.
      V0X=V0
      R0X=R0
      DX=D
      QX=Q2
      ION1=ION
      IF(NOQ2.EQ.1)ION1=4
      XIN=XIN-DXIN
C-------------
      DO 300 MTIM=1,MTIMES
      XIN=XIN+DXIN
      XFI=XIN+DXIN
      XDIFF=(XFI-XIN)/2.
      XAV=(XFI+XIN)/2.
C------------------------------------------
      DO 200 J=JFIRST,JLAST
      IF(Z(J).NE.0.) GOTO  150
      R0X=XAV
      XXX0=XXAV-XAV
      GAU0=GAU(XXX0,WIDTH,XXDIFF)
      LMAX=0
      IF(LIMBAR.EQ.1) LMAX=-1
      CALL FUSION(ELAB,AP,ZP,AT,ZT,QX,V0X,R0X,DX,XTH,APUSH,FPUSH,
     1            CRED,NORUTH,NOPROX,
     2            NOCURV,ION1,SIGL,FLGRAZ,LMAX,SIGF,ITEST,IOPT)
      SUM=SUM+WEIGHT(J)*SIGF*GAU0
      LMAX2=MAX0(LMAX2,LMAX)
      DO 140 L=1,LMAX2
      SIGML(L)=SIGML(L)+
     1  WEIGHT(J)*SIGL(L)*GAU0
 140  CONTINUE
      GOTO 200
C------------------------
 150  XX1=Z(J)*XDIFF+XAV
      XXX1=XXAV-XX1
      GAU1=GAU(XXX1,WIDTH,XXDIFF)
      R0X=XX1
 160  LMAX=0
      IF(LIMBAR.EQ.1) LMAX=-1
      CALL FUSION(ELAB,AP,ZP,AT,ZT,QX,V0X,R0X,DX,XTH,APUSH,FPUSH,
     1            CRED,NORUTH,NOPROX,
     2            NOCURV,ION1,SIGL,FLGRAZ,LMAX,SIGF1,ITEST,IOPT)
      LMAX2=MAX0(LMAX2,LMAX)
      XX2=-Z(J)*XDIFF+XAV
      XXX2=XXAV-XX2
      GAU2=GAU(XXX2,WIDTH,XXDIFF)
      R0X=XX2
 170  LMAX=0
      IF(LIMBAR.EQ.1) LMAX=-1
      CALL FUSION(ELAB,AP,ZP,AT,ZT,QX,V0X,R0X,DX,XTH,APUSH,FPUSH,
     1            CRED,NORUTH,NOPROX,
     2            NOCURV,ION1,SIGL1,FLGRAZ,LMAX,SIGF2,ITEST,IOPT)
      SUM=SUM+WEIGHT(J)*(SIGF1*GAU1+SIGF2*GAU2)
      LMAX2=MAX0(LMAX2,LMAX)
      DO 180 L=1,LMAX2
      SIGML(L)=SIGML(L)+
     1  WEIGHT(J)*(SIGL(L)*GAU1+SIGL1(L)*GAU2)
 180  CONTINUE
 200  CONTINUE
C------------------------------------------
 300  CONTINUE
      SIGF=XDIFF*SUM
      LMAX2=MIN0(200,LMAX2)
      DO 210 L=1,LMAX2
 210  SIGML(L)=SIGML(L) *XDIFF
C--------------------------------------------------
      LMAX=LMAX2
      RETURN
      END
C---------------------------------------------------------------------
      FUNCTION GAU(X,SI,XLIM)
C     SPECIAL GAUSS
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(ABS(X).GT.XLIM .OR. XLIM.LE.0.) GOTO 100
      IF(ABS(SI).LT.0.01*XLIM)GOTO 80
      SI2=2.*SI
      DUM=X*X/(SI2  *SI)
      IF(DUM.GT.100.) GOTO 100
      T=XLIM/(1.41421*SI)
!-----The following changed by Ning Wang -----
!      T=DERF(T)
      T=1.-ERFC(T)
!---------------------------------------------
      DUM=EXP(-DUM)
      GAU=DUM/(2.5066*SI*T)
C     SQRT(2*PI)
      RETURN
 80   GAU=0.5/XLIM
 100  GAU=0.
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DOUBLE PRECISION FUNCTION ERFC (X)
      implicit real*8 (a-h,o-z)
      DOUBLE PRECISION S15ADF
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
      DOUBLE PRECISION       T,XHI,XLO,Y
      INTRINSIC                        ABS, EXP
      DATA XLO/-6.50D0/
      DATA XHI/9.5D0/
      IFAIL = 0
      IF (X.GE.XHI) GOTO 20
      IF (X.LE.XLO) GOTO 40
      T = 1.0D0 - 7.5D0/(DABS(X)+3.75D0)
      AY = ((((((+3.328130055126039D-10
     *    *T-5.718639670776992D-10)*T-4.066088879757269D-9)
     *    *T+7.532536116142436D-9)*T+3.026547320064576D-8)
     *    *T-7.043998994397452D-8)*T-1.822565715362025D-7)
      BY=((((((AY*T+6.575825478226343D-7)*T+7.478317101785790D-7)
     *    *T-6.182369348098529D-6)*T+3.584014089915968D-6)
     *    *T+4.789838226695987D-5)*T-1.524627476123466D-4)
      Y=(((BY*T-2.553523453642242D-5)*T+1.802962431316418D-3)
     *    *T-8.220621168415435D-3)*T+2.414322397093253D-2
      Y = (((((Y*T-5.480232669380236D-2)*T+1.026043120322792D-1)
     *    *T-1.635718955239687D-1)*T+2.260080669166197D-1)
     *    *T-2.734219314954260D-1)*T + 1.455897212750385D-1
C
      S15ADF =dexp(-X*X)*Y
      IF (X.LT.0.0D0) S15ADF = 2.0D0 - S15ADF
      ERFC = S15ADF
      RETURN
   20 S15ADF = 0.0D0
      ERFC = S15ADF
      RETURN
   40 S15ADF = 2.0D0
      ERFC = S15ADF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC