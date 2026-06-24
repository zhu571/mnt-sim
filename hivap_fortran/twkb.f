      FUNCTION TWKB(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,V0,R0,D,
     1           XTH,APUSH,FPUSH,DR,OCOS,ECM,L,RIN,ROUT,RFUS,CRED,IOPT)
C
C  Berechnet Transmissionskoeffizienten fuer Fusion mittels WKB -
C  naeherung ( Eingang WKB )
C  letzte Modifizierung 10.8.1987
C  Stand: 3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NPOINT(7),KEY(8),Z(24),WEIGHT(24)
      DATA NPOINT/2,3,4,5,6,10,15/,KEY/1,2,4,6,9,12,17,25/,
     1     Z/0.5773502  ,0.0000000  ,0.7745966  ,
     2       0.3399810  ,0.8611363  ,0.0000000  ,0.5384693  ,
     3       0.9061798  ,0.2386191  ,0.6612093  ,0.9324695  ,
     4       0.1488743  ,0.4333953  ,0.6794095  ,0.8650633  ,
     5       0.9739065  ,0.0000000  ,0.2011940  ,0.3941513  ,
     6       0.5709721  ,0.7244177  ,0.8482065  ,0.9372733  ,
     7       0.9879925  /,
     8WEIGHT/1.0000000  ,0.8888888  ,0.5555555  ,
     9       0.6521451  ,0.3478548  ,0.5688888  ,0.4786286  ,
     A       0.2369268  ,0.4679139  ,0.3607615  ,0.1713244  ,
     B       0.2955242  ,0.2692667  ,0.2190863  ,0.1494513  ,
     C       0.0666713  ,0.2025782  ,0.1984314  ,0.1861610  ,
     D       0.1662692  ,0.1395706  ,0.1071592  ,0.0703660  ,
     E       0.0307532  /
      DATA HBARC,AMU/197.329,931.502/
C-----------------------------------------------------------------
      AMXEXP=100.
      MTIMES=2
      MPNT=5
      DO  5 I=1,7
      IF(MPNT.EQ.NPOINT(I)) GOTO 6
 5    CONTINUE
      MPNT=5
      I=4
C     DEFAULT
 6    JFIRST=KEY(I)
      JLAST=KEY(I+1)-1
C------------------------------------------------------------------
      IF(DR.LE.0.)DR=0.2
      IUP=4./DR
      IUPP=8./DR
      IF(CRED.LE.0)CRED=1.
      ARED=AP*AT*CRED/(AP+AT)
      F0=2.*ARED*AMU/(HBARC*HBARC)
      F0=SQRT(F0)
      R=RIN
      V=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,
     1       XTH,APUSH,FPUSH,ECM,L,R,IOPT)
      GOTO 8
C-----------------------------------------------------------------
C      ENTRY WKB(OCOS,ECM,L,RIN,ROUT)
      R=RIN
      V=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,
     1       XTH,APUSH,FPUSH,ECM,L,RIN,IOPT)
 8    DO 10 I=1,IUP
      R1=R
      V1=V
      R=R-DR
      V=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,
     1       XTH,APUSH,FPUSH,ECM,L,R,IOPT)
      IF(V.LT.ECM) GOTO 12
 10   CONTINUE
      IF(L.LT.3)GOTO 200
      RI=DMAX1(RFUS,RIN)
      GOTO 25
 12   R3=0.5*(R+R1)
      V2=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,
     1       XTH,APUSH,FPUSH,ECM,L,R3,IOPT)
      IF(V2.LT.ECM)GO TO 15
      V1=V2
      R1=R3
      GOTO 20
 15   V=V2
      R=R3
 20   RI=((V1-ECM)*R-(V-ECM)*R1)/(V1-V)
 25   R=ROUT
      V=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,
     1       XTH,APUSH,FPUSH,ECM,L,ROUT,IOPT)
      DO 30 I=1,IUPP
      R1=R
      V1=V
      R=R+DR
      V=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,
     1       XTH,APUSH,FPUSH,ECM,L,R,IOPT)
      IF(V.LT.ECM) GOTO 35
 30   CONTINUE
      GOTO 201
 35   R3=0.5*(R+R1)
      V2=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,
     1       XTH,APUSH,FPUSH,ECM,L,R3,IOPT)
      IF(V2.LT.ECM)GO TO 45
      V1=V2
      R1=R3
      GOTO 50
 45   V=V2
      R=R3
 50   RO=((V1-ECM)*R-(V-ECM)*R1)/(V1-V)
      RIN  =RI
      ROUT =RO
      IF(RFUS.LT.RO)RI=DMAX1(RI,RFUS)
C----------------------------------------------------------
      DRR=(RO-RI)/FLOAT(MTIMES)
      R1=RI-DRR
      SY=0.
      DO 95 I=1,MTIMES
      R1=R1+DRR
      R2=R1+DRR
      C=(R2-R1)/2.
      D1=(R2+R1)/2.
      SUM=0.
      DO 90 J=JFIRST,JLAST
      IF(Z(J).NE.0.) GOTO 80
      DUM=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,
     1       XTH,APUSH,FPUSH,ECM,L,D1,IOPT)-ECM
      IF(DUM.GT.0.)DUM=SQRT(DUM)
      IF(DUM.LT.0.)DUM=0.
      SUM=SUM+WEIGHT(J)*DUM
      GOTO 90
 80   DUM=Z(J)*C+D1
      DUM1=-Z(J)*C+D1
      DUM=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,
     1       XTH,APUSH,FPUSH,ECM,L,DUM,IOPT)-ECM
      DUM1=VLRO(AP,AT,ZP,ZT,Q2,NORUTH,NOPROX,NOCURV,OCOS,V0,R0,D,
     1       XTH,APUSH,FPUSH,ECM,L,DUM1,IOPT)-ECM
      IF(DUM.GT.0.)DUM=SQRT(DUM)
      IF(DUM.LE.0.)DUM=0.
      IF(DUM1.GT.0.)DUM1=SQRT(DUM1)
      IF(DUM1.LE.0.)DUM1=0.
                     SUM=SUM+WEIGHT(J)*(DUM+DUM1)
 90   CONTINUE
      GAUSS=C*SUM*F0
      SY=SY+GAUSS
 95   CONTINUE
      Y=SY
      TWKB=0.
C----------------------------------------------------------
      DUM=2.*Y
      TWKB=0.
      IF(DUM.LT.AMXEXP)TWKB=EXP(-DUM)
      RETURN
C----------------------------------------------------------
 200  IF(L.LE.2) WRITE(6,210)  L,OCOS,RIN,R,V,ECM,IUP,DR
 210  FORMAT(' TWKB TROUBLE RI  L,OCOS=',I4,F8.3/' RIN,R,V,ECM',4F9.3/
     1       ' IUP,DR',I6,F8.2)
      TWKB=0.
      RETURN
 201  IF(L.LE.2) WRITE(6,220)  L,OCOS,ROUT,R,V,ECM
 220  FORMAT(' TWKB TROUBLE RO  L,OCOS=',I4,F8.3/' ROUT,R,V,ECM',4F9.3)
      TWKB=0.
      RETURN
      END