      SUBROUTINE TLD(H4,N1,J1,H1,Z,W1,STPLTH)
C
C  Berechnet Transmissionskoeffizienten in 'Double Precision' nach dem
C  optischen Modell, T(2*S+1,LMAX)
C  wird von OVER2 und OM aufgerufen
C  Stand: 3.5.1994
C
C-------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SG(2),S(16,4),F(101,2),G(101,3),R(2,3),H(502),A(502),B(5
     102),Q(502),D(300),E(300),Y(300),U(300),Z1(5),U1(7),Y1(7)
      COMMON/TLJ/ T(2,31),V(15),V1(3)
      DOUBLE PRECISION T,V,V1,H4,H1,Z,W1,STPLTH
C     H4=0.0478326*AREDUCED*ELAB*AT/(AP+AT)=0.04783*AREDUCED*ECM=K**2
C     W1=0.0478326*AREDUCED
C     N1 OUTPUT,ESSENTIALLY ROUNDED UP VALUE OF V1(3)
C     J1 OUTPUT (2S+1 OF PROJECTILE)
C     H1=SQRT(H4)=K=WAVENUMBER=1/LAMBDA_SLASH
C     Z=0.03478*ZP*ZT*ARED/K
C     STPLTH  INTEGRATION STEPS IN FERMI (0.1 IS TYPICAL)
C     V1(1)=  INTRINSIC SPIN OF PROJECTILE (S)
C     V1(2)=  INTRINSIC SPIN OF TARGET PRESUMABLY
C     V1(3)=  LMAX MAX ORBITAL ANG MOM+1
C     V(1 )=  REAL POT RADIUS
C     V(2 )=  REAL POT DIFFUSENESS
C     V(3 )=  IMAG POT RADIUS
C     V(4 )=  IMAG POT DIFFUSENESS
C     V(5 )=  1 SURFACE IMAG, 0 VOLUME IMAG
C     V(6 )=  DEPTH SPIN ORBIT POT
C     V(7 )=  DEPTH REAL POT
C     V(8 )=  DEPTH IMAG POT
C     V(9 )=  COULOMB RADIUS
C     V(10)=  SPIN ORBIT POT RADIUS
C     V(11)=  SPIN ORBIT POT DIFFUSENESS
C-------------------------------------------------------------------
C      DOUBLE PRECISION A,B,D,E,F,G,H,Q,R,S,U,Y,A1,A2,B1,B7,DT,D1,HY,HZ,
C     1                 H2,P5,P6,P7,P8,P9,RB,SG,TM,T1,T2,T3,T4,T5,T6,T7,
C     2                 T8,T9,U1,U2,V4,V5,W2,Y1,Y2,Z1,Z2,DTT,DT3,HZ3,
C     3                 T101,T102,TT1,TT2
C-----------------------------------------------------------------------
      B1 = 0.
 1    I1=2.0*V1(1) + 2
      I3=I1+2
      I4=I1*(I1-2)
      J1=I1-1
      B7=V(1)
      V(1)=DABS(V(1))
      RB=V(1)+9.*V(2)
      RB=DMAX1(RB,V(10)+9.0*V(11))
      T1=V(3)+9.*V(4)
      IF(RB-T1)2,3,3
 2    RB=T1
 3    IRB=10.*RB
      RB=IRB
      RB=RB/10.
      M4=RB/STPLTH
      IF(M4-299)5,5,4
 4    M4=299
 5    D1=RB/FLOAT(M4)
      M4=M4+1
      M1=M4+1
      M3=M4+2
      N1=(V1(3)+1.0)
      V1(3)=FLOAT(N1-1)
      Z2=Z*Z
      DT=2.*D1
      DTT=D1*D1/12.
      W2=W1
      W1=W1*DTT
      DT3=-2.*DTT
      H2=1.+DTT*H4
      HY=DT3*H1*Z
      S(16,1)=1.
      S(16,2)=1.
      S(16,3)=0.
      S(16,4)=0.
      T1=1.
      T2=2.
      T3=0.
C------------------------
      DO 6 I=1,15
      J=16-I
      T7=T1*Z/T2
      T6=-T3*(T3+1.)+Z2
      T4=T6/T2
      T5=(2.+T6)/T2
      S(J,1)=T7*S(J+1,1)-T4*S(J+1,3)
      S(J,2)=T7*S(J+1,2)-T5*S(J+1,4)
      S(J,3)=T7*S(J+1,3)+T4*S(J+1,1)
      S(J,4)=T7*S(J+1,4)+T5*S(J+1,2)
      T1=T1+2.
      T2=T2+2.
 6    T3=T3+1.
C------------------------
      TM=RB
      IF(Z)7,10,7
 7    DO 9 I=1,4
      T1=( DEXP((13.816+DLOG( DABS(S(1,I))))/15.0))/H1
      IF(T1-TM)9,9,8
 8    TM=T1
 9    CONTINUE
 10   M2=(TM-RB)/D1 + 1
      TM=RB+D1*FLOAT(M2)
      T1=Z2+16.
      SG(1)=-Z+Z*(DLOG(T1))/2.0+3.5*DATAN(Z/4.0)-DATAN(Z)
     1 -DATAN(Z/2.0)-DATAN(Z/3.)
     2      -Z*(1.+(Z2-48.)/(30.*T1*T1)+(Z2*Z2-160.*Z2+1280.)/(105.*T1*
     3T1*T1*T1))/(12.*T1)
      SG(2)=SG(1)-1.5707963+DATAN(Z)
      T1=TM+DT
      DO 13 I=1,2
      T1=T1-D1
      T2=T1*H1
      T3=T2-Z*DLOG(2.0*T2)
      DO 12 J=1,2
      T7=0.
      T8=0.
      DO 11 K=1,15
      T7=(T7+S(K,J))/T2
 11   T8=(T8+S(K,J+2))/T2
      T7=T7+1.
      T4=T3+SG(J)
      T5=DCOS(T4)
      T6=DSIN(T4)
 12   G(J,I)=T7*T5-T8*T6
      R(1,I)=H2+HY/T1
 13   R(2,I)=R(1,I)+DT3/(T1*T1)
C-----------------------------------
      DO 14 I=1,M2
      T1=T1-D1
      R(1,3)=H2+HY/T1
      R(2,3)=R(1,3)+DT3/(T1*T1)
      DO 14 J=1,2
      G(J,3)=((12.-10.*R(J,2))*G(J,2)-R(J,1)*G(J,1))/R(J,3)
      DO 14 K=1,2
      R(J,K)=R(J,K+1)
 14   G(J,K)=G(J,K+1)
C-----------------------------------
      T1=T1+D1
      DO 26 I=1,2
      T2=T1*H1
      T3=1./T2
      T1=T1-D1
      T7=0.
      N4=0
      N2=(T2/1.4142)*DSQRT(25.-2.*Z*T3+10.*DSQRT((Z*T3-.5)**2+6.))
      IF(N2-N1-8)15,16,16
 15   N2=N1+8
 16   IF (N2-500) 17,17,46
 17   T6=FLOAT(N2)
      N3=N2+1
      T5=T7
      H(N3+1)=0.
      H(N3)=1.0E-20
      A(N3)=DSQRT(Z2+(T6+1.0)**2)/(T6+1.0)
      DO 21 K=1,N2
      M=N3-K
      IF(N4-M)18,19,19
 18   A(M)=DSQRT(Z2+T6*T6)/T6
      B(M)=(2.*T6+1.)*(Z/(T6*(T6+1.))+T3)
      T6=T6-1.
 19   H(M)=(B(M)*H(M+1)-A(M+1)*H(M+2))/A(M)
      IF(DABS(H(M+2))-(10.0**30)) 21,21,20
 20   H(M)=H(M)/(10.0**25)
      H(M+1)=H(M+1)/(10.0**25)
 21   CONTINUE
      N4=N2
      N2=N2+10
      T7=H(2)/H(1)
      IF(DABS((T5-T7)/T7)-0.0001) 23,23,22
 22   IF(N2.LE.500) GO TO 16
      N2 = 500
 23   T5=1./(A(1)*(H(1)*G(2,I)-H(2)*G(1,I)))
      F(1,I)=T5*H(1)
      F(2,I)=T5*H(2)
      DO 26 K=3,N1
      IF(DABS(H(K)/H(K-1))-10.0**15) 25,25,24
 24   T5=T5/(10.0**25)
 25   F(K,I)=T5*H(K)
 26   G(K,I)=(B(K-2)*G(K-1,I)-A(K-2)*G(K-2,I))/A(K-1)
C-----------------------------------
      Q(1)=0.
      Q(2)=0.
      H(1)=0.
      H(2)=0.
      A(1)=0.
      A(2)=1.0E-20
      B(1)=0.
      B(2)=1.0E-20
      HZ3=H2+1.5*HY/V(9)
      HZ=-HY/(2.*(V(9)**3))
      T101=1.0/DEXP(V(10)/V(11))
      T1=1.0/DEXP(V(1)/V(2))
      T102=DEXP(D1/V(11))
      T2=DEXP(D1/V(2))
      IF(B7)27,28,28
 27   A2=V(5)*V(5)/16.
      V(5)=1.
      A1=4.*W2*A2
      B1=DEXP(4.0*A2*H4)
 28   V5=V(5)
      IF(V5)29,30,30
 29   V(5)=-V5
      V4=V(4)
      V(4)=0.69*V(4)
 30   T9=V(8)*(1.-V(5))
      T3=1.0/DEXP(V(3)/V(4))
      T4=DEXP(D1/V(4))
      T6=0.
C-----------------------------------------------------
      DO 40 K=1,M4
      T6=T6+D1
      Y(K)=T6*T6
      T1=T1*T2
      T101=T101*T102
      T3=T3*T4
      T5=V(7)/(1.+T1)
      T7=1./(1.+T3)
      IF(V5)31,34,34
 31   T8=((T6-V(3))/V4)**2
      IF(T8-10.)33,32,32
 32   T8=0.
      GO TO 37
 33   T8=V(8)*DEXP(-T8)
      GO TO 37
 34   T8=4.*T3*T7*T7*V(8)
      IF(B7)35,37,37
 35   P5=A2/(T5*T5+T8*T8)
      P8=T1/(1.+T1)
      P9=T3*T7
      P6=-T5*P8/V(2)
      P7=T8*(1.-2.*P9)/V(4)
      P8=P6*(1.-2.*P8)/V(2)
      P9=T8*(1.-6.*P9*(1.-P9))/(V(4)*V(4))
      U2=P5*((T5*P6+T8*P7)*2./T6+T5*P8+T8*P9)
      Y2=P5*((T5*P7-T8*P6)*2./T6+T5*P9-T8*P8)
      Y1(1)=T8/(B1*B1+2.*B1*A1*T5)
      U1(1)=(T5+T8*A1*Y1(1))/(B1+A1*T5)
      DO 36 J=1,6
      P5=A1*Y1(J)-Y2
      P6=DSIN(P5)
      P5=DCOS(P5)
      P7=1.0/(B1*DEXP(A1*U1(J)-U2))
      U1(J+1)=(T5*P5+T8*P6)*P7
 36   Y1(J+1)=(T8*P5-T5*P6)*P7
      T5=U1(7)-((U1(7)-U1(6))**2)/(U1(7)-2.*U1(6)+U1(5))
      T8=Y1(7)-((Y1(7)-Y1(6))**2)/(Y1(7)-2.*Y1(6)+Y1(5))
 37   H(K+2)=W1*(T9*T7+V(5)*T8)
      U(K)=2.*W1*V(6)*T101/(V(11)*T6*((1.+T101)**2))
      IF(T6-V(1))38,39,39
 38   E(K)=W1*T5+HZ3+HZ*Y(K)
      GO TO 40
 39   E(K)=W1*T5+H2+HY/T6
 40   CONTINUE
      T8=0.
C------------------------------------------------------------------
C                                            I ORBIT AM LOOP
      DO 45 I=1,N1
      I2=2*I
      T1=DTT*T8*(T8+1.)
      T8=T8+1.
      DO 41 K=1,M4
 41   D(K)=E(K)-T1/Y(K)
C--------------------------------------------
C                                J INTRINSIC SPIN LOOP
      DO 45 J=1,J1
      L=I2-I3+2*J
      IF(IABS(I2-I1)-L)43,43,42
 42   T(J,I)=0.
      GO TO 45
 43   T9=(FLOAT(L*(L+2)-I2*(I2-2)-I4))/4.
      DO 44 K=1,M4
      Q(K+2)=D(K)+T9*U(K)
      T3=12.-10.*Q(K+1)
      T4=10.*H(K+1)
      T1=T3*A(K+1)-Q(K)*A(K)+T4*B(K+1)+H(K)*B(K)
      T2=T3*B(K+1)-Q(K)*B(K)-T4*A(K+1)-H(K)*A(K)
      T3=Q(K+2)**2+H(K+2)**2
      A(K+2)=(Q(K+2)*T1+H(K+2)*T2)/T3
 44   B(K+2)=(Q(K+2)*T2-H(K+2)*T1)/T3
      T3=A(M3)**2+B(M3)**2
      T1=(A(M1)*A(M3)+B(M1)*B(M3))/T3
      T2=(A(M3)*B(M1)-A(M1)*B(M3))/T3
      T5=F(I,2)-F(I,1)*T1
      T6=F(I,1)*T2
      T3=T5-G(I,1)*T2
      T4=G(I,2)-G(I,1)*T1+T6
      T7 = (T3/T4)*T3+ T4
      T1=(T3*T5+T4*T6)/T7
      T2=(T4*T5-T3*T6)/T7
      T7 = T4
      T1 = T1/T7
      TT1=T1
      IF(T1.LT.0.1D-15)TT1=0.
      T2 = T2/T7
      TT2=T2
      IF(T2.LT.0.1D-15)TT2=0.
C     T(J,I)=4.*(T1-T1**2-T2**2)
      T(J,I)=4.*(T1-TT1*TT1-TT2*TT2)
 45   CONTINUE
C----------------------------------------     END ORBIT AM LOOP
C--------------------------------------------------------------------
 46   CONTINUE
      RETURN
      END