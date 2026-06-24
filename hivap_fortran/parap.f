      SUBROUTINE PARAP(ELAB,AP,AT,ZP,ZT,FINERT,RATIOS,IFUS,LIMBAR,JCAL,
     1                NOPRNT,  MX,SIGML,CRSROT,FLGRAZ)
C
C  Berechnet Fusionsquerschnitte mit Parabelmethode (Hill-Wheeler)
C  Stand: 3.5.1994
C
C     QUICK CALCULATION OF NUCLEAR CROSS SECTIONS VIA OPTICAL MODEL
C     WITH PARABOLIC BARRIER - SEE THOMAS, PHYS REV 116,703 (1959) -
C
C     STANDARD V0,R0,D = 67.,1.17,0.574 (BLANN-PLASIL)
C     IFUS=5   V0,R0,D = 70.,1.25,0.44 (SIKKELAND)
C     IFUS=7     READS V0,R0,D FROM CARDS
C     LIMBAR=1        LMAX LIMITED BY      FISSION BARRIER>0.1MEV
C     OUTPUT MX IS EROT+0.5 FOR MAX AM
C     FLGRAZ OUTPUT:GRAZING AM,
C     NOPRNT=1 REDUCED PRINT OUT
C     JCAL=0 STANDARD
C     JCAL=2 SIMPLIFIED VERSION :NO AM BUT EFFECTIVE EXCIT (EXCIT-EROT)
C            WHERE EROT=YRAST FROM MODIFIED LDR
C     JCAL=3 SAME AS '2' BUT EROT=SPHERIC NUCLEUS(RIGID)
C     FINERT,RATIOS(1-6) PARAMETERS FOR MODIFIED LDR
C     SIGML(L) PARTIAL CROSS SECTIONS
C     CRSROT(EROT+0.5) SUM OF SIGML CORESPONDING TO ROT.ENERGY AROUND
C                      EROT(1MEV INTERVALS),NEEDED FOR JCAL=2,3
C-----------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIGML(1),CRSROT(1),RATIOS(6),ARRAY(4)
      COMMON/FISRT/ROT0,A2MS,RKAPPA,X0,Y0
C
      LMAX=200
      DO1 I=1,LMAX
      CRSROT(I)=0.
 1    SIGML(I)=0.
 2    CONTINUE
C------------------------------------- POTENTIAL PARAMETERS
      F=67.
      IF(AP.LE.4.)F=1100.
      R00=1.17
      D=0.574
      IF(IFUS.NE.5) GOTO 3
      F=70.
      R00=1.25
      D=0.44
C3    IF(IFUS.EQ.7) READ(5,1015) F,R00,D
 3    IF(IFUS.NE.7) GOTO 4
      CALL MYRD(ARRAY,3,23,5)
      F=ARRAY(1)
      R00=ARRAY(2)
      D=ARRAY(3)
 4    WRITE(6,1016) F,R00,D,ELAB
C---------------------------------------------------------------
      ACN=AP+AT
      E  =ELAB*AT/ACN
      AR=AT**0.333333+AP**0.333333
      IF(AP.LE.4.) AR=AT**0.333333
      COUL= 1.4393 *ZT*ZP
      RA=R00*AR
      U=AT*AP/(AT+AP)
      H2=41.814
      CROTL=H2/(1.16*ACN**1.666666)
      CR2=D*COUL/F
      CR3P=H2*D/(F*U)
      PWAVE=SQRT(2.*931.16*U*E)
      WAVEL=19.732/PWAVE
      ARWAVE=3141.59*(WAVEL**2)
C-------------------------------------------------------------
 8    IF(NOPRNT.EQ.1) GOTO 80
      WRITE(6,1002)
 1002 FORMAT(3H0  )
      WRITE(6,1003)
      WRITE(6,1004)
 1003 FORMAT(61H0 PROJECTILE  TARGET   PROJECTILE  TARGET    ENERGY    E
     1NERGY  )
 1004 FORMAT('    MASS       MASS      CHARGE    CHARGE     LAB       '
     1	   //' CM    ')
      WRITE(6,1005)
 1005 FORMAT(3H   )
      WRITE(6,1006) AP,AT,ZP,ZT,ELAB,E
 1006 FORMAT(4X,F7.3,5(3X,F7.3))
      WRITE(6,1002)
      WRITE(6,1007)
 1007 FORMAT(124H ANGULAR TRANSMISSION BARRIER BAR-RADIUS  BARRIER    BA
     1RRIER   RELATIVE CURVATURE  PROBABLE SIGMA(L) ROTOR-E EFFROT EFFSI
     2GMA    )
      WRITE(6,1008)
 1008 FORMAT(124H MOMENTA COEFFICIENTS  RADII   CONSTANT    (MEV)    CUR
     1VATURE  SIGMA(L)   SIGN     L-ENERGY  (MB)     (MEV)   (MEV)  (MB)
     2       )
      WRITE(6,1005)
      WRITE(6,1009)
 1009 FORMAT('     L         TL        RB        R0       VLB        '//
     1'WL       RSIGL    SD2VLB     SIGXL2   SIGL     EROTL  EROTLI '//
     2'SUMDXL    ')
      WRITE(6,1005)
C---------------------------------------------------------------------
C
 80   SUMAL2=0.0
      SUMTL=0.0
      SUMDXL=0.0
      EROTLI=0.0
      IREP=0
      TL2=0.
C                                      PARTIAL WAVE LOOP
      DO 29 K=1,LMAX
      L=K-1
      R=RA
      AL=FLOAT(L)
      CR3=CR3P*AL*(AL+1)
 9    RB=RA-D*DLOG(CR2/R**2+CR3/R**3)
      DELTR=ABS(RB-R)
      IF(DELTR-0.001) 11,10,10
 10   R=RB
      GO TO 9
 11   CONTINUE
      R0=RB/AR
      VLB=COUL/RB+(H2/2./U)*(AL*(AL+1)/RB**2)-F*EXP(-(RB-RA)/D)
      D2VLB=2.*COUL/RB**3+(H2/U)*3.*AL*(AL+1)/RB**4
     1 -(F/D**2)*EXP(-(RB-RA)/D)
      IF(D2VLB) 12,13,14
 12   SD2VLB=-1.0
      GO TO 15
 13   SD2VLB=0.0
      GO TO 15
 14   SD2VLB=1.0
 15   CONTINUE
      WL2=ABS((H2/U)*D2VLB)
      WL=SQRT(WL2)
      DUM=2.*3.1459*(VLB-E)/WL
      IF(DUM.LT.-20.) GOTO 152
      TL=1./(1.+EXP(DUM))
      GOTO 155
 152  TL=1.
 155  CONTINUE
      RSIGL=TL*(2.*AL+1.)
      SIGL=ARWAVE*RSIGL
      SIGML(K)=SIGL
      IF(TL.LT.0.5) IREP=IREP+1
      IF(IREP.GE.2) GOTO 160
      FLGRAZ=L-1
      TL1=TL2
      TL2=TL
      IF(L.EQ.0) TL0=TL
 160  CONTINUE
C-----------------------------------------------
      A=ACN
      Z=ZP+ZT
      AN=A-Z
      ONLYGS=0.
      CALL FISROT(A,Z,AN,AL,DELR,DELSP,ERO,ONLYGS)
      BAR=DELSP-DELR
C     WRITE(6,161)BAR,DELSP,DELR
C 161 FORMAT('Bar=',F8.3,' DELSP=',F8.3,' DELR=',F8.3)
      IF(LIMBAR.EQ.0) GOTO 170
      IF(BAR.GT.0.1) GOTO 170
      FLGRAZ=AL
      GOTO 30
 170  CONTINUE
C----------------------------------------------EROTL
      EROTL=ERO
      IF(JCAL.EQ.3) GOTO 19
      RIGITY=FINERT
C     IF(AL.GE.RATIOS(6))GOTO 18
C     CALL IP2(RATIOS(2),RATIOS(4),RATIOS(6),AL,RATIOS(1),RATIOS(3),
C    1         RATIOS(5),RIGITY)
 18   EROTL=DELR*(AL+1.)/RIGITY
      IF(AL.GT.0.)EROTL=EROTL/AL
C----------------------------------------------
 19   SIGXL2=AL*(AL+1.)*RSIGL
      IF(EROTL-0.5-EROTLI) 20,20,22
 20   SUMDXL=SUMDXL+SIGL
      IF(TL-0.0001    ) 22,21,21
 21   CONTINUE
      GO TO 25
 22   IF(L) 24,24,23
 23   CONTINUE
C23   WRITE(6,1010) EROTLI,SUMDXL
      LX=EROTL+.5
      CRSROT(LX)=SUMDXL
      MX=LX
 1010 FORMAT(3H+  ,107X,F5.1,1X,F7.2)
 24   CONTINUE
      EROTLI=EROTLI+1.0
      SUMDXL=SIGL
 25   CONTINUE
      IF(TL.GT.0.99) GOTO 28
C----------------------------------------
      IF(NOPRNT.EQ.1) GOTO 28
      WRITE(6,1011) L,TL,RB,R0,VLB,WL,RSIGL,SD2VLB,SIGXL2,SIGL,EROTL
 1011 FORMAT(3X,I3,4X,F10.8,2X,F6.3,4X,F6.3,4X,F6.1,4X,F6.3,4X,F6.2,4X,F
     16.3,1X,F10.2,2X,F7.3,3X,F7.3)
      IF(L) 26,26,28
 26   IF(TL-.0001)27,28,28
 27   WRITE(6,1010) EROTLI,SUMDXL
C---------------------------------------
 28   CONTINUE
      SUMTL=SUMTL+RSIGL
      SUMAL2=SUMAL2+SIGXL2
      IF(TL-.0001)30,29,29
 29   CONTINUE
C
C                                      END PARTIAL WAVE LOOP
C---------------------------------------------------------------------
 30   SIGMA=ARWAVE*SUMTL
      IF(LIMBAR.EQ.1) GOTO 40
      IF(TL0.LT.0.5) GOTO 32
      FLGRAZ=FLGRAZ+(TL1-0.5)/(TL1-TL2)
 32   IF(FLGRAZ.LT.0.)FLGRAZ=0.
      WRITE(6,1012) SIGMA,FLGRAZ
 1012 FORMAT(' PARAP X-SECTION(MBARN) =',  F8.1,4X,'LGRAZE',F8.1 /)
      AVEL2=SUMAL2/SUMTL
      AVROE=AVEL2*CROTL
      RMSL2=SQRT(AVEL2)
      IF(NOPRNT.EQ.0)WRITE(6,1013) RMSL2,AVEL2,AVROE
      GOTO 50
 40   WRITE(6,1017) SIGMA,FLGRAZ
 50   CONTINUE
 1013 FORMAT(23H0ROOT MEAN SQUARE L =  ,F6.2,10X,16H RMSL SQUARED = ,F10
     1.3,10X,33H AVERAGE ROTATION ENERGY (MEV) = ,F6.3)
 1014 FORMAT(6E10.3)
 1015 FORMAT(10F8.3)
 1016 FORMAT(/' PARAP',6X,'V0=',F8.3,6X,'R0=',F8.3,6X,'D=',F8.3,6X,
     1 'ELAB=',F8.1)
 1017 FORMAT(' ZERO FISSION BARRIER LIMITED  SIGFUS,LCRIT',F10.3,F10.0)
      RETURN
      END