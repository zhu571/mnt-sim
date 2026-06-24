      SUBROUTINE FUSE(ELAB,AP,AT,ZP,ZT,SIGML,JDIM,IFUS,FLGRAZ)
C
C  Berechnet Fusionsquerschnitte, verschiedene Optionen, einschliess-
C  lich Bass-Modell,  modifiziert am 26.1.1980
C  Stand: 3.5.1994
C
C     PERFORMS SMOOTH CUT OFF OF UPPER PARTIAL WAVES (ENTRANCE)
C
C     IFUS=0  NO MODIFICATION OF SIGML
C     IFUS=1  FUSION GIVEN LCRIT AND DELTAL(DELTAL=0 IS SHARP CUTOFF)
C             WILL READ LCRIT,DELTAL
C     IFUS=2  READ RATIO SIGFUS TO SIG REACT
C     IFUS=3  READ LCRIT
C     IFUS=4  BASS1 CUT OFF
C     IFUS=5  FUSION GIVEN SIGFUS  AND DELTAL
C             (DELTAL=0 IS SHARP CUTOFF)
C             WILL READ SIGFUS,DELTAL
C     IFUS>5  NO MODIFICATION
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SIGML(JDIM)
      DIMENSION SIGFS (200) ,ARRAY(10)
C ------------------------------------------------------------------
      IF(IFUS.LE.0 .OR. IFUS.GT.5) RETURN
C --------------------------------------------------------------------
C
      IF(IFUS.EQ.2 .OR. IFUS .EQ. 3) THEN
         IF(FLGRAZ .LT.1.) RETURN
         LGRAZ=FLGRAZ+0.5
         DO 30 I=1,JDIM
         IF(SIGML(I))30,30,20
 20      LMAX=I-1
 30      CONTINUE
         WRITE(6,35) LMAX,FLGRAZ
 35      FORMAT(/' FUSE: VALUES FROM PARAP LMAX,LGRAZE',I10,F10.1)
         LMAX=MIN0(JDIM-1,LMAX)
         LMAX1=LMAX+1
      ENDIF
C --------------------------------------------------------------------
C
      ECM=ELAB*AT/(AT+AP)
      U=AT*AP/(AT+AP)
      PWAVE=SQRT(2.*931.16*U*ECM)
      WAVEL=19.732/PWAVE
      ARWAVE=3141.59*(WAVEL**2.)
      GOTO(50,80,160,180,60),IFUS
C----------------------------------------------------------------
C     FUSION GIVEN LCRIT,DELTAL(WELL ABOVE CB)                 IFUS=1
C
 50   CALL MYRD(ARRAY,2,23,5)
      FLCRIT=ARRAY(1)
      DELTAL=ARRAY(2)
      LUP =FLCRIT+5.*DELTAL+1.5
      LUP =MIN(JDIM,LUP )
      SIGFUS=0.
      DO 55 L=1,LUP
      AL=L-1
      TL=1.
      IF(DELTAL.LE.0)GOTO 53
      DUM=(AL-FLCRIT)/DELTAL
      TL=1./(1.+EXP( DUM))
 53   SIGML(L)=ARWAVE*TL*(2.*AL+1.)
      SIGFUS=SIGFUS+SIGML(L)
 55   CONTINUE
      IF (LUP.LT.JDIM) THEN
          DO 57 L=LUP+1,JDIM
 57       SIGML(L)=0.
      ENDIF
      LMAX=LUP-1
      WRITE(6,56) ELAB,FLCRIT,DELTAL,LMAX,SIGFUS
 56   FORMAT(' FUSE OPTION 1',3X,'ELAB=',F8.2,3X,'FLCRIT,DELTAL',2F8.2,
     1       3X,'LMAX=',I3,3X,'SIGFUS(MB)=',E12.4/)
      RETURN
C ------------------------------------------------------------------
C                             GIVEN SIGFUS,DELTAL               IFUS=5
 60   CALL MYRD(ARRAY,2,23,5)
      SIGFUS=ARRAY(1)
      DELTAL=ARRAY(2)
C
      FLCRIT=0.
      DO 61 L=1,JDIM
 61   SIGML(L)=0.
      DO 62 I=1,3
 62   FLCRIT=SQRT(SIGFUS/ARWAVE-FLCRIT)
      IREP=0
 63   SIGF=0.
      LUP =FLCRIT+5.*DELTAL+1.5
      LUP =MIN(JDIM,LUP )
      DO 65 L=1,LUP
      AL=L-1
      TL=1.
      IF(DELTAL.GT.0.01) THEN
          DUM=(AL-FLCRIT)/DELTAL
          IF(DUM.GE.150.) TL=0.
          IF(ABS(DUM).LT.150.) TL=1./(1.+EXP( DUM))
      ENDIF
      SIGML(L)=ARWAVE*TL*(2.*AL+1.)
      SIGF=SIGF+SIGML(L)
 65   CONTINUE
      IF(DELTAL.LE.0.) GOTO 68
      IF(IREP.GT.1) GOTO 74
      DSIG=SIGFUS-SIGF
      DUM=0.01*SIGFUS
      IF(ABS(DSIG).LT.DUM)GOTO 74
      IREP=IREP+1
      DUM=     SIGF /(FLCRIT*FLCRIT)
      DLCRIT=DSIG/(2.*DUM*FLCRIT)
      FLCRIT=FLCRIT+DLCRIT
      GOTO 63
 68   SIGF=SIGFUS-DSIG
      IF(DSIG.GT.0.) GOTO 69
      SIGML(LUP)=SIGML(LUP)+DSIG
      GOTO 74
 69   LUP=LUP+1
      SIGML(LUP)=DSIG
      FLCRIT=LUP-1
 74   LMAX=LUP-1
      WRITE(6,75) ELAB,SIGFUS,DELTAL,SIGF,FLCRIT,LMAX
 75   FORMAT(/' FUSE OPTION 5  ELAB(MEV)=',F8.2/
     1       ' ENTERED SIGFUS,DELTAL:',E12.4,F8.3,5X,
     2       ' OUT SIGFUS,FLCRIT:',E12.4,F8.2,3X,'LMAX=',I3)
      RETURN
C ------------------------------------------------------------------
C            GIVEN RATIO( 0 T0 1) COMPARED TO PARAP         IFUS=2
 80   CALL MYRD(ARRAY,1,23,5)
      RATIO=ARRAY(1)
 81   FORMAT( 9F8.3)
 90   IF(RATIO.LE.0. .OR. RATIO .GE.1.) RETURN
      FLCRIT=SQRT(RATIO)*FLGRAZ
      GOTO 195
C ------------------------------------------------------------------
C            SMALLER LCRIT        COMPARED TO PARAP         IFUS=3
 160  CALL MYRD(ARRAY,1,23,5)
      FLCRIT=ARRAY(1)
      GOTO 195
C-------------------------------------------------------------------
C     BASS MODEL  NPA231(1974)45  EQUATION NR IN COL 72-80  IFUS=4
C
 180  AS=17.
      RZERO=1.07
      FF=5./7.
      D=1.35
      AP3=AP**0.33333
      AT3=AT**0.33333
      SUMA3=AP3+AT3
      R12=RZERO*SUMA3
      DUM=AP3*AT3*(AP3+AT3)
      X=0.0792*ZP*ZT/DUM
      DUM=SUMA3*SUMA3*AP*AT*AT3*AP3
      Y=1.065*(AP+AT)/DUM
      VR=AS*AP3*AT3*D/R12
      F=VR
      VCOUL=1.44*ZT*ZP/R12
      E1=VCOUL*(1.+0.5*(1.-X)/X-D/(X*R12))
      E1LAB=E1*(AT+AP)/AT
      CFUE1=SQRT((1.-X)/(2.*Y))
      E2=VCOUL*(1.+0.5*(1.-X)/(X*FF*FF)-D/(X*R12))
      E2LAB=E2*(AT+AP)/AT
      CFUE2=CFUE1/FF
      PRINT 170
 170  FORMAT(/ 2X,'BASS MODEL' )
      PRINT 181,E1,E1LAB,CFUE1,E2,E2LAB,CFUE2
 181  FORMAT(1X ,'E1,E1LAB,LFUE1',3F10.3/1H0,'E2,E2LAB,LFUE2',3F10.3 )
C
C     FOR ECM.LT.E1 USE INVERTED PARABOLA(INVPAR) FOR LCRIT
C     FOR E1.LT.ECM.LT.E2   LCRIT IS LINEAR INTERPOLATION BETWEEN CFUE1
C     FOR ECM.GT.E2  LCRIT=CFUE2=CONSTANT
      IF(ELAB-E1LAB)183,183,185
 183  CALL INVPAR(ELAB,F,RZERO,D,AP,ZP,AT,ZT,SIGFS ,FLCRIT)
      DO 184 L=1,200
 184  SIGML(L)=0.
      DO 187 L=1,199
      SIGML(L)=SIGFS (L)
      IF(SIGFS (L+1).LE..1E-3) GOTO 188
 187  CONTINUE
 188  LUP=L
      LUP1=L+1
      GOTO 218
 185  IF(ELAB-E2LAB) 186,186,190
 186  ECM=ELAB*AT/(AT+AP)
      FLCRIT=CFUE1*CFUE1      +(CFUE2*CFUE2-CFUE1*CFUE1)*(ECM-E1)
     1 /(E2-E1)
      FLCRIT=SQRT(FLCRIT)
      GOTO 195
 190  FLCRIT=CFUE2
 195  LCRIT=FLCRIT+0.5
 200  FDL=FLGRAZ-FLCRIT
      IDL=FDL
      IF(FDL.GT.0.) GOTO 201
      IDL=-FDL+1.
      IDL=-IDL
 201  IF(IFUS.EQ.3) GOTO 202
      IF(FDL.LE.0.) GOTO 400
 202  CONTINUE
C
C     SMOOTH CUT OFF FUSION PARTIAL X-SECTIONS
C
      DO 205 L=1,LMAX1,10
      L1=MIN0(L+9,LMAX)
      WRITE(6,225) (SIGML(L2),L2=L,L1)
 205  CONTINUE
C
      LUP=LMAX-IDL
      LUP1=LUP+1
      LDN1=MAX0(1,1-IDL)
      DO 210 L=LDN1,LUP1
      AL=L
      L1=AL+FDL
      FL1=L1
      X=AL+FDL-FL1
      FL=L1+L1-1
      FL1=L1+L1+1
      DUM=(1.-X)*SIGML(L1)/FL+X*SIGML(L1+1)/FL1
      FL=L+L-1
      SIGFS (L)=DUM*FL
 210  CONTINUE
      DO 213 L=LDN1,LUP1
 213  SIGML(L)=SIGFS (L)
        LUP1=LUP1+1
      IF(LUP1.GT.LMAX1) GOTO 217
 214  DO 215 L=LUP1,LMAX1
 215  SIGML(L)=0.
 217  CONTINUE
C
 218  CONTINUE
      WRITE(6,220) LMAX,FLGRAZ,LUP,FLCRIT,ELAB,RATIO
 220  FORMAT(//2X,'LMAX,LGRAZ,LUP,LCRIT',I8,F8.1,I8,F8.1/
     1 1H0,'ELAB,RATIO',2F10.3/2X,'FUSION PARTIAL X-SECTIONS')
      DO 230 L=1,LUP1,10
      L1=MIN0(L+9,LUP1)
      WRITE(6,225) (SIGML(L2),L2=L,L1)
 225  FORMAT(1X,10F11.3)
      FLGRAZ=FLCRIT
 230  CONTINUE
 400  CONTINUE
      RETURN
      END