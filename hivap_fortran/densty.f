      SUBROUTINE  DENSTY (NROW, NSTEP,LDBM,LEP1,LEVPRN,KHS)
C
C  Berechnet die Zustandsdichten, wird von HIVAP und EVA aufgerufen
C  Stand: 3.5.1994
C
      INCLUDE 'common.f'
      DIMENSION ARRAY(300),YJK(300),XLDM(9),THMIN(9)
C     DIMENSION OF ARRAY MUST BE >= JDIM
      DATA XLDM,THMIN/.3 ,.4  , .5, .6 ,.67, .7 , .8, .9 , 1.,
     1                .57,.555,.55,.495,.47,.485,.62,.795,1.0/
C--------------------------------------------------------------------
C     MINIMUM MOMENT OF INERTIA FROM CPS74
      AMPAR1=0.0
      AMPAR2=30.
      AMPAR3=2.0
      AMPAR4=0.3
      AMPAR5=10.
      KJC1=4*JDIM +1
      PINERF=(YRSMO(KJC1+30)-YRSMO(KJC1))/(31.*30.)
      PINERT=(YRSMO(31)-YRSMO(1))/(31.*30.)
      IF (PINERF.LT.0. .OR. PINERT.LT.0.)WRITE(6,210)YRSMO(1),YRSMO(2),
     1YRSMO(KJC1),YRSMO(KJC1+1),NROW,NSTEP
 210  FORMAT(' DENSTY YRSMO ERROR',4E12.4,2I6)
C     ABOVE IS HBAR-SQUARE OVER INERTIA FOR SADDLE AND GS
      BETAF=(PINERT/PINERF-1.)/0.28
      LEP=LEP1
      LEP=MAX0(1,LEP)
      IDUM=0
      DO 10 K=1,5
 10   IDUM=IDUM+LADJ(K)
      IF(IDUM.NE.0)  GO TO 316
C      LADJ IS USED IN ENTRANCE FROM EVAP, TO CONTROL THE TABLE THAT
C           IS BEING ADJUSTED
C---------------------------------------------------------------------
C                                        RESET
      SOR=0.
      IF(EXMAX1.LT.200.) GOTO 15
      SOR=SQRT(ACN*EXMAX1*0.1)
  15  DO 20  K=1,IDEN
      OMEGF(K)=0.
      OMEGG(K)=0.
      OMEGN(K)=0.
  20  OMEGP(K)=0.
      DO 30  K=1,IDNA
  30  OMEGA(K)=0.
      DO 40 K=1,IEDIM
      ENHANC(K)=0.
  40  ENHANS(K)=0.
C----------------------------------------------------------------------
C                  SET DIMENSIONS         EXMAX(N,P,A,G,F)
      EXMAXG=EXMAX1-DEL(1)
      EXMAXN=EXMAX1+Q(2)-DEL(2)
      EXMAXP=EXMAX1+Q(3)-DEL(3)
      EXMAXA=EXMAX1+Q(4)-DEL(4)
      EXMAXF=EXMAX1-DEL(5)
C----------------------------------------------------------------------
C                                           NW(G,N,P,A,F)
      IT=1
      EXMAX=EXMAXG
      KLAST=MIN0(N1+2,JDIM)
      KMAX=KLAST+1
    4 DO 302 K=1,KLAST
      KI=KMAX-K
      IF(EJAY(KI).LE.EXMAX) GO TO 303
  302 CONTINUE
  303 GO TO (5,6,7,8,9),IT
    5 IT=2
      EXMAX=EXMAXN
      KLAST=MIN0(N1+LN-1,JDIM)
      KMAX=KLAST+JDIM+1
      NWG=KI
      GO TO 4
    6 IT=3
      EXMAX=EXMAXP
      KLAST=MIN0(N1+LP-1,JDIM)
      KMAX=KLAST+2*JDIM+1
      NWN=KI-JDIM
      GO TO 4
    7 IT=4
      EXMAX=EXMAXA
      KLAST=MIN0(N1+LA-1,JDIM)
      KMAX=KLAST+3*JDIM+1
      NWP=KI-2*JDIM
      GO TO 4
    8 IT=5
      EXMAX=EXMAXF
      KLAST=MIN0(N1,JDIM)
      KMAX=KLAST+4*JDIM+1
      NWA=KI-3*JDIM
      GOTO 4
 9    NWF=KI-4*JDIM
C-----------------------------------------------------------------------
C                                              EXMIN(G,N,P,A,F)
C                                         DECIMAL PARTS OF EXMAX
      EXMING=EXMAXG-AINT(EXMAXG)
      EXMINN=EXMAXN-AINT(EXMAXN)
      EXMINP=EXMAXP-AINT(EXMAXP)
      EXMINA=EXMAXA-AINT(EXMAXA)
      EXMINF=EXMAXF-AINT(EXMAXF)
C-----------------------------------------------------------------------
C                                              MW(G,N,P,A,F)
      MWG=EXMAXG-EXMING+1.0001
      MWG=MIN0(MWG,IDEN/NWG)
      MWG=MAX0(MWG,1)
      MWN=EXMAXN-EXMINN+1.0001
      MWN=MIN0(MWN,IDEN/NWN)
      MWN=MAX0(MWN,1)
      MWP=EXMAXP-EXMINP+1.0001
      MWP=MIN0(MWP,IDEN/NWP)
      MWP=MAX0(MWP,1)
      MWA=EXMAXA-EXMINA+1.0001
      MWA=MIN0(MWA,IDNA/NWA)
      MWA=MAX0(MWA,1)
      MWF=EXMAXF-EXMINF+1.0001
      MWF=MIN0(MWF,IDEN/NWF)
      MWF=MAX0(MWF,1)
C-----------------------------------------------------------------------
C                                   AM DEPENDENCE AF/AN AND SHELLS
      AMASS=MTGT+MPROJ+2-NSTEP-NROW
      ZEE=IZT+IZP+1-NROW
      AN=AMASS-ZEE
      JSHAP=97
      JBAR=99
      DLIM=120/AMASS
      KJC=4*JDIM
      DO 3040 KJ=1,JDIM
      J=KJC+KJ
      AA=KJ-1
      CALL FISROT(AMASS,ZEE,AN,AA,DELR,SPT,EROT,0.D0)
      IF(KJ.EQ.1)SPT0=SPT
      IF(KJ.EQ.31)SPT30=SPT
      IF(SPT0.LT.5.) GOTO 3030
      IF(KJ.EQ.2) HHF=SPT-SPT0
      IF(KJ.EQ.2) HHG=DELR
      SPT1=AA*(AA+1)*HHF+SPT0 -SPT
      DELR1=AA*(AA+1)*HHG -DELR
      IF((SPT1.GT.DLIM .OR. DELR1.GT.DLIM) .AND. JSHAP.EQ.97)JSHAP=KJ-1
      IF(SPT-DELR .LT. DLIM .AND.JBAR.EQ.99) JBAR=KJ-1
 3030 DUM=(SPT-DELR)/SPT0
      DUM=DMAX1(AMPAR1,DUM)
 3040 ARRAY(KJ)=DUM
      PINERF=(SPT30-SPT0)/(30.*30.)
      RJLIM=0.5*(JBAR+JSHAP)
      RJLIM=DMAX1(AMPAR2,RJLIM)
      RJDIF=0.2*(JBAR-JSHAP)
      RJDIF=DMAX1(AMPAR3,RJDIF)
      RJDIF=DMIN1(AMPAR5,RJDIF)
      IF(LPRINT.LT.3) WRITE(6,3041) JSHAP,JBAR,RJLIM,RJDIF
 3041 FORMAT(' DENSTY: JSHAP=',I3,4X,'JBAR=',I3,4X,2F8.1)
C-----------------------------------------------------------------------
C                                   GAMMA PARAMETERS
      LE=1
      IT=1
      IT1=IT
  304 CR=R(1)
      CAL=AL(1)
      MA=KZ(NSTEP)+IN(NSTEP)
      IZ=KZ(NSTEP)
      IMAX=MWG
      KE0=0
      NJ=NWG
      DELTAC=DELTA(1)
      SHELL=SHELK(1)
      EXMAX=EXMAXG
      KNDEX=0
      KJC=0
      PINERT=(YRSMO(KJC+31)-YRSMO(KJC+1))/(31.*30.)
      PINERT=PINERT**1.5 /24.
      IAX=0
C     KNDEX IS OMEG OFFSET,KJC IS EJAY OFFSET
C***********************************************************************
C                              CALCULATE LEVEL DENSITY
  305 RAL=CR*CAL
      IF(EXMAX.LE.0.)GOTO 3111
      IF(CR.LE.0.)      RAL=1.
      SC=CAL
C-----------------------------------------------------------------------
C                                  ENHANCEMENT FUNCTION
C                 USE ONLY FOR EXCIT < 100 MEV AND LOW ANG.MOMENTA
C                 USES SAME FOR GAMMA AND N CHANNELS
      IF(IENH.EQ.0) GOTO 360
      IF(LE.NE.1) GOTO 360
      KEMAX=EXMAX+5.
      KEMAX=MIN0(100,KEMAX)
      DO 350 KE=1,KEMAX
      E=KE
      KKE=KE+KE0
      IF(IT1.NE.5) GOTO 345
      BETAMX=BETAF
      UEFF=E
      IF(ESHELL.NE.0.) UEFF=E+SHELL*(1.-EXP(-E/ESHELL))
      GOTO 348
 345  CALL BETA(IZ,MA,SHELL,DUM,E,ESHELL,BETAMX,FROT,UEFF,AMPAR1,AMPAR2)
 348  IF(UEFF.LE.0.) GOTO 350
      TEMP=SQRT(UEFF/CAL)
      SIG2=TEMP/PINERF
      IF(IT1.NE.5) SIG2=(1.+0.28*BETAMX)*TEMP/PINERT
      TCROT=37.8*BETAMX/FLOAT(MA)**0.33333
      EDROT=CAL*TCROT*TCROT
      IF(EDROT.LE.0.)GOTO 350
      DUM=UEFF/EDROT
      IF(EDCOLL.GT.0.) DUM=DUM+UEFF/EDCOLL
      ENHANC(KKE)=0.
      IF(DUM.GT.150.) GOTO 350
      ENHANC(KKE)=SIG2*EXP(-DUM)
      ENHANS(KKE)=UEFF
 350  CONTINUE
C     IF(IENHP.EQ.0 ) GOTO 360
C     WRITE(6,352) IT1,SHELL,BETA0,BETAF
 352  FORMAT(/' DENSTY  ENHANC-FACTOR  IT,SHELL,BETA0=',I3,F8.2,F8.3,
     1       2X,'BETAF=',F8.3/)
C     WRITE(6,354) (ENHANC(KE+KE0),KE=1,KEMAX)
 354  FORMAT(10F8.2)
C-----------------------------------------------------------------------
C                                        LEVDENS AS IN BM69
C--------------------------------------------------
 360  DO 490 KJ=1,NJ
      RJ=KJ-1
      IROW=(KJ-1)*IMAX
      J=KJC+KJ
      YJ=EJAY(J)
      IF(LDBM.EQ.3 .OR. LDBM.EQ.4) YJ=YRSMO(J)
      PREEX=1.
      IF(IAX.EQ.0) PREEX=KJ+KJ-1
      IF(IPREEX.EQ.1) PREEX= PREEX*PINERT
      SHELJ=YJ-YRSMO(J)+SHELL
 405  IF(IT.NE.5) GOTO 410
C      IT=5 IS FISSION CHANNEL
      SSS=SHELK(1)-SHELK(5)
      IF(ISHELL.EQ.0) SSS=0.
      IF(LDBM.EQ.3 .OR. LDBM.EQ.4) YJ=YRSMO(J)-SSS
      ALJF=(AL(5)-AL(1))*(1./(1.+EXP((RJ-RJLIM)/RJDIF) ))
      IF(LDBM.EQ.2 .OR.LDBM.EQ.4) ALJF=AL(5)-AL(1)
      CAL=AL(1)+ALJF
      SHELJ=SHELL*ARRAY(KJ)
 410  DUM=1.05*ABS(SHELJ)
      ESHELJ=DMAX1(ESHELL,DUM)
      KUP=1
      IF(IAX.EQ.1) KUP=KJ
      IU=1
C--------------------------------------------
      DO 470 KE=LE,IMAX
      INDEX=KNDEX+IROW+KE
      E=EXMAX+1.-FLOAT(KE)
      OM=0.
      PREEX1=PREEX
      DO 455 KVAL=1,KUP
      U=E-YJ
      IF(IAX.EQ.1) U=U-YJK(KVAL)
      IF(U.LE.0.) GOTO 455
      HE=1.
      IF(U.LT.UCRIT .AND. KHS.EQ.1) HE=(1.-U*U/(UCRIT*UCRIT))
      U=U -DELTAC*HE
      IF(U.LE.0.) GOTO 455
      IF(IENH.EQ.0) GOTO 450
      IU=U+0.5
      IU=MAX0(1,IU)
      IU=MIN0(100,IU)
      RU=0.
      IF(U.LT.100.) RU=U-IU
      U=ENHANS(IU+KE0)+RU
      IF(U.LE.0.) GOTO 455
 450  CALS=CAL
      U1=DMAX1(AMPAR4,U)
      FF=1.
      IF(IENH.EQ.0) FF=1.+SHELJ*(1.-EXP(-U1/ESHELJ))/U1
      IF(ISHELL.NE.0)CALS=CAL*FF
      SQ=SQRT(CALS*U)
      SQ=SQ+SQ-SOR
      OM= OM+       SQRT(CALS)*PREEX1*EXP(SQ)/(U1*U1)
      IF(KVAL.EQ.2) PREEX1=2.*PREEX
 455  CONTINUE
      OMEGG(INDEX) =OM * (1.+ENHANC(IU+KE0))
C     GOTO 470
C460  OMEGG(INDEX)=0.
 470  CONTINUE
C---------------------------------------------
 490  CONTINUE
C
      IF(LEVPRN.NE.1) GOTO 3111
      WRITE(6,493) MA,KZ(NSTEP),EXMAX,CAL,UCRIT,ESHELL,SHELL,
     1             DELTAC
 493  FORMAT(' DENSTY MASS,Z',2I4/' EXMAX AL UCRIT ESHELL',4F10.3/
     1         ' SHELL,DELTA',2F12.5,2F10.3/)
      DO 495 KJ=1,NJ
      IROW=(KJ-1)*IMAX
      WRITE(6,491) KJ,IT,LE,IMAX,NJ
 491  FORMAT(' DENSTY KJ,IT,LE,IMAX,NJ',5I6)
      DO 492 KE=LE,IMAX
      POPA(KE)=0.
      INDEX=KNDEX+IROW+KE
      IF(OMEGG(INDEX).LE.0.) GOTO 492
      POPA(KE)=DLOG (OMEGG(INDEX))
 492  CONTINUE
      WRITE(6,496) (POPA(KE),KE=LE,IMAX)
 495  CONTINUE
 496  FORMAT(10F10.3)
      DO 497 KE=LE,IMAX
 497  POPA(KE)=0.
C
C-----------------------------------------------------------------------
C                                                  SWITCH
 3111 GO TO (312,313,314,515,315,320,322,324,330),IT
C-----------------------------------------------------------------------
C                                     NEUTRON PARAMETERS
  312 IT=2
      IT1=IT
 3121 CR=R(2)
      CAL=AL(2)
      MA=KZ(NSTEP)+IN(NSTEP)-1
      IZ=KZ(NSTEP)
      DELTAC=DELTA(2)
      SHELL=SHELK(2)
      IMAX=MWN
      KE0=0
      NJ=NWN
      EXMAX=EXMAXN
      KNDEX=IDEN
      KJC=JDIM
      PINERT=(YRSMO(KJC+31)-YRSMO(KJC+1))/(31.*30.)
      PINERT=PINERT**1.5 /24.
      IAX=0
      GO TO 305
C-----------------------------------------------------------------------
C                                       PROTON PARAMETERS
  313 IT=3
      IT1=IT
 3131 CR=R(3)
      IU0=0
      CAL=AL(3)
      MA=KZ(NSTEP)+IN(NSTEP)-1
      IZ=KZ(NSTEP)-1
      DELTAC=DELTA(3)
      SHELL=SHELK(3)
      IMAX=MWP
      KE0=100
      NJ=NWP
      EXMAX=EXMAXP
      KNDEX=2*IDEN
      KJC=2*JDIM
      PINERT=(YRSMO(KJC+31)-YRSMO(KJC+1))/(31.*30.)
      PINERT=PINERT**1.5 /24.
      GO TO 305
C-----------------------------------------------------------------------
C                                      ALPHA  PARAMETERS
  314 IT=4
      IT1=IT
 3141 CR=R(4)
      CAL=AL(4)
      MA=KZ(NSTEP)+IN(NSTEP) -4
      IZ=KZ(NSTEP)-2
      DELTAC=DELTA(4)
      SHELL=SHELK(4)
      IMAX=MWA
      KE0=200
      NJ=NWA
      EXMAX=EXMAXA
      KNDEX=3*IDEN
      KJC=3*JDIM
      PINERT=(YRSMO(KJC+31)-YRSMO(KJC+1))/(31.*30.)
      PINERT=PINERT**1.5 /24.
      IAX=0
      GO TO 305
C-----------------------------------------------------------------------
C                                    FISSION  PARAMETERS
  515 IT=5
      IT1=IT
 5151 CR=R(5)
      CAL=AL(5)
      MA=KZ(NSTEP)+IN(NSTEP)
      IZ=KZ(NSTEP)
      DELTAC=DELTA(5)
      SHELL=SHELK(5)
      IMAX=MWF
      KE0=300
      NJ=NWF
      EXMAX=EXMAXF
      KNDEX=3*IDEN+IDNA
      KJC=4*JDIM
      PINERT=(YRSMO(KJC+31)-YRSMO(KJC+1))/(31.*30.)
      IF(AX.LE.0.) GOTO 520
C     APPROXIMATE DETERMINATION OF INERTIA PARAMETER PARALLEL TO K-AXIS
      IAX=1
      XLD=FLOAT(IZ*IZ)/FLOAT(MA)
      DUM=(MA-2*IZ)/FLOAT(MA)
      DUM=50.883*(1.-1.7826*DUM*DUM)
      XLD=XLD/DUM
      PINERK=TABIP(XLDM,THMIN,9,2,0,XLD,0)
      DUM=(YRSMO(31)-YRSMO(1))/(31.*30.)
      THMN=PINERK
      PINERK=DUM/PINERK
      PINER0=PINERK-PINERT
      IF(I7.NE.99)
     1WRITE(6,517) PINERT,PINERK,DUM,XLD,THMN
 517  FORMAT(' DENSTY  INERF, INERK, INERGS',3F9.4,4X,'XLD THMIN',2F9.4)
      I7=99
      DO 518 KVAL=1,NJ
      DUM=KVAL-1
 518  YJK(KVAL)= PINER0*(DUM+1.)*DUM
C
      PINERT=PINERT*SQRT(PINERK)/24.
      GOTO 305
 520  PINERT=PINERT**1.5/24.
      GO TO 305
C***********************************************************************
C                                            EXIT
  315 CONTINUE
      LEVPRN=0
      RETURN
C***********************************************************************
C                                     ENTRY POINT FROM EVA
C-----------------------------------------------------------------------
C                                         READJUST GAMMAS
  316 IF(LADJ(1).NE.1)     GO TO 320
      IT=6
      IT1=1
      LE=MWG-LEP
      NJ=NWG
      IMAX=MWG
      KNDEX=0
C
  317 DO 318 KJ=1,NJ
      INDEX=KNDEX+(KJ-1)*IMAX
      DO 318 K=1,LE
      INDEX=INDEX+1
      ITM=INDEX+LEP
      OMEGG(INDEX)=OMEGG(ITM)
  318 OMEGG(ITM)=0.
      GO TO (319,321,323,325, 335),IT1
C--------------------------------------------
  319 LE=MWG+1-LEP
      LADJ(1)=0
      EXMAXG=EXMAXG-FLOAT(LEP)
      GO TO 304
C-----------------------------------------------------------------------
C                                          READJUST NEUTRONS
  320 IF(LADJ(2).NE.1)      GO TO 322
      IT=7
      IT1=2
      LE=MWN-LEP
      NJ=NWN
      IMAX=MWN
      KNDEX=IDEN
      GO TO 317
C-------------------------------------------
  321 LE=MWN+1-LEP
      LADJ(2)=0
      EXMAXN=EXMAXN-FLOAT(LEP)
      GO TO 3121
C-----------------------------------------------------------------------
C                                         READJUST PROTONS
  322 IF(LADJ(3).NE.1)     GO TO 324
      IT=8
      IT1=3
      LE=MWP-LEP
      NJ=NWP
      IMAX=MWP
      KNDEX=2*IDEN
      GO TO 317
C--------------------------------------------
  323 LE=MWP+1-LEP
      LADJ(3)=0
      EXMAXP=EXMAXP-FLOAT(LEP)
      GO TO 3131
C-----------------------------------------------------------------------
C                                         READJUST ALPHAS
  324 IF(LADJ(4).NE.1) GOTO 330
      IT=9
      IT1=4
      LE=MWA-LEP
      NJ=NWA
      IMAX=MWA
      KNDEX=3*IDEN
      GO TO 317
C------------------------------------------
  325 LE=MWA+1-LEP
      LADJ(4)=0
      EXMAXA=EXMAXA-FLOAT(LEP)
      GO TO 3141
C----------------------------------------------------------------------
C                                         READJUST FISSION
  330 IF(LADJ(5).NE.1) GOTO 315
      IT1=5
      LE=MWF-LEP
      NJ=NWF
      IMAX=MWF
      KNDEX=3*IDEN+IDNA
      GO TO 317
C------------------------------------------
  335 LE=MWF+1-LEP
      LADJ(5)=0
      EXMAXF=EXMAXF-FLOAT(LEP)
      GO TO 515
      END
