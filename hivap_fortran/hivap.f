      PROGRAM HIVAP
C
C  Program 'HIVAP' for oalculations of fusion- and evaporation residue
C  cross sections
C  Author: W. Reisdorf, GSI
C  modified for use at 'rzri6f' by F.P.Hessberger
C  Informations on modifications on
C     rzri6f:/u/hess/fortpro/hivapn/aainfo.dat (in german)
C  Option for fission hindrance include; working up to about 130 MeV
C  ecitation energy, for higher values eventually problems with
C  level densities ( >3x10**39)
C  Status: June, 7th, 1994
C
      INCLUDE 'common.f'
      COMMON/DELFIS/EDELFIS,ADELFIS,INOF
      COMMON/PUSHPA/SPUSH,SCREIPA
      DIMENSION ARRAY(40),SHEL(5)
      CHARACTER*1 AW(72)
      CHARACTER*4 BLNK
      CHARACTER*10 BLNK2
      CHARACTER*1 SZ
      CHARACTER*6 AUT1
      CHARACTER*6 AUT2
      CHARACTER*35 AUTOR
      CHARACTER*12 BEARB
      DATA BLNK/'    '/ BLNK2/'          '/
      DATA AUT1/'F.P.Hessberger'/,BEARB/'modified by '/
C
      SZ=CHAR(223)
      AUTOR=BLNK2//BEARB//AUT1
C-----------------------------------------------------------------------
C
      SCREIPA=0.
      AMPAR1=0.
      INOF=0
      LEVPRN=0
      CALL ALLOCH
      EPS=0.01
      IPOPS=4096
C     DIMENSION OF POP,POPN,POPP,POPA
      MASSES=9
C     LOGICAL UNIT FOR MASS EXCESS TABLE
      IYR=14
C     LOGICAL UNIT FOR TRANSMISSION COEFF IS IUNIT=21
C     LOGICAL UNIT FOR INPUT SUMMARY IS 23
      IDEN=2500
      IDNA=5000
      NAUP=25
      NZUP=25
      NAUP2=NAUP+2
      NZUP2=NZUP+2
      NUMIS=5
      JDIM=300
C     ALSO DEFINED IN 'TOT'
      IEDIM=400
C     ALSO DEFINED IN 'TOT'  (SPECE)
      KEN=25
      KEP=25
      KEP1=20+KEP
      KEG=20
      KEA=35
      KEA1=30+KEA
      KEF=20
      LN=12
      LP=12
      LA=20
      CUT=0.
      FRACT2=0.
      ABSMIN=0.
      PRCN=0.
      SIGLOW=0.
      VRN=0.
      R0RN=0.
      ADIFRN=0.
      VIN=0.
      R0IN=0.
      ADIFIN=0.
      VRP=0.
      R0RP=0.
      ADIFRP=0.
      VIP=0.
      R0IP=0.
      ADIFIP=0.
      RCLMBP=0.
      CBFACP=1.
      VRA=0.
      R0RA=0.
      ADIFRA=0.
      VIA=0.
      R0IA=0.
      ADIFIA=0.
      RCLMBA=0.
      CBFACA=1.
      SIGR0=0.
      NOLEP=1
      LEP=1
      NOLJI=1
      NOLJF=1
      IPTR=0
      WRITE(6,10)
 10   FORMAT(10X,'Program HIVAP,  Author  W.Reisdorf  ',
     1' Version for "rzri6f",  7. 6. 1994'/)
      WRITE(6,*) AUTOR
      WRITE(6,*) BLNK2
C
C******************************************************************
C                                                 NEW CASE
 99   CONTINUE
 100  FORMAT(18A4)
 101  FORMAT(18I4)
C 102  FORMAT(9F8.3)
      DO 105 K=1,IPOPS
      POP (K)=0.
      POPN(K)=0.
      POPP(K)=0.
      POPA(K)=0.
      PENTRY(K)=0.
 105  CONTINUE
      DO 107 IE=1,NENMX
      DO 107 NSTEP=1,NAUP
      DO 107 NROW=1,NZUP
 107  XSECTE(NROW,NSTEP,IE)=0.
      READ(5, 100,END=9999)(TITLE(I),I=1,18)
      WRITE(6,100)(TITLE(I),I=1,18)
      WRITE( 23  ,100)(TITLE(I),I=1,18)
      WRITE(6,*) BLNK2
      NOMORE=0
      DO 110 I=1,18
      IF(TITLE(I).NE.BLNK) GOTO 120
 110  CONTINUE
      GOTO 9999
C---------------------------------------------------------------
 120  CALL MYRD(ARRAY,18,23,5)
      MPROJ = ARRAY(1 )+EPS
      IZP   = ARRAY(2 )+EPS
      MTGT  = ARRAY(3 )+EPS
      IZT   = ARRAY(4 )+EPS
      ISHELL= ARRAY(5 )+EPS
      ISHELF= ARRAY(6 )+EPS
      IPAIR = ARRAY(7 )+EPS
      IF(ISHELL.EQ.0) IPAIR=0
      MC    = ARRAY(8 )+EPS
      IF(ISHELL.EQ.0) MC=1
      MP    = ARRAY(9 )+EPS
      IF(ISHELL.EQ.0) MP=1
      IBF   = ARRAY(10)+EPS
      IFISRT= ARRAY(11)+EPS
      NOF   = ARRAY(12)+EPS
      NON   = ARRAY(13)+EPS
      NOP   = ARRAY(14)+EPS
      NOA   = ARRAY(15)+EPS
      NOG   = ARRAY(16)+EPS
      IDISC = ARRAY(17)+EPS
      IGAM  = ARRAY(18)+EPS
!---The following added by N.Wang ----------------
	READ(101,*)MPROJ,IZP,MTGT,IZT 
!	PRINT*,"INPUT  A1, Z1, A2, Z2:"  !BY K.ZHAO
!	READ(*,*)MPROJ,IZP,MTGT,IZT
C--------------------------------------------------
C  Includef for fission hindrance  NOF=2, 2.5.94 FPH
C-----------------------------
      IF(NOF.NE.2) GOTO 1201
      INOF=1
      NOF=0
      CALL MYRD(ARRAY,2,23,5)
      EDELFIS=ARRAY(1)+EPS
      ADELFIS=ARRAY(2)+EPS
C-----------------------------
 1201 IF(IDISC.NE.0) WRITE(25, 9993)(TITLE(I),I=1,18)
      NZ1=0
      ACN=MTGT+MPROJ
      AP=MPROJ
      ZP=IZP
      AT=MTGT
      ZT=IZT
      IACN=ACN
      ZCN=IZT+IZP
      IZCN=ZCN
      NOPUT=0
      IF(NZ1.LT.0) NOPUT=1
      NZ1=IABS(NZ1)
      NZ1=MAX0(1,NZ1)
C-----------------------------
      CALL MYRD(ARRAY,18,23,5)
C-----------------------------
      NA    = ARRAY(1 ) + EPS
      NZ    = ARRAY(2 ) + EPS
      MASSES= ARRAY(3 ) + EPS
      NUMB  = ARRAY(4 ) + EPS
      IOVER = ABS(ARRAY(5))+EPS
      IF(ARRAY(5).LT.0.) IOVER=-IOVER
      INERT = ARRAY(6 ) + EPS
      INERF = ARRAY(7 ) + EPS
      FINERT= ARRAY(8 )
      IF(FINERT.LE.0.)FINERT=1.
      ILIM  = ARRAY(9 ) + EPS
      LPRINT= ARRAY(10) + EPS
      LOGUN = ARRAY(11) + EPS
      ICOR  = ARRAY(12) + EPS
      NUMISO= ARRAY(13) + EPS
      LDBM  = ABS(ARRAY(14)) + EPS
      IF(ARRAY(14).LT.0.) LDBM=-LDBM
      LDBM5=0
      IF(LDBM.EQ.5) LDBM5=1
      IRAST = ARRAY(15) + EPS
      IOWKB = ARRAY(16) + EPS
      ITRANS= ARRAY(17) + EPS
      JFJI  = ARRAY(18) + EPS
      IF(NA.EQ.0) NA=14
      NZ=MAX0(1,NZ)
      NZ=MIN0(14,NZ)
      LOGUN=0
      IF(MASSES.NE.0)MASSES=9
      DEL(5)=0.
      DO 130 K=1,5
 130  DJ(K)=0.
      IF(JFJI.NE.1) GOTO 132
C------------------------------------
      CALL MYRD(ARRAY,3,23,5)
C------------------------------------
      DJ(2)=ARRAY(1)
      DJ(3)=ARRAY(2)
      DJ(4)=ARRAY(3)
C      IF(JFJI.EQ.1) WRITE(6,102) DJ(2),DJ(3),DJ(4)
C-----------------------------------------------------------------------
C                                                  GAMMA NORM
 132  IOUT=LPRINT
      CALL GAMMAS(ACN,IOUT,DELG,IGAM,NOE2,STRIPE,IOPT)
      DEL(1)=DELG
C 1320 WRITE(6,1321) DEL(1)
C 1321 FORMAT(' DEL(1) AFTER GAMMAS',F8.1)
C-----------------------------------------------------------------------
C                                       LEVEL DENSITY PARAMETERS
      CALL MYRD(ARRAY,8,23,5)
C------------------------------
      CST   =ARRAY(1 )
      CLD   =ARRAY(2 )
      BARFAC=ARRAY(3 )
      ESHELL=ARRAY(4 )
      BAR0  =ARRAY(5 )
      SHELL0=ARRAY(6 )
      DELT  =ARRAY(7 )
      QVALUE=ARRAY(8 )
C------------------------------
      IPREEX=1
      AX=0.
      BETA0=0.
      UCRIT=0.
      EDCOLL=0.
      IENH=0
      IF(LDBM5.EQ.1) THEN
        CALL MYRD(ARRAY,7,23,5)
        LDBM=ABS(ARRAY(1))+EPS
        IF(ARRAY(1).LT.0.) LDBM=-LDBM
        IPREEX=ARRAY(2)+EPS
        AX=ARRAY(3)
        IENH=ARRAY(4)+EPS
        BETA0=ARRAY(5)
        EDCOLL=ARRAY(6)
        UCRIT=ARRAY(7)
      ENDIF
 133  IF(BARFAC.LE.0.) BARFAC=1.
      IF(ESHELL.LE.0.) ESHELL=18.5
      IF(DELT.LE.0.)DELT=11.
 135  FORMAT(4F8.3,6I4)
      IF(UCRIT.LE.0.)UCRIT=5.
      KHS=0
C     KHS=1 IS PAIRING ACCORDING TO KHS, NOT YET TESTED
      IF(IENH.NE.0) AX=1.
      IF(TZERO.EQ.0.)TZERO=0.2
      IF(LDBM.NE.(-1)) GOTO 220
      IRC=NA*NZ
      CALL MYRD(ALS,IRC,23,5)
      IRC=NA*NZ
      CALL MYRD(DELTAS,IRC,23,5)
C-------------------------------------------------------------
C                                         CUT-OFF PARAMETERS
 220  IF(ILIM.NE.1)GOTO 226
      CALL MYRD(ARRAY,9,23,5)
      CUT   =ARRAY(1 )
      FRACT2=ARRAY(2 )
      ABSMIN=ARRAY(3 )
      PRCN  =ARRAY(4 )
      SIGLOW=ARRAY(5 )
C     DEL(1)=ARRAY(6 )
      DUM   =ARRAY(6 )
      NOLEP =ARRAY(7 )+EPS
      NOLJI =ARRAY(8 )+EPS
      NOLJF =ARRAY(9 )+EPS
      IF(JFJI.EQ.1)NOLEP=1
      IF(JFJI.EQ.1)NOLJI=1
      IF(JFJI.EQ.1)NOLJF=1
C      IF(ILIM.EQ.1)WRITE(6,225) CUT,FRACT2,ABSMIN,PRCN ,SIGLOW,DEL(1),
C     1                          NOLEP,NOLJI,NOLJF
C 225  FORMAT(' CUT,FRACT2,ABSMIN=',3E12.4/' PRCN,SIGLOW,DEL1',3F8.3,4X,
C     1       'NOLEP,NOLJI,NOLJF=',3I4)
 226  IF(CUT.LE.0. .OR. CUT.GT.0.15)CUT=0.1E-2
      PRCNT=1.-PRCN*0.01
      IF(PRCN.LE.0. .OR. PRCN.GT.20.) PRCNT=1.-CUT
      IF(FRACT2.LE.0.)FRACT2=0.001
      IF(DEL(1).LE.0.)DEL(1)=1.
C-------------------------------------------------------------
C                                       TRANSMISSION
      TLOW=0.
      DEL(2)=0.1
C      WRITE(6,2261) DEL(2)
C 2261 FORMAT(' DEL(2)=',F8.2)
      QT=0.
      QQ2=0.
      IOMN=0
      IOMP=0
      IOMA=0
      IPTR=0
      IF(IOVER.GE.0) GOTO 229
       WRITE(6,2260)
 2260 FORMAT(' IOVER nicht Standart OM')
      CALL MYRD(ARRAY,11,23,5)
      IOVER =ARRAY(1 ) + EPS
      IOWKB =ARRAY(2 ) + EPS
      ITRANS=ARRAY(3 ) + EPS
      SIGLOW=ARRAY(4 )
      TLOW  =ARRAY(5 )
      DEL(2)=ARRAY(6)
      QT    =ARRAY(7 )
      QQ2=QT
      IOMN=ARRAY(8)+EPS
      IOMP=ARRAY(9)+EPS
      IOMA=ARRAY(10)+EPS
      IPTR  =ARRAY(11) + EPS
C------------------------------  OPT MOD N
      IF(IOMN.EQ.0) GOTO 2280
      CALL MYRD(ARRAY,6,23,5)
      VRN   = ARRAY(1)
      R0RN  = ARRAY(2)
      ADIFRN= ARRAY(3)
      VIN   = ARRAY(4)
      R0IN  = ARRAY(5)
      ADIFIN= ARRAY(6)
 2280 IF(IOMP.EQ.0) GOTO 2281
      CALL MYRD(ARRAY,8,23,5)
      VRP   = ARRAY(1)
      R0RP  = ARRAY(2)
      ADIFRP= ARRAY(3)
      VIP   = ARRAY(4)
      R0IP  = ARRAY(5)
      ADIFIP= ARRAY(6)
      RCLMBP= ARRAY(7)
      CBFACP= ARRAY(8)
 2281 IF(IOMA.EQ.0) GOTO 229
      CALL MYRD(ARRAY,8,23,5)
      VRA   = ARRAY(1)
      R0RA  = ARRAY(2)
      ADIFRA= ARRAY(3)
      VIA   = ARRAY(4)
      R0IA  = ARRAY(5)
      ADIFIA= ARRAY(6)
      RCLMBA= ARRAY(7)
      CBFACA= ARRAY(8)
 229  IF(SIGLOW.EQ.0.)SIGLOW=0.001
C-------------------------------------------------------------
C                                       ISOMERS
      IF(NUMISO.EQ.0) GOTO  230
      NUMISO=MIN0(NUMIS,NUMISO)
      IRC=NUMISO*3
      ISO=0
      CALL MYRD(ARRAY,IRC,23,5)
      DO 2290 I=1,IRC,3
      ISO=ISO+1
      ISOZ(ISO)=ARRAY(I  )+EPS
      ISOA(ISO)=ARRAY(I+1)+EPS
 2290 ISOJ(ISO)=ARRAY(I+2)+EPS
      IF(LPRINT.EQ.0) WRITE(6,2291) NUMISO,(ISOZ(I),ISOA(I),ISOJ(I),
     1                I=1,NUMISO)
 2291 FORMAT(' NUMISO',I3,4X,12I4)
C-------------------------------------------------------------
C                                        INERTIAS
 230  IF(INERT.EQ.0) GOTO 255
      CALL MYRD(RATIOS,6,23,5)
      IF(LPRINT.EQ.0) WRITE(6,2301) RATIOS
 2301 FORMAT(' RATIOS',6F8.3)
C
C-----------------------------------------------------------------------
C                              BARRIERS AND    SEPARATION ENERGIES
  255 CALL SEPENF(MTGT,IZT,MPROJ,IZP,QC,NA,NZ,POP,3,ISHELL,SHELLS,
     1            SHELL0,DELTAS,IPAIR,KHS)
      IF(ABS(QVALUE).LT.0.1) GOTO 270
      QC=QVALUE
 270  CONTINUE
      QEFF=QC
      IF(ISHELL.EQ.0)QEFF=QC+SHELLS(1,1)
      IF(ISHELL.EQ.0)SHELLS(1,1)=0.
      WRITE(6,269) QC,QEFF
 269  FORMAT(/' benutzter Q-Wert',F10.3,1X,'MeV',4X,
     1'effektiver Q-Wert',F10.3,1X,'MeV')
C
C***********************************************************************
C                                                  ELAB LOOP
 400  CONTINUE
      IRC=-10
      IEOF=0
      CALL MYRD(ARRAY,IRC,23,5)
      IF(IRC.EQ.99) IEOF=1
      IF(IEOF.EQ.1) GOTO 3000
      IF(ARRAY(1).LE.0.) GOTO 3000
      EQ    =ARRAY(1)
      IEXC  =ARRAY(2 ) +EPS
      IFUS  =ARRAY(3 ) +EPS
      LIMBAR=ARRAY(4 ) +EPS
      JLOWER=ARRAY(5 ) +EPS
      JUPPER=ARRAY(6 ) +EPS
      NEWFIS=ARRAY(7 ) +EPS
      ITSTRT=ARRAY(8 ) +EPS
      JFIS  =ARRAY(9 ) +EPS
      EQ1   =ARRAY(10)
!--------- added by N.Wang 2007-6-15-------------------
	READ(101,*)ECM  
	IF(ECM.LE.0)GOTO 3000
	EQ=ECM+QC
	write(*,'(1x,"Ecm=",f12.4)')Ecm
!------------------------------------------------------
      IFUS1=0
      IF(IFUS.EQ.10) IFUS1=1
      IF(IFUS.EQ.10) IFUS=8
      IF(IFUS.EQ.11) IFUS1=2
      IF(IFUS.EQ.11) IFUS=8
      IF(LPRINT.LT.3)WRITE(6,403) EQ,EQ1
 403  FORMAT(1H1,4X,' ENERGY',2F10.3)
      IF(EQ1.GT.0.) IFUS=8
      IF(EQ1.LT.EQ) GOTO 404
      H=EQ1
      EQ1=EQ
      EQ =H
 404  JLOWER=MAX0(JLOWER,1)
      JUPPER=MIN0(JUPPER,JDIM)
      IF(JFIS.GT.0 .AND. NOMORE.EQ.0) WRITE(19,100)(TITLE(I),I=1,18)
      NOMORE=1
      IF(NEWFIS.NE.1) GOTO 407
      CALL MYRD(ARRAY,2,23,5)
      CLD1=ARRAY(1)
      BFC=ARRAY(2)
      IF(BFC .GT.0.) BARFAC=BFC
      IF(CLD1.GT.0.)CLD=CLD1
 407  DO 408 I=1,IPOPS
      POPN(I)=0.
      POPP(I)=0.
      POPA(I)=0.
 408  POP(I)=0.
      DO 410 ISO=1,NUMISO
 410  SIGISO(ISO)=0.
      DO 411 I=1,JDIM
 411  FLAND(I)=0.
      DO 412 NROW=1,NZUP
      DO 412 NSTEP=1,NAUP
 412  XSECT(NROW,NSTEP)=0.
      DO 413 J=1,JDIM
      PER(J)=1.
      PFIS(J)=0.
 413  YFISS(J)=0.
      DO 4131 K=1,KEN
 4131 SPECNT(K)=0.
      DO 4132 K=1,KEP1
 4132 SPECPT(K)=0.
      DO 4133 K=1,KEA1
 4133 SPECAT(K)=0.
      DO 4134 K=1,KEF
 4134 SPECFT(K)=0.
      DO 4135 K=1,KEG
 4135 SPECGT(K)=0.
      DO 4136 K=1,KEG
 4136 SPECQT(K)=0.
C--------------------------------------------------------------------
C             GET COMPOUND NUCLEUS POPULATION
C
      AP=MPROJ
      AT=MTGT
      ZT=IZT
      ZP=IZP
      CM=AT/(AP+AT)
      IF(IEXC.EQ.1) GOTO 414
      ELAB=EQ
      ELAB1=EQ1
      EXCIT =ELAB*CM+QC
      EXCIT1=ELAB1*CM+QC
      IF(EQ1.LE.0.)EXCIT1=0.
      GOTO 415
 414  EXCIT=EQ
      EXCIT1=EQ1
      ELAB=(EXCIT-QC)/CM
      ELAB1=(EXCIT1-QC)/CM
      IF(EQ1.LE.0.)ELAB1=0.
C415  EXMAX2=EXCIT
 415  H4=0.04783*AP*CM*CM*ELAB
      C6=31.42/H4
C--------------------------------------
      IF(IFUS.EQ.9 ) GOTO 420
      IF(IFUS.EQ.1 .OR. IFUS.EQ.5) GOTO 417
      IF(IFUS.EQ.8) GOTO 430
      IF(AP.GT.4.) GOTO 416
      KPART=0
      IDUM=AP+0.01
      IDUM1=ZP+0.01
      IF(IDUM.EQ.4 .AND. IDUM1.EQ.2)KPART=3
      IF(IDUM.EQ.1 .AND. IDUM1.EQ.0)KPART=1
      IF(IDUM.EQ.1 .AND. IDUM1.EQ.1)KPART=2
      IF(KPART.EQ.0) GOTO 416
      LMAX=15
      IF(KPART.EQ.3) LMAX=31
      DO 418 L=1,50
 418  POPN(L)=0.
      IPOTS=0
      IF(IFUS.EQ.7) IPOTS=1
      IOMPR=1
      IF(LPRINT.GT.3) IOMPR=0
      CALL OM(KPART,ELAB,LMAX,AT,ZT,SIGABS,POPN(1),IPOTS,IOMPR)
      FLGRAZ=SIGABS/C6
      FLGRAZ=SQRT(FLGRAZ)-0.5
      FLGRAZ=DMAX1(AMPAR1,FLGRAZ)
      GOTO 419
 416  CALL PARAP(ELAB,AP,AT,ZP,ZT,FINERT,RATIOS,IFUS,LIMBAR,0,1,IDUM,
     1           POPN(1),POPN(401),FLGRAZ)
 419  IF(LIMBAR.EQ.1) GOTO 424
 417  CALL FUSE (ELAB,AP,AT,ZP,ZT,POPN(1),JDIM,IFUS,FLGRAZ)
      GOTO 424
C-------------------------------------------------- IFUS=9
 420  IRC=JUPPER-JLOWER+1
      CALL MYRD(POPN(JLOWER),IRC,23,5)
 421  FORMAT(9F8.3)
C--------------------------------------------------
 424  JDIM1=JDIM+1
      DO 425 J=1,JDIM
      JJ=JDIM1 -J
      IF(POPN(JJ).GT.0.) GOTO 426
 425  CONTINUE
 426  JMAX=JJ
      IF(JUPPER.EQ.0)JUPPER=JMAX
      SIGFUS=0.
      DO 427 J=1,JDIM
      IF(J.LT.JLOWER)POPN(J)=0.
      IF(J.GT.JUPPER)POPN(J)=0.
      SIGFUS=SIGFUS+POPN(J)
 427  CONTINUE
      DO 428 J=JDIM1,IPOPS
 428  POPN(J)=0.
      M2=1
      GOTO 460
C-------------------------------------- IFUS=8
 430  IF(IFUS1.EQ.0) CALL MYRD(ARRAY,9,23,5)
      IF(IFUS1.EQ.1) CALL MYRD(ARRAY,10,23,5)
      IF(IFUS1.EQ.2) CALL MYRD(ARRAY,14,23,5)
      V0=ARRAY(1)
      R0=ARRAY(2)
      D=ARRAY(3)
      Q2   =ARRAY(4)
      IF(Q2.LT.1.) Q2=1.09*ZT*(AT**0.666667)*Q2*(1.+0.15577*Q2)
C     IF Q2 IS LESS THAN ONE IT IS INTERPRETED AS DEFO BETA
      CRED=ARRAY(5)
      NOCURV=ARRAY(6)+0.01
      NOPROX=ARRAY(7)+0.01
      IOPT1=ABS(ARRAY(8))+0.01
      IF(ARRAY(8).LT.0.) IOPT1=-IOPT1
      IF(IOPT1.EQ.5) D=1.
      ITEST1=ABS(ARRAY(9))+0.01
      IF(ARRAY(9).LT.0.) ITEST1=-ITEST1
      IF(IFUS1.GT.0) SIGR0=ARRAY(10)
      CUTOF=2.5
      XTH=0.
      IF(IFUS1.LT.2) GOTO 436
      CUTOF=ARRAY(11)
      IF(CUTOF.LE.0.) CUTOF=2.5
      DUM=CUTOF*SIGR0
      IF(DUM.LT.50.) GOTO 434
      WRITE(6,432) SIGR0,CUTOF
 432  FORMAT(' UNREASONABLE PARAMETERS SIGR0 AND/OR CUTOF',2F8.3)
      GOTO 9991
 434  XTH=ARRAY(12)
      APUSH=ARRAY(13)
      IF(APUSH.LE.0.) APUSH=12.
      FPUSH=ARRAY(14)
      SPUSH=ARRAY(15)
      IF(FPUSH.LE.0.) FPUSH=0.75
 436  IF(CRED.LE.0.)CRED=1.
      IF(LIMBAR.EQ.1) LMAX=-1
C---------------------------------
      ELB=ELAB
      IENG=1
      M2=1
      RNORM=1.
      SIGFUS=0.
      FLGRAZ=0.
      IESTEP=1
      JMAX=0
      IF(EQ1.LE.0.) GOTO 444
      DUM=EXCIT-EXCIT1
      M2=DUM+1.5
 442  IDUM=M2/IESTEP
      IF(IDUM.LT.20) GOTO 443
      IESTEP=IESTEP+1
      GOTO 442
 443  DELB=IESTEP/CM
      RNORM=M2
C
 444  DO 452 IE=1,M2,IESTEP
      LMAX=0
      IF(LIMBAR.EQ.1)LMAX=-1
      FLGRZ=0.
      R01=R0
      IF(IOPT1.EQ.5) GOTO 4446
      IF(R01.GT.0.) GOTO 4448
 4446 AP3=AP**0.333333
      AT3=AT**0.333333
      RPS=1.28*AP3-0.76+0.8/AP3
      RTS=1.28*AT3-0.76+0.8/AT3
      IF(IOPT1.EQ.1 .OR. IOPT1.EQ.2 .OR.IOPT1.GT.3) GOTO 4445
      RPS=RPS-1./RPS
      RTS=RTS-1./RTS
 4445 R01=(RTS+RPS)/(AP3+AT3)
 4448 MPT=20
      CALL FUSIO (ELB ,AP,ZP,AT,ZT,Q2,V0,R01,D,XTH,APUSH,FPUSH,
     1            SIGR0,CRED,1,NOPROX,NOCURV,
     2        1  ,POP ,  FLGRZ,LMAX,SIGFU ,ITEST1,IOPT1,MPT,CUTOF)
      LMAX=MIN0(LMAX,JDIM)
      JMAX=MAX0(JMAX,LMAX)
      JMAX1=LMAX
      JMAX1=MIN0(JMAX1,JDIM)
      ICHECK=JMAX1*M2
      IF(ICHECK.LE.IPOPS)GOTO 445
      WRITE(6,446) M2,JMAX1,IPOPS
 446  FORMAT(' DIMENSIONS TOO LARGE  M2*JMAX1>IPOPS',3I6)
      GOTO 9991
 445  IF(JUPPER.GT.0)JMAX1=MIN0(LMAX,JUPPER)
      FLGRAZ=DMAX1(FLGRAZ,FLGRZ)
      DO 448 J=JLOWER,JMAX1
      DO 447 IEE=1,IESTEP
      INDEX=IE-1+IEE
      IF(INDEX.GT.M2) GOTO 447
      INDEX=INDEX+(J-1)*M2
      POPN(INDEX)=POP(J)/RNORM
      SIGFUS=SIGFUS+POPN(INDEX)
 447  CONTINUE
 448  CONTINUE
      ELB=ELB-DELB
 452  CONTINUE
C
      IF(EQ1.LE.0.) GOTO 450
      WRITE(6, 449) EXCIT,EXCIT1,M2,IESTEP,RNORM,DELB,SIGFUS
 449  FORMAT(' CROSS SECTION AVERAGED',2F6.1,2I6,3F8.3/)
 450  IF(JUPPER.EQ.0)JUPPER=JMAX
      DO 453 J=1,IPOPS
 453  POP(J)=0.
C---------------------------------
 460  IF(SIGFUS.GT.0.1E-15)GOTO 463
      WRITE(6,461) ELAB,EXCIT
 461  FORMAT(//'E_lab =',F8.3,' MeV',4X,'E* =',F8.3,' MeV',3X,
     1'nicht gerechnet, SIGFUS kleiner als 1 nb'//)
      GO TO 400
C---------------------------------
 463  N2=JUPPER
      JUPYR=JUPPER+6
      JUPYR=MIN0(JUPYR,JDIM-1)
      SIGMXZ =0.
      SIGEVA=0.
      DO 464 NROW=1,NZUP
 464  SIGZ(NROW)=0.
      SIGFIS=0.
 462  IU1=1-MOD(NZ1,2)
      IU2=0
      IDUM=(NZ1+1)/2
      IF(MOD(IDUM,2).NE.0) IU2=1
      DO 465 I=1,NAUP2
      DO 465 K=1,NZUP2
 465  STORE(I,K)=0.
      REWIND 10
      REWIND 11
      REWIND 12
      REWIND 13
      IF(NZ1.EQ.1) GOTO 469
 1001 FORMAT(18I4)
 1002 FORMAT(3F8.2,4E12.4)
 1003 FORMAT(9F6.1)
 1004 FORMAT(6E12.4)
      NZ1=NROW+1
      IF(MOD(NZ1,2).EQ.0) GOTO 469
      IUNIT=13-IU2
      IREAD=1
      DO 467 NSTEP=1,NA
      IF(STORE(NZ1+1,NSTEP).LT.2) GOTO 467
      CALL GETPUT(POP,IPOPS,IUNIT,IREAD,IDUM,IDUM,DUM,IDUM,IDUM,DUM,DUM)
 467  CONTINUE
 469  CONTINUE
C     WRITE(6,860) IPAIR,CST,DELT,DELTA
C
 468  IDUM=EXCIT+1.
      IF(IDUM.LT.IEDIM) GOTO 474
      WRITE(6,472) IEDIM
      GOTO 400
 472  FORMAT(' die hoechste erlaubte Anregungsenergie ist ',I5,' MEV')
 474  CONTINUE
C
C******************************* ROW LOOP (NZ) *************************
      EXCIT0=EXCIT
      IF(ISHELL.EQ.0) EXCIT=EXCIT-QC+QEFF
      EXMAX2=EXCIT
      DO 1000 NROW=NZ1,NZ
C
      IU1=1-IU1
      IF(MOD(NROW,2).EQ.0)IU2=1-IU2
      KZ(1)=IZCN+1-NROW
      ZPAR=KZ(1)
      APAR =IACN+1-NROW
      IN(1)=APAR-ZPAR
      DO 476 NSTEP=2,NA
      KZ(NSTEP)=KZ(1)
 476  IN(NSTEP)=IN(1)+1-NSTEP
C
 480  DO 490 NSTEP=1,NAUP
      TOTAL(NSTEP)=0.
      SUMEG(NSTEP)=0.
      SUMJG(NSTEP)=0.
      AVGJ(NSTEP)=0.
      AVGE(NSTEP)=0.
      AVGEF(NSTEP)=0.
      AVGJF(NSTEP)=0.
      SIG  (NSTEP)=0.
      SIGNN(NSTEP)=0.
      SIGP(NSTEP)=0.
      SIGA(NSTEP)=0.
      TRIM(NSTEP)=0.
      RP (NSTEP)=0.
      RLOSTG(NSTEP)=0.
      FLOST (NSTEP)=0.
      RLOSTN(NSTEP)=0.
      RLOSTP(NSTEP)=0.
      RLOSTA(NSTEP)=0.
      RLOSTF(NSTEP)=0.
      IND(NSTEP)=0
      ISEQ(NSTEP)=0
 490  CONTINUE
      DO 492 NSTEP = 1,NA
 492  ISEQ(NSTEP)=1
      DO 493 K=1,KEN
 493  SPECN(K)=0.
      DO 494 K=1,KEP1
 494  SPECP(K)=0.
      DO 495 K=1,KEA1
 495  SPECA(K)=0.
      DO 496 K=1,KEF
 496  SPECF(K)=0.
      DO 497 K=1,KEG
 497  SPECG(K)=0.
      DO 498 K=1,KEG
 498  SPECQ(K)=0.
      IF(NROW.EQ.NZ1) GOTO 510
      DO 499 K=1,IPOPS
      POP (K)=0.
      POPN(K)=0.
      POPP(K)=0.
      POPA(K)=0.
      PENTRY(K)=0.
 499  CONTINUE
 510  SIGMAX=0.
C
C*********************************************************************
C                                    SEQUENCE LOOP (NA)
      DO 800 NSTEP = 1,NA
      DO 630 ISO=1,NUMISO
      IF(ISOA(ISO).EQ.NSTEP .AND. ISOZ(ISO).EQ.NROW) GOTO 631
 630  CONTINUE
      ISO=0
 631  CONTINUE
C
      APAR=ZPAR+IN(NSTEP)
      IF(NSTEP.EQ.1) GOTO 635
      IF(MOD(NSTEP,3).NE.0 .OR. ITRANS.EQ.0) GOTO 438
 635                            IUNIT=21
      IF(LPRINT.LT.3)
     1WRITE(6,636) APAR,ZPAR,SIGLOW,IUNIT,IOWKB,KEN,KEP,KEA,LN,LP,LA,DEL
 636  FORMAT(2F6.0,F9.4,8I6/' DEL',5F8.1)
             CALL TRANSM(APAR,ZPAR,SIGLOW,TLOW,IPTR,IUNIT,IOWKB,QT)
      IF(LPRINT.LT.3)
     1WRITE(6,475) APAR,ZPAR,SIGLOW,IUNIT,DEL(2),DEL(3),DEL(4)
 475  FORMAT(' TRANSM  APAR,ZPAR,SIGLOW,UNIT,DEL',2F6.1,E10.2,I4,3F8.1)
C                         TRANSM CALLS OVER2(WHICH CALLS TLD) OR OWKB
C
 438                                   CALL SEQIN(NROW,NSTEP,LDBM,
     1                                         CST,CLD)
C
C
 640  APAR=ZPAR+IN(NSTEP)
      DO 641 K=1,5
 641  SHEL(K)=SHELK(K)
                         CALL ROT (NSTEP,NROW,APAR,ZPAR,SHEL,ITSTRT,NOF)
C
 660  IF(NSTEP.GT.1)
     1SIGMAX=DMAX1(SIGMAX,SIG(NSTEP-1))
      DUM2=FRACT2*SIGMAX
      SUMLOW=DMAX1(ABSMIN,DUM2)
                                       CALL GRIDS(NROW,NSTEP,CUTPOP)
C                                  GRIDS CALLS SUMPOP,GETPUT,CUTOFF,TOT
      IF(ISEQ(NSTEP).EQ.0) GO TO 700
      CALL DENSTY(NROW,NSTEP,LDBM,LEP,LEVPRN,KHS)
C
      IF(NOLEP.EQ.1)IAMX=10000
      CALL EVA(NROW,NSTEP,CUTPOP,NOLEP,NOLJI,NOLJF,LEVPRN,LDBM,KHS)
      IF(LPRINT.LT.3)
     1WRITE(6,1006) NSTEP
 1006 FORMAT( ' EVA-,STEP',I3)
C                                      EVA CALLS POPUL
C
                                       CALL SEQOUT(NROW,NSTEP,NOPUT)
                                       GOTO 800
 700  IF(NSTEP.EQ.NA) GOTO 800
      NSTEP1=NSTEP+1
      DO 710 NST=NSTEP1,NA
      IF(STORE(NROW,NST).GT.0.) GOTO 800
 710  CONTINUE
      GOTO 850
 800  CONTINUE
C
C******************** END SEQUENCE LOOP(NA) ****************************
C
 850  CALL ROWOUT(NROW,LOGUN,SIGMXZ,NOPUT,EXCIT0)
C 860  FORMAT(' IPAIR,CST,DELT,DELTA',I4,7F8.3)
      IF(NOPUT.EQ.1) GOTO 990
      DUM=0.
      IDUM=0
      CALL GETPUT(POP,IPOPS,10,IDUM,IDUM,IDUM,DUM,IDUM,IDUM,DUM,DUM)
      REWIND 10
      CALL GETPUT(POP,IPOPS,11,IDUM,IDUM,IDUM,DUM,IDUM,IDUM,DUM,DUM)
      REWIND 11
      IF(MOD(NROW,2).EQ.0)GOTO 990
      CALL GETPUT(POP,IPOPS,12,IDUM,IDUM,IDUM,DUM,IDUM,IDUM,DUM,DUM)
      REWIND 12
      CALL GETPUT(POP,IPOPS,13,IDUM,IDUM,IDUM,DUM,IDUM,IDUM,DUM,DUM)
      REWIND 13
 990  SIGMXZ=DMAX1(SIGMXZ,SIGZ(NROW))
      DUM=FRACT2*SIGMXZ
      IF(SIGZ(NROW).LT.DUM) GOTO 1500
 1000 CONTINUE
C                                     END ROW LOOP
C***********************************************************************
 1500 IF(LPRINT .GT. 4) GOTO 1510
      WRITE(6,403)EQ,EQ1
 1510 CALL ENGOUT(LOGUN,IDISC,EXCIT0)
 2000 CONTINUE
      GOTO 400
C                                     END ELAB LOOP
C***********************************************************************
 3000 CALL CASOUT(LOGUN,IDISC)
	CALL LINKFUS(MPROJ,IZP,MTGT,IZT,NZ,NA)
      IF(IEOF.EQ.1) GOTO 9999
      GOTO 99
C***********************************************************************
C                                       EXIT
 9999 WRITE(6,9992)
 9992 FORMAT(///' Ende HIVAP ')
 9993 FORMAT('C ',18A4)
 9991 CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(22)
      CLOSE(23)
      CLOSE(5)
      CLOSE(6)
      CLOSE(9)
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(19)
      CLOSE(21)
      CLOSE(25)
      STOP
      END
      BLOCK DATA
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ELEM/ELEMNT/S/BLANK,SYMB,NU
      CHARACTER*4 ELEMNT(128)/
     1 ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     2 'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',
     3 'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     4 'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr',
     5 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     6 'Sb','Te',' J','Xe','Cs','Ba','La','Ce','Pr','Nd',
     7 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     8 'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     9 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     X 'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     Y 'Md','No','Lr','Rf','Ha','Sg','Ns','Hs','Mt','110/',
     Z '111/','112/','113/','114/','115/','116/','117/','118/','119/',
     A '120/','121/','122/','123/','124/','125/','126/','127/','128/'/
C
      CHARACTER*1 BLANK/' '/,       NU(15)/'1','2','3','4','5','6','7',
     1 '8','9','A','B','C','D','E','F'/,SYMB(7)/'G','N','P','A','F',
     2 '*','X'/
      END