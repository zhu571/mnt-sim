      SUBROUTINE SEPENF(MTGT,IZT,MPROJ,IZP,QVALUE,NA,NZ,POP,KMAX1,
     1                 ISHELL,SHELLS,SHELL0,DELTAS,IPAIR,KHS)
C
C Berechnet Q-Werte, Separationsenergien fuer n,p,alphas, Spalt-
C barrieren, Schalenkorrektur-Energien
C ruft MSBEN, FISROT, liest Massentafel ( von logischer Einheit 9 )
C Einfuehrung experimenteller Paarungskorrekturen IPAIR=4 (6.10.1981)
C Einfuehrung gg-Kerne als Referenz fuer Schalenkorrekturen
C IPAIR=4 (18.12.1981)
C Stand: 5.5.1994
C
C     NEEDED EXTERNAL BE,MP,MC,NA,NZ,EXC,MTGT,IZT,MPROJ,IZP,QVALUE,
C     POP(DIM ABOUT 100),IND(DIM=15),ZCN,ACN,NUMB,MASSES,ELEMNT
C     'MASSES'  LOG UNIT MASS TABLE
C
C******* GET SEPARATION ENERGIES N,P,ALPHA,D,T *************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IND(27),POP(100),SHELLS(1),DELTAS(1),SMAS(27,27)
      DIMENSION MAGIC(6),MAGICN(6)
      COMMON/YR/EJAY(1500),RATIOS(6),RATIOF(6),FINERT,FINERF,
     &          JDIM,JUPYR,IYR,IRAST,INERF,NOROTF
      COMMON/BIND/BE(25,25,5),EXC(27,27),DELT,MC,MP,NUMB,MASSES
      COMMON/ELEM/ELEMNT(128)
      COMMON/FISRT/ROT0,A2MS,RKAPPA,X0,Y0,RKAPA4
      COMMON/BARR/BFLDM(25,25),BF(25,25),BAR0,BARFAC,IBF,IFISRT
      CHARACTER*4 ELEMNT
      DATA MAGIC/114,82,64,50,28,20/
      DATA MAGICN/184,126,82,50,28,20/
      NZMX=25
      NAMX=25
      KAKZMX=25*25
      NZMX2=NZMX+2
      NAMX2=NAMX+2
      NA=MIN0(NAMX-1,NA)
      NZ=MIN0(NZMX-1,NZ)
C*********************************************************************
C
      IF(INIT.EQ.99) GOTO 15
      A2MS=17.9439
      R0MS=1.2249
      C3MS=0.863987/R0MS
      RKAPPA=1.7826
      NOROTF=0
      FINERF=1.
      RKAPA4=0.
 15   WRITE(6,106)
 106  FORMAT(/' Unterprogramm SEPENF')
      IF(IPAIR.NE.4) GOTO 84
      DO 80 KAKZ=1,KAKZMX
 80   DELTAS(KAKZ)=0.
C-------------------------------------
 84   IF(IFISRT.NE.1) GOTO 100
      CALL MYRD(POP,7,23,5)
      A2MS=POP(1)
      R0MS=POP(2)
      C3MS=POP(3)
      RKAPPA=POP(4)
      NOROTF=POP(5)
      FINERF=POP(6)
      RKAPA4=POP(7)
C-------------------------------------
 100  IF(A2MS.LE.0.)A2MS=17.9439
      IF(R0MS.LE.0.)R0MS=1.2249
      IF(C3MS.LE.0)C3MS=0.863987/R0MS
      IF(RKAPPA.LE.0.) RKAPPA=1.7826
      ROT0=51.827/(R0MS*R0MS)
      Y0=ROT0/A2MS
      X0=C3MS/(2.*A2MS)
      IF(FINERF.LE.0.01) FINERF=1.
      IF(IFISRT.EQ.2) GOTO 110
C      WRITE(6,105) A2MS,R0MS,C3MS,RKAPPA,ROT0,X0,Y0,FINERF,NOROTF,
C     1             RKAPA4
C 105  FORMAT(' FISROT PARAMETERS'/' A2MS=',F9.4,4X,'R0MS=',F8.4,4X,
C     1 'C3MS=',F8.4,4X,'KAPPA=',F8.4,4X,'ROT0,X0,Y0=',F 9.3,2F9.5/
C     2 'FINERF,NOROTF=',F8.3,I4,4X,'KAPPA4=',F9.4/)
 107  FORMAT(18I4)
      GOTO 115
 110  WRITE(6,111)
 111  FORMAT(/' benutzt Sierks Programm BARFIT ')
C-------------------------------------- BE'S WITH MS67
 115  IZCN=IZT+IZP
      IACN=MTGT+MPROJ
      INCN=IACN-IZCN
      ZCN=IZCN
      ACN=IACN
      ZP=IZP
      AP=MPROJ
      ZT=IZT
      AT=MTGT
      KMAX=KMAX1
      IF(KMAX.EQ.0)KMAX=3
      DO 120 KK=1,NZMX2
      DO 120 JJ=1,NAMX2
 120  EXC(KK,JJ)=0.
      CALL  MSBEN(ZCN,ACN,NA,NZ,DUM,DUM)
      DO 122 KK=1,NZMX2
      DO 122 JJ=1,NAMX2
      SMAS(KK,JJ)=0.
 122  EXC(KK,JJ)=0.
C-------------------------------------------- USER BE'S --
      IF(NUMB.LE.0)  GOTO 130
      CALL MYRD(POP,1,23,5)
      QVALUE=POP(1)
      NUMB1=4*NUMB
      CALL MYRD(POP,NUMB1,23,5)
      DO 123 L=1,NUMB1,4
      IZ=IZCN+1-POP(L)
      IA=IACN+2-IZ-POP(L+1)
      K=POP(L+2)+0.01
      IF(IZ.LE.NZMX .AND. IA.LE.NAMX .AND. K.LE.KMAX)
     1             BE(IZ,IA,K)=POP(L+3)
 123  CONTINUE
 125  FORMAT(3(3I4,F8.3))
      GOTO 25
C----------------------------------- QVAL WITH MS67 +SHELLS AND PAIRING
 130  MCC=MC
      MPP=MP
      MP=0
      MC=0
      CALL MSBEN(ZP,AP,0,0,EXCP ,SHLL)
      CALL MSBEN(ZT,AT,0,0,EXCT ,SHLL)
      CALL MSBEN(ZCN,ACN,0,0,EXCCN,SHLLCN)
      MP=MPP
      MC=MCC
      QVALUE=EXCP+EXCT-EXCCN
      IF(MASSES.EQ.0)
     1 WRITE(6,127) MP,MC,ZP,AP,EXCP,ZT,AT,EXCT,ZCN,ACN,EXCCN,QVALUE
 127  FORMAT(' MSBEN MP,MC=',2I2,3(2F6.0,F9.3,3X),'QVALUE=',F9.3)
C-------------------------------------------------------------------
C                 EXP MASSES FROM DISC(EXC),LD MASSES MS67(SMAS)
      IF(MASSES.EQ.0) GOTO 25
 132  MCC=MC
      MPP=MP
      MP=1
      MC=1
      DO 135 KK=1,NZMX2
      DO 135 JJ=1,NAMX2
 135  EXC(KK,JJ)=0.
      NA4=NA+4
      NZ4=NZ+4
C--------------------------
 140  READ(MASSES,*) IZZ,IADN,IAUP
      IF(IZZ.EQ.0) GOTO 152
      NPOW=IAUP-IADN+1
      READ(MASSES,*)(POP(I),I=1,NPOW)
 145  FORMAT(10F8.3)
      IF(IZZ.EQ.IZP) EXCP=POP(MPROJ-IADN+1)
      IF(IZZ.EQ.IZT) EXCT=POP(MTGT-IADN+1)
      IF(IZZ.EQ.IZCN) EXCCNT=POP(IACN-IADN+1)
      IZ=IZCN-IZZ+2
      RZZ=IZZ
      IF(IZ.GT.NZ4) GOTO 140
      IF(IZ.LT.1) GOTO 152
C--------------------------------
      DO 150  IA=1,NA4
      IAA=IZZ+INCN+1 + 1-IA
      RAA=IAA
      CALL MSBEN(RZZ,RAA,0,0,SMAS(IZ,IA),DUM)
C                       AMONIX HERE
      I=IAA-IADN+1
      IF(I.LT.1  .OR.  I.GT.NPOW) GOTO 150
      EXC(IZ,IA)=POP(I)
 150  CONTINUE
C--------------------------------
      GOTO 140
C---------------------------------- TEST OUTPUT FOR TESTSE
 152  IF(JDIM.NE.1) GOTO 160
      WRITE (6,153 )
 153  FORMAT(10X,'MASS EXCESSES TABLE AND LDM')
      NA0=NA+4
      NA2=MIN0(12,NA4)
      NA1=MIN0(NA0,13)
 154  DO 159 IZ = 1,NZ4
      IIZ=IZCN+1   +1-IZ
      IDUM=IIZ+INCN+1   +1
      DO 155   IA=1,NA4
 155  IND (IA)=IDUM-IA
      WRITE(6,290 ) (ELEMNT(IIZ),IND (NA1-IA),IA=1,NA2)
      WRITE(6,300 ) (EXC(IZ,NA1-IA),IA=1,NA2)
      WRITE(6,300 ) (SMAS(IZ,NA1-IA),IA=1,NA2)
 159  CONTINUE
      IF(NA1.EQ.NA0) GOTO 160
      NA2=MIN0(12,NA0-NA1)
      NA1=NA1+12
      NA1=MIN0(NA0,NA1)
      GOTO 154
C-----------------------------------   GET EXP Q-VALUE
 160  MP=MPP
      MC=MCC
      DUM=ABS(EXCCNT)
      ACNLD=EXCCN-SHLLCN
      IF(DUM.GT.0.0001)EXCCN=EXCCNT
      SHLLCN=EXCCN-ACNLD
      IF(MASSES.NE.5) REWIND MASSES
      QVALUE=EXCP+EXCT-EXCCN
      WRITE (6,164) QVALUE
 164  FORMAT(' Q-Wert',F10.3,1X,'MeV')
      IF(ISHELL.EQ.0) GOTO 25
C-----------------------------------   GET EXP BE'S
      DO 220  IZ=1,NZ
      DO 210  IA=1,NA
      DUM=ABS( EXC(IZ+1,IA+1) )
      IF(DUM.LT.0.0001) GOTO 210
      DUM=ABS(EXC(IZ+1,IA+2))
      IF(DUM.LT.0.0001) GOTO 170
      BE(IZ,IA,1)=8.072 + EXC(IZ+1,IA+2) - EXC(IZ+1,IA+1)
 170  DUM=ABS(EXC(IZ+2,IA+1))
      IF(DUM.LT.0.0001) GOTO 180
      BE(IZ,IA,2)=7.289 + EXC(IZ+2,IA+1) - EXC(IZ+1,IA+1)
 180  DUM=ABS(EXC(IZ+3,IA+3))
      IF(DUM.LT.0.0001) GOTO 190
      BE(IZ,IA,3)=2.425 + EXC(IZ+3,IA+3) - EXC(IZ+1,IA+1)
 190  DUM=ABS(EXC(IZ+2,IA+2))
      IF(DUM.LT.0.0001) GOTO 200
      BE(IZ,IA,4)=13.136 + EXC(IZ+2,IA+2) - EXC(IZ+1,IA+1)
 200  DUM=ABS(EXC(IZ+2,IA+3))
      IF(DUM.LT.0.0001) GOTO 210
      BE(IZ,IA,5)=14.950 + EXC(IZ+2,IA+3) - EXC(IZ+1,IA+1)
 210  CONTINUE
 220  CONTINUE
 230  WRITE(6,234)
 234  FORMAT(' Massen aus der Datei "mlz.dat" ')
      DO 240  I=1,100
 240  POP(I)=0.
C
C---------------------------------------------------
C                                    GET SHELLS AND DELTAS
 25   IFLAG=0
      IF(ISHELL.EQ.0 .AND. IPAIR.EQ.0) IFLAG=1
      IF(IFLAG.EQ.0) GOTO 2505
      ISHELL=2
      IPAIR =4
      IF( MASSES.EQ.0) ISHELL=1
      IF( MASSES.EQ.0) IPAIR=2
      GOTO 2510
 2505 IF(MASSES.EQ.0 .AND. ISHELL.EQ.2) ISHELL=1
      IF(IPAIR.EQ.4  .AND. ISHELL.NE.2) IPAIR=2
 2510 INZ=0
      INA=0
      MCC=MC
      MPP=MP
      MC=0
      MP=0
C-------------------------------
      NZ2=NZ+2
      DO 40 KZ=1,NZ2
      KZ1=(KZ-1)*NAMX2
      IZ=IZCN+1-KZ
      IZOD=MOD(IZ,2)
      Z=IZ
C---------------------
      NA2=NA+2
      DO 30 KA=1,NA2
      KAKZ=KZ1+KA
      IA=IACN+2-KZ-KA
      IN=IA-IZ
      INOD=MOD(IN,2)
      A=IA
      CALL MSBEN(Z,A,INA,INZ,CMASS,SHELL)
      SHELLS(KAKZ)=SHELL
      SMASS=CMASS-SHELL
      IF(IPAIR.NE.4) GOTO 26
      ODDEV=-DELT*(1-IZOD-INOD)/SQRT(A)
      SMASS=SMASS-ODDEV
26    IF(ISHELL.EQ.2 .AND. EXC(KZ+1,KA+1).NE.0.)SHELLS(KAKZ)=
     1                     EXC(KZ+1,KA+1)-SMASS

C----------- Added by C.W.Shen, use Shell from FRDM(1995) -------
!      CALL SHELL_FILE(Z,A,SHELL)
!      SHELLS(KAKZ)=SHELL
!      SMASS=CMASS-SHELL
!      IF(IPAIR.NE.4) GOTO 26
!      ODDEV=-DELT*(1-IZOD-INOD)/SQRT(A)
!      SMASS=SMASS-ODDEV
!26    CONTINUE
C---------------------------------------------------------------

      IF(IPAIR.NE.4) GOTO 27
      IDUM=IZOD+INOD
      IF(IDUM.EQ.KHS+KHS)DELTAS(KAKZ)=0.
      DUM1=EXC(KZ,KA+1)
      DUM2=EXC(KZ+2,KA+1)
      DUM3=EXC(KZ+1,KA)
      DUM4=EXC(KZ+1,KA+2)
      IF(ABS(DUM1).LT.0.0001) GOTO 2620
      IF(ABS(DUM2).LT.0.0001) GOTO 2620
      IF(ABS(DUM3).LT.0.0001) GOTO 2620
      IF(ABS(DUM4).LT.0.0001) GOTO 2620
      DUM1=DUM1-SMAS(KZ,KA+1)
      DUM2=DUM2-SMAS(KZ+2,KA+1)
      DUM3=DUM3-SMAS(KZ+1,KA)
      DUM4=DUM4-SMAS(KZ+1,KA+2)
      DUM5=0.25*(DUM1+DUM2+DUM3+DUM4)-SHELLS(KAKZ)
      IF(IZOD.EQ.KHS) GOTO 2610
      DELTAS(KAKZ)=0.5*(DUM1+DUM2)-SHELLS(KAKZ)
      IF(KHS.EQ.0) GOTO 2610
      DO 2605 MAG=1,6
      IF(IZ.EQ.MAGIC(MAG)) DELTAS(KAKZ)=0.
 2605 CONTINUE
 2610 IF(INOD.EQ.KHS) GOTO 2630
      DELTAN=0.5*(DUM3+DUM4)-SHELLS(KAKZ)
      DO 2615 MAG=1,6
      IF(IN.EQ.MAGICN(MAG)) DELTAN=0.
 2615 CONTINUE
      DELTAS(KAKZ)=DELTAS(KAKZ)+DELTAN
      GOTO 2630
 2620 DELTAS(KAKZ)=-ODDEV-(1-2*KHS)*DELT/SQRT(A)
C------------------------------------------------------------C
 2630 SHELLS(KAKZ)=SHELLS(KAKZ)+DELTAS(KAKZ)-KHS*DELT/SQRT(A)
C-----FOLLOWING FOR OVERALL SHIFT OF REFERENCE SURFACE ------C
 27   SHELLS(KAKZ)=SHELLS(KAKZ)+SHELL0
 30   CONTINUE
C---------------------
 40   CONTINUE
      IF(IFLAG.EQ.0) GOTO 50
      IF(ABS(SHELLS(1)).GT.0.001) SHLLCN=SHELLS(1)
      DO 44 KAKZ=1,KAKZMX
      DELTAS(KAKZ)=0.
 44   SHELLS(KAKZ)=0.
      ISHELL=0
      IPAIR=0
      IFLAG=0
C--------------------------------------
 50   MP=MPP
      MC=MCC
C------------------------------------------------------------------
      IF(ISHELL.EQ.1)WRITE(6,60) SHELL0
      IF(ISHELL.EQ.2)WRITE(6,61) SHELL0
 60   FORMAT(/9X,'MS SHELLS  SHELL0=',F8.2)
 61   FORMAT(/9X,'Exp. Schaleneffekte,  SHELL0=',F8.2)
      NA2=MIN0(12,NA)
      NA0=NA+1
      NA1=MIN0(NA0,13)
 65   DO 90 IZ = 1,NZ
      IIZ=IZCN+1-IZ
      IDUM=IACN-IZ+2
      IZ1=(IZ-1)*NAMX2+NA1
      DO 70    IA=1,NA
 70   IND (IA)=IDUM-IA
      WRITE(6,75  ) (ELEMNT(IIZ),IND (NA1-IA),IA=1,NA2)
 75   FORMAT(8X,12(2X,A4,I3,1X))
      WRITE(6,85  ) (SHELLS(IZ1-IA) ,IA=1,NA2)
 85   FORMAT(8X,12F10.2/)
 90   CONTINUE
      IF(NA1.EQ.NA0) GOTO 410
      NA2=MIN0(12,NA0-NA1)
      NA1=NA1+12
      NA1=MIN0(NA0,NA1)
      GOTO 65
C----------------------------------- GET F BARRIERS BF,BFLDM----
 410  AL=0.
      DUM=0.
      DO 420 IZ=1,NZ
      Z=ZCN+1-IZ
      IZEE=Z+0.01
      KZ1=(IZ-1)*NAMX2
      DO 415 IA=1,NA
      A=ACN+2-IZ-IA
      MASS=A+0.01
      KAKZ=KZ1+IA
      AN=A-Z
      IF(IFISRT.EQ.0) CALL FISROT(A,Z,AN,AL,DELR,DELSP,ERO,DUM)
      IF(IFISRT.EQ.0) BFLDM(IA,IZ)=DELSP-DELR
      IF(IFISRT.EQ.2) CALL BARFIT(IZEE,MASS,0,BFIS,DELR,ELMAX)
      IF(IFISRT.EQ.2) BFLDM(IA,IZ)=BFIS

C     ------------------------------ The following added by Ning Wang --------
      IF(IFISRT.EQ.3) BFLDM(IA,IZ) = BFMWS(A,Z)       !Modified Woods-Saxon
!      IF(IFISRT.EQ.3) BFLDM(IA,IZ) = BFEMPIRIC(A,Z)  !added by C.W.Shen
!      IF(IFISRT.EQ.3) BFLDM(IA,IZ) = BFSierk(A,Z)    !test Sierk's barrier
!	 IF(IFISRT.EQ.3) BFLDM(IA,IZ) = BFLiquid(A,Z)   !test Liquid Drop barrier

      IF(IBF.EQ.0)BF(IA,IZ)=BARFAC*BFLDM(IA,IZ)+BAR0
      IF(IBF.EQ.1)BF(IA,IZ)=BARFAC*BFLDM(IA,IZ)+BAR0 - SHELLS(KAKZ)
!      WRITE(*,'(X,A,2F7.1,2x,2F11.5)') "Z, A, Bfldm, Shell:",
!     &        Z,A,Bfldm(ia,iz),shells(kakz)
!     Here, BF(IA,IZ)=BFS(IA,IZ) in COMMON/BARR/ ------ N. WANG -------------
 415  CONTINUE
 420  CONTINUE
      IF(IBF.LT.2) GOTO 249
C     IF(IBF.EQ.2) READ(5,1070) BF(1,1)
      IF(IBF.EQ.2) CALL MYRD(BF,1,23,5)
 1070 FORMAT(9F8.3)
      DO 425 IZ=1,NZ
      DO 424 IA=1,NA
 424  BF(IA,IZ)=BF(1,1)
 425  CONTINUE
C     IF(IBF.EQ.3) READ(5,1070) ((BF(IA,IZ),IA=1,NA),IZ=1,NZ)
      IRC=NA*NZ
      IF(IBF.EQ.3) CALL MYRD(BF,IRC,23,5)
C---------------------------------------------------------------
 249  IF(MASSES.NE.0) GOTO 280
 250  WRITE(6,251)
 251  FORMAT(' MS LYSEKIL 1967')
      IF(MC)260,260,255
 255  WRITE (6,256 )
 256  FORMAT(' LIQUID DROP MASSES')
      GO TO 265
 260  WRITE (6,261 )
 261  FORMAT(' SHELL CORRECTED MASSES')
 265  IF (MP) 275,275,270
 270  WRITE (6,271 )
 271  FORMAT(' ZERO PAIRING')
      GO TO 280
 275  WRITE (6,276 )
 276  FORMAT(' WITH PAIRING')
C-------------------------------------
 280  WRITE(6,282) BARFAC,BAR0
 282  FORMAT(/10X,' BARFAC=',F8.3,6X,'BAR0=',F8.3)
      WRITE (6,281 )
 281  FORMAT(/10X,'Separationsenergien (BN,BP,BA) und Spaltbarrieren')
      NA0=NA+1
      NA2=MIN0(12,NA)
      NA1=MIN0(NA0,13)
 283  DO 310 IZ = 1,NZ
      IIZ=IZCN+1-IZ
      IDUM=IACN-IZ+2
      KZ1=(IZ-1)*NAMX2
      DO 285   IA=1,NA
 285  IND (IA)=IDUM-IA
      WRITE(6,290 ) (ELEMNT(IIZ),IND (NA1-IA),IA=1,NA2)
 290  FORMAT(8X,12(2X,A4,I3,1X))
      DO 295  K=1,KMAX
 295  WRITE(6,300 ) (BE(IZ,NA1-IA,K),IA=1,NA2)
 300  FORMAT(8X,12F10.2)
      WRITE(6,300)(BF(NA1-IA,IZ),IA=1,NA2)
      IF(IPAIR.EQ.4) WRITE(6,300)(DELTAS(KZ1+NA1-IA),IA=1,NA2)
 310  CONTINUE
      IF(NA1.EQ.NA0) GOTO 313
      NA2=MIN0(12,NA0-NA1)
      NA1=NA1+12
      NA1=MIN0(NA0,NA1)
      GOTO 283
C----------------------------------------------------
 313  QVALUE=EXCP+EXCT-EXCCN
      IF(ISHELL.EQ.0) SHELLS(1)=SHLLCN
      INIT=99
      RETURN
      END

C--------------------------------------------------------------------

      REAL*8 FUNCTION BFLiquid(A,Z)
      REAL*8 A,Z,I,X,U,N
 
	Ec0=0.7053*Z*Z/A**0.3333
	Es0=17.9439*(1-1.7826*((A-2*Z)/A)**2)*A**0.6666
	X=Ec0/(2*Es0)
	if(X.gt.0.3333.and.X.lt.0.6666) U = 0.38*(0.75-X)*Es0
	if(X.ge.0.6666.and.X.lt.1)      U = 0.83*(1-X)**3*Es0
      BFLiquid = U

	END 

C---------------------------------------------------------------------
      function BfSierk(A,Z)
!	A. Sierk, Phys. Rev. C33 (1986) 2039.

      IMPLICIT REAL*8 (A-H,O-Z)
	il=0
	IZEE=int(Z)
	MASS=int(A)
	CALL BARFIT(IZEE,MASS,0,BFIS,DELR,ELMAX)
	Bfsierk=bfis
	return
	end


C--------------------------------------------------------------------
C     Empiric Liquid Drop fission barrier, added by C.W.Shen
C     from: M. Dahlinger and D. Vermeulen, Nucl.Phys.A376(1982)94
C-----
      REAL*8 FUNCTION BFEMPIRIC(A,Z)
      REAL*8 A,Z,I,X,U,N
      N = A-Z
      I = (N-Z)/A
      X = Z*Z/49.22D0/A/(1-0.3803*I*I-20.489*I**4)
      U = 0.368-5.057*X+8.93*X*X-8.71*X**3
      BFEMPIRIC = 0.7322*Z*Z*10**(U)*A**(-1D0/3D0)

	END 

C--------------------------------------------------------------------
C------------ SHELL CORRECTION FROM FILE (Koura OR Moller) ----------
      SUBROUTINE SHELL_FILE(Z,A,SHELL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8   SHELL,TSHELL
      IZ = Z+0.1
      IA = A+0.1
      OPEN(93,FILE='shell.dat',STATUS='OLD')
      IEND = 0
      IOS = 0
      DO WHILE(IOS==0)
          READ(93,*,IOSTAT=IOS) NZ,NA,TSHELL
          IF(IOS==0) THEN
              IF(NZ.EQ.IZ .AND. NA.EQ.IA) THEN
                  SHELL = TSHELL
                  GOTO 10
              END IF
          END IF
      END DO
10    CLOSE(93)
      END
              
C     Modified Woods-Saxon potential for fission barrier, added by Ning Wang
 
      REAL*8 FUNCTION BFMWS(A,Z)
      IMPLICIT REAL*8 (A-H,O-Z)

	A1=INT(A/2.)
	Z1=INT(Z/2.)
	A2=A-A1
	Z2=Z-Z1

	VMIN=100000
	DO DIS=7,16,0.1
	V=VMWS(A1,Z1,A2,Z2,DIS)
      IF(V.GT.VMIN)THEN
	EXIT
	ELSE
	VMIN=V
	ENDIF
	ENDDO
	Bs=VMIN

	IF(DIS.GE.15)THEN
	BFMWS=0.0
	RETURN
	ENDIF
	
	VMAX=0
	DO DIS=16,7,-0.1
	V=VMWS(A1,Z1,A2,Z2,DIS)
      IF(V.LT.VMAX)THEN
	EXIT
	ELSE
	VMAX=V
	ENDIF
	ENDDO
	B0=VMAX
	
	BFMWS=B0-Bs
	IF(BFMWS.LT.0.1)BFMWS=0.0   

	RETURN

	END 

