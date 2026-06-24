      SUBROUTINE GRIDS(NROW,NSTEP,CUTPOP)
C
C  Bereitet die E-J - Bevoelkerung vor der Hauptrechnung (in EVA) zu;
C  ruft GETPUT, SUMPOP, TOT, CUTOFF auf
C  Stand 3.5.1994
C
      INCLUDE 'common.f'
!------- The following added by Kai Zhao (’‘ø≠) -----
	M=0 
	N=0 
	EXMAX=0.
!----------------------------------------------------
      IPRINT=LPRINT
      IF(NSTEP.EQ.1 .AND. LPRINT.LT.3) WRITE(6,5)
 5    FORMAT(1H1,' NEW ROW')
      IF(LPRINT.LT.2)
     1WRITE (6,10) NROW,NSTEP
 10   FORMAT(/' GRIDS  NROW,NSTEP',2I4)
C
C                                      GET P
C
 250  MN=M2*N2
      SUMN=0.
      DO 260 I=1,MN
 260  SUMN=SUMN+POPN(I)
      IF(LPRINT.LT.2)
     1WRITE(6,265) SUMN,SUMLOW
 265  FORMAT(' SUMN,SUMLOW',2E14.4)
      IDUM= STORE(NROW,NSTEP)
      IF(IDUM.EQ.0 .OR. IDUM.EQ.2) GOTO 275
      IUNIT=10+IU1
      IREAD=1
      CALL GETPUT(POPP,IPOPS,IUNIT,IREAD,M,N,EXMAX,MMZ,MMA,AACN,ZZCN)
      MN=N*M
      SUMP=0.
      DO 270 I=1,MN
 270  SUMP=SUMP+POPP(I)
      I=AACN+0.001
      II=ZZCN+0.001
      IF(LPRINT.LT.2)
     1WRITE(6,272) IUNIT,SUMP,M,N,EXMAX,MMA,ELEMNT(MMZ),I,ELEMNT(II)
 272  FORMAT(' GOT POPP FROM UNIT',I4,6X,'SUMP',E12.4,3X,'M,N,EXMAX',
     1        2I4,F8.1,4X,I4,A4,4X,'CN=',I4,A4)
      GOTO 280
 275  DO 276 I=1,IPOPS
 276  POPP(I)=0.
      SUMP=0.
 280  DO 282 I=1,IPOPS
 282  POPA(I)=0.
      CALL SUMPOP(N ,M ,EXMAX ,POPP,SUMP,IPOPS,N2,M2,EXMAX2,POPN,
     1                 SUMN  ,POPA,FLOST(NSTEP),IPRINT)
C
C                                      GET A
C
 300  IDUM= STORE(NROW,NSTEP)
      IF(IDUM.EQ.0 .OR. IDUM.EQ.1) GOTO 350
      IDIM=-1
      IUNIT=12+IU2
      K=4
      IREAD=1
      CALL GETPUT(POPN,IPOPS,IUNIT,IREAD,M,N,EXMAX,MMZ,MMA,AACN,ZZCN)
      MN=M*N
      SUMA=0.
      DO 310 I=1,MN
 310  SUMA=SUMA+POPN(I)
      IF(LPRINT.LT.2)
     1WRITE(6,312) IUNIT,SUMA,M,N,EXMAX,MMA,ELEMNT(MMZ)
 312  FORMAT(' GOT POPA FROM UNIT',I4,6X,'SUMA',E12.4,3X,'M,N,EXMAX',
     1        2I4,F8.1,4X,I4,A4)
      GOTO 360
 350  DO 355 I=1,IPOPS
 355  POPN(I)=0.
      SUMA=0.
 360  DO 362 I=1,IPOPS
 362  POP(I)=0.
      CALL SUMPOP(N2,M2,EXMAX2,POPA,SUMN,IPOPS,N,M ,EXMAX, POPN,
     1            SUMA       ,POP ,FLOST(NSTEP),IPRINT)
      N1=N
      M1=M
      EXMAX1=EXMAX
C
C                           TRIM NEW PARENT
C
      CALL TOT(N1,M1,POP,SPECJ,AMX,AVGJ(NSTEP),TOTAL(NSTEP),SPECE,
     1         AVGE(NSTEP),EXMAX1)
      IF(TOTAL(NSTEP).GT.SUMLOW) GOTO 450
      ISEQ(NSTEP)=0
      RETURN
 450  ALEVEL=0.
      CUT1=CUT
      PRCNT1=PRCNT
      IF(NROW.EQ.1 .AND. NSTEP.EQ.1) CUT1=CUT*0.1
      IF(NROW.EQ.1 .AND. NSTEP.EQ.1) PRCNT1=1.-(1.-PRCNT)*0.1
      CALL CUTOFF(N1,M1,POP,SPECJ,TOTAL,NSTEP,AMX,CUT1,PRCNT1,ICUT,
     1            TRIM(NSTEP),IPOPS,EXMAX1,ALEVEL,IPRINT)
      CALL TOT(N1,M1,POP,SPECJ,AMX,AVGJ(NSTEP),TOTAL(NSTEP),SPECE,
     1         AVGE(NSTEP),EXMAX1)
      IDUM=M1*N1
      DO 455 I=1,IDUM
 455  PENTRY(I)=POP(I)
      CUTPOP=ALEVEL/2.
      DUM=0.05*CUT*TOTAL(NSTEP)
      IF(ALEVEL.LE.0. .OR. CUTPOP.GT.DUM) CUTPOP=DUM
      DUM=0.
      MM1=MIN0(IEDIM,M1)
      DO 460 I=1,MM1
 460  DUM=DMAX1(DUM,SPECE(I))
      IAMX=DUM*10000./TOTAL(NSTEP)  +0.5
C
C                           SET LIMITS OF NEW DAUGHTERS
C                                                       N
 505  K=2
      N2=1
      M2=1
      EXMAX2=EXMAX1+Q(K)   -DEL(K)
      IF(EXMAX2.LE.0.)GOTO 518
      N22=MIN0(N1+5,JDIM)
      JZERO=JDIM*(K-1)+1
      DO 509 J1=1,N22
      JYR=N22-J1+JZERO
      IF(EJAY(JYR).LE.EXMAX2) GO TO 510
  509 CONTINUE
  510 N2=JYR-JZERO+1
C
  511 M22=EXMAX2+1.
      M2=M1+KEN
      M2=MIN0(M22,M2)
      IF((M2*N2).LE.IPOPS) GOTO 518
      M2=IPOPS/N2
 518  MN=M2*N2
      DO  530  I=1,IPOPS
 530  POPN(I)=0.
C                                                       P
C
 605  K=3
      N3=1
      M3=1
      EXMAX3=EXMAX1+Q(K)   -DEL(K)
      IF(EXMAX3.LE.0.)GOTO 618
      N33=MIN0(N1+5,JDIM)
      JZERO=JDIM*(K-1)+1
      DO 609 J1=1,N33
      JYR=N33-J1+JZERO
      IF(EJAY(JYR).LE.EXMAX3) GO TO 610
  609 CONTINUE
  610 N3=JYR-JZERO+1
C
  611 M33=EXMAX3+1.
      M3=M1+KEP
      M3=MIN0(M33,M3)
      IF((M3*N3).LE.IPOPS) GOTO 618
      M3=IPOPS/N3
 618  MN=M3*N3
      DO  630  I=1,IPOPS
 630  POPP(I)=0.
C
C                                                        A
  705 K=4
      N4=1
      M4=1
      EXMAX4=EXMAX1+Q(K)   -DEL(K)
      IF(EXMAX4.LE.0.) GOTO 718
      N44=MIN0(N1+5,JDIM)
      JZERO=JDIM*(K-1)+1
      DO 709 J1=1,N44
      JYR=N44-J1+JZERO
      IF(EJAY(JYR).LE.EXMAX4) GO TO 710
  709 CONTINUE
  710 N4=JYR-JZERO+1
C
  711 M44=EXMAX4+1.
      M4=M1+KEP
      M4=MIN0(M44,M4)
      IF((M4*N4).LE.IPOPS) GOTO 718
      M4=IPOPS/N4
 718  MN=M4*N4
      DO  730  I=1,IPOPS
 730  POPA(I)=0.
      IF(LPRINT.LT.2)
     1WRITE(6,740) M1,N1,EXMAX1,M2,N2,EXMAX2,M3,N3,EXMAX3,M4,N4,EXMAX4
 740  FORMAT(' DIMENSIONS',4(2I4,F8.1))
      RETURN
      END