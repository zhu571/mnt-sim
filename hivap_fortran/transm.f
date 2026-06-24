      SUBROUTINE TRANSM(APAR,ZPAR,SIGLOW,TLOW,IPRNT,IUNIT,IOWKB,Q2)
C
C  Liest Transmissionskoeffizienten fuer n,p,alphas aus Datei
C  (logische Einheit 21) oder nimmt Werte aus Rechnungen ( ruft OVER2
C  oder WKB auf) und schreibt sie wahlweise auf Datei)
C  Allokierung von Einheit 21 fuer rzri6f geandert (25.4.94 FPH)
C  Stand: 4.5.1994
C
C --------------------------------------------------------------------
C     IF IPRNT=2 WRITE SIGINV ON LOG UNIT 11 (FOR EVAP4)
C     IOWKB=0  OM FOR N,P,A ,WILL WRITE ON DISC IF IOVER>0
C          =1  WKB FOR N,P,A ,WILL NOT WRITE ON DISC OR READ FROM IT
C          =2  OM FOR N,P   WKB FOR A ,WILL WRITE ON DISC IF IOVER>0
C     IOVER=0  CALCULATE, NO READ, NO WRITE
C           1  READ IF ON DISC,ELSE CALC AND WRITE ON DISC,EXCEPT IF
C              IOWKB=1
C           2  NO READ,CALC AND WRITE ON DISC
C     DEL(2)  INPUT:<0  WILL ASK FOR POT PARAMETERS FOR N
C     DEL(3)  INPUT:<0  WILL ASK FOR POT PARAMETERS FOR P
C     DEL(4)  INPUT:<0  WILL ASK FOR POT PARAMETERS FOR A
C --------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/TRNS/TCOFN(25,13),TCOFP(25,13),TCOFA(50,21),DEL(5),
     1            KEN,KEP,KEA,LN,LP,LA,IOVER
      COMMON/TLJ/T(2,31),V(15),V1(3)
      COMMON/POTS/ VRN,R0RN,ADIFRN,VIN,R0IN,ADIFIN,
     1             VRP,R0RP,ADIFRP,VIP,R0IP,ADIFIP,RCLMBP,CBFACP,
     2             VRA,R0RA,ADIFRA,VIA,R0IA,ADIFIA,RCLMBA,CBFACA,QQ2
C --------------------------------------------------------------------
      DIMENSION SIG(50)
C      DOUBLE PRECISION T,V,V1
C --------------------------------------------------------------------
C     Aenderungen 25.4.1994  FPH
      DIMENSION IX(3),JX(3)
      CHARACTER*12  DATNAM
      CHARACTER*1   MASSNAM
      CHARACTER*3   MASSWERT
      CHARACTER*1   ZNAM
      CHARACTER*3   ZWERT
      CHARACTER*4   DATSPEZ
      CHARACTER*1   MZX(3)
      CHARACTER*1   ZZX(3)
      DATA MASSNAM/'m'/,ZNAM/'z'/,DATSPEZ/'.dat'/
C
      AMPAR1=0.1
      IUNIT=21
      KDIM=25
      KDIMA=50
      JEMAX=48
C     ----------------------------------
C     Eingefuegt 25.4.1994   FPH
C       Allokierung Dateien "mAPARzZPAR.dat" fuer Einheit 21
C       Umrechnung IX(I) --> JX(I)  CHAR(48,48,50...) = '0,1,2...
      NUFA=1
      IAPAR=APAR+0.01
      IZPAR=ZPAR+0.01
  11  IF(NUFA.EQ.1) IL=IAPAR
      IF(NUFA.EQ.2) IL=IZPAR
      DO 10 J=3,1,-1
      I=J-1
      IX(J)=IL/10**I
      JX(J)=48+IX(J)
      IF(NUFA.EQ.1) MZX(J)=CHAR(JX(J))
      IF(NUFA.EQ.2) ZZX(J)=CHAR(JX(J))
      IL=IL-IX(J)*10**I
  10  CONTINUE
      NUFA=NUFA+1
      IF(NUFA.LT.3) GOTO 11
      MASSWERT=MZX(3)//MZX(2)//MZX(1)
      ZWERT=ZZX(3)//ZZX(2)//ZZX(1)
      DATNAM=MASSNAM//MASSWERT//ZNAM//ZWERT//DATSPEZ
C     WRITE(6,*) DATNAM
C
      IREAD=0
      IWRIT=0
C
      IF(IOWKB.NE.1) THEN
C       -------------------
        IF(IOVER.EQ.1) THEN
           CALL HALLOC(DATNAM)
           IREAD=1
           IWRIT=0
        ENDIF  ! IOVER=1
C       -------------------
      IF(IOVER.EQ.2) THEN
           CALL HALLOC(DATNAM)
           IWRIT=1
           IREAD=0
        ENDIF   ! IOVER=2
C       ------------------
      ENDIF  ! IOWKB NOT 1
C
C************************************  NEUTRONS  **********************
      FACTOR=2.
      IF(IUNIT.GT.4)GOTO 1000
      GOTO(1000,2000,3000,9999),IUNIT
 1000 K=1
      IF(DEL(2).LT.0.)K=-1
      IF(KEN.EQ.0)KEN=25
      IF(LN.EQ.0)LN=12
      IF(IREAD.EQ.1) GOTO 22
      DUMF=DMAX1(DEL(2),AMPAR1)
      IF(IOWKB.NE.1)
     1CALL OVER2(K,KEN,LN,TCOFN,APAR,ZPAR,DUMF,SIGLOW,TLOW,KDIM,IPRNT)
      IF(IOWKB.EQ.1)
     1CALL OWKB (K,KEN,LN,TCOFN,APAR,ZPAR,DUMF,SIGLOW,TLOW,KDIM,Q2,
     2           IPRNT)
      DEL(2)=DUMF
      IF(IWRIT.EQ.0 .OR. IOWKB.EQ.1) GOTO 310
      WRITE (6,1001) K
 1001 FORMAT(' WRITE ON DISC  K=',I3)
      WRITE (IUNIT,212) TCOFN,DEL(2)
      GOTO 310
C
 22   IUN=IUNIT
      READ (IUNIT,212) TCOFN,DEL(2)
 212  FORMAT(6E12.4)
C
 310  IF(IPRNT.EQ.0) GOTO 3100
      WRITE (6,315) APAR,ZPAR,DUMF
 315  FORMAT(' TCOF ',6X,3F10.0/)
      DO 330 L=1,LN
      L0=L-1
      WRITE (6,320) L0
 320  FORMAT(I12)
      WRITE (6,325) (TCOFN(IE,L),IE=1,KEN)
 325  FORMAT(10E12.4)
 330  CONTINUE
      CALL OUT2(TCOFN,KEN,LN,0,0,0,FACTOR,10,0,0)
C
 3100 LN1=LN+1
      LN0=LN-1
      XMP=1.
      XM=(APAR-XMP)/APAR
      CONST=31.42/(0.04783*XMP*(XM**2))
      DO 335 I=1,JEMAX
 335  SIG(I)=0.
      DO 1210 I=1,KEN
      L1=LN
      FL=LN+LN-1
      SUM=TCOFN(I,L1)*FL
      E=I
      E=E+DEL(2)-1.
      DO 1110  L=1,LN0
      FL=FL-2.
      L1=L1-1
      SUM=SUM+TCOFN(I,L1)*FL
      TCOFN(I,L1  )=TCOFN(I,L1)  +TCOFN(I,L1+1 )
 1110 CONTINUE
      TCOFN(I,LN+1)=2.*SUM
      SIG(I)=SUM*CONST/E
 1210 CONTINUE
C
C     SPIN FACTOR 2
      DO 1220 I=1,KEN
      DO 1220 L=1,LN
 1220 TCOFN(I,L)=2.*TCOFN(I,L)
C
      IF(IPRNT.EQ.0) GOTO 2000
      WRITE(6,340) K,DEL(2)
 340  FORMAT(/' TOTAL X-SECTIONS',I4,4X,'STARTING WITH',F8.1,' MEV')
      WRITE(6,345) (SIG(I),I=1,KEN)
 345  FORMAT(10F10.4)
      IF(IPRNT.NE.2) GOTO 2000
      WRITE(11,1228)
 1228 FORMAT(' OUTPUT TRANSM FORMAT 1 0 0 '/
     1 ' INVERSE CROSS SECTIONS NEUTRONS')
      WRITE(11,1215) (SIG(I),I=1,JEMAX)
      DUM=1.
      WRITE(11,1231) DEL(2),DUM
C
C**************************************  PROTONS  **********************
 2000 IF(IUNIT.EQ.1) GOTO 9999
      K=2
      IF(DEL(3).LT.0.) K=-2
      IF(KEP.EQ.0)KEP=25
      IF(LP.EQ.0) LP=12
      IF(IREAD.EQ.1) GOTO 251
      DUMF=1.
      IF(IOWKB.NE.1)
     1CALL OVER2(K,KEP,LP,TCOFP,APAR,ZPAR,DUMF,SIGLOW,TLOW,KDIM,IPRNT)
      IF(IOWKB.EQ.1)
     1CALL OWKB (K,KEP,LP,TCOFP,APAR,ZPAR,DUMF,SIGLOW,TLOW,KDIM,Q2,
     2           IPRNT)
      DEL(3)=DUMF
      IF(IWRIT.EQ.0 .OR. IOWKB.EQ.1) GOTO 420
      WRITE (6,1001) K
      WRITE (IUNIT,212) TCOFP,DEL(3)
      GOTO 420
C
  251 READ(IUNIT,212) TCOFP,DEL(3)
C
 420  IF(IPRNT.EQ.0) GOTO 3200
      WRITE (6,315) APAR,ZPAR,DUMF
      DO 430 L=1,LP
      L0=L-1
      WRITE (6,320) L0
      WRITE (6,325) (TCOFP(IE,L),IE=1,KEP)
 430  CONTINUE
      CALL OUT2(TCOFP,KEP,LP,0,0,0,FACTOR,10,0,0)
 3200 LP1=LP+1
      LP0=LP-1
      XMP=1.
      XM=(APAR-XMP)/APAR
      CONST=31.42/(0.04783*XMP*(XM**2))
      DO 435 I=1,JEMAX
 435  SIG(I)=0.
      DO 1211 I=1,KEP
      FL=LP+LP-1
      L1=LP
      SUM=TCOFP(I,L1)*FL
      E=I
      E=E+DEL(3)-1.
      DO 1111  L=1,LP0
      FL=FL-2.
      L1=L1-1
      SUM=SUM+TCOFP(I,L1)*FL
 1111 TCOFP(I,L1)  =TCOFP(I,L1)  +TCOFP(I,L1+1 )
      TCOFP(I,LP+1)=2.*SUM
      SIG(I)=SUM*CONST/E
 1211 CONTINUE
C
C     SPIN FACTOR 2
      DO 1221 I=1,KEP
      DO 1221 L=1,LP
 1221 TCOFP(I,L)=2.*TCOFP(I,L)
C
 1215 FORMAT(6E12.4)
      IF(IPRNT.EQ.0) GOTO 3000
      WRITE(6,340) K,DEL(3)
      WRITE(6,345) (SIG(I),I=1,KEP)
      IF(IPRNT.NE.2) GOTO 3000
      WRITE(11,1229)
 1229 FORMAT(' INVERSE CROSS SECTIONS PROTONS')
      WRITE(11,1215) (SIG(I),I=1,JEMAX)
      DUM=1.
      WRITE(11,1231) DEL(3),DUM
C
C***************************************  ALPHAS  ********************
 3000 IF(IUNIT.EQ.2) GOTO 9999
      K=3
      IF(DEL(4).LT.0.)K=-3
      IF(KEA.EQ.0)KEA=50
      IF(LA.EQ.0)LA=20
      IF(IREAD.EQ.1) GOTO 253
      DUMF=2.
      IF(IOWKB.EQ.0)
     1CALL OVER2(K,KEA,LA,TCOFA,APAR,ZPAR,DUMF,SIGLOW,TLOW,KDIMA,IPRNT)
      IF(IOWKB.GT.0)
     1CALL OWKB (K,KEA,LA,TCOFA,APAR,ZPAR,DUMF,SIGLOW,TLOW,KDIMA,Q2,
     2           IPRNT)
      DEL(4)=DUMF
      IF(IWRIT.EQ.0 .OR. IOWKB.EQ.1) GOTO 530
      WRITE (6,1001) K
      WRITE (IUNIT,212) TCOFA,DEL(4)
      GOTO 530
C
 253  READ (IUNIT,212) TCOFA,DEL(4)
C
 530  IF(IPRNT.EQ.0) GOTO 3400
      WRITE (6,315) APAR,ZPAR,DUMF
      DO 540 L=1,LA
      L0=L-1
      WRITE (6,320) L0
      WRITE (6,325) (TCOFA(IE,L),IE=1,KEA)
 540  CONTINUE
      CALL OUT2(TCOFA,KEA,LA,0,0,0,FACTOR,10,0,0)
 3400 LA1=LA+1
      LA0=LA-1
      XMP=4.
      XM=(APAR-XMP)/APAR
      CONST=31.42/(0.04783*XMP*(XM**2))
      DO 535 I=1,JEMAX
 535  SIG(I)=0.
      DO 1212 I=1,KEA
      FL=LA+LA-1
      L1=LA
      SUM=TCOFA(I,L1)*FL
      E=I
      E=E+DEL(4)-1.
      DO 1112  L=1,LA0
      FL=FL-2.
      L1=L1-1
      SUM=SUM+TCOFA(I,L1)*FL
 1112 TCOFA(I,L1  )=TCOFA(I,L1)  +TCOFA(I,L1+1 )
      TCOFA(I,LA+1)=   SUM
      SIG(I)=SUM*CONST/E
 1212 CONTINUE
      IF(IPRNT.EQ.0) GOTO 9999
      WRITE(6,340) K,DEL(4)
      WRITE(6,345) (SIG(I),I=1,KEA)
      IF(IPRNT.NE.2) GOTO 9999
      WRITE(11,1230)
 1230 FORMAT(' INVERSE CROSS SECTION ALPHAS')
      WRITE(11,1215) (SIG(I),I=1,JEMAX)
      DUM=1.
      WRITE(11,1231) DEL(4),DUM
 1231 FORMAT(/' X0,DX=',2F8.3/)
 1232 FORMAT(/)
 9999 IF(IWRIT.EQ.1 .OR. IREAD.EQ.1)  THEN
      CLOSE(21)
      ENDIF
      IF(IPRNT.EQ.2) WRITE(11,1232)
      RETURN
      END