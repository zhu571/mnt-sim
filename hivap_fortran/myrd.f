      SUBROUTINE MYRD(ARRAY,IR  ,IPR,W)
C
C  Unterprogramm zum Einlesen der Eingabedaten
C  Stand: 3.5.1994
C
C     IRC1 NUMBER OF VARIABLES TO BE READ INTO 'ARRAY'
C     IPR  LOG UNIT FOR PRINT OUT  (IF<6, NO PRINT OUT)
C     W    LOG UNIT TO READ FROM
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ARRAY(36)
      INTEGER*4 W5,W
      CHARACTER*1 AW(72),AST,COMM,BLANK
      DATA COMM,AST,BLANK/'C','*',' '/
C--------------------------------------
      IRC1=IR
      W5=W
      W5=MAX0(5,W5)
      IF(IRC1.EQ.0)RETURN
      IDUM=1
      ITEST=0
      IF(IRC1.LT.0)ITEST=1
      IRC1=IABS(IRC1)
      I=1
      IRC2=0
 10   IRC=0
C--------------------------------------
      CALL DECODN(ARRAY(I),AW,IRC,W5,72)
C--------------------------------------
      IF(IRC.EQ.99) GOTO 99
      IF(IPR.GT.5) WRITE(IPR,15) AW
 15   FORMAT(1X,72A1)
      IF(AW(1).EQ.AST .OR. (AW(1).EQ.COMM .AND. AW(2).EQ.BLANK))GOTO 10
      IF(ITEST.EQ.0) GOTO 18
      IF(IRC.GT.0) IDUM=ABS(ARRAY(1))+0.01
      IF(IDUM.EQ.0) RETURN
 18   IRC2=IRC2+IRC
      IF(IRC2.GE.IRC1) GOTO 20
      I=IRC2+1
      GOTO 10
 20   RETURN
 99   IRC1=99
      RETURN
      END