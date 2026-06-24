      SUBROUTINE DECODN (ARRAY,A,IRC,LOGUN,NN)
C
C  Unterprogramm zum Einlesen der Eingabedaten
C  Stand: 3.5.1994
C
C     'A'      IS READ IN A-FORMAT AND/OR DECODED
C     'IRC'    OUTPUT:NUMBER OF NUMBERS FOUND
C              OUTPUT:99 END OF FILE FOUND
C              OUTPUT:-1 DECIMAL EXPONENT EXCEEDS 99
C     'LOGUN'  LOG.UNIT FROM WHICH 'A' IS READ
C              0  NO READING
C              1-4  =5
C     'NN'     'A' WILL BE READ AND/OR DECODED FROM A(1) TO A(N)
C              DEFAULT FOR ZERO IS 72
C---------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ARRAY(1)
      CHARACTER B,A(1),SY(15),MINUS,ESIGN,PERIOD,BLANK,COMMA,EQUAL,DOPP,
     1          SEMI
      CHARACTER SLASH
      EQUIVALENCE (SY(11),MINUS),(SY(14),ESIGN),(SY(12),PERIOD),(SY(15),
     1 COMMA)
      DATA SY,BLANK/'0','1','2','3','4','5','6','7','8','9','-','.','+',
     1 'E',',',' '/,EQUAL,DOPP  /'=',':'/,SLASH/'/'/,SEMI/';'/
C---------------------------------------------------------------------
C                                                         INITIALIZE
      N=NN
      IF(NN.LE.0) N=72
      N1=N+1
      IFLAG=0
      ISEP=1
      DEC=0.
      NUMB=0
      IE=0
      IEST=0
      S=1.
      IS=1
      IW=1
      IBLANK=1
      IW1=1
      L=1
      F=1.
      IF(LOGUN.LT.1) GOTO 17
      LOGUN1=MAX0(LOGUN,5)
C----------------------------------------------------------------------
C                                           READ UP TO N  SYMBOLS
      READ(LOGUN1,15,END=99) (A(I),I=1,N)
 15   FORMAT(72A1,50A1)
C----------------------------------------------------------------------
C                                                SYMBOL LOOP
 17   IRC=0
      DO 80 I=1,N1
      IF(I.EQ.N1) GOTO 25
      B=A(I)
C----------------------------- IS IT ONE OF THE SYMBOLS 'SY' (0 TO +)
      DO 20 K=1,13
      IF(B.EQ.SY(K) .AND.ISEP.EQ.1) GOTO 30
 20   CONTINUE
C----------------------------- NO IT'S A SEPARATION
 25   IF(B.EQ.ESIGN .AND. IBLANK.EQ.0) GOTO 60
      ISEP=0
      IF(B.EQ.COMMA .OR. B.EQ.BLANK .OR. B.EQ.EQUAL .OR. B.EQ.DOPP
     1   .OR. B.EQ.SLASH .OR. B.EQ.SEMI .OR. I.EQ.N1) ISEP=1
      IBLANK=IBLANK+1
      IEST=0
      IF(IBLANK.NE.1) GOTO 50
      IF(IFLAG.EQ.0) GOTO 28
      IF(ISEP.NE.1) GOTO 55
      IW=IW+1
      IFLAG=0
      GO TO 50
 28   IS=1
      S=1.
      L=1
      GOTO 50
C----------------------------- IS IT A NUMBER OR IS IT + - E .
 30   IBLANK=0
      IF(IEST.EQ.0) GOTO 32
      IEST=0
      L=3
 32   IF(K.GT.10) GOTO 70
C----------------------------- IT'S A NUMBER
      IFLAG=1
      GOTO (35,40,45),L
 35   NUMB=NUMB*10+K-1
      GOTO 80
 40   F=F*0.1
      DEC=DEC+(K-1)*F
      GOTO 80
 45   IE=IE*10 + K-1
      IF(IE.LT.99) GOTO 80
      IRC=-1
      RETURN
C---------------------------- IS IT END OF NUMBER-WORD
 50   IF(IW.EQ.IW1) GOTO 80
C---------------------------- YES IT IS
      R=NUMB
      D=10.**FLOAT(IE)
      IF(IS.LT.0)ARRAY(IW-1)=S*(R+DEC)/D
      IF(IS.GT.0)ARRAY(IW-1)=S*(R+DEC)*D
 55   IW1=IW
      DEC=0.
      F=1.
      NUMB=0
      IS=1
      S=1.
      IE=0
      L=1
      GOTO 80
C---------------------------- E-SIGN FLAG
 60   IEST=1
      GOTO 80
C---------------------------- IT'S . OR E OR MINUS
 70   IF(B.EQ.PERIOD) L=2
      IF(B.EQ.MINUS.AND.L.EQ.1) S=-1.
      IF(B.EQ.MINUS .AND. L.EQ.3) IS=-1
C----------------------------
 80   CONTINUE
C------------------------------------------------------------------
C                                         END OF SYMBOL LOOP
      IRC=IW-1
      RETURN
 99   IRC=99
      RETURN
      END