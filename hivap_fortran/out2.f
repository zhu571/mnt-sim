      SUBROUTINE OUT2(POP,M,N,M0,N0,ISUPP,DIG1,MXDD,INV,LIN)
C
C  Druck-Routine fuer zweidimensionales Array, modifiziert am 18.11.88
C  Stand: 3.5.1994
C
C     MODIFIED TO ALLOW CALL BY VALUE OF DIG AND MXD 18 NOV 86

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION POP(1),AM(15)
      CHARACTER*1 ARR(128)
      COMMON/S/BLANK,SYMB,NU
      CHARACTER*1 BLANK,NU(15),SYMB(7)
      CHARACTER*4 DASH(2)
      DATA DASH/' M> ',' N> '/
C---------------------------------------------------------------------
      AMX=0.
      MN=M*N
      DO 10 I=1,MN
 10   AMX=MAX(AMX,POP(I))
      IF(AMX.LE.0.) RETURN
C
      DIG=DIG1
      IF(DIG.LE.1.)DIG=2.
      MXD=MXDD
      IF(MXD.LE.0)MXD=5
      MXD=MIN0(15,MXD)
C
      AM(1)=AMX
      IF(LIN.EQ.1) GOTO 30
      DO 20 I=2,15
 20   AM(I)=AM(I-1)/DIG
      GOTO 35
 30   DIG=MXD
      DUM=AMX/DIG
      DO 32 I=2,15
 32   AM(I)=AM(I-1)-DUM
C
 35   IH=MAX0(M,N)
      IV=MIN0(M,N)
      IDH=1
      IHH=IH
      IF(IH.LT.128) GOTO 50
 40   IDH=IDH+1
      IHH=IH /IDH
      IF(IHH.GT.128)GOTO 40
 50   IDV=1
      IVV=IV
      IF(IV.LT.60) GOTO 60
 55   IDV=IDV+1
      IVV=IV /IDV
      IF(IVV.GT.60) GOTO 55
C
 60   IF(IH.EQ.IV) GOTO 65
      IF(IH.EQ.N) GOTO 70
C---------------------------------------- M IS HORIZONTAL
 65   IPARH=1
      IF(INV.EQ.1 .OR. INV.EQ.3) IPARH=2
      IPARV=3
      IF(INV.GT.2) IPARV=4
      IDASH=1
      J00=N0-1
      GOTO 80
C---------------------------------------- N IS HORIZONTAL
 70   IPARH=3
      IF(INV.GE.2) IPARH=4
      IPARV=2
      IF(INV.EQ.1 .OR. INV.EQ.3) IPARV=1
      IDASH=2
      J00=M0-1
C-----------------------------------------------------------------
 80   DO 90 I=1,128
 90   ARR(I)=BLANK
      IF(IV.LT.3) GOTO 110
      DO 100 I=1,IHH,5
 100  ARR(I)=SYMB(6)
      WRITE(6,103) LIN,IDH,IDV,DIG,ISUPP
      IF(IHH.EQ.128) WRITE(6,101) DASH(IDASH),(ARR(I),I=1,IHH)
      IF(IHH.LT.128) WRITE(6,102) DASH(IDASH),(ARR(I),I=1,IHH)
 101  FORMAT(A4,128A1)
 102  FORMAT(A4,1X,127A1)
 103  FORMAT(/' OUT2:  LINFLAG',I2,4X,'STEPS H,V',2I2,4X,'CUT-PAR',
     1       F8.2,4X,'ZERO-LINE SUPPRESS FLAG',I2)
C-----------------------------------------------------------------
 110  MXD1=MXD+1
      DO 200 J1=1,IV,IDV
      IF(IPARV.EQ.1) J=J1
      IF(IPARV.EQ.2) J=IV+1-J1
      IF(IPARV.EQ.3) J=(J1-1)*IH
      IF(IPARV.EQ.4) J=(IV-J1)*IH
      L=0
      I0=0
      DO 180 I1=1,IH,IDH
      I0=I0+1
      IF(IPARH.EQ.1) I=I1
      IF(IPARH.EQ.2) I=IH+1-I1
      IF(IPARH.EQ.3) I=(I1-1)*IV
      IF(IPARH.EQ.4) I=(IH-I1)*IV
      JI=J+I
      ARR(I0)=BLANK
      DUM=POP(JI)
      IF(DUM.LE.0.)GOTO 180
      DO 170 MD=1,MXD
      IF(DUM.GE.AM(MD)) GOTO 175
 170  CONTINUE
      GOTO 180
 175  L=L+1
      ARR(I0)=NU(MXD1-MD)
 180  CONTINUE
      IF(L.EQ.0 .AND. ISUPP.EQ.1) GO TO 200
      J2=IABS(J00+J1)
      IF(IHH.LT.128) GOTO 192
      WRITE(6,190) J2,(ARR(I),I=1,I0)
 190  FORMAT(I4,128A1)
      GOTO 200
 192  WRITE(6,193) J2,(ARR(I),I=1,I0)
 193  FORMAT(I4,1X,127A1)
 200  CONTINUE
      RETURN
      END