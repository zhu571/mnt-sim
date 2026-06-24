      SUBROUTINE SEQIN(NROW,NSTEP,LDBM, CST,CLD)
C
C  Unterprogramm zur Initialisierung spezifischer Variabler eines
C  gegebenen, aktuellen Mutterkerns fuer Verdampfung (Separations-
C  energien, Zustandsdichteparameter etc.)
C  Stand: 3.5.1994
C
      INCLUDE 'common.f'
!      data r0wr,avwr,aswr,akwr/1.153,0.04543,0.1246,0.15231/  !cancelled by N.Wang
!      data r0ts,avts,asts,akts/1.160,0.04385,0.1584,0.33200/
      data r0wr,avwr,aswr,akwr/1.153,0.04543,0.1355,0.1426/  !W. Reisdorf, Z. Phys. A 300, 227-238 (1981)
C                          INITIALISATIONS FOR CURRENT STEP
      IA=NSTEP
      IZ=NROW
      AVGF=0.
      DO    J=1,JDIM
        YPOPJ(J)=0.
        YERJ(J)=0.
        YFISJ(J)=0.
      ENDDO
      IF(IOPT.EQ.1) THEN
        DO    J=1,JDIM
         FLAND(J)=0.
        ENDDO
      ENDIF   ! IOPT=1
      DO 45 K=1,5
      AVESP(K)=0.
 45   SSP(K)=0.
      DO 50 I=1,IPOPS
 50   TABLE(I)=BLANK
      IDUM=2*LN+1
      DO 55 I=1,IDUM
 55   DJN(I)=0.
      IDUM=2*LP+1
      DO 60 I=1,IDUM
 60   DJP(I)=0.
      IDUM=2*LA+1
      DO 65 I=1,IDUM
 65   DJA(I)=0.
      IF(LPRINT.lt.3) then
        DO 70 K=1,KEN
 70     SPECN(K)=0.
        DO 72 K=1,KEP1
 72     SPECP(K)=0.
        DO 74 K=1,KEA1
 74     SPECA(K)=0.
        DO 76 K=1,KEF
 76     SPECF(K)=0.
        DO 78 K=1,KEG
 78     SPECG(K)=0.
        DO 79 K=1,KEG
 79     SPECQ(K)=0.
      endif  ! lprint<3
C-----------------------------------------------------------------------
C                                 Q-VALUES FOR CURRENT STEP
      Q(1)=0.
      Q(5)=0.
      DO 100 K=2,4
 100  Q(K)=-BE(NROW,NSTEP,K-1)
C-----------------------------------------------------------------------
C                                  LEVEL DENSITY PARAMETERS
      IF(BARFAC.LE.0) BARFAC=1.
      IF(DELT.LE.0.) DELT=12.
      IF(IPAIR.EQ.0)DELT=0.
      IF(TZERO.LE.0.) TZERO=0.2
C ------------------------
      IAF=KZ(NSTEP)+IN(NSTEP)
      AF=IAF
      ZF=KZ(NSTEP)
C ---------------------------------------------------
      IF(LDBM.EQ.(-1)) then
        if(cld.le.0.) cld=1.
        AL(1)=ALS(IA,IZ)
        AL(2)=ALS(IA+1,IZ)
        AL(3)=ALS(IA,IZ+1)
        AL(4)=ALS(IA+2,IZ+2)
        AL(5)=CLD*AL(1)
      endif   ! Ldbm=-1
c     -----------------------
      IF(LDBM.eq.0) then
        if(cld.le.0.) cld=1.
        if(cst.le.0.) cst=10.
        AL(1)=AF/CST
        AL(2)=(AF-1.)/CST
        AL(3)=AL(2)
        AL(4)=(AF-4.)/CST
        AL(5)=CLD*AL(1)
      endif   ! ldbm=0
c     -----------------------
      if(ldbm.gt.0) then
          r0ld=r0wr
          avld=avwr
          asld=aswr
          akld=akwr

        if(abs(CST-1.15).lt.0.25) r0ld=cst  
c
        AL(1)=avld*r0ld**3*AF+asld*r0ld**2*(AF**0.666667)+
     1        akld*r0ld*(AF**0.33333)
        AF1=AF-1.
        AL(2)=avld*r0ld**3*AF1+asld*r0ld**2*(AF1**0.666667)+
     1        akld*r0ld*(AF1**0.33333)
        AL(3)=AL(2)
        AF1=AF-4.
        AL(4)=avld*r0ld**3*AF1+asld*r0ld**2*(AF1**0.666667)+
     1        akld*r0ld*(AF1**0.33333)
        CALL DROPB (ZF,AF,INDEX,DPBAR,Y,1,0)
        AL(5)=avld*r0ld**3*AF+asld*r0ld**2*(AF**0.666667)*BS+
     1        akld*r0ld*(AF**0.33333)*BK
 
      endif   ! ldbm>0

c     -----------------------
 144  IF(IPAIR.ne.4) then
        ODD1=MOD(KZ(NSTEP),2)
        ODD2=MOD(IN(NSTEP),2)
        DE=ABS(DELT    )/SQRT(AF)
        DELTA(1)=(1.-ODD1-ODD2)*DE
        DELTA(3)=(ODD1-ODD2)*DE
        DELTA(2)=(ODD2-ODD1)*DE
        DELTA(4)=DELTA(1)
        DELTA(5)=DELTA(1)
        IF(IPAIR   .eq.2 ) then
          DO     I=1,5
          DELTA(I)=DELTA(I)-DE
          enddo
        endif ! ipair=2
      endif  ! ipair not 4
C -------------------------------
      IF(ipair.eq.4 .or. LDBM.eq.(-1)) then
        DELTA(1)=DELTAS(IA,IZ)
        DELTA(2)=DELTAS(IA+1,IZ)
        DELTA(3)=DELTAS(IA,IZ+1)
        DELTA(4)=DELTAS(IA+2,IZ+2)
        DELTA(5)=DELTAS(IA,IZ)
      endif  ! ipair=4
c     ---------------------------------------------------------
 240  IF(LPRINT.LT.2) WRITE(6,245) IPAIR,DELT,CST,AL,DELTA
 245  FORMAT(' SEQIN:PAIR,DELT,CST',I4,F7.2,F8.3,
     1' AL=',5F7.2,' DELTA=',5F6.2)
      IF(LPRINT.LT.2 .AND. LDBM.GE.1)WRITE(6,250)CLD,BS,BK,DPBAR,AF,ZF,Y
 250  FORMAT(' AF/AN=',F8.3,2X,' BS,BK=',2F9.4,' DPBAR=',F7.2,
     1   ' A=',F5.0,' Z=',F5.0,' Y=',F8.3)
c     ---------------------------------------------------------------
c                                     shell corrections
      IF(ISHELL.gt.0) then
        SHELK(1)=SHELLS(IA,IZ)
        SHELK(2)=SHELLS(IA+1,IZ)
        SHELK(3)=SHELLS(IA,IZ+1)
        SHELK(4)=SHELLS(IA+2,IZ+2)
        SHELK(5)=0.
        IF(ISHELF.EQ.1) SHELK(5)=SHELK(1)
      IF(ISHELF.EQ.2) SHELK(5)=BFS(NSTEP,NROW)-BFLDM(NSTEP,NROW)*BARFAC+
     1                         SHELLS(NSTEP,NROW)
        IF(LPRINT.LT.3)
     1    WRITE(6,265) SHELK
 265    FORMAT(' SHELL CORRECTIONS',4X,5F8.1)
      endif ! shell corections
c     --------------------------------------------------------------
 300  RETURN
      END