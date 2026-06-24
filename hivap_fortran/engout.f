      SUBROUTINE ENGOUT(LOGUN,IDISC,EXCIT0)
C
C  Steuert  Ausgabe
C  Stand: 3.5.1994
C
      INCLUDE 'common.f'
C
      DIMENSION SIGMAS(50)
C     DIMENSION IS NAUP+NZUP
      NEXC=NEXC+1
                    LOGUN=6
      IF(LPRINT. GT.4)LOGUN=22
      IACN =ACN
      IZCN =ZCN
      I0=IACN-(IACN/10)*10+1
      DO 20 I=1,10
      I1=I0-I
      IF(I1.LT.0)I1=I1+10
 20   IND(I)=I1
      NANZ=NA+NZ
      DO 30 N=1,NANZ
 30   SIGMAS(N)=0.
      DO 40  NROW=1,NZ
      DO 35  NSTEP=1,NA
      N=NSTEP+NROW-1
 35   SIGMAS(N)=SIGMAS(N)+XSECT(NROW,NSTEP)
 40   CONTINUE
C---------------------------------------------------------
      WRITE (LOGUN,50)
 50   FORMAT(/   ' HIVAP  ENGOUT'/)
      WRITE(LOGUN,55)(TITLE(I),I=1,18)
 55   FORMAT(18A4)
      IF(ELAB1.LE.0.) GOTO 58
      WRITE(LOGUN,57) ELAB1,ELAB,EXCIT1,EXCIT0,IACN,ELEMNT(IZCN)
 57   FORMAT(' AVERAGED CALCULATION  RANGE(LAB):',2F8.1,4X,'OR(EXCIT)',
     1        2F8.1/' CN=',I4,A4)
      GOTO 59
 58   WRITE (LOGUN,60) ELAB,EXCIT0,IACN,ELEMNT(IZCN)
 59   WRITE (LOGUN,65) SIGFUS,FLGRAZ,SIGEVA,SIGFIS
 60   FORMAT(4X,'ELAB=',F8.1,6X,'EXCIT=',F8.1,6X,'CN=',I4,A4)
 65   FORMAT(4X,'FUSION X-SECT=',E12.4,6X,'LCRIT=',F6.1/
     1 4X,'EVAP SO FAR',E12.4,6X,'FISSION',E12.4)
C-----------------------------------------------------
      ELABE(NEXC)=ELAB
      ELABE1(NEXC)=ELAB1
      SIGFUE(NEXC)=SIGFUS
      SIGEVE(NEXC)=SIGEVA
      SIGFIE(NEXC)=SIGFIS
      FLGRAE(NEXC)=FLGRAZ
C--------------------------------------------------------
      WRITE (LOGUN,105)
 105  FORMAT(/' MASS DISTRIBUTION')
      WRITE(LOGUN,106)(IND(I),I=1,10)
 106  FORMAT(I9,4I12,6X,5I12)
      DO 120 N=1,NANZ,10
      N1=MIN0(N+9,NANZ)
      MAS=IACN -N +1
      WRITE (LOGUN,115) MAS,(SIGMAS(N2),N2=N,N1)
 115  FORMAT(I4,5E12.4,6X,5E12.4)
 120  CONTINUE
C--------------------------------------------------------
      WRITE(LOGUN,125)
 125  FORMAT(/' Z-DISTIBUTION')
      NZ0=MIN0(NZ,12)
      MZ0=ZCN-NZ
      WRITE(LOGUN,130) (ELEMNT(MZ0+NROW),NROW=1,NZ0)
 130  FORMAT(12(4X,A4,3X))
      WRITE (LOGUN,135) (SIGZ(NZ+1-NROW),NROW=1,NZ0)
 135  FORMAT(12E11.3)
C--------------------------------------------------------
      WRITE (LOGUN,150)
 150  FORMAT(/'  CROSS SECTIONS (MB)')
      NA0=MIN0(NA,12)
      ZAV=0.
      AAV=0.
      ZDUM=0.
      ADUM=0.
      SUM=0.
      DO 170 NROW=1,NZ
      MZ=IZCN -NROW+1
      ZZ=MZ
      IDUM=IACN -NROW+2
      DO 155 NSTEP=1,NA
      XSECTE(NROW,NSTEP,NEXC)=XSECT(NROW,NSTEP)
      IND(NSTEP)=IDUM-NSTEP
      AA=IND(NSTEP)
      XS=XSECT(NROW,NSTEP)
      ZAV=ZAV+ZZ*XS
      AAV=AAV+AA*XS
      ZDUM=ZDUM+ZZ*ZZ*XS
      ADUM=ADUM+AA*AA*XS
      SUM=SUM+XS
 155  CONTINUE
      WRITE(LOGUN,160)(IND(NA+1-NSTEP),ELEMNT(MZ),NSTEP=1,NA0)
      WRITE(LOGUN,165)(XSECT(NROW,NA+1-NSTEP),NSTEP=1,NA0)
 160  FORMAT(12(I5,A4,2X))
 165  FORMAT(12E11.3)
 170  CONTINUE
C----------------------------------------------------------
      IF(SUM.LE.0.) GOTO 176
      ZAV=ZAV/SUM
      AAV=AAV/SUM
      ZDUM=ZDUM/SUM
      ADUM=ADUM/SUM
      DELTAZ=IZCN-ZAV
      DELTAA=IACN-AAV
      FWHMZ=0.
      FWHMA=0.
      DUMZ=ZDUM-ZAV*ZAV
      DUMA=ADUM-AAV*AAV
      IF(DUMZ.LT.0. .OR. DUMA.LT.0.) GOTO 300
      IF(DELTAZ.GT.0.01) FWHMZ=2.36*SQRT(ZDUM-ZAV*ZAV)
      IF(DELTAA.GT.0.01) FWHMA=2.36*SQRT(ADUM-AAV*AAV)
      WRITE(LOGUN,172) ZAV,DELTAZ,FWHMZ,AAV,DELTAA,FWHMA
 172  FORMAT(/' AVERAGE Z,DELTAZ,FWHMZ SO FAR', 3F8.2/
     1        ' AVERAGE A,DELTAA,FWHMA SO FAR', 3F8.2/)
C------------------------------------------------------------
      ROUND=DEL(1)-AINT(DEL(1))
      WRITE(LOGUN,311) SYMB(1)
      CALL OUT1(SPECGT,20,ROUND,1.D0,LOGUN)
      IF(NOE2.EQ.1) GOTO 304
      ROUND=DEL(1)-AINT(DEL(1))
      WRITE(LOGUN,310) SYMB(1)
      CALL OUT1(SPECQT,20,ROUND,1.D0,LOGUN)
 304  ROUND=DEL(2)-AINT(DEL(2))
      WRITE(LOGUN,310) SYMB(2)
      CALL OUT1(SPECNT,25,ROUND,1.D0,LOGUN)
      ROUND=DEL(3)-AINT(DEL(3))
      WRITE(LOGUN,310) SYMB(3)
      CALL OUT1(SPECPT,45,ROUND,1.D0,LOGUN)
      ROUND=DEL(4)-AINT(DEL(4))
      WRITE(LOGUN,310) SYMB(4)
      CALL OUT1(SPECAT,80,ROUND,1.D0,LOGUN)
      IF(NOF.EQ.1) GOTO 176
      ROUND=DEL(5)-AINT(DEL(5))
      WRITE(LOGUN,310) SYMB(5)
      CALL OUT1(SPECFT,20,ROUND,1.D0,LOGUN)
 310  FORMAT(' SPECTRUM ',A1)
 311  FORMAT(/' SPECTRUM ',A1)
C------------------------------------------------------------
 176  IF(IOPT.NE.0) GOTO 180
      WRITE(LOGUN,173)
 173  FORMAT(/' HALFLIFE FIRST STEP IN UNITS OF 10**(-22) SEC AS FUNTION
     1 OF J')
      WRITE(LOGUN,175)(FLAND(I),I=1,JUPYR)
 175  FORMAT(10E12.3)
C----------------------------------------------------------
 180  IF(NUMISO.EQ.0) GOTO 200
      WRITE(LOGUN,185)
 185  FORMAT(/' ISOMER CROSS SECTIONS')
      DO 190 ISO=1,NUMISO
      MZ=IZCN-ISOZ(ISO)+1
      IDUM=IACN-ISOZ(ISO)-ISOA(ISO)+2
      WRITE(LOGUN,187) IDUM,ELEMNT(MZ),SIGISO(ISO)
 187  FORMAT(I5,A4,2X,E12.4)
 190  CONTINUE
C----------------------------------------------------------
 200  IF(JFIS.NE.0) THEN
         DO 205 J=1,JDIM
         J1=JDIM+1-J
         IF(YFISS(J1).GT.0.) GOTO 206
 205     CONTINUE
 206     N1MX=J1
         IF(JFIS.EQ.5) THEN
            DO 208 J=1,N1MX
 208        YFISS(J)=1.00000-YFISS(J)
            WRITE(19,210) ELAB,N1MX,JFIS
 210        FORMAT('C 1-YFISS',4X,'ELAB=',F8.3,I6,I4/'H: X-1    Y')
         ELSE
            WRITE(19,211) ELAB,N1MX,JFIS
 211        FORMAT('C YFISS',4X,'ELAB=',F8.3,I6,I4/'H: X-1    Y')
         ENDIF
         WRITE(19,215) (J,YFISS(J),J=1,N1MX)
 215     FORMAT(4(I4,E12.4,';'))
      ENDIF
C----------------------------------------------------------
 250  IF(NEXC.GE.NENMX) CALL CASOUT(LOGUN,IDISC)
      RETURN
C----------------------------------------------------------
 300  WRITE(LOGUN,301) ZDUM,ADUM,ZAV,AAV ,SUM
 301  FORMAT(' ENGOUT NEG SQRT  ZDUM,ADUM,ZAV,AAV=',4E12.4,' SUM=',E12.4
     1      )
      GOTO 176
      END