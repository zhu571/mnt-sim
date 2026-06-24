      SUBROUTINE CASOUT(LOGUN,IDISC)
C
C  Steuerung der Datenausgabe
C  Stand: 5. 5. 1994
C
      INCLUDE 'common.f'
      COMMON/DELFIS/EDELFIS,ADELFIS,INOF
      DIMENSION XSUM(30),EXCITE(NENMX),EXCITE1(NENMX)
      DIMENSION ZAV(NENMX),AAV(NENMX),SUMZN(NENMX),ZSTDEV(NENMX)
      DIMENSION ASTDEV(NENMX),SUMN(NENMX,27)
C ------------------------------------------------------------
      AMPAR1=0.
      LOGUN=6
      IF(NEXC.EQ.0) RETURN
      AP=MPROJ
      AT=MTGT
      WRITE(LOGUN,50)
 50   FORMAT(/' HIVAP  CASOUT')
      WRITE(LOGUN,55)(TITLE(I),I=1,18)
 55   FORMAT(/18A4/)
      IF(INOF.NE.1) GOTO 57
      WRITE(LOGUN,56) EDELFIS,ADELFIS
 56   FORMAT(/2X,'Rechnung mit Behinderung der Spaltung; Parameter:',
     11X,'E_Beh=',F5.1,2X,'DE_Beh=',F5.1)
 57   WRITE(LOGUN,58)
 58   FORMAT(/5X,'E_lab  A*MeV   E*/MeV',5X,'SI_fus',6X,'SI_VR ',4X,
     1       'SI_Spalt  L_krit',5X,'b_Spalt'/)
C
      DO 70 I=1,NEXC
      ZAV(I)=0.
      AAV(I)=0.
      SUMZN(I)=0.
      DO NROW=1,NZ
        SUMN(I,NROW)=0.
      ENDDO
      ZSTDEV(I)=0.
      ASTDEV(I)=0.
      DUM=ELABE(I)*FLOAT(MTGT )/ACN +  QC
      DUM0=ELABE1(I)*FLOAT(MTGT )/ACN +  QC
      DUM1=ELABE(I)/FLOAT(MPROJ)
      DUM11=ELABE1(I)/FLOAT(MPROJ)
      DUM2=0.
      IF(SIGFUE(I).GT.0.) DUM2=SIGFIE(I)*100./SIGFUE(I)
      WRITE(LOGUN,60)ELABE(I),DUM1,DUM,SIGFUE(I),SIGEVE(I),SIGFIE(I),
     1               FLGRAE(I),DUM2
      IF(ELABE1(I).GT.0.)WRITE(LOGUN,60) ELABE1(I),DUM11,DUM0
      EXCITE(I)=DUM
      IF(ELABE1(I).GT.0.) EXCITE1(I)=DUM0
 60   FORMAT(F10.3,F8.2,F8.1,3E12.4,F8.1,F12.2)
 70   CONTINUE
C---------------------------------------------
      IZCN=ZCN+0.01
      IACN=ACN+0.01
      WRITE (LOGUN,745)
745   Format(/'     Querschnitte / mbarn')
      NA0=MIN0(NA,12)
C
      DO 170 NROW=1,NZ
      MZ=IZCN -NROW+1
      IDUM=IACN -NROW+2
      ZZ=MZ-IZCN
C
      DO 155 NSTEP=1,NA
      IND(NSTEP)=IDUM-NSTEP
      AA=IND(NSTEP)-IACN
      DO IE =1,NEXC
        XS=XSECTE(NROW,NSTEP,IE)
        ZAV(IE)=ZAV(IE)+ZZ*XS
        AAV(IE)=AAV(IE)+AA*XS
        ZSTDEV(IE)=ZSTDEV(IE)+ZZ*ZZ*XS
        ASTDEV(IE)=ASTDEV(IE)+AA*AA*XS
        SUMZN(IE)=SUMZN(IE)+XS
        SUMN(IE,NROW)=SUMN(IE,NROW)+XS
      ENDDO                     ! IE
 155  CONTINUE                  ! NSTEP
C
      WRITE(LOGUN,'(/''  E*/MeV '',12(I5,A4,1X))')
     1            (IND(NA+1-NSTEP),ELEMNT(MZ),NSTEP=1,NA0)
      DO     IE=1,NEXC
         WRITE(LOGUN,'(F8.1,12E10.3)')
     1            EXCITE(IE),(XSECTE(NROW,NA+1-NSTEP,IE),NSTEP=1,NA0)
         IF(EXCITE1(IE).GT.0.) WRITE(LOGUN,'(F8.1)') EXCITE1(IE)
      ENDDO   ! IE
 170  CONTINUE
C
      IF(NZ.GT.2) THEN
        WRITE(LOGUN,'(/'' Z-Verteilung  '')')
        NZ0=MIN(NZ,12)
        MZ0=ZCN-NZ
        WRITE(LOGUN,'('' E*/MeV '',12(4X,A4,3X))')
     1            (ELEMNT(MZ0+NROW),NROW=1,NZ0)
        DO IE=1,NEXC
         IF(SUMZN(IE).GT.0.) THEN
            WRITE(LOGUN,'(F7.1,E11.3,11E10.3)') EXCITE(IE),
     1           (SUMN(IE,NZ+1-NROW),NROW=1,NZ0)
            IF(EXCITE1(IE).GT.0.) WRITE(LOGUN,'(F7.1)') EXCITE1(IE)
            ZAV(IE)=ZAV(IE)/SUMZN(IE)
            AAV(IE)=AAV(IE)/SUMZN(IE)
            ZSTDEV(IE)=ZSTDEV(IE)/SUMZN(IE)
            DUM=DMAX1(AMPAR1,ZSTDEV(IE)-ZAV(IE)**2)
            ZSTDEV(IE)=SQRT(DUM)
            ZAV(IE)=ZAV(IE)+ZCN
            ASTDEV(IE)=ASTDEV(IE)/SUMZN(IE)
            DUM=DMAX1(AMPAR1,ASTDEV(IE)-AAV(IE)**2)
            ASTDEV(IE)=SQRT(DUM)
            AAV(IE)=AAV(IE)+ACN
         ELSE
            ZAV(IE)=0.
            AAV(IE)=0.
            ZSTDEV(IE)=0.
            ASTDEV(IE)=0.
         ENDIF
        ENDDO             ! IE
C
        WRITE(LOGUN,'(/'' Z und A Mittelwerte    ''/
     1''   E*/MeV     MW.Z  St.Abw.Z      MW.A  St.Abw.A'') ')
        DO IE=1,NEXC
          WRITE(LOGUN,'(F8.1,4F10.2)')
     1       EXCITE(IE),ZAV(IE),ZSTDEV(IE),AAV(IE),ASTDEV(IE)
        ENDDO
      ENDIF   ! NZ>3
C ----------------------------------------------------------------
      IF(IDISC.EQ.0) GOTO 300
      IF(IDISC.EQ.2) GOTO 210
      DO 200 NROW=1,NZ
      MZ=IZCN-NROW+1
      DO 185 NSTEP=1,NA
      DO 171 I=1,NEXC
      IF(XSECTE(NROW,NSTEP,I).GT.0.) GOTO 174
 171  CONTINUE
      GOTO 185
 174  MA=IACN-NROW-NSTEP+2
      IF(MA.GT.99) WRITE(25,172) MA,ELEMNT(MZ),MA,ELEMNT(MZ)
      IF(MA.LE.99) WRITE(25,173) MA,ELEMNT(MZ),MA,ELEMNT(MZ)
 172  FORMAT('H:  X  Y  ''',I3,A4,''''/'C ',I3,A4)
 173  FORMAT('H:  X  Y  ''',I2,A4,''''/'C ',I2,A4)
      WRITE(25,175) (EXCITE(I),XSECTE(NROW,NSTEP,I),I=1,NEXC)
 175  FORMAT(4(F6.1,E11.3,';'))
 185  CONTINUE
 200  CONTINUE
C-----------------------------
      IF(IDISC.EQ.3) GOTO 210
      GOTO 230
 210  WRITE(25,212)
 212  FORMAT('H:   X Y  ''SIGFUS'''/'C EXCIT VS SIGFUS(MB)')
      WRITE(25,175) (EXCITE(I), SIGFUE(I),I=1,NEXC)
      WRITE(25,214)
 214  FORMAT('H:   X Y  ''SIGER'''/'C EXCIT VS SIGEVA(MB)')
      WRITE(25,175)(EXCITE(I),SIGEVE(I),I=1,NEXC)
      WRITE(25,217)
 217  FORMAT('H:   X Y ''SIGFISS'''/'C EXCIT VS SIGFISS(MB)')
      WRITE(25,175)(EXCITE(I),SIGFIE(I),I=1,NEXC)
C --------
      DO 226 NROW=1,NZ
      MZ=IZCN-NROW+1
      WRITE(25,220) ELEMNT(MZ)
 220  FORMAT('H: X Y  ''',A4,'''')
      DO 224 I=1,NEXC
         XSUM(I)=0.
         DO 222 NSTEP=1,NA
         XSUM(I)=XSUM(I)+XSECTE(NROW,NSTEP,I)
 222     CONTINUE
 224  CONTINUE
      WRITE(25,175) (EXCITE(I),XSUM(I),I=1,NEXC)
 226  CONTINUE
C-------------
 230  DUM=FLOAT(MTGT)/ACN
      WRITE(25 ,231) QC,DUM,NEXC
 231  FORMAT('C   QC=',F8.2,4X,'AT/ACN=',F9.4,4X,'NEXC=',I4)
C
      DUM0=-QC
      DUM1=1.
      WRITE(25,232) DUM0,DUM1
 232  FORMAT('C RESC  X0',F9.3,2X,'X1',F4.0,6X,' FROM 1',
     1       4X,'SHIFT TO ECM')
C
      DUM0=-QC*ACN/AT
      DUM1=ACN/AT
      WRITE(25,234) DUM0,DUM1
 234  FORMAT('C RESC  X0',F9.3,2X,'X1',F10.5,' FROM 1',
     1      4X,'SHIFT TO ELAB')
C
      DUM1=ACN/(AP*AT)
      DUM0=DUM1*(-QC)
      WRITE(25,236) DUM0,DUM1
 236  FORMAT('C RESC  X0',F9.3,2X,'X1',F10.6,' FROM 1',4X,'SHIFT TO MEV/
     1U')
 300  NEXC=0
      RETURN
      END