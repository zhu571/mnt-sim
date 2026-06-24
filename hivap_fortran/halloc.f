      SUBROUTINE HALLOC(DATNAM)
C
C  Allokiert Ein- und Ausgabe - Datensaetze fuer HIVAP
C  Stand:  3.5.1994
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*12 DATNAM
      OPEN(21,FILE=DATNAM,STATUS='UNKNOWN',IOSTAT=IST)
      RETURN
      END