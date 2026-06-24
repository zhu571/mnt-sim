      SUBROUTINE ALLOCH
C
C  Allokiert Ein- und Ausgabe - Datensaetze fuer HIVAP
C  (Einheit 21 wird in halloc.f allokiert)
C  Stand:  3.5.1994
C
      INTEGER*4 IST
      INTEGER*4 KST
      OPEN(5,FILE='hivapein.dat',STATUS='OLD',IOSTAT=KST)
      OPEN(6,FILE='hivaperg.dat',STATUS='UNKNOWN')
      OPEN(9,FILE='mlz.dat',STATUS='OLD')
      OPEN(10, ACCESS='SEQUENTIAL', STATUS='SCRATCH' )
      OPEN(11, ACCESS='SEQUENTIAL', STATUS='SCRATCH' )
      OPEN(12, ACCESS='SEQUENTIAL', STATUS='SCRATCH' )
      OPEN(13, ACCESS='SEQUENTIAL', STATUS='SCRATCH' )
      OPEN(14,FILE='dummi1.dat',STATUS='UNKNOWN')
      OPEN(15,FILE='dummi2.dat',STATUS='UNKNOWN')
      OPEN(16,FILE='dummi3.dat',STATUS='UNKNOWN')
      OPEN(19,FILE='dummi4.dat',STATUS='UNKNOWN')
      OPEN(22,FILE='dummi6.dat',STATUS='UNKNOWN')
      OPEN(23,FILE='hivapaus.dat',STATUS='UNKNOWN')
      OPEN(25,FILE='dummi7.dat',STATUS='UNKNOWN')
!---------- The following files are opened by N.Wang ------
	OPEN(1236,FILE='TEMP0.DAT')  
	OPEN(3,   FILE='TEMP1.DAT')  
	OPEN(101, FILE='INPUT.DAT')
!----------------------------------------------------------
      RETURN
      END