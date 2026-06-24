      SUBROUTINE DROPB(Z,A,INDEX,DPBAR,Y,LDM,IPRNT)
C
C  Droplet-Modell (frei nach Myers)
C  dient zur Berechnung des Oberflaechenverhaeltnisses von Sattelpunkt-
C  und Gleichgewichtskonfiguration, was a_f/a_n beeinflusst;
C  ruft DROPL auf.
C  Stand: 3.5.1994
C
C       DROPLET BARRIER
C    /* INDEX = 0 NUCLEUS TOO LIGHT FOR Y-FAMILY                   */
C    /*       = 1 SOLN WITH Y > .44                                */
C    /*       = 2 SOLN WITH Y = .44                                */
C    /*       = 3 SOLN WITH Y < .44                                */
C    /*       = 4 SOLN FOR HEAVY NUCLEI USING FULL Y_FAMILY SEARCH */
C    /*       = 5 FISSION UNSTABLE                                 */
C    /*       = 6 SOMETHING WRONG                                  */
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(36),BCA(36),BSA(36),BKA(36),BRA(36),BVA(36),BWA(36)
      DIMENSION YA(3 ), AA(3 )
      COMMON/BEES/BS,BC,BK,BR,BV,BW
      DATA
     1  XA
     2   /.30, .32, .34, .36, .38,
     3    .40, .42, .44, .46, .48,
     4    .50, .52, .54, .56, .58,
     5    .60, .62, .64, .66, .68,
     6    .70, .72, .74, .76, .78,
     7    .80, .82, .84, .86, .88,
     8    .90, .92, .94, .96, .98,
     9   1.00/
       DATA
     1  BCA
     2   /.82500, .82332, .82173, .82023, .81881,
     3    .81747, .81621, .81503, .81393, .81293,
     4    .81206, .81133, .81081, .81056, .81073,
     5    .81155, .81347, .81749, .82584, .84250,
     6    .86664, .89017, .91008, .92667, .94060,
     7    .95238, .96239, .97088, .97807, .98409,
     8    .98905, .99303, .99609, .99827, .99957,
     9   1.00000/
      DATA
     1  BSA
     2  /1.27314, 1.27418, 1.27522, 1.27627, 1.27732,
     3   1.27837, 1.27941, 1.28042, 1.28141, 1.28235,
     4   1.28320, 1.28394, 1.28450, 1.28477, 1.28458,
     5   1.28352, 1.28126, 1.27619, 1.26532, 1.24296,
     6   1.20963, 1.17623, 1.14717, 1.12229, 1.10085,
     7   1.08224, 1.06604, 1.05195, 1.03974, 1.02927,
     8   1.02044, 1.01319, 1.00750, 1.00338, 1.00086,
     9   1.00000/
      DATA
     1BKA
     2  /1.5968, 1.5917, 1.5867, 1.5817, 1.5768,
     3   1.5718, 1.5668, 1.5617, 1.5564, 1.5509,
     4   1.5452, 1.5389, 1.5321, 1.5243, 1.5151,
     5   1.5036, 1.4886, 1.4668, 1.4216, 1.3727,
     6   1.3014, 1.2418, 1.1951, 1.1576, 1.1267,
     7   1.1009, 1.0790, 1.0609, 1.0458, 1.0328,
     8   1.0224, 1.0141, 1.0079, 1.0035, 1.0009,
     9   1.0000/
      DATA
     1BRA
     2  / .4548,  .4538,  .4527,  .4516,  .4506,
     3    .4495,  .4486,  .4476,  .4468,  .4462,
     4    .4458,  .4458,  .4463,  .4477,  .4506,
     5    .4560,  .4661,  .4860,  .5295,  .6247,
     6    .7578,  .8660,  .9386,  .9845, 1.0119,
     7   1.0264, 1.0323, 1.0326, 1.0292, 1.0241,
     8   1.0181, 1.0123, 1.0071, 1.0032, 1.0007,
     9   1.0000/
      DATA
     1BVA
     2  / .7832,  .7810,  .7787,  .7765,  .7741,
     3    .7718,  .7694,  .7669,  .7643,  .7617,
     4    .7589,  .7561,  .7532,  .7502,  .7472,
     5    .7444,  .7422,  .7419,  .7473,  .7682,
     6    .8071,  .8462,  .8786,  .9048,  .9257,
     7    .9426,  .9564,  .9673,  .9761,  .9832,
     8    .9887,  .9930,  .9961,  .9983,  .9996,
     9   1.0000/
      DATA
     1BWA
     2  / .5782,  .5732,  .5685,  .5641,  .5601,
     3    .5564,  .5531,  .5500,  .5474,  .5451,
     4    .5434,  .5422,  .5410,  .5425,  .5447,
     5    .5493,  .5583,  .5752,  .6096,  .6781,
     6    .7673,  .8392,  .8899,  .9255,  .9503,
     7    .9677,  .9797,  .9877,  .9930,  .9963,
     8    .9983,  .9993,  .9998, 1.0000, 1.0000,
     9   1.0000/
C---------------------------------------------------------
      DPBAR = 3.0E+30
      Y = 3.0E+30

C        /* CHECK TO SEE IF THE SPHERE IS STABLE */
      BS = BSA(35)
      BC = BCA(35)
      BK = BKA(35)
      BR = BRA(35)
      BV = BVA(35)
      BW = BWA(35)
      IDUM=0
      CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IDUM )
      CALL DROPL (Z,A,1,DROP,DEL,EPS,LDM,IPRNT)
      IF((F2-DROP).GT.0.)  GO TO      2
      INDEX = 5
C     INDEX = 5                       /* SPHERE IS UNSTABLE     */
      RETURN
C        /* IF ICODE = 0 LOOK FOR MAX., IF ICODE = 1 LOOK FOR MIN. */
 2    ICODE = 0
C        /* FOR HEAVY NUCLEI USE ALL BEES IPERT = 0             */
C        /* WHEN IPERT = 1 A PERTURBATION APPROACH IS USED      */
C        /* IF Z > 68 A NON-PERT CALCULATION IS ATTEMPTED       */
      IPERT = 0
      IF(Z.GT.68)    GO TO      31
 30   IPERT = 1
      BS = BSA(35)
      BC = BCA(35)
      BK = 1
      BR = 1
      BV = 1
      BW = 1
      CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IDUM )
C
 31   DO 50 I=3,36
         J = 37 - I
         F1 = F2
         BS = BSA(J)
         BC = BCA(J)
         IF(IPERT.EQ.1)    GO TO      3
         BK = BKA(J)
         BR = BRA(J)
         BV = BVA(J)
         BW = BWA(J)
 3       CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IDUM )
         IF(ICODE.EQ.1)    GO TO      11
C        /* LOOKING FOR A MAXIMUM */
         IF(I.EQ.25 .AND. IPERT.EQ.0) GOTO     30
         IF((F2-F1).GT.0.)   GO TO      10
         IF((I-24).LT.0)    GO TO      12
         IF((I-24).EQ.0)    GO TO      13
         IF((I-24).GT.0)    GO TO      14
C        /* GOOD MAXIMUM */
 12      INDEX = 3
         IF(IPERT.EQ.0)    INDEX = 4
         GO TO      15
C        /* CONTINUE SEARCH LOOKING FOR MINIMUM */
 13      ICODE = 1
         GO TO      10
C        /* THIS SHOULD NOT HAPPEN */
 14      INDEX = 6
         RETURN
C        /* LOOKING FOR MINIMUM */
 11      IF((F2-F1).LT.0.)   GO TO      10
         IF(I.LT.14)    GO TO      14
         INDEX = 1
         GO TO      15
 10   CONTINUE
 50   CONTINUE
 17   IF(ICODE.EQ.0)    GO TO      16
C     /* TOO LIGHT, NO MINIMUM */
      INDEX = 0
      RETURN
C        /* NO MAX FOUND SO IT MUS BE AT Y = .44 */
 16   INDEX = 2
      J = 13
 15   J = J + 1
      BS = BSA(J)
      BC = BCA(J)
      BR = BRA(J)
      BK = BKA(J)
      BV = BVA(J)
      BW = BWA(J)
      RI2=(A-2.*Z)/A
      RI2=RI2*RI2
      XLDM= .7053 *Z*Z/(2.*A*17.9439*(1.-1.78  *RI2))
      IF(IPRNT.EQ.1)WRITE(6,60) XA(J),XLDM
 60   FORMAT(/' SADDLEPOINT X  =',F10.3,4X,'XLDM=',F10.3)
      CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IPRNT)
      DPBAR = F2 - DROP
      Y = 1.E0 - XA(J)
      IF(IPRNT.EQ.1)WRITE(6,70)XA(J),DPBAR
      YA(2)=DPBAR
      IREP=0
      JJ=J-1
      IF(JJ.LT.1)RETURN
 65   BS = BSA(JJ)
      BC = BCA(JJ)
      BR = BRA(JJ)
      BK = BKA(JJ)
      BV = BVA(JJ)
      BW = BWA(JJ)
      CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IDUM )
      DPBAR=F2-DROP
      IF(IREP.EQ.0)YA(1)=DPBAR
      IF(IPRNT.EQ.1)WRITE(6,70)XA(JJ),DPBAR
 70   FORMAT(' X=',F10.3,4X,'BARRIER=',F10.2)
      IREP=IREP+1
      IF(IREP.GT.1) GOTO 80
      JJ=J+1
      IF(JJ.GT.36)RETURN
      GOTO 65
 80   YA(3)=DPBAR
      J1=J-1
      CALL PARABO(XA(J1 ),YA,AA,X0,BAR)
      IF(IPRNT.EQ.1) WRITE(6,90)X0,BAR
 90   FORMAT(' PARABOLA X=',F10.3,4X,'BARRIER=',F10.2)
      J2=J+1
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BSA(J1),BSA(J),BSA(J2),BS)
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BCA(J1),BCA(J),BCA(J2),BC)
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BRA(J1),BRA(J),BRA(J2),BR)
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BKA(J1),BKA(J),BKA(J2),BK)
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BVA(J1),BVA(J),BVA(J2),BV)
      CALL IP2(XA(J1),XA(J),XA(J2),X0,BWA(J1),BWA(J),BWA(J2),BW)
      CALL DROPL (Z,A,2,F2,DEL,EPS,LDM,IPRNT)
      DPBAR = F2 - DROP
      Y = 1.E0 - X0
      RETURN
C
      END