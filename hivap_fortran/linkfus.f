	SUBROUTINE LINKFUS(MPROJ,IZP,MTGT,IZT,NZ,NA)
C     This code is to calculate fusion (capture) excitation function
C     based on the fusion barrier obtained by ETF2+Skyrme energy-density funtional
C     together with the parameterized barrier distribution. 
C
C     Ning Wang 2005-12-2     Reference: Nuclear Physics A 768 (2006) 80�C98
C                                        ArXiv: Nucl-th/0509069
C     
!	INCLUDE 'alloch.f'

	implicit real*8 (a-h,o-z)
      CHARACTER*2 TEXT,PROG,TARG
	INTEGER*4 NPROG, NTARG 
	dimension Text(110)
	dimension Mwt(110)  !reference system
	data mwt/1,4,7,9,11,12,14,16,19,21,23,25,27,29,31,33,36,40,
     &        40,41,45,48,51,52,55,56,59,59,64,66,70,73,75,79,80,84, 
     & 86,88,89,92,93,96,98,102,103,107,108,113,115,119,122,128,128,132,
     & 133,138,
     & 139,141,141,145,145,147,152,158,159,163,165,168,169,173,175,
     & 179,181,184,187,191,193,196,197,201,205,208,209,210,210,222,
     & 223,226,
     & 227,232,231,238,237,244,243,247,247,251,252,257,258,259,260,
     & 261,262,263,264,265,268,269/

      DATA TEXT/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     &'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ',
     &'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
     &'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     &'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',
     &'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re',
     &'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
     &'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md',
     &'No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds'/
 


!	OPEN(2,FILE='BarETF2.SkM')  !Fusion barrier with SkM* Skyrme force

!	WRITE(*,*)"INPUT REACTION: A1, Z1, A2, Z2"
!	READ(*,*)A1,Z1,A2,Z2
	A1=FLOAT(MPROJ)
	Z1=FLOAT(IZP)
	A2=FLOAT(MTGT)
	Z2=FLOAT(IZT)

	PROG=TEXT(Z1)
      TARG=TEXT(Z2)
	Nprog=int(A1)
	Ntarg=int(A2)
	WRITE(*,*)Nprog,PROG," + ",Ntarg,TARG, "   Fusion reaction"

	OPEN(212,FILE='SIGMA.DAT')
!	WRITE(212,*)"===================================================="
!	WRITE(212,*)Nprog,PROG," + ",Ntarg,TARG, "   Fusion reaction"
!	WRITE(212,*)"===================================================="
!	WRITE(212,*)"      ECM,            SIGFUS,
!     &        SIGEvR,         SIGFIS"

	OPEN(313,FILE='SIGXPN.DAT')
!	WRITE(313,*)"===================================================="
!	WRITE(313,*)Nprog,PROG," + ",Ntarg,TARG, "   Fusion reaction"
!	WRITE(313,*)"===================================================="
!	WRITE(313,*)"    Z,    N,      ECM,          SIGFUS, 
!     &       SIGEvR,       SIGFIS"

!	do while(.not.EOF(2))
!	read(2,*)A1i,Z1i,A2i,Z2i,R0i,B0i,hw0i,Rsi,Bsi
!	if(  A1i.eq.A1.and.Z1i.eq.Z1.and.A2i.eq.A2.and.Z2i.eq.Z2
!     & .or.A1i.eq.A2.and.Z1i.eq.Z2.and.A2i.eq.A1.and.Z2i.eq.Z1)then
!	R0=R0i
!	B0=B0i
!	hw0=hw0i
!	Bs=Bsi
!	write(*,'(1x,"R0, B0, hw=",3f12.6)')B0,R0,hw0
!	if(Bsi.gt.B0i)write(*,*)"ATTENTION: fusion pocket disappears!"
!	exit
!	endif
!	enddo
!	rewind(2)
!	if(R0.eq.0.or.B0.eq.0)then
!	Write(*,*)"Sorry, the data for this system are not available!"
!	write(*,*)"Input B0, R0, hw0"
!	read(*,*)B0, R0, hw0
!	Stop
!	endif

	call Barrier(A1,Z1,A2,Z2,B0,R0,hw0) !Chin.Phys.Lett.24 (2007) 905
	write(*,'(1x,"FUSION BARRIER: B0, R0, hw0=",3f12.4)')B0,R0,hw0

101	EcmL=B0*(1-0.2)
	EcmH=B0*(1+0.3)

	Q=Qvalue(A1,Z1,A2,Z2)
	write(*,'(1X,"Q=",F12.4)')Q


	f=0.926   !Nucl. Phys. A 768 (2006) 80
	gamma=1.  


	N1=int(A1-Z1)
	N2=int(A2-Z2)
	fs1=0   !flag for neutron shell-closure nuclei
	fs2=0
	if(N1.eq.20.or.N1.eq.28.or.N1.eq.50.or.N1.eq.82.or.N1.eq.126)fs1=1
	if(N2.eq.20.or.N2.eq.28.or.N2.eq.50.or.N2.eq.82.or.N2.eq.126)fs2=1
	if(A1.gt.(Mwt(Z1)+3).or.fs1.eq.1.or.
     &   A2.gt.(Mwt(Z2)+3).or.fs2.eq.1)then

	A10=float(Mwt(Z1))
	A20=float(Mwt(Z2))
	write(*,'(1x,"The reference system is:",i8,f8.0,i8,f8.0)')
     &         Mwt(Z1),Z1,Mwt(Z2),Z2

	Q0=Qvalue(A10,Z1,A20,Z2)
	write(*,'(1X,"Q0=",F12.4)')Q0

	C0=0.5
	if(Q-Q0.gt.0)C0=0.1
	gamma=1-C0*(Q-Q0)+0.5*(fs1+fs2)
	if(gamma.lt.0.5)gamma=0.5
	write(*,'(1X,"Gamma=",F12.2)')gamma
	endif


	call mwkb(R0,B0,hw0,EcmL,EcmH,f,gamma,NZ,NA)

	stop
	end
!===================================================

	FUNCTION Qvalue(A1,Z1,A2,Z2)
      implicit double precision(a-h,o-z)
	DIMENSION POP(100)
	MASSES=9
	IZP=INT(Z1)
	IZT=INT(Z2)
	IZCN=INT(Z1+Z2)
	MPROJ=INT(A1)
	MTGT=INT(A2)
	IACN=INT(A1+A2)
	REWIND(MASSES)
C     DO WHILE (.NOT.EOF(MASSES))
	DO
140   READ(MASSES,*,IOSTAT=IOS) IZZ,IADN,IAUP
      IF(IOS.NE.0) EXIT
      NPOW=IAUP-IADN+1
      READ(MASSES,*,IOSTAT=IOS)(POP(I),I=1,NPOW)
      IF(IOS.NE.0) EXIT
      IF(IZZ.EQ.IZP) EXCP=POP(MPROJ-IADN+1)
      IF(IZZ.EQ.IZT) EXCT=POP(MTGT-IADN+1)
      IF(IZZ.EQ.IZCN) EXCCN=POP(IACN-IADN+1)
      QVALUE=EXCP+EXCT-EXCCN
	IF(IZZ.GT.IZCN)EXIT
	ENDDO
	END

!========================================================
	subroutine mwkb(R0,B0,hw,EcmL,EcmH,f,gamma,NZ,NA)
	implicit real*8 (a-h,o-z)
!     to calculate the fusion excitation function by multi-dimentional WKB method
!     wangning
!	implicit real*8 (a-h,o-z)

	PI=3.1415926
	Bc=B0*F
	width1=(B0-Bc)/4.
	width2=(B0-Bc)/2.

	B1=Bc+width1
	B2=Bc+width2

	F1=2*Sqrt(PI)*width1/sqrt(gamma)
	F2=2*Sqrt(PI)*width2

	SIGZNS=0.
	RLOSTFZNS=0.
	ECM0=0.
	rewind(1236)
10	READ(1236,'(3E20.12,2I8,2E20.12)',ERR=100,END=1000)
     +ELAB,ECM,SIGFUS,NROW,NSTEP,SIGZN,RLOSTFZN

	if(Ecm.ne.Ecm0)then
	sig1=0
	sig2=0
	sig1up=0
	sig2up=0
	sig1dn=0
	sig2dn=0
	DB=0.1
	do B=B2-6*width2,B2+6*width2,DB
  !according to nucl-th/0302025 and C.Y.Wong, Phys.Rev.Lett. 31,(1973)766
	if(Ecm/B.gt.1.3)then
	sig=R0*R0*hw/2./Ecm*(2.*PI/hw*(Ecm-B))*10 
	else
	sig=R0*R0*hw/2./Ecm*log(1+exp(2.*PI/hw*(Ecm-B)))*10 
	endif
	Wb1=Exp(-gamma*(B-B1)**2/(2*width1)**2)/F1
	Wb2=Exp(-(B-B2)**2/(2*width2)**2)/F2

	sig1=sig1+sig*Wb1*DB
	sig2=sig2+sig*Wb2*DB


	Wb1up=Exp(-0.5*(B-B1)**2/(2*width1)**2)
     &     /(2*Sqrt(PI)*width1/sqrt(0.5))
	Wb1dn=Exp(-20.*(B-B1)**2/(2*width1)**2)
     &     /(2*Sqrt(PI)*width1/sqrt(20.))

	sig1up=sig1up+sig*Wb1up*DB
	sig1dn=sig1dn+sig*Wb1dn*DB


	enddo

	avsig=(sig1+sig2)/2
	sigma=min(sig1,avsig)

	avsigup=(sig1up+sig2)/2
	avsigdn=(sig1dn+sig2)/2
	sigup=min(sig1up,avsigup)
	sigdn=min(sig1dn,avsigdn)


	endif
	
	write(3,'(1x,3f15.9,2i6,2f15.9)')
     +Ecm,sigma,SIGFUS,NROW,NSTEP,
     +SIGZN/SIGFUS*sigma,RLOSTFZN/SIGFUS*sigma

	IF(ECM.NE.ECM0.AND.ECM0.GT.0)THEN
	WRITE(212,'(6F16.9)')ECM0,sigma0,SIGZNS,RLOSTFZNS
     &    ,max(sigup0,sigma0*(1+0.18)),min(sigdn0,sigma0*(1-0.18))  !error bar
	SIGZNS=0.
	RLOSTFZNS=0.
	ENDIF

	SIGZNS=SIGZNS+SIGZN/SIGFUS*sigma
	RLOSTFZNS=RLOSTFZNS+RLOSTFZN/SIGFUS*sigma

	ECM0=ECM
	sigma0=sigma
	sigup0=sigup
	sigdn0=sigdn
	GOTO 10

100	PRINT*,"ERROR:",ELAB,ECM,SIGFUS,NROW,NSTEP,SIGZN,RLOSTFZN

1000	PRINT*,"THE PROGRAM STOPS! GOOD LUCK!"


	DO NZTPS=1,NZ
	DO NSTPS=1,NA

	REWIND(3)

	DO

	READ(3,'(1x,3f15.9,2i6,2f15.9)',IOSTAT=IOS)
     +Ecm,sigma,SIGFUS,NROW,NSTEP,
     +SIGZN,RLOSTFZN
	IF(IOS.NE.0) EXIT

	IF(NSTPS.EQ.NSTEP.AND.NZTPS.EQ.NROW)THEN
	WRITE(313,'(2I6,4F15.9)')
     +NZTPS-1,NSTPS-1,ECM,SIGMA,SIGZN,RLOSTFZN
	ENDIF

	ENDDO

	ENDDO
	ENDDO

!---------- mean barrier height and most probable barrier height, 2008-08-21, N. Wang ----------
	avbar=0
	Dmax=0
	do B=B2-6*width2,B2+6*width2,DB
	Wb1=Exp(-gamma*(B-B1)**2/(2*width1)**2)/F1
	Wb2=Exp(-(B-B2)**2/(2*width2)**2)/F2
	Wbav=(Wb1+Wb2)/2.
	Wb=Wbav
	if(Wb1.le.Wb2.and.B.lt.Bc)Wb=Wb1
	avbar=avbar+B*wb*DB
	avD=avD+wb*DB
		if(wb.gt.Dmax)then
		Dmax=wb
		Bpeak=B
		endif
	enddo
	write(*,*)
 	write(*,'(1x,"Mean barrier height (in MeV) : ",f12.4)') avbar/avD
	write(*,'(1x,"Most probable barrier height : ",f12.4)') Bpeak 
!-----------------------------------------------------------------------------------------------
	end		

	subroutine Barrier(A1,Z1,A2,Z2,B0,R0,hw0)
	implicit real*8 (a-h,o-z)
	
	common/mydata/XX(300),YY(300),Npoints
	dimension BAR(300),DIST(300)

	I=0
	do DIS=16,7,-0.25
	VB=VMWS(A1,Z1,A2,Z2,DIS)
	I=I+1
	BAR(I)=VB
	DIST(I)=DIS
	enddo

!=========== FIND R0, B0, hw0 and Bs  ============  
	KEY=0
	DO I=2,300
	IF(BAR(I).LT.BAR(I-1).AND.KEY.EQ.0)THEN
	KEY=1
	B0=BAR(I-1)
	R0=DIST(I-1)
	ENDIF
	IF(BAR(I).GT.BAR(I-1).AND.KEY.EQ.1)THEN
	Bs=BAR(I-1)
	Rs=DIST(I-1)
	KEY=2
	ENDIF
	IF(KEY.EQ.2)EXIT
	ENDDO

!=== To obtain the curvature of barrier through data fitting using an inverted parabola ===
	NFit=0
	DO I=1,300
	IF(DIST(I).GE.R0-1.25.AND.DIST(I).LE.R0+1.25
     &                     .AND.DIST(I).GE.Rs)THEN
	Nfit=Nfit+1
	XX(Nfit)=DIST(I)
	YY(Nfit)=BAR(I)
	ENDIF
	ENDDO
	Npoints=Nfit
	CALL FITTING(hw0)

	return
	end


      function VMWS(A1,Z1,A2,Z2,R) 
!     Modified Woods-Saxon potential fit to the potential barrier 
!     obtained by Skyrme energy-density functional plus ETF2 [Nucl.phys.A768(2006)80]
!     Ning Wang 2007-5-22  Chin.Phys.Lett.24 (2007) 905
	implicit real*8 (a-h,o-z)

	r0=1.27
	c=-1.37
	as=-44.16
	a=0.75
	cs=-0.40


	d0=r0*(A1**0.3333+A2**0.3333)+c
	V0=as*(1+cs*((A1-2*Z1)/A1+(A2-2*Z2)/A2))
     &   *(A1**0.3333*A2**0.3333)/(A1**0.3333+A2**0.3333)
	VN=V0/(1+Exp((R-d0)/a))

	VMWS=Vn+Z1*Z2*1.44/R
	return
	end


      subroutine FITTING(hw0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	implicit real*8 (a-h,o-z)


	EXTERNAL fun
	common/mydata/XX(300),YY(300),Npoints

      dimension Xmin(10),Xo(10),Eo(10,10),X00(10),E00(10,10)

      DATA((E00(I,J),I=1,10),J=1,10)/1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     1	                         0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,
     2                             0.,0.,1.,0.,0.,0.,0.,0.,0.,0.,
     3                             0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,
     4	                         0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,
     5	                         0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,
     6	                         0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,
     7	                         0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,
     8                             0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,
     9  	                         0.,0.,0.,0.,0.,0.,0.,0.,0.,1./


	x1=XX(1)
	y1=YY(1)
	x2=XX((1+Npoints)/2)
	y2=YY((1+Npoints)/2)
	x3=XX(Npoints)
	y3=YY(Npoints)

	factor=(x1-x2)*(x1-x3)*(x2-x3)
	a0=x1*x3*(x3-x1)*y2+x2*x2*(x3*y1-x1*y3)+x2*(x1*x1*y3-x3*x3*y1)
	b0=x3*x3*(y1-y2)+x1*x1*(y2-y3)+x2*x2*(y3-y1)
	c0=x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2)

	X00=0
	X00(1)=a0/factor
	X00(2)=b0/factor
	X00(3)=c0/factor

!the parameter of DSC method
      n=3                                                  !!!!!!!!!!!!!!n=?
	h=0.01
	err=1.e-9
	Xo=X00
	Eo=E00

	   call DSC(n,h,err,Xo,Eo,Xmin)

!	   fvalu=fun(Xmin)

	a=xmin(1)
	b=xmin(2)
	c=xmin(3)

!	print*,"Parabola Fitting... V=a+b*R+c*R^2"
!	print*,"a,b,c=",a,b,c
!	PRINT*,"chi^2=",fvalu

	hw0=-c      
	return

	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!the method of DSC  
C ====  WRITTEN BY MIN LIU  ===
 
      SUBROUTINE DSC(n,h,err,Xo,Eo,Xmin)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	dimension X0(10),X1(10),Z(10),d(10),fmop(10),pp(10),Xmin(10)
	dimension Xo(10)
      dimension E(10,10),Eo(10,10),P(10,10),Q(10,10)
      E=Eo
	X0=Xo
10    Z=X0
	K=1
	h0=h
20	do J=1,n
	   pp(J)=E(K,J)
	enddo

      call search(t,valu,h0,Z,pp)

      d(k)=t
	Z=Z+t*pp

	if (K.lt.n) then
	   K=K+1
	   goto 20
	endif
      pp=Z-X0
	call search(t,valu,h0,Z,pp)
	X1=Z+t*(Z-X0)
      I=1
	do J=1,n
	   p(n,J)=d(n)*E(n,J)
	enddo
30	if (I.lt.n)then
	   do J=1,n
	      L=n-I
	      M=n-I+1
	      p(L,J)=P(M,J)+d(L)*E(L,J)
	   enddo
	   I=I+1
         goto 30
	endif
      fmop=0
	do I=1,n
         do K=I,n
            fmop(I)=fmop(I)+d(K)**2
	   enddo
	   fmop(I)=sqrt(fmop(I))
      enddo
	if(fmop(1).le.h) then
	   if(fmop(1).lt.err)then
	      goto 100
	   else
	      X0=X1
	      h=0.1*h
	      goto 10
	   endif
      else
          I=n
40        if(fmop(I).lt.1.00000e-15)then
	      goto 50
	    else
	      L=I-1
	      do J=1,n
	         E(I,J)=(d(L)*P(I,J)-fmop(I)**2*E(L,J))/fmop(L)/fmop(I)
	      enddo
          endif
50        if(I.eq.2)then
             do J=1,n
	          E(1,J)=P(1,J)/fmop(1)
	       enddo
	       X0=Z
	       Z=X1
	       d(1)=t
	       K=2
	       goto  20
	    else
	       I=I-1
	       goto 40
          endif
	endif

100   Xmin=X1
      return
	end
	       

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!one dimension search for the minimal value
      Subroutine search(t,valu,h0,X,pp)
      implicit double precision(a-h,o-z)
	dimension X(10),X1(10),X2(10),X3(10),X4(10),pp(10)
	err=1.0e-5
	errl=0.1
	t1=0
	h=h0
	X1=X+t1*pp
	f1=fun(X1)
	t2=h
	X2=X+t2*pp
	f2=fun(X2)
	if(f2.gt.f1)then
	   h=-h
	   t3=t1
	   f3=f1
	else
	   goto 2000
	endif

1000	   t1=t2
         f1=f2
	   t2=t3
	   f2=f3
2000	t3=t2+h
	X3=X+t3*pp
	f3=fun(X3)
	if(f2.gt.f3)then                                         
	   h=h+h
	   goto 1000
      else                                                      
3000	   c1=(f3-f1)/(t3-t1)
	   c2=((f2-f1)/(t2-t1)-c1)/(t2-t3)
	   if(c2.eq.0.)then                                       
	     tmin=t2
	     fmin=f2
	     goto 4000
	   else                                                    
	     t4=0.5*(t1+t3-c1/c2)
	     sign1=(t4-t1)*(t3-t4)
	     if(sign1.gt.0.)then                                   
	        X4=X+t4*pp
	        f4=fun(X4)
	        if(abs(f2).lt.errl)then                            
	           AA=1.0
	        else                                               
	           AA=abs(f2)
	        endif                                              
	        if(abs(f2-f4)/AA.lt.err)then                       
	           if(f4.lt.f2)then                                
	              tmin=t4
	              fmin=f4
	              goto 4000
                 else                                            
	              tmin=t2
	              fmin=f2
	              goto 4000
	           endif                                           
	         else                                              
	            if((t4-t2)*h.gt.0.)then                        
	               if(f2.gt.f4)then                            
	                  t1=t2
	                  f1=f2
	                  t2=t4
	                  f2=f4
	               else                                        
	                  t3=t4
	                  f3=f4
	               endif                                       
	            else                                           
	               if(f2.gt.f4)then                            
	                   t3=t2
	                   f3=f2
	                   t2=t4
	                   f2=f4
	               else                                        
	                   t1=t4
	                   f1=f4
	               endif                                       
	            endif                                          
	            goto 3000
	        endif                                              
	     else                                                                                                 
	        tmin=t2
	        fmin=f2
	        goto 4000
	     endif                                                 
         endif                                                   
	endif                                                      
4000  t=tmin
      valu=fmin
	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!calculate the energy
      function fun(x)
      implicit double precision(a-h,o-z)
	common/mydata/XX(300),YY(300),Npoints

	dimension x(10)

!	f=a+bx+cx^2	
	a=x(1)
	b=x(2)
	c=x(3)
	sum=0
	do i=1,Npoints
	f=a+b*XX(i)+c*XX(i)*XX(i)
	sum=sum+(f-YY(i))*(f-YY(i))
	enddo
      fun=sum

!	print*,a,b,c,fun

 	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



