program main

!&&& part 1

! 1.0 include library

implicit none


! 1.1 control area &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! 1.1.1 indicators

integer, parameter :: isys=1		! =1:	Windows
									! =2:	Mac
									
integer, parameter :: ict=1			! =1:	Ct+Kt+1=(1-delta)Kt+Kt^alpha	(consumption with CD production)
									! =2:	Ct+Kt+1=Kt	(cake eating)
					
integer, parameter :: iut=1			! =1:	Ut=log(Ct)
									! =2:	Ut=Ct^(1-gamma)/(1-gamma)
									
integer, parameter :: ialgo=2		! =1:	sequential search
									! =2:	binary search(concavity & monotonicity)
									! =3:	binary + gss & linear interpolation
									! =4:	binary + gss & cubic spline

integer, parameter :: irandom=0		! =0:	deterministic value fn iteraion
									! =1:	AR(1) to Markov
									! =2:	simple stochastic (given zspace and p)
									
integer, parameter :: ishow=0		! =1:	print stepwise result
									! =0:	do not show result


! 1.1.2 parameters


real,parameter :: beta=0.99
real,parameter :: alpha=0.3
real,parameter :: delta=0.1
real,parameter :: gamma=0.5
real,parameter :: theta=0.0
integer,parameter :: maxage=3
integer,parameter :: retage=3

real,parameter :: kmax=14.5
real,parameter :: kmin=0.0
integer,parameter :: kgrid=151

real :: gradkm(9)=(/ 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 /)
real,parameter :: tol=0.00001
real :: tau=theta/(2.0+theta)	!tau=0, so pen=0 too

real :: L=real(retage-1)/real(maxage)
real,parameter :: kinitmax=9.9
real,parameter :: kinitmin=0.1
integer,parameter :: kinitgrid=99

integer,parameter :: maxiter=10000

! declare variabls

! basic
real(8) v(kgrid,maxage),K,w,r,tK,d1(kgrid,maxage)
integer d(kgrid,maxage)

! feature
real(8) pen

! parameter experiment
integer gm,initm
real(8) kinitstep,kinit,kinitspace(kinitgrid)

! indexing spaces
real(8) kspace(kgrid),fktr(kgrid),fktw(kgrid),fkt(kgrid),kstep
integer kmaxindr(kgrid),kmaxindw(kgrid)

! loop
real(8) kdiff,gradk
integer i,j,iter,t

! aggregation
real(8) kcross(maxage),ktp1(maxage),sum,temp01
integer temp,ktp(maxage)

! sequential method
real(8) vtemp,vmax

! binary search
real(8) v3(3)
integer js,jmin,jmax,jl,ju

! gss
real(8) vx

! cspline
real(8) sig,qn,un,uu(kgrid,maxage),p,vsod(kgrid,maxage)

! else
integer dummy







! open file

if (isys==1) open(unit=7,file='C:\Users\NREM\Desktop\Dropbox\cversion\ckv')
if (isys==2) open(unit=7,file='ckv')

! kspace

kstep=(kmax-kmin)/real(kgrid-1)
kspace=(/ ( kmin+real(i-1)*kstep, i=1,kgrid ) /)


! print *, kspace
! print *, kinitspace
! pause


! 88 format (4(I5,2X))

! parameter loop begin &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

kinitstep=(kinitmax-kinitmin)/real(kinitgrid-1)
kinitspace=(/ ( kinitmin+real(i-1)*kinitstep, i=1,kinitgrid ) /)

do initm=29,29
	kinit=kinitspace(initm)

do gm=5,5
	gradk=gradkm(gm)





!&&& part 2
! initialize for loop

	iter=0
	kdiff=10.0
	K=kinit
	
	if (ialgo==4) call cspline
	
	
! K loop begin &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&& part 3 value fn loop

do while((kdiff>tol).and.(iter<=maxiter))
	sum=0.0
	iter=iter+1
	w=(1.0-alpha)*(K**alpha)*(L**(-alpha))
	r=alpha*(K**(alpha-1.0))*(L**(1.0-alpha))-delta
	
	pen=theta*(1.0-tau)*w
	
	! upper limit for Kt+1 known after rt, wt and pen are known
	! Ct+Kt+1=(1+rt)Kt+(1-tau)wt 		workers upper limit fktw
	! Ct+Kt+1=(1+rt)Kt+pen				retirees upper limit fktr
	! all workers upper limit for Kt+1 are same
	! all retirees upper limit for Kt+1 are same
	
	do i=1,kgrid
		fktr(i)=(1.0+r)*kspace(i)+pen
		fktw(i)=(1.0+r)*kspace(i)+(1.0-tau)*w
		kmaxindr(i)=count(fktr(i)>kspace)
		kmaxindw(i)=count(fktw(i)>kspace)
		if (kmaxindr(i)==0) kmaxindr(i)=1
		if (kmaxindr(i)==0) print *, 'kmaxindr(i)=0'
		if (kmaxindw(i)==0) kmaxindw(i)=1
		if (kmaxindw(i)==0) print *, 'kmaxindw(i)=0'
	end do
	

	! value fn and decision rule for maxage
	
	d(:,maxage)=1
	d1(:,maxage)=0.0
	v(:,maxage)=(/( ((1.0+r)*kspace(i)+pen)**(1.0-gamma)/(1.0-gamma),i=1,kgrid )/)
	
	
	
	! age loop
	do t=maxage-1,1,-1
		if (t>retage-1) fkt(:)=fktr(:)
		if (t<retage)	fkt(:)=fktw(:)
		
		
		! sequential method
		if (ialgo==1) then
			do i=1,kgrid
				if (t>retage-1) jmax=kmaxindr(i)
				if (t<retage)	jmax=kmaxindw(i)
				vmax=-1000.0
				js=1				
				do j=1,jmax
					vtemp=util(i,kspace(j))+beta*v(j,t+1)
					if (vtemp>=vmax) then
						vmax=vtemp
						js=j
	                end if
				end do
				v(i,t)=vmax
				d(i,t)=js
			end do
			go to 119
		end if
			
! 			print *,'prelinary optimal kt+1 is on grid js=',js
!			print *,'for kt on grid i=',i, 'vmax=',vmax
! 			pause

		! binary search (default)
		js=1
		do i=1,kgrid
			jmin=js
			if (t>retage-1) jmax=kmaxindr(i)
			if (t<retage)	jmax=kmaxindw(i)
			do while ((jmax-jmin)>2)
				jl=floor(real(jmin+jmax)/2.0)
				ju=jl+1
				v3(1)=util(i,kspace(jl))+beta*v(jl,t+1)
				v3(2)=util(i,kspace(ju))+beta*v(ju,t+1)
				if (v3(2)>v3(1)) jmin=jl
				if (v3(2)<=v3(1)) jmax=ju 
			end do
			v3(1)=util(i,kspace(jmin))+beta*v(jmin,t+1)
			v3(3)=util(i,kspace(jmax))+beta*v(jmax,t+1)
			if (jmax>jmin) v3(2)=util(i,kspace(jmin+1))+beta*v(jmin+1,t+1)
			if (jmax<=jmin) v3(2)=v3(1)
			js=jmin+(maxloc(v3,dim=1))-1
			


			
			! only binary, no interpolation
			if (ialgo==2) then
				v(i,t)=maxval(v3)
				d(i,t)=js
				go to 118
			end if
! 			print *,'preliminary optimal kt+1 is on grid js=',js
! 			print *,'for kt on grid i=',i
! ! 			pause
	
			

			! gss + linear interp & cubic spline (default)
			if (js==1) then	! left boundary
				if (fkt(i)<=kspace(1)) then	! there are cases where all ct<0
					d1(i,t)=kspace(1)
 					v(i,t)=maxval(v3)
					print *,'fkt(i)<=kspace(1)'
				else
					d1(i,t)=gss(1,2)	! m, i cannot be used in gss
					v(i,t)=vx
				end if
			else if (js==kgrid) then	! right boundary
				d1(i,t)=gss(kgrid-1,kgrid)
				v(i,t)=vx
			else	! middle
				d1(i,t)=gss(js-1,js+1)
				v(i,t)=vx
			end if
			
			
			
118			dummy=1
		end do ! kt loop
119		dummy=1
	end do ! t loop

	
	! aggregation
	ktp(1)=1
	ktp1(1)=1.0	! ktp(t) is position of starting kt at each age
	kcross(1)=0.0	! kcross(t) is exact value of starting kt at each age
	do t=2,maxage
		if (ialgo<3) then
			ktp(t)=d(ktp(t-1),t-1)
			kcross(t)=kspace(ktp(t))
		else
			temp=floor(ktp1(t-1))
			temp01=ktp1(t-1)-temp
			kcross(t)=(1.0-temp01)*d1(temp,t-1)+temp01*d1(temp+1,t-1)	! d1 is the value of the optimal kt+1
			ktp1(t)=(kcross(t)-kmin)/kstep+1.0
		end if
		sum=sum+kcross(t)
	end do
	
	
	! update
	tK=sum/real(maxage)
	kdiff=abs(K-tK)/K
	
 		print *, 'iteration :', iter
 		print *, 'K is:', K
 		print *, 'tK is:', tK
 		print *, 'kdiff is:', kdiff
	
	K=gradk*K+(1.0-gradk)*tK
    if (ialgo==4) call cspline
end do ! K loop

! K loop end &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&







!&&& part 4 screen print between each parameter set



 print *, 'gradk= ', gradk, 'kinit= ', kinit
! 
! if (kdiff<=tol) then
 	print *, 'Converged, iter= ', iter, 'kdiff= ', kdiff
 	print *, 'K= ', K, 'kgrid=', kgrid
! 	print *, ''
! 	print *, ''
! 
 	if (any(d==kgrid)) print *, 'Kt+1 has reached upper bound of state space'
 	if (any(d(:,1:maxage-1)==1)) print *, 'Kt+1 has reached lower bound of state space'
! 	print *, 'd(kgrid,maxage) is:'
! 	do i=1,kgrid
! 		write (7,888) d(i,:)
! 		print *, d(i,:)
! 		write (7,*) d(i,:)
! 	end do
! 
! end if
! 
! print *, ''
! print *, ''
! 

end do ! end gradk loop
	
end do ! end kinit loop

! parameter loop end &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&






!&&& part 5 internal process &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

contains

! utility

real(8) function util(i1,ktp11)	!i1 is integer, ktp11 is real

implicit none
integer i1
real(8) ct1,ktp11
if (t>retage-1) ct1=fktr(i)-ktp11
if (t<retage)	ct1=fktw(i)-ktp11

if (ct1>0) then
	util=ct1**(1.d0-gamma)/(1.d0-gamma)
else	! this occurs when kmaxind=1 & all ct<=0
	util=-10000.0
end if

end function







! golden section search

real(8) function gss(ia,ib)	! [ia,ib] is range of grid position to be searched

implicit none
integer ia,ib
real(8) ia1,ib1,a,b,c,a1,b1,fa1,fb1,ppp,via,vib,va,vb,vjs,va1,vb1
real(8) :: tol_1=0.0000001
! real kspace1(kgrid)
! kspace1=(/ ( kmin+real(i-1)*kstep, i=1,kgrid ) /)

888	format (F20.15,F20.15)

ppp=(sqrt(5.d0)-1.d0)/2.d0	! p=.618


ia1=kspace(ia)
a=ia1
via=v(ia,t+1)
va=via
ib1=kspace(ib)
b=ib1
vib=v(ib,t+1)
vb=vib
c=kspace(js)
vjs=v(js,t+1)	! vjs=va if js=1

! print *,'ia=',ia,'ib=',ib

! 888	format (F20.15,F20.15)


if (b>fkt(i)) then	! reset fkt(i) as new right bound
	if (ialgo==3) then
		vb=((b-fkt(i))*vjs+(fkt(i)-c)*vib)/kstep
	else
		vb=splint(fkt(i))
	end if
	b=fkt(i)
	print *,'b>fkt(i)'
end if

	a1=ppp*a+(1.d0-ppp)*b	!by gss
	b1=(1.d0-ppp)*a+ppp*b
	
if ((js>1).and.(js<kgrid)) then
	if (ialgo==3) then
		if (a1<c) va1=((a1-ia1)*vjs+(c-a1)*via)/kstep
		if (a1>=c) va1=((a1-c)*vib+(ib1-a1)*vjs)/kstep
		if (b1<c) vb1=((b1-ia1)*vjs+(c-b1)*via)/kstep
		if (b1>=c) vb1=((b1-c)*vib+(ib1-b1)*vjs)/kstep
	else
		va1=splint(a1)
		vb1=splint(b1)
	end if
else
	if (ialgo==3) then	
		va1=ppp*va+(1.d0-ppp)*vb
		vb1=(1.d0-ppp)*va+ppp*vb
	else
		va1=splint(a1)
		vb1=splint(b1)
	end if
end if

if (iter==10 .and. t==2) then
	if (i==2) then
		write(7,888) ia1,via
		write(7,888) ib1,vib
		write(7,888) c,vjs
		write(7,888) a1,va1
		write(7,888) b1,vb1
	end if
end if

fa1=util(i,a1)+beta*va1	!calculate util & lip value fn
fb1=util(i,b1)+beta*vb1

! print *,'a1=',a1,'b1=',b1


do while (abs(b-a)>tol_1)
	if ((js==1).or.(js==kgrid)) then
		if (fa1>fb1) then	! move towards left bound a, construct new point a1
			b=b1
			vb=vb1
			b1=a1
			vb1=va1
			fb1=fa1
			a1=ppp*a+(1.d0-ppp)*b
			if (ialgo==3) va1=ppp*va+(1.d0-ppp)*vb
			if (ialgo==4) va1=splint(a1)
			
			fa1=util(i,a1)+beta*va1	
! 			print *,'fa1>fb1'
! 			print *,'new a1=',a1,'b1=',b1
! 			pause
			if (iter==10 .and. t==2) then
				if (i==2) then
					write(7,888) a1,va1
				end if
			end if
		else	! move towards right bound b, construct new point b1
			a=a1
			va=va1
			a1=b1
			va1=vb1
			fa1=fb1
			b1=ppp*b+(1.d0-ppp)*a
			if (ialgo==3) vb1=ppp*vb+(1.d0-ppp)*va
			if (ialgo==4) vb1=splint(b1)
			fb1=util(i,b1)+beta*vb1
! 			print *,'fa1<=fb1'
! 			print *,'a1=',a1,'new b1=',b1
! 			pause
			if (iter==10 .and. t==2) then
				if (i==2) then
					write(7,888) b1,vb1
				end if
			end if
		end if
    else	! js<kgrid & js>1
		if (fa1>fb1) then
			b=b1
			b1=a1
			a1=ppp*a+(1.d0-ppp)*b
			if (ialgo==3) then
				if (a1<c) va1=((a1-ia1)*vjs+(c-a1)*via)/kstep
				if (a1>=c) va1=((a1-c)*vib+(ib1-a1)*vjs)/kstep
			else
				va1=splint(a1)
			end if
			fa1=util(i,a1)+beta*va1
! 			print *,'fa1>fb1'
! 			print *,'new a1=',a1,'b1=',b1
! 			pause
			if (iter==10 .and. t==2) then
				if (i==2) then
					write(7,888) a1,va1
				end if
			end if
		else
			a=a1
			a1=b1
			b1=ppp*b+(1.d0-ppp)*a
			if (ialgo==3) then
				if (b1<c) vb1=((b1-ia1)*vjs+(c-b1)*via)/kstep
				if (b1>=c) vb1=((b1-c)*vib+(ib1-b1)*vjs)/kstep
			else
				vb1=splint(b1)
			end if
			fb1=util(i,b1)+beta*vb1
! 			print *,'fa1<=fb1'
! 			print *,'a1=',a1,'new b1=',b1
! 			pause
			if (iter==10 .and. t==2) then
				if (i==2) then
					write(7,888) b1,vb1
				end if
			end if
		end if
	end if
end do
gss=(a1+b1)/2.d0
vx=(fa1+fb1)/2.d0

! print *,'gss=',gss,'vx=',vx
! print *,'************************************************************'


end function



! subroutine cspline

subroutine cspline
implicit none
integer i12


vsod=0.d0


vsod(1,:)=0.d0	! each age has a unique valeu fn, so each valeu fn line needs to be approximated with seperate 2nd order derivative vector
uu(1,:)=0.d0
qn=0.d0
un=0.d0

do t=1,maxage
	do i12=2,kgrid-1
		sig=(kspace(i12)-kspace(i12-1))/(kspace(i12+1)-kspace(i12-1))	! sig=hi/(hi+hi+1) on lecture note
		p=sig*vsod(i12-1,t)+2.d0
		vsod(i12,t)=(sig-1.d0)/p
		uu(i12,t)=(6.d0*((v(i12+1,t)-v(i12,t))/(kspace(i12+1)-kspace(i12))-(v(i12,t)-v(i12-1,t))/(kspace(i12)-kspace(i12-1)))/(kspace(i12+1)-kspace(i12-1))-sig*uu(i12-1,t))/p
	end do

	vsod(kgrid,t)=(un-qn*uu(kgrid-1,t))/(qn*vsod(kgrid-1,t)+1.d0)

	do i12=kgrid-1,1,-1
		vsod(i12,t)=vsod(i12,t)*vsod(i12+1,t)+uu(i12,t)
	end do
end do

end subroutine



! function splint
real(8) function splint(a1_1)

implicit none
real(8) a1_1,pa,pb
integer ia_1
ia_1=floor((a1_1-kspace(1))/kstep)+1
pa=(kspace(ia_1+1)-a1_1)/kstep
pb=1.d0-pa
splint=pa*v(ia_1,t+1)+pb*v(ia_1+1,t+1)+((pa**3-pa)*vsod(ia_1,t+1)+(pb**3-pb)*vsod(ia_1+1,t+1))*(kstep**2)/6.d0

end function











end program