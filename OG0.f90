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

! utility fn
real :: betam(4)=(/ 0.978,0.95,0.995,1.011 /)
real :: gammam(4)=(/ 1.25,2.0,5.0,8.0 /)

! production fn
real :: alpham(4)=(/ 0.5,0.63,0.69,0.729 /)	! labor's income share
real,parameter :: tfp=1.005
real,parameter :: psi=0.232		! capital's income share 0.277 or 0.232

! policy
real :: thetam(10)=(/ 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 /)
real :: phim(10)=(/ 0.0,0.1,0.2,0.30,0.4,0.5,0.6,0.7,0.8,0.9 /)

! data
real :: deltam(4)=(/ 0.064,0.069,0.08,0.125 /)
real :: growthm(4)=(/ 0.01,0.0165,0.02,0.05 /)
real :: rhom(4)=(/ 0.008,0.012,0.015,0.021 /)


! age structure
integer,parameter :: maxage=65
integer,parameter :: retage=45

! state space
real,parameter :: kmax=32.0
real,parameter :: kmin=0.0
integer,parameter :: kgrid=107

! loop

integer,parameter :: maxiter=51
real,parameter :: tol=0.001

! initial guess
real,parameter :: kinitmax=5.0
real,parameter :: kinitmin=0.1
integer,parameter :: kinitgrid=50

! real,parameter :: binit=0.042

real,parameter :: binitmax=0.5
real,parameter :: binitmin=0.01
integer,parameter :: binitgrid=50


! parameters



real :: gradm(9)=(/ 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 /)

! real :: tau=theta/(2.0+theta)	!tau=0, so pen=0 too
! real :: L=real(retage-1)/real(maxage)




! basic
real(8) vr(kgrid,retage:maxage),vw(2,kgrid,retage-1),dr1(kgrid,retage:maxage),dw1(2,kgrid,retage-1)
real(8) K,w,r,tK,d1(kgrid,maxage),L,beq,beqm(maxage),tbeq
real(8) pr(kgrid,retage:maxage),pw(2,kgrid,retage-1)
integer dr(kgrid,retage:maxage),dw(2,kgrid,retage-1)

! parameter experiment 1
real(8) beta,gamma,alpha,theta,phi,delta,g,rho,utax,ag
integer bm,gam,am,tm,pm,dm,gm,rm

! feature
real(8) pen,sumeff,temp1,pentax
real(8) effcross(retage-1),efflong(retage-1),ptrans(2,2)

! population
real(8) consurv(maxage),abssurv(maxage),sumpop,mu(maxage)


! parameter experiment 2
integer gkm,gbm,kinitm,binitm
real(8) kinitstep,kinit,kinitspace(kinitgrid),gradk
real(8) binitstep,binit,binitspace(binitgrid)


! indexing spaces
real(8) fktr(kgrid,retage:maxage),fktw(2,kgrid,retage-1),fkt(kgrid),kstep,kspace(kgrid)
integer kmaxindr(kgrid,retage:maxage),kmaxindw(2,kgrid,retage-1),maxgrid


! loop
real(8) kdiff,gradb,beqdiff
integer i,j,iter,t,m,iterkmaxind,itermaxgrid,iterlowinterest,iterp

! aggregation
real(8) kcross(maxage),ktp1(maxage),sumk,temp01,klong(maxage),sumbeq,output,price
integer temp,ktp(maxage)

! sequential method
real(8) vtemp,vmax

! binary search
real(8) v3(3)
integer js,jmin,jmax,jl,ju,jsm(2)

! gss
real(8) vx

! cspline
real(8) sig,qn,un,uu(kgrid,maxage),p,vsod(kgrid,maxage)

! else
integer dummy
real(8) cons







! open file

if (isys==1) open(unit=7,file='C:\Users\NREM\Desktop\Dropbox\cversion\comboeff.txt')
if (isys==1) open(unit=8,file='C:\Users\NREM\Desktop\Dropbox\cversion\psurv.txt')
if (isys==1) open(unit=10,file='C:\Users\NREM\Desktop\Dropbox\cversion\result.txt')

if (isys==2) open(unit=7,file='comboeff.txt')
if (isys==2) open(unit=8,file='psurv.txt')
if (isys==2) open(unit=10,file='result.txt')

777 format( 4(F5.2,1X) 2(F7.4,1X) 4(I6,1X) )
888 format( 8(F6.3,1X) 2(F9.6,1X) 4(I6,1X) )
! kspace





read(7,*) ( effcross(t),t=1,retage-1 )


ptrans=reshape((/0.94,0.94,0.06,0.06/),(/2,2/))
read(8,*) ( consurv(t), t=1,maxage )

abssurv(1) = 1.0
do t=2,maxage
	abssurv(t) = abssurv(t-1)*consurv(t)
end do










! print *, kspace
! print *, kinitspace
! pause


! 88 format (4(I5,2X))

! parameter loop begin &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

kinitstep=(kinitmax-kinitmin)/real(kinitgrid-1)
kinitspace=(/ ( kinitmin+real(i-1)*kinitstep, i=1,kinitgrid ) /)


binitstep=(binitmax-binitmin)/real(binitgrid-1)
binitspace=(/ ( binitmin+real(i-1)*binitstep, i=1,binitgrid ) /)
	
kstep=(kmax-kmin)/real(kgrid-1)
kspace=(/ ( kmin+real(i-1)*kstep, i=1,kgrid ) /)

iterp=0

! data
do gm=2,2 ! default =2
	g=growthm(gm)

    efflong=(/( ((1.0+g)**(t-1))*effcross(t), t=1,retage-1 )/)
    
do rm=2,2 ! default =2
	rho=rhom(rm)

	sumpop=0.0
	do t=1,maxage
		sumpop=sumpop+abssurv(t)/((1.0+rho)**(t-1))
	end do
	mu(1)=1.0/sumpop
	do t=2,maxage
		mu(t)=consurv(t)*mu(t-1)/(1.0+rho)
	end do

	L=0.0
	do t=1,retage-1
		L=L+0.94*effcross(t)*mu(t)
	end do

    
! parameter
do gam=2,2 ! default =2
	gamma=gammam(gam)
	
do dm=1,1 ! default =2
	delta=deltam(dm)	
	
do bm=1,1 ! default =1
	beta=betam(bm)
	
do am=4,4 ! default =3
	alpha=alpham(am)


! tax
do tm=4,4 ! default =4
	theta=thetam(tm)
	
do pm=4,4 ! default =4
	phi=phim(pm)


! algorithm
do binitm=4,4 ! default =4
	binit=binitspace(binitm)

do gbm=4,4 ! default =4
	gradb=gradm(gbm)
	
do kinitm=18,18 ! default =18
	kinit=kinitspace(kinitm)

do gkm=8,8 ! default =8
	gradk=gradm(gkm)



    
    

iterp=iterp+1
utax=(0.06/0.94)*phi
ag=(1.0+g)*(1.0+rho)-1.0
	
!&&& part 2
! initialize for loop

	iter=0
    iterkmaxind=0
    itermaxgrid=0
	iterlowinterest=0
	kdiff=10.0
	
	K=kinit
	beq=binit
	
!	if (ialgo==4) call cspline
	
	
! K loop begin &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&& part 3 value fn loop

do while((kdiff>tol).and.(iter<=maxiter))
	beqm=(/( beq*(1.0+g)**(t-1),t=1,maxage )/)
	iter=iter+1
! 	w=(1.0-alpha)*(K**alpha)*(L**(-alpha))
! 	r=alpha*(K**(alpha-1.0))*(L**(1.0-alpha))-delta
	
	w=alpha*tfp*(L**(alpha-1.0))*(K**psi)
	r=psi*tfp*(L**alpha)*(K**(psi-1.0))-delta
	
	
	if (r<ag) then
! 		print *, "Interest rate is less than g rate."
! 		write (10,*), 'Interest rate is less than g rate.'
		iterlowinterest=iterlowinterest+1
	end if

	sumeff=0.0
	do t=1,retage-1
		sumeff=sumeff+efflong(t)
	end do
	pen=sumeff/(real(retage-1))*theta*w
	
	temp1=0.0
	do t=retage,maxage
		temp1=temp1+mu(t)*(1.0+g)**(1-t)
	end do
	pentax=pen*temp1
	temp1=0.0
	do t=1,retage-1
		temp1=temp1+mu(t)*effcross(t)
	end do
	pentax=pentax/(0.94*w*temp1)
	
! 	pen=theta*(1.0-tau)*w
	
	! upper limit for Kt+1 known after rt, wt and pen are known
	! Ct+Kt+1=(1+rt)Kt+(1-tau)wt 		workers upper limit fktw
	! Ct+Kt+1=(1+rt)Kt+pen				retirees upper limit fktr
	! all workers upper limit for Kt+1 are same
	! all retirees upper limit for Kt+1 are same
	
	
	do t=1,maxage
		do i=1,kgrid
			if (t>=retage) then
				fktr(i,t)=(1.0+r)*( kspace(i) + beqm(t) ) + pen
				kmaxindr(i,t)=count(fktr(i,t)>kspace)
				if (kmaxindr(i,t)==0) then
!                     print *, 'kmaxindr(i,t)=0'
				    !if (kmaxindr(i,t)==0) write (10,*), 'kmaxindr(i,t)=0'
				    kmaxindr(i,t)=1
                    iterkmaxind=iterkmaxind+1
                end if
			else
				fktw(1,i,t)=(1.0+r)*( kspace(i) + beqm(t) ) + (1.0-pentax-utax)*w*efflong(t)
				fktw(2,i,t)=(1.0+r)*( kspace(i) + beqm(t) ) + phi*w*efflong(t)
				kmaxindw(1,i,t)=count(fktw(1,i,t)>kspace)
				kmaxindw(2,i,t)=count(fktw(2,i,t)>kspace)
				if (kmaxindw(1,i,t)==0) then
!                     print *, 'kmaxindw(1,i,t)=0'
				    !if (kmaxindw(1,i,t)==0) write (10,*), 'kmaxindw(1,i,t)=0'
				    kmaxindw(1,i,t)=1
                    iterkmaxind=iterkmaxind+1
                end if
                if (kmaxindw(2,i,t)==0) then
!                     print *, 'kmaxindw(2,i,t)=0'
				    !if (kmaxindw(2,i,t)==0) write (10,*), 'kmaxindw(2,i,t)=0'
				    kmaxindw(2,i,t)=1
                    iterkmaxind=iterkmaxind+1
                end if
			end if
		end do
	end do
	

	! value fn and decision rule for maxage
	
	dr(:,maxage)=1
	dr1(:,maxage)=0.0
	
	vr(:,maxage)=(/( fktr(i,maxage)**(1.0-gamma)/(1.0-gamma), i=1,kgrid )/)
	
	
	
	! age loop
	do t=maxage-1,1,-1
! 		if (t>=retage) fkt(:)=fktr(:)
! 		if (t<retage)	fkt(:)=fktw(:)
		
		
		! sequential method
! 		if (ialgo==1) then
! 			do i=1,kgrid
! 				if (t>=retage) jmax=kmaxindr(i)
! 				if (t<retage)	jmax=kmaxindw(i)
! 				vmax=-1000.0
! 				js=1				
! 				do j=1,jmax
! 					vtemp=util(i,kspace(j))+beta*v(j,t+1)
! 					if (vtemp>=vmax) then
! 						vmax=vtemp
! 						js=j
! 	                end if
! 				end do
! 				v(i,t)=vmax
! 				d(i,t)=js
! 			end do
! 			go to 119
! 		end if
! 			
! 			print *,'prelinary optimal kt+1 is on grid js=',js
!			print *,'for kt on grid i=',i, 'vmax=',vmax
! 			pause

		! binary search (default)
		js=1
        jsm=1
		do i=1,kgrid
			
			
			if (t>=retage) then ! all retirees
				jmin=js
				jmax=kmaxindr(i,t)
				do while ((jmax-jmin)>2)
					jl=floor(real(jmin+jmax)/2.0)
					ju=jl+1
					v3(1)=util(i,kspace(jl))+beta*consurv(t+1)*vr(jl,t+1)
					v3(2)=util(i,kspace(ju))+beta*consurv(t+1)*vr(ju,t+1)
					if (v3(2)>v3(1)) jmin=jl
					if (v3(2)<=v3(1)) jmax=ju
				end do
				v3(1)=util(i,kspace(jmin))+beta*consurv(t+1)*vr(jmin,t+1)
				v3(3)=util(i,kspace(jmax))+beta*consurv(t+1)*vr(jmax,t+1)				
				if (jmax>jmin) v3(2)=util(i,kspace(jmin+1))+beta*consurv(t+1)*vr(jmin+1,t+1)
				if (jmax<=jmin) v3(2)=v3(1)
				js=jmin+(maxloc(v3,dim=1))-1
                if (js==0) print *, 'js=0'
                if (js==0) pause
				
				! only binary, no interpolation
				if (ialgo==2) vr(i,t)=maxval(v3)
				if (ialgo==2) dr(i,t)=js
				if (ialgo==2) go to 118
				
			else if (t==retage-1) then	! agent is 44, retire next period
				
				do m=1,2
					jmin=jsm(m)
					jmax=kmaxindw(m,i,t)
					do while ((jmax-jmin)>2)
						jl=floor(real(jmin+jmax)/2.0)
						ju=jl+1
						v3(1)=util(i,kspace(jl))+beta*consurv(t+1)*vr(jl,t+1)
						v3(2)=util(i,kspace(ju))+beta*consurv(t+1)*vr(ju,t+1)
						if (v3(2)>v3(1)) jmin=jl
						if (v3(2)<=v3(1)) jmax=ju
					end do
					v3(1)=util(i,kspace(jmin))+beta*consurv(t+1)*vr(jmin,t+1)
					v3(3)=util(i,kspace(jmax))+beta*consurv(t+1)*vr(jmax,t+1)				
					if (jmax>jmin) v3(2)=util(i,kspace(jmin+1))+beta*consurv(t+1)*vr(jmin+1,t+1)
					if (jmax<=jmin) v3(2)=v3(1)
					jsm(m)=jmin+(maxloc(v3,dim=1))-1
                    if (jsm(m)==0) print *, 'jsm(m)=0'
                    if (jsm(m)==0) pause
					
					! only binary, no interpolation
					if (ialgo==2) vw(m,i,t)=maxval(v3)
					if (ialgo==2) dw(m,i,t)=jsm(m)
                end do
                if (ialgo==2) go to 118
				
			else if (t<retage-1) then	! all working
				
				do m=1,2
					jmin=jsm(m)
					jmax=kmaxindw(m,i,t)
					do while ((jmax-jmin)>2)
						jl=floor(real(jmin+jmax)/2.0)
						ju=jl+1
						v3(1)=util(i,kspace(jl))+beta*consurv(t+1)*( ptrans(m,1)*vw(1,jl,t+1) + ptrans(m,2)*vw(2,jl,t+1) )
						v3(2)=util(i,kspace(ju))+beta*consurv(t+1)*( ptrans(m,1)*vw(1,ju,t+1) + ptrans(m,2)*vw(2,ju,t+1) )
						if (v3(2)>v3(1)) jmin=jl
						if (v3(2)<=v3(1)) jmax=ju
					end do
					v3(1)=util(i,kspace(jmin))+beta*consurv(t+1)*( ptrans(m,1)*vw(1,jmin,t+1) + ptrans(m,2)*vw(2,jmin,t+1) )
					v3(3)=util(i,kspace(jmax))+beta*consurv(t+1)*( ptrans(m,1)*vw(1,jmax,t+1) + ptrans(m,2)*vw(2,jmax,t+1) )
					if (jmax>jmin) v3(2)=util(i,kspace(jmin+1))+beta*consurv(t+1)*( ptrans(m,1)*vw(1,jmin+1,t+1) + ptrans(m,2)*vw(2,jmin+1,t+1) )
					if (jmax<=jmin) v3(2)=v3(1)
					jsm(m)=jmin+(maxloc(v3,dim=1))-1
                    if (jsm(m)==0) print *, 'jsm(m)=0'
                    if (jsm(m)==0) pause
                    
					! only binary, no interpolation
					if (ialgo==2) vw(m,i,t)=maxval(v3)
					if (ialgo==2) dw(m,i,t)=jsm(m)
                end do
                if (ialgo==2) go to 118
				
			end if
			

! 			print *,'preliminary optimal kt+1 is on grid js=',js
! 			print *,'for kt on grid i=',i
! ! 			pause
	
			

			! gss + linear interp & cubic spline (default)
! 			if (js==1) then	! left boundary
! 				if (fkt(i)<=kspace(1)) then	! there are cases where all ct<0
! 					d1(i,t)=kspace(1)
!  					v(i,t)=maxval(v3)
! 					print *,'fkt(i)<=kspace(1)'
! 				else
! 					d1(i,t)=gss(1,2)	! m, i cannot be used in gss
! 					v(i,t)=vx
! 				end if
! 			else if (js==kgrid) then	! right boundary
! 				d1(i,t)=gss(kgrid-1,kgrid)
! 				v(i,t)=vx
! 			else	! middle
! 				d1(i,t)=gss(js-1,js+1)
! 				v(i,t)=vx
! 			end if
			
			
			
118			dummy=1
		end do ! kt loop
119		dummy=1
	end do ! t loop
	
	
! 	print *, 'max asset obtained in this loop', maxval( maxval( dw(1,:,:),dim=1 ) )
	
	pw=0.0
	pr=0.0
	pw(1,1,1)=0.94
	pw(2,1,1)=0.06
	
	! working agents
	do t=2,retage-1
		do i=1,kgrid
			do m=1,2
				j=dw(m,i,t-1)
				pw(1,j,t)=pw(1,j,t)+pw(m,i,t-1)*ptrans(m,1)
				pw(2,j,t)=pw(2,j,t)+pw(m,i,t-1)*ptrans(m,2)
			end do
		end do
	end do
	
	! new retiree
	
	do i=1,kgrid
		do m=1,2
			j=dw(m,i,retage-1)
			pr(j,retage)=pr(j,retage)+pw(m,i,retage-1)
		end do
	end do
	
	! old retiree
	do t=retage+1,maxage
		do i=1,kgrid
			j=dr(i,t-1)
			pr(j,t)=pr(j,t)+pr(i,t-1)
		end do
	end do
	
	
	! check
! 	if ( any(dr==kgrid) .or. any(dw==kgrid) ) print *, 'Kt+1 has reached upper bound of state space'
! 	if ( any(dw==1) .or. any( dr(:,retage:maxage-1)==1 ) ) print *, 'Kt+1 has reached lower bound of state space'
	
	maxgrid=0
	
	do t=1,retage-1
		do i=1,kgrid
			if ( dw(1,i,t)>maxgrid .and. pw(1,i,t)>0 ) maxgrid=dw(1,i,t)
			if ( dw(1,i,t)>=kgrid .and. pw(1,i,t)>0 ) then
! 				print *, 'Kt+1 reached upper bound'
				!write (10,*), 'Kt+1 reached upper bound'
                itermaxgrid=itermaxgrid+1
! 				print *, 'age =',t
! 				print *, 'Kt =',i
!				print *, 'proportion is',pw(1,i,t)
			end if
		end do
	end do
! 	
! 	print *, 'max asset grid =',maxgrid
	

	! compute longitudinal profiles for a given cohort
	klong=0.0
	
	do t=1,retage-1
		do i=1,kgrid
			do m=1,2
				j=dw(m,i,t)
				klong(t)=klong(t)+kspace(j)*pw(m,i,t)
			end do
		end do
	end do
	
	do t=retage,maxage
		do i=1,kgrid
			j=dr(i,t)
			klong(t)=klong(t)+kspace(j)*pr(i,t)
		end do
	end do
	
	kcross=(/( klong(t)*(1.0+g)**(1-t), t=1,maxage )/)
	
	
	! aggregation
! 	ktp(1)=1
! 	ktp1(1)=1.0	! ktp(t) is position of starting kt at each age
! 	kcross(1)=0.0	! kcross(t) is exact value of starting kt at each age
! 	do t=2,maxage
! 		if (ialgo<3) then
! 			ktp(t)=d(ktp(t-1),t-1)
! 			kcross(t)=kspace(ktp(t))
! 		else
! 			temp=floor(ktp1(t-1))
! 			temp01=ktp1(t-1)-temp
! 			kcross(t)=(1.0-temp01)*d1(temp,t-1)+temp01*d1(temp+1,t-1)	! d1 is the value of the optimal kt+1
! 			ktp1(t)=(kcross(t)-kmin)/kstep+1.0
! 		end if
! 		sumk=sumk+kcross(t)
! 	end do
	
	sumk=0.0
	sumbeq=0.0
	
	do t=1,maxage-1
		sumk=sumk+kcross(t)*mu(t)
		sumbeq=sumbeq+kcross(t)*mu(t)*(1.0-consurv(t+1))
	end do
	
	
	
	! update
! 	tK=sumk/real(maxage)
	output=tfp*(L**alpha)*(K**psi)
	price=(1.0-alpha-psi)*((1.0+ag)/(r-ag))*output
	tk=(sumk-price)/(1.0+ag)
	tbeq=sumbeq/(1.0+ag)
	
	kdiff=abs(K-tK)/K
	beqdiff=abs(beq-tbeq)/beq
! 	
!  		print *, 'iteration :', iter
!  		print *, 'K is:', K, 'tK is:', tK
! 		print *, 'kidff is:',kdiff
!  		print *, 'beq is:', beq, 'tbeq is:', tbeq
!  		print *, 'beqdiff is:', beqdiff
! 		print *, '                     '
	K=gradk*K+(1.0-gradk)*tK
	beq=gradb*beq+(1.0-gradb)*tbeq
	
!    if (ialgo==4) call cspline
end do ! K loop


! K loop end &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&







!&&& part 4 screen print between each parameter set



!  print *, 'gradk= ', gradk, 'kinit= ', kinit
! ! 
! ! if (kdiff<=tol) then
!  	print *, 'Converged, iter= ', iter, 'kdiff= ', kdiff
!  	print *, 'K= ', K, 'kgrid=', kgrid
! 	print *, ''
! 	print *, ''
! 
!  	if (any(d==kgrid)) print *, 'Kt+1 has reached upper bound of state space'
!  	if (any(d(:,1:maxage-1)==1)) print *, 'Kt+1 has reached lower bound of state space'
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
! print *, '                             '
!print *, 'gradk =',gradk
print *, 'iterp =',iterp
print *, 'gamma =',gamma
print *, 'delta =',delta
print *, 'beta =',beta
print *, 'alpha =',alpha
! print *, 'theta =',theta
! print *, 'phi =',phi
! print *, '                             '
!print *, "Initial K =",kinit
print *, 'Final K =',K
!print *, 'Initial bequest =',binit
print *, 'Final bequest',beq
print *, 'Iteration',iter
print *, '_____________________________'
! print *, '                             '
! print *, '                             '
! write (10,*) Kinit,K

!write (10,*), '                             '
! write (10,*), 'gradk =',gradk
! write (10,*), '                             '
! write (10,*), "Initial K =",kinit
! write (10,*), 'Final K =',K
! write (10,*), 'Initial bequest =',binit
! write (10,*), 'Final bequest',beq
! write (10,*), 'Iteration',iter
! write (10,*), 'binit,beq,K,iter'
!write (10,777) kinit,gradk,K,beq,iter,iterkmaxind,itermaxgrid,iterlowinterest
write (10,888) g,rho,gamma,delta,beta,alpha,theta,phi,K,beq,iter,iterkmaxind,itermaxgrid,iterlowinterest

!write (10,*), '_____________________________'
!write (10,*), '                             '
!write (10,*), '                             '

end do ! gradk loop
	
end do ! kinit loop

end do ! gradb loop

end do ! binit loop





end do ! phi loop

end do ! theta loop



end do ! alpha loop

end do ! beta loop

end do ! delta loop

end do ! gama loop


end do ! rho loop

end do ! g loop

! parameter loop end &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&






!&&& part 5 internal process &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

contains

! utility

real(8) function util(i1,ktp11)	!i1 is integer, ktp11 is real

implicit none
integer i1
real(8) ct1,ktp11
if (t>=retage) ct1=fktr(i,t)-ktp11
if (t<retage) ct1=fktw(m,i,t)-ktp11

if (ct1>0) then
	util=ct1**(1.d0-gamma)/(1.d0-gamma)
else	! this occurs when kmaxind=1 & all ct<=0
	util=-10000000.0
end if

! if (i==64 .and. i==1) print *, 'util =', util
end function






! 
! ! golden section search
! 
! real(8) function gss(ia,ib)	! [ia,ib] is range of grid position to be searched
! 
! implicit none
! integer ia,ib
! real(8) ia1,ib1,a,b,c,a1,b1,fa1,fb1,ppp,via,vib,va,vb,vjs,va1,vb1
! real(8) :: tol_1=0.0000001
! ! real kspace1(kgrid)
! ! kspace1=(/ ( kmin+real(i-1)*kstep, i=1,kgrid ) /)
! 
! 888	format (F20.15,F20.15)
! 
! ppp=(sqrt(5.d0)-1.d0)/2.d0	! p=.618
! 
! 
! ia1=kspace(ia)
! a=ia1
! via=v(ia,t+1)
! va=via
! ib1=kspace(ib)
! b=ib1
! vib=v(ib,t+1)
! vb=vib
! c=kspace(js)
! vjs=v(js,t+1)	! vjs=va if js=1
! 
! ! print *,'ia=',ia,'ib=',ib
! 
! ! 888	format (F20.15,F20.15)
! 
! 
! if (b>fkt(i)) then	! reset fkt(i) as new right bound
! 	if (ialgo==3) then
! 		vb=((b-fkt(i))*vjs+(fkt(i)-c)*vib)/kstep
! 	else
! 		vb=splint(fkt(i))
! 	end if
! 	b=fkt(i)
! 	print *,'b>fkt(i)'
! end if
! 
! 	a1=ppp*a+(1.d0-ppp)*b	!by gss
! 	b1=(1.d0-ppp)*a+ppp*b
! 	
! if ((js>1).and.(js<kgrid)) then
! 	if (ialgo==3) then
! 		if (a1<c) va1=((a1-ia1)*vjs+(c-a1)*via)/kstep
! 		if (a1>=c) va1=((a1-c)*vib+(ib1-a1)*vjs)/kstep
! 		if (b1<c) vb1=((b1-ia1)*vjs+(c-b1)*via)/kstep
! 		if (b1>=c) vb1=((b1-c)*vib+(ib1-b1)*vjs)/kstep
! 	else
! 		va1=splint(a1)
! 		vb1=splint(b1)
! 	end if
! else
! 	if (ialgo==3) then	
! 		va1=ppp*va+(1.d0-ppp)*vb
! 		vb1=(1.d0-ppp)*va+ppp*vb
! 	else
! 		va1=splint(a1)
! 		vb1=splint(b1)
! 	end if
! end if
! 
! if (iter==10 .and. t==2) then
! 	if (i==2) then
! 		write(7,888) ia1,via
! 		write(7,888) ib1,vib
! 		write(7,888) c,vjs
! 		write(7,888) a1,va1
! 		write(7,888) b1,vb1
! 	end if
! end if
! 
! fa1=util(i,a1)+beta*va1	!calculate util & lip value fn
! fb1=util(i,b1)+beta*vb1
! 
! ! print *,'a1=',a1,'b1=',b1
! 
! 
! do while (abs(b-a)>tol_1)
! 	if ((js==1).or.(js==kgrid)) then
! 		if (fa1>fb1) then	! move towards left bound a, construct new point a1
! 			b=b1
! 			vb=vb1
! 			b1=a1
! 			vb1=va1
! 			fb1=fa1
! 			a1=ppp*a+(1.d0-ppp)*b
! 			if (ialgo==3) va1=ppp*va+(1.d0-ppp)*vb
! 			if (ialgo==4) va1=splint(a1)
! 			
! 			fa1=util(i,a1)+beta*va1	
! ! 			print *,'fa1>fb1'
! ! 			print *,'new a1=',a1,'b1=',b1
! ! 			pause
! 			if (iter==10 .and. t==2) then
! 				if (i==2) then
! 					write(7,888) a1,va1
! 				end if
! 			end if
! 		else	! move towards right bound b, construct new point b1
! 			a=a1
! 			va=va1
! 			a1=b1
! 			va1=vb1
! 			fa1=fb1
! 			b1=ppp*b+(1.d0-ppp)*a
! 			if (ialgo==3) vb1=ppp*vb+(1.d0-ppp)*va
! 			if (ialgo==4) vb1=splint(b1)
! 			fb1=util(i,b1)+beta*vb1
! ! 			print *,'fa1<=fb1'
! ! 			print *,'a1=',a1,'new b1=',b1
! ! 			pause
! 			if (iter==10 .and. t==2) then
! 				if (i==2) then
! 					write(7,888) b1,vb1
! 				end if
! 			end if
! 		end if
!     else	! js<kgrid & js>1
! 		if (fa1>fb1) then
! 			b=b1
! 			b1=a1
! 			a1=ppp*a+(1.d0-ppp)*b
! 			if (ialgo==3) then
! 				if (a1<c) va1=((a1-ia1)*vjs+(c-a1)*via)/kstep
! 				if (a1>=c) va1=((a1-c)*vib+(ib1-a1)*vjs)/kstep
! 			else
! 				va1=splint(a1)
! 			end if
! 			fa1=util(i,a1)+beta*va1
! ! 			print *,'fa1>fb1'
! ! 			print *,'new a1=',a1,'b1=',b1
! ! 			pause
! 			if (iter==10 .and. t==2) then
! 				if (i==2) then
! 					write(7,888) a1,va1
! 				end if
! 			end if
! 		else
! 			a=a1
! 			a1=b1
! 			b1=ppp*b+(1.d0-ppp)*a
! 			if (ialgo==3) then
! 				if (b1<c) vb1=((b1-ia1)*vjs+(c-b1)*via)/kstep
! 				if (b1>=c) vb1=((b1-c)*vib+(ib1-b1)*vjs)/kstep
! 			else
! 				vb1=splint(b1)
! 			end if
! 			fb1=util(i,b1)+beta*vb1
! ! 			print *,'fa1<=fb1'
! ! 			print *,'a1=',a1,'new b1=',b1
! ! 			pause
! 			if (iter==10 .and. t==2) then
! 				if (i==2) then
! 					write(7,888) b1,vb1
! 				end if
! 			end if
! 		end if
! 	end if
! end do
! gss=(a1+b1)/2.d0
! vx=(fa1+fb1)/2.d0
! 
! ! print *,'gss=',gss,'vx=',vx
! ! print *,'************************************************************'
! 
! 
! end function
! 
! 
! 
! ! subroutine cspline
! 
! subroutine cspline
! implicit none
! integer i12
! 
! 
! vsod=0.d0
! 
! 
! vsod(1,:)=0.d0	! each age has a unique valeu fn, so each valeu fn line needs to be approximated with seperate 2nd order derivative vector
! uu(1,:)=0.d0
! qn=0.d0
! un=0.d0
! 
! do t=1,maxage
! 	do i12=2,kgrid-1
! 		sig=(kspace(i12)-kspace(i12-1))/(kspace(i12+1)-kspace(i12-1))	! sig=hi/(hi+hi+1) on lecture note
! 		p=sig*vsod(i12-1,t)+2.d0
! 		vsod(i12,t)=(sig-1.d0)/p
! 		uu(i12,t)=(6.d0*((v(i12+1,t)-v(i12,t))/(kspace(i12+1)-kspace(i12))-(v(i12,t)-v(i12-1,t))/(kspace(i12)-kspace(i12-1)))/(kspace(i12+1)-kspace(i12-1))-sig*uu(i12-1,t))/p
! 	end do
! 
! 	vsod(kgrid,t)=(un-qn*uu(kgrid-1,t))/(qn*vsod(kgrid-1,t)+1.d0)
! 
! 	do i12=kgrid-1,1,-1
! 		vsod(i12,t)=vsod(i12,t)*vsod(i12+1,t)+uu(i12,t)
! 	end do
! end do
! 
! end subroutine
! 
! 
! 
! ! function splint
! real(8) function splint(a1_1)
! 
! implicit none
! real(8) a1_1,pa,pb
! integer ia_1
! ia_1=floor((a1_1-kspace(1))/kstep)+1
! pa=(kspace(ia_1+1)-a1_1)/kstep
! pb=1.d0-pa
! splint=pa*v(ia_1,t+1)+pb*v(ia_1+1,t+1)+((pa**3-pa)*vsod(ia_1,t+1)+(pb**3-pb)*vsod(ia_1+1,t+1))*(kstep**2)/6.d0
! 
! end function
! 
! 










end program