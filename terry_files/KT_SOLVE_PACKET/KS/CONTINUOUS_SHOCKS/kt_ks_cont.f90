!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kt_ks_cont.f90
!
! Fortran code for the KS solution of the Khan and Thomas (2008) 
! model, allowing for continuous aggregate shocks.
!
! 'Alternative Methods for Solving Heterogeneous Firm Models'
! Stephen Terry (2015)
!
! This Version : 12/18/15
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module modparams
implicit none
    
!This module contains fixed model parameter values, available
!to the rest of the code. It also declares constants which are made available
!to subprograms when set elsewhere. It also declares some variables
!threadprivate, which makes them unique to OpenMP parallel threads.
integer, parameter :: knum = 10 !grid size for idio capital
integer, parameter :: znum = 5  !grid size for idio productivity
integer, parameter :: anum = 5  !grid size for agg productivity
integer, parameter :: kbarnum = 10 !grid size for agg capital
integer, parameter :: kdensenum = 50 !grid size for idio capital on dist
integer, parameter :: numper = 2500 !number of periods in unconditional simulation
integer, parameter :: numdiscard = 500 !number of periods discarded from unconditional simulation 
integer, parameter :: seedint=2503 !random number seed
integer, parameter :: maxvfit = 1000 !max number of vf iters
integer, parameter :: maxaccelit = 50 !number of howard accelerations
integer, parameter :: maxpit=100 !max number of price iterations when clearing market
integer, parameter :: maxfcstit=30 !max number of GE forecast rule iterations
integer, parameter :: maxergit = 5000 !max number of iters to compute ergodic dist of agg prod
integer, parameter :: pcutoff = 17 !iterations at which to restart price clearing algorithms
integer, parameter :: nummicro = 7 !number of micro moments to compute
integer, parameter :: doVFI = 1 !complete VF iteration or read from files?
integer, parameter :: numbeta = 4 !number of coeffs in fcst rule

double precision, parameter :: ainit = 1.0 !initial agg prod value
double precision, parameter :: alpha = 0.256 !capital elasticity
double precision, parameter :: nu = 0.640 !labor elasticity
double precision, parameter :: phi = 2.4 !disutility of labor
double precision, parameter :: xibar = 0.0083 !upper bound of capital AC dist
double precision, parameter :: beta = 0.977 !discount rate
double precision, parameter :: delta = 0.069 !capital depreciation rate
double precision, parameter :: kmin= 0.1 !min of idio capital grid
double precision, parameter :: kmax = 8.0 !max of idio capital grid
double precision, parameter :: kbarmin =1.25 !min of agg capital grid
double precision, parameter :: kbarmax = 2.0 !max of agg capital grid
double precision, parameter :: nstdevz = 2.0 !number of st dev's to span for idio prod discretization
double precision, parameter :: rhoz = 0.859  !persistence of idio prod shock
double precision, parameter :: sigmaz = 0.022 !st dev of shock to idio prod
double precision, parameter :: nstdeva = 2.0 !number of st dev's to span for agg prod discretization
double precision, parameter :: rhoa = 0.859  !persistence of agg prod shock
double precision, parameter :: sigmaa=0.014 !st dev of shock to agg prod
double precision, parameter :: vferrortol=1e-4 !tolerance of VF iteration
double precision, parameter :: kprimeerrortol = 1e-4 !tolerance on policy convergence
double precision, parameter :: xistarerrortol = 1e-4 !tolerance on adjustment threshold convergence
double precision, parameter :: perrortol = 1e-4 !tolerance on mkt clearing/price
double precision, parameter :: brenttol = 1e-6 !tolerance on Brent optimization of adjusting capital
double precision, parameter :: fcsterrortol = 1e-3 !tolerance on fcst rule convergence
double precision, parameter :: fcstgain=0.5 !dampening parameter for fcst rule update
double precision, parameter :: plb=2.0 !lb for price clearing
double precision, parameter :: pub=2.55 !ub for price clearing
double precision, parameter :: maxextrap = 5.0 !max percentage extrapolation in A

!some stuff to be available globally
double precision :: k0(knum),z0(znum),a0(anum),kbar0(kbarnum),kdense0(kdensenum),pr_mat_z(znum,znum),pr_mat_a(anum,anum),&
					kfcstmat(numbeta),pfcstmat(numbeta),V(znum,anum,knum,kbarnum),&
					Vold(znum,anum,knum,kbarnum),Vna(znum,anum,knum,kbarnum),Va(znum,anum,knum,kbarnum),&
					kprime(znum,anum,knum,kbarnum),xistar(znum,anum,knum,kbarnum),V2old(znum,anum,knum,kbarnum),&
					kprime2(znum,anum,knum,kbarnum),xistar2(znum,anum,knum,kbarnum),kprimeold(znum,anum,knum,kbarnum),&
					xistarold(znum,anum,knum,kbarnum),distkz(znum,kdensenum,numper),kprimep(znum,kdensenum),&
					xistarp(znum,kdensenum),padjustp(znum,kdensenum),perrorsim(numper),kbarfcstsim(numper),&
					kbarsim(numper),psim(numper),asimshock(numper),pfcstsim(numper),ysim(numper),isim(numper),&
					kfcststore(numbeta,maxfcstit),pfcststore(numbeta,maxfcstit),kfcstmatnew(numbeta),pfcstmatnew(numbeta),&
					ergz0(znum),ergz0old(znum),pr_mat_z_dum(znum,znum),pstore(anum,kbarnum),Kprimestore(anum,kbarnum),&
					perrorstore(anum,kbarnum),Kbarnoagg(anum),pnoagg(anum),Kbaraggrule(numper),paggrule(numper),Nsim(numper),&
					MICROsim(nummicro,numper),kfcsterror,pfcsterror,XOLS(numper-numdiscard,4),pYOLS(numper-numdiscard),&
					KYOLS(numper-numdiscard),WORK(2*4*(numper-numdiscard)),asim(numper)

integer :: asimpos(numper),kbarfcstinds(anum,kbarnum),info
integer, allocatable :: seedarray(:)

!stuff to be available to subfunctions and which may need to be varied across OpenMP threads
double precision :: pval,wval,aval,weight,kbarfcstval
integer :: kbarfcstind,RHSact,RHSzct
!$omp threadprivate(pval,wval,aval,kbarfcstind,weight,RHSact,RHSzct)
    
end module modparams

program kt_ks_cont
use omp_lib
use base_lib
use modparams
implicit none

integer :: zct,act,debugind,kct,zprimect,aprimect,kbarct,statect,stateprimect,&
    vfct,ct,zvalpos,avalpos,kvalpos,pfcstvalpos,kbarfcstvalpos,zctpos,&
    actpos,t,piter,kprimeind,kprimeindnoadj,fcstiter,perct,&
    adjustbias,accelit,pinitflag,numstates,seeddim,shockct,simct,ind,&
    kbarfcstindstore,momct
       
double precision :: start,finish,vferror,zval,kval,kbarval,kprimeval,&
    pfcstval,Vnaval,Vaval,wfcstval,Vnextval,xival,Vval,kprimeerror,xistarerror,&
    Vnextval1,Vnextval2,Yactualp,Iactualp,Cactualp,yval,perror,Kprimeactualp,&
    kprimevalnoadj,weightnoadj,xmean,ymean,x2mean,xymean,xval,&
    ergerror,kfcstbias,pfcstbias,p1val,p2val,p1error,p2error,shockprobIRF,&
    shockprobdenom,Nactualp,wgt,pvalstore,wvalstore,avalstore,weightstore,&
    fa,fb,fc,pvala,pvalb,pvalc,ival,padjust,EVval,kbarfcstvalstore

start = omp_get_wtime()

open(13,file="kt_ks_cont.txt")

!$omp parallel
write(*,*) "parallel hello to you."
!$omp end parallel

!!!!!!! INSERT PARAMETERS
numstates = znum*anum*knum*kbarnum

!write constants
open(8,file="constants.txt")
write(8,*) xibar,delta,numper,numdiscard
close(8)

!capital grids
call linspace(k0,log(kmin),log(kmax),knum); k0 = exp(k0);
call linspace(kbar0,kbarmin,kbarmax,kbarnum);
call linspace(kdense0,log(kmin),log(kmax),kdensenum); kdense0=exp(kdense0);

!fcst rule initialization
kfcstmat = (/  0.0874, 0.8282, 0.4868, 0.4469 /)
kfcststore(:,:)=0.0; kfcststore(:,1)=kfcstmat;

pfcstmat = (/ 0.9898, -0.3937, -0.6308, 0.7175 /)
pfcststore(:,:)=0.0; pfcststore(:,1)=pfcstmat;

!discretize exogenous processes A, z, and simulate A
call discretize_simulate()

!HAVE THIS IN HERE BECAUSE YOU INITIALLY HAD WRONG UPDATE
open(8,file="kbarsim.txt")
open(9,file="psim.txt")
do t=1,numper
read(8,*) kbarsim(t)
read(9,*) psim(t)
end do !t
close(8)
close(9)
call update_fcst()
kfcstmat = kfcstmat + fcstgain * (kfcstmatnew - kfcstmat)
pfcstmat = pfcstmat + fcstgain * (pfcstmatnew - pfcstmat)
!REMOVE THIS PORTION OF CODE ONCE YOU'VE GOT A SOLUTION

!intialize the value function
Vold(:,:,:,:) = 0.0
V2old(:,:,:,:) = 0.0
kprimeold(:,:,:,:) = 0.0
xistarold(:,:,:,:) = 0.5

if (doVFI==0) then ;!read in the Vold and V2old matrices from data files
open(8,file="kprimeold.txt")
do zct=1,znum
do act=1,anum
do kct=1,knum
do kbarct=1,knum
read(8,*) kprimeold(zct,act,kct,kbarct)
end do
end do
end do
end do
close(8)
end if

!initialize kprime for Howard acceleration
if (doVFI==1) then
!$omp parallel private(zct,act,kct,kbarct)
!$omp do collapse(4)
do zct=1,znum
do act=1,anum
do kct=1,knum
do kbarct=1,kbarnum
 	   
    kprimeold(zct,act,kct,kbarct) = k0(kct)

end do !kbarct
end do !kct
end do !act
end do !zct
!$omp end do nowait
!$omp end parallel
end if

!!!INSERT FCST RULE LOOPING APPARATUS HERE
do fcstiter=1,maxfcstit

!!!!!!!!!!!!!!!!! GIVEN A FCST RULE, PERFORM VFI
do vfct=1,maxvfit
    
    !note that with howard iteration, it is kprime that's important, not Vold
    Vold(:,:,:,:) = 0.0; V2old(:,:,:,:) = 0.0;
    
    !here is howard acceleration step
    do accelit=1,maxaccelit
    
      do zct=1,znum
        do act=1,anum
        do kct=1,knum
        do kbarct=1,kbarnum
            
            !determine states
            zval = z0(zct); aval = a0(act); kval = k0(kct); kbarval = kbar0(kbarct);
            
            !determine fcsts
            kbarfcstval = fkfcst(aval,kbarval)
            pval = fpfcst(aval,kbarval)
            wval = phi/pval;
            
            !! determine Vna(z,a,k,K)
            Vnaval = pval * freduced(zval,kval)
            kprimeval = max((1.0-delta)*kval,kmin)
                                   
            !add expectations
            EVval = EVfunc(zct,kprimeval,aval,kbarfcstval)
            Vnaval = Vnaval + beta * EVval
            
            !! determine Va(z,a,k,K)
            kprimeval = kprimeold(zct,act,kct,kbarct)
            Vaval = pval *  freduced(zval,kval)
            Vaval = Vaval - pval*(kprimeval - (1.0-delta)*kval)

            !add expectations
            EVval = EVfunc(zct,kprimeval,aval,kbarfcstval)
            Vaval = Vaval + beta * EVval
			
			!compute implications for xi and V
            xival = (Vaval - Vnaval)/phi
            Vval = -1.0 * phi * expecxi(xival)&
            + cdfxi(xival) * Vaval + (1.0 - cdfxi(xival)) * Vnaval;      
            V(zct,act,kct,kbarct) = Vval

            
        end do
        end do
        end do
        end do
        Vold = V

        !given Vold from last iteration, spline the function along k dimension
        do zct=1,znum
        do act=1,anum
        do kbarct=1,kbarnum
            call spline(k0,Vold(zct,act,:,kbarct),knum,dble(1.0e30),dble(1.0e30),V2old(zct,act,:,kbarct))
        end do !kbarct
        end do !act
        end do !zct
        
    end do !accelct    

    do zct=1,znum
	do act=1,anum
	do kct=1,knum
	do kbarct=1,kbarnum
       
        !determine states
        zval = z0(zct); aval = a0(act); kval = k0(kct); kbarval = kbar0(kbarct);
        
        !determine fcsts
		kbarfcstval = fkfcst(aval,kbarval)
		pval = fpfcst(aval,kbarval)
        wval = phi/pval;
        
        !! determine Vna(z,a,k,K)
        Vnaval = pval * freduced(zval,kval)
        kprimeval = max((1.0-delta)*kval,kmin)
        EVval = EVfunc(zct,kprimeval,aval,kbarfcstval)
	    Vnaval = Vnaval + beta * EVval
        Vna(zct,act,kct,kbarct) = Vnaval
        !end of block computing Vnaval
      
        !now, this block determines Va(z,a,k,K)
        !first, you need to optimize the function wrt kprimeval
        RHSzct =zct
        RHSact = act
        Vaval = brent(k0(1),k0(knum/2),k0(knum),fnRHS,brenttol,kprimeval) !note this sets correct kprimeval, not correct vaval
        Vaval = -1.0 * Vaval
      
        !at this point, Vaval = -p k' + beta E V'
        !now, construct current period return remaining portions
        Vaval = Vaval + pval * (freduced(zval,kval)+(1.0-delta)*kval)
        Va(zct,act,kct,kbarct) = Vaval
        kprime(zct,act,kct,kbarct) = kprimeval
        
        xival = (Vaval - Vnaval)/phi
        xistar(zct,act,kct,kbarct) = xival
        !!end of block computuing Vaval
        
        !now, process the info to compute the value function
        Vval = -1.0 * phi* expecxi(xival)&
            + cdfxi(xival) * Vaval + (1.0 - cdfxi(xival)) * Vnaval;
        V(zct,act,kct,kbarct) = Vval
        
    end do !kbarct
    end do !kct
    end do !act
    end do !zct

    vferror = maxval(abs((log(V)-log(Vold))))    
    kprimeerror = maxval(abs(log(kprime)-log(kprimeold)))
    xistarerror = maxval(abs(log(xistar)-log(xistarold)))
	
    if (mod(vfct,1)==0) then
    write(13,*) "VF iter = ",vfct,"VF error = ",vferror
    write(13,*) "VF iter = ",vfct,"Kprime error = ",kprimeerror
    write(13,*) "VF iter = ",vfct,"Xistar error = ",xistarerror
    write(13,*) " "
    
    write(*,*) "VF iter = ",vfct,"VF error = ",vferror
    write(*,*) "VF iter = ",vfct,"Kprime error = ",kprimeerror
    write(*,*) "VF iter = ",vfct,"Xistar error = ",xistarerror
    write(*,*) " "

    end if
    if (doVFI==0.and.fcstiter==1) exit; !in this case, only need one run through to initialize the V2 and param matrices
    if (kprimeerror<kprimeerrortol .and. xistarerror<xistarerrortol) exit
    Vold=V
    kprimeold = kprime
    xistarold = xistar
	
end do !vfct

finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished VFI in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished VFI in ",finish-start," seconds."


if (minval(kprime)<1.05*k0(1)) then
    write(13,*) "Lowest idio capital close to bottom of grid."
    write(*,*) "Lowest idio capital close to bottom of grid."
end if
if (maxval(kprime)>0.95*k0(knum)) then
    write(13,*) "Highest idio capital close to top of grid."
    write(*,*) "Highest idio capital close to top of grid."
end if

open(8,file="k0.txt")
do kct=1,knum; write(8,*) k0(kct); end do
close(8)

open(8,file="kbar0.txt")
do kbarct=1,kbarnum; write(8,*) kbar0(kbarct); end do
close(8)

open(8,file="V.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) V(zct,act,kct,kbarct)
end do; end do;end do; end do
close(8)

open(8,file="Vna.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) Vna(zct,act,kct,kbarct)
end do; end do;end do; end do
close(8)

open(8,file="Va.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) Va(zct,act,kct,kbarct)
end do; end do;end do; end do
close(8)

open(8,file="xistar.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) xistar(zct,act,kct,kbarct)
end do; end do;end do; end do
close(8)

open(8,file="kprime.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) kprime(zct,act,kct,kbarct)
end do; end do;end do; end do
close(8)

open(8,file="padjust.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) cdfxi(xistar(zct,act,kct,kbarct))
end do; end do;end do; end do
close(8)
!!!!END OF VF BLOCK

!NOW THAT THE VFI IS DONE, SIMULATE THE MODEL
distkz(:,:,:) = 0.0

!initialize
kbarsim(1) = 1.55
distkz(:,:,1) = 1.0
distkz(:,:,1) = distkz(:,:,1) / sum(distkz(:,:,1))
kbaraggrule(:) = 0.0
paggrule(:) = 0.0

!MICROsim(nummicro,numper)
MICROsim(:,:) = 0.0
! 1 = i/k
! 2 = stdev(i/k)
! 3 = P(inaction)
! 4 = P(i/k>=0.2)
! 5 = P(i/k<=-0.2)
! 6 = P(i/k > 0)
! 7 = P(i/k < 0)

do t=1,numper-1

    kprimep(:,:) = 0.0
    xistarp(:,:) = 0.0
    padjustp(:,:) = 0.0
    
    aval = asim(t); 
    kbarval = kbarsim(t)

    kbarfcstval = fkfcst(aval,kbarval)
    pfcstval = fpfcst(aval,kbarval)
	
    !note that we will in future iterate over p
    pvala = plb
    pvalb = pub
    pvalc = pvala + dble(0.67) * (pvalb-pvala) 
  
    do piter=1,maxpit
    
    !brent method for price optimization    
    if (piter==1) pval = pvala
    if (piter==2) pval = pvalb
    if (piter==3) pval = pvalc
    if (piter>3) then
            
	    !first, try inverse quadratic interpolation of the excess demand function
	    pval = ( pvala * fb * fc ) / ( (fa - fb) * (fa - fc) ) &
	        + ( pvalb * fa * fc ) / ( (fb - fa) * (fb - fc ) ) &
	        + ( pvalc * fa * fb ) / ( (fc - fa) * (fc - fb ) )

	    !if it lies within bounds, and isn't too close to the bounds, then done
	    !o/w, take bisection step
	    if ((minval( (/ abs(pvala - pval), abs(pvalb-pval) /) )<&
	            abs( (pvalb-pvala)/dble(9.0) ) ).or.(pval<pvala).or.(pval>pvalb))   then
	        pval = (pvala + pvalb) / dble(2.0)
	    end if
		
    end if
    
    wval = phi/pval;

    do zct=1,znum
	do kct=1,kdensenum          
        if (distkz(zct,kct,t)>0.0) then
        
        zval = z0(zct)
        kval = kdense0(kct)
        
        !block to determine Vnaval
        Vnaval = pval * freduced(zval,kval)
        kprimeval = max(kmin,(1.0-delta)*kval)
        Vnaval = Vnaval + beta * EVfunc(zct,kprimeval,aval,kbarfcstval)    
        
        !block to construct Vaval
        RHSzct =zct
        Vaval = brent(k0(1),k0(knum/2),k0(knum),fnRHS,brenttol,kprimeval) !note this sets correct kprimeval, not correct vaval
        Vaval = -1.0 * Vaval
        !at this point, Vaval = -p k' + beta E V'
        !now, construct current period return remaining portions
        Vaval = Vaval + pval * (freduced(zval,kval)+(1.0-delta)*kval)
        
        !now process this info
        xival = (Vaval - Vnaval)/phi
        Vval = -1.0 * phi * expecxi(xival) + cdfxi(xival)*Vaval &
            + (1.0-cdfxi(xival))*Vnaval
            
        kprimep(zct,kct)=kprimeval
        xistarp(zct,kct)=xival
        padjustp(zct,kct) = cdfxi(xival)
		
    end if !distkz tolerance
	end do !kct
	end do !zct

    Yactualp = 0.0
    Iactualp = 0.0
    Kprimeactualp = 0.0
    Cactualp = 0.0
    Nactualp = 0.0
    
    do zct=1,znum
    do kct=1,kdensenum
         if (distkz(zct,kct,t)>0.0) then
        zval = z0(zct); 
        kval = kdense0(kct)
        kprimeval=kprimep(zct,kct)
        xival=xistarp(zct,kct)
        
        yval = yreduced(zval,kval)
        
        Yactualp = Yactualp  + distkz(zct,kct,t) * yval
        Iactualp = Iactualp + distkz(zct,kct,t) * cdfxi(xival) &
            * (kprimeval - (1.0-delta)*kval)
    
        Kprimeactualp = Kprimeactualp + distkz(zct,kct,t) * cdfxi(xival) * kprimeval &
            +  distkz(zct,kct,t) * (1.0 - cdfxi(xival)) * (1.0-delta)*kval
            
        Nactualp = Nactualp +  distkz(zct,kct,t) * (nreduced(zval,kval) + expecxi(xival))              
            
        end if
    end do !kct
    end do  !zct
    Cactualp = Yactualp - Iactualp

    !are you initializing the brent?
    if (piter==1) fa = (1.0/pval) - Cactualp
    if (piter==2) fb = (1.0/pval) - Cactualp
    if (piter==3) fc = (1.0/pval) - Cactualp
    perror = 1.0/pval - Cactualp
    
    !if not restarting or initializing
   if (piter>3) then 
        if (perror<0) then
            pvalc = pvalb; fc = fb;
            pvalb = pval; fb = perror;
            !pval a doesn't change
        else if (perror>=0) then
            pvalc = pvala; fc = fa;
            pvala = pval; fa = perror;
            !pval b doesn't change
        end if
    end if

    if (abs(perror)<perrortol) exit

    end do !piter
	
    write(13,"(A,I5,A,I5,A,F7.5,A,F7.5,A,F7.5)") "t= ",t,", pit = ",piter,", err = ",perror,", p = ",pval,", aval = ",aval
	write(*,"(A,I5,A,I5,A,F7.5,A,F7.5,A,F7.5)") "t= ",t,", pit = ",piter,", err = ",perror,", p = ",pval,", aval = ",aval
    
    !now, the value of p has been determined, in a market clearing fashion
    kbarsim(t+1)=Kprimeactualp
    psim(t) = pval
    perrorsim(t) = perror
    kbarfcstsim(t+1) = kbarfcstval
    pfcstsim(t) = pfcstval
    ysim(t) = Yactualp
    isim(t) = Iactualp
    Nsim(t) = Nactualp
    
    !DH stat series calculations
    if (t==(numdiscard)) then
        
        kbaraggrule(t) = kbarsim(t)
        paggrule(t) = pval
        kbaraggrule(t+1)=kbarfcstval
    
    else if (t>numdiscard) then
		
        kbaraggrule(t+1)=fkfcst(aval,kbaraggrule(t))
        paggrule(t) = fpfcst(aval,kbaraggrule(t))
        
    end if
    
    !now that market clearing has occured, insert distribution into next period's
    do zct=1,znum
    do kct=1,kdensenum
        if (distkz(zct,kct,t)>0.0) then
        zval = z0(zct);
        kval = kdense0(kct);
        kprimeval = kprimep(zct,kct)
        kprimevalnoadj = max(kmin,(1.0-delta)*kval)
        xival = xistarp(zct,kct)
        
        !bracket kprimeval (which is kprime given that you adjust)
        kprimeind=kct
        call hunt(kdense0,kdensenum,kprimeval,kprimeind)
        
        if (kprimeind<=0) then
            kprimeind = 1; weight = 0.0
        else if (kprimeind>=1.and.kprimeind<=(kdensenum-1)) then
            weight = (kprimeval - kdense0(kprimeind))/(kdense0(kprimeind+1) - kdense0(kprimeind));
        else if (kprimeind>=kdensenum) then
            kprimeind = kdensenum-1; weight = 1.0
        end if
        
        !bracket kprimevalnoadj (which is kprime given that you don't adjust)
        kprimeindnoadj=kct
        call hunt(kdense0,kdensenum,kprimevalnoadj,kprimeindnoadj)
        
        if (kprimeindnoadj<=0) then
            kprimeindnoadj = 1; weightnoadj = 0.0
        else if (kprimeindnoadj>=1.and.kprimeindnoadj<=(kdensenum-1)) then
            weightnoadj = (kprimevalnoadj - kdense0(kprimeindnoadj))/(kdense0(kprimeindnoadj+1) - kdense0(kprimeindnoadj));
        else if (kprimeindnoadj>=kdensenum) then
            kprimeindnoadj = kdensenum-1; weightnoadj = 1.0
        end if
        
        do zprimect=1,znum
            
            !transfer the weight to kprimeind that does adjust
            distkz(zprimect,kprimeind,t+1) = distkz(zprimect,kprimeind,t+1) &
                + pr_mat_z(zct,zprimect) * padjustp(zct,kct) * (1.0 - weight) * distkz(zct,kct,t)
            
            !transfer the weight to kprimeind+1 that does adjust
            distkz(zprimect,kprimeind+1,t+1) = distkz(zprimect,kprimeind+1,t+1) &
                + pr_mat_z(zct,zprimect) * padjustp(zct,kct) * weight * distkz(zct,kct,t)
    
            !transfer the weight to kprimeindnoadj that doesn't adjust
            distkz(zprimect,kprimeindnoadj,t+1) = distkz(zprimect,kprimeindnoadj,t+1) &
                + pr_mat_z(zct,zprimect) * (1.0 - padjustp(zct,kct)) * (1.0 - weightnoadj) * distkz(zct,kct,t)
            
            
            !transfer the weight to kprimeindnoadj+1 that doesn't adjust
            distkz(zprimect,kprimeindnoadj+1,t+1) = distkz(zprimect,kprimeindnoadj+1,t+1) &
                + pr_mat_z(zct,zprimect) * (1.0 - padjustp(zct,kct)) * weightnoadj * distkz(zct,kct,t)
            
        end do !zprimect
        end if !distkz tolerance
    end do !kct
    end do !zct
    
    !now that other simulation apparatus is complete, input correct data into moment storage
    !MICROsim(nummicro,numper)
    ! 1 = i/k
    ! 2 = stdev(i/k)
    ! 3 = P(inaction)
    ! 4 = P(i/k>=0.2)
    ! 5 = P(i/k<=-0.2)
    ! 6 = P(i/k > 0)
    ! 7 = P(i/k < 0)
    
    do zct=1,znum
    do kct=1,kdensenum
        if (distkz(zct,kct,t)>0.0) then
            zval = z0(zct);
            kval = kdense0(kct);
            kprimeval = kprimep(zct,kct)
            ival = kprimeval - (1.0-delta)*kval !investment conditional upon investment
            kprimevalnoadj = max(kmin,(1.0-delta)*kval)
            xival = xistarp(zct,kct)
            padjust = cdfxi(xival)
            
            !investment rate
            MICROsim(1,t) = MICROsim(1,t) + distkz(zct,kct,t) * padjust * (ival / kval)
            
            !investment rate squared - for stdev construction
            MICROsim(2,t) = MICROsim(2,t) + distkz(zct,kct,t) * padjust * ( (ival / kval)**2.0 )
            
            !P(inaction)
            MICROsim(3,t) = MICROsim(3,t) + distkz(zct,kct,t) * (1.0 - padjust)
            
            !P(pos. spike)
            if ((ival/kval)>=0.2) then
                MICROsim(4,t) = MICROsim(4,t) + distkz(zct,kct,t) * padjust 
            end if
            
            !P(neg. spike)
            if ((ival/kval)<=-0.2) then
                MICROsim(5,t) = MICROsim(5,t) + distkz(zct,kct,t) * padjust 
            end if
            
            !P(pos. invest)
            if ((ival/kval)>0.0) then
                MICROsim(6,t) = MICROsim(6,t) + distkz(zct,kct,t) * padjust 
            end if
            
            !P(neg. invest)
            if ((ival/kval)<0.0) then
                MICROsim(7,t) = MICROsim(7,t) + distkz(zct,kct,t) * padjust 
            end if
            
        end if
    end do !kct
    end do !zct

    !now, convert squared investment moment to stdev
    MICROsim(2,t) = sqrt(MICROsim(2,t) - (MICROsim(1,t)**2.0))
    
end do !t

finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished simulation in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished simulation in ",finish-start," seconds."

!!!NOW THAT SIMULATION IS DONE, DO THE FCST RULE UPDATE

!now, from the simulated values of p and K', update the fcst arrays
call update_fcst()

if (doVFI==1) then; !store the arrays for future use
open(8,file="kprimeold.txt")
do zct=1,znum
do act=1,anum
do kct=1,knum
do kbarct=1,knum
write(8,*) kprimeold(zct,act,kct,kbarct)
end do
end do
end do
end do
close(8)
end if

open(8,file="kdense0.txt")
do kct=1,kdensenum; write(8,*) kdense0(kct); end do
close(8)

open(8,file="psim.txt")
do t=1,numper; write(8,*) psim(t); end do !t
close(8)

open(8,file="kbarsim.txt")
do t=1,numper; write(8,*) kbarsim(t); end do !t
close(8)

open(8,file="perrorsim.txt")
do t=1,numper; write(8,*) perrorsim(t); end do !t
close(8)

open(8,file="kbarfcstsim.txt")
do t=1,numper; write(8,*) kbarfcstsim(t); end do !t
close(8)

open(8,file="pfcstsim.txt")
do t=1,numper; write(8,*) pfcstsim(t); end do !t
close(8)

open(8,file="distkzsim.txt");
do zct=1,znum; do kct=1,kdensenum; do t=1,numper
write(8,*) distkz(zct,kct,t)
end do; end do; end do
close(8)

open(8,file="ysim.txt")
do t=1,numper; write(8,*) ysim(t); end do !t
close(8)

open(8,file="isim.txt")
do t=1,numper; write(8,*) isim(t); end do !t
close(8)

open(8,file="Nsim.txt")
do t=1,numper; write(8,*) Nsim(t); end do !t
close(8)

open(8,file="paggrule.txt")
do t=1,numper; write(8,*) paggrule(t); end do !t
close(8)

open(8,file="kbaraggrule.txt")
do t=1,numper; write(8,*) kbaraggrule(t); end do !t
close(8)

open(8,file="kprimep.txt")
do zct=1,znum; do kct=1,kdensenum; write(8,*) kprimep(zct,kct); end do; end do 
close(8)

open(8,file="xistarp.txt")
do zct=1,znum; do kct=1,kdensenum; write(8,*) xistarp(zct,kct); end do; end do 
close(8)

open(8,file="padjustp.txt")
do zct=1,znum; do kct=1,kdensenum; write(8,*) padjustp(zct,kct); end do; end do 
close(8)

open(8,file="MICROsim.txt")
do t=1,numper
do momct=1,nummicro
write(8,*) MICROsim(momct,t)
end do !momct
end do !t
close(8)

open(8,file="kfcststore.txt")
do ct=1,maxfcstit; write(8,*) kfcststore(:,ct); end do; 
close(8)

open(8,file="pfcststore.txt")
do ct=1,maxfcstit; write(8,*) pfcststore(:,ct); end do; 
close(8)

write(*,*) " "
if (doVFI==0) exit
if (pfcsterror<fcsterrortol.and.kfcsterror<fcsterrortol) exit    

!if forecast rule not yet converged, update
kfcstmatnew = kfcstmat + fcstgain * (kfcstmatnew - kfcstmat)
kfcststore(:,fcstiter+1) = kfcstmatnew
kfcstmat = kfcstmatnew

pfcstmatnew = pfcstmat + fcstgain * (pfcstmatnew - pfcstmat)
pfcststore(:,fcstiter+1) = pfcstmatnew
pfcstmat = pfcstmatnew


write(*,*) " "
end do !fcstiter

write(*,*) "Done with GE at ",omp_get_wtime()-start," seconds."
write(13,*) "Done with GE at ",omp_get_wtime()-start," seconds."

write(*,*) "kfcstmat = ",kfcstmat

write(*,*) "pfcstmat = ",pfcstmat

finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished in ",finish-start," seconds."
close(13); !closing log file

contains

subroutine update_fcst()
implicit none

!takes simulated values of p, K', and runs OLS to get new coefficient matrices
!trying to go from simulated data to values of pfcstmatnew & kfcstmatnew
!data range: numdiscard,numper-1

integer :: M,N,t,ct

!do general setup

M = numper-numdiscard
N = numbeta

XOLS(:,:) = 0.0

ct = 0
do t=numdiscard,numper-1
	ct = ct +1 
	XOLS(ct,1) = dble(1.0)
	XOLS(ct,2) = log(kbarsim(t))
	XOLS(ct,3) = log(asim(t))
	XOLS(ct,4) = log(asim(t)) * log(kbarsim(t))
	
	KYOLS(ct) = log(kbarsim(t+1))
	
end do !t

!do kfcstmatnew
call dgels('N', M, N, 1, XOLS, M, KYOLS, M, WORK, 2*4*(numper-numdiscard),info)
kfcstmatnew = KYOLS(1:numbeta)

!then, do pfcstmatnew
!first, have to re-initialize some arrays because DGELS overwrites data with other stuff

WORK(:) = 0.0
ct = 0
do t=numdiscard,numper-1
	ct = ct +1 
	XOLS(ct,1) = dble(1.0)
	XOLS(ct,2) = log(kbarsim(t))
	XOLS(ct,3) = log(asim(t))
	XOLS(ct,4) = log(asim(t)) * log(kbarsim(t))
	
	pYOLS(ct) = log(psim(t))

end do !t

!call LAPACK to perform OLS 
call dgels('N', M, N, 1, XOLS, M, pYOLS, M, WORK, 2*4*(numper-numdiscard),info)
pfcstmatnew = pYOLS(1:numbeta)


write(13,*) " "; write(13,*) "kfcstmat (old, new) "
write(*,*) " "; write(*,*) "kfcstmat (old, new) "

write(13,*) kfcstmat
write(*,*) kfcstmat

write(13,*) kfcstmatnew
write(*,*) kfcstmatnew

write(13,*) " "
write(*,*) " "

write(13,*) "pfcstmat (old, new) "
write(*,*) "pfcstmat (old, new) "

write(13,*) pfcstmat
write(*,*) pfcstmat

write(13,*) pfcstmatnew
write(*,*) pfcstmatnew

write(13,*) " "
write(*,*) " "

kfcsterror = maxval(abs(kfcstmat-kfcstmatnew))
pfcsterror = maxval(abs(pfcstmat-pfcstmatnew))

write(13,*) " "; write(13,*) "Fcst Ct = ",fcstiter," K Fcst error = ",kfcsterror;
write(*,*) " "; write(*,*) "Fcst Ct = ",fcstiter," K Fcst error = ",kfcsterror;
write(13,*) "Fcst Ct = ",fcstiter," P Fcst error = ",pfcsterror; write(13,*) " "
write(*,*) "Fcst Ct = ",fcstiter," P Fcst error = ",pfcsterror; write(*,*) " "


end subroutine update_fcst

subroutine discretize_simulate()
implicit none

double precision :: asimgrid(anum)
integer :: act,zct,statect,ct,t,aprimect,shockct


!draw random seeds
call random_seed(size=seeddim)

!insert random seed into seedarray
allocate(seedarray(seeddim))
do ct=1,seeddim
    seedarray(ct) = seedint + ct
end do !ct
call random_seed(put=seedarray)

!BEGIN DISCRETIZATION PORTION
!discretize idio prod process
call tauchen(znum,rhoz,sigmaz,nstdevz,pr_mat_z,z0)

!discretize agg prod process
call tauchen(anum,rhoa,sigmaa,nstdeva,pr_mat_a,a0)

open(8,file="pr_mat_z.txt")
do zct=1,znum
write(8,*) pr_mat_z(zct,:)
end do !zct
close(8)

open(8,file="z0.txt")
do zct=1,znum
write(8,*) z0(zct)
end do !zct
close(8)

open(8,file="pr_mat_a.txt")
do act=1,anum
write(8,*) pr_mat_a(act,:)
end do !azct
close(8)

open(8,file="a0.txt")
do act=1,anum
write(8,*) a0(act)
end do !act
close(8)
!!!END DISCRETIZATION PORTION

asim(1) = ainit

do t=2,numper
  
  asim(t) = rhoa * log(asim(t-1)) + sigmaa * rnormdraw()
  asim(t) = exp(asim(t))
  
end do

open(8,file="asim.txt")
do t=1,numper
write(8,*) asim(t)
end do !t
close(8)

!END UNCONDITIONAL SIMULATION

end subroutine discretize_simulate

double precision function fkfcst(aggprod,kbarval)

integer :: ct
double precision :: aggprod,kbarval

fkfcst = 0.0

fkfcst = fkfcst + kfcstmat(1)
fkfcst = fkfcst + kfcstmat(2) * log(kbarval)
fkfcst = fkfcst + kfcstmat(3) * log(aggprod)
fkfcst = fkfcst + kfcstmat(4) * log(kbarval) * log(aggprod)

fkfcst = exp(fkfcst)

end function fkfcst

double precision function fpfcst(aggprod,kbarval)

integer :: ct
double precision :: aggprod,kbarval

fpfcst = 0.0

fpfcst = fpfcst + pfcstmat(1)
fpfcst = fpfcst + pfcstmat(2) * log(kbarval)
fpfcst = fpfcst + pfcstmat(3) * log(aggprod)
fpfcst = fpfcst + pfcstmat(4) * log(kbarval) * log(aggprod)

fpfcst = exp(fpfcst)

end function fpfcst

subroutine ghquad(aggprod,aprimenodes,aprimewgts)
implicit none

!this function takes as input the value of agg. productivity today (aggprod)
!and outputs Guass-Hermite quadrature weights and nodes for agg prod tomorrow

double precision, intent(in) :: aggprod
double precision, intent(out):: aprimenodes(anum),aprimewgts(anum)

!input 5-point Gauss-Hermite quadrature weights
aprimewgts = (/0.01125741,0.22207592,0.53333333,0.22207592,0.01125741/)
aprimewgts = aprimewgts / sum(aprimewgts)

!input grid for integration based on Gauss-Hermite nodes
aprimenodes = (/-2.020182870,-0.9585724646,0.0,0.9585724646,2.020182870/)
aprimenodes = sqrt(2.0)*sigmaa*aprimenodes + rhoa * log(aggprod)
aprimenodes = exp(aprimenodes)

end subroutine ghquad

double precision function EVfunc(zct,kprimeval,aggprod,kbarprimeval)
implicit none

!!!!!
!this function takes as input (zct,k',A,K') & evaluates
!E_{z',A'|z,A} V(z',k';A',K') with:
!
!Tauchen discretization of z --> z'
!Cubic spline evaluation in k'
!Gauss-Hermite quadrature in A'
!Bilinear interpolation in (A',K')
!!!!!

!input declarations
integer :: zct
double precision :: kprimeval,aggprod,kbarprimeval

!other declarations
integer :: zprimect,aprimect,kbarind,aprimeind
double precision :: kbarwgt,aprimewgt,ghnodes(anum),ghwgts(anum),aprimeval,Vnextval

	!first, get the quad nodes for integration of log(A') ~ N(rhoa * log(A), sigmaa^2) 
	call ghquad(aggprod,ghnodes,ghwgts)

	!second, bracket aggregate capital and get weights
	kbarind = kbarnum/2
	call hunt(kbar0,kbarnum,kbarprimeval,kbarind)
	if (kbarind<=0) then
	kbarind =  1
	else if (kbarind>=kbarnum) then
	kbarind = kbarnum-1
	end if
	kbarwgt = ( kbarprimeval - kbar0(kbarind) ) / ( kbar0(kbarind+1) - kbar0(kbarind) )
	kbarwgt = min(max(kbarwgt,-1.0*maxextrap),1.0+maxextrap)
	
	!initialize the sum of interest
	EVfunc = 0.0

	!now, loop over A'
	do aprimect=1,anum
		
		!extract agg prod next period
		aprimeval = ghnodes(aprimect)
		
		!now, bracket the aggregate productivity next period
		aprimeind = anum/2
		call hunt(a0,anum,aprimeval,aprimeind)
		if (aprimeind<=0) then
		aprimeind =  1
		else if (aprimeind>=anum) then
		aprimeind = anum-1
		end if
		aprimewgt = ( aprimeval - a0(aprimeind) ) / ( a0(aprimeind+1) - a0(aprimeind) )
		aprimewgt = min(max(aprimewgt,-1.0*maxextrap),1.0+maxextrap)
		

		!now, loop over z'
		do zprimect=1,znum

			!evaluate & increment for (1,1)
			call splint(k0,Vold(zprimect,aprimect,:,kbarind),V2old(zprimect,aprimect,:,kbarind),&
            knum,kprimeval,Vnextval)
		
			EVfunc = EVfunc + pr_mat_z(zct,zprimect) * ghwgts(aprimect) * ( &
				(1.0-aprimewgt)*(1.0-kbarwgt)*Vnextval )
			
			!evaluate & increment for (1,2)
			call splint(k0,Vold(zprimect,aprimect,:,kbarind+1),V2old(zprimect,aprimect,:,kbarind+1),&
            knum,kprimeval,Vnextval)
	
			EVfunc = EVfunc + pr_mat_z(zct,zprimect) * ghwgts(aprimect) * ( &
				(1.0-aprimewgt)*kbarwgt*Vnextval )
				
			!evaluate & increment for (2,1)
			call splint(k0,Vold(zprimect,aprimect+1,:,kbarind),V2old(zprimect,aprimect+1,:,kbarind),&
            knum,kprimeval,Vnextval)
	
			EVfunc = EVfunc + pr_mat_z(zct,zprimect) * ghwgts(aprimect) * ( &
				aprimewgt*(1.0-kbarwgt)*Vnextval )
				
			!evaluate & increment for (2,2)
			call splint(k0,Vold(zprimect,aprimect+1,:,kbarind+1),V2old(zprimect,aprimect+1,:,kbarind+1),&
            knum,kprimeval,Vnextval)
	
			EVfunc = EVfunc + pr_mat_z(zct,zprimect) * ghwgts(aprimect) * ( &
				aprimewgt*kbarwgt*Vnextval )
			

			
		end do !zprimect
	end do !aprimect
	!at this point, EVfunc = E_{z',A'|z,A}V(z',k';A',K')
end function EVfunc

double precision function cdfxi(xi)
implicit none
double precision :: xi
if (xi<0.0) cdfxi = 0.0
if (xi>0.0.and.xi<=xibar) cdfxi = xi/xibar
if (xi>xibar) cdfxi = 1.0
end function cdfxi
    
double precision function expecxi(xi)
implicit none
double precision :: xi
if(xi<0.0) expecxi = 0.0
if(xi<=xibar.and.xi>=0.0)expecxi = ( xi ** 2.0 ) / (2.0*xibar)
if(xi>xibar)expecxi = ( xibar ** 2.0 ) / (2.0*xibar)
end function expecxi

double precision function nreduced(z,k)
implicit none
double precision :: z,k

double precision :: exponentnu

exponentnu = 1.0 / (1.0 - nu)

nreduced = ( nu **  exponentnu ) * ( wval ** ( -1.0 * exponentnu ) )
nreduced = nreduced * (z ** exponentnu) * (aval ** exponentnu) * ( k ** (alpha * exponentnu) )
end function nreduced

double precision function freduced(z,k)
implicit none
double precision :: z,k

double precision :: exponentnu

exponentnu = 1.0 / (1.0 - nu)

freduced = ( nu ** (nu * exponentnu) ) * (1.0 - nu) * ( wval ** ( -1.0 * nu * exponentnu ) )
freduced = freduced * (z ** exponentnu) * (aval ** exponentnu) * ( k ** (alpha * exponentnu) )
end function freduced

double precision function yreduced(z,k)
implicit none
double precision :: z,k

double precision :: exponentnu

exponentnu = 1.0 / (1.0 - nu)

yreduced = ( (nu/wval) ** (nu * exponentnu) ) * ( (aval * z ) ** exponentnu  )
yreduced = yreduced * ( k ** ( alpha * exponentnu ) )

end function yreduced

double precision function fnRHS(kprimeval)
implicit none

!returns -1 (-pk' + beta E V(z',k',a',K'^F) ) given k'

double precision :: kprimeval,Vnextval1,Vnextval2,EnextVval
integer :: zprimect,aprimect,stateprimect,statect

!current return    
fnRHS = -pval * kprimeval

!expected future return
EnextVval = EVfunc(RHSzct,kprimeval,aval,kbarfcstval)

!set up the Bellman Eqn
fnRHS = fnRHS + beta * EnextVval

!you want to MINIMIZE this
fnRHS = -1.0 * fnRHS
end function fnRHS

double precision function rnormdraw()
implicit none

double precision, parameter :: PI = 3.141592653589793238462643383279 !define PI constant

!this function draws N(0,1) using Box-Muller tranform and 2 uniform random draws

double precision :: u1,u2

double precision :: theta,R

call random_number(u1)
call random_number(u2)

theta = 2.0 * PI * u1
R = sqrt(-2.0 * log(u2))

rnormdraw = R * sin(theta)

end function rnormdraw

end program kt_ks_cont