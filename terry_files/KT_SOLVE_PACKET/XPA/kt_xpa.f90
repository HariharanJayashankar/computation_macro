!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kt_xpa.f90
!
! Fortran code for the XPA solution of the Khan and Thomas (2008) 
! model, in the baseline version.
!
! 'Alternative Methods for Solving Heterogeneous Firm Models'
! Stephen Terry (2015)
!
! This Version : 12/21/15
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module modparams
use omp_lib
implicit none

integer, parameter :: knum = 10 !size of idio capital grid
integer, parameter :: kbarnum = 10 !size of agg capital grid
integer, parameter :: znum = 5;  !size of idio prod grid
integer, parameter :: anum = 5;  !size of agg prod grid
integer, parameter :: nummicro = 7 !number of micro moments to store
integer, parameter :: numstates = znum*anum*knum*kbarnum !total number of states
integer, parameter :: kdensenum = 50 !number of idio capital points on simulation grid
integer, parameter :: numper = 2500 !number of periods in unconditional simulation
integer, parameter :: numdiscard = 500 !number of periods to discard from unconditional simulation
integer, parameter :: seedint=2503 !seed for random shocks
integer, parameter :: ainit = 3 !initial agg prod state in simulations
integer, parameter :: numperIRF=50 !number of periods in each simulated IRF economy
integer, parameter :: numsimIRF=2000 !total number of simulated IRF economies
integer, parameter :: shockperIRF=25 !period at which to shock the economy in IRF
integer, parameter :: shockanumIRF=anum !agg prod destination upon shock
integer, parameter :: maxvfit = 1000 !max number of VF iterations
integer, parameter :: maxaccelit = 50 !number of Howard accelerations
integer, parameter :: maxpit=100 !max number of price iterations for clearing
integer, parameter :: maxfcstit=30 !max number of iterations on forecasting system
integer, parameter :: maxergit = 5000 !max number of iterations to determine ergodic dists of a,z
integer, parameter :: doVFI = 1 !do VFI or read in files from storage?
integer, parameter :: adjustbias = 1 !adjust for Jensen's inequality bias via SS solution?
integer, parameter :: doIRF = 1 !do IRF simulation?

double precision, parameter :: alpha = 0.256 !capital elasticity
double precision, parameter :: nu = 0.640 !labor elasticity
double precision, parameter :: phi = 2.4 !disutility of labor
double precision, parameter :: xibar = 0.0083 !upper bound of capital AC dist
double precision, parameter :: beta = 0.977 !discount rate
double precision, parameter :: delta = 0.069 !capital depreciation rate
double precision, parameter :: kmin= 0.1 !min of idio capital grid
double precision, parameter :: kmax = 8.0 !max of idio capital grid
double precision, parameter ::kbarmin =1.25 !min of agg capital grid
double precision, parameter ::kbarmax = 2.0 !max of agg capital grid
double precision, parameter ::nstdevz = 2.0 !number of st dev's to span for idio prod discretization
double precision, parameter ::rhoz = 0.859; !persistence of idio prod shock
double precision, parameter ::sigmaz = 0.022 !st dev of shock to idio prod
double precision, parameter ::nstdeva = 2.0 !number of st dev's to span for agg prod discretization
double precision, parameter ::rhoa = 0.859; !persistence of agg prod shock
double precision, parameter ::sigmaa=0.014 !st dev of shock to agg prod
double precision, parameter ::ergdistatol=1.0e-5 !tolerance on computing ergodic distribution
double precision, parameter ::shocksizeIRF=sigmaa !shock size in IRF simulations
double precision, parameter ::vferrortol=1.0e-4 !tolerance of VF iteration
double precision, parameter ::kprimeerrortol = 1.0e-4 !tolerance on policy convergence
double precision, parameter ::xistarerrortol = 1.0e-4 !tolerance on adjustment threshold convergence
double precision, parameter ::perrortol = 1.0e-4 !tolerance on mkt clearing/price
double precision, parameter ::brenttol = 1.0e-6 !tolerance on Brent optimization for capital
double precision, parameter ::fcsterrortol = 1.0e-3 !tolerance on fcst rule convergence
double precision, parameter ::fcstgain=0.5 !dampening parameter for fcst rule update 
double precision, parameter ::ergerrortol = 1.0e-10; !tolerance on finding error distribution

double precision :: k0(knum),z0(znum),a0(anum),kbar0(kbarnum),pr_mat(anum*znum,anum*znum),pr_mat_a(anum,anum),&
	pr_mat_z(znum,znum),V(znum,anum,knum,kbarnum),Vold(znum,anum,knum,kbarnum),&
	Vna(znum,anum,knum,kbarnum),Va(znum,anum,knum,kbarnum),kprime(znum,anum,knum,kbarnum),&
	xistar(znum,anum,knum,kbarnum),V2old(znum,anum,knum,kbarnum),kprime2(znum,anum,knum,kbarnum),&
	xistar2(znum,anum,knum,kbarnum),kprimeold(znum,anum,knum,kbarnum),&
	xistarold(znum,anum,knum,kbarnum),asimshock(numper),kbarsim(numper),psim(numper),&
	kdense0(kdensenum),distkz(znum,kdensenum,numper),kprimep(znum,kdensenum),xistarp(znum,kdensenum),&
	padjustp(znum,kdensenum),perrorsim(numper),kbarfcstsim(numper),pfcstsim(numper),ysim(numper),isim(numper),&
	kfcststore(anum,2,maxfcstit),pfcststore(anum,2,maxfcstit),Kss(anum),pss(anum),kfcstmat(anum,2),pfcstmat(anum,2),&
	kfcstmatnew(anum,2),pfcstmatnew(anum,2),asimgrid(anum),ergz0(znum),ergz0old(znum),&
	pr_mat_z_dum(znum,znum),pstore(anum,kbarnum),Kprimestore(anum,kbarnum),perrorstore(anum,kbarnum),&
	Kbarnoagg(anum),pnoagg(anum),Kbaraggrule(numper),paggrule(numper),Nsim(numper),ergdista(anum),ergdistaold(anum),&
	asimshockIRF(numperIRF,numsimIRF),arriveshockIRF(numsimIRF),kbarsimIRF(numperIRF,numsimIRF,2),&
	psimIRF(numperIRF,numsimIRF,2),perrorsimIRF(numperIRF,numsimIRF,2),kbarfcstsimIRF(numperIRF,numsimIRF,2),&
	pfcstsimIRF(numperIRF,numsimIRF,2),ysimIRF(numperIRF,numsimIRF,2),isimIRF(numperIRF,numsimIRF,2),&
	distIRF(znum,kdensenum,numperIRF,numsimIRF,2),NsimIRF(numperIRF,numsimIRF,2),MICROsim(nummicro,numper)

integer :: asimpos(numper),asimposIRF(numperIRF,numsimIRF,2),seeddim

integer, allocatable :: seedarray(:)

!stuff to be available to subfunctions and which may need to be varied across OpenMP threads
double precision :: pval,wval,aval,weight,RHSkval
integer :: kbarfcstind,RHSact,RHSzct
!$omp threadprivate(pval,wval,aval,kbarfcstind,weight,RHSact,RHSzct,RHSkval)

end module modparams

program kt_xpa
use omp_lib
use base_lib
use modparams
implicit none

integer :: debugind,zct,kct,act,zprimect,aprimect,kbarct,&
    statect,stateprimect,vfct,ct,zvalpos,avalpos,kvalpos,pfcstvalpos,&
    kbarfcstvalpos,zctpos,actpos,kbardum,t,&
    piter,kprimeind,kprimeindnoadj,fcstiter,perct,&
	ergit,accelit,&
    pinitflag,loopct,kbarst,kbarend,&
    simct,shockct,momct
    
double precision :: start,finish,vferror,zval,kval,kbarval,&
    kbarfcstval,kprimeval,Vnaval,Vaval,Vnextval,xival,Vval,&
    kprimeerror,xistarerror,Vnextval1,Vnextval2,&
    Yactualp,Iactualp,Cactualp,yval,perror,&
    Kprimeactualp,kprimevalnoadj,weightnoadj,xmean,ymean,x2mean,&
    xymean,kfcsterror,pfcsterror,xval,&
    ergerror,kfcstbias,pfcstbias,p1val,p2val,&
    p1error,p2error,mval,Nactualp,&
    shockprobdenom,shockprobIRF,ival,padjust,pfcstval

start = omp_get_wtime()


open(13,file="kt_xpa.txt")


write(*,*) "available procs = ",omp_get_num_procs()
write(*,*) "available threads = ",omp_get_max_threads()
write(*,*) "using threads = ",omp_get_num_threads()

!$omp parallel

write(*,*) "parallel hello to you."

!$omp end parallel

!!!!!!! INSERT PARAMETERS

!read in the no agg unc capital value for later use
open(8,file="Kbarnoagg.txt")
do act=1,anum
read(8,*) Kbarnoagg(act)
end do
close(8)

open(8,file="pnoagg.txt")
do act=1,anum
read(8,*) pnoagg(act)
end do
close(8)

!fcst rule initialization
kfcstmat(1,:) = (/   5.5703123796310788E-002, 0.81914773272391017 /)
kfcstmat(2,:) = (/   6.9665036089883670E-002, 0.81538518691927342 /)
kfcstmat(3,:) = (/   8.1809515956025064E-002, 0.81478097273136640 /)
kfcstmat(4,:) = (/   9.4412606856224934E-002, 0.81379804303600745 /)
kfcstmat(5,:) = (/  0.10996621006719427     , 0.80876644960499511 /)

kfcststore(:,:,:)=0.0; kfcststore(:,:,1)=kfcstmat;

pfcstmat(1,:) = (/  1.0156048668373361  ,   -0.40906585902130099 /)
pfcstmat(2,:) = (/ 0.99893966987725991  ,   -0.40622933253252513 /)
pfcstmat(3,:) = (/ 0.98065875147066461  ,   -0.40188307388436484 /)
pfcstmat(4,:) = (/ 0.96172566916304070  ,   -0.39616735766609229 /)
pfcstmat(5,:) = (/ 0.94552394727216549  ,   -0.39495042951610010 /)


pfcststore(:,:,:)=0.0; pfcststore(:,:,1)=pfcstmat;

!find the fixed point of the forecasting rules
do act=1,anum
    Kss(act) = kfcstmat(act,1)/(1.0 - kfcstmat(act,2))
    pss(act) = pfcstmat(act,1) + pfcstmat(act,2) * Kss(act)
    Kss(act) = exp(Kss(act))
    pss(act) = exp(pss(act))
end do !act


call discretize_simulate()


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
do kct=1,knum
do act=1,anum
do kbarct=1,kbarnum
 	   
    kprimeold(zct,act,kct,kbarct) = k0(kct)

end do
end do
end do
end do
!$omp end do nowait
!$omp end parallel
end if




!!!INSERT FCST RULE LOOPING APPARATUS HERE
do fcstiter=1,maxfcstit


!!!!!!!!!!!!!!!!! GIVEN A FCST RULE, PERFORM VFI

do vfct=1,maxvfit
    
    !note that with howard iteration, it is kprime that's important, not Vold
    Vold(:,:,:,:) = 0.0; V2old(:,:,:,:) = 0.0;
    !write(*,*) "made it here 0"
    !here is howard acceleration step
    do accelit=1,maxaccelit
        !write(*,*) "accelit = ",accelit
        !$omp parallel private(zct,act,kct,kbarct,zval,kval,kbarval,&
        !$omp& kbarfcstval,Vnaval,kprimeval,statect,&
        !$omp& zprimect,aprimect,stateprimect,Vnextval,Vaval,&
        !$omp& xival,Vval,Vnextval1,Vnextval2)
        
        !$omp do collapse(4)        
        do zct=1,znum
        do act=1,anum
        do kct=1,knum
        do kbarct=1,kbarnum
            
            !determine states
            zval = z0(zct); aval = a0(act); kval = k0(kct); kbarval = kbar0(kbarct);
            
            !determine fcsts
            kbarfcstval = kfcstmat(act,1) + kfcstmat(act,2)*log(kbarval); kbarfcstval = exp(kbarfcstval);
            pval = pfcstmat(act,1) + pfcstmat(act,2)*log(kbarval); pval = exp(pval);
            wval = phi/pval;
            
            !! determine Vna(z,a,k,K)
            Vnaval = pval * freduced(zval,kval)
            kprimeval = max((1.0-delta)*kval,kmin)
            statect=(act-1)*znum+zct
            do zprimect=1,znum
            do aprimect=1,anum
                stateprimect = (aprimect-1)*znum+zprimect
                
                !now, need to evaluate Vold at (zprimect,aprimect,kprimeval,kbarfcstval)
                kbarfcstind = kbarct
                call hunt(kbar0,kbarnum,kbarfcstval,kbarfcstind); !kbarfcstval is in interval kbarfcstind
                
                if (kbarfcstind<=0) then 
                    weight = 0.0
                    kbarfcstind=1
                else if (kbarfcstind>=1.and.kbarfcstind<=(kbarnum-1)) then
                    weight = (kbarfcstval - kbar0(kbarfcstind))/( kbar0(kbarfcstind+1)-kbar0(kbarfcstind) )
                else if (kbarfcstind>=kbarnum) then
                    weight = 1.0
                    kbarfcstind = kbarnum-1
                end if
                
                call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
                knum,kprimeval,Vnextval1)
                call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
                knum,kprimeval,Vnextval2)
                
                Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)
                
                
                !now that next period's value function is evaluated, add to expectation
                Vnaval = Vnaval + beta * pr_mat(statect,stateprimect) * Vnextval
            end do
            end do
            
            
            
            !! determine Va(z,a,k,K)
            kprimeval = kprimeold(zct,act,kct,kbarct)
            Vaval = pval *  freduced(zval,kval)
            Vaval = Vaval - pval*(kprimeval - (1.0-delta)*kval)
            
            statect=(act-1)*znum+zct
            do zprimect=1,znum
            do aprimect=1,anum
                stateprimect = (aprimect-1)*znum+zprimect
                
                !now, need to evaluate Vold at (zprimect,aprimect,kprimeval,kbarfcstval)
                kbarfcstind = kbarct
                call hunt(kbar0,kbarnum,kbarfcstval,kbarfcstind); !kbarfcstval is in interval kbarfcstind
                
                if (kbarfcstind<=0) then 
                    weight = 0.0
                    kbarfcstind=1
                else if (kbarfcstind>=1.and.kbarfcstind<=(kbarnum-1)) then
                    weight = (kbarfcstval - kbar0(kbarfcstind))/( kbar0(kbarfcstind+1)-kbar0(kbarfcstind) )
                else if (kbarfcstind>=kbarnum) then
                    weight = 1.0
                    kbarfcstind = kbarnum-1
                end if
                
                call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
                knum,kprimeval,Vnextval1)
                call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
                knum,kprimeval,Vnextval2)
                
                Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)
                
                
                !now that next period's value function is evaluated, add to expectation
                Vaval = Vaval + beta * pr_mat(statect,stateprimect) * Vnextval
            end do
            end do
            
            xival = (Vaval - Vnaval)/phi
            Vval = -1.0 * phi * expecxi(xival)&
            + cdfxi(xival) * Vaval + (1.0 - cdfxi(xival)) * Vnaval;
            
            V(zct,act,kct,kbarct) = Vval
            
        end do
        end do
        end do
        end do
        !$omp end do nowait
        !$omp end parallel
        
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
    
    
    
    
    !$omp parallel private(zct,act,kct,kbarct,zval,kval,kbarval,&
    !$omp& kbarfcstval,Vnaval,kprimeval,statect,&
    !$omp& zprimect,aprimect,stateprimect,Vnextval,Vaval,&
    !$omp& xival,Vval,Vnextval1,Vnextval2)
        
    !$omp do collapse(4)
    do zct=1,znum
    do act=1,anum
    do kct=1,knum
    do kbarct=1,kbarnum

        !determine states
        zval = z0(zct); aval = a0(act); kval = k0(kct); kbarval = kbar0(kbarct);
        
        !determine fcsts
        kbarfcstval = kfcstmat(act,1) + kfcstmat(act,2)*log(kbarval); kbarfcstval = exp(kbarfcstval);
        pval = pfcstmat(act,1) + pfcstmat(act,2)*log(kbarval); pval = exp(pval);
        wval = phi/pval;
        
        !! determine Vna(z,a,k,K)
        Vnaval = pval * freduced(zval,kval)
        kprimeval = max((1.0-delta)*kval,kmin)
        statect=(act-1)*znum+zct
        do zprimect=1,znum
        do aprimect=1,anum
            stateprimect = (aprimect-1)*znum+zprimect
            
            !now, need to evaluate Vold at (zprimect,aprimect,kprimeval,kbarfcstval)
            kbarfcstind = kbarct
            call hunt(kbar0,kbarnum,kbarfcstval,kbarfcstind); !kbarfcstval is in interval kbarfcstind
            
            if (kbarfcstind<=0) then 
                weight = 0.0
                kbarfcstind=1
            else if (kbarfcstind>=1.and.kbarfcstind<=(kbarnum-1)) then
                weight = (kbarfcstval - kbar0(kbarfcstind))/( kbar0(kbarfcstind+1)-kbar0(kbarfcstind) )
            else if (kbarfcstind>=kbarnum) then
                weight = 1.0
                kbarfcstind = kbarnum-1
            end if
            
            call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
            knum,kprimeval,Vnextval1)
            call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
            knum,kprimeval,Vnextval2)
            
            Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)
            
            
            !now that next period's value function is evaluated, add to expectation
            Vnaval = Vnaval + beta * pr_mat(statect,stateprimect) * Vnextval
        end do
        end do
        
        Vna(zct,act,kct,kbarct) = Vnaval
        !end of block computing Vnaval
        
        !now, this block determines Va(z,a,k,K)
     	RHSact = act
     	RHSzct = zct
     	RHSkval = kval
        Vaval = brent(k0(1),k0(knum/2),k0(knum),fVa,brenttol,kprimeval)
        Vaval = -1.0 * Vaval; !this is to account for Brent routine's minimization
        
        Va(zct,act,kct,kbarct) = Vaval
        kprime(zct,act,kct,kbarct) = kprimeval
        
        xival = (Vaval - Vnaval)/phi
        xistar(zct,act,kct,kbarct) = xival
        !!end of block computuing Vaval
        
        !now, process the info to compute the value function
        Vval = -1.0 * phi *  expecxi(xival)&
            + cdfxi(xival) * Vaval + (1.0 - cdfxi(xival)) * Vnaval;
  
        V(zct,act,kct,kbarct) = Vval
        
    end do !kbarct
    end do !kct
    end do !act
    end do !zct
    !$omp end do nowait
    !$omp end parallel


    vferror = maxval(abs( (log(V)-log(Vold) )) )    
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
    if (doVFI==0.and.fcstiter==1) exit; !in this case, only need one run through to initialize the V2 matrix
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

!!!!END OF VF BLOCK

!!!!NOW THAT VF IS DONE, DO THE AGGREGATION ON THE (A,K) GRID

do act=1,anum
do kbarct=1,kbarnum
    
    aval = a0(act)
    kval = kbar0(kbarct)
    kbarval = kval
    
    kbarfcstval = kfcstmat(act,1) + kfcstmat(act,2)*log(kbarval);
    kbarfcstval = exp(kbarfcstval)
    pval = pfcstmat(act,1)+pfcstmat(act,2)*log(kbarval); 
    pval = exp(pval);
    
    
    !bisection algorithm over price
    p1val = 1.0
    p2val = 4.0
    
    do piter=1,maxpit
        
    pval = (p1val+p2val)/2.0        
        
    wval = phi/pval
    Yactualp=0.0; Iactualp = 0.0; Kprimeactualp=0.0;
    
    !given pval, do the optimization    
    do zct=1,znum
        
        zval = z0(zct)
        
        !!!block to form Vna at the particular point (z,a,k,K)
        Vnaval = pval * freduced(zval,kval)
        kprimeval = max((1.0-delta)*kval,kmin)
        statect=(act-1)*znum+zct
        do zprimect=1,znum
        do aprimect=1,anum
            stateprimect = (aprimect-1)*znum+zprimect
            
            !now, need to evaluate Vold at (zprimect,aprimect,kprimeval,kbarfcstval)
            kbarfcstind = kbarnum/2
            call hunt(kbar0,kbarnum,kbarfcstval,kbarfcstind); !kbarfcstval is in interval kbarfcstind
            
            if (kbarfcstind<=0) then 
                weight = 0.0
                kbarfcstind=1
            else if (kbarfcstind>=1.and.kbarfcstind<=(kbarnum-1)) then
                weight = (kbarfcstval - kbar0(kbarfcstind))/( kbar0(kbarfcstind+1)-kbar0(kbarfcstind) )
            else if (kbarfcstind>=kbarnum) then
                weight = 1.0
                kbarfcstind = kbarnum-1
            end if
            
            call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
            knum,kprimeval,Vnextval1)
            call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
            knum,kprimeval,Vnextval2)
            
            Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)
            
            
            !now that next period's value function is evaluated, add to expectation
            Vnaval = Vnaval + beta * pr_mat(statect,stateprimect) * Vnextval
        end do
        end do
        !!!end of block computing Vnaval
        
        !!!block computing Vaval
        
        RHSact=act
        RHSzct=zct        
        RHSkval=kval
        Vaval = brent(k0(1),k0(knum/2),k0(knum),fVa,brenttol,kprimeval)
        Vaval = -1.0 * Vaval; !this is to account for Brent routine's minimization
        
        xival = (Vaval - Vnaval)/phi
        !!end of block computing Vaval
        
        !now, process these decisions
        
        yval = yreduced(zval,kval)
        
        Yactualp = Yactualp + ergz0(zct)*yval
        Iactualp = Iactualp + ergz0(zct)*( cdfxi(xival) * (kprimeval - (1.0 - delta)*kval) )
        Kprimeactualp = Kprimeactualp + ergz0(zct) *  cdfxi(xival)*kprimeval &
        + ergz0(zct) * (1.0 - cdfxi(xival))*(1.0-delta)*kval
        
    end do !zct
    
    !now, compute the implied consumption, etc, via the array kprimep and xistarp from above
   
    
    Cactualp = Yactualp - Iactualp
    perror = 1.0/pval - Cactualp
    !end of bisection loop over pval --> wval    
    if (abs(perror)<perrortol) exit
    if (perror<0.0) then
        p2val = pval
    else if (perror>0.0) then
        p1val = pval
    end if
    
    
    end do !piter
    
    Kprimestore(act,kbarct) = Kprimeactualp
    pstore(act,kbarct) = pval
    perrorstore(act,kbarct) = perror

    write(*,"(A,I3,A,I3,A,F10.5,A,F8.5)") "a = ",act,", kbar = ",kbarct,", perror = ",perror,", pval = ",pval
    
end do !kbarct
end do !act


!!!NOW THAT THE DATA (A,K) --> (P,K') IS OBTAINED, UPDATE THE FCST RULES


do act=1,anum
    
    !do kprime first
    xmean = 0.0
    ymean = 0.0
    x2mean = 0.0
    xymean = 0.0
    
    kbarst = 1
    kbarend = kbarnum
    
    do kbarct = kbarst,kbarend
        
        xval = log(kbar0(kbarct))
        yval = log(Kprimestore(act,kbarct))

        xmean = xmean + xval
        ymean = ymean + yval
        x2mean = x2mean + xval ** 2.0
        xymean = xymean + xval * yval

    end do !period
    
    xmean = xmean / dble(kbarend-kbarst+1)
    ymean = ymean / dble(kbarend-kbarst+1)
    x2mean = x2mean / dble(kbarend-kbarst+1)
    xymean = xymean / dble(kbarend-kbarst+1)

    kfcstmatnew(act,2) = ( xymean - xmean * ymean ) / ( x2mean - ( xmean ** 2.0 ) )
    kfcstmatnew(act,1) = ymean - kfcstmatnew(act,2) * xmean
    
    if (adjustbias==1) then
        !now, adjust for fcstbias
        kfcstbias = kfcstmatnew(act,1) + kfcstmatnew(act,2)* log(Kbarnoagg(act)) - log(Kbarnoagg(act))
        kfcstmatnew(act,1) = kfcstmatnew(act,1) - kfcstbias
    
    end if
    
    
    !then do price rule
    xmean = 0.0
    ymean = 0.0
    x2mean = 0.0
    xymean = 0.0
    do kbarct = kbarst,kbarend
        
        xval = log(kbar0(kbarct))
        yval = log(pstore(act,kbarct))

        xmean = xmean + xval
        ymean = ymean + yval
        x2mean = x2mean + xval ** 2.0
        xymean = xymean + xval * yval
    
    end do !period
    
    xmean = xmean / dble(kbarend-kbarst+1)
    ymean = ymean / dble(kbarend-kbarst+1)
    x2mean = x2mean / dble(kbarend-kbarst+1)
    xymean = xymean / dble(kbarend-kbarst+1)
    
    pfcstmatnew(act,2) = ( xymean - xmean * ymean ) / ( x2mean - ( xmean ** 2.0 ) )
    pfcstmatnew(act,1) = ymean - pfcstmatnew(act,2) * xmean
    
    if (adjustbias==1) then
        !now, adjust for fcstbias
        pfcstbias = pfcstmatnew(act,1) + pfcstmatnew(act,2)* log(Kbarnoagg(act)) - log(pnoagg(act))
        pfcstmatnew(act,1) = pfcstmatnew(act,1) - pfcstbias
    
    end if
end do !act




write(13,*) " "; write(13,*) "kfcstmat = "
write(*,*) " "; write(*,*) "kfcstmat = "
do act=1,anum
write(13,*) kfcstmat(act,:)
write(*,*) kfcstmat(act,:)
end do !act
write(13,*) " "
write(*,*) " "

write(13,*) "kfcstmatnew = "
write(*,*) "kfcstmatnew = "
do act=1,anum
write(13,*) kfcstmatnew(act,:)
write(*,*) kfcstmatnew(act,:)
end do !act
write(13,*) " "
write(*,*) " "

write(13,*) "pfcstmat = "
write(*,*) "pfcstmat = "
do act=1,anum
write(13,*) pfcstmat(act,:)
write(*,*) pfcstmat(act,:)
end do !act
write(13,*) " "
write(*,*) " "

write(13,*) "pfcstmatnew = "
write(*,*) "pfcstmatnew = "
do act=1,anum
write(13,*) pfcstmatnew(act,:)
write(*,*) pfcstmatnew(act,:)
end do !act
write(13,*) " "
write(*,*) " "

kfcsterror = maxval(abs(kfcstmat-kfcstmatnew))
pfcsterror = maxval(abs(pfcstmat-pfcstmatnew))


write(13,*) " "; write(13,*) "Fcst Ct = ",fcstiter," K Fcst error = ",kfcsterror;
write(*,*) " "; write(*,*) "Fcst Ct = ",fcstiter," K Fcst error = ",kfcsterror;
write(13,*) "Fcst Ct = ",fcstiter," P Fcst error = ",pfcsterror; write(13,*) " "
write(*,*) "Fcst Ct = ",fcstiter," P Fcst error = ",pfcsterror; write(*,*) " "

finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished aggregation in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished aggregation in ",finish-start," seconds."

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



open(8,file="ergz0.txt")
do zct=1,znum; write(8,*) ergz0(zct); end do
close(8)

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

!write constants
open(8,file="constants.txt")
write(8,*) xibar,delta,numper,numdiscard,numperIRF,numsimIRF,shockperIRF,doIRF
close(8)

open(8,file="Kprimestore.txt")
do act=1,anum; do kbarct=1,kbarnum
write(8,*) Kprimestore(act,kbarct)
end do; end do;
close(8)

open(8,file="pstore.txt")
do act=1,anum; do kbarct=1,kbarnum
write(8,*) pstore(act,kbarct)
end do; end do;
close(8)

open(8,file="perrorstore.txt")
do act=1,anum; do kbarct=1,kbarnum
write(8,*) perrorstore(act,kbarct)
end do; end do;
close(8)

write(*,*) " "
if (pfcsterror<fcsterrortol.and.kfcsterror<fcsterrortol) exit    
if (maxfcstit>1) then
kfcstmatnew = kfcstmat + fcstgain * (kfcstmatnew - kfcstmat)
kfcststore(:,:,fcstiter+1) = kfcstmatnew
kfcstmat = kfcstmatnew

pfcstmatnew = pfcstmat + fcstgain * (pfcstmatnew - pfcstmat)
pfcststore(:,:,fcstiter+1) = pfcstmatnew
pfcstmat = pfcstmatnew
end if




end do !fcstiter

finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished GE solution in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished GE solution in ",finish-start," seconds."

open(8,file="kfcststore.txt")
do fcstiter=1,maxfcstit; do act=1,anum; write(8,*) kfcststore(act,:,fcstiter); end do; end do
close(8)

open(8,file="pfcststore.txt")
do fcstiter=1,maxfcstit; do act=1,anum; write(8,*) pfcststore(act,:,fcstiter); end do; end do
close(8)


!NOW THAT THE FCST RULES ARE FOUND, SIMULATE THE MODEL

!initialize
kbarsim(1) = 1.55
MICROsim(:,:) = 0.0
distkz(:,:,:) = 0.0
distkz(:,:,1) = 1.0
distkz(:,:,1) = distkz(:,:,1) / sum(distkz(:,:,1))
kbaraggrule(:) = 0.0
paggrule(:) = 0.0

do t=1,numper-1

	!extract aggregate prod
	act = asimpos(t)
	aval = a0(act)
	RHSact=act
	
	!extract agg capital & determine grid point/weight of next period forecast
	kbarval = kbarsim(t)
	
	kbarfcstval = kfcstmat(act,1) + kfcstmat(act,2)*log(kbarval);
    kbarfcstval = exp(kbarfcstval)
    
    pfcstval = pfcstmat(act,1) + pfcstmat(act,2)*log(kbarval);
    pfcstval = exp(pfcstval);
    
	kbarfcstind = kbarnum/2
	call hunt(kbar0,kbarnum,kbarfcstval,kbarfcstind); !kbarfcstval is in interval kbarfcstind

	if (kbarfcstind<=0) then 
		weight = 0.0
		kbarfcstind=1
	else if (kbarfcstind>=1.and.kbarfcstind<=(kbarnum-1)) then
		weight = (kbarfcstval - kbar0(kbarfcstind))/( kbar0(kbarfcstind+1)-kbar0(kbarfcstind) )
	else if (kbarfcstind>=kbarnum) then
		weight = 1.0
		kbarfcstind = kbarnum-1
	end if
				
	
	!perform bisection over price
	p1val = 2.0
	p2val = 2.5
	
	do piter=1,maxpit
		
		pval = (p1val+p2val)/2.0        
		wval = phi/pval

		
		!given pval, do the optimization    
		!$omp parallel private(zval,kval,Vnaval,kprimeval,statect,zprimect,&
		!$omp& aprimect,stateprimect,Vnextval1,Vnextval2,Vnextval,Vaval,xival) &
		!$omp& copyin(aval,RHSact,pval,wval,kbarfcstind,weight)
		!$omp do collapse(2)
		do zct=1,znum
		do kct=1,kdensenum
			if (distkz(zct,kct,t)>0.0) then
								
				zval = z0(zct)
				kval = kdense0(kct)
		
				!!!block to form Vna at the particular point (z,a,k,K)
				Vnaval = pval * freduced(zval,kval)
				kprimeval = max((1.0-delta)*kval,kmin)
				statect=(act-1)*znum+zct
				do zprimect=1,znum
				do aprimect=1,anum
					stateprimect = (aprimect-1)*znum+zprimect
			
					!now, need to evaluate Vold at (zprimect,aprimect,kprimeval,kbarfcstval)
					call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
					knum,kprimeval,Vnextval1)
					call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
					knum,kprimeval,Vnextval2)
			
					Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)
			
					!now that next period's value function is evaluated, add to expectation
					Vnaval = Vnaval + beta * pr_mat(statect,stateprimect) * Vnextval
				end do
				end do
				!!!end of block computing Vnaval
		
				!!!block computing Vaval
				RHSzct=zct        
				RHSkval=kval
				Vaval = brent(k0(1),k0(knum/2),k0(knum),fVa,brenttol,kprimeval)
				Vaval = -1.0 * Vaval; !this is to account for Brent routine's minimization
		
				xival = (Vaval - Vnaval)/phi
				!!end of block computing Vaval

				!store this info for distributional pushforward
				kprimep(zct,kct)=kprimeval
				xistarp(zct,kct)=xival
				padjustp(zct,kct) = cdfxi(xival)
		
			end if !dist weight > 0
		end do !kct
		end do !zct
		!$omp end do nowait
		!$omp end parallel
		
		
		!now, process this info into aggregates
		Yactualp=0.0; Iactualp = 0.0; Kprimeactualp=0.0; Nactualp = 0.0
		!$omp parallel private(zct,kct,zval,kval,xival,padjust) &
		!$omp& reduction(+:Yactualp,Iactualp,Kprimeactualp,Nactualp)
		!$omp do collapse(2)
		do zct=1,znum
		do kct=1,kdensenum
			if (distkz(zct,kct,t)>0.0) then
				zval = z0(zct)
				kval = kdense0(kct)
				kprimeval = kprimep(zct,kct)
				xival = xistarp(zct,kct)
				padjust = padjustp(zct,kct)
					
				Yactualp = Yactualp + distkz(zct,kct,t)*yreduced(zval,kval)
				Iactualp = Iactualp + distkz(zct,kct,t)*( padjust * (kprimeval - (1.0 - delta)*kval) )
				Kprimeactualp = Kprimeactualp + distkz(zct,kct,t) *  padjust*kprimeval &
				+ distkz(zct,kct,t) * (1.0 - padjust)*(1.0-delta)*kval
				Nactualp = Nactualp + distkz(zct,kct,t) * (nreduced(zval,kval) + expecxi(xival))  
			end if !dist tolerance
		end do !kct
		end do !zct
		!$omp end do nowait
		!$omp end parallel
		
		!now, compute the implied consumption, etc, via the array kprimep and xistarp from above
   		Cactualp = Yactualp - Iactualp
		perror = 1.0/pval - Cactualp

		!end of bisection loop over pval --> wval    
		if (abs(perror)<perrortol) exit
		if (perror<0.0) then
			p2val = pval
		else if (perror>0.0) then
			p1val = pval
		end if
		
    end do !piter
	

	
	if (mod(t,50)==1) then 
        write(13,"(A,I5,A,I5,A,F7.5,A,F7.5)") "t= ",t,", pit = ",piter,", err = ",perror,", p = ",pval
        write(*,"(A,I5,A,I5,A,F7.5,A,F7.5)") "t= ",t,", pit = ",piter,", err = ",perror,", p = ",pval
    end if 
    
    
    !now the value of p has been determined store aggregates
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
        kbaraggrule(t+1)=kfcstmat(act,1)+kfcstmat(act,2)*log(kbaraggrule(t))
        kbaraggrule(t+1) = exp(kbaraggrule(t+1))
        
        paggrule(t)=pfcstmat(act,1)+pfcstmat(act,2)*log(kbaraggrule(t))
        paggrule(t) = exp(paggrule(t))
        
    end if
    
    !now that market clearing has occured, insert distribution into next period's
    do zct=1,znum
    do kct=1,kdensenum
        if ( distkz(zct,kct,t)>0.0) then
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
        end if !dist>0
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
write(13,"(A,F15.1,A)") "Finished uncond simulation in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished uncond simulation in ",finish-start," seconds."

open(8,file="kdense0.txt")
do kct=1,kdensenum; write(8,*) kdense0(kct); end do
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

open(8,file="Nsim.txt")
do t=1,numper; write(8,*) Nsim(t); end do !t
close(8)

open(8,file="isim.txt")
do t=1,numper; write(8,*) isim(t); end do !t
close(8)

open(8,file="paggrule.txt")
do t=1,numper; write(8,*) paggrule(t); end do !t
close(8)

open(8,file="kbaraggrule.txt")
do t=1,numper; write(8,*) kbaraggrule(t); end do !t
close(8)


open(8,file="MICROsim.txt")
do t=1,numper
do momct=1,nummicro
write(8,*) MICROsim(momct,t)
end do !momct
end do !t
close(8)

if (doIRF==1) then

!initialize the distribution  
distIRF(:,:,:,:,:) = 0.0
do shockct=1,2
do simct=1,numsimIRF
    kbarsimIRF(1,simct,shockct) = 1.55
    distIRF(:,:,1,simct,shockct) = 1.0
    distIRF(:,:,1,simct,shockct) = distIRF(:,:,1,simct,shockct) / sum(distIRF(:,:,1,simct,shockct))
end do !simct
end do !shockct

do shockct=1,2
do simct=1,numsimIRF
do t=1,numperIRF-1

 	act=asimposIRF(t,simct,shockct)
    aval = a0(act)
	RHSact=act
	
    kbarval = kbarsimIRF(t,simct,shockct)
    kbarfcstval = kfcstmat(act,1) + kfcstmat(act,2)*log(kbarval); kbarfcstval = exp(kbarfcstval);
    
    pfcstval = pfcstmat(act,1) + pfcstmat(act,2)*log(kbarval); pfcstval = exp(pfcstval);
    
	kbarfcstind = kbarnum/2
	call hunt(kbar0,kbarnum,kbarfcstval,kbarfcstind); !kbarfcstval is in interval kbarfcstind

	if (kbarfcstind<=0) then 
		weight = 0.0
		kbarfcstind=1
	else if (kbarfcstind>=1.and.kbarfcstind<=(kbarnum-1)) then
		weight = (kbarfcstval - kbar0(kbarfcstind))/( kbar0(kbarfcstind+1)-kbar0(kbarfcstind) )
	else if (kbarfcstind>=kbarnum) then
		weight = 1.0
		kbarfcstind = kbarnum-1
	end if
				
					
	!perform bisection over price
	p1val = 2.0
	p2val = 2.5
	
	do piter=1,maxpit
		
		pval = (p1val+p2val)/2.0        
		wval = phi/pval

		
		!given pval, do the optimization    
		!$omp parallel private(zval,kval,Vnaval,kprimeval,statect,zprimect,&
		!$omp& aprimect,stateprimect,Vnextval1,Vnextval2,Vnextval,Vaval,xival) &
		!$omp& copyin(aval,RHSact,pval,wval,kbarfcstind,weight)
		!$omp do collapse(2)
		do zct=1,znum
		do kct=1,kdensenum
			if (distIRF(zct,kct,t,simct,shockct)>0.0) then
								
				zval = z0(zct)
				kval = kdense0(kct)
		
				!!!block to form Vna at the particular point (z,a,k,K)
				Vnaval = pval * freduced(zval,kval)
				kprimeval = max((1.0-delta)*kval,kmin)
				statect=(act-1)*znum+zct
				do zprimect=1,znum
				do aprimect=1,anum
					stateprimect = (aprimect-1)*znum+zprimect
			
					!now, need to evaluate Vold at (zprimect,aprimect,kprimeval,kbarfcstval)
					call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
					knum,kprimeval,Vnextval1)
					call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
					knum,kprimeval,Vnextval2)
			
					Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)
			
					!now that next period's value function is evaluated, add to expectation
					Vnaval = Vnaval + beta * pr_mat(statect,stateprimect) * Vnextval
				end do
				end do
				!!!end of block computing Vnaval
		
				!!!block computing Vaval
				RHSzct=zct        
				RHSkval=kval
				Vaval = brent(k0(1),k0(knum/2),k0(knum),fVa,brenttol,kprimeval)
				Vaval = -1.0 * Vaval; !this is to account for Brent routine's minimization
		
				xival = (Vaval - Vnaval)/phi
				!!end of block computing Vaval

				!store this info for distributional pushforward
				kprimep(zct,kct)=kprimeval
				xistarp(zct,kct)=xival
				padjustp(zct,kct) = cdfxi(xival)
		
			end if !dist weight > 0
		end do !kct
		end do !zct
		!$omp end do nowait
		!$omp end parallel
		
		
		!now, process this info into aggregates
		Yactualp=0.0; Iactualp = 0.0; Kprimeactualp=0.0; Nactualp = 0.0
		!$omp parallel private(zct,kct,zval,kval,xival,padjust) &
		!$omp& reduction(+:Yactualp,Iactualp,Kprimeactualp,Nactualp)
		!$omp do collapse(2)
		do zct=1,znum
		do kct=1,kdensenum
			if (distIRF(zct,kct,t,simct,shockct)>0.0) then
				zval = z0(zct)
				kval = kdense0(kct)
				kprimeval = kprimep(zct,kct)
				xival = xistarp(zct,kct)
				padjust = padjustp(zct,kct)
					
				Yactualp = Yactualp + distIRF(zct,kct,t,simct,shockct)*yreduced(zval,kval)
				Iactualp = Iactualp + distIRF(zct,kct,t,simct,shockct)*( padjust * (kprimeval - (1.0 - delta)*kval) )
				Kprimeactualp = Kprimeactualp + distIRF(zct,kct,t,simct,shockct) *  padjust*kprimeval &
				+ distIRF(zct,kct,t,simct,shockct) * (1.0 - padjust)*(1.0-delta)*kval
				Nactualp = Nactualp + distIRF(zct,kct,t,simct,shockct) * (nreduced(zval,kval) + expecxi(xival))  
			end if !dist tolerance
		end do !kct
		end do !zct
		!$omp end do nowait
		!$omp end parallel
		
		!now, compute the implied consumption, etc, via the array kprimep and xistarp from above
   		Cactualp = Yactualp - Iactualp
		perror = 1.0/pval - Cactualp

		!end of bisection loop over pval --> wval    
		if (abs(perror)<perrortol) exit
		if (perror<0.0) then
			p2val = pval
		else if (perror>0.0) then
			p1val = pval
		end if
		
    end do !piter
	
	if (mod(t,numperIRF/2)==1) then 
        write(13,"(A,I5,A,I5,A,I5,A,I5,A,F7.5,A,F7.5)") "t= ",t,", sim = ",simct,", ct = ",shockct,&
            ", pit = ",piter,", err = ",perror,", p = ",pval
        write(*,"(A,I5,A,I5,A,I5,A,I5,A,F7.5,A,F7.5)") "t= ",t,", sim = ",simct,", ct = ",shockct,&
            ", pit = ",piter,", err = ",perror,", p = ",pval
	end if 
    
    !now, the value of p has been determined, in a market clearing fashion
    kbarsimIRF(t+1,simct,shockct)=Kprimeactualp
    psimIRF(t,simct,shockct) = pval
    perrorsimIRF(t,simct,shockct) = perror
    kbarfcstsimIRF(t+1,simct,shockct) = kbarfcstval
    pfcstsimIRF(t,simct,shockct) = pfcstval
    ysimIRF(t,simct,shockct) = Yactualp
    isimIRF(t,simct,shockct) = Iactualp
    NsimIRF(t,simct,shockct) = Nactualp
        
    !now that market clearing has occured, insert distribution into next period's
    do zct=1,znum
    do kct=1,kdensenum
        if ( distIRF(zct,kct,t,simct,shockct)>0.0) then
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
            distIRF(zprimect,kprimeind,t+1,simct,shockct) = distIRF(zprimect,kprimeind,t+1,simct,shockct) &
                + pr_mat_z(zct,zprimect) * padjustp(zct,kct) * (1.0 - weight) * distIRF(zct,kct,t,simct,shockct)
            
            !transfer the weight to kprimeind+1 that does adjust
            distIRF(zprimect,kprimeind+1,t+1,simct,shockct) = distIRF(zprimect,kprimeind+1,t+1,simct,shockct) &
                + pr_mat_z(zct,zprimect) * padjustp(zct,kct) * weight * distIRF(zct,kct,t,simct,shockct)
    
            !transfer the weight to kprimeindnoadj that doesn't adjust
            distIRF(zprimect,kprimeindnoadj,t+1,simct,shockct) = distIRF(zprimect,kprimeindnoadj,t+1,simct,shockct) &
                + pr_mat_z(zct,zprimect) * (1.0 - padjustp(zct,kct)) * (1.0 - weightnoadj) * distIRF(zct,kct,t,simct,shockct)
            
            
            !transfer the weight to kprimeindnoadj+1 that doesn't adjust
            distIRF(zprimect,kprimeindnoadj+1,t+1,simct,shockct) = distIRF(zprimect,kprimeindnoadj+1,t+1,simct,shockct) &
                + pr_mat_z(zct,zprimect) * (1.0 - padjustp(zct,kct)) * weightnoadj * distIRF(zct,kct,t,simct,shockct)
            
        end do !zprimect
        end if !dist>0
    end do !kct
    end do !zct
    
    
end do !t
end do !simct
end do !shockct

finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished IRF in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished IRF in ",finish-start," seconds."  
  
!psimIRF 
open(8,file="psimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) psimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)

!kbarsimIRF 
open(8,file="kbarsimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) kbarsimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)    
    
!perrorsimIRF 
open(8,file="perrorsimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) perrorsimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)      
    
!pfcstsimIRF 
open(8,file="pfcstsimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) pfcstsimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)          
    
!kbarfcstsimIRF 
open(8,file="kbarfcstsimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) kbarfcstsimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)       
    
!ysimIRF 
open(8,file="ysimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) ysimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)     
    
!isimIRF 
open(8,file="isimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) isimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)    
    
!NsimIRF 
open(8,file="NsimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) NsimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)    
    
          
end if !doIRF




finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished in ",finish-start," seconds."
close(13); !closing log file


contains


subroutine discretize_simulate()
implicit none

integer :: zct,ergit,act,statect,stateprimect,ct,t,aprimect
double precision :: ergerror

!capital grids
call linspace(k0,log(kmin),log(kmax),knum); k0 = exp(k0);
call linspace(kbar0,kbarmin,kbarmax,kbarnum);
call linspace(kdense0,log(kmin),log(kmax),kdensenum); kdense0=exp(kdense0);

!discretize idio prod process
call tauchen(znum,rhoz,sigmaz,nstdevz,pr_mat_z,z0)

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

!determine the ergodic distribution of z0
write(*,*) "Finding ergodic distribution of z0.";
pr_mat_z_dum = pr_mat_z
ergz0old(:)=1.0/dble(znum)
do ergit=1,maxergit
    pr_mat_z_dum = matmul(pr_mat_z,pr_mat_z_dum)
    ergz0(:) = pr_mat_z_dum(1,:)
    ergerror = maxval(abs(ergz0 - ergz0old))
    !if (mod(ergit,100)==0) write(*,*) "ergit = ",ergit,"ergerror = ",ergerror
    if (ergerror<ergerrortol) exit
    ergz0old = ergz0
end do
ergz0 = ergz0/sum(ergz0)
if (ergit<maxergit) write(*,*) "Found ergodic distribution for z. Error = ",ergerror
if (ergit>=maxergit) write(*,*) "Problem with ergodic distribution for z."


!discretize agg prod process
call tauchen(anum,rhoa,sigmaa,nstdeva,pr_mat_a,a0)

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

!create unified transition matrix
do zct=1,znum
do act=1,anum
    statect = (act-1)*znum + zct
    do zprimect=1,znum
    do aprimect=1,anum
        stateprimect = (aprimect-1)*znum+zprimect
        pr_mat(statect,stateprimect)=pr_mat_a(act,aprimect)*pr_mat_z(zct,zprimect)
    end do
    end do
end do
end do

do statect=1,znum*anum
    pr_mat(statect,:) = pr_mat(statect,:) / sum(pr_mat(statect,:))
end do !statect




!draw random seed array size
call random_seed(size=seeddim)
allocate(seedarray(seeddim))


!call the simulation random numbers once, so they don't change with iteration
do ct=1,seeddim
    seedarray(ct) = seedint + ct
end do !ct
call random_seed(put=seedarray)


call random_number(asimshock)
asimpos(1) = ainit

do t=2,numper
    asimgrid(1)=pr_mat_a(asimpos(t-1),1)
    do act=2,anum
        asimgrid(act) = asimgrid(act-1) + pr_mat_a(asimpos(t-1),act)
    end do !act
    
      
    if (asimshock(t)<asimgrid(1)) then 
        aprimect=1; 
    end if
    do act=2,anum
        if (asimgrid(act-1)<=asimshock(t).and.asimshock(t)<asimgrid(act)) then
            aprimect=act
        end if
    end do !act
    asimpos(t) = aprimect
end do


open(8,file="asimpos.txt")
do t=1,numper
write(8,*) asimpos(t)
end do !t
close(8)

open(8,file="asimshock.txt")
do t=1,numper
write(8,*) asimshock(t)
end do !t
close(8)





!if doing IRF, determine ergodic distribution of agg prod, then IRF shock prob, then perform IRF agg prod simulation
if (doIRF==1) then

call random_number(asimshockIRF)
call random_number(arriveshockIRF)

!initialize ergodic distribution
ergdistaold(:) = 0.0
ergdistaold(anum/2) = 1.0

do ct=1,maxergit
    do act=1,anum
    do aprimect=1,anum
        ergdista(aprimect) = ergdista(aprimect) + pr_mat_a(act,aprimect)*ergdistaold(act)
    end do !aprimect
    end do !act
    
    if (maxval(abs(ergdista-ergdistaold))<ergdistatol) exit
    
    ergdistaold = ergdista
    ergdista(:) = 0.0
    
end do !ct

!now that the ergodic distribution is done, compute size of IRF shock probability
shockprobdenom = 0.0
shockprobdenom = log(a0(shockanumIRF))
do act=1,anum
    shockprobdenom = shockprobdenom  - ergdista(act) * log(a0(act))
end do !act
shockprobIRF = shocksizeIRF / shockprobdenom
write(*,*) "shockprobIRF = ",shockprobIRF
!now, perform the simulation
!asimshockIRF(numperIRF,numsimIRF),asimposIRF(numperIRF,numsimIRF,2)

asimposIRF(1,:,:) = ainit

do ct=1,numsimIRF !counting IRF simulations
do t=2,numperIRF !counting IRF
    
    if (t<shockperIRF) then
        
        !both transit as normal
        act = asimposIRF(t-1,ct,1)
        
        !bins to use for comparison
        asimgrid(:) = 0.0
        asimgrid(1) = pr_mat_a(act,1)
        do aprimect=2,anum
            asimgrid(aprimect) = asimgrid(aprimect-1) + pr_mat_a(act,aprimect)
        end do !aprimect
        
        if (asimshockIRF(t,ct)<asimgrid(1)) then
            asimposIRF(t,ct,1) = 1
            asimposIRF(t,ct,2) = 1
        else if (asimshockIRF(t,ct)>=asimgrid(1)) then
            do aprimect=2,anum
                if (asimshockIRF(t,ct)>=asimgrid(aprimect-1).and.asimshockIRF(t,ct)<asimgrid(aprimect)) then
                    asimposIRF(t,ct,1) = aprimect
                    asimposIRF(t,ct,2) = aprimect
                end if
            end do !aprimect
        end if
        
    else if (t==shockperIRF) then
        
        !first, decide if you're transiting normally or not
        if (arriveshockIRF(ct)>shockprobIRF) then
                
                !in this case, transit normally as above
                act = asimposIRF(t-1,ct,1)
                
                !bins to use for comparison
                asimgrid(:) = 0.0
                asimgrid(1) = pr_mat_a(act,1)
                do aprimect=2,anum
                    asimgrid(aprimect) = asimgrid(aprimect-1) + pr_mat_a(act,aprimect)
                end do !aprimect
                
                if (asimshockIRF(t,ct)<asimgrid(1)) then
                    asimposIRF(t,ct,1) = 1
                    asimposIRF(t,ct,2) = 1
                else if (asimshockIRF(t,ct)>=asimgrid(1)) then
                    do aprimect=2,anum
                        if (asimshockIRF(t,ct)>=asimgrid(aprimect-1).and.asimshockIRF(t,ct)<asimgrid(aprimect)) then
                            asimposIRF(t,ct,1) = aprimect
                            asimposIRF(t,ct,2) = aprimect
                        end if
                    end do !aprimect
                end if
            
            
            
            
        else if (arriveshockIRF(ct)<=shockprobIRF) then
            !in this case, transit normally only for version 2, not version 1
            
            !version 2 is shocked, goes to shockasnumIRF position
            asimposIRF(t,ct,2) = shockanumIRF
            
            !version 1 transits normally
            act = asimposIRF(t-1,ct,1)
            
            !bins to use for comparison
            asimgrid(:) = 0.0
            asimgrid(1) = pr_mat_a(act,1)
            do aprimect=2,anum
                asimgrid(aprimect) = asimgrid(aprimect-1) + pr_mat_a(act,aprimect)
            end do !aprimect
            
            if (asimshockIRF(t,ct)<asimgrid(1)) then
                asimposIRF(t,ct,1) = 1
            else if (asimshockIRF(t,ct)>=asimgrid(1)) then
                do aprimect=2,anum
                    if (asimshockIRF(t,ct)>=asimgrid(aprimect-1).and.asimshockIRF(t,ct)<asimgrid(aprimect)) then
                        asimposIRF(t,ct,1) = aprimect
                    end if
                end do !aprimect
            end if
        
            
            
        end if 
    
    else if (t>shockperIRF) then
        
        !both transit as normal, but separately
    

        !version 1 transits normally
        act = asimposIRF(t-1,ct,1)
        
        !bins to use for comparison
        asimgrid(:) = 0.0
        asimgrid(1) = pr_mat_a(act,1)
        do aprimect=2,anum
            asimgrid(aprimect) = asimgrid(aprimect-1) + pr_mat_a(act,aprimect)
        end do !aprimect
        
        if (asimshockIRF(t,ct)<asimgrid(1)) then
            asimposIRF(t,ct,1) = 1
        else if (asimshockIRF(t,ct)>=asimgrid(1)) then
            do aprimect=2,anum
                if (asimshockIRF(t,ct)>=asimgrid(aprimect-1).and.asimshockIRF(t,ct)<asimgrid(aprimect)) then
                    asimposIRF(t,ct,1) = aprimect
                end if
            end do !aprimect
        end if
        
            
        !version 2 transits normally
        act = asimposIRF(t-1,ct,2)
        
        !bins to use for comparison
        asimgrid(:) = 0.0
        asimgrid(1) = pr_mat_a(act,1)
        do aprimect=2,anum
            asimgrid(aprimect) = asimgrid(aprimect-1) + pr_mat_a(act,aprimect)
        end do !aprimect
        
        if (asimshockIRF(t,ct)<asimgrid(1)) then
            asimposIRF(t,ct,2) = 1
        else if (asimshockIRF(t,ct)>=asimgrid(1)) then
            do aprimect=2,anum
                if (asimshockIRF(t,ct)>=asimgrid(aprimect-1).and.asimshockIRF(t,ct)<asimgrid(aprimect)) then
                    asimposIRF(t,ct,2) = aprimect
                end if
            end do !aprimect
        end if
        
        
    
    end if
    
end do !t
end do !ct


open(8,file="asimposIRF.txt")
do ct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) asimposIRF(t,ct,shockct)
end do !shockct
end do !t
end do !ct
close(8)



end if !doIRF==1

end subroutine discretize_simulate


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

double precision function freduced(z,k)
implicit none
double precision :: z,k

double precision :: exponentnu

exponentnu = 1.0 / (1.0 - nu)

freduced = ( nu ** (nu * exponentnu) ) * (1.0 - nu) * ( wval ** ( -1.0 * nu * exponentnu ) )
freduced = freduced * (z ** exponentnu) * (aval ** exponentnu) * ( k ** (alpha * exponentnu) )
end function freduced


double precision function nreduced(z,k)
implicit none
double precision :: z,k

double precision :: exponentnu

exponentnu = 1.0 / (1.0 - nu)

nreduced = ( nu **  exponentnu ) * ( wval ** ( -1.0 * exponentnu ) )
nreduced = nreduced * (z ** exponentnu) * (aval ** exponentnu) * ( k ** (alpha * exponentnu) )
end function nreduced

double precision function yreduced(z,k)
implicit none
double precision :: z,k

double precision :: exponentnu

exponentnu = 1.0 / (1.0 - nu)

yreduced = ( (nu/wval) ** (nu * exponentnu) ) * ( (aval * z ) ** exponentnu  )
yreduced = yreduced * ( k ** ( alpha * exponentnu ) )

end function yreduced

double precision function fVa(kprime)
implicit none

double precision :: kprime

integer :: ct,act,zct,aprimect,zprimect,kbarct,zdumct,adumct,&
    statect,stateprimect

double precision :: zval,kval,Vnextval,Vnextval1,Vnextval2

zct = RHSzct
act = RHSact
kval = RHSkval

zval = z0(zct)

fVa = freduced(zval,kval) - kprime + (1.0 - delta) * kval
fVa = fVa * pval

statect=(act-1)*znum + zct
do zprimect=1,znum
do aprimect=1,anum
    stateprimect=(aprimect-1)*znum+zprimect
    
    call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
    knum,kprime,Vnextval1)
    call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
    knum,kprime,Vnextval2)
    
    Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)        

    !now that Vnextval is evaluated, add to expectation
    fVa = fVa + beta * pr_mat(statect,stateprimect) * Vnextval
    
end do 
end do

fVa  = -1.0 * fVa; !take into account the minimization in the Brent routine
end function fVa

end program kt_xpa