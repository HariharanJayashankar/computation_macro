!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kt_param.f90
!
! Fortran code for the PARAM solution of the Khan and Thomas (2008) 
! model.
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
!to subprograms when set elsewhere.

!real parameters
double precision, parameter :: alpha = 0.256 !capital elasticity
double precision, parameter :: nu = 0.640 !labor elasticity
double precision, parameter :: phi = 2.4 !labor disutility
double precision, parameter :: xibar = 0.0083 !upper bound of capital AC distribution
double precision, parameter :: beta = 0.977 !discount rate
double precision, parameter :: delta = 0.069 !capital depreciation rate
double precision, parameter :: kmin= 0.1 !min value on idio capital grid
double precision, parameter :: kmax = 8.0 !max value on idio capital grid
double precision, parameter :: kbarmin= 1.25 !min value on agg capital grid
double precision, parameter :: kbarmax = 2.0 !max value on agg capital grid
double precision, parameter :: rhoz = 0.859 !persistence of idio prod
double precision, parameter :: sigmaz = 0.022 !stdev of shock to log idio prod
double precision, parameter :: nstdevz = 2.0 !multiple of st dev of idio prod process to discretize
double precision, parameter :: rhoa = 0.859 !persistence of idio prod
double precision, parameter :: sigmaa = 0.014 !stdev of shock to log idio prod
double precision, parameter :: nstdeva = 2.0 !multiple of st dev of idio prod process to discretize
double precision, parameter :: vferrortol=1e-4 !tolerance for vf error 
double precision, parameter :: kprimeerrortol = 1e-4 !tolerance for error on kprime
double precision, parameter :: xistarerrortol = 1e-4 !tolerance for error on xistar
double precision, parameter :: plb = 2.0 !lower bound for bisection on price
double precision, parameter :: pub = 2.4 !upper bound for bisection on price
double precision, parameter :: perrortol = 1e-5 !tolerance for bisection on price
double precision, parameter :: poltol = 1e-6 !tolerance for optimization of idio capital policies
double precision, parameter :: disterrortol = 1e-5!tolerance on distribution fixed point
double precision, parameter :: broydenrhotol = 1e-6!tolerance on change in rho values in Broyden optimization
double precision, parameter :: broydengradtol = 0.01 !tolerance on size of gradient in Broyden optimization
double precision, parameter :: broydenfunctol = 1e-7 !tolerance on size of function eval diff
double precision, parameter :: stepsize = 0.001 !Broyden step size
double precision, parameter :: reftol = 1e-2 !tolerance for convergence of reference moments
double precision, parameter :: refdamp = 0.5 !dampening parameter for update on reference moments
double precision, parameter :: pfixeddamp = 0.025 !dampening parameter on p in fixed-point iteration
double precision, parameter :: Kfixeddamp = 0.05 !dampening parameter on K' in fixed-point iteration
double precision, parameter :: fixedtol = 1e-3 !relative tolerance fixed-point iteration
double precision, parameter :: shocksizeIRF=sigmaa !percentage shock size
double precision, parameter :: ergdistatol = 1e-5 !tolerance for computing aggregate ergodic distribution


!integer parameters
integer, parameter :: doIRF = 1 !do IRF calculations?
integer, parameter :: knum = 10 !number of idio capital nodes for spline solution
integer, parameter :: znum = 5 !number of idio prod states in discretization
integer, parameter :: anum = 5 !number of agg prod states in discretization
integer, parameter :: kbarnum = 10 !number of agg capital states in discretization
integer, parameter :: kdensenum = 50 !number of idio capital points in distribution
integer, parameter :: maxvfit = 100 !max number of policy iterations
integer, parameter :: maxpit = 100 !max number of iterations on price bisection
integer, parameter :: maxdistit = 5000 !max number of iterations for distribution fixed point
integer, parameter :: numconstants = 14 !number of constants to pass to output
integer, parameter :: momnum = 25 !number of moments to make available for parametrized distributions
integer, parameter :: momuse = 4 !number of moments to use for parametrized distributions
integer, parameter :: nsimp = 36 !number of Simpson quadrature nodes (needs to be even)
integer, parameter :: maxbroydit = 150000 !max number of Broyden iters
integer, parameter :: maxrefit = 1 !max number of iterations on reference moments
integer, parameter :: seedint = 2503 !RNG seed
integer, parameter :: maxfixedit = 500 !max number of iterations for fixed point on (K',p)
integer, parameter :: numper = 2500 !# of periods to simulate
integer, parameter :: numdiscard = 500 !# of periods to discard
integer, parameter :: ainit = 3 !initial agg prod state
integer, parameter :: kbarinit = kbarnum/2 !initial agg capital state
integer, parameter :: OPTMETHOD = 3 !1 is broyden, 2 is PSO, 3 is SR1, referring to solution of dist coeffs
integer, parameter :: nummicro = 7 !micro moments to compute
integer, parameter :: numsimIRF=2000
integer, parameter :: numperIRF=50
integer, parameter :: shockperIRF=25
integer, parameter :: shockanumIRF = anum !grid point that you impose in IRF period

!parameters for SR1 quasi-newton with step size line search
integer, parameter :: maxstepnum = 100
double precision, parameter :: stepfrac = 0.9
double precision, parameter :: stepinit = 0.09
double precision, parameter :: c1 = 1e-3
double precision, parameter :: c2 = 0.9
double precision, parameter :: stepfallback = 0.00001

!insert the pso parameters for rho parameter minimization
double precision, parameter :: bound = 30.0
integer, parameter :: npart = 500
integer, parameter :: nvar = znum*momuse
double precision, parameter :: xtol = 1.0e-3
double precision, parameter :: ftol = 1.0e-3
double precision, parameter :: xquicktol = 1.0e-3
integer, parameter :: xquicknum = 25
double precision :: phipso(2)
integer, parameter :: maxpsoit =5000

!moment stuff to be available to subfunctions
double precision :: rhomat(znum,momuse),intvec(znum),mommat(znum,momuse),momstoremat(znum,momnum)
double precision :: simpnodes(nsimp+1),simpweights(nsimp+1)

!other stuff to be available
double precision, allocatable :: Vold(:,:,:,:),V2old(:,:,:,:),pr_mat(:,:),k0(:),pr_mat_z(:,:),pr_mat_a(:,:),&
    Va(:,:,:,:),Vna(:,:,:,:),z0(:),a0(:),ergdistz(:),momREF(:,:,:),rhoREF(:,:,:),intREF(:,:),V(:,:,:,:),kprime(:,:,:,:),&
    xistar(:,:,:,:),pVF(:,:),KprimeVF(:,:),kbar0(:),rhoSS(:,:),intSS(:),kdense0(:),ergoldz(:),ergnewz(:),&
    asimshockIRF(:,:),ergdista(:),ergdistaold(:),arriveshockIRF(:),asimshock(:),distkzsim(:,:,:)
    
double precision :: pval,wval,aval,kbarval,kbarfcstwgt

integer :: seeddim

integer, allocatable :: seedarray(:),asimpos(:),asimposIRF(:,:,:)

integer :: kbarfcstind,RHSzct,RHSact,RHSkbarct
!$omp threadprivate(RHSzct,RHSact,RHSkbarct)

end module modparams

program kt_param
use base_lib
use omp_lib
use modparams
implicit none

integer :: ct,momct,zct,act,refct,vfct,kbarct,kct,zprimect,statect,stateprimect,aprimect,&
    ind,fixediter,t,shockct,simct
    
   
double precision :: start,finish,pSS,KSS,kbarprimeval,ergerrz,wgt,kprimeval,Kprimeactual,Cactual,&
    perror,Kprimeerror,fixederror,vferror,xistarerror,kpolerror,Yactual,Iactual,Nactual,kval,Pintval,&
    pfcstval,kbarfcstval,momREFerror,shockprobdenom,shockprobIRF,zval,ival,xival,padjust

double precision, allocatable :: constantvec(:),momSS(:,:),kprimeSS(:,:),&
    kprimeold(:,:,:,:),xistarold(:,:,:,:),VSS(:,:),kprimevec(:),RHSvec(:),&
    xistarsimp(:,:),kprimesimp(:,:),xistarSS(:,:),momSIM(:,:,:),rhoSIM(:,:,:),intSIM(:,:),&
    kbarsim(:),psim(:),ysim(:),isim(:),nsim(:),densitySIM(:,:,:),lb(:),ub(:),rhovec(:),&
    gradvec(:),momstoreSS(:,:),pfcstsim(:),kbarfcstsim(:),momREFstore(:,:,:,:),momREFnew(:,:,:),&
    pDHsim(:),kbarDHsim(:),ysimIRF(:,:,:),isimIRF(:,:,:),nsimIRF(:,:,:),psimIRF(:,:,:),&
    kbarsimIRF(:,:,:),rhoIRF(:,:,:,:,:),momIRF(:,:,:,:,:),intIRF(:,:,:,:),MICROsim(:,:)

!!!!!!PROGRAM PRELIMINARIES

!open log file
start = omp_get_wtime()
open(13,file="kt_param.txt")

!parallelization check

!$omp parallel
write(*,*) "Parallel hello to you!"
write(13,*) "Parallel hello to you!"
!$omp end parallel

write(*,*) " "
write(13,*) " "

!do allocate memory based on dimensions set in parameter module
allocate(V(znum,knum,anum,kbarnum),Vold(znum,knum,anum,kbarnum),V2old(znum,knum,anum,kbarnum),pr_mat_a(anum,anum),&
    pr_mat_z(znum,znum),pr_mat(anum*znum,anum*znum),k0(knum),kdense0(kdensenum),z0(znum),a0(anum),kbar0(kbarnum),&
    constantvec(numconstants),momSS(znum,momuse),rhoSS(znum,momuse),intSS(znum),momREF(znum,momuse,anum),&
    rhoREF(znum,momuse,anum),intREF(znum,anum),kprimeSS(znum,knum),kprimeold(znum,knum,anum,kbarnum),&
    kprime(znum,knum,anum,kbarnum),Va(znum,knum,anum,kbarnum),Vna(znum,knum,anum,kbarnum),ergdistz(znum),&
    ergnewz(znum),ergoldz(znum),VSS(znum,knum),pVF(anum,kbarnum),KprimeVF(anum,kbarnum),&
    kprimevec(znum),RHSvec(znum),xistarsimp(znum,nsimp+1),kprimesimp(znum,nsimp+1),&
    xistar(znum,knum,anum,kbarnum),xistarold(znum,knum,anum,kbarnum),xistarSS(znum,knum),momSIM(znum,momuse,numper),&
    rhoSIM(znum,momuse,numper),intSIM(znum,numper),asimpos(numper),asimshock(numper),kbarsim(numper),&
    psim(numper),ysim(numper),isim(numper),nsim(numper),densitySIM(znum,kdensenum,numper),lb(nvar),ub(nvar),&
    rhovec(znum*momuse),gradvec(znum*momuse),momstoreSS(znum,momnum),pfcstsim(numper),kbarfcstsim(numper),&
    momREFstore(znum,momuse,anum,maxrefit),momREFnew(znum,momuse,anum),pDHsim(numper),kbarDHsim(numper),&
    asimposIRF(numperIRF,numsimIRF,2),asimshockIRF(numperIRF,numsimIRF),ergdista(anum),ergdistaold(anum),&
    arriveshockIRF(numsimIRF),ysimIRF(numperIRF,numsimIRF,2),isimIRF(numperIRF,numsimIRF,2),nsimIRF(numperIRF,numsimIRF,2),&
    psimIRF(numperIRF,numsimIRF,2),kbarsimIRF(numperIRF,numsimIRF,2),rhoIRF(znum,momuse,numperIRF,numsimIRF,2),&
    momIRF(znum,momuse,numperIRF,numsimIRF,2),intIRF(znum,numperIRF,numsimIRF,2),MICROsim(nummicro,numper),&
	distkzsim(znum,kdensenum,numper))

!!!!!!!SET UP GRIDS AND DISCRETIZE IDIO PROD

!set up pso
phipso = (/2.05,2.05/)
lb(:) = -1.0*bound
ub(:) = 1.0 * bound

!constants to store
open(8,file="constantvec.txt")
write(8,*) znum,anum,knum,kbarnum,kdensenum,numper,numdiscard,nsimp,kmin,kmax,momuse,&
	numperIRF,numsimIRF,shockperIRF
close(8)

call discretize_simulate()

!call Simpson quadrature weights and nodes for moment integration
call qsimpweightsnodes(kmin,kmax,nsimp,simpweights,simpnodes)

!!!!!!!!!READ IN SS COEFFICIENTS AND MOMENTS
open(8,file="momSS.txt")
do zct=1,znum
do momct=1,momnum
    read(8,*) momstoreSS(zct,momct)
end do !momct
end do !zct
close(8)

momSS = momstoreSS(:,1:momuse)

!solve for the coefficients and integrals of the SS distribution
mommat = momSS
call findrhoSR1wolfe(0)
rhoSS = rhomat
intSS = intvec

write(13,*) "Found coefficients for SS distribution."
write(*,*) "Found coefficients for SS distribution."

write(13,*) " "
write(*,*) " "

!insert these SS moments into the array storing current value of reference moments, 
!coefficients, and integrals
do act=1,anum
    momREF(:,:,act) = momSS
    rhoREF(:,:,act) = rhoSS
    intREF(:,act) = intSS
end do !act

!read in SS capital policy, for initialization
open(8,file="kprimeSS.txt")
do zct=1,znum
do kct=1,knum
    read(8,*) kprimeSS(zct,kct)
end do !kct
end do !zct
close(8)

!initialize the policy function with SS values
do act=1,anum
do kbarct=1,kbarnum
    kprimeold(:,:,act,kbarct) = kprimeSS
end do !kbarct
end do !act

!read in SS cutoff values, for initialization
open(8,file="xistarSS.txt")
do zct=1,znum
do kct=1,knum
    read(8,*) xistarSS(zct,kct)
end do !kct
end do !zct
close(8)

!initialize the cutoff values with SS values
do act=1,anum
do kbarct=1,kbarnum
    xistarold(:,:,act,kbarct) = xistarSS
end do !kbarct
end do !act

!read in SS values of the VF
open(8,file="VSS.txt")
do zct=1,znum
do kct=1,knum
read(8,*) VSS(zct,kct)
end do !kct
end do !zct
close(8)

!initialize the VF with SS values, and spline the resulting VF
do act=1,anum
do kbarct=1,kbarnum
    Vold(:,:,act,kbarct) = VSS
end do !kbarct
end do !act

!spline the resulting VF
V2old(:,:,:,:) = 0.0
do act=1,anum
do kbarct=1,kbarnum
    do zct=1,znum
        call spline(k0,Vold(zct,:,act,kbarct),knum,dble(1.0e30),dble(1.0e30),V2old(zct,:,act,kbarct))
    end do !zct
end do !kbarnum
end do !act

!read in SS price and agg capital
open(8,file="pandKSS.txt")
read(8,*) pSS
read(8,*) KSS
close(8)

!set up agg capital grid
if (kbarnum>1) then
call linspace(kbar0,log(kbarmin),log(kbarmax),kbarnum);
kbar0=exp(kbar0);
else if (kbarnum==1) then
kbar0 = KSS
end if

!initialize guesses for the mapping (a,K) --> (p,K')
do act=1,anum
do kbarct=1,kbarnum
    pVF(act,kbarct) = pSS
    KprimeVF(act,kbarct) = KSS
end do !kbarct
end do !act

!!!!!BIG LOOP ON REFERENCE MOMENTS - THIS IS MODEL SOLUTION LOOP
do refct=1,maxrefit

!!!!!VFI, GIVEN CURRENT VALUES OF REFERENCE MOMENTS
write(*,*) "Performing VFI for reference loop ",refct,"."
write(13,*) "Performing VFI for reference loop ",refct,"."
write(*,*) " "
write(13,*) " "
do vfct=1,maxvfit
    write(*,*) " "
    write(13,*) " "
    
    write(*,*) "Now doing VF iter number ",vfct,"."
    write(13,*) "Now doing VF iter number ",vfct,"."
    
    !in each VF iteration, loop over aggregate states
    do act=1,anum
    do kbarct=1,kbarnum
        
        !what are aggregate states?
        aval = a0(act)
        kbarval = kbar0(kbarct)
        RHSact = act !this needs to be available within RHS
        RHSkbarct = kbarct !this can be available in RHS functions as well
        
        !given Vold as well as (a,K), solve for fixed point on (p,K')
        
        !start apparatus for solution on (p,K')
        pval = pVF(act,kbarct)
        kbarprimeval = KprimeVF(act,kbarct)
        
        do fixediter=1,maxfixedit

            wval = phi / pval !available inside RHS
            
            !what are linear interpolation values for next period VF in agg K'?
            call kbarprimelinterp(kbarprimeval,kbarfcstind,kbarfcstwgt,kbar0,kbarnum)
            
            !this loop obtains the policies k'(z), as well as (part of) the RHS at these values
            do zct=1,znum
                RHSzct = zct
                RHSvec(zct) = -1.0*brent(k0(1),k0(knum/2),k0(knum),fnRHS,poltol,kprimeval)
                kprimevec(zct) = kprimeval
            end do !zct
            
            !construct the policies by calling subroutine (only constructs them for RHSact,RHSkbarct)
            call constructPolicySimp(kprimevec,xistarsimp,kprimesimp)

            !given policies at integration points, can construct implied K',C, which requires moment, coeff,
            ! and integral arrays
            mommat = momREF(:,:,RHSact)
            rhomat = rhoREF(:,:,RHSact)
            intvec = intREF(:,RHSact)
            call impliedaggs(xistarsimp,kprimesimp,Cactual,Kprimeactual,Yactual,Iactual,Nactual)

            !compute FPI errors on agg capital and price
            perror = (1.0/pval - Cactual) 
            Kprimeerror = (Kprimeactual - kbarprimeval)/abs(kbarprimeval)
            fixederror = maxval(abs((/perror,Kprimeerror/)))
            
            !exit FPI if converged
            if (fixederror < fixedtol ) exit
            
            !if not converged, update
            pval = pval + pfixeddamp * perror
            kbarprimeval = kbarprimeval + Kfixeddamp * Kprimeerror
            
        end do !fixediter
        
        write(*,"(A,F5.3,A,F5.3,A,F7.5,A,I6,A,F7.5,A,F7.5)") "K' = ",kbarprimeval,", p = ",pval,&
            " err = ",fixederror,", it = ",fixediter-1," a = ",aval," K = ",kbarval
        write(13,"(A,F5.3,A,F5.3,A,F7.5,A,I6,A,F7.5,A,F7.5)") "K' = ",kbarprimeval,", p = ",pval,&
            " err = ",fixederror,", it = ",fixediter-1," a = ",aval," K = ",kbarval
    
        !store price and agg capital next period associated with (a,K)
        pVF(act,kbarct) = pval
        KprimeVF(act,kbarct) = kbarprimeval
        !close apparatus for fixed point on (p,K')
        
        !then, once converged, compute Bellman RHS for this value of (a,K)
        call constructBellman(kprimevec)
         
    end do !act
    end do !kbarct
    
    !compute Bellman equation error and test for convergence
    vferror = maxval(abs(log(V) - log(Vold)))
    xistarerror = maxval(abs(xistar - xistarold))
    kpolerror = maxval(abs(kprime - kprimeold))
    
    !test for convergence and exit if needed
    write(*,"(A,F7.5,A,F7.5,A,F7.5)") " VFerr = ",vferror," XIerr = ",xistarerror," k'err = ",kpolerror
    write(13,"(A,F7.5,A,F7.5,A,F7.5)") " VFerr = ",vferror," XIerr = ",xistarerror," k'err = ",kpolerror
    if (vferror<vferrortol) exit
    
    !if not yet converged, iterate stuff forward and respline the VF
    Vold = V
    kprimeold = kprime
    xistarold = xistar
    
    !spline the VF for the next iteration
    V2old(:,:,:,:) = 0.0
    do act=1,anum
    do kbarct=1,kbarnum
    do zct=1,znum
        call spline(k0,Vold(zct,:,act,kbarct),knum,dble(1.0e30),dble(1.0e30),V2old(zct,:,act,kbarct))
    end do !zct
    end do !kbarnum
    end do !act    
    
end do !vfct

write(*,*) "Done with VFI for ref loop ",refct," at ",omp_get_wtime()-start," seconds."
write(13,*) "Done with VFI for ref loop ",refct," at ",omp_get_wtime()-start," seconds."

write(*,*) " "
write(13,*) " "

!write output files

open(8,file="a0.txt")
do act=1,anum
write(8,*) a0(act)
end do 
close(8)

open(8,file="z0.txt")
do zct=1,znum
write(8,*) z0(zct)
end do 
close(8)

open(8,file="k0.txt")
do kct=1,knum
write(8,*) k0(kct)
end do 
close(8)

open(8,file="kdense0.txt")
do kct=1,kdensenum
write(8,*) kdense0(kct)
end do 
close(8)

open(8,file="kbar0.txt")
do kct=1,kbarnum
write(8,*) kbar0(kct)
end do 
close(8)

open(8,file="simpweights.txt")
do kct=1,nsimp+1
write(8,*) simpweights(kct)
end do !kct
close(8)

open(8,file="simpnodes.txt")
do kct=1,nsimp+1
write(8,*) simpnodes(kct)
end do !kct
close(8)

open(8,file="ergdistz.txt")
do zct=1,znum
write(8,*) ergdistz(zct)
end do !zct
close(8)

open(8,file="pr_mat_z.txt")
do zct=1,znum
write(8,*) pr_mat_z(zct,:)
end do !zct
close(8)

open(8,file="pr_mat_a.txt")
do act=1,anum
write(8,*) pr_mat_a(act,:)
end do !zct
close(8)

open(8,file="V.txt")
do zct=1,znum
do kct=1,knum
do act=1,anum
do kbarct=1,kbarnum
write(8,*) V(zct,kct,act,kbarct)
end do 
end do 
end do 
end do 
close(8)

open(8,file="kprime.txt")
do zct=1,znum
do kct=1,knum
do act=1,anum
do kbarct=1,kbarnum
write(8,*) kprime(zct,kct,act,kbarct)
end do 
end do 
end do 
end do 
close(8)

open(8,file="xistar.txt")
do zct=1,znum
do kct=1,knum
do act=1,anum
do kbarct=1,kbarnum
write(8,*) xistar(zct,kct,act,kbarct)
end do 
end do 
end do 
end do 
close(8)

open(8,file="pVF.txt")
do act=1,anum
do kbarct=1,kbarnum
write(8,*) pVF(act,kbarct)
end do 
end do 
close(8)

open(8,file="KprimeVF.txt")
do act=1,anum
do kbarct=1,kbarnum
write(8,*) KprimeVF(act,kbarct)
end do 
end do 
close(8)

!!!!!END VFI

!!!!!NOW, SIMULATE THE MODEL

!first, need to initialize the moments
kbarsim(:) = 0.0
kbarval = 1.55 !what is agg capital at assumed starting point?
kbarsim(1) = kbarval

!store moments, coefficients, and integrals
momSIM(:,:,1) = momREF(:,:,asimpos(1))
momSIM(:,1,1) = momREF(:,1,asimpos(1)) * (kbarval / sum(ergdistz * momREF(:,1,asimpos(1))))
intSIM(:,1) = intREF(:,asimpos(1))
rhoSIM(:,:,1) = rhoREF(:,:,asimpos(1))

pfcstsim(:) = 0.0
kbarfcstsim(:) = 0.0
pDHsim(:) = 0.0
kbarDHsim(:) = 0.0
psim(:) = 0.0
ysim(:) = 0.0
isim(:) = 0.0
nsim(:) = 0.0

MICROsim(:,:) = 0.0

open(8,file="distkzsim.txt")
close(8,status='DELETE')	


do t=1,numper-1
    
    !for this period, what is the aggregate productivity?
    act = asimpos(t)
    aval = a0(act)
    RHSact = act
    
    kbarval = kbarsim(t)
    
    !obtain forecasted values of (p,K')        
    call kbarprimelinterp(kbarval,ind,wgt,kbar0,kbarnum)        
    pfcstval = (1.0 - wgt) * pVF(act,ind) + wgt * pVF(act,ind+1)
    kbarfcstval = (1.0 - wgt) * KprimeVF(act,ind) + wgt * KprimeVF(act,ind+1)
    
    !store forecasts, note timing difference for (p,K')
    pfcstsim(t) = pfcstval
    kbarfcstsim(t+1) = kbarfcstval
    
    !start apparatus for solution on (p,K'), at initial guesses from SS
    pval = pfcstval
    kbarprimeval = kbarfcstval
        
    !insert this period's distributional info (moments, coefficient, and integration weights), into global arrays
    mommat = momSIM(:,:,t)
    rhomat = rhoSIM(:,:,t)
    intvec = intSIM(:,t)
    
    
        do fixediter=1,maxfixedit

        wval = phi / pval !available inside RHS
        
        !what are linear interpolation values for next period VF in agg K'?
        call kbarprimelinterp(kbarprimeval,kbarfcstind,kbarfcstwgt,kbar0,kbarnum)
        
        !this loop obtains the policies k'(z), as well as (part of) the RHS at these values
        do zct=1,znum
            RHSzct = zct
            RHSvec(zct) = -1.0*brent(k0(1),k0(knum/2),k0(knum),fnRHS,poltol,kprimeval)
            kprimevec(zct) = kprimeval
        end do !zct
            
        !construct the policies by calling subroutine (only constructs them for RHSact,RHSkbarct)
        call constructPolicySimp(kprimevec,xistarsimp,kprimesimp)
    
        !given policies at integration points, can construct implied K',C, which requires moment, coeff,
        ! and integral arrays
        call impliedaggs(xistarsimp,kprimesimp,Cactual,Kprimeactual,Yactual,Iactual,Nactual)

        !compute FPI errors on agg capital and price
        perror = (1.0/pval - Cactual) 
        Kprimeerror = (Kprimeactual - kbarprimeval)/abs(kbarprimeval)
        fixederror = maxval(abs((/perror,Kprimeerror/)))
        
        !exit FPI if converged
        if (fixederror < fixedtol ) exit
        
        !if not converged, update
        pval = pval + pfixeddamp * perror
        kbarprimeval = kbarprimeval + Kfixeddamp * Kprimeerror
        
        end do !fixediter
        write(*,*) " "
        write(*,"(A,I4,A,F4.2,A,F4.2,A,F7.5,A,I4)") "T = ",t," K' = ",kbarprimeval,", p = ",pval,&
            " err = ",fixederror,", it = ",fixediter-1
            
        write(13,"(A,I4,A,F4.2,A,F4.2,A,F7.5,A,I4)") "T = ",t," K' = ",kbarprimeval,", p = ",pval,&
            " err = ",fixederror,", it = ",fixediter-1            
        
    !close apparatus for fixed point on (p,K')
    
    !store price and agg capital next period associated with (a,K)
    psim(t) = pval
    kbarprimeval = max(kbarmin,min(kbarmax,kbarprimeval))
    kbarsim(t+1) = kbarprimeval    
    ysim(t) = Yactual
    isim(t) = Iactual
    nsim(t) = Nactual
    
    !DH stat series for p and K' 
    if (t==numdiscard) then
    
        pDHsim(t) = pval
        kbarDHsim(t) = kbarval
        
        call kbarprimelinterp(kbarDHsim(t),ind,wgt,kbar0,kbarnum)        
        kbarDHsim(t+1) = (1.0 - wgt) * KprimeVF(act,ind) + wgt * KprimeVF(act,ind+1)

    else if (t>numdiscard) then        
    
        call kbarprimelinterp(kbarDHsim(t),ind,wgt,kbar0,kbarnum)        
        pDHsim(t) = (1.0 - wgt) * pVF(act,ind) + wgt * pVF(act,ind+1)
        kbarDHsim(t+1) = (1.0 - wgt) * KprimeVF(act,ind) + wgt * KprimeVF(act,ind+1)
        
    end if
    
    
    !now, compute MICRO moments, which needs to occur with this period's rhomat, mommat, intvec, etc.   
     !MICROsim(nummicro,numper)
    ! 1 = i/k
    ! 2 = stdev(i/k)
    ! 3 = P(inaction)
    ! 4 = P(i/k>=0.2)
    ! 5 = P(i/k<=-0.2)
    ! 6 = P(i/k > 0)
    ! 7 = P(i/k < 0)
    
    do zct=1,znum
    do kct=1,nsimp+1
        
        zval = z0(zct);
        kval = simpnodes(kct);
        kprimeval = kprimesimp(zct,kct)
        ival = kprimeval - (1.0-delta)*kval !investment conditional upon investment
        xival = xistarsimp(zct,kct)
        padjust = cdfxi(xival)
        
        !investment rate
        MICROsim(1,t) = MICROsim(1,t) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust * (ival / kval)
        
        !investment rate squared - for stdev construction
        MICROsim(2,t) = MICROsim(2,t) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust * ( (ival / kval)**2.0 )
        
        !P(inaction)
        MICROsim(3,t) = MICROsim(3,t) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * (1.0 - padjust)
        
        !P(pos. spike)
        if ((ival/kval)>=0.2) then
            MICROsim(4,t) = MICROsim(4,t) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust
        end if
        
        !P(neg. spike)
        if ((ival/kval)<=-0.2) then
            MICROsim(5,t) = MICROsim(5,t) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust
        end if
        
        !P(pos. invest)
        if ((ival/kval)>0.0) then
            MICROsim(6,t) = MICROsim(6,t) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust
        end if
        
        !P(neg. invest)
        if ((ival/kval)<0.0) then
            MICROsim(7,t) = MICROsim(7,t) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust 
        end if
        
    
    end do !kct
    end do !zct
        
        
    !now, convert squared investment moment to stdev
    MICROsim(2,t) = sqrt(MICROsim(2,t) - (MICROsim(1,t)**2.0))
        
	!now, store the distribution in distkzsim
	open(8,file="distkzsim.txt",position='APPEND')
	do zct=1,znum
	do kct=1,kdensenum
		kval = kdense0(kct)
		write(8,*) Fk(kval,zct)
	end do !kct
	end do !zct
	close(8)
	
    !now that you have the price today and agg capital tomorrow, need to compute all of the moments for tomorrow
    call nextperiodmoments(xistarsimp,kprimesimp,momSIM(:,:,t+1))
            
    !now, find coefficients and integral weights for tomorrow's moments
    mommat = momSIM(:,:,t+1)    
    
    if (OPTMETHOD==1) then
        
        call findrhoBroyden()
        
    else if (OPTMETHOD==2) then
    
        call pso(rhovec,Pintval,Pintrhovec,lb,ub,nvar,npart,xtol,xquicktol,xquicknum,ftol,maxpsoit,phipso)
        gradvec = gradPint()
        write(*,*) "max grad error = ",maxval(abs(gradvec))
        
        !store coefficients and integral array
        ct=0
        do zct=1,znum
        do momct=1,momuse
            ct=ct+1
            rhomat(zct,momct) = rhovec(ct)
        end do !momct
        end do !zct
        
    else if (OPTMETHOD==3) then
        
        !call findrhoSR1()
        call findrhoSR1wolfe(1)
        
    end if
    
    rhoSIM(:,:,t+1) = rhomat
    intSIM(:,t+1) = intvec
            
    !now, record this period's distribution as a density
    do zct=1,znum
    do kct=1,kdensenum
        kval = kdense0(kct)
        densitySIM(zct,kct,t) = (ergdistz(zct) * Fk(kval,zct)) / intSIM(zct,t)
    end do !kct
    end do !zct
        
end do !t

write(*,*) "Done with simulation for ref loop ",refct," at ",omp_get_wtime()-start," seconds."
write(13,*) "Done with simulation for ref loop ",refct," at ",omp_get_wtime()-start," seconds."

!!!!!END SIMULATION

open(8,file="momSIM.txt")
do t=1,numper
do momct=1,momuse
do zct=1,znum
write(8,*) momSIM(zct,momct,t)
end do !zct
end do !momct
end do !t
close(8)

open(8,file="rhoSIM.txt")
do t=1,numper
do momct=1,momuse
do zct=1,znum
write(8,*) rhoSIM(zct,momct,t)
end do !zct
end do !momct
end do !t
close(8)

open(8,file="intSIM.txt")
do t=1,numper
do zct=1,znum
write(8,*) intSIM(zct,t)
end do !zct
end do !t
close(8)

open(8,file="densitySIM.txt")
do t=1,numper
do zct=1,znum
do kct=1,kdensenum
write(8,*) densitySIM(zct,kct,t)
end do 
end do 
end do 
close(8)

open(8,file="psim.txt")
do t=1,numper
write(8,*) psim(t)
end do !t
close(8)

open(8,file="pfcstsim.txt")
do t=1,numper
write(8,*) pfcstsim(t)
end do !t
close(8)

open(8,file="kbarfcstsim.txt")
do t=1,numper
write(8,*) kbarfcstsim(t)
end do !t
close(8)

open(8,file="pDHsim.txt")
do t=1,numper
write(8,*) pDHsim(t)
end do !t
close(8)

open(8,file="kbarDHsim.txt")
do t=1,numper
write(8,*) kbarDHsim(t)
end do !t
close(8)

open(8,file="kbarsim.txt")
do t=1,numper
write(8,*) kbarsim(t)
end do !t
close(8)

open(8,file="ysim.txt")
do t=1,numper
write(8,*) ysim(t)
end do !t
close(8)

open(8,file="isim.txt")
do t=1,numper
write(8,*) isim(t)
end do !t
close(8)

open(8,file="nsim.txt")
do t=1,numper
write(8,*) nsim(t)
end do !t
close(8)


open(8,file="MICROsim.txt")
do t=1,numper
do momct=1,nummicro
write(8,*) MICROsim(momct,t)
end do !momct
end do !t
close(8)

!!!!!NOW, DO REFERENCE MOMENT UPDATE, STORAGE, AND ERROR CALCULATION
write(*,*) " "
write(13,*) " "

write(*,*) "Now, updating reference moments based on simulation."
write(13,*) "Now, updating reference moments based on simulation."

!first, store old reference moments
if (refct==1) momREFstore(:,:,:,:) = 0.0
momREFstore(:,:,:,refct) = momREF

!now, update
momREFnew(:,:,:) = 0.0

!loop over agg prod states
do act=1,anum

    ct = 0
    do t=numdiscard,numper-1
        if (asimpos(t)==act) then
            ct=ct+1
            momREFnew(:,:,act) = momREFnew(:,:,act) + momSIM(:,:,t)
        end if
    end do !t
    
    !now compute mean from sum
    momREFnew(:,:,act) = (1.0/dble(ct)) * momREFnew(:,:,act)

end do !act

!what is the error in the reference moments?
momREFerror = maxval(abs(momREFnew - momREF))


!now, output the old and new moments
write(*,*) " "
write(13,*) " "

write(*,*) "   Old Moments, New Moments"
write(13,*) "   Old Moments, New Moments"
do act=1,anum
write(*,*) "Agg Prod Pos = ",act
write(13,*) "Agg Prod Pos = ",act
do zct=1,znum
do momct=1,momuse
    write(*,*) momREF(zct,momct,act)," ",momREFnew(zct,momct,act)
    write(13,*) momREF(zct,momct,act)," ",momREFnew(zct,momct,act)
end do !momct
end do !zct
end do !act

write(*,*) "Error for moment update in refct ",refct," is ",momREFerror,"."
write(13,*) "Error for moment update in refct ",refct," is ",momREFerror,"."

!if reference moments have converged, exit
if (momREFerror<reftol.or.refct==maxrefit) exit

!if reference moments have not yet converged, then update
momREF = momREF + refdamp * (momREFnew - momREF)

!now, compute new coefficients and integrals, then restart
write(*,*) " "
write(13,*) " "

do act=1,anum
    write(*,*) "Finding coeffs and ints for ref moments of act ",act,"."
    write(13,*) "Finding coeffs and ints for ref moments of act ",act,"."

    mommat = momREF(:,:,act)
    call findrhoSR1wolfe(1)
    rhoREF(:,:,act) = rhomat
    intREF(:,act) = intvec
    
end do !act

write(*,*) " "
write(13,*) " "

end do !refct
!END OF BIG MODEL SOLUTION LOOP


!!!!!!DO IRF SIMULATION
if (doIRF==1) then 
write(*,*) "Done with solution at ",omp_get_wtime() - start," seconds. Now doing IRF."
write(13,*) "Done with solution at ",omp_get_wtime() - start," seconds.  Now doing IRF."
!first, need to initialize the moments
!(t,simct,shockct)
kbarsimIRF(:,:,:) = 0.0
kbarval = 1.55 !what is agg capital at assumed starting point?
kbarsimIRF(1,:,:) = kbarval

!store moments, coefficients, and integrals
!(zct,momct,t,simct,shockct) for mom and rho, (zct,t,simct,shockct) for int
do simct=1,numsimIRF
do shockct=1,2
momIRF(:,:,1,simct,shockct) = momREF(:,:,asimposIRF(1,simct,shockct))
momIRF(:,1,1,simct,shockct) = momREF(:,1,asimposIRF(1,simct,shockct)) * &
    (kbarval / sum(ergdistz * momREF(:,1,asimposIRF(1,simct,shockct))))
intIRF(:,1,simct,shockct) = intREF(:,asimposIRF(1,simct,shockct))
rhoIRF(:,:,1,simct,shockct) = rhoREF(:,:,asimposIRF(1,simct,shockct))
end do 
end do 

psimIRF(:,:,:) = 0.0
ysimIRF(:,:,:) = 0.0
isimIRF(:,:,:) = 0.0
nsimIRF(:,:,:) = 0.0

do shockct=1,2
do simct=1,numsimIRF
do t=1,numperIRF-1
    
    !for this period, what is the aggregate productivity?
    act = asimposIRF(t,simct,shockct)
    aval = a0(act)
    RHSact = act
    
    kbarval = kbarsimIRF(t,simct,shockct)
    
    !obtain forecasted values of (p,K')        
    call kbarprimelinterp(kbarval,ind,wgt,kbar0,kbarnum)        
    pfcstval = (1.0 - wgt) * pVF(act,ind) + wgt * pVF(act,ind+1)
    kbarfcstval = (1.0 - wgt) * KprimeVF(act,ind) + wgt * KprimeVF(act,ind+1)
            
    !start apparatus for solution on (p,K'), at initial guesses from SS
    pval = pfcstval
    kbarprimeval = kbarfcstval
        
    !insert this period's distributional info (moments, coefficient, and integration weights), into global arrays
    mommat(:,:) = 0.0
    rhomat(:,:) = 0.0
    intvec(:) = 0.0
    mommat = momIRF(:,:,t,simct,shockct)
    rhomat = rhoIRF(:,:,t,simct,shockct)
    intvec = intIRF(:,t,simct,shockct)
    
        do fixediter=1,maxfixedit

        wval = phi / pval !available inside RHS
        
        !what are linear interpolation values for next period VF in agg K'?
        call kbarprimelinterp(kbarprimeval,kbarfcstind,kbarfcstwgt,kbar0,kbarnum)
        
        !this loop obtains the policies k'(z), as well as (part of) the RHS at these values
        do zct=1,znum
            RHSzct = zct
            RHSvec(zct) = -1.0*brent(k0(1),k0(knum/2),k0(knum),fnRHS,poltol,kprimeval)
            kprimevec(zct) = kprimeval
        end do !zct
            
        !construct the policies by calling subroutine (only constructs them for RHSact,RHSkbarct)
        call constructPolicySimp(kprimevec,xistarsimp,kprimesimp)
    
        !given policies at integration points, can construct implied K',C, which requires moment, coeff,
        ! and integral arrays
        call impliedaggs(xistarsimp,kprimesimp,Cactual,Kprimeactual,Yactual,Iactual,Nactual)

        !compute FPI errors on agg capital and price
        perror = (1.0/pval - Cactual) 
        Kprimeerror = (Kprimeactual - kbarprimeval)/abs(kbarprimeval)
        fixederror = maxval(abs((/perror,Kprimeerror/)))
        
        !exit FPI if converged
        if (fixederror < fixedtol ) exit
        
        !if not converged, update
        pval = pval + pfixeddamp * perror
        kbarprimeval = kbarprimeval + Kfixeddamp * Kprimeerror
    
        end do !fixediter
        
        if (t==1) then
            write(*,*) " "
            write(*,"(A,I4,A,I5,A,I3,A,F7.5,A,I4)") "T = ",t," sim = ",simct," shock = ",shockct," err = ",&
                fixederror," it = ",fixediter-1
            
            write(13,"(A,I4,A,I5,A,I3,A,F7.5,A,I4)") "T = ",t," sim = ",simct," shock = ",shockct," err = ",&
                fixederror," it = ",fixediter-1
        end if
    !close apparatus for fixed point on (p,K')
    
    !store price and agg capital next period associated with (a,K)
    psimIRF(t,simct,shockct) = pval
    kbarprimeval = max(kbarmin,min(kbarmax,kbarprimeval))
    kbarsimIRF(t+1,simct,shockct) = kbarprimeval    
    ysimIRF(t,simct,shockct) = Yactual
    isimIRF(t,simct,shockct) = Iactual
    nsimIRF(t,simct,shockct) = Nactual
    
    !now that you have the price today and agg capital tomorrow, need to compute all of the moments for tomorrow
    call nextperiodmoments(xistarsimp,kprimesimp,momIRF(:,:,t+1,simct,shockct))
            
    !now, find coefficients and integral weights for tomorrow's moments
    mommat = momIRF(:,:,t+1,simct,shockct)    
    
    if (OPTMETHOD==1) then
        
        call findrhoBroyden()
        
    else if (OPTMETHOD==2) then
    
        call pso(rhovec,Pintval,Pintrhovec,lb,ub,nvar,npart,xtol,xquicktol,xquicknum,ftol,maxpsoit,phipso)
        gradvec = gradPint()
        write(*,*) "max grad error = ",maxval(abs(gradvec))
        
        !store coefficients and integral array
        ct=0
        do zct=1,znum
        do momct=1,momuse
            ct=ct+1
            rhomat(zct,momct) = rhovec(ct)
        end do !momct
        end do !zct
        
    else if (OPTMETHOD==3) then
        
        !call findrhoSR1()
        call findrhoSR1wolfe(1)
        
    end if
    
    rhoIRF(:,:,t+1,simct,shockct) = rhomat
    intIRF(:,t+1,simct,shockct) = intvec
                    
        
end do !t
end do !simct
end do !shockct

write(*,*) "Done with IRF sims at ",omp_get_wtime() - start," seconds."
write(13,*) "Done with IRF sims at ",omp_get_wtime() - start," seconds."

!now, output IRF data
open(8,file="ysimIRF.txt")
do shockct=1,2
do simct=1,numsimIRF
do t=1,numperIRF
write(8,*) ysimIRF(t,simct,shockct)    
end do
end do 
end do
close(8)

open(8,file="isimIRF.txt")
do shockct=1,2
do simct=1,numsimIRF
do t=1,numperIRF
write(8,*) isimIRF(t,simct,shockct)    
end do
end do 
end do
close(8)

open(8,file="nsimIRF.txt")
do shockct=1,2
do simct=1,numsimIRF
do t=1,numperIRF
write(8,*) nsimIRF(t,simct,shockct)    
end do
end do 
end do
close(8)

open(8,file="psimIRF.txt")
do shockct=1,2
do simct=1,numsimIRF
do t=1,numperIRF
write(8,*) psimIRF(t,simct,shockct)    
end do
end do 
end do
close(8)

open(8,file="kbarsimIRF.txt")
do shockct=1,2
do simct=1,numsimIRF
do t=1,numperIRF
write(8,*) kbarsimIRF(t,simct,shockct)    
end do
end do 
end do
close(8)

end if !doIRF
!!!!!!!END IRF SIMULATION

finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished all calculations in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished all calculations ",finish-start," seconds."
close(13); !closing log file

contains

subroutine discretize_simulate()
implicit none

integer :: zct,act,statect,stateprimect,ct,t

double precision :: asimgrid(anum)

!call dimension of random seed
call random_seed(size=seeddim)
allocate(seedarray(seeddim))
!set up idio capital grids
call linspace(k0,log(kmin),log(kmax),knum);
k0=exp(k0);

call linspace(kdense0,log(kmin),log(kmax),kdensenum);
kdense0=exp(kdense0);

!do discretization of the z process to set up z0 grid and transitions
call tauchen(znum,rhoz,sigmaz,nstdevz,pr_mat_z,z0)

!do discretization of the a process to set up a0 grid and transitions
call tauchen(anum,rhoa,sigmaa,nstdeva,pr_mat_a,a0)

!create unified transition matrix
pr_mat(:,:) = 0.0
do zct=1,znum
do act=1,anum
    statect = (act-1)*znum+zct
    do zprimect=1,znum
    do aprimect=1,anum
        stateprimect=(aprimect-1)*znum+zprimect
        pr_mat(statect,stateprimect) = pr_mat_z(zct,zprimect)*pr_mat_a(act,aprimect)
    end do !aprimect
    end do !zprimect
    pr_mat(statect,:) = pr_mat(statect,:)/sum(pr_mat(statect,:))
end do !act
end do !zct

!actually seed the RNG
do ct=1,seeddim
    seedarray(ct) = seedint + ct
end do !ct
call random_seed(put=seedarray)

!draw the random shocks for the simulation and perform simulation of agg productivity
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

do ct=1,maxdistit
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

!!!!!!!!FIND ERGODIC DISTRIBUTION OF IDIO PROD
ergdistz(:) = 0.0
ergdistz(znum/2) = 1.0
ergoldz = ergdistz

do ct=1,maxdistit
    
    !iterate micro dist forward
    ergnewz(:) = 0.0
    do zct=1,znum
    do zprimect=1,znum
        ergnewz(zprimect) = ergnewz(zprimect) + pr_mat_z(zct,zprimect)*ergoldz(zct)
    end do !zprimect
    end do !zct
    
    !check error condition and exit if needed
    ergerrz = maxval(abs(ergnewz-ergoldz))
    if (ergerrz<disterrortol) exit
    !if not converged, update and round
    ergoldz = ergnewz
    ergoldz = ergoldz / sum(ergoldz)
end do !ct
ergdistz = ergnewz

end subroutine discretize_simulate

subroutine nextperiodmoments(xistarsimp,kprimesimp,momnextmat)

!given a policy function, and this period's moments in mommat, rhomat, and intvec, 
!this function returns next period's moments in "momnextmat"

double precision, intent(in) :: xistarsimp(znum,nsimp+1),kprimesimp(znum,nsimp+1)
double precision, intent(out) :: momnextmat(znum,momuse)

integer :: zct,kct,zprimect,momct
double precision :: kval,adjustmat(znum,nsimp+1,2),kprimeval,xival,padjust

!construct adjustmat
!adjustmat(zct,kct,1) = this period weight on z0(zct),simpnodes(kct),and non adjustment
!adjustmat(zct,kct,2) = this period weight on z0(zct),simpnodes(kct),and adjustment

adjustmat(:,:,:) = 0.0
do zct=1,znum
do kct=1,nsimp+1
    
    !what are states today as well as policy?
    kval = simpnodes(kct)
    kprimeval = kprimesimp(zct,kct)
    xival = xistarsimp(zct,kct)
    padjust = cdfxi(xival)
    
    adjustmat(zct,kct,1) = (1.0-padjust) *  ( ergdistz(zct) * simpweights(kct) * Fk(kval,zct) ) / intvec(zct)
    adjustmat(zct,kct,2) = padjust *  ( ergdistz(zct) * simpweights(kct) * Fk(kval,zct) ) / intvec(zct)
    
end do !kct
end do !zct


!now, fill in the first moments for each zprimect
momnextmat(:,:) = 0.0
do zprimect=1,znum
    
    do zct=1,znum
    do kct=1,nsimp+1
        kval = simpnodes(kct)
        kprimeval = kprimesimp(zct,kct)
        momnextmat(zprimect,1) = momnextmat(zprimect,1) + ( pr_mat_z(zct,zprimect) / ergdistz(zprimect) ) * &
            ( adjustmat(zct,kct,1) * max((1.0-delta)*kval,kmin) + adjustmat(zct,kct,2) * kprimeval )
    end do !kct
    end do !zct
    
end do !zprimect

!now, fill in the higher moments for each zprimect
do momct=2,momuse
    do zprimect=1,znum
        
        do zct=1,znum
        do kct=1,nsimp+1
            
            kval = simpnodes(kct)
            kprimeval = kprimesimp(zct,kct)
            
            momnextmat(zprimect,momct) = momnextmat(zprimect,momct) + ( pr_mat_z(zct,zprimect) / ergdistz(zprimect) ) * &
                ( adjustmat(zct,kct,1) * ( max((1.0-delta)*kval,kmin) -  momnextmat(zprimect,1) )**dble(momct) +  &
                adjustmat(zct,kct,2) * ( kprimeval - momnextmat(zprimect,1)  )**dble(momct) )
            
        end do !kct
        end do !zct
        
    end do !zprimect
end do !momct

return    
end subroutine nextperiodmoments

subroutine impliedaggs(xistarsimp,kprimesimp,Cactual,Kprimeactual,Yactual,Iactual,Nactual)
    
!need kbarval, RHSact, RHSkbarct, which are global, then
!can use rhomat, mommat, and intvec to construct implied moments

double precision, intent(in) :: xistarsimp(znum,nsimp+1),&
    kprimesimp(znum,nsimp+1) 
double precision, intent(out) :: Cactual,Kprimeactual,Yactual,Iactual,Nactual

integer :: zct,kct
double precision :: Kref,Kshift,wgt,&
    kval,zval,xival,padjust,kprimeval

!what is the implied capital from reference moments?
Kref = sum(ergdistz * mommat(:,1))
Kshift = kbarval / Kref !note that kbarval is available to the subroutine

!modify mommat by multiplication with this value    
mommat(:,1) = Kshift * mommat(:,1)

Yactual = 0.0
Iactual = 0.0    
Kprimeactual = 0.0
Nactual = 0.0
do zct=1,znum    
do kct=1,nsimp+1
    
    !what are states here?
    zval = z0(zct)
    kval = simpnodes(kct)
    
    !what are adjustment probs and capital policies?
    xival = xistarsimp(zct,kct)
    padjust = cdfxi(xival)
    kprimeval = kprimesimp(zct,kct)
    
    !what is the weight at this point, for integration against density?
    wgt = ( ergdistz(zct) * simpweights(kct) * Fk(kval,zct) ) / intvec(zct)
    
    !what are aggregates?
    Yactual = Yactual + yreduced(zval,kval) * wgt
    Iactual = Iactual + padjust*(kprimeval - (1.0-delta)*kval) * wgt
    Kprimeactual = Kprimeactual + ( padjust * kprimeval + (1.0 - padjust)*(1.0-delta)*kval ) * wgt
    Nactual = Nactual + (nreduced(zval,kval) + expecxi(xival) ) * wgt
end do !zct
end do !kct

Cactual = Yactual - Iactual
    
end subroutine impliedaggs

subroutine constructBellman(kprimevec)
implicit none

!note that Va,Vna,V,and Vold are all computable from here
double precision, intent(in) :: kprimevec(znum)

integer :: zct,kct,zprimect,aprimect
double precision :: Vaval,Vnaval,kprimeval,Vnextval,zval,kval,xival,Vval

do zct=1,znum
do kct=1,knum
    
    RHSzct = zct
    
    !what states are there here?
    zval = z0(zct)
    kval = k0(kct)
    
    !construct value when not adjusting
    Vnaval = pval * freduced(zval,kval)
    kprimeval = max((1.0-delta)*kval,kmin) !no adjustment in capital
    
    !add continuation value when not adjusting
    do zprimect=1,znum
    do aprimect=1,anum

        !evaluate continuation at kbarfcstind
        call splint(k0,Vold(zprimect,:,aprimect,kbarfcstind),V2old(zprimect,:,aprimect,kbarfcstind),knum,kprimeval,Vnextval)
        Vnaval = Vnaval + beta * pr_mat_z(RHSzct,zprimect) * pr_mat_a(RHSact,aprimect) * (1.0 - kbarfcstwgt) *  Vnextval   
        
        !then, evaluate continuation at kbarfcstind+1
        call splint(k0,Vold(zprimect,:,aprimect,kbarfcstind+1),V2old(zprimect,:,aprimect,kbarfcstind+1),knum,kprimeval,Vnextval)
        Vnaval = Vnaval + beta * pr_mat_z(RHSzct,zprimect) * pr_mat_a(RHSact,aprimect) * kbarfcstwgt *  Vnextval   
        
    end do !aprimect
    end do !zprimect
    
    !construct value when adjusting
    kprimeval = kprimevec(zct)
    Vaval = pval*freduced(zval,kval)!current period return
    Vaval = Vaval - pval * ( kprimeval - ( 1.0 - delta ) * kval )
    
    !add continuation value when not adjusting
    do zprimect=1,znum
    do aprimect=1,anum

        !evaluate continuation at kbarfcstind
        call splint(k0,Vold(zprimect,:,aprimect,kbarfcstind),V2old(zprimect,:,aprimect,kbarfcstind),knum,kprimeval,Vnextval)
        Vaval = Vaval + beta * pr_mat_z(RHSzct,zprimect) * pr_mat_a(RHSact,aprimect) * (1.0 - kbarfcstwgt) *  Vnextval   
        
        !then, evaluate continuation at kbarfcstind+1
        call splint(k0,Vold(zprimect,:,aprimect,kbarfcstind+1),V2old(zprimect,:,aprimect,kbarfcstind+1),knum,kprimeval,Vnextval)
        Vaval = Vaval + beta * pr_mat_z(RHSzct,zprimect) * pr_mat_a(RHSact,aprimect) * kbarfcstwgt *  Vnextval   
        
    end do !aprimect
    end do !zprimect        
    
    !now that the value when adjusting and when not adjusting are determined, create and store total value
    xival = (Vaval-Vnaval)/phi  !cutoff AC value
    
    !what is the implied total value?
    Vval = -1.0*phi*expecxi(xival) + cdfxi(xival) * Vaval + (1.0 - cdfxi(xival))*Vnaval
    
    xistar(zct,kct,RHSact,RHSkbarct) = xival
    kprime(zct,kct,RHSact,RHSkbarct) = kprimevec(zct)
    Va(zct,kct,RHSact,RHSkbarct) = Vaval
    Vna(zct,kct,RHSact,RHSkbarct) = Vnaval
    V(zct,kct,RHSact,RHSkbarct) = Vval
   
end do !kct
end do !zct


end subroutine constructBellman

subroutine constructPolicySimp(kprimevec,xistarsimp,kprimesimp)
implicit none

!xistarsimp(znum,nsimp+1,anum,kbarnum),kprimesimp(znum,nsimp+1,anum,kbarnum)
double precision, intent(in) :: kprimevec(znum)
double precision, intent(inout) :: xistarsimp(znum,nsimp+1),&
    kprimesimp(znum,nsimp+1)
    
integer :: zct,kct,zprimect,aprimect
double precision :: Vaval,Vnaval,kprimeval,Vnextval,zval,kval,xival

do zct=1,znum
do kct=1,nsimp+1
    
    RHSzct = zct
    
    !what states are there here?
    zval = z0(zct)
    kval = simpnodes(kct)
    
    !construct value when not adjusting
    Vnaval = pval * freduced(zval,kval)
    kprimeval = max((1.0-delta)*kval,kmin) !no adjustment in capital
    
    
    !add continuation value when not adjusting
    do zprimect=1,znum
    do aprimect=1,anum

        !evaluate continuation at kbarfcstind
        call splint(k0,Vold(zprimect,:,aprimect,kbarfcstind),V2old(zprimect,:,aprimect,kbarfcstind),knum,kprimeval,Vnextval)
        Vnaval = Vnaval + beta * pr_mat_z(RHSzct,zprimect) * pr_mat_a(RHSact,aprimect) * (1.0 - kbarfcstwgt) *  Vnextval   
        
        !then, evaluate continuation at kbarfcstind+1
        call splint(k0,Vold(zprimect,:,aprimect,kbarfcstind+1),V2old(zprimect,:,aprimect,kbarfcstind+1),knum,kprimeval,Vnextval)
        Vnaval = Vnaval + beta * pr_mat_z(RHSzct,zprimect) * pr_mat_a(RHSact,aprimect) * kbarfcstwgt *  Vnextval   
        
    end do !aprimect
    end do !zprimect
    
    !construct value when adjusting
    kprimeval = kprimevec(zct)
    Vaval = pval*freduced(zval,kval)!current period return
    Vaval = Vaval - pval * ( kprimeval - ( 1.0 - delta ) * kval )
    
        
    !add continuation value when not adjusting
    do zprimect=1,znum
    do aprimect=1,anum

        !evaluate continuation at kbarfcstind
        call splint(k0,Vold(zprimect,:,aprimect,kbarfcstind),V2old(zprimect,:,aprimect,kbarfcstind),knum,kprimeval,Vnextval)
        Vaval = Vaval + beta * pr_mat_z(RHSzct,zprimect) * pr_mat_a(RHSact,aprimect) * (1.0 - kbarfcstwgt) *  Vnextval   
        
        !then, evaluate continuation at kbarfcstind+1
        call splint(k0,Vold(zprimect,:,aprimect,kbarfcstind+1),V2old(zprimect,:,aprimect,kbarfcstind+1),knum,kprimeval,Vnextval)
        Vaval = Vaval + beta * pr_mat_z(RHSzct,zprimect) * pr_mat_a(RHSact,aprimect) * kbarfcstwgt *  Vnextval   
        
    end do !aprimect
    end do !zprimect        
    
    !now that the value when adjusting and when not adjusting are determined, create and store total value
    xival = (Vaval-Vnaval)/phi  !cutoff AC value
    xistarsimp(zct,kct) = xival
    kprimesimp(zct,kct) = kprimevec(zct)
    
end do !kct
end do !zct
    
    
end subroutine constructPolicySimp

subroutine kbarprimelinterp(kbarprimeval,ind,wgt,kbar0,kbarnum)
implicit none

!goes from K' --> (wgt,ind) for linear interpolation on agg capital state in continuation VF

integer, intent(in) :: kbarnum
double precision, intent(in) :: kbarprimeval
double precision, intent(in) :: kbar0(kbarnum)
double precision, intent(out) :: wgt
integer, intent(out) :: ind

ind = kbarnum/2
call hunt(kbar0,kbarnum,kbarprimeval,ind); !kbarprimeval is in interval kbarfcstind
if (ind<=0) then 
    wgt = 0.0
    ind=1
else if (ind>=1.and.ind<=(kbarnum-1)) then
    wgt = (kbarprimeval - kbar0(ind))/( kbar0(ind+1)-kbar0(ind) )
else if (ind>=kbarnum) then
    wgt = 1.0
    ind = kbarnum-1
end if

end subroutine

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
if(xi<=xibar.and.xi>=0.0)expecxi = ( xi ** 2.0 ) / (2.0 * xibar)
if(xi>xibar)expecxi = ( xibar ** 2.0 ) / (2.0 * xibar)
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

double precision function fnRHS(kprimeval)
implicit none

!this function returns -1 * ( -p * k' + beta * E V(z',k',a',K') ), for minimization wrt k'

double precision :: kprimeval,Vnextval

integer :: zprimect,aprimect

!note those things available to the function:
!k0, Vold, V2old, knum, anum, znum, pr_mat, kbarfcstind,kbarfcstwgt (all threads)
!RHSzct, RHSact (one thread only via "threadprivate" status)

!current period return (the only terms relevant to k optimization)
fnRHS = -pval * kprimeval

!add the continuation value
do zprimect=1,znum
do aprimect=1,anum    
    
    
    !evaluate continuation at kbarfcstind
    call splint(k0,Vold(zprimect,:,aprimect,kbarfcstind),V2old(zprimect,:,aprimect,kbarfcstind),knum,kprimeval,Vnextval)
    fnRHS = fnRHS + beta * pr_mat_z(RHSzct,zprimect) * pr_mat_a(RHSact,aprimect) * (1.0 - kbarfcstwgt) *  Vnextval   
    
    !then, evaluate continuation at kbarfcstind+1
    call splint(k0,Vold(zprimect,:,aprimect,kbarfcstind+1),V2old(zprimect,:,aprimect,kbarfcstind+1),knum,kprimeval,Vnextval)
    fnRHS = fnRHS + beta * pr_mat_z(RHSzct,zprimect) * pr_mat_a(RHSact,aprimect) * kbarfcstwgt *  Vnextval   
    
end do !aprimect
end do !zprimect

fnRHS = -1.0 * fnRHS

end function fnRHS

subroutine convertpolicy(kprimedense,kprimedensewgt,kprimedenseind)
implicit none    

double precision, intent(in) :: kprimedense(znum,kdensenum)
double precision, intent(out) :: kprimedensewgt(znum,kdensenum)
integer, intent(out) :: kprimedenseind(znum,kdensenum)

integer :: zct,kct,ind
double precision :: wgt,kprimeval

!loop over states
do zct=1,znum
do kct=1,kdensenum
    
    !what is continuous policy when adjusting?
    kprimeval = kprimedense(zct,kct)
    
    !which interval is the policy in?
    ind = kct
    call hunt(kdense0,kdensenum,kprimeval,ind)
    
    !modify for endpoints if needed
    if (ind<=0) then
        wgt = 0.0
        ind = 1
    else if (ind>=1.and.ind<=(kdensenum-1)) then
        wgt = (kprimeval - kdense0(ind))/(kdense0(ind+1)-kdense0(ind))
    else if (ind>=kdensenum) then
        wgt = 1.0
        ind = kdensenum-1
    end if
    
    kprimedensewgt(zct,kct) = wgt
    kprimedenseind(zct,kct) = ind
    
end do !kct
end do !zct
    
end subroutine convertpolicy

double precision function Fk(kval,zct)
implicit none

!this function evaluates the density-proportional function
!Fk = exp(rho_1 * (k-m^zct_1)+....+rho_momuse * ((k - m^zct_1)^momuse - m^zct_momuse) )

double precision :: kval
integer :: zct

double precision :: rhovec(momuse),momvec(momuse)
integer :: momct

!extract correct value of rho from global variable
rhovec = rhomat(zct,:)
momvec = mommat(zct,:)

Fk = rhovec(1)*(kval - momvec(1))

do momct=2,momuse
    Fk = Fk + rhovec(momct) * ( (kval-momvec(1))**dble(momct) - momvec(momct) )
end do !momct

Fk = exp(Fk)

end function Fk

function Pint()
implicit none

!Pint is \sum_{z=1}^n_z \int_kmin^kmax Fk(k,z) dk, where Fk is above, i.e. Pint is the objective that
!is minimized in the process of finding distributions matching the indicated moments

double precision :: Pint

integer :: zct,kct
double precision :: kval,wgt,addval

Pint = 0.0
intvec(:) = 0.0
do zct=1,znum
    do kct=1,nsimp+1
        
        kval = simpnodes(kct)
        wgt = simpweights(kct)
        addval = wgt * Fk(kval,zct)
        
        Pint = Pint + addval
        intvec(zct) = intvec(zct) + addval
    end do !kct
end do !zct

end function Pint

function gradPint()
implicit none

!gradPint is the znum * momuse x 1 gradient of Pint() from above, wrt to {rho^z_1,....,rho^1_momuse }_{z=1}^znum
    
double precision :: gradPint(znum*momuse)    

integer :: momct,zct,kct,gradct
double precision :: rhovec(momuse),momvec(momuse),kval,wgt

gradPint(:) = 0.0

do zct=1,znum

!extract moments and rho's for zct    
rhovec = rhomat(zct,:)
momvec = mommat(zct,:)
    
do momct=1,momuse
    !track location in the gradient
    gradct = (zct-1)*momuse + momct
    
    !perform integration for each entry
    do kct=1,nsimp+1
        
        kval = simpnodes(kct)
        wgt = simpweights(kct)
        if (momct>1) then
            gradPint(gradct) = gradPint(gradct) + wgt * ( (kval - momvec(1))**dble(momct) - momvec(momct) ) * Fk(kval,zct)
        else if (momct==1) then
            gradPint(gradct) = gradPint(gradct) + wgt * (kval - momvec(1)) * Fk(kval,zct)
        end if
        
    end do !kct
    
    
end do !momct
end do !zct
    
end function gradPint

double precision function Pintrhovec(rhovec)
implicit none

double precision, intent(in) :: rhovec(znum*momuse)

integer :: zct,momct,ct
ct=0
do zct=1,znum
do momct=1,momct
    ct=ct+1
    rhomat(zct,momct) = rhovec(ct)
end do !momct
end do !zct

Pintrhovec = Pint()
    
end function Pintrhovec

subroutine findrhoBroyden()
implicit none

double precision :: newJinv(znum*momuse,znum*momuse),oldJinv(znum*momuse,znum*momuse),oldbigrho(znum*momuse),&
    newbigrho(znum*momuse),oldgrad(znum*momuse),newgrad(znum*momuse),diffrho(znum*momuse),diffgrad(znum*momuse),&
    numval(znum*momuse),denomval,multval(znum*momuse)
double precision :: rhoerror,graderror,funcerror,oldfunc,newfunc
integer :: zct,momct,gradct,iter,ct1,ct2
   
!make initial guesses for rho and evaluate gradient at those guesses

do zct=1,znum
do momct=1,momuse
    if (mod(momct,2)==0) then
        rhomat(zct,momct) = -1.0
    else
        rhomat(zct,momct) = 0.0
    end if
end do !momct
end do

do zct=1,znum
do momct=1,momuse
    gradct=(zct-1)*momuse + momct
    oldbigrho(gradct) = rhomat(zct,momct)
end do !zct
end do 

do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    rhomat(zct,momct) = oldbigrho(gradct)
end do !momct
end do !zct
oldgrad = gradPint()
oldfunc = Pint()

do zct=1,znum
do momct=1,momuse
    if (mod(momct,2)==0) then
        rhomat(zct,momct) = -5.0
    else
        rhomat(zct,momct) = 0.0
    end if
end do !momct
end do

do zct=1,znum
do momct=1,momuse
    gradct=(zct-1)*momuse + momct
    newbigrho(gradct) = rhomat(zct,momct)
end do !zct
end do 

do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    rhomat(zct,momct) = newbigrho(gradct)
end do !momct
end do !zct
newgrad = gradPint()  
newfunc = Pint()

!what are the implied initial differences in rho and the gradient?
diffrho = newbigrho - oldbigrho
diffgrad = newgrad - oldgrad

!make initial guesses for inv Hessian approx - identity matrix
oldJinv(:,:) = 0.0
do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    oldJinv(gradct,gradct) = 1.0
end do !momct
end do !zct

!actually do the Broyden iterations
do iter=1,maxbroydit
    
    !set some stuff to zero
    numval(:) = 0.0
    denomval = 0.0
    
    !determine numerator
    numval = matmul(oldJinv,diffgrad)
    numval = diffrho - numval

    !determine denominator
    multval = matmul(oldJinv,diffgrad)
    denomval = sum(diffrho*multval)
    
    !determinemultval
    multval(:) = 0.0
    multval = matmul(diffrho,oldJinv)
    
    !what is new guess for inv Hessian?
    do ct1=1,znum*momuse
    do ct2=1,znum*momuse
        newJinv(ct1,ct2) = oldJinv(ct1,ct2) + (1.0 / denomval) * numval(ct1)* multval(ct2)
    end do !ct1
    end do !ct2
    
    !what is new guess for rho?
    oldbigrho = newbigrho
    oldgrad = newgrad
    oldJinv = newJinv
    oldfunc = newfunc
    
    !take new step, at fixed stepsize length
    newbigrho = newbigrho - stepsize*matmul(newJinv,newgrad)
       
    !what is the new value for the gradient, i.e. evaluate gradient at newbigrho
    do zct=1,znum
    do momct=1,momuse
        gradct = (zct-1)*momuse + momct
        rhomat(zct,momct) = newbigrho(gradct)
    end do !momct
    end do !zct
    newgrad = gradPint()
    newfunc = Pint()

    rhoerror = maxval(abs(newbigrho - oldbigrho))/maxval(abs(oldbigrho))
    graderror = maxval(abs(newgrad))
    
    if (iter>5000) then
    funcerror = abs(1.0 - newfunc/oldfunc)
    else if (iter<=5000) then
    funcerror = 1.0
    end if

    if (mod(iter,1000)==1) then
        write(*,"(A,I5,A,F7.4,A,F7.4,A,F7.4,A,F7.4)") "it = ",iter,", rho = ",rhoerror,", grad = ",&
            graderror,", fun = ",funcerror,", val = ",Pint()
    end if
    
    if (iter>10) then
    if (graderror<broydengradtol.or.funcerror<broydenfunctol.or.rhoerror<broydenrhotol) exit
    end if
    
    diffgrad = newgrad - oldgrad
    diffrho = newbigrho - oldbigrho
    
end do !iter

write(*,"(A,I5,A,F7.4,A,F7.4,A,F7.4,A,F7.4)") "it = ",iter-1,", rho = ",rhoerror,", grad = ",&
            graderror,", fun = ",funcerror,", val = ",Pint()

end subroutine findrhoBroyden

subroutine findrhoSR1wolfe(initflag)
implicit none

double precision :: newJinv(znum*momuse,znum*momuse),oldJinv(znum*momuse,znum*momuse),oldbigrho(znum*momuse),&
    newbigrho(znum*momuse),oldgrad(znum*momuse),newgrad(znum*momuse),diffrho(znum*momuse),diffgrad(znum*momuse),&
    numval(znum*momuse),denomval,nummat(znum*momuse,znum*momuse),stepdirec(znum*momuse)
double precision :: rhoerror,graderror,funcerror,oldfunc,newfunc,stepval
integer :: zct,momct,gradct,iter,ct1,ct2,stepct,initflag
logical :: wolfe1,wolfe2,nannewgrad,nannewfunc

!make initial guesses for rho and evaluate gradient at those guesses

if (initflag==0) then
do zct=1,znum
do momct=1,momuse
    if (mod(momct,2)==0) then
        rhomat(zct,momct) = -1.0
    else
        rhomat(zct,momct) = 0.0
    end if
end do !momct
end do

else if (initflag==1) then

    rhomat = rhoSS

end if

do zct=1,znum
do momct=1,momuse
    gradct=(zct-1)*momuse + momct
    oldbigrho(gradct) = rhomat(zct,momct)
end do !zct
end do 

do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    rhomat(zct,momct) = oldbigrho(gradct)
end do !momct
end do !zct
oldgrad = gradPint()
oldfunc = Pint()

if (initflag==0) then
do zct=1,znum
do momct=1,momuse
    if (mod(momct,2)==0) then
        rhomat(zct,momct) = -1.3
    else
        rhomat(zct,momct) = 0.0
    end if
end do !momct
end do

else if (initflag==1) then

    rhomat = rhoSS*1.1

end if

do zct=1,znum
do momct=1,momuse
    gradct=(zct-1)*momuse + momct
    newbigrho(gradct) = rhomat(zct,momct)
end do !zct
end do 

do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    rhomat(zct,momct) = newbigrho(gradct)
end do !momct
end do !zct
newgrad = gradPint()  
newfunc = Pint()

!what are the implied initial differences in rho and the gradient?
diffrho = newbigrho - oldbigrho
diffgrad = newgrad - oldgrad

!make initial guesses for inv Hessian approx - identity matrix
oldJinv(:,:) = 0.0
do zct=1,znum
do momct=1,momuse
    gradct = (zct-1)*momuse + momct
    oldJinv(gradct,gradct) = 1.0
end do !momct
end do !zct

!actually do the SR1 iterations
do iter=1,maxbroydit
    
    !set some stuff to zero
    numval(:) = 0.0
    nummat(:,:) = 0.0
    denomval = 0.0
    
    !determine numerator vector
    numval = matmul(oldJinv,diffgrad)
    numval = diffrho - numval
    
    !determine numerator matrix
    do ct1=1,znum*momuse
    do ct2=1,znum*momuse
        nummat(ct1,ct2) = numval(ct1) * numval(ct2)
    end do
    end do
    
    !determine denominator value
    do ct1=1,znum*momuse
        denomval = denomval + numval(ct1)*diffgrad(ct1)
    end do 
    
    !what is new guess for inv Hessian?
    newJinv = oldJinv + (1.0 / denomval) * nummat
    
    !what is new guess for rho?
    oldbigrho = newbigrho
    oldgrad = newgrad
    oldJinv = newJinv
    oldfunc = newfunc
    
    !choose a step length, according to Wolfe conditions
    stepdirec = matmul(newJinv,newgrad)
    stepval = stepinit
    do stepct=1,maxstepnum
        
        newbigrho = oldbigrho - stepval*stepdirec
           
        !what is the new value for the gradient, i.e. evaluate gradient at newbigrho
        do zct=1,znum
        do momct=1,momuse
            gradct = (zct-1)*momuse + momct
            rhomat(zct,momct) = newbigrho(gradct)
        end do !momct
        end do !zct
        newgrad = gradPint()
        newfunc = Pint()
        
        !check for nonexistence of new function or gradient value
        nannewfunc = isnan(newfunc)
        nannewgrad = .FALSE.
        do ct1=1,znum*momuse
            if (isnan(newgrad(ct1))) nannewgrad = .TRUE.
        end do 
        
        !check wolfe sufficient decrease condition
        wolfe1 = ( newfunc <= (oldfunc + c1 * stepval * sum(stepdirec*oldgrad) ) )
        
        !check wolfe slope condition
        wolfe2 = ( sum(stepdirec * newgrad) >= c2 * sum(stepdirec * oldgrad ) )
        
        !if existence and wolfe conditions hold, then done 
        if ((nannewfunc.eqv..FALSE.).and.(nannewgrad.eqv..FALSE.).and.&
            (wolfe1.eqv..TRUE.).and.(wolfe2.eqv..TRUE.)) exit
        
        !if made it past, then you failed, so update
        if (stepct<(maxstepnum-1)) then
            stepval = stepval * stepfrac
        else if (stepct==(maxstepnum-1)) then
            stepval = stepfallback
        end if
        

    end do !stepct
    
    rhoerror = maxval(abs(newbigrho - oldbigrho))/maxval(abs(oldbigrho))
    graderror = maxval(abs(newgrad))
    funcerror = abs(newfunc-oldfunc)/abs(oldfunc)

    
    !call sleep(10)
    if (iter>10) then
    if (graderror<broydengradtol) exit
    end if
    
    diffgrad = newgrad - oldgrad
    diffrho = newbigrho - oldbigrho
    
end do !iter

write(*,"(A,I5,A,F7.4,A,F7.4,A,I5,A,F7.4)") "it = ",iter-1,", stepval = ",stepval,", grad = ",&
            graderror,", stepct = ",stepct-1,", val = ",Pint()

!check for issues and return SS coefficients if there are optimization failures
if ((nannewfunc.eqv..TRUE.).or.(nannewgrad.eqv..TRUE.)) then
    write(*,*) "ISSUE WITH OPTIMIZATION, RETURNING SS VALUES"
    rhomat = rhoSS
    intvec = intSS
end if


end subroutine findrhoSR1wolfe

end program kt_param
