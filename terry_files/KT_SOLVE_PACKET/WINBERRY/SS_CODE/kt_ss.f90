!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kt_ss.f90
!
! Fortran code for the KS solution of the Khan and Thomas (2008) 
! model in the steady-state for use in the WINBERRY solution.
!
! 'Alternative Methods for Solving Heterogeneous Firm Models'
! Stephen Terry (2017)
!
! This Version : 01/16/17
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module modparams
implicit none

!This module contains fixed model parameter values, available
!to the rest of the code. It also declares constants which are made available
!to subprograms when set elsewhere. It also declares some variables
!threadprivate, which makes them unique to OpenMP parallel threads.

!real parameters
double precision, parameter :: alpha = 0.256 !capital elasticity
double precision, parameter :: nu = 0.640 !labor elasticity
double precision, parameter :: phi = 2.4 !labor disutility
double precision, parameter :: xibar = 0.0083 !upper bound of capital AC distribution
double precision, parameter :: beta = 0.977 !discount rate
double precision, parameter :: delta = 0.069 !capital depreciation rate
double precision, parameter :: aval = 1.0 !fixed value of aggregate productivity
double precision, parameter :: kmin= 0.1 !min value on idio capital grid
double precision, parameter :: kmax = 8.0 !max value on idio capital grid
double precision, parameter :: rhoz = 0.859 !persistence of idio prod
double precision, parameter :: sigmaz = 0.022 !stdev of shock to log idio prod
double precision, parameter :: nstdevz = 2.0 !multiple of st dev of idio prod process to discretize
double precision, parameter :: vferrortol=1e-4 !tolerance for vf error 
double precision, parameter :: kprimeerrortol = 1e-6 !tolerance for error on kprime
double precision, parameter :: xistarerrortol = 1e-6 !tolerance for error on xistar
double precision, parameter :: plb = 2.15 !lower bound for bisection on price
double precision, parameter :: pub = 2.3 !upper bound for bisection on price
double precision, parameter :: perrortol = 1e-5 !tolerance for bisection on price
double precision, parameter :: poltol = 1e-6 !tolerance for optimization of idio capital policies
double precision, parameter :: disterrortol = 1e-6!tolerance on distribution fixed point
double precision, parameter :: broydenrhotol = 1e-6!tolerance on change in rho values in Broyden optimization
double precision, parameter :: broydengradtol = 1e-2 !tolerance on size of gradient in Broyden optimization
double precision, parameter :: broydenfunctol = 1e-7 !tolerance on size of function eval diff
double precision, parameter :: stepsize = 0.001 !Broyden step size

!integer parameters
integer, parameter :: knum = 10 !number of idio capital nodes for spline solution
integer, parameter :: znum = 5 !number of idio prod states in discretization
integer, parameter :: kdensenum = 50 !number of idio capital points in distribution
integer, parameter :: maxpolit = 1000 !max number of policy iterations
integer, parameter :: maxaccelit = 50 !max number of Howard accelerations per VF iteration
integer, parameter :: maxpit = 100 !max number of iterations on price bisection
integer, parameter :: maxdistit = 5000 !max number of iterations for distribution fixed point
integer, parameter :: numconstants = 16 !number of constants to pass to output
integer, parameter :: momnum = 25 !number of moments to make available for parametrized distributions
integer, parameter :: momuse = 4 !number of moments to use for parametrized distributions
integer, parameter :: nsimp = 100 !number of Simpson quadrature nodes (needs to be even)
integer, parameter :: maxbroydit = 100000 !max number of Broyden iters


!parameters for SR1 quasi-newton with step size line search
integer, parameter :: maxstepnum = 100
double precision, parameter :: stepfrac = 0.9
double precision, parameter :: stepinit = 0.09
double precision, parameter :: c1 = 1e-3
double precision, parameter :: c2 = 0.9
double precision, parameter :: stepfallback = 0.00001


!insert the pso parameters for rho parameter minimization
double precision, parameter :: bound = 50.0
integer, parameter :: nvar = znum*momuse
integer, parameter :: npart = 500
double precision, parameter :: xtol = 1e-3
double precision, parameter :: ftol = 1e-3
double precision, parameter :: xquicktol = 1e-3
double precision, parameter :: phipso(2) = (/2.05,2.05/)
integer, parameter :: xquicknum = 5
integer, parameter :: maxpsoit = 1000

!real values available elsewhere
double precision :: pval,wval

!real allocatable arrays available elsewhere
double precision, allocatable :: Vold(:,:),V2old(:,:),pr_mat_z(:,:),k0(:),kdense0(:)

!idio prod variable which will be threadprivate below
integer :: RHSzct
!$omp threadprivate(RHSzct)

!moment stuff to be available to subfunctions
double precision :: rhomat(znum,momuse),intvec(znum),mommat(znum,momuse),momstoremat(znum,momnum)
double precision :: simpnodes(nsimp+1),simpweights(nsimp+1),rhoSS(znum,momuse),intSS(znum)


end module modparams

program kt_ss
use base_lib
use omp_lib
use modparams
implicit none

integer :: zct,kct,zprimect,polct,accelit,ct,distct,piter,ind,momct
    
integer, allocatable :: kprimedenseind(:,:),knoadjdenseind(:,:)
    
double precision :: start,finish,kprimeval,Vnaval,Vaval,Vval,xival,Vnextval,vferror,kprimeerror,&
    xistarerror,disterror,wgt,Yvalp,Ivalp,Cvalp,yval,ival,perror,zval,kval,padjust,pvala,pvalb,Cerror,&
    Kvalp,Nvalp,rhovec(znum*momuse),lb(nvar),ub(nvar),Pvalpso

double precision, allocatable :: V(:,:),kprime(:,:),kprimeold(:,:),xistar(:,:),xistarold(:,:),&
    Va(:,:),Vna(:,:),z0(:),dist(:,:),distold(:,:),kprimedense(:,:),xistardense(:,:),constantvec(:),&
    kprimedensewgt(:,:),knoadjdense(:,:),knoadjdensewgt(:,:)

!!!!!!PROGRAM PRELIMINARIES

!open log file
start = omp_get_wtime()
open(13,file="kt_ss.txt")

!parallelization check

!$omp parallel
write(*,*) "Parallel hello to you!"
write(13,*) "Parallel hello to you!"
!$omp end parallel

write(*,*) " "
write(13,*) " "

!do allocate memory based on dimensions set in parameter module
allocate(k0(knum),z0(znum),pr_mat_z(znum,znum),V(znum,knum),Vold(znum,knum),V2old(znum,knum),kprime(znum,knum),&
    kprimeold(znum,knum),xistar(znum,knum),xistarold(znum,knum),Va(znum,knum),Vna(znum,knum),dist(znum,kdensenum),&
    distold(znum,kdensenum),kdense0(kdensenum),kprimedense(znum,kdensenum),xistardense(znum,kdensenum),&
    kprimedensewgt(znum,kdensenum),kprimedenseind(znum,kdensenum),knoadjdense(znum,kdensenum),&
    knoadjdensewgt(znum,kdensenum),knoadjdenseind(znum,kdensenum))

!!!!!!!SET UP GRIDS AND DISCRETIZE IDIO PROD

!pso optimization bounds
lb(:) = -1.0*bound
ub(:) = bound

!set up idio capital grids
call linspace(k0,log(kmin),log(kmax),knum);
k0=exp(k0);
call linspace(kdense0,log(kmin),log(kmax),kdensenum);
kdense0=exp(kdense0);

!do discretization of the z process to set up z0 grid and transitions
call tauchen(znum,rhoz,sigmaz,nstdevz,pr_mat_z,z0)

!call Simpson quadrature weights and nodes for moment integration
call qsimpweightsnodes(kmin,kmax,nsimp,simpweights,simpnodes)

!!!!!!!SOLVE THE MODEL

!price bisection loop
pvala = plb
pvalb = pub

do piter = 1,maxpit

!bisect current price interval
pval = 0.5 * (pvala + pvalb)

!given price, what is implied wage?
wval = phi/pval

!now, initialize the value and policy functions
xistarold(:,:) = 0.5
do zct=1,znum
    kprimeold(zct,:) = k0    
end do !zct

!initialize value function and second derivatives to 0
Vold(:,:) = 0.0
V2old(:,:) = 0.0

!!!!!!!POLICY ITERATION LOOP
write(*,*) "Doing policy iterations."
write(13,*) "Doing policy iterations."

do polct=1,maxpolit
    
    !do Howard acceleration to create RHS of Bellman equation given kprimeold
    do accelit=1,maxaccelit
        
        !loop over the state space
        do zct=1,znum
        do kct=1,knum
            
            !extract states
            zval = z0(zct)
            kval = k0(kct)
            
            !determine value when not adjusting
            Vnaval = pval * freduced(zval,kval) !current period return
            kprimeval = max((1.0-delta)*kval,kmin) !transition without adjustment
            
            !add continuation value when not adjusting
            do zprimect=1,znum
                call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprimeval,Vnextval)
                Vnaval = Vnaval + beta * pr_mat_z(zct,zprimect) * Vnextval
            end do !zprimect
            
            !determine value when adjusting
            Vaval = pval * freduced(zval,kval) !current period return
            kprimeval = kprimeold(zct,kct) !what is the policy, based on last pol array?
            Vaval = Vaval - pval * ( kprimeval - ( 1.0 - delta ) * kval )
            
            !add continuation value when adjusting
            do zprimect=1,znum
                call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprimeval,Vnextval)
                Vaval = Vaval + beta * pr_mat_z(zct,zprimect) * Vnextval
            end do !zprimect
            
            !now that Va and Vna are determined, what is V?
            xival = (Vaval - Vnaval)/phi !cutoff AC draw value for adjustment
            Vval = - 1.0 * phi * expecxi(xival) + cdfxi(xival) * Vaval + &
                (1.0 - cdfxi(xival)) * Vnaval
            V(zct,kct) = Vval
            
            
        end do !kct
        end do !zct
        
        !now, convert Vold to V, and spline the value function, 
        !for increment in Howard acceleration
        Vold = V
        V(:,:) = 0.0
        
        do zct=1,znum
            call spline(k0,Vold(zct,:),knum,dble(1.0e30),dble(1.0e30),V2old(zct,:))
        end do !zct
        
    end do !accelit
    
    !now that Howard acceleration is done, update policies
    
    !loop over states
    do zct=1,znum
        
        !find optimal capital policy when adjusting, which only depends on zct
        RHSzct = zct !this is made available to the fnRHS function
        Vaval = brent(k0(1),k0(knum/2),k0(knum),fnRHS,poltol,kprimeval) !this returns kprimeval as optimal capital policy
        kprime(zct,:) = kprimeval !store optimal policy in kprime array
                
    do kct=1,knum
        
        !extract state values
        zval = z0(zct)
        kval = k0(kct)
        
        !what is value when not adjusting?
        Vnaval = pval * freduced(zval,kval) !current period return
        kprimeval = max((1.0-delta)*kval,kmin) !no adjustment in capital
        
        !add continuation value
        do zprimect=1,znum
            call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprimeval,Vnextval)
            Vnaval = Vnaval + beta * pr_mat_z(zct,zprimect) * Vnextval
        end do !zprimect
        Vna(zct,kct) = Vnaval !store in Vna array
        
        !what is value when adjusting?
        Vaval = pval*freduced(zval,kval)!current period return
        kprimeval = kprime(zct,kct) !extract policy determined above        
        Vaval = Vaval - pval * ( kprimeval - ( 1.0 - delta ) * kval )
        
        !add continuation value
        do zprimect=1,znum
            call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprimeval,Vnextval)
            Vaval = Vaval + beta * pr_mat_z(zct,zprimect) * Vnextval
        end do !zprimect
        
        Va(zct,kct) = Vaval !store in Va array
        
        !now that the value when adjusting and when not adjusting are determined, create and store total value
        xival = (Vaval-Vnaval)/phi  !cutoff AC value
        xistar(zct,kct) = xival !store in xistar array
        Vval = -1.0*phi*expecxi(xival) + cdfxi(xival) * Vaval &
            + (1.0 - cdfxi(xival))*Vnaval
            
        V(zct,kct) = Vval !store total value
        
    end do !kct
    end do !zct
    
    !now, determine VF error and exit if needed
    vferror = maxval(abs(V-Vold))
    kprimeerror = maxval(abs(kprime-kprimeold))
    xistarerror = maxval(abs(xistar-xistarold))
    
    if (mod(polct,1)==0.and.mod(piter,1)==0) then
    write(13,*) "Pol iter = ",polct,"VF error = ",vferror
    write(13,*) "Pol iter = ",polct,"Kprime error = ",kprimeerror
    write(13,*) "Pol iter = ",polct,"Xistar error = ",xistarerror
    write(13,*) " "
    
    write(*,*) "Pol iter = ",polct,"VF error = ",vferror
    write(*,*) "Pol iter = ",polct,"Kprime error = ",kprimeerror
    write(*,*) "Pol iter = ",polct,"Xistar error = ",xistarerror
    write(*,*) " "
    end if
    
    !if policies and cutoffs or the VF itself have converged, stop
    if (kprimeerror<kprimeerrortol .and. xistarerror<xistarerrortol) exit
    
    !if not converged, update
    kprimeold = kprime
    xistarold = xistar
    
end do !polct
write(*,*) "Done with policy iterations."
write(13,*) "Done with policy iterations."

write(*,*) " "
write(13,*) " "

!!!!!COMPUTE ERGODIC DISTRIBUTION AND IMPLICED PRICES

!what are policies and AC cutoffs at kdense0 rather than coarse k0?
kprimedense(:,:) = 0.0
xistardense(:,:) = 0.0

do zct=1,znum
    
    !capital policies pinned down by z, xistar comes from spline evaluation
    kprimedense(zct,:) = kprime(zct,1)
    
    do kct=1,kdensenum
    
        !extract state values
        zval = z0(zct)
        kval = kdense0(kct)
        
        !what is value when not adjusting?
        Vnaval = pval * freduced(zval,kval) !current period return
        kprimeval = max((1.0-delta)*kval,kmin) !no adjustment in capital
        
        !add continuation value
        do zprimect=1,znum
            call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprimeval,Vnextval)
            Vnaval = Vnaval + beta * pr_mat_z(zct,zprimect) * Vnextval
        end do !zprimect
        
        !what is value when adjusting?
        Vaval = pval*freduced(zval,kval)!current period return
        kprimeval = kprime(zct,1) !extract policy determined above, picking first element b/c only depends on z
        Vaval = Vaval - pval * ( kprimeval - ( 1.0 - delta ) * kval )
        
        !add continuation value
        do zprimect=1,znum
            call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprimeval,Vnextval)
            Vaval = Vaval + beta * pr_mat_z(zct,zprimect) * Vnextval
        end do !zprimect
        
        !now that the value when adjusting and when not adjusting are determined, create and store total value
        xival = (Vaval-Vnaval)/phi  !cutoff AC value
        
        xistardense(zct,kct) = xival
        
    end do !kct
end do !zct

!now that the dense continuous policy is available, convert to grid indexes/weights
kprimedensewgt(:,:) = 0.0
kprimedenseind(:,:) = 0
call convertpolicy(kprimedense,kprimedensewgt,kprimedenseind)

!what are the noadjustment next period values, weights, and indexes?
knoadjdense(:,:) = 0.0
do zct=1,znum
do kct=1,kdensenum
    knoadjdense(zct,kct) = max((1.0-delta)*kdense0(kct),kmin)
end do !kct    
end do !zct

knoadjdensewgt(:,:)=0.0
knoadjdenseind(:,:) = 0
call convertpolicy(knoadjdense,knoadjdensewgt,knoadjdenseind)

!now, iterate to find ergodic distribution and find implied price

!initialize distributions
distold(:,:) = 0.0
distold(znum/2,kdensenum/2) = 1.0
distold = distold / sum(distold)

dist(:,:) = 0.0

!loop over distribution pushforwards
write(*,*) "Now doing distribution iterations."
write(13,*) "Now doing distribution iterations."
do distct = 1,maxdistit
    
    !loop over states
    do zct=1,znum
    do kct=1,kdensenum
        
        !what is the AC cutoff from above, and the associated prob of adjustment
        xival = xistardense(zct,kct)
        padjust = cdfxi(xival)
        
        !if you do NOT adjust, what are index/weight?
        ind = knoadjdenseind(zct,kct)
        wgt = knoadjdensewgt(zct,kct)
        
        !iterate non-adjusting mass forward
        do zprimect=1,znum
            
            dist(zprimect,ind) = dist(zprimect,ind) + distold(zct,kct) * pr_mat_z(zct,zprimect) * &
                (1.0 - padjust) * (1.0 - wgt)
            
            dist(zprimect,ind+1) = dist(zprimect,ind+1) + distold(zct,kct) * pr_mat_z(zct,zprimect) * &
                (1.0 - padjust) * wgt

        end do !zprimect
        
        !if you DO adjust, what are index/weight?
        ind = kprimedenseind(zct,kct)
        wgt = kprimedensewgt(zct,kct)
        
        !iterate adjusting mass forward
        do zprimect=1,znum
            
            dist(zprimect,ind) = dist(zprimect,ind) + distold(zct,kct) * pr_mat_z(zct,zprimect) * &
                padjust * (1.0 - wgt)
            
            dist(zprimect,ind+1) = dist(zprimect,ind+1) + distold(zct,kct) * pr_mat_z(zct,zprimect) * &
                padjust * wgt

        end do !zprimect
        
    end do !kct
    end do !zct
    
    disterror = maxval(abs(dist - distold))
    dist = dist/sum(dist)
    
    if (mod(distct,10)==0.and.mod(piter,1)==0) then
        write(13,*) "Dist iter = ",distct,"Dist error = ",disterror
        write(*,*) "Dist iter = ",distct,"Dist error = ",disterror
    end if
    
    !has the distribution converged?
    if (disterror<disterrortol) exit
    
    !if not, iterate forward
    distold = dist
    dist(:,:) = 0.0
    
end do !distct
write(*,*) "Finished with distribution iterations."
write(13,*) "Finished with distribution iterations."

write(*,*) " "
write(13,*) " "

!now, compute implied prices from the distribution above
Yvalp = 0.0
Ivalp = 0.0
Cvalp = 0.0
Kvalp = 0.0
Nvalp = 0.0

!$omp parallel private(zct,kct,zval,kval,kprimeval,xival,padjust,&
!$omp& yval,ival) reduction(+:Yvalp,Ivalp,Kvalp,Nvalp)
!$omp do collapse(2)
do zct=1,znum
do kct=1,kdensenum
    
    !extract states
    zval = z0(zct)
    kval = kdense0(kct)
    
    !policy and cutoff?
    kprimeval = kprimedense(zct,kct)
    xival = xistardense(zct,kct)
    padjust = cdfxi(xival)
    
    !output and investment when adjusting?
    yval = yreduced(zval,kval)
    ival = kprimeval - (1.0 - delta)*kval
    
    !increment aggregates
    Yvalp = Yvalp + dist(zct,kct) * yval
    Ivalp = Ivalp + dist(zct,kct) * padjust * ival
    Kvalp = Kvalp + dist(zct,kct) * kval
    Nvalp = Nvalp + dist(zct,kct) * (nreduced(zval,kval) + expecxi(xival))
    
end do !kct
end do !zct
!$omp end do nowait
!$omp end parallel

!implied consumption and clearing error?
Cvalp = Yvalp - Ivalp
Cerror  = 1/pval - Cvalp
perror = pvalb - pvala

!report aggregates
if (mod(piter,1)==0) then
write(*,"(A,I5,A,F7.5,A,F7.5)") "P Iter = ",piter,", p = ",pval,", perror = ",perror
write(*,"(A,F7.5,A,F7.5,A,F7.5)") "C(p) = ",Cvalp,", 1/p = ",1/pval,", diff = ",Cerror
write(*,"(A,F7.5,A,F7.5,A,F7.5,A,F7.5)") "Y(p) = ",Yvalp,", I(p) = ",Ivalp,", K(p) = ",Kvalp,", N(p) = ",Nvalp
write(*,*) "#######################################################"
write(*,*) " "

write(13,"(A,I5,A,F7.5,A,F7.5)") "P Iter = ",piter,", p = ",pval,", perror = ",perror
write(13,"(A,F7.5,A,F7.5,A,F7.5)") "C(p) = ",Cvalp,", 1/p = ",1/pval,", diff = ",Cerror
write(13,"(A,F7.5,A,F7.5,A,F7.5,A,F7.5)") "Y(p) = ",Yvalp,", I(p) = ",Ivalp,", K(p) = ",Kvalp,", N(p) = ",Nvalp
write(13,*) "#######################################################"
write(13,*) " "

end if
!if price interval is small enough, exit clearing bisection
if (perror<perrortol) exit

!if prices not yet converged, update bisection endpoints
if (Cerror<0.0) then
    pvalb = pval
else if (Cerror>=0.0) then
    pvala = pval
end if


end do !piter 
!this ends loop over price bisection


!now that the SS equilibrium has been computed, determine the moments of the SS distribution, and store them

!rhomat(znum,momuse) (zct,momct) = rho^zct_momct
!intvec(znum) (zct) = 1 / rho^zct_0
!mommat(znum,momuse) (zct,momct) = m^zct_momct, momct=1,# of moments used in solution
!momstoremat(znum,momnum) (zct,momct) = m^zct_momct, momct=1,# of moments to store in case needed for agg solution (momnum>=momuse)

rhomat(:,:) = 0.0
intvec(:) = 0.0
mommat(:,:) = 0.0
momstoremat(:,:) = 0.0

!compute means of SS dist, conditional on zct, noting dist(zct=1,znum,kct=1,kdensenum)
do zct=1,znum
    momstoremat(zct,1) = sum(dist(zct,:)*kdense0)/sum(dist(zct,:))
end do !zct

!now compute the other centered moments
do momct=2,momnum
    do zct=1,znum
        momstoremat(zct,momct) = sum( dist(zct,:) * ( ( kdense0 - momstoremat(zct,1) )**dble(momct) ) )/sum(dist(zct,:))
    end do !zct
end do !momct

mommat = momstoremat(:,1:momuse)




!!!!!!OUTPUT DATA TO TEXT FILES

!constants to pass to MATLAB
constantvec = (/ pval,wval,Yvalp,Ivalp,Nvalp,dble(znum),dble(knum),dble(kdensenum),alpha,nu,&
    phi,xibar,beta,delta,aval,Kvalp/)
open(8,file="constantvec.txt")
do ct=1,numconstants
write(8,*) constantvec(ct)
end do !ct
close(8)

!idio capital solution nodes
open(8,file="k0.txt")
do kct=1,knum
write(8,*) k0(kct)
end do !kct
close(8)

!idio capital distribution grid
open(8,file="kdense0.txt")
do kct=1,kdensenum
write(8,*) kdense0(kct)
end do !kct
close(8)

!idio prod grid
open(8,file="z0.txt")
do ct=1,znum
write(8,*) z0(ct)    
end do !ct
close(8)

!transition matrix for idio prod
open(8,file="pr_mat_z.txt")
do zct=1,znum
write(8,*) pr_mat_z(zct,:)
end do !zct
close(8)

!idio capital policy
open(8,file="kprimeSS.txt")
do zct=1,znum
do kct=1,knum
write(8,*) kprime(zct,kct)    
end do !kct
end do !zct
close(8)

!adjustment thresholds
open(8,file="xistarSS.txt")
do zct=1,znum
do kct=1,knum
write(8,*) xistar(zct,kct)    
end do !kct
end do !zct
close(8)

!firm value
open(8,file="VSS.txt")
do zct=1,znum
do kct=1,knum
write(8,*) V(zct,kct)    
end do !kct
end do !zct
close(8)

!firm value
open(8,file="V2old.txt")
!do zct=1,znum
do kct=1,knum
write(8,*) V2old(:,kct)    
end do !kct
!end do !zct
close(8)

!firm value when adjusting
open(8,file="VaSS.txt")
do zct=1,znum
do kct=1,knum
write(8,*) Va(zct,kct)    
end do !kct
end do !zct
close(8)

!firm value when not adjusting
open(8,file="VnaSS.txt")
do zct=1,znum
do kct=1,knum
write(8,*) Vna(zct,kct)    
end do !kct
end do !zct
close(8)

!capital policy on dense grid
open(8,file="kprimedense.txt")
do zct=1,znum
do kct=1,kdensenum
write(8,*) kprimedense(zct,kct)    
end do !kct
end do !zct
close(8)

!xistar on densegrid
open(8,file="xistardense.txt")
do zct=1,znum
do kct=1,kdensenum
write(8,*) xistardense(zct,kct)    
end do !kct
end do !zct
close(8)

!ergodic distribution
open(8,file="dist.txt")
do zct=1,znum
do kct=1,kdensenum
write(8,*) dist(zct,kct)    
end do !kct
end do !zct
close(8)


write(*,*) "started at ",omp_get_wtime()-start," seconds"
!call findrhobroyden()
!call findrhoSR1()
call findrhoSR1wolfe(0)

write(*,*) "Pint = ",Pint()
write(*,*) "gradPint = ",gradPint()
open(8,file="Densities.txt")
do zct=1,znum
    !write(*,*) "density for zct = ",zct
do kct=1,kdensenum
    kval = kdense0(kct)
write(8,*) Fk(kval,zct)
end do !kct
end do !zct
close(8)
do zct=1,znum
    do momct=1,momuse
        write(*,*) rhomat(zct,momct)
    end do
end do

write(*,*) intvec

!write some aggregates
open(8,file="aggsSS.txt")
write(8,*) pval
write(8,*) Kvalp
write(8,*) Yvalp
write(8,*) Ivalp
write(8,*) Nvalp
write(8,*) aval
close(8)

!write moments
open(8,file="momSS.txt")
do zct=1,znum
do momct=1,momnum
    write(8,*) momstoremat(zct,momct)
end do !momct
end do !zct
close(8)

!write simp weights and nodes
open(8,file="simpnodes.txt")
do ct=1,nsimp+1
    write(8,*) simpnodes(ct)
end do !zct
close(8)

!write simp weights and nodes
open(8,file="simpweights.txt")
do ct=1,nsimp+1
    write(8,*) simpweights(ct)
end do !zct
close(8)

!write coefficients
open(8,file="rhoSS.txt")
do zct=1,znum
do momct=1,momuse
    write(8,*) rhomat(zct,momct)
end do !momct
end do !zct
close(8)

!write integral array
open(8,file="intSS.txt")
do zct=1,znum
    write(8,*) intvec(zct)
end do !zct
close(8)


finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished in ",finish-start," seconds."
close(13); !closing log file

contains

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

!this function returns -1 * ( -p * k' + beta * E V(z',k') ), for minimization wrt k'

double precision :: kprimeval,Vnextval

integer :: zprimect

!note those things available to the function:
!k0, Vold, V2old, knum, pr_mat_z (all threads)
!RHSzct (one thread only via "threadprivate" status)

!current period return (the only terms relevant to k optimization)
fnRHS = -pval * kprimeval

!add the continuation value
do zprimect=1,znum
    call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprimeval,Vnextval)
    fnRHS = fnRHS + beta * pr_mat_z(RHSzct,zprimect) *  Vnextval   
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

!make initial guesses for inv Jacobian approx - identity matrix
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
    
    !what is new guess for inv Jacobian?
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
    funcerror = abs(newfunc-oldfunc)/abs(oldfunc)

    !if ((isnan(rhoerror).eqv..FALSE.).and.(isnan(graderror).eqv..FALSE.)) then
    if (mod(iter,100)==1) then
        write(*,"(A,I5,A,F7.4,A,F7.4,A,F7.4,A,F7.4)") "it = ",iter,", rho = ",rhoerror,", grad = ",&
            graderror,", fun = ",funcerror,", val = ",Pint()
    end if
     
    !call sleep(10)
    if (iter>10) then
    if (graderror<broydengradtol.or.funcerror<broydenfunctol) exit
    if (rhoerror<broydenrhotol) then
        write(*,*) "Error finding minimizer. Exiting Brodyen."    
        exit
    end if
    end if
    
    diffgrad = newgrad - oldgrad
    diffrho = newbigrho - oldbigrho
    
end do !iter


end subroutine findrhoBroyden


subroutine findrhoSR1()
implicit none

double precision :: newJinv(znum*momuse,znum*momuse),oldJinv(znum*momuse,znum*momuse),oldbigrho(znum*momuse),&
    newbigrho(znum*momuse),oldgrad(znum*momuse),newgrad(znum*momuse),diffrho(znum*momuse),diffgrad(znum*momuse),&
    numval(znum*momuse),denomval,nummat(znum*momuse,znum*momuse)
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

!make initial guesses for inv Jacobian approx - identity matrix
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
    
    !what is new guess for inv Jacobian?
    newJinv = oldJinv + (1.0 / denomval) * nummat
    
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
    funcerror = abs(newfunc-oldfunc)/abs(oldfunc)

    if (mod(iter,2500)==1) then
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

end subroutine findrhoSR1



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

end program kt_ss