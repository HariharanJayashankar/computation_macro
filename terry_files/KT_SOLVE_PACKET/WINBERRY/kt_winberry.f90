!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kt_winberry.f90
!
! Fortran code to solve the Khan and Thomas (2008)
! model using the Winberry (2016) method.
!
! 'Alternative Methods for Solving Heterogeneous Firm Models'
! Stephen Terry (2017)
!
! This Version : 01/13/17
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
double precision, parameter :: kmin= 0.1 !min value on idio capital grid
double precision, parameter :: kmax = 8.0 !max value on idio capital grid
double precision, parameter :: nstdevz = 2.0 !number of st dev's to span for idio prod discretization
double precision, parameter :: rhoz = 0.859  !persistence of idio prod shock
double precision, parameter :: sigmaz = 0.022 !st dev of shock to idio prod
double precision, parameter :: nstdeva = 2.0 !number of st dev's to span for agg prod discretization
double precision, parameter :: rhoa = 0.859  !persistence of agg prod shock
double precision, parameter :: sigmaa=0.014 !st dev of shock to agg prod
double precision, parameter :: shocksizeIRF=sigmaa !shock size in IRF simulations
double precision, parameter :: vferrortol=1e-4 !tolerance for vf error 
double precision, parameter :: kprimeerrortol = 1e-6 !tolerance for error on kprime
double precision, parameter :: xistarerrortol = 1e-6 !tolerance for error on xistar
double precision, parameter :: plb = 2.0 !lower bound for bisection on price
double precision, parameter :: pub = 2.4 !upper bound for bisection on price
double precision, parameter :: perrortol = 1e-5 !tolerance for bisection on price
double precision, parameter :: poltol = 1e-6 !tolerance for optimization of idio capital policies
double precision, parameter :: disterrortol = 1e-5!tolerance on distribution fixed point
double precision, parameter :: broydenrhotol = 1e-6!tolerance on change in rho values in Broyden optimization
double precision, parameter :: broydengradtol = 1e-2 !tolerance on size of gradient in Broyden optimization
double precision, parameter :: broydenfunctol = 1e-7 !tolerance on size of function eval diff
double precision, parameter :: stepsize = 0.001 !Broyden step size
double precision, parameter :: diffstep  = 1.0e-6 !step size on numerical differentiation

!integer parameters
integer, parameter :: knum = 10 !number of idio capital nodes for spline solution
integer, parameter :: znum = 5 !number of idio prod states in discretization
integer, parameter :: anum = 5 !number of agg prod states in discretization
integer, parameter :: kdensenum = 50 !number of idio capital points in distribution
integer, parameter :: maxpolit = 1000 !max number of policy iterations
integer, parameter :: maxaccelit = 50 !max number of Howard accelerations per VF iteration
integer, parameter :: maxpit = 100 !max number of iterations on price bisection
integer, parameter :: maxergit = 5000 !max number of iters to compute ergodic dist of agg prod
integer, parameter :: momnum = 25 !number of moments to make available for parametrized distributions
integer, parameter :: momuse = 4 !number of moments to use for parametrized distributions
integer, parameter :: nsimp = 100 !number of Simpson quadrature nodes (needs to be even)
integer, parameter :: maxbroydit = 100000 !max number of Broyden iters
integer, parameter :: numper = 2500 !number of periods in unconditional simulation
integer, parameter :: numdiscard = 500 !number of periods discarded from unconditional simulation 
integer, parameter :: seedint=2503 !random number seed
integer, parameter :: ainit=3 !initial gridpoint for agg prod in simulations
integer, parameter :: numsimIRF=2000 !number of IRF economies to simulate
integer, parameter :: numperIRF=50 !number of periods per IRF economy
integer, parameter :: shockperIRF=25 !period at which to shock each IRF economy
integer, parameter :: shockanumIRF = anum !grid point that you impose in IRF period
integer, parameter :: numX = 2*znum*knum+znum*momuse+5+znum !number of elements in the X vector
integer, parameter :: numeta = 2*znum*knum+znum !number of expectational errors in the eta vector
integer, parameter :: numeps = 1 !number of aggregate shocks
integer, parameter :: nummicro = 7 !number of micro moments to compute

integer :: doMICROsim !if 1, then compute micro moments in Fsys, if 0, then don't
integer, parameter :: doMICROsimstore = 0 !if 1, then compute micro moments in Fsys, if 0, then don't
integer, parameter :: doIRF = 1 !if 1, then do IRF simulation

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

!real arrays available elsewhere
double precision :: momstoreSS(znum,momnum),momSS(znum,momuse),rhoSS(znum,momuse),intSS(znum),&
	kprimeSS(znum,knum),xistarSS(znum,knum),Vss(znum,knum),VaSS(znum,knum),VnaSS(znum,knum),&
	pSS,KSS,YSS,ISS,NSS,aSS,k0(knum),z0(znum),pr_mat_z(znum,znum),pmin1,p,Ymin1,Y,Imin1,I,Nmin1,N,&
	amin1,a,V(znum,knum),V2(znum,knum),wmin1,w,ergdistz(znum),F1(numX,numX),F2(numX,numX),&
	F3(numX,numeta),F4(numX,numeps),a0(anum),asimshock(numper),asimshockIRF(numperIRF,numsimIRF),&
	arriveshockIRF(numsimIRF),ergdista(anum),ergdistaold(anum),pr_mat_a(anum,anum),kdense0(kdensenum)

!integer arrays available elsewhere
integer :: asimpos(numper),asimposIRF(numperIRF,numsimIRF,2)
integer, allocatable :: seedarray(:)
integer :: seeddim

!idio prod variable which will be threadprivate below
integer :: RHSzct
!$omp threadprivate(RHSzct)

!moment stuff to be available to subfunctions
double precision :: rhomat(znum,momuse),intvec(znum),mommat(znum,momuse),momstoremat(znum,momnum)
double precision :: simpnodes(nsimp+1),simpweights(nsimp+1)


end module modparams

program kt_winberry
use base_lib
use omp_lib
use modparams
implicit none

double precision :: start,finish,Xvec(numX),Xmin1vec(numX),etavec(numeta),epsvec(numeps),&
	XSS(numX),etaSS(numeta),epsSS(numeps),FSS(numX),ergerrz,ergoldz(znum),ergnewz(znum),h,&
	Fvec(numX),MICROSS(nummicro),MICROvec(nummicro),Xsim(numX,numper),MICROsimstore(nummicro,numper),&
	epssim(numper)
integer :: ct,zct,momct,kct,zprimect,ct2,t

!!!!!!PROGRAM PRELIMINARIES

!open log file
start = omp_get_wtime()
open(13,file="kt_winberry.txt")

!parallelization check

!$omp parallel
write(*,*) "Parallel hello to you!"
write(13,*) "Parallel hello to you!"
!$omp end parallel

write(*,*) " "
write(13,*) " "

open(8,file="constants.txt")
write(8,*) numX
write(8,*) numeta
write(8,*) numeps
write(8,*) znum
write(8,*) knum
write(8,*) anum
write(8,*) numper
write(8,*) numperIRF
write(8,*) numsimIRF
write(8,*) shockperIRF
write(8,*) sigmaa
write(8,*) numdiscard
write(8,*) momuse
close(8)


!!!!!!!READ IN DATA FROM SS SOLUTION AND DO SOME BASIC SETUP

!discretize exogenous processes A, z, and simulate A
call discretize_simulate()

!set up idio capital grids
call linspace(k0,log(kmin),log(kmax),knum);
k0=exp(k0);

call linspace(kdense0,log(kmin),log(kmax),kdensenum);
kdense0=exp(kdense0);


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

!read in SS capital policy, for initialization
open(8,file="kprimeSS.txt")
open(9,file="xistarSS.txt")
open(10,file="VSS.txt")
open(11,file="VaSS.txt")
open(12,file="VnaSS.txt")
do zct=1,znum
do kct=1,knum
    read(8,*) kprimeSS(zct,kct)
    read(9,*) xistarSS(zct,kct)
    read(10,*) VSS(zct,kct)
    read(11,*) VaSS(zct,kct)
    read(12,*) VnaSS(zct,kct)
end do !kct
end do !zct
close(8)
close(9)
close(10)
close(11)
close(12)

!read some aggregates
open(8,file="aggsSS.txt")
read(8,*) pSS
read(8,*) KSS
read(8,*) YSS
read(8,*) ISS
read(8,*) NSS
read(8,*) aSS
close(8)

!!!!!!!!FIND ERGODIC DISTRIBUTION OF IDIO PROD
ergdistz(:) = 0.0
ergdistz(znum/2) = 1.0
ergoldz = ergdistz

do ct=1,maxergit
    
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

!!!!!!CONSTRUCT XSS,etaSS,epsSS vectors

!!!Construct Xss = (Va,Vna,moments,logp,logY,logI,logN,logA,kprime)
XSS(:) = 0.0

!insert Va
ct=0
do zct=1,znum
do kct=1,knum
	ct = ct + 1
	Xss(ct) = VaSS(zct,kct)
end do !zct
end do !kct

!insert Vna
do zct=1,znum
do kct=1,knum
	ct = ct + 1
	Xss(ct) = VnaSS(zct,kct)
end do !zct
end do !kct

!insert moments 
do zct=1,znum
do momct=1,momuse
	ct = ct + 1
	Xss(ct) = momSS(zct,momct)
end do !zct
end do !momct

!insert aggregates
ct = ct +1; Xss(ct) = log(pSS)
ct = ct +1; Xss(ct) = log(Yss)
ct = ct +1; Xss(ct) = log(Iss)
ct = ct +1; Xss(ct) = log(Nss)
ct = ct +1; Xss(ct) = log(aSS)

!insert kprime policies
do zct=1,znum
ct = ct +1; Xss(ct) = kprimeSS(zct,1) !note that kprimeSS doesn't vary by k, so this just selects the first entry
end do !zct

!insert zeros into expectational errors and agg shock series
etaSS(:) = 0.0
epsSS(:) = 0.0

!call Fsys for SS vals
doMICROsim = 1
FSS(:) = 0.0
MICROSS(:) = 0.0
call Fsys(Xss,Xss,etaSS,epsSS,FSS,MICROSS)

!!!!NOW, DO NUMERICAL LINEARIZATION AROUND SS
doMICROsim = 0
!Jacobian wrt Xvec
F1(:,:) = 0.0
do ct=1,numX
	if (mod(ct,20)==0) then
		write(*,*) "Computing derivative for X state ",ct,"."
	end if	
	Xvec = Xss
	h = Xvec(ct)*diffstep
	if (abs(h)<1.0e-10) then
		h = diffstep
	end if
	Xvec(ct) = Xvec(ct)+h
	call Fsys(Xvec,Xss,etaSS,epsSS,Fvec,MICROvec)
	F1(:,ct) = F1(:,ct)+(Fvec-FSS)/h
end do !ct

!Jacobian wrt Xmin1vec
F2(:,:) = 0.0
do ct=1,numX
	if (mod(ct,20)==0) then
		write(*,*) "Computing derivative for Xmin1 state ",ct,"."
	end if	
	Xmin1vec = Xss
	h = Xmin1vec(ct)*diffstep
	if (abs(h)<1.0e-10) then
		h = diffstep
	end if
	Xmin1vec(ct) = Xmin1vec(ct)+h
	call Fsys(Xss,Xmin1vec,etaSS,epsSS,Fvec,MICROvec)
	F2(:,ct) = F2(:,ct)+(Fvec-FSS)/h
end do !ct

!Jacobian wrt eta
F3(:,:) = 0.0
do ct=1,numeta
	if (mod(ct,20)==0) then
		write(*,*) "Computing derivative for eta state ",ct,"."
	end if	
	etavec = etaSS
	h = etavec(ct)*diffstep
	if (abs(h)<1.0e-10) then
		h = diffstep
	end if
	etavec(ct) = etavec(ct)+h
	call Fsys(Xss,Xss,etavec,epsSS,Fvec,MICROvec)
	F3(:,ct) = F3(:,ct)+(Fvec-FSS)/h
end do !ct


!Jacobian wrt eps
F4(:,:) = 0.0
write(*,*) "Computing derivative for eps shock."
epsvec = epsSS
h = epsvec(1)*diffstep
if (abs(h)<1.0e-10) then
	h = diffstep
end if
epsvec = epsvec+h;
call Fsys(Xss,Xss,etaSS,epsvec,Fvec,MICROvec)
F4(:,1) = F4(:,1)+(Fvec-FSS)/h

!XSS.txt
open(8,file="XSS.txt")
do ct=1,numX
write(8,*) XSS(ct)
end do !ct
close(8)

!FSS.txt
open(8,file="FSS.txt")
do ct=1,numX
write(8,*) FSS(ct)
end do !ct
close(8)

!MICROSS.txt
open(8,file="MICROSS.txt")
do ct=1,nummicro
write(8,*) MICROSS(ct)
end do !ct
close(8)

!F1.txt
open(8,file="F1.txt")
do ct=1,numX
do ct2=1,numX
write(8,*) F1(ct,ct2)
end do !ct2
end do !ct
close(8)

!F2.txt
open(8,file="F2.txt")
do ct=1,numX
do ct2=1,numX
write(8,*) F2(ct,ct2)
end do !ct2
end do !ct
close(8)

!F3.txt
open(8,file="F3.txt")
do ct=1,numX
do ct2=1,numeta
write(8,*) F3(ct,ct2)
end do !ct2
end do !ct
close(8)

!F4.txt
open(8,file="F4.txt")
do ct=1,numX
do ct2=1,numeps
write(8,*) F4(ct,ct2)
end do !ct2
end do !ct
close(8)

!ergdistz.txt
open(8,file="ergdistz.txt")
do ct=1,znum
write(8,*) ergdistz(ct)
end do !ct
close(8)

!now, call the gensys solver in MATLAB (which itself uses a mex file)
call system('/Applications/MATLAB_R2016b.app/bin/matlab -nosplash -nodisplay -r "call_gensys"')

finish = omp_get_wtime()
write(*,*) "Finished model solution at ",finish-start," seconds."
write(13,*) "Finished model solution at ",finish-start," seconds."


!now, call the MATLAB file which performs the simulation (which itself uses a mex file)
call system('/Applications/MATLAB_R2016b.app/bin/matlab -nosplash -nodisplay -r "call_simulate"')

finish = omp_get_wtime()
write(*,*) "Finished model unconditional simulation at ",finish-start," seconds."
write(13,*) "Finished model unconditional simulation at ",finish-start," seconds."


if (doIRF==1) then
!now, perform the IRF simulation
call system('/Applications/MATLAB_R2016b.app/bin/matlab -nosplash -nodisplay -r "call_simulate_IRF"')

finish = omp_get_wtime()
write(*,*) "Finished IRF simulation at ",finish-start," seconds."
write(13,*) "Finished IRF simulation at ",finish-start," seconds."
end if


if (doMICROsimstore==1) then
open(8,file="distkzsim.txt")
close(8,status='DELETE')	
	
!now, extract the micro moments
epssim(:) = 0.0
Xsim(:,:) = 0.0
open(8,file="Xsim.txt")
open(9,file="epssim.txt")
do t=1,numper
	read(9,*) epssim(t)
do ct=1,numX
read(8,*) Xsim(ct,t)
end do !ct
end do !t
close(8)
close(9)
write(*,*) "Finished reading in Xsim and epssim."

doMICROsim=1
MICROsimstore(:,:) = 0.0
do t=2,numper
	write(*,*) "t = ",t
	call Fsys(Xsim(:,t),Xsim(:,t-1),etaSS,epssim(t),Fvec,MICROsimstore(:,t-1))
end do !t

finish = omp_get_wtime()
write(*,*) "Finished micro moment simulation at ",finish-start," seconds."
write(13,*) "Finished micro moment simulation at ",finish-start," seconds."
end if 

!!!WRITE SOME FILES
!z0.txt
open(8,file="z0.txt")
open(9,file="pr_mat_z.txt")
do zct=1,znum
write(8,*) z0(zct)
write(9,*) pr_mat_z(zct,:)
end do !zct
close(8)
close(9)

!k0.txt
open(8,file="k0.txt")
do kct=1,knum
write(8,*) k0(kct)
end do !kct
close(8)

if (doMICROsimstore==1) then
open(8,file="MICROsim.txt")
do t=1,numper
do momct=1,nummicro
write(8,*) MICROsimstore(momct,t)
end do !momct
end do !t
close(8)
end if

finish = omp_get_wtime()
write(*,*) "Finished program at ",finish-start," seconds."
write(13,*) "Finished program at ",finish-start," seconds."
close(13)



contains


subroutine discretize_simulate()
implicit none

double precision :: asimgrid(anum),shockprobdenom,shockprobIRF
integer :: act,zct,statect,ct,t,aprimect,shockct

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
!draw random seeds
call random_seed(size=seeddim)

!insert random seed into seedarray
allocate(seedarray(seeddim))
do ct=1,seeddim
    seedarray(ct) = seedint + ct
end do !ct
call random_seed(put=seedarray)

!draw random shocks
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
!END UNCONDITIONAL SIMULATION

!START IRF SIMULATION
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
    
    if (maxval(abs(ergdista-ergdistaold))<disterrortol) exit
    
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
!!END IRF SIMULATION


end subroutine discretize_simulate


subroutine Fsys(X,Xmin1,eta,eps,F,MICROsim)
implicit none

!inputs & outputs
double precision :: X(numX),Xmin1(numX),eta(numeta),eps(numeps)

double precision :: F(numX),MICROsim(nummicro)

!other variables
integer :: ct,zct,kct,momct,zprimect,insertct
double precision :: Vmin1(znum,knum),Vamin1(znum,knum),Va(znum,knum),Vnamin1(znum,knum),&
	Vna(znum,knum),mommin1(znum,momuse),mom(znum,momuse),kprimevec(znum),kprimeval,Vaval,Vnaval,&
	Vnextval,zval,kval,xival,xistarsimp(znum,nsimp+1),kprimesimp(znum,nsimp+1),RHSVa(znum,knum),&
	RHSVna(znum,knum),RHSV(znum,knum),momnextmat(znum,momuse),padjust,adjustmat(znum,nsimp+1,2),&
	Yactual,Iactual,Kprimeactual,Nactual,Cactual,wgt,etaA(znum,knum),etaNA(znum,knum),epsA,ival,&
	etakprime(znum),RHSkprime(znum),Vprimeval,Va2min1(znum,knum),Vna2min1(znum,knum)


double precision :: FRHS(numX),FLHS(numX)

FLHS(:) = -100000.0
FRHS(:) = 100000.0
F(:) = 100000.0
insertct = 0

!!!!extract values from X,Xmin1 vectors

!Va
ct=0
do zct=1,znum
do kct=1,knum
	ct = ct + 1
	Va(zct,kct) = X(ct)
	Vamin1(zct,kct) = Xmin1(ct)	
end do !zct
end do !kct

!Vna
do zct=1,znum
do kct=1,knum
	ct = ct + 1
	Vna(zct,kct) = X(ct)
	Vnamin1(zct,kct) = Xmin1(ct)	
end do !zct
end do !kct

!moments
do zct=1,znum
do momct=1,momuse
	ct = ct + 1
	mom(zct,momct) = X(ct)
	mommin1(zct,momct) = Xmin1(ct)	
end do !zct
end do !momct

!aggregates
ct=ct+1; p = exp(X(ct)); pmin1 = exp(Xmin1(ct)); 
ct=ct+1; Y= exp(X(ct));  Ymin1 = exp(Xmin1(ct))
ct=ct+1; I = exp(X(ct)); Imin1 = exp(Xmin1(ct))
ct=ct+1; N = exp(X(ct)); Nmin1 = exp(Xmin1(ct))
ct=ct+1; a = exp(X(ct)); amin1 = exp(Xmin1(ct))

!capital policy
do zct=1,znum
	ct=ct+1; kprimevec(zct)=Xmin1(ct);
end do !zct


!note that prices imply wages
w = phi/p; wmin1 = phi/pmin1

!!!!!extract eta values

!etaA
ct=0
do zct=1,znum
do kct=1,knum
	ct=ct+1
	etaA(zct,kct) = eta(ct)
end do !kct
end do !zct

!etaNA
do zct=1,znum
do kct=1,knum
	ct=ct+1
	etaNA(zct,kct) = eta(ct)
end do !kct
end do !zct

!eta kprime
do zct=1,znum
	ct = ct+1
	etakprime(zct) = eta(ct)
end do !zct

!extract eps value
epsA = eps(1)
!!!!!Goal #1: Compute implied RHS of Bellman eqns

!first, construct the continuation value function
do zct=1,znum
do kct=1,knum
	xival = (Va(zct,kct) - Vna(zct,kct))/phi
	V(zct,kct) = -1.0 * phi*expecxi(xival)+ &
		cdfxi(xival)*Va(zct,kct) +&
		(1.0-cdfxi(xival))*Vna(zct,kct)
end do !kct
end do !zct

!then, spline the continuation value function
V2(:,:) = 0.0
do zct=1,znum
	call spline(k0,V(zct,:),knum,dble(1.0e30),dble(1.0e30),V2(zct,:))
end do !zct

!now, construct the kprime FOC's
RHSkprime(:) = 0.0

do zct=1,znum
do zprimect=1,znum
	RHSzct = zct
	kprimeval = kprimevec(zct)
	call splintprime(k0,V(zprimect,:),V2(zprimect,:),knum,kprimeval,Vprimeval)
	RHSkprime(zct) = RHSkprime(zct) + beta * pr_mat_z(RHSzct,zprimect) * Vprimeval
end do !zprimect
end do !zct


RHSVa(:,:) = 0.0
RHSVna(:,:) = 0.0

!must now construct RHS of Bellman eqn's given the optimal policies
do zct=1,znum
do kct=1,knum
    
    RHSzct = zct
    
    !what states are there here?
    zval = z0(zct)
    kval = k0(kct)
    
    !construct value when not adjusting
    Vnaval = p * freduced(zval,kval)
    kprimeval = max((1.0-delta)*kval,kmin) !no adjustment in capital
        
    !add continuation value when not adjusting
    do zprimect=1,znum

        !evaluate continuation at kbarfcstind
        call splint(k0,V(zprimect,:),V2(zprimect,:),knum,kprimeval,Vnextval)
        Vnaval = Vnaval + beta * pr_mat_z(RHSzct,zprimect) * Vnextval   
        
    end do !zprimect
    
    !store RHS of Vna
    RHSVna(zct,kct) = Vnaval
    
    !construct value when adjusting
    kprimeval = kprimevec(zct)
    Vaval = p*freduced(zval,kval)!current period return
    Vaval = Vaval - p * ( kprimeval - ( 1.0 - delta ) * kval )
    
        
    !add continuation value when not adjusting
    do zprimect=1,znum

        !evaluate continuation at kbarfcstind
        call splint(k0,V(zprimect,:),V2(zprimect,:),knum,kprimeval,Vnextval)
        Vaval = Vaval + beta * pr_mat_z(RHSzct,zprimect) *  Vnextval   
        
    end do !zprimect        

    !store RHS of Va    
    RHSVa(zct,kct) = Vaval
    
end do !kct
end do !zct

!Fill in the Bellman equation portions of FLHS and FRHS

!Va block
do zct=1,znum
do kct=1,knum
	insertct = insertct + 1 
	FLHS(insertct) = Vamin1(zct,kct)
	FRHS(insertct) = RHSVa(zct,kct) + etaA(zct,kct)
end do !kct
end do !zct

!Vna block
do zct=1,znum
do kct=1,knum
	insertct = insertct + 1 
	FLHS(insertct) = Vnamin1(zct,kct)
	FRHS(insertct) = RHSVna(zct,kct) + etaNA(zct,kct)
end do !kct
end do !zct


!!!!!Goal #3: Construct next period moments

!must construct the optimal next-period capital and adjustment probs at simpson nodes
!in order to conduct the required quadrature for the next period

!first, spline the Va and Vna values, then construct xival
Va2min1(:,:) = 0.0
Vna2min1(:,:) = 0.0
do zct=1,znum
	call spline(k0,Vamin1(zct,:),knum,dble(1.0e30),dble(1.0e30),Va2min1(zct,:))
	call spline(k0,Vnamin1(zct,:),knum,dble(1.0e30),dble(1.0e30),Vna2min1(zct,:))
end do !zct

xistarsimp(:,:) = 0.0
kprimesimp(:,:) = 0.0
do zct=1,znum
    do kct=1,nsimp+1
    
        !extract state values
        zval = z0(zct)
        kval = simpnodes(kct)
        
        kprimeval = kprimevec(zct) !extract policy determined above, picking first element b/c only depends on z
        
		call splint(k0,Vamin1(zct,:),Va2min1(zct,:),knum,kval,Vaval)
		call splint(k0,Vnamin1(zct,:),Vna2min1(zct,:),knum,kval,Vnaval)
        
        !now that the value when adjusting and when not adjusting are determined, create and store total value
        xival = (Vaval-Vnaval)/phi  !cutoff AC value
        
	xistarsimp(zct,kct) = xival
    kprimesimp(zct,kct) = kprimevec(zct)
        
    end do !kct
end do !zct


!find the correct values of rhomat, intvec based on "min1" moments
mommat = mommin1 !integrate against current distribution
call findrhoSR1wolfe(1) !this sets rhomat and intvec correctly, to match values of mommin1

!now, compute next period moments

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


!now, insert the moments into the FRHS and FLHS system
do zct=1,znum
do momct=1,momuse
	insertct = insertct + 1
	FLHS(insertct) = mom(zct,momct)
	FRHS(insertct) = momnextmat(zct,momct)
end do !momct
end do !zct

!!!!!Goal #4: Construct this period aggregates

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

!insert aggregates into the system

insertct = insertct + 1; FLHS(insertct) = 1.0/p; FRHS(insertct) = Cactual
insertct = insertct + 1; FLHS(insertct) = Y;     FRHS(insertct) = Yactual
insertct = insertct + 1; FLHS(insertct) = I;     FRHS(insertct) = Iactual
insertct = insertct + 1; FLHS(insertct) = N;     FRHS(insertct) = Nactual

!!!!!Goal #5: Construct agg shock eqn
insertct = insertct + 1
FLHS(insertct) = log(a)
FRHS(insertct) = rhoa * log(amin1) + epsA

!!!!!Goal #6: insert policy FOC's into the system
do zct=1,znum
	insertct=insertct+1
	FLHS(insertct) = p
	FRHS(insertct) = RHSkprime(zct) + etakprime(zct)
end do !zct

!!!NOW, TAKE DIFFERENCE
F = FLHS-FRHS

!!!NOW, IF YOU HAVE INDICATED THAT MICRO MOMENTS MUST BE COMPUTED, DO THOSE NOW
if (doMICROsim==0) then
	
	MICROsim(:) = 1000000.0
	
else if (doMICROsim==1) then
		
    !now, compute MICRO moments, which needs to occur with this period's rhomat, mommat, intvec, etc.   
     !MICROsim(nummicro,numper)
    ! 1 = i/k
    ! 2 = stdev(i/k)
    ! 3 = P(inaction)
    ! 4 = P(i/k>=0.2)
    ! 5 = P(i/k<=-0.2)
    ! 6 = P(i/k > 0)
    ! 7 = P(i/k < 0)
    
    MICROsim(:) = 0.0
    
    do zct=1,znum
    do kct=1,nsimp+1
        
        zval = z0(zct)
        kval = simpnodes(kct)
        kprimeval = kprimesimp(zct,kct)
        ival = kprimeval - (1.0-delta)*kval !investment conditional upon investment
        xival = xistarsimp(zct,kct)
        padjust = cdfxi(xival)
        
        !investment rate
        MICROsim(1) = MICROsim(1) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust * (ival / kval)
        
        !investment rate squared - for stdev construction
        MICROsim(2) = MICROsim(2) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust * ( (ival / kval)**2.0 )
        
        !P(inaction)
        MICROsim(3) = MICROsim(3) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * (1.0 - padjust)
        
        !P(pos. spike)
        if ((ival/kval)>=0.2) then
            MICROsim(4) = MICROsim(4) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust
        end if
        
        !P(neg. spike)
        if ((ival/kval)<=-0.2) then
            MICROsim(5) = MICROsim(5) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust
        end if
        
        !P(pos. invest)
        if ((ival/kval)>0.0) then
            MICROsim(6) = MICROsim(6) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust
        end if
        
        !P(neg. invest)
        if ((ival/kval)<0.0) then
            MICROsim(7) = MICROsim(7) + (simpweights(kct)*ergdistz(zct)*Fk(kval,zct)/intvec(zct)) * padjust 
        end if
        
    
    end do !kct
    end do !zct
        
        
    !now, convert squared investment moment to stdev
    MICROsim(2) = sqrt(MICROsim(2) - (MICROsim(1)**2.0))

	open(8,file="distkzsim.txt",position='APPEND')
	do zct=1,znum
	do kct=1,kdensenum
		kval=kdense0(kct)
		write(8,*) Fk(kval,zct)
	end do !kct
	end do !zct
	close(8)
	
end if

end subroutine Fsys


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

!NOTE THE USE OF amin1,w here!

exponentnu = 1.0 / (1.0 - nu)

freduced = ( nu ** (nu * exponentnu) ) * (1.0 - nu) * ( w ** ( -1.0 * nu * exponentnu ) )
freduced = freduced * (z ** exponentnu) * (amin1 ** exponentnu) * ( k ** (alpha * exponentnu) )
end function freduced

double precision function nreduced(z,k)
implicit none
double precision :: z,k

double precision :: exponentnu

!NOTE THE USE OF amin1,w here!

exponentnu = 1.0 / (1.0 - nu)

nreduced = ( nu **  exponentnu ) * ( w ** ( -1.0 * exponentnu ) )
nreduced = nreduced * (z ** exponentnu) * (amin1 ** exponentnu) * ( k ** (alpha * exponentnu) )
end function nreduced

double precision function yreduced(z,k)
implicit none
double precision :: z,k

double precision :: exponentnu

!NOTE THE USE OF amin1,w here!

exponentnu = 1.0 / (1.0 - nu)

yreduced = ( (nu/w) ** (nu * exponentnu) ) * ( (amin1 * z ) ** exponentnu  )
yreduced = yreduced * ( k ** ( alpha * exponentnu ) )

end function yreduced



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

!check for issues and return SS coefficients if there are optimization failures
if ((nannewfunc.eqv..TRUE.).or.(nannewgrad.eqv..TRUE.)) then
    write(*,*) "ISSUE WITH OPTIMIZATION, RETURNING SS VALUES"
    rhomat = rhoSS
    intvec = intSS
end if


end subroutine findrhoSR1wolfe



end program kt_winberry