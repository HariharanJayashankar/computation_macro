!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kt_ss.f90
!
! Fortran code for the KS solution of the Khan and Thomas (2008) 
! model in the steady-state for use in the XPA solution.
!
! 'Alternative Methods for Solving Heterogeneous Firm Models'
! Stephen Terry (2015)
!
! This Version : 12/18/15
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program kt_ss
use base_lib
implicit none

character(len=*), parameter :: basedir="."

integer :: znum,knum,kdensenum,waitsec,maxvfit,maxaccelit,maxpit,maxergit,zct,kct,&
    zprimect,vfct,accelit,ct,zdumct,zvalpos,avalpos,kvalpos,pvalpos,&
    zctpos,nparam,numconstants,ergct,kprimeind,kprimeindnoadj,piter,anum,act
    
double precision :: start,finish,alpha,nu,phi,xibar,beta,delta,kmin,kmax,rhoz,sigmaz,&
    vferrortol,kprimeerrortol,xistarerrortol,perrortol,pgain,pcorrect,brenttol,ergerrortol,&
    pval,wval,aval,zval,kval,kprimeval,Vnaval,Vaval,Vval,xival,Vnextval,vferror,kprimeerror,&
    xistarerror,ergerror,weight,weightnoadj,Kprimeactual,Yactual,Cactual,Iactual,yval,ival,&
    perror,pinit,rhoa,sigmaa,nstdeva,nstdevz

double precision, allocatable :: k0(:),z0(:),pr_mat_z(:,:),V(:,:),Vold(:,:),V2old(:,:),&
    kprime(:,:),kprimeold(:,:),xistar(:,:),xistarold(:,:),param(:),Va(:,:),Vna(:,:),constantvec(:),&
    distkz(:,:),distkzold(:,:),kdense0(:),kprimedense(:,:),xistardense(:,:),Kbarnoaggstore(:),&
    pnoaggstore(:),a0(:),pr_mat_a(:,:)


!!program control

call cpu_time(start)
open(13,file="kt_ss.txt")

!!!!INSERT PARAMETERS
alpha = 0.256
nu = 0.640
phi = 2.4
xibar = 0.0083
beta = 0.977
delta = 0.069

aval = 1.0; 

knum = 10
kmin= 0.1
kmax = 8.0

nstdevz = 2.0 !number of stdevs to cover in discretization of idio prod
znum = 5; 
rhoz = 0.859; 
sigmaz = 0.022; 

nstdeva = 2.0 !number of stdevs to cover in discretization of agg prod
anum = 5; 
rhoa = 0.859;
sigmaa = 0.014;

kdensenum = 50

waitsec=10

maxvfit = 1000
vferrortol=1e-4
kprimeerrortol = 1e-4
xistarerrortol = 1e-4

maxaccelit = 50

maxpit=2000
perrortol = 1e-4
pgain=0.05
pcorrect=1.5

pinit = 2.4; !initial guess for p

brenttol = 1e-6

maxergit = 5000
ergerrortol = 1e-6; !tolerance on finding error distribution

nparam = 2*znum*knum + znum*znum + knum + 12

numconstants = 2

!!do allocations
allocate(k0(knum),z0(znum),pr_mat_z(znum,znum),V(znum,knum),Vold(znum,knum),&
V2old(znum,knum),kprime(znum,knum),kprimeold(znum,knum),xistar(znum,knum),&
xistarold(znum,knum),param(nparam),Va(znum,knum),Vna(znum,knum),constantvec(numconstants),&
distkz(znum,kdensenum),distkzold(znum,kdensenum),kdense0(kdensenum),kprimedense(znum,kdensenum),&
xistardense(znum,kdensenum),Kbarnoaggstore(anum),pnoaggstore(anum),a0(anum),pr_mat_a(anum,anum))

!save constants for plotting via MATLAB
constantvec = (/ xibar,delta /)

!!set up idio capital grids
call linspace(k0,log(kmin),log(kmax),knum); k0=exp(k0);
call linspace(kdense0,log(kmin),log(kmax),kdensenum); kdense0=exp(kdense0);


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



!!LOOP OVER A
do act=1,anum
aval = a0(act)
write(*,*) " NOW DOING A = ",aval



!!!now, initialize the value and policy functions
Vold(:,:) = 0.0
V2old(:,:) = 0.0
xistarold(:,:) = 0.5

do zct=1,znum
    kprimeold(zct,:) = k0    
end do !zct


!!!!!!!!!!!!!INSERT LOOP OVER P HERE
pval = pinit

do piter=1,maxpit

wval = phi/pval


!!!!!!!!!!!!!GIVEN P, DO HOWARD POLICY ACCELERATION/VFI


do vfct=1,maxvfit
    Vold(:,:) = 0.0; V2old(:,:) = 0.0
    
    do accelit=1,maxaccelit
        
        do zct=1,znum
        do kct=1,knum
            
            zval = z0(zct); kval = k0(kct);
            
            !determine Vnaval
            Vnaval = pval * freduced(zval,kval,aval,wval,alpha,nu)
            kprimeval = max((1.0-delta)*kval,kmin)
            
            do zprimect=1,znum
                call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprimeval,Vnextval)
                Vnaval = Vnaval + beta * pr_mat_z(zct,zprimect) * Vnextval
            end do !zprimect
            
            !determine Vaval
            Vaval = pval * freduced(zval,kval,aval,wval,alpha,nu)
            kprimeval = kprimeold(zct,kct)
            Vaval = Vaval - pval*(kprimeval - (1.0-delta)*kval)
            
            do zprimect=1,znum
                call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprimeval,Vnextval)
                Vaval = Vaval + beta * pr_mat_z(zct,zprimect) * Vnextval
            end do !zprimect
            
            !determined Vval
            xival = (Vaval - Vnaval)/phi
            Vval = -1.0 * phi * expecxi(xival,xibar)&
                + cdfxi(xival,xibar) * Vaval + (1.0-cdfxi(xival,xibar)) * Vnaval
                
            V(zct,kct) = Vval
            
        end do !kct
        end do !zct
        
        Vold = V
        
        do zct=1,znum
        do kct=1,knum
            call spline(k0,Vold(zct,:),knum,dble(1.0e30),dble(1.0e30),V2old(zct,:))
        end do !zct
        end do !kct
        
    end do !accelit
    
    
    !!!!now that the acceleration is done, do traditional VFI, with actual optimization via Brent
    
    !first,insert values into "param" vector for use in optimization
    ct=0
    do zdumct=1,znum; !note that we use zdumct here rather than zct because zct (and act) are reserved
    do kct=1,knum
        ct=ct+1
        param(ct) = Vold(zdumct,kct)
    end do    
    end do    

    do zdumct=1,znum
    do kct=1,knum
        ct=ct+1
        param(ct) = V2old(zdumct,kct)
    end do    
    end do    

    do zdumct=1,znum
    do zprimect=1,znum
        ct=ct+1
        param(ct) = pr_mat_z(zdumct,zprimect)
    end do    
    end do    

    do kct=1,knum
        ct=ct+1
        param(ct) = k0(kct)
    end do 

    ct=ct+1; param(ct) = alpha;
    ct=ct+1; param(ct) = beta;
    ct=ct+1; param(ct) = delta;
    ct=ct+1; param(ct) = nu;
    ct=ct+1; param(ct) = phi;
    ct=ct+1; param(ct) = zval; zvalpos = ct;
    ct=ct+1; param(ct) = aval; avalpos = ct;
    ct=ct+1; param(ct) = kval; kvalpos = ct;
    ct=ct+1; param(ct) = pval; pvalpos = ct;
    ct=ct+1; param(ct) = dble(zct); zctpos = ct;
    
    ct=ct+1; param(ct) = dble(znum);
    ct=ct+1; param(ct) = dble(knum);
    
    
    !now, loop over states
    
    do zct=1,znum
    do kct=1,knum
        
        !determine states
        zval = z0(zct); kval = k0(kct);
        
        !this block determines Vnaval
        Vnaval = pval * freduced(zval,kval,aval,wval,alpha,nu)
        kprimeval = max((1.0-delta)*kval,kmin)
        
        do zprimect=1,znum
            call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprimeval,Vnextval)
            Vnaval = Vnaval + beta * pr_mat_z(zct,zprimect) * Vnextval
        end do !zprimect
        Vna(zct,kct) = Vnaval
        
        
        !this block determines Vaval
        param(zvalpos) = zval;
        param(avalpos) = aval;
        param(kvalpos) = kval;
        param(pvalpos) = pval;
        param(zctpos) = zct;
        
        Vaval = brentparam(k0(1),k0(knum/2),k0(knum),fVa,brenttol,kprimeval,param,nparam)
        Vaval = -1.0 * Vaval
        
        Va(zct,kct)=Vaval
        kprime(zct,kct) = kprimeval
        
        xival = (Vaval-Vnaval)/phi
        xistar(zct,kct) = xival
        
        !now, process this info to find the value function
        Vval = -1.0*phi*expecxi(xival,xibar)&
            + cdfxi(xival,xibar) * Vaval + (1.0 - cdfxi(xival,xibar))*Vnaval
            
        V(zct,kct) = Vval
        
        
    end do !kct
    end do !znum
    
    vferror = maxval(abs( (log(V)-log(Vold) )) )    
    kprimeerror = maxval(abs(log(kprime)-log(kprimeold)))
    xistarerror = maxval(abs(log(xistar)-log(xistarold)))
    if (mod(vfct,5)==0.and.mod(piter,10)==0) then
    write(13,*) "VF iter = ",vfct,"VF error = ",vferror
    write(13,*) "VF iter = ",vfct,"Kprime error = ",kprimeerror
    write(13,*) "VF iter = ",vfct,"Xistar error = ",xistarerror
    write(13,*) " "
    
    write(*,*) "VF iter = ",vfct,"VF error = ",vferror
    write(*,*) "VF iter = ",vfct,"Kprime error = ",kprimeerror
    write(*,*) "VF iter = ",vfct,"Xistar error = ",xistarerror
    write(*,*) " "
    end if
    
    if ((kprimeerror<kprimeerrortol .and. xistarerror<xistarerrortol).or.(vferror<vferrortol)) exit
    
    Vold=V
    kprimeold = kprime
    xistarold = xistar
end do !vfct



!!!!!!!!!GIVEN VALUE FUNCTION AND POLICIES, NOW COMPUTE IMPLIED ERGODIC AND PRICES

!evaluate, just once, the optimal policies implied via above at kdense0
do zct=1,znum
do kct=1,kdensenum
    
    !determine states
    zval = z0(zct)
    kval = kdense0(kct)
    
    !this block determines Vnaval
    Vnaval = pval * freduced(zval,kval,aval,wval,alpha,nu)
    kprimeval = max((1.0-delta)*kval,kmin)
    
    do zprimect=1,znum
        call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprimeval,Vnextval)
        Vnaval = Vnaval + beta * pr_mat_z(zct,zprimect) * Vnextval
    end do !zprimect

    !this block determines Vaval
    param(zvalpos) = zval;
    param(avalpos) = aval;
    param(kvalpos) = kval;
    param(pvalpos) = pval;
    param(zctpos) = zct;
    
    Vaval = brentparam(k0(1),k0(knum/2),k0(knum),fVa,brenttol,kprimeval,param,nparam)
    Vaval = -1.0 * Vaval
    

    kprimedense(zct,kct) = kprimeval
    xistardense(zct,kct) = (Vaval-Vnaval)/phi
    

end do !kct
end do !zct



!initialize distributions
distkzold(:,:) = 0.0
distkzold(3,kdensenum/2) = 1.0; distkzold = distkzold / sum(distkzold)

distkz(:,:) = 0.0


do ergct=1,maxergit
    
    do zct=1,znum
    do kct=1,kdensenum
        
        !extract policies from above & determine states
        kval = kdense0(kct)
        kprimeval = kprimedense(zct,kct)
        xival = xistardense(zct,kct)
        
        !determine which interval the policy is in, set up weight
        kprimeind = kct
        call hunt(kdense0,kdensenum,kprimeval,kprimeind)
            
        if (kprimeind<=0) then
            weight = 0.0
            kprimeind = 1
        else if (kprimeind>=1.and.kprimeind<=(kdensenum-1)) then
            weight = (kprimeval-kdense0(kprimeind))/(kdense0(kprimeind+1)-kdense0(kprimeind))
        else if (kprimeind>=kdensenum) then
            weight = 1.0
            kprimeind = kdensenum-1
        end if
        
        !determined which interval noadjustment is in, set up that weight
        kprimeindnoadj = kct
        kprimeval = max((1.0-delta)*kval,kmin)
        call hunt(kdense0,kdensenum,kprimeval,kprimeindnoadj)
            
        if (kprimeindnoadj<=0) then
            weightnoadj = 0.0
            kprimeindnoadj = 1
        else if (kprimeindnoadj>=1.and.kprimeindnoadj<=(kdensenum-1)) then
            weightnoadj = (kprimeval-kdense0(kprimeindnoadj))&
            /(kdense0(kprimeindnoadj+1)-kdense0(kprimeindnoadj))
        else if (kprimeindnoadj>=kdensenum) then
            weightnoadj = 1.0
            kprimeindnoadj = kdensenum-1
        end if
        
        
        do zprimect=1,znum
            
            !those that adjust and go to kprimeind
            distkz(zprimect,kprimeind) = distkz(zprimect,kprimeind) &
                + pr_mat_z(zct,zprimect) * cdfxi(xival,xibar) * (1.0-weight) * distkzold(zct,kct)
            
            !those that adjust and go to kprimeind+1
            distkz(zprimect,kprimeind+1) = distkz(zprimect,kprimeind+1) &
                + pr_mat_z(zct,zprimect) * cdfxi(xival,xibar) * weight * distkzold(zct,kct)
            
            !those that don't adjust and go to kprimeindnoadj
            distkz(zprimect,kprimeindnoadj) = distkz(zprimect,kprimeindnoadj) &
                + pr_mat_z(zct,zprimect) * (1.0 - cdfxi(xival,xibar)) * (1.0-weightnoadj) * distkzold(zct,kct)
                
            !those that don't adjust and go to kprimeindnoadj+1
            distkz(zprimect,kprimeindnoadj+1) = distkz(zprimect,kprimeindnoadj+1) &
                + pr_mat_z(zct,zprimect) * (1.0 - cdfxi(xival,xibar)) * weightnoadj * distkzold(zct,kct)
            
        end do !zprimect
        
        
    end do 
    end do !zct
    
    
    ergerror = maxval(abs(distkz-distkzold))
    distkz=distkz/sum(distkz)
    
    if (mod(ergct,15)==0.and.mod(piter,10)==0) then
        write(13,*) "Erg iter = ",ergct,"Erg error = ",ergerror
        write(*,*) "Erg iter = ",ergct,"Erg error = ",ergerror
    end if
    
    if (ergerror<ergerrortol) exit
    
    distkzold = distkz
    distkz(:,:) = 0.0

end do !ergct


!!!!!!!!!!!!!!END ONE ITERATION WITH ONE P, COMPUTE ERRORS AND UPDATE P

Kprimeactual = 0.0
Yactual = 0.0
Cactual = 0.0
Iactual = 0.0

do zct=1,znum
do kct=1,kdensenum
    
    zval = z0(zct); kval = kdense0(kct);
    yval = yreduced(zval,kval,aval,wval,alpha,nu)
    kprimeval = kprimedense(zct,kct)
    ival = kprimeval - (1.0-delta)*kval
    xival = xistardense(zct,kct)
    
    Yactual = Yactual + distkz(zct,kct) * yval
    
    Kprimeactual = Kprimeactual + distkz(zct,kct) * cdfxi(xival,xibar) * kprimeval &
        + distkz(zct,kct) * (1.0 - cdfxi(xival,xibar)) * (1.0 - delta) * kval
        
    Iactual = Iactual + distkz(zct,kct) * cdfxi(xival,xibar) * (kprimeval - (1.0 -delta)*kval) 
    
end do !kct
end do !zct

Cactual = Yactual - Iactual

perror = pval - (1.0/Cactual)

if (mod(piter,10)==0) then
    write(*,*) "P iter = ",piter,"P error = ",perror; write(*,*) " "
    write(13,*) "P iter = ",piter,"P error = ",perror; write(13,*) " "
    
end if

if (abs(perror)<perrortol) exit
if (Cactual > 0.0) pval = pval + pgain * (1.0/Cactual - pval)
if (Cactual <= 0.0) pval = pval * pcorrect


end do !piter

write(*,*) "Equilibrium K = ",Kprimeactual
write(13,*) "Equilibrium K = ",Kprimeactual

write(*,*) "Equilibrium p = ",pval
write(13,*) "Equilibrium p = ",pval

write(*,*) " "
write(13,*) " "

Kbarnoaggstore(act) = Kprimeactual
pnoaggstore(act) = pval

write(*,*) "END OF A = ",aval
write(*,*)
end do !act

open(8,file="Kbarnoagg.txt")
do act=1,anum
write(8,*) Kbarnoaggstore(act)
end do
close(8)

open(8,file="pnoagg.txt")
do act=1,anum
write(8,*) pnoaggstore(act)
end do
close(8)

open(8,file="k0.txt")
do kct=1,knum
write(8,*) k0(kct)
end do !kct
close(8)

open(8,file="V.txt")
do zct=1,znum; do kct=1,knum
write(8,*) V(zct,kct)
end do; end do
close(8)

open(8,file="Va.txt")
do zct=1,znum; do kct=1,knum
write(8,*) Va(zct,kct)
end do; end do
close(8)

open(8,file="Vna.txt")
do zct=1,znum; do kct=1,knum
write(8,*) Vna(zct,kct)
end do; end do
close(8)

open(8,file="kprime.txt")
do zct=1,znum; do kct=1,knum
write(8,*) kprime(zct,kct)
end do; end do
close(8)

open(8,file="xistar.txt")
do zct=1,znum; do kct=1,knum
write(8,*) xistar(zct,kct)
end do; end do
close(8)

open(8,file="padjust.txt")
do zct=1,znum; do kct=1,knum
write(8,*) cdfxi(xistar(zct,kct),xibar)
end do; end do
close(8)

open(8,file="expeckprime.txt")
do zct=1,znum; do kct=1,knum
write(8,*) cdfxi(xistar(zct,kct),xibar)*kprime(zct,kct) &
    + (1.0-cdfxi(xistar(zct,kct),xibar))*(1.0-delta)*k0(kct)
end do; end do
close(8)

open(8,file="constants.txt")
do ct=1,numconstants
write(8,*) constantvec(ct)
end do 
close(8)

open(8,file="expeckprimedense.txt")
do zct=1,znum; do kct=1,kdensenum
write(8,*) cdfxi(xistardense(zct,kct),xibar)*kprimedense(zct,kct) &
    + (1.0-cdfxi(xistardense(zct,kct),xibar))*(1.0-delta)*kdense0(kct)
end do; end do
close(8)

open(8,file="xistardense.txt")
do zct=1,znum; do kct=1,kdensenum
write(8,*) xistardense(zct,kct)
end do; end do
close(8)

open(8,file="padjustdense.txt")
do zct=1,znum; do kct=1,kdensenum
write(8,*) cdfxi(xistardense(zct,kct),xibar)
end do; end do
close(8)

open(8,file="kprimedense.txt")
do zct=1,znum; do kct=1,kdensenum
write(8,*) kprimedense(zct,kct)
end do; end do
close(8)

open(8,file="kdense0.txt")
do kct=1,kdensenum
write(8,*) kdense0(kct)
end do
close(8)

open(8,file="distkz.txt")
do zct=1,znum; do kct=1,kdensenum
write(8,*) distkz(zct,kct)
end do; end do
close(8)

call cpu_time(finish)
write(13,"(A,F15.1,A)") "Finished in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished in ",finish-start," seconds."
close(13); !closing log file

contains


double precision function cdfxi(xi,xibar)
implicit none
double precision :: xi,xibar
if (xi<0.0) cdfxi = 0.0
if (xi>0.0.and.xi<=xibar) cdfxi = xi/xibar
if (xi>xibar) cdfxi = 1.0
end function cdfxi
    
double precision function expecxi(xi,xibar)
implicit none
double precision :: xi,xibar
if(xi<0.0) expecxi = 0.0
if(xi<=xibar.and.xi>=0.0)expecxi = ( xi ** 2.0 ) / (2.0*xibar)
if(xi>xibar)expecxi = ( xibar ** 2.0 ) / (2.0*xibar)
end function expecxi

double precision function freduced(z,k,a,w,alpha,nu)
implicit none
double precision :: z,k,a,w,alpha,nu

double precision :: exponentnu

exponentnu = 1.0 / (1.0 - nu)

freduced = ( nu ** (nu * exponentnu) ) * (1.0 - nu) * ( w ** ( -1.0 * nu * exponentnu ) )
freduced = freduced * (z ** exponentnu) * (a ** exponentnu) * ( k ** (alpha * exponentnu) )
end function freduced

double precision function yreduced(z,k,a,w,alpha,nu)
implicit none
double precision :: z,k,a,w,alpha,nu

double precision :: exponentnu

exponentnu = 1.0 / (1.0 - nu)

yreduced = ( (nu/w) ** (nu * exponentnu) ) * ( (a * z ) ** exponentnu  )
yreduced = yreduced * ( k ** ( alpha * exponentnu ) )

end function yreduced





double precision function fVa(kprime,param,nparam)
implicit none

integer :: nparam
double precision :: param(nparam),kprime

integer :: ct,zct,znum,knum,zprimect,kct,zdumct

double precision :: alpha,beta,delta,nu,zval,aval,kval,pval,phi,wval,Vnextval

double precision, allocatable :: Vold(:,:),V2old(:,:),pr_mat_z(:,:),k0(:)

!now, read in the stuff from the vector param

knum = int(param(nparam))
znum = int(param(nparam-1))

allocate(Vold(znum,knum),V2old(znum,knum),k0(knum),pr_mat_z(znum,znum))

ct=0
do zdumct=1,znum; !note that we use zdumct here rather than zct because zct (and act) are reserved
do kct=1,knum
    ct=ct+1
    Vold(zdumct,kct) = param(ct)
end do    
end do    

do zdumct=1,znum
do kct=1,knum
    ct=ct+1
    V2old(zdumct,kct) = param(ct)
end do    
end do    


do zdumct=1,znum
do zprimect=1,znum
    ct=ct+1
    pr_mat_z(zdumct,zprimect) = param(ct)
end do    
end do    

do kct=1,knum
    ct=ct+1
    k0(kct) = param(ct)
end do 

ct=ct+1; alpha = param(ct);
ct=ct+1; beta = param(ct);
ct=ct+1; delta = param(ct);
ct=ct+1; nu = param(ct);
ct=ct+1; phi = param(ct);
ct=ct+1; zval = param(ct);
ct=ct+1; aval = param(ct);
ct=ct+1; kval = param(ct);
ct=ct+1; pval = param(ct);
ct=ct+1; zct = int(param(ct));

!!now that parameters are read in
wval = phi/pval
fVa = freduced(zval,kval,aval,wval,alpha,nu) - kprime + (1.0 - delta) * kval
fVa = fVa * pval

do zprimect=1,znum
    call splint(k0,Vold(zprimect,:),V2old(zprimect,:),knum,kprime,Vnextval)

    !now that Vnextval is evaluated, add to expectation
    fVa = fVa + beta * pr_mat_z(zct,zprimect) * Vnextval
end do

fVa  = -1.0 * fVa; !take into account the minimization in the Brent routine
end function fVa

end program kt_ss