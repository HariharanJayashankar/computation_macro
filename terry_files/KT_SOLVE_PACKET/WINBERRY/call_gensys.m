%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call_gensys.m
%
% MATLAB code called by kt_reiter.f90 to 
% implement a call to the Sims gensys solver
% for linear rational expectations models in 
% the Winberry (2016) solution of the Khan and Thomas (2008)
% model.
%
% 'Alternative Methods for Solving Heterogeneous Firm Models'
% Stephen Terry (2017)
%
% This Version : 01/13/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;clc;

irfhorz= 20;

constants = importdata('constants.txt');
numX = constants(1);
numeta = constants(2);
numeps = constants(3);
znum = constants(4);
knum = constants(5);
anum = constants(6);
numper = constants(7);
numperIRF = constants(8);
numsimIRF = constants(9);
shockperIRF = constants(10);
sigmaa = constants(11);
numdiscard = constants(12);
momuse = constants(13);

%set up the indexing of the state
Vaind = 1:(znum*knum);
Vnaind = (znum*knum+1):(2*znum*knum);
momind = (Vnaind(end)+1):(Vnaind(end)+znum*momuse);
pind = momind(end)+1;
Yind = pind+1;
Iind = Yind+1;
Nind = Iind + 1;
aind = Nind+1;


XSS = importdata('XSS.txt');
a0 = importdata('a0.txt');
asimpos = importdata('asimpos.txt');
asimposIRFvec = importdata('asimposIRF.txt');

F1vec = importdata('F1.txt');
F2vec = importdata('F2.txt');
F3vec = importdata('F3.txt');
F4vec = importdata('F4.txt');

F1=zeros(numX,numX);
F2 = F1;
F3 = zeros(numX,numeta);
F4 = zeros(numX,numeps);

ct=0;
for ct1=1:numX;
    for ct2=1:numX;
        ct=ct+1;
        F1(ct1,ct2) = F1vec(ct);
        F2(ct1,ct2) = F2vec(ct);
    end;
end;

ct=0;
for ct1=1:numX;
    for ct2=1:numeta;
        ct=ct+1;
        F3(ct1,ct2) = F3vec(ct);
    end;
end;

ct=0;
for ct1=1:numX;
    for ct2=1:numeps;
        ct=ct+1;
        F4(ct1,ct2) = F4vec(ct);
    end;
end;

%put in gensys notation
g0 = -F1; g1 = F2; c = zeros(numX,1); 
psi = F4; pi = F3;
[G1,C,impact,fmat,fwt,ywt,gev,eu]=gensys(g0,g1,c,psi,pi);

%Model is now solved
%output is matrices for system X_t - X_SS= (Asol) (X_{t-1}-X_SS) + (Bsol) epsilon_t,
%where Asol is numX x numX and Bsol is numX x 1
Asol = G1; Bsol=impact;

%save output for later simulation
save SIMDATAMAT
quit;
