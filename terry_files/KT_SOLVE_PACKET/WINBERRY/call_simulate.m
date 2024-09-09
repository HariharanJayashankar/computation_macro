%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call_simulate.m
%
% MATLAB code to called by kt_reiter.f90 to 
% implement unconditional simulation for the
% the Winberry (2016) solution of the Khan and Thomas (2008)
% model.
%
% 'Alternative Methods for Solving Heterogeneous Firm Models'
% Stephen Terry (2017)
%
% This Version : 01/13/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%load in data from call_gensys, including model solution
load SIMDATAMAT.mat

%%%perform unconditional simulation of the model

%output is matrices for system X_t - X_SS= (Asol) (X_{t-1}-X_SS) + (Bsol) epsilon_t,
logasim = log(a0(asimpos));

Xsim = zeros(numX,numper);
epssim = zeros(numper,1);
for t=2:numper;
    if (mod(t,250)==1); disp(['Winberry Sim for Period t = ' num2str(t)]); end;
    Xmin1 = Xsim(:,t-1);
    epsval = (logasim(t) - Asol(aind,:)*Xmin1)/Bsol(aind);
    Xsim(:,t) = Asol * Xmin1 + Bsol*epsval;
    epssim(t) = epsval;
end;
Xsim = Xsim + repmat(XSS,1,numper);

psim = Xsim(pind,:)';
Ysim = Xsim(Yind,:)';
Isim = Xsim(Iind,:)';
Nsim = Xsim(Nind,:)';


dlmwrite('psim.txt',exp(psim));
dlmwrite('ysim.txt',exp(Ysim));
dlmwrite('isim.txt',exp(Isim));
dlmwrite('Nsim.txt',exp(Nsim));
dlmwrite('epssim.txt',epssim);

Xsimvec = reshape(Xsim,numel(Xsim),1);
dlmwrite('Xsim.txt',Xsimvec);

quit;