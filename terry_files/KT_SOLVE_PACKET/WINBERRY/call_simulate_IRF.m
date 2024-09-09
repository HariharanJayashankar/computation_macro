%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call_simulate_IRF.m
%
% MATLAB code to called by kt_reiter.f90 to 
% implement IRF simulation for the
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

%%%perform IRF simulation of the model

%convert vectorized exog IRF agg prod process to matrix form
%(I have verified this delivers the same as loop-based read-in)
asimposIRF = reshape(asimposIRFvec,[2 numperIRF numsimIRF]);
asimposIRF = permute(asimposIRF,[2 3 1]);

%recall that solution is matrices for system (X_t - X_SS)= (Asol) (X_{t-1}-X_SS) + (Bsol) epsilon_t,
logasimIRF = log(a0(asimposIRF));

XsimIRF = zeros(numX,numperIRF,numsimIRF,2);
epssimIRF = zeros(numperIRF,numsimIRF,2);
for simct=1:numsimIRF;
for shockct=1:2;
for t=2:numperIRF;
    Xmin1 = XsimIRF(:,t-1,simct,shockct);
    epsval = (logasimIRF(t,simct,shockct) - Asol(aind,:)*Xmin1)/Bsol(aind);
    XsimIRF(:,t,simct,shockct) = Asol * Xmin1 + Bsol*epsval;
    epssimIRF(t,simct,shockct) = epsval;
end;
end;
end;
XsimIRF = XsimIRF + repmat(XSS,[1 numperIRF numsimIRF 2]);

psimIRF = squeeze(XsimIRF(pind,:,:,:));
YsimIRF = squeeze(XsimIRF(Yind,:,:,:));
IsimIRF = squeeze(XsimIRF(Iind,:,:,:));
NsimIRF = squeeze(XsimIRF(Nind,:,:,:));

psimIRF = reshape(psimIRF,numel(psimIRF),1);
YsimIRF = reshape(YsimIRF,numel(YsimIRF),1);
IsimIRF = reshape(IsimIRF,numel(IsimIRF),1);
NsimIRF = reshape(NsimIRF,numel(NsimIRF),1);

dlmwrite('ysimIRF.txt',exp(YsimIRF))
dlmwrite('isimIRF.txt',exp(IsimIRF))
dlmwrite('psimIRF.txt',exp(psimIRF))
dlmwrite('NsimIRF.txt',exp(NsimIRF))
dlmwrite('epssimIRF.txt',epssimIRF)

quit;