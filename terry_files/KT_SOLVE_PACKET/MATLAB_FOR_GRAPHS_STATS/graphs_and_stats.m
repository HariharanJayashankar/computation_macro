%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphs_and_stats.m
%
% MATLAB code to read in simulated data and produce
% tables and figures.
%
% 'Alternative Methods for Solving Heterogeneous Firm Models'
% Stephen Terry (2017)
%
% This Version : 01/16/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

plotsamp = 1500:1550; %sample for plot of unconditional simulation
plotsamp_perturb = plotsamp+1;

tplot = 1300; %period for distributional plots

lwidnum = 1.33;
titlesizenum = 20;
labelsizenum = 16;

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
kdensenum = constants(13);
momuse = constants(14);
nummicro = constants(15);

%sample for analysis
samp = (numdiscard+1):(numper-1);
samp_perturb = samp+1;

IRFsamp = (shockperIRF-1):(shockperIRF+20);
IRFsamp_perturb = IRFsamp+1;

numX_REITER = 2*znum*knum+znum*kdensenum+5+znum;
distind_REITER = (2*znum*knum+1):(2*znum*knum+znum*kdensenum);
numX_WINBERRY = 3*znum*knum+znum*momuse+5;

a0 = importdata('a0.txt');
k0 = importdata('k0.txt');
kdense0 = importdata('kdense0.txt');
kbar0 = importdata('kbar0.txt');
z0 = importdata('z0.txt');
ergdistz = importdata('ergdistz.txt');

%%%%%%%%%%%%%%%%%%
%%%%%READ IN AND PLOT UNCONDITIONAL SIMULATION DATA
%%%%%%%%%%%%%%%%%%

%exog shock a
asimpos = importdata('asimpos.txt');
asim = a0(asimpos);

%Y
Ysim_KS = importdata('./KS/ysim.txt');
Ysim_XPA = importdata('./XPA/ysim.txt');
Ysim_PARAM = importdata('./PARAM/ysim.txt');
Ysim_REITER = importdata('./REITER/ysim.txt');
Ysim_WINBERRY = importdata('./WINBERRY/ysim.txt');

%I
Isim_KS = importdata('./KS/isim.txt');
Isim_XPA = importdata('./XPA/isim.txt');
Isim_PARAM = importdata('./PARAM/isim.txt');
Isim_REITER = importdata('./REITER/isim.txt');
Isim_WINBERRY = importdata('./WINBERRY/isim.txt');

%N
Nsim_KS = importdata('./KS/Nsim.txt');
Nsim_XPA = importdata('./XPA/Nsim.txt');
Nsim_PARAM = importdata('./PARAM/nsim.txt');
Nsim_REITER = importdata('./REITER/Nsim.txt');
Nsim_WINBERRY = importdata('./WINBERRY/Nsim.txt');

%p
psim_KS = importdata('./KS/psim.txt');
psim_XPA = importdata('./XPA/psim.txt');
psim_PARAM = importdata('./PARAM/psim.txt');
psim_REITER = importdata('./REITER/psim.txt');
psim_WINBERRY = importdata('./WINBERRY/psim.txt');

%make unconditional simulation figure
figure;
subplot(2,2,1);
plot(plotsamp,log(Ysim_REITER(plotsamp_perturb)),'b',...
    plotsamp,log(Ysim_WINBERRY(plotsamp_perturb)),'m',...
    plotsamp,log(Ysim_PARAM(plotsamp)),'g',...
    plotsamp,log(Ysim_XPA(plotsamp)),'r',...
    plotsamp,log(Ysim_KS(plotsamp)),'k',...
        'LineWidth',lwidnum); 
title('Output','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
axis([-Inf Inf -0.75 -0.1])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;
legend('REITER','WINBERRY','PARAM','XPA','KS','Location','NorthEast');
legend boxoff;

subplot(2,2,2);
plot(plotsamp,log(Isim_REITER(plotsamp_perturb)),'b',...
    plotsamp,log(Isim_WINBERRY(plotsamp_perturb)),'m',...
    plotsamp,log(Isim_PARAM(plotsamp)),'g',...
    plotsamp,log(Isim_XPA(plotsamp)),'r',...
    plotsamp,log(Isim_KS(plotsamp)),'k',...
    'LineWidth',lwidnum);
title('Investment','FontSize',titlesizenum)
axis([-Inf Inf -3.0 -1.6])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(2,2,3);
plot(plotsamp,log(Nsim_REITER(plotsamp_perturb)),'b',...
    plotsamp,log(Nsim_WINBERRY(plotsamp_perturb)),'m',...
    plotsamp,log(Nsim_PARAM(plotsamp)),'g',...
    plotsamp,log(Nsim_XPA(plotsamp)),'r',...
    plotsamp,log(Nsim_KS(plotsamp)),'k',...
    'LineWidth',lwidnum);
title('Labor','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -1.225 -0.975])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(2,2,4);
plot(plotsamp,-log(psim_REITER(plotsamp_perturb)),'b',...
    plotsamp,-log(psim_WINBERRY(plotsamp_perturb)),'m',...
    plotsamp,-log(psim_PARAM(plotsamp)),'g',...
    plotsamp,-log(psim_XPA(plotsamp)),'r',...
    plotsamp,-log(psim_KS(plotsamp)),'k',...
    'LineWidth',lwidnum);
title('Consumption','FontSize',titlesizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -0.9 -0.7])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

saveas(gcf,'UNCOND_2BY2.pdf')


%%%%%%%%%%%%%%%%%%
%%%%%READ IN AND PLOT DISTRIBUTIONAL DATA
%%%%%%%%%%%%%%%%%%

%note: each of the methods stores the info slightly differently, so
%processing isn't simultaneous

%KS, XPA
distkz_KSvec = importdata('./KS/distkzsim.txt');
distkz_XPAvec = importdata('./XPA/distkzsim.txt');

distkz_KS=zeros(znum,kdensenum,numper);
distkz_XPA=zeros(znum,kdensenum,numper);
ct=0;
for zct=1:znum; for kct=1:kdensenum; for t=1:numper;
     ct=ct+1;
     distkz_KS(zct,kct,t)=distkz_KSvec(ct);
     distkz_XPA(zct,kct,t)=distkz_XPAvec(ct);
end;end;end;

%PARAM
distkz_PARAMvec = importdata('./PARAM/distkzsim.txt');
distkz_PARAM=zeros(znum,kdensenum,numper);
ct=0;
for t=1:(numper-1);
    for zct=1:znum;for kct=1:kdensenum;
            ct=ct+1;
        distkz_PARAM(zct,kct,t) = distkz_PARAMvec(ct);
    end;end;

    for zct=1:znum;
        distkz_PARAM(zct,:,t) = distkz_PARAM(zct,:,t)/sum(distkz_PARAM(zct,:,t));
        distkz_PARAM(zct,:,t) = distkz_PARAM(zct,:,t)*ergdistz(zct);
    end;
end;

%REITER
Xsim_REITERvec = importdata('./REITER/Xsim.txt');
Xsim_REITER = reshape(Xsim_REITERvec,[numX_REITER numper]);

distkz_REITERvec = Xsim_REITER(distind_REITER,:);
distkz_REITER = zeros(znum,kdensenum,numper);
for t=1:numper;
ct=0;
for zct=1:znum;
    for kct=1:kdensenum;
    ct=ct+1;        
    distkz_REITER(zct,kct,t) = distkz_REITERvec(ct,t);
    end;
end;
end;

%WINBERRY
distkz_WINBERRYvec = importdata('./WINBERRY/distkzsim.txt');
distkz_WINBERRY=zeros(znum,kdensenum,numper);
ct=0;
for t=1:(numper-1);
    for zct=1:znum;for kct=1:kdensenum;
            ct=ct+1;
        distkz_WINBERRY(zct,kct,t) = distkz_WINBERRYvec(ct);
    end;end;

    for zct=1:znum;
        distkz_WINBERRY(zct,:,t) = distkz_WINBERRY(zct,:,t)/sum(distkz_WINBERRY(zct,:,t));
        distkz_WINBERRY(zct,:,t) = distkz_WINBERRY(zct,:,t)*ergdistz(zct);
    end;
    
end;

%actually create distributional simulation figure
distKS_plot = squeeze(distkz_KS(:,:,tplot));
distXPA_plot = squeeze(distkz_XPA(:,:,tplot));
distPARAM_plot = squeeze(distkz_PARAM(:,:,tplot));
distREITER_plot = squeeze(distkz_REITER(:,:,tplot));
distWINBERRY_plot = squeeze(distkz_WINBERRY(:,:,tplot));


figure;
subplot(3,2,1);
plot(kdense0,distKS_plot(1,:),'k',...
    kdense0,distKS_plot(2,:),'r',...
    kdense0,distKS_plot(3,:),'g',...
    kdense0,distKS_plot(4,:),'b',...
    kdense0,distKS_plot(5,:),'m','LineWidth',lwidnum);
title('KS','FontSize',titlesizenum)
ylabel('Density','FontSize',labelsizenum)
ax=gca;
ax.FontSize = labelsizenum;
axis([0.2 3.5 0 0.075])

subplot(3,2,3);
plot(kdense0,distXPA_plot(1,:),'k',...
    kdense0,distXPA_plot(2,:),'r',...
    kdense0,distXPA_plot(3,:),'g',...
    kdense0,distXPA_plot(4,:),'b',...
    kdense0,distXPA_plot(5,:),'m','LineWidth',lwidnum);
title('XPA','FontSize',titlesizenum)
ylabel('Density','FontSize',labelsizenum)
ax=gca;
ax.FontSize = labelsizenum;
axis([0.2 3.5 0 0.075])

subplot(3,2,5);
plot(kdense0,distREITER_plot(1,:),'k',...
    kdense0,distREITER_plot(2,:),'r',...
    kdense0,distREITER_plot(3,:),'g',...
    kdense0,distREITER_plot(4,:),'b',...
    kdense0,distREITER_plot(5,:),'m','LineWidth',lwidnum);
title('REITER','FontSize',titlesizenum)
xlabel('Capital','FontSize',labelsizenum)
ylabel('Density','FontSize',labelsizenum)
ax=gca;
ax.FontSize = labelsizenum;
axis([0.2 3.5 0 0.075])

subplot(3,2,2);
plot(kdense0,distPARAM_plot(1,:),'k',...
    kdense0,distPARAM_plot(2,:),'r',...
    kdense0,distPARAM_plot(3,:),'g',...
    kdense0,distPARAM_plot(4,:),'b',...
    kdense0,distPARAM_plot(5,:),'m','LineWidth',lwidnum);
title('PARAM','FontSize',titlesizenum)
ax=gca;
ax.FontSize = labelsizenum;
axis([0.2 3.5 0 0.075])


subplot(3,2,4);
plot(kdense0,distWINBERRY_plot(1,:),'k',...
    kdense0,distWINBERRY_plot(2,:),'r',...
    kdense0,distWINBERRY_plot(3,:),'g',...
    kdense0,distWINBERRY_plot(4,:),'b',...
    kdense0,distWINBERRY_plot(5,:),'m','LineWidth',lwidnum);
title('WINBERRY','FontSize',titlesizenum)
xlabel('Capital','FontSize',labelsizenum)
ax=gca;
ax.FontSize = labelsizenum;
axis([0.2 3.5 0 0.075])

saveas(gcf,'DIST_3BY2.pdf')


%%%%%%%%%%%%%%%%%%
%%%%%READ IN MICRO MOMENTS FROM UNCONDTIONAL SIMULATION
%%%%%%%%%%%%%%%%%%

%read in raw data on micro moments over uncondtional sim
MICROsim_KS = importdata('./KS/MICROsim.txt');
MICROsim_XPA = importdata('./XPA/MICROsim.txt');
MICROsim_PARAM = importdata('./PARAM/MICROsim.txt');
MICROsim_REITER = importdata('./REITER/MICROsim.txt');
MICROsim_WINBERRY = importdata('./WINBERRY/MICROsim.txt');

MICROsim_KS = reshape(MICROsim_KS,nummicro,numper);
MICROsim_XPA = reshape(MICROsim_XPA,nummicro,numper);
MICROsim_PARAM = reshape(MICROsim_PARAM,nummicro,numper);
MICROsim_REITER = reshape(MICROsim_REITER,nummicro,numper);
MICROsim_WINBERRY = reshape(MICROsim_WINBERRY,nummicro,numper);

%compute mean statistics
MICRO_TABLE = zeros(nummicro,5);
MICRO_TABLE(:,1) = mean(MICROsim_KS(:,samp),2);
MICRO_TABLE(:,2) = mean(MICROsim_XPA(:,samp),2);
MICRO_TABLE(:,3) = mean(MICROsim_PARAM(:,samp),2);
MICRO_TABLE(:,4) = mean(MICROsim_REITER(:,samp),2);
MICRO_TABLE(:,5) = mean(MICROsim_WINBERRY(:,samp),2);

%%%%%%%%%%%%%%%%%%
%%%%%COMPUTE HP-FILTERED BC STATS
%%%%%%%%%%%%%%%%%%

HPdata_KS = log([Ysim_KS(samp) Isim_KS(samp) Nsim_KS(samp) asim(samp) psim_KS(samp)]);
HPdata_XPA = log([Ysim_XPA(samp) Isim_XPA(samp) Nsim_XPA(samp) asim(samp) psim_XPA(samp)]);
HPdata_PARAM = log([Ysim_PARAM(samp) Isim_PARAM(samp) Nsim_PARAM(samp) asim(samp) psim_PARAM(samp)]);
HPdata_REITER = log([Ysim_REITER(samp_perturb) Isim_REITER(samp_perturb) Nsim_REITER(samp_perturb) ...
    asim(samp) psim_REITER(samp_perturb)]);
HPdata_WINBERRY = log([Ysim_WINBERRY(samp_perturb) Isim_WINBERRY(samp_perturb) Nsim_WINBERRY(samp_perturb) asim(samp) psim_WINBERRY(samp_perturb)]);

[HPtrend_KS,HPdetrend_KS] = hptrend(HPdata_KS,100);
[HPtrend_XPA,HPdetrend_XPA] = hptrend(HPdata_XPA,100);
[HPtrend_PARAM,HPdetrend_PARAM] = hptrend(HPdata_PARAM,100);
[HPtrend_REITER,HPdetrend_REITER] = hptrend(HPdata_REITER,100);
[HPtrend_WINBERRY,HPdetrend_WINBERRY] = hptrend(HPdata_WINBERRY,100);

HP_VOL_TABLE = zeros(5,5);
HP_VOL_TABLE(1,:) = 100*var(HPdetrend_KS).^0.5; HP_VOL_TABLE(1,2:end)=HP_VOL_TABLE(1,2:end)./HP_VOL_TABLE(1,1);
HP_VOL_TABLE(2,:) = 100*var(HPdetrend_XPA).^0.5; HP_VOL_TABLE(2,2:end)=HP_VOL_TABLE(2,2:end)./HP_VOL_TABLE(2,1);
HP_VOL_TABLE(3,:) = 100*var(HPdetrend_PARAM).^0.5; HP_VOL_TABLE(3,2:end)=HP_VOL_TABLE(3,2:end)./HP_VOL_TABLE(3,1);
HP_VOL_TABLE(4,:) = 100*var(HPdetrend_REITER).^0.5; HP_VOL_TABLE(4,2:end)=HP_VOL_TABLE(4,2:end)./HP_VOL_TABLE(4,1);
HP_VOL_TABLE(5,:) = 100*var(HPdetrend_WINBERRY).^0.5; HP_VOL_TABLE(5,2:end)=HP_VOL_TABLE(5,2:end)./HP_VOL_TABLE(5,1);

HP_CORR_TABLE = zeros(5,5);

corrmat = corrcoef(HPdetrend_KS); HP_CORR_TABLE(1,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_XPA); HP_CORR_TABLE(2,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_PARAM); HP_CORR_TABLE(3,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_REITER); HP_CORR_TABLE(4,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_WINBERRY); HP_CORR_TABLE(5,:) = corrmat(1,:);

%%%%%%%%%%%%%%%%%%
%%%%%READ IN AND COMPUTE IRFs
%%%%%%%%%%%%%%%%%%

%read in exogenous process
asimposIRFvec = importdata('asimposIRF.txt');
asimposIRF = zeros(numperIRF,numsimIRF,2);

ct=0;
for simct=1:numsimIRF;
for t=1:numperIRF;
for shockct=1:2;
    ct = ct+1;
    asimposIRF(t,simct,shockct) = asimposIRFvec(ct);
end;
end;
end;

asimIRF = a0(asimposIRF);

%Read in simulated series for each method. Different storage conventions
%imply that this isn't entirely joint

%KS
YsimIRF_KSvec = importdata('./KS/ysimIRF.txt');
IsimIRF_KSvec = importdata('./KS/isimIRF.txt');
NsimIRF_KSvec = importdata('./KS/NsimIRF.txt');
psimIRF_KSvec = importdata('./KS/psimIRF.txt');

YsimIRF_KS = zeros(numperIRF,numsimIRF,2);
IsimIRF_KS = zeros(numperIRF,numsimIRF,2);
NsimIRF_KS = zeros(numperIRF,numsimIRF,2);
psimIRF_KS = zeros(numperIRF,numsimIRF,2);

ct=0;
for simct=1:numsimIRF;
for t=1:numperIRF;
for shockct=1:2;
ct=ct+1;

YsimIRF_KS(t,simct,shockct) = YsimIRF_KSvec(ct);    
IsimIRF_KS(t,simct,shockct) = IsimIRF_KSvec(ct);
NsimIRF_KS(t,simct,shockct) = NsimIRF_KSvec(ct);
psimIRF_KS(t,simct,shockct) = psimIRF_KSvec(ct);

end;
end;
end;


%XPA
YsimIRF_XPAvec = importdata('./XPA/ysimIRF.txt');
IsimIRF_XPAvec = importdata('./XPA/isimIRF.txt');
NsimIRF_XPAvec = importdata('./XPA/NsimIRF.txt');
psimIRF_XPAvec = importdata('./XPA/psimIRF.txt');

YsimIRF_XPA = zeros(numperIRF,numsimIRF,2);
IsimIRF_XPA = zeros(numperIRF,numsimIRF,2);
NsimIRF_XPA = zeros(numperIRF,numsimIRF,2);
psimIRF_XPA = zeros(numperIRF,numsimIRF,2);

ct=0;
for simct=1:numsimIRF;
for t=1:numperIRF;
for shockct=1:2;
ct=ct+1;

YsimIRF_XPA(t,simct,shockct) = YsimIRF_XPAvec(ct);    
IsimIRF_XPA(t,simct,shockct) = IsimIRF_XPAvec(ct);
NsimIRF_XPA(t,simct,shockct) = NsimIRF_XPAvec(ct);
psimIRF_XPA(t,simct,shockct) = psimIRF_XPAvec(ct);

end;
end;
end;

%PARAM
YsimIRF_PARAMvec = importdata('./PARAM/ysimIRF.txt');
IsimIRF_PARAMvec = importdata('./PARAM/isimIRF.txt');
NsimIRF_PARAMvec = importdata('./PARAM/nsimIRF.txt');
psimIRF_PARAMvec = importdata('./PARAM/psimIRF.txt');

YsimIRF_PARAM = zeros(numperIRF,numsimIRF,2);
IsimIRF_PARAM = zeros(numperIRF,numsimIRF,2);
NsimIRF_PARAM = zeros(numperIRF,numsimIRF,2);
psimIRF_PARAM = zeros(numperIRF,numsimIRF,2);

ct=0;
for shockct=1:2;
for simct=1:numsimIRF;
for t=1:numperIRF;
ct=ct+1;

YsimIRF_PARAM(t,simct,shockct) = YsimIRF_PARAMvec(ct);    
IsimIRF_PARAM(t,simct,shockct) = IsimIRF_PARAMvec(ct);
NsimIRF_PARAM(t,simct,shockct) = NsimIRF_PARAMvec(ct);
psimIRF_PARAM(t,simct,shockct) = psimIRF_PARAMvec(ct);

end;
end;
end;


%REITER
YsimIRF_REITERvec = importdata('./REITER/ysimIRF.txt');
IsimIRF_REITERvec = importdata('./REITER/isimIRF.txt');
NsimIRF_REITERvec = importdata('./REITER/NsimIRF.txt');
psimIRF_REITERvec = importdata('./REITER/psimIRF.txt');

YsimIRF_REITER = zeros(numperIRF,numsimIRF,2);
IsimIRF_REITER = zeros(numperIRF,numsimIRF,2);
NsimIRF_REITER = zeros(numperIRF,numsimIRF,2);
psimIRF_REITER = zeros(numperIRF,numsimIRF,2);

ct=0;
for shockct=1:2;
for simct=1:numsimIRF;
for t=1:numperIRF;
ct=ct+1;

YsimIRF_REITER(t,simct,shockct) = YsimIRF_REITERvec(ct);    
IsimIRF_REITER(t,simct,shockct) = IsimIRF_REITERvec(ct);
NsimIRF_REITER(t,simct,shockct) = NsimIRF_REITERvec(ct);
psimIRF_REITER(t,simct,shockct) = psimIRF_REITERvec(ct);

end;
end;
end;

%WINBERRY
YsimIRF_WINBERRYvec = importdata('./WINBERRY/ysimIRF.txt');
IsimIRF_WINBERRYvec = importdata('./WINBERRY/isimIRF.txt');
NsimIRF_WINBERRYvec = importdata('./WINBERRY/NsimIRF.txt');
psimIRF_WINBERRYvec = importdata('./WINBERRY/psimIRF.txt');

YsimIRF_WINBERRY = zeros(numperIRF,numsimIRF,2);
IsimIRF_WINBERRY = zeros(numperIRF,numsimIRF,2);
NsimIRF_WINBERRY = zeros(numperIRF,numsimIRF,2);
psimIRF_WINBERRY = zeros(numperIRF,numsimIRF,2);

ct=0;
for shockct=1:2;
for simct=1:numsimIRF;
for t=1:numperIRF;
ct=ct+1;

YsimIRF_WINBERRY(t,simct,shockct) = YsimIRF_WINBERRYvec(ct);    
IsimIRF_WINBERRY(t,simct,shockct) = IsimIRF_WINBERRYvec(ct);
NsimIRF_WINBERRY(t,simct,shockct) = NsimIRF_WINBERRYvec(ct);
psimIRF_WINBERRY(t,simct,shockct) = psimIRF_WINBERRYvec(ct);

end;
end;
end;

%now, process and compute IRFs

%exog process first
aIRF = 100*(mean(log(asimIRF(:,:,2)),2) - mean(log(asimIRF(:,:,1)),2));

%output
YIRF_KS = 100*(mean(log(YsimIRF_KS(:,:,2)),2) - mean(log(YsimIRF_KS(:,:,1)),2));
YIRF_XPA = 100*(mean(log(YsimIRF_XPA(:,:,2)),2) - mean(log(YsimIRF_XPA(:,:,1)),2));
YIRF_PARAM = 100*(mean(log(YsimIRF_PARAM(:,:,2)),2) - mean(log(YsimIRF_PARAM(:,:,1)),2));
YIRF_REITER = 100*(mean(log(YsimIRF_REITER(:,:,2)),2) - mean(log(YsimIRF_REITER(:,:,1)),2));
YIRF_WINBERRY = 100*(mean(log(YsimIRF_WINBERRY(:,:,2)),2) - mean(log(YsimIRF_WINBERRY(:,:,1)),2));

%investment
IIRF_KS = 100*(mean(log(IsimIRF_KS(:,:,2)),2) - mean(log(IsimIRF_KS(:,:,1)),2));
IIRF_XPA = 100*(mean(log(IsimIRF_XPA(:,:,2)),2) - mean(log(IsimIRF_XPA(:,:,1)),2));
IIRF_PARAM = 100*(mean(log(IsimIRF_PARAM(:,:,2)),2) - mean(log(IsimIRF_PARAM(:,:,1)),2));
IIRF_REITER = 100*(mean(log(IsimIRF_REITER(:,:,2)),2) - mean(log(IsimIRF_REITER(:,:,1)),2));
IIRF_WINBERRY = 100*(mean(log(IsimIRF_WINBERRY(:,:,2)),2) - mean(log(IsimIRF_WINBERRY(:,:,1)),2));

%labor
NIRF_KS = 100*(mean(log(NsimIRF_KS(:,:,2)),2) - mean(log(NsimIRF_KS(:,:,1)),2));
NIRF_XPA = 100*(mean(log(NsimIRF_XPA(:,:,2)),2) - mean(log(NsimIRF_XPA(:,:,1)),2));
NIRF_PARAM = 100*(mean(log(NsimIRF_PARAM(:,:,2)),2) - mean(log(NsimIRF_PARAM(:,:,1)),2));
NIRF_REITER = 100*(mean(log(NsimIRF_REITER(:,:,2)),2) - mean(log(NsimIRF_REITER(:,:,1)),2));
NIRF_WINBERRY = 100*(mean(log(NsimIRF_WINBERRY(:,:,2)),2) - mean(log(NsimIRF_WINBERRY(:,:,1)),2));

%price
pIRF_KS = 100*(mean(log(psimIRF_KS(:,:,2)),2) - mean(log(psimIRF_KS(:,:,1)),2));
pIRF_XPA = 100*(mean(log(psimIRF_XPA(:,:,2)),2) - mean(log(psimIRF_XPA(:,:,1)),2));
pIRF_PARAM = 100*(mean(log(psimIRF_PARAM(:,:,2)),2) - mean(log(psimIRF_PARAM(:,:,1)),2));
pIRF_REITER = 100*(mean(log(psimIRF_REITER(:,:,2)),2) - mean(log(psimIRF_REITER(:,:,1)),2));
pIRF_WINBERRY = 100*(mean(log(psimIRF_WINBERRY(:,:,2)),2) - mean(log(psimIRF_WINBERRY(:,:,1)),2));


%plot IRF figure
IRFticksamp = IRFsamp - (shockperIRF);


figure;
subplot(2,2,1);
plot(plotsamp,log(asim(plotsamp)),'k',...
        'LineWidth',lwidnum); 
title('Unconditional Simulation','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
axis([-Inf Inf -0.075 0.075])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;
xlabel('Year','FontSize',labelsizenum)

subplot(2,2,2);
plot(IRFsamp,aIRF(IRFsamp),'k',...
    IRFsamp,zeros(size(IRFsamp)),'k',...
        'LineWidth',lwidnum); 
title('Impulse Response','FontSize',titlesizenum)
ylabel('Percent','FontSize',labelsizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -1 2])
ax=gca;
ax.XTick = [IRFsamp(2) IRFsamp(7) IRFsamp(12) IRFsamp(17) IRFsamp(22)];
ax.XTickLabel = {IRFticksamp(2) IRFticksamp(7) IRFticksamp(12) IRFticksamp(17) IRFticksamp(22) };
ax.FontSize = labelsizenum;
saveas(gcf,'EXOG_PROD.pdf')

figure;
subplot(2,2,1);
plot(IRFsamp,YIRF_REITER(IRFsamp_perturb),'b',...
    IRFsamp,YIRF_WINBERRY(IRFsamp_perturb),'m',...
    IRFsamp,YIRF_PARAM(IRFsamp),'g',...
    IRFsamp,YIRF_XPA(IRFsamp),'r',...
    IRFsamp,YIRF_KS(IRFsamp),'k',...
    IRFsamp,zeros(size(IRFsamp)),'k',...
        'LineWidth',lwidnum); 
title('Output','FontSize',titlesizenum)
ylabel('Percent','FontSize',labelsizenum)
axis([-Inf Inf -1 3])
ax=gca;
ax.XTick = [IRFsamp(2) IRFsamp(7) IRFsamp(12) IRFsamp(17) IRFsamp(22)];
ax.XTickLabel = {IRFticksamp(2) IRFticksamp(7) IRFticksamp(12) IRFticksamp(17) IRFticksamp(22) };
ax.FontSize = labelsizenum;


subplot(2,2,2);
plot(IRFsamp,IIRF_REITER(IRFsamp_perturb),'b',...
    IRFsamp,IIRF_WINBERRY(IRFsamp_perturb),'m',...
    IRFsamp,IIRF_PARAM(IRFsamp),'g',...
    IRFsamp,IIRF_XPA(IRFsamp),'r',...
    IRFsamp,IIRF_KS(IRFsamp),'k',...
    IRFsamp,zeros(size(IRFsamp)),'k',...
        'LineWidth',lwidnum); 
title('Investment','FontSize',titlesizenum)
axis([-Inf Inf -5 17.5])
ax=gca;
ax.XTick = [IRFsamp(2) IRFsamp(7) IRFsamp(12) IRFsamp(17) IRFsamp(22)];
ax.XTickLabel = {IRFticksamp(2) IRFticksamp(7) IRFticksamp(12) IRFticksamp(17) IRFticksamp(22) };
ax.FontSize = labelsizenum;
legend('REITER','WINBERRY','PARAM','XPA','KS','Location','NorthEast');
legend boxoff;

subplot(2,2,3);
plot(IRFsamp,NIRF_REITER(IRFsamp_perturb),'b',...
    IRFsamp,NIRF_WINBERRY(IRFsamp_perturb),'m',...
    IRFsamp,NIRF_PARAM(IRFsamp),'g',...
    IRFsamp,NIRF_XPA(IRFsamp),'r',...
    IRFsamp,NIRF_KS(IRFsamp),'k',...
    IRFsamp,zeros(size(IRFsamp)),'k',...
        'LineWidth',lwidnum); 
title('Labor','FontSize',titlesizenum)
ylabel('Percent','FontSize',labelsizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -1 3])
ax=gca;
ax.XTick = [IRFsamp(2) IRFsamp(7) IRFsamp(12) IRFsamp(17) IRFsamp(22)];
ax.XTickLabel = {IRFticksamp(2) IRFticksamp(7) IRFticksamp(12) IRFticksamp(17) IRFticksamp(22) };
ax.FontSize = labelsizenum;

subplot(2,2,4);
plot(IRFsamp,-pIRF_REITER(IRFsamp_perturb),'b',...
    IRFsamp,-pIRF_WINBERRY(IRFsamp_perturb),'m',...
    IRFsamp,-pIRF_PARAM(IRFsamp),'g',...
    IRFsamp,-pIRF_XPA(IRFsamp),'r',...
    IRFsamp,-pIRF_KS(IRFsamp),'k',...
    IRFsamp,zeros(size(IRFsamp)),'k',...
        'LineWidth',lwidnum); 
title('Consumption','FontSize',titlesizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -1 3])
ax=gca;
ax.XTick = [IRFsamp(2) IRFsamp(7) IRFsamp(12) IRFsamp(17) IRFsamp(22)];
ax.XTickLabel = {IRFticksamp(2) IRFticksamp(7) IRFticksamp(12) IRFticksamp(17) IRFticksamp(22) };
ax.FontSize = labelsizenum;
saveas(gcf,'IRF_2BY2.pdf')

%%%%%%%%%%%%%%%%%%
%%%%%COMPUTE ACCURACY STATS AND PRODUCE FCST FIGURES ***FROM SOLUTION
%%%%%SIMULATION***
%%%%%%%%%%%%%%%%%%

%KS
Kbarsim_KS = importdata('./KS/kbarsim.txt');
KbarFCST_KS = importdata('./KS/kbarfcstsim.txt');
KbarDH_KS = importdata('./KS/kbaraggrule.txt');
pFCST_KS = importdata('./KS/pfcstsim.txt');
pDH_KS = importdata('./KS/paggrule.txt');

%XPA
Kbarsim_XPA = importdata('./XPA/kbarsim.txt');
KbarFCST_XPA = importdata('./XPA/kbarfcstsim.txt');
KbarDH_XPA = importdata('./XPA/kbaraggrule.txt');
pFCST_XPA = importdata('./XPA/pfcstsim.txt');
pDH_XPA = importdata('./XPA/paggrule.txt');

%PARAM
Kbarsim_PARAM = importdata('./PARAM/kbarsim.txt');
KbarFCST_PARAM = importdata('./PARAM/kbarfcstsim.txt');
KbarDH_PARAM = importdata('./PARAM/kbarDHsim.txt');
pFCST_PARAM = importdata('./PARAM/pfcstsim.txt');
pDH_PARAM = importdata('./PARAM/pDHsim.txt');


%%COMPUTE DH STATS
DH_TABLE = zeros(3,4);

%KS
DH_TABLE(1,1) = max(100*abs(log(Kbarsim_KS(samp))-log(KbarDH_KS(samp))));
DH_TABLE(1,2) = mean(100*abs(log(Kbarsim_KS(samp))-log(KbarDH_KS(samp))));
DH_TABLE(1,3) = max(100*abs(log(psim_KS(samp))-log(pDH_KS(samp))));
DH_TABLE(1,4) = mean(100*abs(log(psim_KS(samp))-log(pDH_KS(samp))));

%XPA
DH_TABLE(2,1) = max(100*abs(log(Kbarsim_XPA(samp))-log(KbarDH_XPA(samp))));
DH_TABLE(2,2) = mean(100*abs(log(Kbarsim_XPA(samp))-log(KbarDH_XPA(samp))));
DH_TABLE(2,3) = max(100*abs(log(psim_XPA(samp))-log(pDH_XPA(samp))));
DH_TABLE(2,4) = mean(100*abs(log(psim_XPA(samp))-log(pDH_XPA(samp))));

%PARAM
DH_TABLE(3,1) = max(100*abs(log(Kbarsim_PARAM(samp))-log(KbarDH_PARAM(samp))));
DH_TABLE(3,2) = mean(100*abs(log(Kbarsim_PARAM(samp))-log(KbarDH_PARAM(samp))));
DH_TABLE(3,3) = max(100*abs(log(psim_PARAM(samp))-log(pDH_PARAM(samp))));
DH_TABLE(3,4) = mean(100*abs(log(psim_PARAM(samp))-log(pDH_PARAM(samp))));


%COMPUTE R^2s & RMSEs
ESS_KS = zeros(anum,2);
ESS_PARAM = zeros(anum,2);
ESS_XPA = zeros(anum,2);

TSS_KS = zeros(anum,2);
TSS_PARAM = zeros(anum,2);
TSS_XPA = zeros(anum,2);

MEANS_KS = zeros(anum,2);
MEANS_PARAM = zeros(anum,2);
MEANS_XPA = zeros(anum,2);

COUNTS = zeros(anum,1);

%first, compute means
for t = samp(1):samp(end);
    act = asimpos(t);
    COUNTS(act) = COUNTS(act) + 1;
    
    MEANS_KS(act,:) = MEANS_KS(act,:) + [log(psim_KS(t)) log(Kbarsim_KS(t+1))];
    MEANS_PARAM(act,:) = MEANS_PARAM(act,:) + [log(psim_PARAM(t)) log(Kbarsim_PARAM(t+1))];
    MEANS_XPA(act,:) = MEANS_XPA(act,:) + [log(psim_XPA(t)) log(Kbarsim_XPA(t+1))];
    
end;
MEANS_KS = MEANS_KS./repmat(COUNTS,1,2);
MEANS_PARAM = MEANS_PARAM./repmat(COUNTS,1,2);
MEANS_XPA = MEANS_XPA./repmat(COUNTS,1,2);

%then, compute ESS and TSS
for t = samp(1):samp(end);
    act = asimpos(t);
   
    %ESS
    ESS_KS(act,:) = ESS_KS(act,:) + [log(pFCST_KS(t)/psim_KS(t)) log(KbarFCST_KS(t)/Kbarsim_KS(t))].^2;
    ESS_PARAM(act,:) = ESS_PARAM(act,:) + [log(pFCST_PARAM(t)/psim_PARAM(t)) log(KbarFCST_PARAM(t)/Kbarsim_PARAM(t))].^2;
    ESS_XPA(act,:) = ESS_XPA(act,:) + [log(pFCST_XPA(t)/psim_XPA(t)) log(KbarFCST_XPA(t)/Kbarsim_XPA(t))].^2;
    
    %TSS
    TSS_KS(act,:) = TSS_KS(act,:) + [log(MEANS_KS(act,1)/psim_KS(t)) log(MEANS_KS(act,2)/Kbarsim_KS(t))].^2;
    TSS_PARAM(act,:) = TSS_PARAM(act,:) + [log(MEANS_PARAM(act,1)/psim_PARAM(t)) log(MEANS_PARAM(act,2)/Kbarsim_PARAM(t))].^2;
    TSS_XPA(act,:) = TSS_XPA(act,:) + [log(MEANS_XPA(act,1)/psim_XPA(t)) log(MEANS_XPA(act,2)/Kbarsim_XPA(t))].^2;
end;

%now, actually implement the R^2 and RMSE formulas
R2_KS = 1 - ESS_KS./TSS_KS;
R2_PARAM = 1 - ESS_PARAM./TSS_PARAM;
R2_XPA = 1 - ESS_XPA./TSS_XPA;

RMSE_KS = 100*sqrt(ESS_KS./repmat(COUNTS,1,2));
RMSE_PARAM = 100*sqrt(ESS_PARAM./repmat(COUNTS,1,2));
RMSE_XPA = 100*sqrt(ESS_XPA./repmat(COUNTS,1,2));


%%%%%%%%%%%%%%%%%%
%%%%%COMPUTE MEAN DIFFS FOR IRFS & UNCOND SIMS
%%%%%%%%%%%%%%%%%%

UNCOND_DIFFS_TABLE = zeros(4,4);

data_KS = log([Ysim_KS(samp) Isim_KS(samp) Nsim_KS(samp) psim_KS(samp)]);
data_XPA = log([Ysim_XPA(samp) Isim_XPA(samp) Nsim_XPA(samp) psim_XPA(samp)]);
data_PARAM = log([Ysim_PARAM(samp) Isim_PARAM(samp) Nsim_PARAM(samp) psim_PARAM(samp)]);
data_REITER = log([Ysim_REITER(samp_perturb) Isim_REITER(samp_perturb) Nsim_REITER(samp_perturb) psim_REITER(samp_perturb)]);
data_WINBERRY = log([Ysim_WINBERRY(samp_perturb) Isim_WINBERRY(samp_perturb) Nsim_WINBERRY(samp_perturb) psim_WINBERRY(samp_perturb)]);

UNCOND_DIFFS_TABLE(1,:) = 100*mean(data_XPA-data_KS,1);
UNCOND_DIFFS_TABLE(2,:) = 100*mean(data_PARAM-data_KS,1);
UNCOND_DIFFS_TABLE(3,:) = 100*mean(data_REITER-data_KS,1);
UNCOND_DIFFS_TABLE(4,:) = 100*mean(data_WINBERRY-data_KS,1);

IRF_DIFFS_TABLE = zeros(4,4);

IRFmeansamp = shockperIRF:(numperIRF-1);

IRFdata_KS = [YIRF_KS(IRFmeansamp) IIRF_KS(IRFmeansamp) NIRF_KS(IRFmeansamp) pIRF_KS(IRFmeansamp)];
IRFdata_XPA = [YIRF_XPA(IRFmeansamp) IIRF_XPA(IRFmeansamp) NIRF_XPA(IRFmeansamp) pIRF_XPA(IRFmeansamp)];
IRFdata_PARAM = [YIRF_PARAM(IRFmeansamp) IIRF_PARAM(IRFmeansamp) NIRF_PARAM(IRFmeansamp) pIRF_PARAM(IRFmeansamp)];
IRFdata_REITER = [YIRF_REITER(IRFmeansamp) IIRF_REITER(IRFmeansamp) NIRF_REITER(IRFmeansamp) pIRF_REITER(IRFmeansamp)];
IRFdata_WINBERRY = [YIRF_WINBERRY(IRFmeansamp) IIRF_WINBERRY(IRFmeansamp) NIRF_WINBERRY(IRFmeansamp) pIRF_WINBERRY(IRFmeansamp)];

IRF_DIFFS_TABLE(1,:) = mean(IRFdata_XPA-IRFdata_KS,1);
IRF_DIFFS_TABLE(2,:) = mean(IRFdata_PARAM-IRFdata_KS,1);
IRF_DIFFS_TABLE(3,:) = mean(IRFdata_REITER-IRFdata_KS,1);
IRF_DIFFS_TABLE(4,:) = mean(IRFdata_WINBERRY-IRFdata_KS,1);


%%%%%%%%%%%%%%%%%%
%%%%%COMPUTE ACCURACY STATS AND PRODUCE FCST FIGURES WITH A SEPARATE SIM
%%%%%%%%%%%%%%%%%%

asimpos_NEWSEED=importdata('./NEW_SEED_DH/asimpos.txt');

%KS
Kbarsim_KS_NEWSEED = importdata('./NEW_SEED_DH/KS/kbarsim.txt');
KbarFCST_KS_NEWSEED = importdata('./NEW_SEED_DH/KS/kbarfcstsim.txt');
KbarDH_KS_NEWSEED = importdata('./NEW_SEED_DH/KS/kbaraggrule.txt');
pFCST_KS_NEWSEED = importdata('./NEW_SEED_DH/KS/pfcstsim.txt');
pDH_KS_NEWSEED = importdata('./NEW_SEED_DH/KS/paggrule.txt');
psim_KS_NEWSEED = importdata('./NEW_SEED_DH/KS/psim.txt');

%XPA
Kbarsim_XPA_NEWSEED = importdata('./NEW_SEED_DH/XPA/kbarsim.txt');
KbarFCST_XPA_NEWSEED = importdata('./NEW_SEED_DH/XPA/kbarfcstsim.txt');
KbarDH_XPA_NEWSEED = importdata('./NEW_SEED_DH/XPA/kbaraggrule.txt');
pFCST_XPA_NEWSEED = importdata('./NEW_SEED_DH/XPA/pfcstsim.txt');
pDH_XPA_NEWSEED = importdata('./NEW_SEED_DH/XPA/paggrule.txt');
psim_XPA_NEWSEED = importdata('./NEW_SEED_DH/XPA/psim.txt');

%PARAM
Kbarsim_PARAM_NEWSEED = importdata('./NEW_SEED_DH/PARAM/kbarsim.txt');
KbarFCST_PARAM_NEWSEED = importdata('./NEW_SEED_DH/PARAM/kbarfcstsim.txt');
KbarDH_PARAM_NEWSEED = importdata('./NEW_SEED_DH/PARAM/kbarDHsim.txt');
pFCST_PARAM_NEWSEED = importdata('./NEW_SEED_DH/PARAM/pfcstsim.txt');
pDH_PARAM_NEWSEED = importdata('./NEW_SEED_DH/PARAM/pDHsim.txt');
psim_PARAM_NEWSEED = importdata('./NEW_SEED_DH/PARAM/psim.txt');

figure;
subplot(3,2,1);
plot(plotsamp,-log(pDH_KS_NEWSEED(plotsamp)),'r',...
    plotsamp,-log(pFCST_KS_NEWSEED(plotsamp)),'g',...
    plotsamp,-log(psim_KS_NEWSEED(plotsamp)),'k',...
    'LineWidth',lwidnum)
title('KS: Consumption','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
axis([-Inf Inf -1 -0.3])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;
legend({'Dynamic','Static','Actual'},'Location','NorthEast','FontSize',10);
legend boxoff;

subplot(3,2,3);
plot(plotsamp,-log(pDH_XPA_NEWSEED(plotsamp)),'r',...
    plotsamp,-log(pFCST_XPA_NEWSEED(plotsamp)),'g',...
    plotsamp,-log(psim_XPA_NEWSEED(plotsamp)),'k',...
    'LineWidth',lwidnum)
title('XPA: Consumption','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
axis([-Inf Inf -0.85 -0.725])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(3,2,5);
plot(plotsamp,-log(pDH_PARAM_NEWSEED(plotsamp)),'r',...
    plotsamp,-log(pFCST_PARAM_NEWSEED(plotsamp)),'g',...
    plotsamp,-log(psim_PARAM_NEWSEED(plotsamp)),'k',...
    'LineWidth',lwidnum)
title('PARAM: Consumption','FontSize',titlesizenum)
xlabel('Year','FontSize',labelsizenum)
ylabel('Log','FontSize',labelsizenum)
axis([-Inf Inf -0.85 -0.725])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;


subplot(3,2,2);
plot(plotsamp,log(KbarDH_KS_NEWSEED(plotsamp)),'r',...
    plotsamp,log(KbarFCST_KS_NEWSEED(plotsamp)),'g',...
    plotsamp,log(Kbarsim_KS_NEWSEED(plotsamp)),'k',...
    'LineWidth',lwidnum)
title('KS: Capital','FontSize',titlesizenum)
axis([-Inf Inf 0.35 0.55])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(3,2,4);
plot(plotsamp,log(KbarDH_XPA_NEWSEED(plotsamp)),'r',...
    plotsamp,log(KbarFCST_XPA_NEWSEED(plotsamp)),'g',...
    plotsamp,log(Kbarsim_XPA_NEWSEED(plotsamp)),'k',...
    'LineWidth',lwidnum)
title('XPA: Capital','FontSize',titlesizenum)
axis([-Inf Inf 0.35 0.55])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(3,2,6);
plot(plotsamp,log(KbarDH_PARAM_NEWSEED(plotsamp)),'r',...
    plotsamp,log(KbarFCST_PARAM_NEWSEED(plotsamp)),'g',...
    plotsamp,log(Kbarsim_PARAM_NEWSEED(plotsamp)),'k',...
    'LineWidth',lwidnum)
title('PARAM: Capital','FontSize',titlesizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf 0.35 0.55])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

saveas(gcf,'ACCURACY_FIG_3BY2.pdf');


%%COMPUTE DH STATS
DH_TABLE_NEWSEED = zeros(3,4);

%KS
DH_TABLE_NEWSEED(1,1) = max(100*abs(log(Kbarsim_KS_NEWSEED(samp))-log(KbarDH_KS_NEWSEED(samp))));
DH_TABLE_NEWSEED(1,2) = mean(100*abs(log(Kbarsim_KS_NEWSEED(samp))-log(KbarDH_KS_NEWSEED(samp))));
DH_TABLE_NEWSEED(1,3) = max(100*abs(log(psim_KS_NEWSEED(samp))-log(pDH_KS_NEWSEED(samp))));
DH_TABLE_NEWSEED(1,4) = mean(100*abs(log(psim_KS_NEWSEED(samp))-log(pDH_KS_NEWSEED(samp))));


%XPA
DH_TABLE_NEWSEED(2,1) = max(100*abs(log(Kbarsim_XPA_NEWSEED(samp))-log(KbarDH_XPA_NEWSEED(samp))));
DH_TABLE_NEWSEED(2,2) = mean(100*abs(log(Kbarsim_XPA_NEWSEED(samp))-log(KbarDH_XPA_NEWSEED(samp))));
DH_TABLE_NEWSEED(2,3) = max(100*abs(log(psim_XPA_NEWSEED(samp))-log(pDH_XPA_NEWSEED(samp))));
DH_TABLE_NEWSEED(2,4) = mean(100*abs(log(psim_XPA_NEWSEED(samp))-log(pDH_XPA_NEWSEED(samp))));

%PARAM
DH_TABLE_NEWSEED(3,1) = max(100*abs(log(Kbarsim_PARAM_NEWSEED(samp))-log(KbarDH_PARAM_NEWSEED(samp))));
DH_TABLE_NEWSEED(3,2) = mean(100*abs(log(Kbarsim_PARAM_NEWSEED(samp))-log(KbarDH_PARAM_NEWSEED(samp))));
DH_TABLE_NEWSEED(3,3) = max(100*abs(log(psim_PARAM_NEWSEED(samp))-log(pDH_PARAM_NEWSEED(samp))));
DH_TABLE_NEWSEED(3,4) = mean(100*abs(log(psim_PARAM_NEWSEED(samp))-log(pDH_PARAM_NEWSEED(samp))));

%COMPUTE R^2s & RMSEs
ESS_KS_NEWSEED = zeros(anum,2);
ESS_PARAM_NEWSEED = zeros(anum,2);
ESS_XPA_NEWSEED = zeros(anum,2);

TSS_KS_NEWSEED = zeros(anum,2);
TSS_PARAM_NEWSEED = zeros(anum,2);
TSS_XPA_NEWSEED = zeros(anum,2);

MEANS_KS_NEWSEED = zeros(anum,2);
MEANS_PARAM_NEWSEED = zeros(anum,2);
MEANS_XPA_NEWSEED = zeros(anum,2);

COUNTS_NEWSEED = zeros(anum,1);

%first, compute means
for t = samp(1):samp(end);
    act = asimpos_NEWSEED(t);
    COUNTS_NEWSEED(act) = COUNTS_NEWSEED(act) + 1;
    
    MEANS_KS_NEWSEED(act,:) = MEANS_KS_NEWSEED(act,:) + [log(psim_KS_NEWSEED(t)) log(Kbarsim_KS_NEWSEED(t+1))];
    MEANS_PARAM_NEWSEED(act,:) = MEANS_PARAM_NEWSEED(act,:) + [log(psim_PARAM_NEWSEED(t)) log(Kbarsim_PARAM_NEWSEED(t+1))];
    MEANS_XPA_NEWSEED(act,:) = MEANS_XPA_NEWSEED(act,:) + [log(psim_XPA_NEWSEED(t)) log(Kbarsim_XPA_NEWSEED(t+1))];
    
end;
MEANS_KS_NEWSEED = MEANS_KS_NEWSEED./repmat(COUNTS_NEWSEED,1,2);
MEANS_PARAM_NEWSEED = MEANS_PARAM_NEWSEED./repmat(COUNTS_NEWSEED,1,2);
MEANS_XPA_NEWSEED = MEANS_XPA_NEWSEED./repmat(COUNTS_NEWSEED,1,2);

%then, compute ESS and TSS
for t = samp(1):samp(end);
    act = asimpos_NEWSEED(t);
   
    %ESS
    ESS_KS_NEWSEED(act,:) = ESS_KS_NEWSEED(act,:) + [log(pFCST_KS_NEWSEED(t)/psim_KS_NEWSEED(t)) log(KbarFCST_KS_NEWSEED(t)/Kbarsim_KS_NEWSEED(t))].^2;
    ESS_PARAM_NEWSEED(act,:) = ESS_PARAM_NEWSEED(act,:) + [log(pFCST_PARAM_NEWSEED(t)/psim_PARAM_NEWSEED(t)) log(KbarFCST_PARAM_NEWSEED(t)/Kbarsim_PARAM_NEWSEED(t))].^2;
    ESS_XPA_NEWSEED(act,:) = ESS_XPA_NEWSEED(act,:) + [log(pFCST_XPA_NEWSEED(t)/psim_XPA_NEWSEED(t)) log(KbarFCST_XPA_NEWSEED(t)/Kbarsim_XPA_NEWSEED(t))].^2;
    
    %TSS
    TSS_KS_NEWSEED(act,:) = TSS_KS_NEWSEED(act,:) + [log(MEANS_KS_NEWSEED(act,1)/psim_KS_NEWSEED(t)) log(MEANS_KS_NEWSEED(act,2)/Kbarsim_KS_NEWSEED(t))].^2;
    TSS_PARAM_NEWSEED(act,:) = TSS_PARAM_NEWSEED(act,:) + [log(MEANS_PARAM_NEWSEED(act,1)/psim_PARAM_NEWSEED(t)) log(MEANS_PARAM_NEWSEED(act,2)/Kbarsim_PARAM_NEWSEED(t))].^2;
    TSS_XPA_NEWSEED(act,:) = TSS_XPA_NEWSEED(act,:) + [log(MEANS_XPA_NEWSEED(act,1)/psim_XPA_NEWSEED(t)) log(MEANS_XPA_NEWSEED(act,2)/Kbarsim_XPA_NEWSEED(t))].^2;
end;

%now, actually implement the R^2 and RMSE formulas
R2_KS_NEWSEED = 1 - ESS_KS_NEWSEED./TSS_KS_NEWSEED;
R2_PARAM_NEWSEED = 1 - ESS_PARAM_NEWSEED./TSS_PARAM_NEWSEED;
R2_XPA_NEWSEED = 1 - ESS_XPA_NEWSEED./TSS_XPA_NEWSEED;

RMSE_KS_NEWSEED = 100*sqrt(ESS_KS_NEWSEED./repmat(COUNTS_NEWSEED,1,2));
RMSE_PARAM_NEWSEED = 100*sqrt(ESS_PARAM_NEWSEED./repmat(COUNTS_NEWSEED,1,2));
RMSE_XPA_NEWSEED = 100*sqrt(ESS_XPA_NEWSEED./repmat(COUNTS_NEWSEED,1,2));

%%%%%%%%%%%%%%%%%%
%%%%%READ IN AND COMPARE SIMS FROM KS WITH MAINTENANCE INVESTMENT
%%%%%%%%%%%%%%%%%%

Ysim_KS_MAINT = importdata('./KS/MAINTENANCE_INV/ysim.txt');
Isim_KS_MAINT = importdata('./KS/MAINTENANCE_INV/isim.txt');
Nsim_KS_MAINT = importdata('./KS/MAINTENANCE_INV/Nsim.txt');
psim_KS_MAINT = importdata('./KS/MAINTENANCE_INV/psim.txt');


%read in raw data on micro moments over uncondtional sim
MICROsim_KS_MAINT = importdata('./KS/MAINTENANCE_INV/MICROsim.txt');
MICROsim_KS_MAINT = reshape(MICROsim_KS_MAINT,nummicro,numper);
MICRO_TABLE_KS_MAINT = mean(MICROsim_KS_MAINT(:,samp),2);

%make unconditional simulation figure

figure;
subplot(2,2,1);
plot(plotsamp,log(Ysim_KS_MAINT(plotsamp)),'r',...
    plotsamp,log(Ysim_KS(plotsamp)),'k',...
    'LineWidth',lwidnum); 
title('Output','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
axis([-Inf Inf -0.75 -0.3])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;
legend('Maintenance Inv','Baseline','Location','NorthEast');
legend boxoff;

subplot(2,2,2);
plot(plotsamp,log(Isim_KS_MAINT(plotsamp)),'r',...
    plotsamp,log(Isim_KS(plotsamp)),'k',...
    'LineWidth',lwidnum);
title('Investment','FontSize',titlesizenum)
axis([-Inf Inf -3.0 -1.6])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(2,2,3);
plot(plotsamp,log(Nsim_KS_MAINT(plotsamp)),'r',...
    plotsamp,log(Nsim_KS(plotsamp)),'k',...
    'LineWidth',lwidnum);
title('Labor','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -1.225 -0.975])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(2,2,4);
plot(plotsamp,-log(psim_KS_MAINT(plotsamp)),'r',...
    plotsamp,-log(psim_KS(plotsamp)),'k',...
    'LineWidth',lwidnum);
title('Consumption','FontSize',titlesizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -0.9 -0.7])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

saveas(gcf,'UNCOND_KS_AND_MAINT.pdf')

%%%%%%%%%%%%%%%%%%
%%%%%READ IN AND COMPARE SERIES WITH DIFFERENT SIZE SHOCKS
%%%%%%%%%%%%%%%%%%

a0_5 = importdata('./DIFF_SIZE_SHOCKS/KS_5/a0.txt');
a0_10 = importdata('./DIFF_SIZE_SHOCKS/KS_10/a0.txt');
a0_25 = importdata('./DIFF_SIZE_SHOCKS/KS_25/a0.txt');
a0_50 = importdata('./DIFF_SIZE_SHOCKS/KS_50/a0.txt');
a0_150 = importdata('./DIFF_SIZE_SHOCKS/KS_150/a0.txt');
a0_200 = importdata('./DIFF_SIZE_SHOCKS/KS_200/a0.txt');

asimpos_5 = importdata('./DIFF_SIZE_SHOCKS/KS_5/asimpos.txt');
asimpos_10 = importdata('./DIFF_SIZE_SHOCKS/KS_10/asimpos.txt');
asimpos_25 = importdata('./DIFF_SIZE_SHOCKS/KS_25/asimpos.txt');
asimpos_50 = importdata('./DIFF_SIZE_SHOCKS/KS_50/asimpos.txt');
asimpos_150 = importdata('./DIFF_SIZE_SHOCKS/KS_150/asimpos.txt');
asimpos_200 = importdata('./DIFF_SIZE_SHOCKS/KS_200/asimpos.txt');

asim_5  = a0_5(asimpos_5);
asim_10  = a0_10(asimpos_10);
asim_25  = a0_25(asimpos_25);
asim_50  = a0_50(asimpos_50);
asim_150  = a0_150(asimpos_150);
asim_200  = a0_200(asimpos_200);

logasim_5  = log(a0_5(asimpos_5));
logasim_10  = log(a0_10(asimpos_10));
logasim_25  = log(a0_25(asimpos_25));
logasim_50  = log(a0_50(asimpos_50));
logasim_150  = log(a0_150(asimpos_150));
logasim_200  = log(a0_200(asimpos_200));

Ysim_KS_5 = importdata('./DIFF_SIZE_SHOCKS/KS_5/ysim.txt');
Ysim_KS_10 = importdata('./DIFF_SIZE_SHOCKS/KS_10/ysim.txt');
Ysim_KS_25 = importdata('./DIFF_SIZE_SHOCKS/KS_25/ysim.txt');
Ysim_KS_50 = importdata('./DIFF_SIZE_SHOCKS/KS_50/ysim.txt');
Ysim_KS_150 = importdata('./DIFF_SIZE_SHOCKS/KS_150/ysim.txt');
Ysim_KS_200 = importdata('./DIFF_SIZE_SHOCKS/KS_200/ysim.txt');

Isim_KS_5 = importdata('./DIFF_SIZE_SHOCKS/KS_5/isim.txt');
Isim_KS_10 = importdata('./DIFF_SIZE_SHOCKS/KS_10/isim.txt');
Isim_KS_25 = importdata('./DIFF_SIZE_SHOCKS/KS_25/isim.txt');
Isim_KS_50 = importdata('./DIFF_SIZE_SHOCKS/KS_50/isim.txt');
Isim_KS_150 = importdata('./DIFF_SIZE_SHOCKS/KS_150/isim.txt');
Isim_KS_200 = importdata('./DIFF_SIZE_SHOCKS/KS_200/isim.txt');

Nsim_KS_5 = importdata('./DIFF_SIZE_SHOCKS/KS_5/Nsim.txt');
Nsim_KS_10 = importdata('./DIFF_SIZE_SHOCKS/KS_10/Nsim.txt');
Nsim_KS_25 = importdata('./DIFF_SIZE_SHOCKS/KS_25/Nsim.txt');
Nsim_KS_50 = importdata('./DIFF_SIZE_SHOCKS/KS_50/Nsim.txt');
Nsim_KS_150 = importdata('./DIFF_SIZE_SHOCKS/KS_150/Nsim.txt');
Nsim_KS_200 = importdata('./DIFF_SIZE_SHOCKS/KS_200/Nsim.txt');

psim_KS_5 = importdata('./DIFF_SIZE_SHOCKS/KS_5/psim.txt');
psim_KS_10 = importdata('./DIFF_SIZE_SHOCKS/KS_10/psim.txt');
psim_KS_25 = importdata('./DIFF_SIZE_SHOCKS/KS_25/psim.txt');
psim_KS_50 = importdata('./DIFF_SIZE_SHOCKS/KS_50/psim.txt');
psim_KS_150 = importdata('./DIFF_SIZE_SHOCKS/KS_150/psim.txt');
psim_KS_200 = importdata('./DIFF_SIZE_SHOCKS/KS_200/psim.txt');

Ysim_REITER_5 = importdata('./DIFF_SIZE_SHOCKS/REITER/ysim5.txt');
Ysim_REITER_10 = importdata('./DIFF_SIZE_SHOCKS/REITER/ysim10.txt');
Ysim_REITER_25 = importdata('./DIFF_SIZE_SHOCKS/REITER/ysim25.txt');
Ysim_REITER_50 = importdata('./DIFF_SIZE_SHOCKS/REITER/ysim50.txt');
Ysim_REITER_150 = importdata('./DIFF_SIZE_SHOCKS/REITER/ysim150.txt');
Ysim_REITER_200 = importdata('./DIFF_SIZE_SHOCKS/REITER/ysim200.txt');

Isim_REITER_5 = importdata('./DIFF_SIZE_SHOCKS/REITER/isim5.txt');
Isim_REITER_10 = importdata('./DIFF_SIZE_SHOCKS/REITER/isim10.txt');
Isim_REITER_25 = importdata('./DIFF_SIZE_SHOCKS/REITER/isim25.txt');
Isim_REITER_50 = importdata('./DIFF_SIZE_SHOCKS/REITER/isim50.txt');
Isim_REITER_150 = importdata('./DIFF_SIZE_SHOCKS/REITER/isim150.txt');
Isim_REITER_200 = importdata('./DIFF_SIZE_SHOCKS/REITER/isim200.txt');

Nsim_REITER_5 = importdata('./DIFF_SIZE_SHOCKS/REITER/Nsim5.txt');
Nsim_REITER_10 = importdata('./DIFF_SIZE_SHOCKS/REITER/Nsim10.txt');
Nsim_REITER_25 = importdata('./DIFF_SIZE_SHOCKS/REITER/Nsim25.txt');
Nsim_REITER_50 = importdata('./DIFF_SIZE_SHOCKS/REITER/Nsim50.txt');
Nsim_REITER_150 = importdata('./DIFF_SIZE_SHOCKS/REITER/Nsim150.txt');
Nsim_REITER_200 = importdata('./DIFF_SIZE_SHOCKS/REITER/Nsim200.txt');

psim_REITER_5 = importdata('./DIFF_SIZE_SHOCKS/REITER/psim5.txt');
psim_REITER_10 = importdata('./DIFF_SIZE_SHOCKS/REITER/psim10.txt');
psim_REITER_25 = importdata('./DIFF_SIZE_SHOCKS/REITER/psim25.txt');
psim_REITER_50 = importdata('./DIFF_SIZE_SHOCKS/REITER/psim50.txt');
psim_REITER_150 = importdata('./DIFF_SIZE_SHOCKS/REITER/psim150.txt');
psim_REITER_200 = importdata('./DIFF_SIZE_SHOCKS/REITER/psim200.txt');

figure;
plot(plotsamp,-log(psim_KS_5(plotsamp)),'r',...
     plotsamp,-log(psim_KS_25(plotsamp)),'g',...
     plotsamp,-log(psim_KS(plotsamp)),'k',...
     plotsamp,-log(psim_KS_150(plotsamp)),'b',...
     plotsamp,-log(psim_KS_200(plotsamp)),'m',...
     plotsamp,-log(psim_REITER_5(plotsamp_perturb)),'r--',...
     plotsamp,-log(psim_REITER_25(plotsamp_perturb)),'g--',...
     plotsamp,-log(psim_REITER(plotsamp_perturb)),'k--',...
     plotsamp,-log(psim_REITER_150(plotsamp_perturb)),'b--',...
     plotsamp,-log(psim_REITER_200(plotsamp_perturb)),'m--',...
     'LineWidth',lwidnum);
title('Consumption','FontSize',titlesizenum)
xlabel('Year','FontSize',labelsizenum)
ylabel('Log','FontSize',labelsizenum)
axis([-Inf Inf -0.95 -0.65])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;
legend('5%','25%','Baseline','150%','200%','Location','NorthEast');
legend boxoff;
saveas(gcf,'DIFF_SIZE_SHOCKS.pdf')
 
%do KS HP filtering
HPdata_KS_5 = log([Ysim_KS_5(samp) Isim_KS_5(samp) Nsim_KS_5(samp) asim_5(samp) psim_KS_5(samp)]);
HPdata_KS_25 = log([Ysim_KS_25(samp) Isim_KS_25(samp) Nsim_KS_25(samp) asim_25(samp) psim_KS_25(samp)]);
HPdata_KS_50 = log([Ysim_KS_50(samp) Isim_KS_50(samp) Nsim_KS_50(samp) asim_50(samp) psim_KS_50(samp)]);
HPdata_KS_150 = log([Ysim_KS_150(samp) Isim_KS_150(samp) Nsim_KS_150(samp) asim_150(samp) psim_KS_150(samp)]);
HPdata_KS_200 = log([Ysim_KS_200(samp) Isim_KS_200(samp) Nsim_KS_200(samp) asim_200(samp) psim_KS_200(samp)]);

[HPtrend_KS_25,HPdetrend_KS_25] = hptrend(HPdata_KS_25,100);
[HPtrend_KS_50,HPdetrend_KS_50] = hptrend(HPdata_KS_50,100);
[HPtrend_KS_150,HPdetrend_KS_150] = hptrend(HPdata_KS_150,100);
[HPtrend_KS_200,HPdetrend_KS_200] = hptrend(HPdata_KS_200,100);

HP_VOL_TABLE_KS = zeros(5,5);
HP_VOL_TABLE_KS(1,:) = 100*var(HPdetrend_KS_25).^0.5; HP_VOL_TABLE_KS(1,2:end)=HP_VOL_TABLE_KS(1,2:end)./HP_VOL_TABLE_KS(1,1);
HP_VOL_TABLE_KS(2,:) = 100*var(HPdetrend_KS_50).^0.5; HP_VOL_TABLE_KS(2,2:end)=HP_VOL_TABLE_KS(2,2:end)./HP_VOL_TABLE_KS(2,1);
HP_VOL_TABLE_KS(3,:) = HP_VOL_TABLE(1,:);
HP_VOL_TABLE_KS(4,:) = 100*var(HPdetrend_KS_150).^0.5; HP_VOL_TABLE_KS(4,2:end)=HP_VOL_TABLE_KS(4,2:end)./HP_VOL_TABLE_KS(4,1);
HP_VOL_TABLE_KS(5,:) = 100*var(HPdetrend_KS_200).^0.5; HP_VOL_TABLE_KS(5,2:end)=HP_VOL_TABLE_KS(5,2:end)./HP_VOL_TABLE_KS(5,1);

HP_CORR_TABLE_KS = zeros(5,5);

corrmat = corrcoef(HPdetrend_KS_25); HP_CORR_TABLE_KS(1,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_KS_50); HP_CORR_TABLE_KS(2,:) = corrmat(1,:);
HP_CORR_TABLE_KS(3,:) = HP_CORR_TABLE(1,:);
corrmat = corrcoef(HPdetrend_KS_150); HP_CORR_TABLE_KS(4,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_KS_200); HP_CORR_TABLE_KS(5,:) = corrmat(1,:);


%do REITER HP filtering
HPdata_REITER_5 = log([Ysim_REITER_5(samp_perturb) Isim_REITER_5(samp_perturb) Nsim_REITER_5(samp_perturb) asim_5(samp) psim_REITER_5(samp_perturb)]);
HPdata_REITER_10 = log([Ysim_REITER_10(samp_perturb) Isim_REITER_10(samp_perturb) Nsim_REITER_10(samp_perturb) asim_10(samp) psim_REITER_10(samp_perturb)]);
HPdata_REITER_25 = log([Ysim_REITER_25(samp_perturb) Isim_REITER_25(samp_perturb) Nsim_REITER_25(samp_perturb) asim_25(samp) psim_REITER_25(samp_perturb)]);
HPdata_REITER_50 = log([Ysim_REITER_50(samp_perturb) Isim_REITER_50(samp_perturb) Nsim_REITER_50(samp_perturb) asim_50(samp) psim_REITER_50(samp_perturb)]);
HPdata_REITER_150 = log([Ysim_REITER_150(samp_perturb) Isim_REITER_150(samp_perturb) Nsim_REITER_150(samp_perturb) asim_150(samp) psim_REITER_150(samp_perturb)]);
HPdata_REITER_200 = log([Ysim_REITER_200(samp_perturb) Isim_REITER_200(samp_perturb) Nsim_REITER_200(samp_perturb) asim_200(samp) psim_REITER_200(samp_perturb)]);

[HPtrend_REITER_5,HPdetrend_REITER_5] = hptrend(HPdata_REITER_5,100);
[HPtrend_REITER_10,HPdetrend_REITER_10] = hptrend(HPdata_REITER_10,100);
[HPtrend_REITER_25,HPdetrend_REITER_25] = hptrend(HPdata_REITER_25,100);
[HPtrend_REITER_50,HPdetrend_REITER_50] = hptrend(HPdata_REITER_50,100);
[HPtrend_REITER_150,HPdetrend_REITER_150] = hptrend(HPdata_REITER_150,100);
[HPtrend_REITER_200,HPdetrend_REITER_200] = hptrend(HPdata_REITER_200,100);

HP_VOL_TABLE_REITER = zeros(7,5);

HP_VOL_TABLE_REITER(1,:) = 100*var(HPdetrend_REITER_5).^0.5;  HP_VOL_TABLE_REITER(1,2:end)=HP_VOL_TABLE_REITER(1,2:end)./HP_VOL_TABLE_REITER(1,1);
HP_VOL_TABLE_REITER(2,:) = 100*var(HPdetrend_REITER_10).^0.5; HP_VOL_TABLE_REITER(2,2:end)=HP_VOL_TABLE_REITER(2,2:end)./HP_VOL_TABLE_REITER(2,1);
HP_VOL_TABLE_REITER(3,:) = 100*var(HPdetrend_REITER_25).^0.5; HP_VOL_TABLE_REITER(3,2:end)=HP_VOL_TABLE_REITER(3,2:end)./HP_VOL_TABLE_REITER(3,1);
HP_VOL_TABLE_REITER(4,:) = 100*var(HPdetrend_REITER_50).^0.5; HP_VOL_TABLE_REITER(4,2:end)=HP_VOL_TABLE_REITER(4,2:end)./HP_VOL_TABLE_REITER(4,1);
HP_VOL_TABLE_REITER(5,:) = HP_VOL_TABLE(4,:);
HP_VOL_TABLE_REITER(6,:) = 100*var(HPdetrend_REITER_150).^0.5; HP_VOL_TABLE_REITER(6,2:end)=HP_VOL_TABLE_REITER(6,2:end)./HP_VOL_TABLE_REITER(6,1);
HP_VOL_TABLE_REITER(7,:) = 100*var(HPdetrend_REITER_200).^0.5; HP_VOL_TABLE_REITER(7,2:end)=HP_VOL_TABLE_REITER(7,2:end)./HP_VOL_TABLE_REITER(7,1);

HP_CORR_TABLE_REITER = zeros(7,5);

corrmat = corrcoef(HPdetrend_REITER_5);  HP_CORR_TABLE_REITER(1,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_REITER_10); HP_CORR_TABLE_REITER(2,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_REITER_25); HP_CORR_TABLE_REITER(3,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_REITER_50); HP_CORR_TABLE_REITER(4,:) = corrmat(1,:);
HP_CORR_TABLE_REITER(5,:) = HP_CORR_TABLE(4,:);
corrmat = corrcoef(HPdetrend_REITER_150); HP_CORR_TABLE_REITER(6,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_REITER_200); HP_CORR_TABLE_REITER(7,:) = corrmat(1,:);

%max shock size differences
MAX_SHOCK_SIZE_DIFFS_TABLE = 100*max(abs(HPdata_KS_5-HPdata_REITER_5));
MAX_SHOCK_SIZE_DIFFS_TABLE = [MAX_SHOCK_SIZE_DIFFS_TABLE; 100*max(abs(HPdata_KS_25-HPdata_REITER_25))];
MAX_SHOCK_SIZE_DIFFS_TABLE = [MAX_SHOCK_SIZE_DIFFS_TABLE; 100*max(abs(HPdata_KS-HPdata_REITER))];
MAX_SHOCK_SIZE_DIFFS_TABLE = [MAX_SHOCK_SIZE_DIFFS_TABLE; 100*max(abs(HPdata_KS_150-HPdata_REITER_150))];
MAX_SHOCK_SIZE_DIFFS_TABLE = [MAX_SHOCK_SIZE_DIFFS_TABLE; 100*max(abs(HPdata_KS_200-HPdata_REITER_200))];

%mean shock size differences
MEAN_SHOCK_SIZE_DIFFS_TABLE = 100*mean(abs(HPdata_KS_5-HPdata_REITER_5));
MEAN_SHOCK_SIZE_DIFFS_TABLE = [MEAN_SHOCK_SIZE_DIFFS_TABLE; 100*mean(abs(HPdata_KS_25-HPdata_REITER_25))];
MEAN_SHOCK_SIZE_DIFFS_TABLE = [MEAN_SHOCK_SIZE_DIFFS_TABLE; 100*mean(abs(HPdata_KS-HPdata_REITER))];
MEAN_SHOCK_SIZE_DIFFS_TABLE = [MEAN_SHOCK_SIZE_DIFFS_TABLE; 100*mean(abs(HPdata_KS_150-HPdata_REITER_150))];
MEAN_SHOCK_SIZE_DIFFS_TABLE = [MEAN_SHOCK_SIZE_DIFFS_TABLE; 100*mean(abs(HPdata_KS_200-HPdata_REITER_200))];

%%%%%%%%%%%%%%%%%%
%%%%%READ IN SUBSIDY DATA
%%%%%%%%%%%%%%%%%%
    
Ysim_KS_SUB_10 = importdata('./SUBSIDY/KS/10/ysim.txt');
Isim_KS_SUB_10 = importdata('./SUBSIDY/KS/10/isim.txt');
Nsim_KS_SUB_10 = importdata('./SUBSIDY/KS/10/Nsim.txt');
psim_KS_SUB_10 = importdata('./SUBSIDY/KS/10/psim.txt');
Kbarsim_KS_SUB_10 = importdata('./SUBSIDY/KS/10/kbarsim.txt');
KbarFCST_KS_SUB_10 = importdata('./SUBSIDY/KS/10/kbarfcstsim.txt');
KbarDH_KS_SUB_10 = importdata('./SUBSIDY/KS/10/kbaraggrule.txt');
pFCST_KS_SUB_10 = importdata('./SUBSIDY/KS/10/pfcstsim.txt');
pDH_KS_SUB_10 = importdata('./SUBSIDY/KS/10/paggrule.txt');

Ysim_KS_SUB_25 = importdata('./SUBSIDY/KS/25/ysim.txt');
Isim_KS_SUB_25 = importdata('./SUBSIDY/KS/25/isim.txt');
Nsim_KS_SUB_25 = importdata('./SUBSIDY/KS/25/Nsim.txt');
psim_KS_SUB_25 = importdata('./SUBSIDY/KS/25/psim.txt');
Kbarsim_KS_SUB_25 = importdata('./SUBSIDY/KS/25/kbarsim.txt');
KbarFCST_KS_SUB_25 = importdata('./SUBSIDY/KS/25/kbarfcstsim.txt');
KbarDH_KS_SUB_25 = importdata('./SUBSIDY/KS/25/kbaraggrule.txt');
pFCST_KS_SUB_25 = importdata('./SUBSIDY/KS/25/pfcstsim.txt');
pDH_KS_SUB_25 = importdata('./SUBSIDY/KS/25/paggrule.txt');


Ysim_KS_SUB_33 = importdata('./SUBSIDY/KS/33/ysim.txt');
Isim_KS_SUB_33 = importdata('./SUBSIDY/KS/33/isim.txt');
Nsim_KS_SUB_33 = importdata('./SUBSIDY/KS/33/Nsim.txt');
psim_KS_SUB_33 = importdata('./SUBSIDY/KS/33/psim.txt');
Kbarsim_KS_SUB_33 = importdata('./SUBSIDY/KS/33/kbarsim.txt');
KbarFCST_KS_SUB_33 = importdata('./SUBSIDY/KS/33/kbarfcstsim.txt');
KbarDH_KS_SUB_33 = importdata('./SUBSIDY/KS/33/kbaraggrule.txt');
pFCST_KS_SUB_33 = importdata('./SUBSIDY/KS/33/pfcstsim.txt');
pDH_KS_SUB_33 = importdata('./SUBSIDY/KS/33/paggrule.txt');

Ysim_REITER_SUB_10 = importdata('./SUBSIDY/REITER/10/ysim.txt');
Isim_REITER_SUB_10 = importdata('./SUBSIDY/REITER/10/isim.txt');
Nsim_REITER_SUB_10 = importdata('./SUBSIDY/REITER/10/Nsim.txt');
psim_REITER_SUB_10 = importdata('./SUBSIDY/REITER/10/psim.txt');

Ysim_REITER_SUB_25 = importdata('./SUBSIDY/REITER/25/ysim.txt');
Isim_REITER_SUB_25 = importdata('./SUBSIDY/REITER/25/isim.txt');
Nsim_REITER_SUB_25 = importdata('./SUBSIDY/REITER/25/Nsim.txt');
psim_REITER_SUB_25 = importdata('./SUBSIDY/REITER/25/psim.txt');

Ysim_REITER_SUB_33 = importdata('./SUBSIDY/REITER/33/ysim.txt');
Isim_REITER_SUB_33 = importdata('./SUBSIDY/REITER/33/isim.txt');
Nsim_REITER_SUB_33 = importdata('./SUBSIDY/REITER/33/Nsim.txt');
psim_REITER_SUB_33 = importdata('./SUBSIDY/REITER/33/psim.txt');



figure;
subplot(2,2,1);
plot(plotsamp,log(Ysim_KS_SUB_33(plotsamp)),'b',...
     plotsamp,log(Ysim_KS_SUB_25(plotsamp)),'g',...
     plotsamp,log(Ysim_KS_SUB_10(plotsamp)),'r',...
     plotsamp,log(Ysim_KS(plotsamp)),'k',...
     plotsamp,log(Ysim_REITER_SUB_33(plotsamp_perturb)),'b--',...
     plotsamp,log(Ysim_REITER_SUB_25(plotsamp_perturb)),'g--',...
     plotsamp,log(Ysim_REITER_SUB_10(plotsamp_perturb)),'r--',...
     plotsamp,log(Ysim_REITER(plotsamp_perturb)),'k--',...
        'LineWidth',lwidnum); 
title('Output','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
axis([-Inf Inf -0.75 -0.15])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;
legend('33%','25%','10%','Baseline','Location','NorthEast');
legend boxoff;

subplot(2,2,2);
plot(plotsamp,log(Isim_KS_SUB_33(plotsamp)),'b',...
     plotsamp,log(Isim_KS_SUB_25(plotsamp)),'g',...
     plotsamp,log(Isim_KS_SUB_10(plotsamp)),'r',...
     plotsamp,log(Isim_KS(plotsamp)),'k',...
     plotsamp,log(Isim_REITER_SUB_33(plotsamp_perturb)),'b--',...
     plotsamp,log(Isim_REITER_SUB_25(plotsamp_perturb)),'g--',...
     plotsamp,log(Isim_REITER_SUB_10(plotsamp_perturb)),'r--',...
     plotsamp,log(Isim_REITER(plotsamp_perturb)),'k--',...
        'LineWidth',lwidnum); 
title('Investment','FontSize',titlesizenum)
axis([-Inf Inf -3.0 -1.5])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(2,2,3);
plot(plotsamp,log(Nsim_KS_SUB_33(plotsamp)),'b',...
     plotsamp,log(Nsim_KS_SUB_25(plotsamp)),'g',...
     plotsamp,log(Nsim_KS_SUB_10(plotsamp)),'r',...
     plotsamp,log(Nsim_KS(plotsamp)),'k',...
     plotsamp,log(Nsim_REITER_SUB_33(plotsamp_perturb)),'b--',...
     plotsamp,log(Nsim_REITER_SUB_25(plotsamp_perturb)),'g--',...
     plotsamp,log(Nsim_REITER_SUB_10(plotsamp_perturb)),'r--',...
     plotsamp,log(Nsim_REITER(plotsamp_perturb)),'k--',...
    'LineWidth',lwidnum);
title('Labor','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -1.225 -0.975])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(2,2,4);
plot(plotsamp,-log(psim_KS_SUB_33(plotsamp)),'b',...
     plotsamp,-log(psim_KS_SUB_25(plotsamp)),'g',...
     plotsamp,-log(psim_KS_SUB_10(plotsamp)),'r',...
     plotsamp,-log(psim_KS(plotsamp)),'k',...
     plotsamp,-log(psim_REITER_SUB_33(plotsamp_perturb)),'b--',...
     plotsamp,-log(psim_REITER_SUB_25(plotsamp_perturb)),'g--',...
     plotsamp,-log(psim_REITER_SUB_10(plotsamp_perturb)),'r--',...
     plotsamp,-log(psim_REITER(plotsamp_perturb)),'k--',...
    'LineWidth',lwidnum);
title('Consumption','FontSize',titlesizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -0.9 -0.7])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

saveas(gcf,'SUBSIDY_KS_VS_REITER.pdf')



%%%NOW, READ IN DATA AND COMPUTE IRF FIGURES FOR THE SUBSIDY EXPERIMENTS


%KS
YsimIRF_KS_SUB_10vec = importdata('./SUBSIDY/KS/10/ysimIRF.txt');
YsimIRF_KS_SUB_25vec = importdata('./SUBSIDY/KS/25/ysimIRF.txt');
YsimIRF_KS_SUB_33vec = importdata('./SUBSIDY/KS/33/ysimIRF.txt');

IsimIRF_KS_SUB_10vec = importdata('./SUBSIDY/KS/10/isimIRF.txt');
IsimIRF_KS_SUB_25vec = importdata('./SUBSIDY/KS/25/isimIRF.txt');
IsimIRF_KS_SUB_33vec = importdata('./SUBSIDY/KS/33/isimIRF.txt');

NsimIRF_KS_SUB_10vec = importdata('./SUBSIDY/KS/10/NsimIRF.txt');
NsimIRF_KS_SUB_25vec = importdata('./SUBSIDY/KS/25/NsimIRF.txt');
NsimIRF_KS_SUB_33vec = importdata('./SUBSIDY/KS/33/NsimIRF.txt');

psimIRF_KS_SUB_10vec = importdata('./SUBSIDY/KS/10/psimIRF.txt');
psimIRF_KS_SUB_25vec = importdata('./SUBSIDY/KS/25/psimIRF.txt');
psimIRF_KS_SUB_33vec = importdata('./SUBSIDY/KS/33/psimIRF.txt');

YsimIRF_KS_SUB_10 = zeros(numperIRF,numsimIRF,2);
YsimIRF_KS_SUB_25 = zeros(numperIRF,numsimIRF,2);
YsimIRF_KS_SUB_33 = zeros(numperIRF,numsimIRF,2);

IsimIRF_KS_SUB_10 = zeros(numperIRF,numsimIRF,2);
IsimIRF_KS_SUB_25 = zeros(numperIRF,numsimIRF,2);
IsimIRF_KS_SUB_33 = zeros(numperIRF,numsimIRF,2);

NsimIRF_KS_SUB_10 = zeros(numperIRF,numsimIRF,2);
NsimIRF_KS_SUB_25 = zeros(numperIRF,numsimIRF,2);
NsimIRF_KS_SUB_33 = zeros(numperIRF,numsimIRF,2);

psimIRF_KS_SUB_10 = zeros(numperIRF,numsimIRF,2);
psimIRF_KS_SUB_25 = zeros(numperIRF,numsimIRF,2);
psimIRF_KS_SUB_33 = zeros(numperIRF,numsimIRF,2);

ct=0;
for simct=1:numsimIRF;
for t=1:numperIRF;
for shockct=1:2;
ct=ct+1;

YsimIRF_KS_SUB_10(t,simct,shockct) = YsimIRF_KS_SUB_10vec(ct);    
YsimIRF_KS_SUB_25(t,simct,shockct) = YsimIRF_KS_SUB_25vec(ct);    
YsimIRF_KS_SUB_33(t,simct,shockct) = YsimIRF_KS_SUB_33vec(ct);    

IsimIRF_KS_SUB_10(t,simct,shockct) = IsimIRF_KS_SUB_10vec(ct);    
IsimIRF_KS_SUB_25(t,simct,shockct) = IsimIRF_KS_SUB_25vec(ct);    
IsimIRF_KS_SUB_33(t,simct,shockct) = IsimIRF_KS_SUB_33vec(ct);    

NsimIRF_KS_SUB_10(t,simct,shockct) = NsimIRF_KS_SUB_10vec(ct);    
NsimIRF_KS_SUB_25(t,simct,shockct) = NsimIRF_KS_SUB_25vec(ct);    
NsimIRF_KS_SUB_33(t,simct,shockct) = NsimIRF_KS_SUB_33vec(ct);    

psimIRF_KS_SUB_10(t,simct,shockct) = psimIRF_KS_SUB_10vec(ct);
psimIRF_KS_SUB_25(t,simct,shockct) = psimIRF_KS_SUB_25vec(ct);
psimIRF_KS_SUB_33(t,simct,shockct) = psimIRF_KS_SUB_33vec(ct);

end;
end;
end;

%REITER
YsimIRF_REITER_SUB_10vec = importdata('./SUBSIDY/REITER/10/ysimIRF.txt');
YsimIRF_REITER_SUB_25vec = importdata('./SUBSIDY/REITER/25/ysimIRF.txt');
YsimIRF_REITER_SUB_33vec = importdata('./SUBSIDY/REITER/33/ysimIRF.txt');

IsimIRF_REITER_SUB_10vec = importdata('./SUBSIDY/REITER/10/isimIRF.txt');
IsimIRF_REITER_SUB_25vec = importdata('./SUBSIDY/REITER/25/isimIRF.txt');
IsimIRF_REITER_SUB_33vec = importdata('./SUBSIDY/REITER/33/isimIRF.txt');

NsimIRF_REITER_SUB_10vec = importdata('./SUBSIDY/REITER/10/NsimIRF.txt');
NsimIRF_REITER_SUB_25vec = importdata('./SUBSIDY/REITER/25/NsimIRF.txt');
NsimIRF_REITER_SUB_33vec = importdata('./SUBSIDY/REITER/33/NsimIRF.txt');

psimIRF_REITER_SUB_10vec = importdata('./SUBSIDY/REITER/10/psimIRF.txt');
psimIRF_REITER_SUB_25vec = importdata('./SUBSIDY/REITER/25/psimIRF.txt');
psimIRF_REITER_SUB_33vec = importdata('./SUBSIDY/REITER/33/psimIRF.txt');

YsimIRF_REITER_SUB_10 = zeros(numperIRF,numsimIRF,2);
YsimIRF_REITER_SUB_25 = zeros(numperIRF,numsimIRF,2);
YsimIRF_REITER_SUB_33 = zeros(numperIRF,numsimIRF,2);

IsimIRF_REITER_SUB_10 = zeros(numperIRF,numsimIRF,2);
IsimIRF_REITER_SUB_25 = zeros(numperIRF,numsimIRF,2);
IsimIRF_REITER_SUB_33 = zeros(numperIRF,numsimIRF,2);

NsimIRF_REITER_SUB_10 = zeros(numperIRF,numsimIRF,2);
NsimIRF_REITER_SUB_25 = zeros(numperIRF,numsimIRF,2);
NsimIRF_REITER_SUB_33 = zeros(numperIRF,numsimIRF,2);

psimIRF_REITER_SUB_10 = zeros(numperIRF,numsimIRF,2);
psimIRF_REITER_SUB_25 = zeros(numperIRF,numsimIRF,2);
psimIRF_REITER_SUB_33 = zeros(numperIRF,numsimIRF,2);

ct=0;
for shockct=1:2;
for simct=1:numsimIRF;
for t=1:numperIRF;
ct=ct+1;

YsimIRF_REITER_SUB_10(t,simct,shockct) = YsimIRF_REITER_SUB_10vec(ct);    
YsimIRF_REITER_SUB_25(t,simct,shockct) = YsimIRF_REITER_SUB_25vec(ct);    
YsimIRF_REITER_SUB_33(t,simct,shockct) = YsimIRF_REITER_SUB_33vec(ct);    

IsimIRF_REITER_SUB_10(t,simct,shockct) = IsimIRF_REITER_SUB_10vec(ct);    
IsimIRF_REITER_SUB_25(t,simct,shockct) = IsimIRF_REITER_SUB_25vec(ct);    
IsimIRF_REITER_SUB_33(t,simct,shockct) = IsimIRF_REITER_SUB_33vec(ct);    

NsimIRF_REITER_SUB_10(t,simct,shockct) = NsimIRF_REITER_SUB_10vec(ct);    
NsimIRF_REITER_SUB_25(t,simct,shockct) = NsimIRF_REITER_SUB_25vec(ct);    
NsimIRF_REITER_SUB_33(t,simct,shockct) = NsimIRF_REITER_SUB_33vec(ct);    

psimIRF_REITER_SUB_10(t,simct,shockct) = psimIRF_REITER_SUB_10vec(ct);    
psimIRF_REITER_SUB_25(t,simct,shockct) = psimIRF_REITER_SUB_25vec(ct);    
psimIRF_REITER_SUB_33(t,simct,shockct) = psimIRF_REITER_SUB_33vec(ct);    

end;
end;
end;

%now, process and compute IRFs

%KS
%output
YIRF_KS_SUB_10 = 100*(mean(log(YsimIRF_KS_SUB_10(:,:,2)),2) - mean(log(YsimIRF_KS_SUB_10(:,:,1)),2));
YIRF_KS_SUB_25 = 100*(mean(log(YsimIRF_KS_SUB_25(:,:,2)),2) - mean(log(YsimIRF_KS_SUB_25(:,:,1)),2));
YIRF_KS_SUB_33 = 100*(mean(log(YsimIRF_KS_SUB_33(:,:,2)),2) - mean(log(YsimIRF_KS_SUB_33(:,:,1)),2));

%investment
IIRF_KS_SUB_10 = 100*(mean(log(IsimIRF_KS_SUB_10(:,:,2)),2) - mean(log(IsimIRF_KS_SUB_10(:,:,1)),2));
IIRF_KS_SUB_25 = 100*(mean(log(IsimIRF_KS_SUB_25(:,:,2)),2) - mean(log(IsimIRF_KS_SUB_25(:,:,1)),2));
IIRF_KS_SUB_33 = 100*(mean(log(IsimIRF_KS_SUB_33(:,:,2)),2) - mean(log(IsimIRF_KS_SUB_33(:,:,1)),2));

%labor
NIRF_KS_SUB_10 = 100*(mean(log(NsimIRF_KS_SUB_10(:,:,2)),2) - mean(log(NsimIRF_KS_SUB_10(:,:,1)),2));
NIRF_KS_SUB_25 = 100*(mean(log(NsimIRF_KS_SUB_25(:,:,2)),2) - mean(log(NsimIRF_KS_SUB_25(:,:,1)),2));
NIRF_KS_SUB_33 = 100*(mean(log(NsimIRF_KS_SUB_33(:,:,2)),2) - mean(log(NsimIRF_KS_SUB_33(:,:,1)),2));

%price
pIRF_KS_SUB_10 = 100*(mean(log(psimIRF_KS_SUB_10(:,:,2)),2) - mean(log(psimIRF_KS_SUB_10(:,:,1)),2));
pIRF_KS_SUB_25 = 100*(mean(log(psimIRF_KS_SUB_25(:,:,2)),2) - mean(log(psimIRF_KS_SUB_25(:,:,1)),2));
pIRF_KS_SUB_33 = 100*(mean(log(psimIRF_KS_SUB_33(:,:,2)),2) - mean(log(psimIRF_KS_SUB_33(:,:,1)),2));


%REITER
%output
YIRF_REITER_SUB_10 = 100*(mean(log(YsimIRF_REITER_SUB_10(:,:,2)),2) - mean(log(YsimIRF_REITER_SUB_10(:,:,1)),2));
YIRF_REITER_SUB_25 = 100*(mean(log(YsimIRF_REITER_SUB_25(:,:,2)),2) - mean(log(YsimIRF_REITER_SUB_25(:,:,1)),2));
YIRF_REITER_SUB_33 = 100*(mean(log(YsimIRF_REITER_SUB_33(:,:,2)),2) - mean(log(YsimIRF_REITER_SUB_33(:,:,1)),2));

%investment
IIRF_REITER_SUB_10 = 100*(mean(log(IsimIRF_REITER_SUB_10(:,:,2)),2) - mean(log(IsimIRF_REITER_SUB_10(:,:,1)),2));
IIRF_REITER_SUB_25 = 100*(mean(log(IsimIRF_REITER_SUB_25(:,:,2)),2) - mean(log(IsimIRF_REITER_SUB_25(:,:,1)),2));
IIRF_REITER_SUB_33 = 100*(mean(log(IsimIRF_REITER_SUB_33(:,:,2)),2) - mean(log(IsimIRF_REITER_SUB_33(:,:,1)),2));

%labor
NIRF_REITER_SUB_10 = 100*(mean(log(NsimIRF_REITER_SUB_10(:,:,2)),2) - mean(log(NsimIRF_REITER_SUB_10(:,:,1)),2));
NIRF_REITER_SUB_25 = 100*(mean(log(NsimIRF_REITER_SUB_25(:,:,2)),2) - mean(log(NsimIRF_REITER_SUB_25(:,:,1)),2));
NIRF_REITER_SUB_33 = 100*(mean(log(NsimIRF_REITER_SUB_33(:,:,2)),2) - mean(log(NsimIRF_REITER_SUB_33(:,:,1)),2));

%price
pIRF_REITER_SUB_10 = 100*(mean(log(psimIRF_REITER_SUB_10(:,:,2)),2) - mean(log(psimIRF_REITER_SUB_10(:,:,1)),2));
pIRF_REITER_SUB_25 = 100*(mean(log(psimIRF_REITER_SUB_25(:,:,2)),2) - mean(log(psimIRF_REITER_SUB_25(:,:,1)),2));
pIRF_REITER_SUB_33 = 100*(mean(log(psimIRF_REITER_SUB_33(:,:,2)),2) - mean(log(psimIRF_REITER_SUB_33(:,:,1)),2));


%now, compute HP-filtering table for subsidy data

%KS
HPdata_KS_SUB_10 = log([Ysim_KS_SUB_10(samp) Isim_KS_SUB_10(samp) Nsim_KS_SUB_10(samp) asim(samp) psim_KS_SUB_10(samp)]);
HPdata_KS_SUB_25 = log([Ysim_KS_SUB_25(samp) Isim_KS_SUB_25(samp) Nsim_KS_SUB_25(samp) asim(samp) psim_KS_SUB_25(samp)]);
HPdata_KS_SUB_33 = log([Ysim_KS_SUB_33(samp) Isim_KS_SUB_33(samp) Nsim_KS_SUB_33(samp) asim(samp) psim_KS_SUB_33(samp)]);

[HPtrend_KS_SUB_10,HPdetrend_KS_SUB_10] = hptrend(HPdata_KS_SUB_10,100);
[HPtrend_KS_SUB_25,HPdetrend_KS_SUB_25] = hptrend(HPdata_KS_SUB_25,100);
[HPtrend_KS_SUB_33,HPdetrend_KS_SUB_33] = hptrend(HPdata_KS_SUB_33,100);

HP_VOL_TABLE_KS_SUB = zeros(4,5);
HP_VOL_TABLE_KS_SUB(1,:) = HP_VOL_TABLE(1,:);
HP_VOL_TABLE_KS_SUB(2,:) = 100*var(HPdetrend_KS_SUB_10).^0.5; HP_VOL_TABLE_KS_SUB(2,2:end)=HP_VOL_TABLE_KS_SUB(2,2:end)./HP_VOL_TABLE_KS_SUB(2,1);
HP_VOL_TABLE_KS_SUB(3,:) = 100*var(HPdetrend_KS_SUB_25).^0.5; HP_VOL_TABLE_KS_SUB(3,2:end)=HP_VOL_TABLE_KS_SUB(3,2:end)./HP_VOL_TABLE_KS_SUB(3,1);
HP_VOL_TABLE_KS_SUB(4,:) = 100*var(HPdetrend_KS_SUB_33).^0.5; HP_VOL_TABLE_KS_SUB(4,2:end)=HP_VOL_TABLE_KS_SUB(4,2:end)./HP_VOL_TABLE_KS_SUB(4,1);

HP_CORR_TABLE_KS_SUB = zeros(4,5);
HP_CORR_TABLE_KS_SUB(1,:) = HP_CORR_TABLE(1,:);
corrmat = corrcoef(HPdetrend_KS_SUB_10); HP_CORR_TABLE_KS_SUB(2,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_KS_SUB_25); HP_CORR_TABLE_KS_SUB(3,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_KS_SUB_33); HP_CORR_TABLE_KS_SUB(4,:) = corrmat(1,:);

%Reiter
HPdata_REITER_SUB_10 = log([Ysim_REITER_SUB_10(samp_perturb) Isim_REITER_SUB_10(samp_perturb) Nsim_REITER_SUB_10(samp_perturb) asim(samp) psim_REITER_SUB_10(samp_perturb)]);
HPdata_REITER_SUB_25 = log([Ysim_REITER_SUB_25(samp_perturb) Isim_REITER_SUB_25(samp_perturb) Nsim_REITER_SUB_25(samp_perturb) asim(samp) psim_REITER_SUB_25(samp_perturb)]);
HPdata_REITER_SUB_33 = log([Ysim_REITER_SUB_33(samp_perturb) Isim_REITER_SUB_33(samp_perturb) Nsim_REITER_SUB_33(samp_perturb) asim(samp) psim_REITER_SUB_33(samp_perturb)]);

[HPtrend_REITER_SUB_10,HPdetrend_REITER_SUB_10] = hptrend(HPdata_REITER_SUB_10,100);
[HPtrend_REITER_SUB_25,HPdetrend_REITER_SUB_25] = hptrend(HPdata_REITER_SUB_25,100);
[HPtrend_REITER_SUB_33,HPdetrend_REITER_SUB_33] = hptrend(HPdata_REITER_SUB_33,100);

HP_VOL_TABLE_REITER_SUB = zeros(4,5);
HP_VOL_TABLE_REITER_SUB(1,:) = HP_VOL_TABLE(4,:);
HP_VOL_TABLE_REITER_SUB(2,:) = 100*var(HPdetrend_REITER_SUB_10).^0.5; HP_VOL_TABLE_REITER_SUB(2,2:end)=HP_VOL_TABLE_REITER_SUB(2,2:end)./HP_VOL_TABLE_REITER_SUB(2,1);
HP_VOL_TABLE_REITER_SUB(3,:) = 100*var(HPdetrend_REITER_SUB_25).^0.5; HP_VOL_TABLE_REITER_SUB(3,2:end)=HP_VOL_TABLE_REITER_SUB(3,2:end)./HP_VOL_TABLE_REITER_SUB(3,1);
HP_VOL_TABLE_REITER_SUB(4,:) = 100*var(HPdetrend_REITER_SUB_33).^0.5; HP_VOL_TABLE_REITER_SUB(4,2:end)=HP_VOL_TABLE_REITER_SUB(4,2:end)./HP_VOL_TABLE_REITER_SUB(4,1);

HP_CORR_TABLE_REITER_SUB = zeros(4,5);
HP_CORR_TABLE_REITER_SUB(1,:) = HP_CORR_TABLE(4,:);
corrmat = corrcoef(HPdetrend_REITER_SUB_10); HP_CORR_TABLE_REITER_SUB(2,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_REITER_SUB_25); HP_CORR_TABLE_REITER_SUB(3,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_REITER_SUB_33); HP_CORR_TABLE_REITER_SUB(4,:) = corrmat(1,:);

%plot the relative volatilies
bardata = [HP_VOL_TABLE_KS_SUB(:,1)/HP_VOL_TABLE_KS_SUB(1,1) HP_VOL_TABLE_REITER_SUB(:,1)/HP_VOL_TABLE_REITER_SUB(1,1)];
figure;
bar(bardata,'BarWidth',lwidnum)
axis([-Inf Inf 0.8 1.4])
colormap([0 0 0 ; 0 0 1]);
legend('KS','REITER','Location','NorthWest')
legend boxoff;
title('Output Volatility')
ylabel('Standard Deviation Relative to Baseline','FontSize',labelsizenum)
set(gca,'XTickLabel',{'Baseline','10%','25%','33%'})
set(gca,'FontSize',labelsizenum)
xlabel('Size-Dependent Tax Cyclicality: \gamma_\tau','FontSize',labelsizenum)
saveas(gcf,'SUBSIDY_VOLS.pdf')

%%NOW, COMPUTE DH STATS FOR THESE SUBSIDY EXPERIMENTS (NEED TO DO THIS WITH
%%A DIFFERENT SEED LATER)

DH_TABLE_SUB = zeros(3,4);

%KS
DH_TABLE_SUB(1,1) = max(100*abs(log(Kbarsim_KS_SUB_10(samp))-log(KbarDH_KS_SUB_10(samp))));
DH_TABLE_SUB(1,2) = mean(100*abs(log(Kbarsim_KS_SUB_10(samp))-log(KbarDH_KS_SUB_10(samp))));
DH_TABLE_SUB(1,3) = max(100*abs(log(psim_KS_SUB_10(samp))-log(pDH_KS_SUB_10(samp))));
DH_TABLE_SUB(1,4) = mean(100*abs(log(psim_KS_SUB_10(samp))-log(pDH_KS_SUB_10(samp))));

DH_TABLE_SUB(2,1) = max(100*abs(log(Kbarsim_KS_SUB_25(samp))-log(KbarDH_KS_SUB_25(samp))));
DH_TABLE_SUB(2,2) = mean(100*abs(log(Kbarsim_KS_SUB_25(samp))-log(KbarDH_KS_SUB_25(samp))));
DH_TABLE_SUB(2,3) = max(100*abs(log(psim_KS_SUB_25(samp))-log(pDH_KS_SUB_25(samp))));
DH_TABLE_SUB(2,4) = mean(100*abs(log(psim_KS_SUB_25(samp))-log(pDH_KS_SUB_25(samp))));

DH_TABLE_SUB(3,1) = max(100*abs(log(Kbarsim_KS_SUB_33(samp))-log(KbarDH_KS_SUB_33(samp))));
DH_TABLE_SUB(3,2) = mean(100*abs(log(Kbarsim_KS_SUB_33(samp))-log(KbarDH_KS_SUB_33(samp))));
DH_TABLE_SUB(3,3) = max(100*abs(log(psim_KS_SUB_33(samp))-log(pDH_KS_SUB_33(samp))));
DH_TABLE_SUB(3,4) = mean(100*abs(log(psim_KS_SUB_33(samp))-log(pDH_KS_SUB_33(samp))));


%COMPUTE R^2s & RMSEs
ESS_KS_SUB_10 = zeros(anum,2);
TSS_KS_SUB_10 = zeros(anum,2);
MEANS_KS_SUB_10 = zeros(anum,2);

ESS_KS_SUB_25 = zeros(anum,2);
TSS_KS_SUB_25 = zeros(anum,2);
MEANS_KS_SUB_25 = zeros(anum,2);

ESS_KS_SUB_33 = zeros(anum,2);
TSS_KS_SUB_33 = zeros(anum,2);
MEANS_KS_SUB_33 = zeros(anum,2);

COUNTS_SUB = zeros(anum,1);

%first, compute means
for t = samp(1):samp(end);
    act = asimpos(t);
    COUNTS_SUB(act) = COUNTS_SUB(act) + 1;
    
    MEANS_KS_SUB_10(act,:) = MEANS_KS_SUB_10(act,:) + [log(psim_KS_SUB_10(t)) log(Kbarsim_KS_SUB_10(t+1))];
    MEANS_KS_SUB_25(act,:) = MEANS_KS_SUB_25(act,:) + [log(psim_KS_SUB_25(t)) log(Kbarsim_KS_SUB_25(t+1))];
    MEANS_KS_SUB_33(act,:) = MEANS_KS_SUB_33(act,:) + [log(psim_KS_SUB_33(t)) log(Kbarsim_KS_SUB_33(t+1))];
    
end;
MEANS_KS_SUB_10 = MEANS_KS_SUB_10./repmat(COUNTS_SUB,1,2);
MEANS_KS_SUB_25 = MEANS_KS_SUB_25./repmat(COUNTS_SUB,1,2);
MEANS_KS_SUB_33 = MEANS_KS_SUB_33./repmat(COUNTS_SUB,1,2);

%then, compute ESS and TSS
for t = samp(1):samp(end);
    act = asimpos(t);
   
    %ESS
    ESS_KS_SUB_10(act,:) = ESS_KS_SUB_10(act,:) + [log(pFCST_KS_SUB_10(t)/psim_KS_SUB_10(t)) log(KbarFCST_KS_SUB_10(t)/Kbarsim_KS_SUB_10(t))].^2;
    ESS_KS_SUB_25(act,:) = ESS_KS_SUB_25(act,:) + [log(pFCST_KS_SUB_25(t)/psim_KS_SUB_25(t)) log(KbarFCST_KS_SUB_25(t)/Kbarsim_KS_SUB_25(t))].^2;
    ESS_KS_SUB_33(act,:) = ESS_KS_SUB_33(act,:) + [log(pFCST_KS_SUB_33(t)/psim_KS_SUB_33(t)) log(KbarFCST_KS_SUB_33(t)/Kbarsim_KS_SUB_33(t))].^2;
    
    %TSS
    TSS_KS_SUB_10(act,:) = TSS_KS_SUB_10(act,:) + [log(MEANS_KS_SUB_10(act,1)/psim_KS_SUB_10(t)) log(MEANS_KS_SUB_10(act,2)/Kbarsim_KS_SUB_10(t))].^2;
    TSS_KS_SUB_25(act,:) = TSS_KS_SUB_25(act,:) + [log(MEANS_KS_SUB_25(act,1)/psim_KS_SUB_25(t)) log(MEANS_KS_SUB_25(act,2)/Kbarsim_KS_SUB_25(t))].^2;
    TSS_KS_SUB_33(act,:) = TSS_KS_SUB_33(act,:) + [log(MEANS_KS_SUB_33(act,1)/psim_KS_SUB_33(t)) log(MEANS_KS_SUB_33(act,2)/Kbarsim_KS_SUB_33(t))].^2;
   
end;

%now, actually implement the R^2 and RMSE formulas
R2_KS_SUB_10 = 1 - ESS_KS_SUB_10./TSS_KS_SUB_10;
R2_KS_SUB_25 = 1 - ESS_KS_SUB_25./TSS_KS_SUB_25;
R2_KS_SUB_33 = 1 - ESS_KS_SUB_33./TSS_KS_SUB_33;


RMSE_KS_SUB_10 = 100*sqrt(ESS_KS_SUB_10./repmat(COUNTS_SUB,1,2));
RMSE_KS_SUB_25 = 100*sqrt(ESS_KS_SUB_25./repmat(COUNTS_SUB,1,2));
RMSE_KS_SUB_33 = 100*sqrt(ESS_KS_SUB_33./repmat(COUNTS_SUB,1,2));


%%%%%%%%%%%%%%%%%%
%%%%%COMPUTE ACCURACY STATS AND PRODUCE FCST FIGURES WITH A SEPARATE SIM
%%%%%FOR THE SUBSIDY EXPERIMENT
%%%%%%%%%%%%%%%%%%

%10 percent
Kbarsim_KS_NEWSEED_SUB_10 = importdata('./SUBSIDY/KS/10/NEWSEED/kbarsim.txt');
KbarFCST_KS_NEWSEED_SUB_10 = importdata('./SUBSIDY/KS/10/NEWSEED/kbarfcstsim.txt');
KbarDH_KS_NEWSEED_SUB_10 = importdata('./SUBSIDY/KS/10/NEWSEED/kbaraggrule.txt');
pFCST_KS_NEWSEED_SUB_10 = importdata('./SUBSIDY/KS/10/NEWSEED/pfcstsim.txt');
pDH_KS_NEWSEED_SUB_10 = importdata('./SUBSIDY/KS/10/NEWSEED/paggrule.txt');
psim_KS_NEWSEED_SUB_10 = importdata('./SUBSIDY/KS/10/NEWSEED/psim.txt');

%25 percent
Kbarsim_KS_NEWSEED_SUB_25 = importdata('./SUBSIDY/KS/25/NEWSEED/kbarsim.txt');
KbarFCST_KS_NEWSEED_SUB_25 = importdata('./SUBSIDY/KS/25/NEWSEED/kbarfcstsim.txt');
KbarDH_KS_NEWSEED_SUB_25 = importdata('./SUBSIDY/KS/25/NEWSEED/kbaraggrule.txt');
pFCST_KS_NEWSEED_SUB_25 = importdata('./SUBSIDY/KS/25/NEWSEED/pfcstsim.txt');
pDH_KS_NEWSEED_SUB_25 = importdata('./SUBSIDY/KS/25/NEWSEED/paggrule.txt');
psim_KS_NEWSEED_SUB_25 = importdata('./SUBSIDY/KS/25/NEWSEED/psim.txt');

%33 percent
Kbarsim_KS_NEWSEED_SUB_33 = importdata('./SUBSIDY/KS/33/NEWSEED/kbarsim.txt');
KbarFCST_KS_NEWSEED_SUB_33 = importdata('./SUBSIDY/KS/33/NEWSEED/kbarfcstsim.txt');
KbarDH_KS_NEWSEED_SUB_33 = importdata('./SUBSIDY/KS/33/NEWSEED/kbaraggrule.txt');
pFCST_KS_NEWSEED_SUB_33 = importdata('./SUBSIDY/KS/33/NEWSEED/pfcstsim.txt');
pDH_KS_NEWSEED_SUB_33 = importdata('./SUBSIDY/KS/33/NEWSEED/paggrule.txt');
psim_KS_NEWSEED_SUB_33 = importdata('./SUBSIDY/KS/33/NEWSEED/psim.txt');

%%COMPUTE DH STATS
DH_TABLE_NEWSEED_SUB = zeros(4,4);

%copy KS DH stats for subsidy = 0%
DH_TABLE_NEWSEED_SUB(1,:) = DH_TABLE_NEWSEED(1,:);

%10 percent
DH_TABLE_NEWSEED_SUB(2,1) = max(100*abs(log(Kbarsim_KS_NEWSEED_SUB_10(samp))-log(KbarDH_KS_NEWSEED_SUB_10(samp))));
DH_TABLE_NEWSEED_SUB(2,2) = mean(100*abs(log(Kbarsim_KS_NEWSEED_SUB_10(samp))-log(KbarDH_KS_NEWSEED_SUB_10(samp))));
DH_TABLE_NEWSEED_SUB(2,3) = max(100*abs(log(psim_KS_NEWSEED_SUB_10(samp))-log(pDH_KS_NEWSEED_SUB_10(samp))));
DH_TABLE_NEWSEED_SUB(2,4) = mean(100*abs(log(psim_KS_NEWSEED_SUB_10(samp))-log(pDH_KS_NEWSEED_SUB_10(samp))));

%25 percent
DH_TABLE_NEWSEED_SUB(3,1) = max(100*abs(log(Kbarsim_KS_NEWSEED_SUB_25(samp))-log(KbarDH_KS_NEWSEED_SUB_25(samp))));
DH_TABLE_NEWSEED_SUB(3,2) = mean(100*abs(log(Kbarsim_KS_NEWSEED_SUB_25(samp))-log(KbarDH_KS_NEWSEED_SUB_25(samp))));
DH_TABLE_NEWSEED_SUB(3,3) = max(100*abs(log(psim_KS_NEWSEED_SUB_25(samp))-log(pDH_KS_NEWSEED_SUB_25(samp))));
DH_TABLE_NEWSEED_SUB(3,4) = mean(100*abs(log(psim_KS_NEWSEED_SUB_25(samp))-log(pDH_KS_NEWSEED_SUB_25(samp))));

%33 percent
DH_TABLE_NEWSEED_SUB(4,1) = max(100*abs(log(Kbarsim_KS_NEWSEED_SUB_33(samp))-log(KbarDH_KS_NEWSEED_SUB_33(samp))));
DH_TABLE_NEWSEED_SUB(4,2) = mean(100*abs(log(Kbarsim_KS_NEWSEED_SUB_33(samp))-log(KbarDH_KS_NEWSEED_SUB_33(samp))));
DH_TABLE_NEWSEED_SUB(4,3) = max(100*abs(log(psim_KS_NEWSEED_SUB_33(samp))-log(pDH_KS_NEWSEED_SUB_33(samp))));
DH_TABLE_NEWSEED_SUB(4,4) = mean(100*abs(log(psim_KS_NEWSEED_SUB_33(samp))-log(pDH_KS_NEWSEED_SUB_33(samp))));


%COMPUTE R^2s & RMSEs
ESS_KS_NEWSEED_SUB_10 = zeros(anum,2);
ESS_KS_NEWSEED_SUB_25 = zeros(anum,2);
ESS_KS_NEWSEED_SUB_33 = zeros(anum,2);

TSS_KS_NEWSEED_SUB_10 = zeros(anum,2);
TSS_KS_NEWSEED_SUB_25 = zeros(anum,2);
TSS_KS_NEWSEED_SUB_33 = zeros(anum,2);

MEANS_KS_NEWSEED_SUB_10 = zeros(anum,2);
MEANS_KS_NEWSEED_SUB_25 = zeros(anum,2);
MEANS_KS_NEWSEED_SUB_33 = zeros(anum,2);

%first, compute means
for t = samp(1):samp(end);
    act = asimpos_NEWSEED(t);

    MEANS_KS_NEWSEED_SUB_10(act,:) = MEANS_KS_NEWSEED_SUB_10(act,:) + [log(psim_KS_NEWSEED_SUB_10(t)) log(Kbarsim_KS_NEWSEED_SUB_10(t+1))];
    MEANS_KS_NEWSEED_SUB_25(act,:) = MEANS_KS_NEWSEED_SUB_25(act,:) + [log(psim_KS_NEWSEED_SUB_25(t)) log(Kbarsim_KS_NEWSEED_SUB_25(t+1))];
    MEANS_KS_NEWSEED_SUB_33(act,:) = MEANS_KS_NEWSEED_SUB_33(act,:) + [log(psim_KS_NEWSEED_SUB_33(t)) log(Kbarsim_KS_NEWSEED_SUB_33(t+1))];
    
end;
MEANS_KS_NEWSEED_SUB_10 = MEANS_KS_NEWSEED_SUB_10./repmat(COUNTS_NEWSEED,1,2);
MEANS_KS_NEWSEED_SUB_25 = MEANS_KS_NEWSEED_SUB_25./repmat(COUNTS_NEWSEED,1,2);
MEANS_KS_NEWSEED_SUB_33 = MEANS_KS_NEWSEED_SUB_33./repmat(COUNTS_NEWSEED,1,2);

%then, compute ESS and TSS
for t = samp(1):samp(end);
    act = asimpos_NEWSEED(t);
   
    %ESS
    ESS_KS_NEWSEED_SUB_10(act,:) = ESS_KS_NEWSEED_SUB_10(act,:) + [log(pFCST_KS_NEWSEED_SUB_10(t)/psim_KS_NEWSEED_SUB_10(t)) log(KbarFCST_KS_NEWSEED_SUB_10(t)/Kbarsim_KS_NEWSEED_SUB_10(t))].^2;
    ESS_KS_NEWSEED_SUB_25(act,:) = ESS_KS_NEWSEED_SUB_25(act,:) + [log(pFCST_KS_NEWSEED_SUB_25(t)/psim_KS_NEWSEED_SUB_25(t)) log(KbarFCST_KS_NEWSEED_SUB_25(t)/Kbarsim_KS_NEWSEED_SUB_25(t))].^2;
    ESS_KS_NEWSEED_SUB_33(act,:) = ESS_KS_NEWSEED_SUB_33(act,:) + [log(pFCST_KS_NEWSEED_SUB_33(t)/psim_KS_NEWSEED_SUB_33(t)) log(KbarFCST_KS_NEWSEED_SUB_33(t)/Kbarsim_KS_NEWSEED_SUB_33(t))].^2;
    
    %TSS
    TSS_KS_NEWSEED_SUB_10(act,:) = TSS_KS_NEWSEED_SUB_10(act,:) + [log(MEANS_KS_NEWSEED_SUB_10(act,1)/psim_KS_NEWSEED_SUB_10(t)) log(MEANS_KS_NEWSEED_SUB_10(act,2)/Kbarsim_KS_NEWSEED_SUB_10(t))].^2;
    TSS_KS_NEWSEED_SUB_25(act,:) = TSS_KS_NEWSEED_SUB_25(act,:) + [log(MEANS_KS_NEWSEED_SUB_25(act,1)/psim_KS_NEWSEED_SUB_25(t)) log(MEANS_KS_NEWSEED_SUB_25(act,2)/Kbarsim_KS_NEWSEED_SUB_25(t))].^2;
    TSS_KS_NEWSEED_SUB_33(act,:) = TSS_KS_NEWSEED_SUB_33(act,:) + [log(MEANS_KS_NEWSEED_SUB_33(act,1)/psim_KS_NEWSEED_SUB_33(t)) log(MEANS_KS_NEWSEED_SUB_33(act,2)/Kbarsim_KS_NEWSEED_SUB_33(t))].^2;
    
end;

%now, actually implement the R^2 and RMSE formulas
R2_KS_NEWSEED_SUB_10 = 1 - ESS_KS_NEWSEED_SUB_10./TSS_KS_NEWSEED_SUB_10;
R2_KS_NEWSEED_SUB_25 = 1 - ESS_KS_NEWSEED_SUB_25./TSS_KS_NEWSEED_SUB_25;
R2_KS_NEWSEED_SUB_33 = 1 - ESS_KS_NEWSEED_SUB_33./TSS_KS_NEWSEED_SUB_33;

RMSE_KS_NEWSEED_SUB_10 = 100*sqrt(ESS_KS_NEWSEED_SUB_10./repmat(COUNTS_NEWSEED,1,2));
RMSE_KS_NEWSEED_SUB_25 = 100*sqrt(ESS_KS_NEWSEED_SUB_25./repmat(COUNTS_NEWSEED,1,2));
RMSE_KS_NEWSEED_SUB_33 = 100*sqrt(ESS_KS_NEWSEED_SUB_33./repmat(COUNTS_NEWSEED,1,2));



%%%%%%%%%%%%%%%%%%
%%%%%READ IN & PROCESS IDIO VOL SHOCK DATA
%%%%%%%%%%%%%%%%%%

%%%%UNCONDITIONAL SIMULATION

%KS
Ysim_KS_VOL_10 = importdata('./IDIO_VOL_SHOCKS/KS/10/ysim.txt');
Isim_KS_VOL_10 = importdata('./IDIO_VOL_SHOCKS/KS/10/isim.txt');
Nsim_KS_VOL_10 = importdata('./IDIO_VOL_SHOCKS/KS/10/Nsim.txt');
psim_KS_VOL_10 = importdata('./IDIO_VOL_SHOCKS/KS/10/psim.txt');
Kbarsim_KS_VOL_10 = importdata('./IDIO_VOL_SHOCKS/KS/10/kbarsim.txt');

Ysim_KS_VOL_25 = importdata('./IDIO_VOL_SHOCKS/KS/25/ysim.txt');
Isim_KS_VOL_25 = importdata('./IDIO_VOL_SHOCKS/KS/25/isim.txt');
Nsim_KS_VOL_25 = importdata('./IDIO_VOL_SHOCKS/KS/25/Nsim.txt');
psim_KS_VOL_25 = importdata('./IDIO_VOL_SHOCKS/KS/25/psim.txt');
Kbarsim_KS_VOL_25 = importdata('./IDIO_VOL_SHOCKS/KS/25/kbarsim.txt');

Ysim_KS_VOL_33 = importdata('./IDIO_VOL_SHOCKS/KS/33/ysim.txt');
Isim_KS_VOL_33 = importdata('./IDIO_VOL_SHOCKS/KS/33/isim.txt');
Nsim_KS_VOL_33 = importdata('./IDIO_VOL_SHOCKS/KS/33/Nsim.txt');
psim_KS_VOL_33 = importdata('./IDIO_VOL_SHOCKS/KS/33/psim.txt');
Kbarsim_KS_VOL_33 = importdata('./IDIO_VOL_SHOCKS/KS/33/kbarsim.txt');

Ysim_KS_VOL_50 = importdata('./IDIO_VOL_SHOCKS/KS/50/ysim.txt');
Isim_KS_VOL_50 = importdata('./IDIO_VOL_SHOCKS/KS/50/isim.txt');
Nsim_KS_VOL_50 = importdata('./IDIO_VOL_SHOCKS/KS/50/Nsim.txt');
psim_KS_VOL_50 = importdata('./IDIO_VOL_SHOCKS/KS/50/psim.txt');
Kbarsim_KS_VOL_50 = importdata('./IDIO_VOL_SHOCKS/KS/50/kbarsim.txt');

Ysim_KS_VOL_100 = importdata('./IDIO_VOL_SHOCKS/KS/100/ysim.txt');
Isim_KS_VOL_100 = importdata('./IDIO_VOL_SHOCKS/KS/100/isim.txt');
Nsim_KS_VOL_100 = importdata('./IDIO_VOL_SHOCKS/KS/100/Nsim.txt');
psim_KS_VOL_100 = importdata('./IDIO_VOL_SHOCKS/KS/100/psim.txt');
Kbarsim_KS_VOL_100 = importdata('./IDIO_VOL_SHOCKS/KS/100/kbarsim.txt');


%REITER
Ysim_REITER_VOL_10 = importdata('./IDIO_VOL_SHOCKS/REITER/10/ysim.txt');
Isim_REITER_VOL_10 = importdata('./IDIO_VOL_SHOCKS/REITER/10/isim.txt');
Nsim_REITER_VOL_10 = importdata('./IDIO_VOL_SHOCKS/REITER/10/Nsim.txt');
psim_REITER_VOL_10 = importdata('./IDIO_VOL_SHOCKS/REITER/10/psim.txt');

Ysim_REITER_VOL_25 = importdata('./IDIO_VOL_SHOCKS/REITER/25/ysim.txt');
Isim_REITER_VOL_25 = importdata('./IDIO_VOL_SHOCKS/REITER/25/isim.txt');
Nsim_REITER_VOL_25 = importdata('./IDIO_VOL_SHOCKS/REITER/25/Nsim.txt');
psim_REITER_VOL_25 = importdata('./IDIO_VOL_SHOCKS/REITER/25/psim.txt');

Ysim_REITER_VOL_33 = importdata('./IDIO_VOL_SHOCKS/REITER/33/ysim.txt');
Isim_REITER_VOL_33 = importdata('./IDIO_VOL_SHOCKS/REITER/33/isim.txt');
Nsim_REITER_VOL_33 = importdata('./IDIO_VOL_SHOCKS/REITER/33/Nsim.txt');
psim_REITER_VOL_33 = importdata('./IDIO_VOL_SHOCKS/REITER/33/psim.txt');

Ysim_REITER_VOL_50 = importdata('./IDIO_VOL_SHOCKS/REITER/50/ysim.txt');
Isim_REITER_VOL_50 = importdata('./IDIO_VOL_SHOCKS/REITER/50/isim.txt');
Nsim_REITER_VOL_50 = importdata('./IDIO_VOL_SHOCKS/REITER/50/Nsim.txt');
psim_REITER_VOL_50 = importdata('./IDIO_VOL_SHOCKS/REITER/50/psim.txt');

Ysim_REITER_VOL_100 = importdata('./IDIO_VOL_SHOCKS/REITER/100/ysim.txt');
Isim_REITER_VOL_100 = importdata('./IDIO_VOL_SHOCKS/REITER/100/isim.txt');
Nsim_REITER_VOL_100 = importdata('./IDIO_VOL_SHOCKS/REITER/100/Nsim.txt');
psim_REITER_VOL_100 = importdata('./IDIO_VOL_SHOCKS/REITER/100/psim.txt');

figure;
subplot(2,2,1);
plot(plotsamp,log(Ysim_KS_VOL_100(plotsamp)),'b',...
     plotsamp,log(Ysim_KS_VOL_50(plotsamp)),'g',...
     plotsamp,log(Ysim_KS_VOL_25(plotsamp)),'r',...
     plotsamp,log(Ysim_KS(plotsamp)),'k',...
     plotsamp,log(Ysim_REITER_VOL_100(plotsamp_perturb)),'b--',...
     plotsamp,log(Ysim_REITER_VOL_50(plotsamp_perturb)),'g--',...
     plotsamp,log(Ysim_REITER_VOL_25(plotsamp_perturb)),'r--',...
     plotsamp,log(Ysim_REITER(plotsamp_perturb)),'k--',...
        'LineWidth',lwidnum); 
title('Output','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
axis([-Inf Inf -0.75 -0.15])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;
legend('100%','50%','25%','Baseline','Location','NorthEast');
legend boxoff;

subplot(2,2,2);
plot(plotsamp,log(Isim_KS_VOL_100(plotsamp)),'b',...
     plotsamp,log(Isim_KS_VOL_50(plotsamp)),'g',...
     plotsamp,log(Isim_KS_VOL_25(plotsamp)),'r',...
     plotsamp,log(Isim_KS(plotsamp)),'k',...
     plotsamp,log(Isim_REITER_VOL_100(plotsamp_perturb)),'b--',...
     plotsamp,log(Isim_REITER_VOL_50(plotsamp_perturb)),'g--',...
     plotsamp,log(Isim_REITER_VOL_25(plotsamp_perturb)),'r--',...
     plotsamp,log(Isim_REITER(plotsamp_perturb)),'k--',...
        'LineWidth',lwidnum); 
title('Investment','FontSize',titlesizenum)
axis([-Inf Inf -3.0 -1.5])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(2,2,3);
plot(plotsamp,log(Nsim_KS_VOL_100(plotsamp)),'b',...
     plotsamp,log(Nsim_KS_VOL_50(plotsamp)),'g',...
     plotsamp,log(Nsim_KS_VOL_25(plotsamp)),'r',...
     plotsamp,log(Nsim_KS(plotsamp)),'k',...
     plotsamp,log(Nsim_REITER_VOL_100(plotsamp_perturb)),'b--',...
     plotsamp,log(Nsim_REITER_VOL_50(plotsamp_perturb)),'g--',...
     plotsamp,log(Nsim_REITER_VOL_25(plotsamp_perturb)),'r--',...
     plotsamp,log(Nsim_REITER(plotsamp_perturb)),'k--',...
    'LineWidth',lwidnum);
title('Labor','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -1.225 -0.975])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(2,2,4);
plot(plotsamp,-log(psim_KS_VOL_100(plotsamp)),'b',...
     plotsamp,-log(psim_KS_VOL_50(plotsamp)),'g',...
     plotsamp,-log(psim_KS_VOL_25(plotsamp)),'r',...
     plotsamp,-log(psim_KS(plotsamp)),'k',...
     plotsamp,-log(psim_REITER_VOL_100(plotsamp_perturb)),'b--',...
     plotsamp,-log(psim_REITER_VOL_50(plotsamp_perturb)),'g--',...
     plotsamp,-log(psim_REITER_VOL_25(plotsamp_perturb)),'r--',...
     plotsamp,-log(psim_REITER(plotsamp_perturb)),'k--',...
    'LineWidth',lwidnum);
title('Consumption','FontSize',titlesizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -0.9 -0.7])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

saveas(gcf,'UNC_KS_VS_REITER.pdf')


%%HP filter data and show vols

%do KS HP filtering
HPdata_KS_VOL_10 = log([Ysim_KS_VOL_10(samp) Isim_KS_VOL_10(samp) Nsim_KS_VOL_10(samp) asim(samp) psim_KS_VOL_10(samp)]);
HPdata_KS_VOL_25 = log([Ysim_KS_VOL_25(samp) Isim_KS_VOL_25(samp) Nsim_KS_VOL_25(samp) asim(samp) psim_KS_VOL_25(samp)]);
HPdata_KS_VOL_33 = log([Ysim_KS_VOL_33(samp) Isim_KS_VOL_33(samp) Nsim_KS_VOL_33(samp) asim(samp) psim_KS_VOL_33(samp)]);
HPdata_KS_VOL_50 = log([Ysim_KS_VOL_50(samp) Isim_KS_VOL_50(samp) Nsim_KS_VOL_50(samp) asim(samp) psim_KS_VOL_50(samp)]);
HPdata_KS_VOL_100 = log([Ysim_KS_VOL_100(samp) Isim_KS_VOL_100(samp) Nsim_KS_VOL_100(samp) asim(samp) psim_KS_VOL_100(samp)]);

[HPtrend_KS_VOL_10,HPdetrend_KS_VOL_10] = hptrend(HPdata_KS_VOL_10,100);
[HPtrend_KS_VOL_25,HPdetrend_KS_VOL_25] = hptrend(HPdata_KS_VOL_25,100);
[HPtrend_KS_VOL_33,HPdetrend_KS_VOL_33] = hptrend(HPdata_KS_VOL_33,100);
[HPtrend_KS_VOL_50,HPdetrend_KS_VOL_50] = hptrend(HPdata_KS_VOL_50,100);
[HPtrend_KS_VOL_100,HPdetrend_KS_VOL_100] = hptrend(HPdata_KS_VOL_100,100);

HP_VOL_TABLE_KS_VOL = zeros(6,5);
HP_VOL_TABLE_KS_VOL(1,:) = HP_VOL_TABLE(1,:);
HP_VOL_TABLE_KS_VOL(2,:) = 100*var(HPdetrend_KS_VOL_10).^0.5; HP_VOL_TABLE_KS_VOL(2,2:end)=HP_VOL_TABLE_KS_VOL(2,2:end)./HP_VOL_TABLE_KS_VOL(2,1);
HP_VOL_TABLE_KS_VOL(3,:) = 100*var(HPdetrend_KS_VOL_25).^0.5; HP_VOL_TABLE_KS_VOL(3,2:end)=HP_VOL_TABLE_KS_VOL(3,2:end)./HP_VOL_TABLE_KS_VOL(3,1);
HP_VOL_TABLE_KS_VOL(4,:) = 100*var(HPdetrend_KS_VOL_33).^0.5; HP_VOL_TABLE_KS_VOL(4,2:end)=HP_VOL_TABLE_KS_VOL(4,2:end)./HP_VOL_TABLE_KS_VOL(4,1);
HP_VOL_TABLE_KS_VOL(5,:) = 100*var(HPdetrend_KS_VOL_50).^0.5; HP_VOL_TABLE_KS_VOL(5,2:end)=HP_VOL_TABLE_KS_VOL(5,2:end)./HP_VOL_TABLE_KS_VOL(5,1);
HP_VOL_TABLE_KS_VOL(6,:) = 100*var(HPdetrend_KS_VOL_100).^0.5; HP_VOL_TABLE_KS_VOL(6,2:end)=HP_VOL_TABLE_KS_VOL(6,2:end)./HP_VOL_TABLE_KS_VOL(6,1);

HP_CORR_TABLE_KS_VOL = zeros(6,5);
HP_CORR_TABLE_KS_VOL(1,:) = HP_CORR_TABLE(1,:);
corrmat = corrcoef(HPdetrend_KS_VOL_10); HP_CORR_TABLE_KS_VOL(2,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_KS_VOL_25); HP_CORR_TABLE_KS_VOL(3,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_KS_VOL_33); HP_CORR_TABLE_KS_VOL(4,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_KS_VOL_50); HP_CORR_TABLE_KS_VOL(5,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_KS_VOL_100); HP_CORR_TABLE_KS_VOL(6,:) = corrmat(1,:);

%do REITER HP filtering
HPdata_REITER_VOL_10 = log([Ysim_REITER_VOL_10(samp_perturb) Isim_REITER_VOL_10(samp_perturb) Nsim_REITER_VOL_10(samp_perturb) asim(samp) psim_REITER_VOL_10(samp_perturb)]);
HPdata_REITER_VOL_25 = log([Ysim_REITER_VOL_25(samp_perturb) Isim_REITER_VOL_25(samp_perturb) Nsim_REITER_VOL_25(samp_perturb) asim(samp) psim_REITER_VOL_25(samp_perturb)]);
HPdata_REITER_VOL_33 = log([Ysim_REITER_VOL_33(samp_perturb) Isim_REITER_VOL_33(samp_perturb) Nsim_REITER_VOL_33(samp_perturb) asim(samp) psim_REITER_VOL_33(samp_perturb)]);
HPdata_REITER_VOL_50 = log([Ysim_REITER_VOL_50(samp_perturb) Isim_REITER_VOL_50(samp_perturb) Nsim_REITER_VOL_50(samp_perturb) asim(samp) psim_REITER_VOL_50(samp_perturb)]);
HPdata_REITER_VOL_100 = log([Ysim_REITER_VOL_100(samp_perturb) Isim_REITER_VOL_100(samp_perturb) Nsim_REITER_VOL_100(samp_perturb) asim(samp) psim_REITER_VOL_100(samp_perturb)]);

[HPtrend_REITER_VOL_10,HPdetrend_REITER_VOL_10] = hptrend(HPdata_REITER_VOL_10,100);
[HPtrend_REITER_VOL_25,HPdetrend_REITER_VOL_25] = hptrend(HPdata_REITER_VOL_25,100);
[HPtrend_REITER_VOL_33,HPdetrend_REITER_VOL_33] = hptrend(HPdata_REITER_VOL_33,100);
[HPtrend_REITER_VOL_50,HPdetrend_REITER_VOL_50] = hptrend(HPdata_REITER_VOL_50,100);
[HPtrend_REITER_VOL_100,HPdetrend_REITER_VOL_100] = hptrend(HPdata_REITER_VOL_100,100);

HP_VOL_TABLE_REITER_VOL = zeros(6,5);
HP_VOL_TABLE_REITER_VOL(1,:) = HP_VOL_TABLE(4,:);
HP_VOL_TABLE_REITER_VOL(2,:) = 100*var(HPdetrend_REITER_VOL_10).^0.5; HP_VOL_TABLE_REITER_VOL(2,2:end)=HP_VOL_TABLE_REITER_VOL(2,2:end)./HP_VOL_TABLE_REITER_VOL(2,1);
HP_VOL_TABLE_REITER_VOL(3,:) = 100*var(HPdetrend_REITER_VOL_25).^0.5; HP_VOL_TABLE_REITER_VOL(3,2:end)=HP_VOL_TABLE_REITER_VOL(3,2:end)./HP_VOL_TABLE_REITER_VOL(3,1);
HP_VOL_TABLE_REITER_VOL(4,:) = 100*var(HPdetrend_REITER_VOL_33).^0.5; HP_VOL_TABLE_REITER_VOL(4,2:end)=HP_VOL_TABLE_REITER_VOL(4,2:end)./HP_VOL_TABLE_REITER_VOL(4,1);
HP_VOL_TABLE_REITER_VOL(5,:) = 100*var(HPdetrend_REITER_VOL_50).^0.5; HP_VOL_TABLE_REITER_VOL(5,2:end)=HP_VOL_TABLE_REITER_VOL(5,2:end)./HP_VOL_TABLE_REITER_VOL(5,1);
HP_VOL_TABLE_REITER_VOL(6,:) = 100*var(HPdetrend_REITER_VOL_100).^0.5; HP_VOL_TABLE_REITER_VOL(6,2:end)=HP_VOL_TABLE_REITER_VOL(6,2:end)./HP_VOL_TABLE_REITER_VOL(6,1);


HP_CORR_TABLE_REITER_VOL = zeros(6,5);
HP_CORR_TABLE_REITER_VOL(1,:) = HP_CORR_TABLE(4,:);
corrmat = corrcoef(HPdetrend_REITER_VOL_10); HP_CORR_TABLE_REITER_VOL(2,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_REITER_VOL_25); HP_CORR_TABLE_REITER_VOL(3,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_REITER_VOL_33); HP_CORR_TABLE_REITER_VOL(4,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_REITER_VOL_50); HP_CORR_TABLE_REITER_VOL(5,:) = corrmat(1,:);
corrmat = corrcoef(HPdetrend_REITER_VOL_100); HP_CORR_TABLE_REITER_VOL(6,:) = corrmat(1,:);


bardata = [HP_VOL_TABLE_KS_VOL([1 3 5 6],1)/HP_VOL_TABLE_KS_VOL(1,1) HP_VOL_TABLE_REITER_VOL([1 3 5 6],1)/HP_VOL_TABLE_REITER_VOL(1,1)];
figure;
bar(bardata,'BarWidth',lwidnum)
axis([-Inf Inf 0.8 1.15])
colormap([0 0 0 ; 0 0 1]);
legend('KS','REITER','Location','NorthWest')
legend boxoff;
title('Output Volatility')
ylabel('Standard Deviation Relative to Baseline','FontSize',labelsizenum)
set(gca,'XTickLabel',{'Baseline','25%','50%','100%'})
set(gca,'FontSize',labelsizenum)
xlabel('Micro-Level Volatility Fluctuations: \gamma_s','FontSize',labelsizenum)
saveas(gcf,'UNC_VOLS.pdf')
%READ IN DATA FROM NEW SEED TO COMPUTE DH STATS, R2'S AND RMSEs


%10
Kbarsim_KS_VOL_10_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/10/NEWSEED/kbarsim.txt');
KbarFCST_KS_VOL_10_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/10/NEWSEED/kbarfcstsim.txt');
KbarDH_KS_VOL_10_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/10/NEWSEED/kbaraggrule.txt');
pFCST_KS_VOL_10_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/10/NEWSEED/pfcstsim.txt');
pDH_KS_VOL_10_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/10/NEWSEED/paggrule.txt');
psim_KS_VOL_10_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/10/NEWSEED/psim.txt');

%25
Kbarsim_KS_VOL_25_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/25/NEWSEED/kbarsim.txt');
KbarFCST_KS_VOL_25_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/25/NEWSEED/kbarfcstsim.txt');
KbarDH_KS_VOL_25_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/25/NEWSEED/kbaraggrule.txt');
pFCST_KS_VOL_25_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/25/NEWSEED/pfcstsim.txt');
pDH_KS_VOL_25_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/25/NEWSEED/paggrule.txt');
psim_KS_VOL_25_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/25/NEWSEED/psim.txt');

%33
Kbarsim_KS_VOL_33_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/33/NEWSEED/kbarsim.txt');
KbarFCST_KS_VOL_33_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/33/NEWSEED/kbarfcstsim.txt');
KbarDH_KS_VOL_33_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/33/NEWSEED/kbaraggrule.txt');
pFCST_KS_VOL_33_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/33/NEWSEED/pfcstsim.txt');
pDH_KS_VOL_33_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/33/NEWSEED/paggrule.txt');
psim_KS_VOL_33_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/33/NEWSEED/psim.txt');

%50
Kbarsim_KS_VOL_50_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/50/NEWSEED/kbarsim.txt');
KbarFCST_KS_VOL_50_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/50/NEWSEED/kbarfcstsim.txt');
KbarDH_KS_VOL_50_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/50/NEWSEED/kbaraggrule.txt');
pFCST_KS_VOL_50_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/50/NEWSEED/pfcstsim.txt');
pDH_KS_VOL_50_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/50/NEWSEED/paggrule.txt');
psim_KS_VOL_50_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/50/NEWSEED/psim.txt');

%100
Kbarsim_KS_VOL_100_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/100/NEWSEED/kbarsim.txt');
KbarFCST_KS_VOL_100_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/100/NEWSEED/kbarfcstsim.txt');
KbarDH_KS_VOL_100_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/100/NEWSEED/kbaraggrule.txt');
pFCST_KS_VOL_100_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/100/NEWSEED/pfcstsim.txt');
pDH_KS_VOL_100_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/100/NEWSEED/paggrule.txt');
psim_KS_VOL_100_NEWSEED = importdata('./IDIO_VOL_SHOCKS/KS/100/NEWSEED/psim.txt');


%%COMPUTE DH STATS
DH_TABLE_NEWSEED_VOL = zeros(6,4);

%copy KS DH stats for vol elasticity = 0%
DH_TABLE_NEWSEED_VOL(1,:) = DH_TABLE_NEWSEED(1,:);

%10 percent
DH_TABLE_NEWSEED_VOL(2,1) = max(100*abs(log(Kbarsim_KS_VOL_10_NEWSEED(samp))-log(KbarDH_KS_VOL_10_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(2,2) = mean(100*abs(log(Kbarsim_KS_VOL_10_NEWSEED(samp))-log(KbarDH_KS_VOL_10_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(2,3) = max(100*abs(log(psim_KS_VOL_10_NEWSEED(samp))-log(pDH_KS_VOL_10_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(2,4) = mean(100*abs(log(psim_KS_VOL_10_NEWSEED(samp))-log(pDH_KS_VOL_10_NEWSEED(samp))));

%25 percent
DH_TABLE_NEWSEED_VOL(3,1) = max(100*abs(log(Kbarsim_KS_VOL_25_NEWSEED(samp))-log(KbarDH_KS_VOL_25_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(3,2) = mean(100*abs(log(Kbarsim_KS_VOL_25_NEWSEED(samp))-log(KbarDH_KS_VOL_25_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(3,3) = max(100*abs(log(psim_KS_VOL_25_NEWSEED(samp))-log(pDH_KS_VOL_25_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(3,4) = mean(100*abs(log(psim_KS_VOL_25_NEWSEED(samp))-log(pDH_KS_VOL_25_NEWSEED(samp))));

%33 percent
DH_TABLE_NEWSEED_VOL(4,1) = max(100*abs(log(Kbarsim_KS_VOL_33_NEWSEED(samp))-log(KbarDH_KS_VOL_33_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(4,2) = mean(100*abs(log(Kbarsim_KS_VOL_33_NEWSEED(samp))-log(KbarDH_KS_VOL_33_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(4,3) = max(100*abs(log(psim_KS_VOL_33_NEWSEED(samp))-log(pDH_KS_VOL_33_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(4,4) = mean(100*abs(log(psim_KS_VOL_33_NEWSEED(samp))-log(pDH_KS_VOL_33_NEWSEED(samp))));

%50 percent
DH_TABLE_NEWSEED_VOL(5,1) = max(100*abs(log(Kbarsim_KS_VOL_50_NEWSEED(samp))-log(KbarDH_KS_VOL_50_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(5,2) = mean(100*abs(log(Kbarsim_KS_VOL_50_NEWSEED(samp))-log(KbarDH_KS_VOL_50_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(5,3) = max(100*abs(log(psim_KS_VOL_50_NEWSEED(samp))-log(pDH_KS_VOL_50_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(5,4) = mean(100*abs(log(psim_KS_VOL_50_NEWSEED(samp))-log(pDH_KS_VOL_50_NEWSEED(samp))));

%100 percent
DH_TABLE_NEWSEED_VOL(6,1) = max(100*abs(log(Kbarsim_KS_VOL_100_NEWSEED(samp))-log(KbarDH_KS_VOL_100_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(6,2) = mean(100*abs(log(Kbarsim_KS_VOL_100_NEWSEED(samp))-log(KbarDH_KS_VOL_100_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(6,3) = max(100*abs(log(psim_KS_VOL_100_NEWSEED(samp))-log(pDH_KS_VOL_100_NEWSEED(samp))));
DH_TABLE_NEWSEED_VOL(6,4) = mean(100*abs(log(psim_KS_VOL_100_NEWSEED(samp))-log(pDH_KS_VOL_100_NEWSEED(samp))));

%COMPUTE R^2s & RMSEs
ESS_KS_NEWSEED_VOL_10 = zeros(anum,2);
ESS_KS_NEWSEED_VOL_25 = zeros(anum,2);
ESS_KS_NEWSEED_VOL_33 = zeros(anum,2);
ESS_KS_NEWSEED_VOL_50 = zeros(anum,2);
ESS_KS_NEWSEED_VOL_100 = zeros(anum,2);

TSS_KS_NEWSEED_VOL_10 = zeros(anum,2);
TSS_KS_NEWSEED_VOL_25 = zeros(anum,2);
TSS_KS_NEWSEED_VOL_33 = zeros(anum,2);
TSS_KS_NEWSEED_VOL_50 = zeros(anum,2);
TSS_KS_NEWSEED_VOL_100 = zeros(anum,2);

MEANS_KS_NEWSEED_VOL_10 = zeros(anum,2);
MEANS_KS_NEWSEED_VOL_25 = zeros(anum,2);
MEANS_KS_NEWSEED_VOL_33 = zeros(anum,2);
MEANS_KS_NEWSEED_VOL_50 = zeros(anum,2);
MEANS_KS_NEWSEED_VOL_100 = zeros(anum,2);

%first, compute means
for t = samp(1):samp(end);
    act = asimpos_NEWSEED(t);

    MEANS_KS_NEWSEED_VOL_10(act,:) = MEANS_KS_NEWSEED_VOL_10(act,:) + [log(psim_KS_VOL_10_NEWSEED(t)) log(Kbarsim_KS_VOL_10_NEWSEED(t+1))];
    MEANS_KS_NEWSEED_VOL_25(act,:) = MEANS_KS_NEWSEED_VOL_25(act,:) + [log(psim_KS_VOL_25_NEWSEED(t)) log(Kbarsim_KS_VOL_25_NEWSEED(t+1))];
    MEANS_KS_NEWSEED_VOL_33(act,:) = MEANS_KS_NEWSEED_VOL_33(act,:) + [log(psim_KS_VOL_33_NEWSEED(t)) log(Kbarsim_KS_VOL_33_NEWSEED(t+1))];
    MEANS_KS_NEWSEED_VOL_50(act,:) = MEANS_KS_NEWSEED_VOL_50(act,:) + [log(psim_KS_VOL_50_NEWSEED(t)) log(Kbarsim_KS_VOL_50_NEWSEED(t+1))];
    MEANS_KS_NEWSEED_VOL_100(act,:) = MEANS_KS_NEWSEED_VOL_100(act,:) + [log(psim_KS_VOL_100_NEWSEED(t)) log(Kbarsim_KS_VOL_100_NEWSEED(t+1))];
    
end;
MEANS_KS_NEWSEED_VOL_10 = MEANS_KS_NEWSEED_VOL_10./repmat(COUNTS_NEWSEED,1,2);
MEANS_KS_NEWSEED_VOL_25 = MEANS_KS_NEWSEED_VOL_25./repmat(COUNTS_NEWSEED,1,2);
MEANS_KS_NEWSEED_VOL_33 = MEANS_KS_NEWSEED_VOL_33./repmat(COUNTS_NEWSEED,1,2);
MEANS_KS_NEWSEED_VOL_50 = MEANS_KS_NEWSEED_VOL_50./repmat(COUNTS_NEWSEED,1,2);
MEANS_KS_NEWSEED_VOL_100 = MEANS_KS_NEWSEED_VOL_100./repmat(COUNTS_NEWSEED,1,2);

%then, compute ESS and TSS
for t = samp(1):samp(end);
    act = asimpos_NEWSEED(t);
   
    %ESS
    ESS_KS_NEWSEED_VOL_10(act,:) = ESS_KS_NEWSEED_VOL_10(act,:) + [log(pFCST_KS_VOL_10_NEWSEED(t)/psim_KS_VOL_10_NEWSEED(t)) log(KbarFCST_KS_VOL_10_NEWSEED(t)/Kbarsim_KS_VOL_10_NEWSEED(t))].^2;
    ESS_KS_NEWSEED_VOL_25(act,:) = ESS_KS_NEWSEED_VOL_25(act,:) + [log(pFCST_KS_VOL_25_NEWSEED(t)/psim_KS_VOL_25_NEWSEED(t)) log(KbarFCST_KS_VOL_25_NEWSEED(t)/Kbarsim_KS_VOL_25_NEWSEED(t))].^2;
    ESS_KS_NEWSEED_VOL_33(act,:) = ESS_KS_NEWSEED_VOL_33(act,:) + [log(pFCST_KS_VOL_33_NEWSEED(t)/psim_KS_VOL_33_NEWSEED(t)) log(KbarFCST_KS_VOL_33_NEWSEED(t)/Kbarsim_KS_VOL_33_NEWSEED(t))].^2;
    ESS_KS_NEWSEED_VOL_50(act,:) = ESS_KS_NEWSEED_VOL_50(act,:) + [log(pFCST_KS_VOL_50_NEWSEED(t)/psim_KS_VOL_50_NEWSEED(t)) log(KbarFCST_KS_VOL_50_NEWSEED(t)/Kbarsim_KS_VOL_50_NEWSEED(t))].^2;
    ESS_KS_NEWSEED_VOL_100(act,:) = ESS_KS_NEWSEED_VOL_100(act,:) + [log(pFCST_KS_VOL_100_NEWSEED(t)/psim_KS_VOL_100_NEWSEED(t)) log(KbarFCST_KS_VOL_100_NEWSEED(t)/Kbarsim_KS_VOL_100_NEWSEED(t))].^2;
    
    %TSS
    TSS_KS_NEWSEED_VOL_10(act,:) = TSS_KS_NEWSEED_VOL_10(act,:) + [log(MEANS_KS_NEWSEED_VOL_10(act,1)/psim_KS_VOL_10_NEWSEED(t)) log(MEANS_KS_NEWSEED_VOL_10(act,2)/Kbarsim_KS_VOL_10_NEWSEED(t))].^2;
    TSS_KS_NEWSEED_VOL_25(act,:) = TSS_KS_NEWSEED_VOL_25(act,:) + [log(MEANS_KS_NEWSEED_VOL_25(act,1)/psim_KS_VOL_25_NEWSEED(t)) log(MEANS_KS_NEWSEED_VOL_25(act,2)/Kbarsim_KS_VOL_25_NEWSEED(t))].^2;
    TSS_KS_NEWSEED_VOL_33(act,:) = TSS_KS_NEWSEED_VOL_33(act,:) + [log(MEANS_KS_NEWSEED_VOL_33(act,1)/psim_KS_VOL_33_NEWSEED(t)) log(MEANS_KS_NEWSEED_VOL_33(act,2)/Kbarsim_KS_VOL_33_NEWSEED(t))].^2;
    TSS_KS_NEWSEED_VOL_50(act,:) = TSS_KS_NEWSEED_VOL_50(act,:) + [log(MEANS_KS_NEWSEED_VOL_50(act,1)/psim_KS_VOL_50_NEWSEED(t)) log(MEANS_KS_NEWSEED_VOL_50(act,2)/Kbarsim_KS_VOL_50_NEWSEED(t))].^2;
    TSS_KS_NEWSEED_VOL_100(act,:) = TSS_KS_NEWSEED_VOL_100(act,:) + [log(MEANS_KS_NEWSEED_VOL_100(act,1)/psim_KS_VOL_100_NEWSEED(t)) log(MEANS_KS_NEWSEED_VOL_100(act,2)/Kbarsim_KS_VOL_100_NEWSEED(t))].^2;
    
end;

%now, actually implement the R^2 and RMSE formulas
R2_KS_NEWSEED_VOL_10 = 1 - ESS_KS_NEWSEED_VOL_10./TSS_KS_NEWSEED_VOL_10;
R2_KS_NEWSEED_VOL_25 = 1 - ESS_KS_NEWSEED_VOL_25./TSS_KS_NEWSEED_VOL_25;
R2_KS_NEWSEED_VOL_33 = 1 - ESS_KS_NEWSEED_VOL_33./TSS_KS_NEWSEED_VOL_33;
R2_KS_NEWSEED_VOL_50 = 1 - ESS_KS_NEWSEED_VOL_50./TSS_KS_NEWSEED_VOL_50;
R2_KS_NEWSEED_VOL_100 = 1 - ESS_KS_NEWSEED_VOL_100./TSS_KS_NEWSEED_VOL_100;

RMSE_KS_NEWSEED_VOL_10 = 100*sqrt(ESS_KS_NEWSEED_VOL_10./repmat(COUNTS_NEWSEED,1,2));
RMSE_KS_NEWSEED_VOL_25 = 100*sqrt(ESS_KS_NEWSEED_VOL_25./repmat(COUNTS_NEWSEED,1,2));
RMSE_KS_NEWSEED_VOL_33 = 100*sqrt(ESS_KS_NEWSEED_VOL_33./repmat(COUNTS_NEWSEED,1,2));
RMSE_KS_NEWSEED_VOL_50 = 100*sqrt(ESS_KS_NEWSEED_VOL_50./repmat(COUNTS_NEWSEED,1,2));
RMSE_KS_NEWSEED_VOL_100 = 100*sqrt(ESS_KS_NEWSEED_VOL_100./repmat(COUNTS_NEWSEED,1,2));

%%%%%%%%%%%%%%%%%%
%%%%%READ IN & PROCESS CONT SHOCK DATA
%%%%%%%%%%%%%%%%%%

%KS cont shocks with discrete sim
Ysim_KS_CONT = importdata('./KS/CONT/ysim.txt');
Isim_KS_CONT = importdata('./KS/CONT/isim.txt');
Nsim_KS_CONT = importdata('./KS/CONT/Nsim.txt');
psim_KS_CONT = importdata('./KS/CONT/psim.txt');
Kbarsim_KS_CONT = importdata('./KS/CONT/kbarsim.txt');

%REITER cont shocks 
Ysim_REITER_CONT = importdata('./REITER/CONT/ysim.txt');
Isim_REITER_CONT = importdata('./REITER/CONT/isim.txt');
Nsim_REITER_CONT = importdata('./REITER/CONT/Nsim.txt');
psim_REITER_CONT = importdata('./REITER/CONT/psim.txt');

asim_CONT = importdata('./KS/CONT/asim.txt');
plotsamp = 900:950;
plotsamp_perturb = plotsamp+1;
%make unconditional simulation figure

figure;
subplot(2,2,1);
plot(plotsamp,log(Ysim_REITER_CONT(plotsamp_perturb)),'b',...
    plotsamp,log(Ysim_KS_CONT(plotsamp)),'k',...
        'LineWidth',lwidnum); 
title('Output','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
axis([-Inf Inf -0.75 -0.3])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;
legend('REITER','KS','Location','NorthEast');
legend boxoff;

subplot(2,2,2);
plot(plotsamp,log(Isim_REITER_CONT(plotsamp_perturb)),'b',...
    plotsamp,log(Isim_KS_CONT(plotsamp)),'k',...
    'LineWidth',lwidnum);
title('Investment','FontSize',titlesizenum)
axis([-Inf Inf -3.0 -1.6])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(2,2,3);
plot(plotsamp,log(Nsim_REITER_CONT(plotsamp_perturb)),'b',...
    plotsamp,log(Nsim_KS_CONT(plotsamp)),'k',...
    'LineWidth',lwidnum);
title('Labor','FontSize',titlesizenum)
ylabel('Log','FontSize',labelsizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -1.225 -0.975])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

subplot(2,2,4);
plot(plotsamp,-log(psim_REITER_CONT(plotsamp_perturb)),'b',...
    plotsamp,-log(psim_KS_CONT(plotsamp)),'k',...
    'LineWidth',lwidnum);
title('Consumption','FontSize',titlesizenum)
xlabel('Year','FontSize',labelsizenum)
axis([-Inf Inf -0.9 -0.7])
ax=gca;
ax.XTick = [plotsamp(1) plotsamp(10) plotsamp(20) plotsamp(30) plotsamp(40) plotsamp(50)];
ax.XTickLabel = {'1' '10' '20' '30' '40' '50'};
ax.FontSize = labelsizenum;

saveas(gcf,'UNCOND_CONT.pdf')

