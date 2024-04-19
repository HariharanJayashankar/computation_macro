%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Econ 702 UMD
% Assignment 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

% fsolve options
options = optimset('Display','off');

% parameters
% model parameters
p.beta = 0.97^(1/12); p.zeta = 1;
p.nu = 1;
p.pistar = 0;
p.epsilon = 5;
p.rho = 0.8;
p.sigma = 0.15;
p.kappa = 0.02;

p.iss = 1/p.beta - 1;

% otehr parameters (numeric mostly)
p.m =  3; % tauchen grid distance
p.na = 50; % number of grids in shock
p.np = 200; % number of price grids
p.gamma = 0.05; % learning rte for equilibrium


% getting shock grid
[aP, agrid] = tauchen(p.m, p.rho, p.sigma, p.na);
p.aP = aP;
p.agrid = agrid;

% getting price grid
pflexprice = flexsoln(1, p);
pss = log(pflexprice((p.na)/2));
p.pgrid = exp(linspace(-3.2*pss, 3.2*pss, p.np)); % log spaced grid
% pss = pflexprice(p.na/2);
% p.pgrid = linspace(0.2*pss, 2*pss, p.np);

agg.Y = 1.5;
w = fsolve(@(x) findEquilibrium(x,  p), 1);
[~,V, polp, pollamb, omegahat, omega, Y, Pi, Ld] = findEquilibrium(w, p);
pdist = sum(omega, 2);
pdistalt = sum(omegahat, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs for pset
%%%%%%%%%%%%%%%%%%%%%%%%%

testgraphs

%% calculating aggregate y manually
