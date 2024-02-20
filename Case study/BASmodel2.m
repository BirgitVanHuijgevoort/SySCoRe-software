% load Building Automation System (BAS) model
% author: Birgit van Huijgevoort

%% This code loads a class object for an LTI model of the BAS tutorial
% 
% Original 4D model:
% -------------------------------------------------------
% x_d[k+1] = Ax_d[k] + Bu[k] + Q_d + SigmaW[k]
% y_d[k]   = [1 0 0 0; 0 1 0 0]
% x_d = [T_z1 T_z1 T_rw,rad1 T_rw,rad2]^T
% u   = T_sa
% -------------------------------------------------------
% -------------------------------------------------------
% author of the 4d model: Nathalie Cauchi
% Original file: https://github.com/natchi92/BASBenchmarks/blob/master/src/CaseStudy1_Ms.m
%
% the model has 7 dimensions, 2 outputs and 1 input
% the input should be designed to
% achieve an objective.
%
% Note that this is an affine term! We will load the LTI part in sysLTI and
% give the affine term 'a' seperately.

% Load parameters needed to build model
Model = load('BASModel2.mat');
Z1m = Model.Z1m;
a = Z1m.F;

% Define system parameters
A = Z1m.A;
B = Z1m.B;
C = Z1m.C;
D = zeros(1,1);
Bw = Z1m.sigma; % last column of Z1m.F is the affine term Q_d
dim = Z1m.dim;

% Specify mean and variance of disturbance w(t)
mu = zeros(4,1); % mean of disturbance
sigma = eye(4);% variance of disturbance

% Set up an LTI model
sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);