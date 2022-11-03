% load Building Automation System (BAS) model
% author: Birgit van Huijgevoort

%% This code loads a class object for an LTI model of the BAS tutorial
% 
% Original 7D model:
% -------------------------------------------------------
% x_c[k+1] = A_cx_c[k] + B_cu_c[k] +F_cd_c[k] + Q_c
% y_c[k]   = [1 0 0 0 0 0 0]
% x_c = [T_z1 T_z2 T_w5 T_w6 T_w2 T_w3 T_w7]^T
% u_c = T_sa
% d_c =[T_out T_hall CO2_1 CO2_2 T_rw,r1 T_rw,r2]^T
% -------------------------------------------------------
% author of the 7d model: Nathalie Cauchi
%
% the model has 7 dimensions, 1 output and 1 input
% the input should be designed to
% achieve an objective.
%
% Note that this is an affine term! We will load the LTI part in sysLTI and
% give the affine term 'a' seperately.

% Load parameters needed to build model
Model = load('BASmodel.mat');
Z1m = Model.Z1m;
a = Model.Qc;

% Define system parameters
A = Z1m.A;
B = Z1m.B;
C = Z1m.C;
D = zeros(1,1);
Bw = Z1m.F(:,1:6);
dim = Z1m.dim;

% Specify mean and variance of disturbance w(t)
mu = [9;15;500;500;35;35]; % mean of disturbance
sigma = [1, zeros(1,5); 0, 1, zeros(1,4); 0, 0 , 100, zeros(1,3);
                0,0,0,100,zeros(1,2); zeros(1,4), 5, 0;
                zeros(1,5), 5];% variance of disturbance

% Set up an LTI model
sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);