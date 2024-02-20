%%% 2D car park case study with option to set the interface function 
% equal to 0 or 1, with 
% 0. (default)    u = uhat
% 1.              u = uhat + K(x-xhat)
% when the interface is set to 0, this script is the same as
% CarPark2D_RunningExample.m
%
% 2D car park is an LTI systems of the form
% x(t+1) = Ax(t) + Bu(t) + Bw w(t)
% y(t) = Cx(t) + Du(t)
%
% Expected runtime = approx 10 seconds

clc
clear
close all

%% Specify system parameters and regions

% Define system parameters
A = 0.9*eye(2);
B = 0.7*eye(2);
C = eye(2);
D = zeros(2);
Bw = eye(2);
dim = length(A);

% Specify mean and variance of disturbance w(t) 
mu = zeros(dim,1); % mean of disturbance
sigma = eye(dim); % variance of disturbance

% Set up an LTI model
sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);
 
% Bounds on state space 
x1l = -10;   % Lowerbound x1
x1u = 10;   % Upperbound x1
x2l = -10;   % Lowerbound x2
x2u = 10;   % Upperbound x2
sysLTI.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');
% Define bounds on input space
sysLTI.U = Polyhedron(combvec([-1,1],[-1,1])');

% Specify regions for the specification
P1 = Polyhedron([4, -4; 4, 0; 10, 0; 10 -4]); % parking region
P2 = Polyhedron([4, 0; 4, 4; 10, 4; 10 0]);  % avoid region

sysLTI.regions = [P1;P2]; % regions that get specific atomic propositions
sysLTI.AP = {'p1', 'p2'}; % with the corresponding atomic propositions

Plot_sysLTI(sysLTI)

% Select an interface function:
%%% 0. (default)    u = uhat
%%% 1.              u = uhat + K(x-xhat)
int_f = 1;

%% Step 1 Translate the specification

% Define the scLTL specification
formula = '(!p2 U p1)';  

% Translate the spec to a DFA
[DFA] = TranslateSpec(formula,sysLTI.AP);

%% Step 2 Finite-state abstraction
tGridStart = tic;

% Construct abstract input space
lu = 3;  % number of abstract inputs in each direction
[uhat,sysLTI.U] = GridInputSpace(lu,sysLTI.U,'interface',int_f,0.6,0.4); % abstract input space

% Construct finite-state abstraction
l = [200, 200];  % number of grid cells in x1- and x2-direction
tol=10^-6;
sysAbs = FSabstraction(sysLTI,uhat,l,tol,DFA,'TensorComputation',true);

tGridEnd = toc(tGridStart);
%% Step 3 Similarity quantification
tSimStart = tic;

% Choose a value for epsilon
epsilon = 1.005;

% Quantify similarity 
[rel, K] = QuantifySim(sysLTI, sysAbs, epsilon, 'interface', int_f);

tSimEnd = toc(tSimStart);

disp(['delta = ', num2str(rel.delta), ', epsilon = ', num2str(rel.epsilon) ])
%% Step 4 Synthesize a robust controller

% Specify threshold
thold = 1e-6; 

% Synthesize an abstract robust controller
[satProb,pol] = SynthesizeRobustController(sysAbs, DFA, rel, thold, true);

% Plot satisfaction probability
plotSatProb(satProb, sysAbs, 'initial', DFA);
%% Step 5 Control refinement

% Refine abstract controller to a continous-state controller
Controller = RefineController(satProb,pol,sysAbs,rel,sysLTI,DFA,int_f,K);

%% Step 6 Implementation
x0 = [4;8];
N = 40;

% Simulate controlled system
xsim = ImplementController(x0,N,Controller);
plotTrajectories(xsim, [x1l, x1u; x2l, x2u], sysLTI);