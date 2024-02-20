%%% Building Automation System (BAS) case study
% Including model order reduction AND KK filtering
% 7D model reduced to a 2D model using model-order reduction
% author: Birgit van Huijgevoort
% Modified: Maico Engelaar
%
% Expected runtime = approx 75 seconds
%
% sysLTI -> original model
% sysLTI_KF -> model with reduced order disturbance
% sysLTIr -> reduced order model

clear;close all;
seed = 4;
rng(seed);

% run Install.m

tStart = tic;
disp('Start building automation system benchmark')

%% Specify system parameters of 7-dimensional model and regions

% Load model into sysLTI, with affine-term a
BASmodel

% Define input bounds
ul = 15;
uu = 30;

% Define state space bounds 
x1l = 19.5; % Lowerbound x1
x1u = 20.5; % Upperbound x1
x2l = 10;   % Lowerbound x2
x2u = 20;   % Upperbound x2
LowerBounds = [x1l;x2l;zeros(5,1)];
UpperBounds = [x1u;x2u;30*ones(5,1)];

% Define Initial distribution of x(0)
mu0 = zeros(sysLTI.dim,1);
Sigma0 = diag([0.1,0.1,3,3,3,3,4]);
InitState = {mu0,Sigma0};
sysLTI.InitState = InitState;

% Transform model, such that w comes from Gaussian distribution with mean 0
% and variance identity
[sysLTI, a] = NormalizeDisturbance(sysLTI,a);

% Since the original dynamics are of an affine system, we will transform
% the system to an LTI system by adjusting the state and output wrt a
% steady state value xss.
% In other words, we remove affine term a by looking at a
% steady state solution xss and obtain a new system with state x-xss and
% output y-C*xss

% Choose steady state values
uss = ul;
xss = -inv(A-eye(dim))*(B*uss+a);

% New bounds on input space (u-uss)
ul = ul-uss;     % Lowerbound input u
uu = uu-uss;     % Upperbound input u

% Define bounded input space
sysLTI.U = Polyhedron(combvec([ul(1),uu(1)])');

% Compute new bounds on state space (x-xss)
x1l = x1l-xss(1);   
x1u = x1u-xss(1);  
x2l = x2l-xss(2);   
x2u = x2u-xss(2);   
% Define bounded state space
sysLTI.X = Polyhedron('lb', LowerBounds-xss, 'ub', UpperBounds-xss);

% Specify region for the specification
P1 = Polyhedron(combvec([x1l,x1u])');

% Region that get specific atomic proposition
sysLTI.regions = [P1];
% Proposition corresponding to the region
sysLTI.AP = {'p1'};

% Select interface function for MOR, u = ur + Qxr + K(x-Pxr)
int_f = 1;
% Interface function between reduced-order model and finite-state
% abstraction is automatically set to default, ur = uhat

%% Step 1 Translate the specification for the original model
t1start = tic;

N = 6; % (finite) time horizon of spec

% Define the scLTL specification
formula = '(p1 & X p1 & X X p1 & X X X p1 & X X X X p1 & X X X X X p1)';

% Translate the spec to a DFA
[DFA] = TranslateSpec(formula, sysLTI.AP);

t1end = toc(t1start);

%% Step 2a Reduce the disturbance on the original model            
t2start = tic;
% Define tuning variables (Make sure that NC=H)
% with H = sysLTI.C = [1 0 0 0 0 0 0]
Cobs = [1 0 0 0 0 0 0; 0 1 0 0 0 0 0]; %N = [1 0];

sysLTI_KF = KKfilter(sysLTI,Cobs);

% Transform model, such that w comes from Gaussian distribution with mean 0
% and variance identity
[sysLTI_KF, atest] = NormalizeDisturbance(sysLTI_KF);

if atest ~=0
    fprintf('Abort program\n');
end

%% Step 2b Model order reduction
f = 0.098;  % tuning parameter for feedback-matrix F
dimr = 2; % desired dimension of reduced-order model
% Construct reduced-order model
[sysLTIr, ~] = ModelReduction(sysLTI_KF,dimr,f);

% Compute projection matrix P and Q for interface function
% u = ur + Qxr + K(x-Pxr) and add them to sysLTIr
sysLTIr = ComputeProjection(sysLTI_KF,sysLTIr);

% Define bounded state space of reduced-order model
sysLTIr.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');
% Define bounded input space of reduced-order model
sysLTIr.U = sysLTI_KF.U;

% Regions that get specific atomic propositions
sysLTIr.regions = [P1];

% Propositions corresponding to the regions
sysLTIr.AP = {'p1'};

%% Step 2 Finite-state abstraction

% Construct abstract input space
lu = 3;  % number of abstract inputs in each direction
[uhat,sysLTIr.U] = GridInputSpace(lu,sysLTIr.U,'interface',int_f,0.6,0.175); % abstract input space

% Reduce the state space to speed up computations [option only available for invariance specs]
[sysLTIr,~] = ReduceX(sysLTIr, sysLTIr.U{2}, P1, 'invariance', 5);

% Construct finite-state abstraction
l = [3000*3000];  % total number of grid cells
tol=10^-6;
sysAbs = FSabstraction(sysLTIr,uhat,l,tol,DFA,'TensorComputation',true);

t2end = toc(t2start);
%% Step 3 Similarity quantification
t3start = tic;

% set values of epsilon 
epsilon_1 = 0.2413; % output deviation for MOR simulation relation
epsilon_2 = 0.1087; % output deviation for gridding simulation relation

% Compute MOR simulation relation
[rel_1, K, kernel] = QuantifySim(sysLTI_KF, sysLTIr, epsilon_1, 'MOR', sysAbs);

% Compute gridding simulation relation
[rel_2] = QuantifySim(sysLTIr, sysAbs, epsilon_2);

% Combine simulation relations
rel = CombineSimRel(rel_1, rel_2, sysLTIr, sysAbs);

t3end = toc(t3start);
%% Step 4 Synthesize a robust controller
t4start = tic;

thold = 1e-6;
% Synthesize an abstract robust controller
[satProb, pol] = SynthesizeRobustController(sysAbs,DFA,rel,thold,false);

t4end = toc(t4start);
%% Step 5 Control refinement
t5start = tic;

% Refine abstract controller to a continous-state controller
Controller = RefineController(satProb,pol,sysAbs,rel,sysLTIr,DFA,int_f,K, 'KKfilter',sysLTI_KF);

t5end = toc(t5start);
%% Step 6 Deployment
t6start = tic;

% Supply an initial state for the LTI system (x_LTI = x_affine - xss)
x0 = mvnrnd(sysLTI.InitState{1}, sysLTI.InitState{2}, 1)';

% Simulate controlled system Ns times
Ns = 1;
xsim = ImplementController(x0, N, Controller, Ns, 'MOR', sysLTIr, kernel, 'KKfilter', sysLTI_KF);

% Plot resulting trajectories
%plotTrajectories(xsim, [LowerBounds UpperBounds], sysLTI_KF, 'shift', xss);

t6end = toc(t6start);
%% Show details on computation time and memory usage
tEnd = toc(tStart);

% Display computation time per step and total.
disp(' ')
disp(['Step 1: ', mat2str(t1end,3), ' seconds'])
disp(['Step 2: ', mat2str(t2end,3), ' seconds'])
disp(['Step 3: ', mat2str(t3end,3), ' seconds'])
disp(['Step 4: ', mat2str(t4end,3), ' seconds'])
disp(['Step 5: ', mat2str(t5end,3), ' seconds'])
disp(['Step 6: ', mat2str(t6end,3), ' seconds'])
disp(' ')
disp(['Total runtime = ', mat2str(tEnd,3), ' seconds']) 
disp(' ')
% Display total memory usage
d = whos();
Mem = sum([d.bytes]);
Mem = Mem*1e-6; % memory usage in Mb

disp(['Memory usage: ', mat2str(Mem), ' Mb'])

%% Show satisfaction probability plot 

% Plot satisfaction probability
plotSatProb(satProb, sysAbs, 'initial', DFA, 'shift', xss, 'MOR');
%%% This plot shows that satifaction probability of the reduced-order model 
% compensated with a steady shift of [xss(1);xss(2)]
% to translate the state in this Figure [xr1;xr2] to the corresponding state of the
% affine system, use x0 = sysLTIr.P*[xr1-xss(1);xr2-xss(2)]+xss

disp('Finished building automation system benchmark')