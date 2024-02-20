%%% Running example, 2D car park case study
% 
% 2D car park is an LTI systems of the form
% x(t+1) = Ax(t) + Bu(t) + Bw w(t)
% y(t) = Cx(t) + Du(t)
%
% Expected runtime = approx 10 seconds

clc
clear
close all

tStart = tic;
disp('Start car park 2D (running example) benchmark')

%% Specify system parameters and regions
tstart = tic;

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
 
% Define bounded state space
sysLTI.X = Polyhedron(combvec([-10,10],[-10,10])');
% Define bounded input space
sysLTI.U = Polyhedron(combvec([-1,1],[-1,1])');

% Specify regions for the specification
P1 = Polyhedron([4, -3.25; 4, 0; 10, 0; 10 -3.25]); % parking region
P2 = Polyhedron([4, 0; 4, 3.25; 10, 3.25; 10 0]);  % avoid region

% Regions that get specific atomic propositions
sysLTI.regions = [P1;P2]; 
% Propositions corresponding to the regions
sysLTI.AP = {'p1', 'p2'}; 

%% Step 1 Translate the specification
t1start = tic;

% Define the scLTL specification
formula = '(!p2 U p1)';  

% Translate the spec to a DFA
[DFA] = TranslateSpec(formula,sysLTI.AP);

t1end = toc(t1start);
%% Step 2 Finite-state abstraction
t2start = tic;

% Construct abstract input space uhat
lu = 7;  % number of abstract inputs in each direction
uhat = GridInputSpace(lu,sysLTI.U); 

% Construct finite-state abstraction
l = [200, 200];  % number of grid cells 
tol=10^-6;  
sysAbs = FSabstraction(sysLTI,uhat,l,tol,DFA,'TensorComputation',true);

t2end = toc(t2start);
%% Step 3 Similarity quantification
t3start = tic;

%[epsilonBounds] = ComputeEpsilonBounds(sysLTI,mu,sigma,sysAbs.beta)

% Choose a value for epsilon
epsilon = 1.005; % interesting values: 0.142, 1.005, 1.4143

% Quantify similarity
simRel = QuantifySim(sysLTI, sysAbs, epsilon);

t3end = toc(t3start);
%% Step 4 Synthesize a robust controller
t4start = tic;

% Specify threshold for convergence error
thold = 1e-6;

% Synthesize an abstract robust controller
[satProb,pol] = SynthesizeRobustController(sysAbs, DFA, simRel, thold, true);

t4end = toc(t4start);
%% Step 5 Control refinement
t5start = tic;

% Refine abstract controller to a continous-state controller
Controller = RefineController(satProb,pol,sysAbs,simRel,sysLTI,DFA);

% Plot satisfaction probability
plotSatProb(satProb, sysAbs, 'initial', DFA);

t5end = toc(t5start);
%% Step 6 Deployment
t6start = tic;

x0 = [-4;-5]; % initial state
N = 40;     % time horizon

% Simulate controlled system
xsim = ImplementController(x0,N,Controller, 10);
plotTrajectories(xsim, [-10, 10; -10, 10], sysLTI);

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

disp('Finished car park 2D (running example) benchmark')