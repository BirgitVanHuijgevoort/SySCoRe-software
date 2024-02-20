%%% 2D nonlinear case study of a van der Pol oscillator 
% of the form
% x(t+1) = f(x,u) + Bw w(t)
% y(t) = Cx(t)
% In this case the input enters in an affine fashion. The full dynamics are
% given in Vanderpol.m in the folder Case study
%
% Be aware: expected runtime = approx 60 minutes (depending on the number
% of cores in your computer)

clear all
clc
close all

run Install.m

tStart = tic;
disp('Start Van der Pol Oscillator benchmark')

%% Specify system parameters and regions

% Load model into sysNonLin
Vanderpol

% Bounds on state space 
x1l = -3;   % Lowerbound x1
x1u = 3;   % Upperbound x1
x2l = -3;   % Lowerbound x2
x2u = 3;   % Upperbound x2
sysNonLin.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');
% Bounds on input space
ul = [-1];   % Lowerbound input u
uu = [1];     % Upperbound input u
sysNonLin.U = Polyhedron([ul(1),uu(1)]');

% Specify regions for the specification
%%% Stay in bounded region
P1 = sysNonLin.X;
%%% Reach in bounded region
P2 = Polyhedron(combvec([-1.4,-0.7],[-2.9,-2])');

sysNonLin.regions = [P1;P2]; % regions that get specific atomic propositions
sysNonLin.AP = {'p1', 'p2'}; % with the corresponding atomic propositions

% Select an interface function:
%%% 0. (default)    u = uhat
%%% 1.              u = uhat + K(x-xhat)
int_f = 1;

%% Step 1 Translate the specification
t1start = tic;

% Define scLTL formula 
formula = '(p1 U p2)';  % p1 = safe region, p2 = target region
 
% Translate formula to deterministic finite automaton 
[DFA] = TranslateSpec(formula,sysNonLin.AP);

t1end = toc(t1start);
%% Step 2 Finite-state abstraction
t2start = tic;

%% Abstraction part 1: PWA approximation   
% Define SQUARE partitions
Np = [41 41]; % number of grid points in each direction (first state space, then if desired input space)

% Perform PWA approximation and quantify difference between original model
% and PWA approximation
[sysPWA] = PWAapproximation(sysNonLin,Np);

%% Abstraction part 2: finite-state abstraction

% Construct abstract input space
lu = 3; % number of abstract inputs in each direction
[uhat,sysPWA.U] = GridInputSpace(lu,sysPWA.U,'interface',int_f,0.6,0.4, 'order'); % abstract input space

% Construct finite-state abstraction
l = [600, 600];  % number of grid cells in x1- and x2-direction
tol = 10^-8;
sysAbs = FSabstraction(sysPWA,uhat,l,tol,DFA,'TensorComputation', true);


t2end = toc(t2start);
%% Step 3 Similarity quantification
t3start = tic;

% Set a value for epsilon
epsilon = 0.1;

% Compute a suitable weighting matrix D for simulation relation 
% ||x-xhat||_D \leq \epsilon
% States that are taken into account when computing D matrix
%%% For the vdPol oscillator, points on the limit cycle are interesting.
States = [1/8*x1l, 6/10*x2u; 5/7*x1u, 5/17*x2u; 2/13*x1u, 5/9*x2l; 3/4*x1l, 1/7*x2l; 0,0]';
[D, ~] = ComputeD(epsilon,sysPWA,sysAbs,'interface',int_f,'states',States);

% Quantify similarity
[rel, sysPWA] = QuantifySim(sysPWA, sysAbs, epsilon, 'interface', int_f, 'weighting', D);


t3end = toc(t3start);
%% Step 4 Synthesize a robust controller
t4start = tic;

% Specify threshold
thold = 1e-6;

% Synthesize an abstract robust controller
[satProb,pol] = SynthesizeRobustController(sysAbs, DFA, rel, thold, false);

t4end = toc(t4start);
%% Step 5 Control refinement
t5start = tic;

% Refine abstract controller to a continous-state controller
Controller = RefineController(satProb,pol,sysAbs,rel,sysPWA,DFA,int_f);

t5end = toc(t5start);
%% Step 6 Deployment
t6start = tic;

x0 = [-1;1]; % initial state
N = 60;     % time horizon

% Simulate controlled system
xsim = ImplementController(x0,N,Controller);
%plotTrajectories(xsim, [x1l,x1u;x2l,x2u], sysPWA);

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

%% Plot satisfaction probability
plotSatProb(satProb, sysAbs, 'initial', DFA);

disp('Finished Van der Pol Oscillator benchmark')
