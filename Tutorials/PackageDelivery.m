%%% Package delivery case study
%
% Package delivery is a 2D LTI systems of the form
% x(t+1) = Ax(t) + Bu(t) + Bw w(t)
% y(t) = Cx(t) + Du(t)
% Be aware: Current implementation requires D = 0!
%
% Expected runtime = approx 10 seconds
% 
% Authors: Oliver Schön, Birgit van Huijgevoort

clc; clear; close all;

% Add toolboxes to path
run Install.m

% Track runtime
tStart = tic;
disp('Start package delivery benchmark')

%% Specify system and regions

% Define system parameters
A = 0.9*eye(2);
B = eye(2);
C = eye(2);
D = zeros(2);
Bw = sqrt(0.2)*eye(2);
dim = length(A);

% Specify mean and variance of disturbance w(t) 
mu = zeros(dim, 1); % Mean of disturbance
sigma = eye(dim); % Variance of disturbance;

% Set up an LTI model
sysLTI = LinModel(A, B, C, D, Bw, mu, sigma);
 
% Bounds on state space 
x1l = -6; % Lowerbound in x1
x1u = 6; % Upperbound in x1
x2l = -6; % Lowerbound in x2
x2u = 6; % Upperbound in x2
% Define bounded state space
sysLTI.X = Polyhedron(combvec([x1l, x1u], [x2l, x2u])');

% Bounds on  and input space
ul = [-1; -1]; % Lowerbound input u
uu = [1; 1]; % Upperbound input u
% Define bounded input space
sysLTI.U = Polyhedron(combvec([ul(1), uu(1)], [ul(2), uu(2)])');

% Specify regions for the specification
% Pick up a parcel at P1 and deliver it to P3. If on this path the agent passes 
% P2, it loses the package and has to pick up a new one (at P1).

% 1) Pick-up region
p1x = [5 5 6 6 5]; % x1-coordinates
p1y = [-1 1 1 -1 -1]; % x2-coordinates
P1 = Polyhedron([p1x; p1y]');

% 2) Region where you lose package
p2x = [0 0 1 1 0]; % x1-coordinates
p2y = [-5 1 1 -5 -5]; % x2-coordinates
P2 = Polyhedron([p2x; p2y]');

% 3) Delivery region
p3x = [-4 -4 -2 -2 -4];    % x1-coordinates
p3y = [-3 -4 -4 -3 -3];    % x2-coordinates  
P3 = Polyhedron([p3x; p3y]');

sysLTI.regions = [P1; P2; P3]; % regions that get specific atomic propositions
sysLTI.AP = {'p1', 'p2', 'p3'}; % with the corresponding atomic propositions

%Plot_sysLTI(sysLTI)
%% Step 1 Translate the specification (or input DFA yourself)
t1start = tic;

% Define the scLTL formula using syntax as in Specification/LTL2BA/README
formula = 'F(p1 & (!p2 U p3))';

% Translate the spec to a DFA
[DFA] = TranslateSpec(formula,sysLTI.AP);

% or specify DFA directly
% act = {' ', 'p3', 'p2', 'p2p3', 'p1', 'p1p3', 'p1p2', 'p1p2p3'};
% 
% trans = [ ...
%     1, 1, 1, 1, 2, 2, 2, 2; ...
%     2, 3, 1, 1, 2, 3, 1, 1; ...
%     0, 0, 0, 0, 0, 0, 0, 0 ...
%     ];
% 
% DFA = struct( ...
%     'S', [1 2 3], ...
%     'S0', 1, ...
%     'F', 3, ...
%     'act', {act}, ...
%     'trans', trans, ...
%     'sink', find([0, 0]) ... % Empty double row vector
%     );

t1end = toc(t1start);
%% Step 2 Finite-state abstraction
t2start = tic;

% Specify granularity of abstraction/gridding
lu = 3; % Division of input space
lx1 = 400; % Division of state space in x1-direction
lx2 = 400;% Division of state space in x2-direction

% Construct abstract input space uhat
uhat = GridInputSpace(lu,sysLTI.U); % abstract input space

% Construct finite-state abstraction
l = [lx1, lx2]; % Number of grid cells in x1- and x2-direction
tol = 10^(-6);
sysAbs = FSabstraction(sysLTI,uhat,l,tol,DFA,'TensorComputation',true);

t2end = toc(t2start);
%% Step 3 Similarity quantification
t3start = tic;

% Choose a value for epsilon
epsilon = 0.075;
 
% Quantify similarity
rel = QuantifySim(sysLTI, sysAbs, epsilon);

t3end = toc(t3start);
%% Step 4 Synthesize a robust controller
t4start = tic;

thold = 1e-6;     % threshold

% Synthesize an abstract robust controller
[satProb, pol] = SynthesizeRobustController(sysAbs, DFA, rel, thold, false);

t4end = toc(t4start);
%% Step 5 Control refinement
t5start = tic;

% Refine abstract controller to a continous-state controller
Controller = RefineController(satProb,pol,sysAbs,rel,sysLTI,DFA);

t5end = toc(t5start);
%% Step 6 Deployment
t6start = tic;

x0 = [-5;-5];
N = 60;

% Simulate controlled system
[xsim, qsim] = ImplementController(x0, N, Controller);

% Show results
%plotTrajectories(xsim, [x1l x1u; x2l x2u], sysLTI, 'DFAinfo', qsim);

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

disp('Finished package delivery benchmark')