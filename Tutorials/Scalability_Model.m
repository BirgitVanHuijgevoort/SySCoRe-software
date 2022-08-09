clear all
close all

%% Specify system parameters and regions
% LTI systems of the form
% x(t+1) = Ax(t) + Bu(t) + Bw w(t)
% y(t) = Cx(t) + Du(t)

% define parameters
Ns = .1;
n = 1:3;
Ns_sequence = [1 Ns.^(n)./factorial(n)];
c_sequence = Ns_sequence.*0;
c_sequence(1)=1;

% define dynamics
dim =length(Ns_sequence)-1;
A = toeplitz(c_sequence(1:end-1),Ns_sequence(1:end-1));
B = Ns_sequence(end-1:-1:1)';
C = eye(dim);
D = zeros(dim,1);
Bw = .1*eye(dim);

% Specify mean and variance of disturbance w(t)
mu = zeros(dim,1);
sigma = eye(dim);

% save all system parameters (incl Bw) into a struct
sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);

%% Spaces and sets

% Take state space equal to safe set
xil = -10;   % Lowerbound x1
xiu = 10;   % Upperbound x1
bounds = repmat({[xil,xiu]}, n(end),1);
sysLTI.X = Polyhedron(combvec(bounds{:})');

% single input 
ul = [-1];   % Lowerbound input u
uu = [1];     % Upperbound input u
sysLTI.U = Polyhedron([ul(1),uu(1)]');

% Specify regions for the specification
% (using same format as used in function inpolygon)
% p1 = combvec(bounds{:});
% P1 = Polyhedron(p1');

pil = -8;   % Lowerbound x1
piu = 8;   % Upperbound x1
boundsp = repmat({[pil,piu]}, n(end),1);
P2 = Polyhedron(combvec(boundsp{:})');

% Definte the regions and the atomic propositions
sysLTI.regions = [P2];
sysLTI.AP = {'p2'};

%% Synthesize scLTL formula  
%%% use LTL2BA and check if determinstic and accepting state with loop with 1.
% input: (sc)LTL formula and atomic propositions (see readme in folder
% LTL2BA)
% output: struct DFA containing (among other) the transitions
N=5; % Horizon of specification

formula = '( X X p2) ';

[DFA] = TranslateSpec(formula, sysLTI.AP);

%% Construct abstract model
disp('start gridding')
ula = 1*ul;   % part of input for actuation (lowerbound)
uua = 1*uu;
ulf = ul-ula;   % part of input for feedback (lowerbound)
uuf = uu-uua;

lu = 20;
uhat = combvec(linspace(ula(1),uua(1),lu));

l = [50,50,50];  % number of grid cells in x1- and x2-direction
tol=10^-4;


% sysAbs = GridSpace_nonlin_tensor(sysLTI,uhat,l,tol);
sysAbs = GridSpace_nonlin_tensor_v2(sysLTI,uhat,l,tol, 'tensortoolbox');


disp(['---> finished gridding in ', num2str(toc), ' seconds.'])
%% Compute delta based on epsilon
disp('start computing eps delta');tic;

% TODO: epsilons delat doesnt work

epsilon = 3;    % should be larger than vector beta!

[delta, D_m, K] = ComputeDelta(epsilon,sysLTI,sysLTI.mu,sysLTI.sigma,sysAbs.beta);

disp(['delta = ', num2str(delta), ', epsilon = ', num2str(epsilon) ])
rel = SimRel(epsilon,delta,D_m);

disp(['---> finished computing eps delta in ', num2str(toc), ' seconds'] )

%% Synthesize controller
tic

rel.NonDetLabels  = NonDeterministicLabelling(sysAbs.states, sysLTI.regions, rel);

toc

tic
disp('Start computing Robust policy')
N = 5;     % time horizon
[V,pol]  = SynthesizeRobustController(sysAbs,DFA, rel, N, false);
disp(['---> finished computing robust policy in ', num2str(toc), ' seconds.'] )