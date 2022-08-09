clc
clear
close all

tStart = tic;

%% Specify system parameters and regions
% LTI systems of the form
% x(t+1) = Ax(t) + Bu(t) + Bw w(t)
% y(t) = Cx(t) + Du(t)

% Define system parameters
A = 0.9;
B = 0.5;
C = 1;
D = 0;
Bw = sqrt(0.5); % corresponds to w from distr with variance 0.5
dim = length(A);

% Specify mean and variance of disturbance w(t) 
mu = zeros(dim,1); % mean of disturbance
sigma = eye(dim); % variance of disturbance

sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);
 
% Bounds on state space 
x1l = -10;   % Lowerbound x1
x1u = 10;   % Upperbound x1
sysLTI.X = Polyhedron('lb', x1l, 'ub', x1u);

% Bounds on  and input space
ul = [-1];   % Lowerbound input u
uu = [1];     % Upperbound input u
%sysLTI.U = Polyhedron(combvec([ul(1),uu(1)],[ul(2),uu(2)])');
sysLTI.U = [ul,uu];

% Specify regions for the specification
P1 = [5 6];    % x1-coordinates
P1 = Polyhedron('lb', P1(1), 'ub', P1(2));

P2 = [6 10];    % x1-coordinates
P2 = Polyhedron('lb', P2(1), 'ub', P2(2));

sysLTI.regions = [P1;P2];
sysLTI.AP = {'p1', 'p2'}; % with the corresponding atomic propositions

Plot_sysLTI(sysLTI)

%% Synthesize scLTL formula (or input DFA yourself)
%%% use LTL2BA and check if determinstic and accepting state with loop with 1.
% input: (sc)LTL formula and atomic propositions (see readme in folder
% LTL2BA)
% output: struct DFA containing (among other) the transitions

formula = '(!p2 U p1)';  % p1 = parking, p2 = avoid region
% formula should use atomic propositions in sysLTI.AP. 

% Make sure your current folder is the main SySCoRe folder
[DFA] = TranslateSpec(formula,sysLTI.AP);

%% Construct abstract model by gridding it
disp('start gridding');
tGridStart = tic;
% input: sysLTI, sigma, space bounds for input (controller) and state space

% Specify division of input space for actuation and feedback
ula = 1*ul;   % part of input for actuation (lowerbound)
uua = 1*uu;
ulf = ul-ula;   % part of input for feedback (lowerbound)
uuf = uu-uua;

lu = 7;
uhat = combvec(linspace(ula,uua,lu));

l = 200;  % number of grid cells in x1- and x2-direction
tol=10^-6;
sysAbs = Gridding(sysLTI,uhat,l,tol);
sysAbs = Gridding(sysLTI,uhat,l,tol,'TensorComputation', true);

% Save some extra system parameters into struct
sysAbs.orig = sysLTI;

label = zeros(1,prod(l));
[label(1:prod(l))] = deal(1);
inP1 = (sysAbs.states>min(P1.V))& (sysAbs.states<max(P1.V));
inP2 = (sysAbs.states>min(P2.V))& (sysAbs.states<max(P2.V));
[label(inP1)] = deal(3);
[label(inP2)] = deal(2);
sysAbs.labels = label;

tGridEnd = toc(tGridStart);
disp(['---> finished gridding in ', num2str(tGridEnd), ' seconds.'])

%% Quantify similarity
disp('Quantify similarity');
tSimStart = tic;

epsilon = 0.2;

[rel, K] = QuantifySim(sysLTI, sysAbs, epsilon, sysLTI.mu,sysLTI.sigma,sysAbs.beta);

tSimEnd = toc(tSimStart);

disp(['delta = ', num2str(rel.delta), ', epsilon = ', num2str(rel.epsilon) ])
disp(['---> finished quantifying similarity in ', num2str(tSimEnd), ' seconds'] )

%% Synthesize controller
disp('start computing robust controller')

N = 60;     % time horizon
[ satProp,pol] = SynthesizeRobustController(sysAbs,DFA, rel, N, true);

tEnd = toc(tStart);
disp(['Total runtime = ', mat2str(tEnd)])

%% Start simulation
x0 = -2;
xsim = x0;
indexing = 1:length(sysAbs.states);

% find initial abstract state
diff = abs(x0.*ones(size(sysAbs.states))-sysAbs.states);
inR = sqrt(diff.*rel.R.*diff)<=epsilon;
indices_valid = indexing(inR);
[value_max, index_aux] = max(satProp(inR)); 
j = indices_valid(index_aux); % find maximizing index of abstract state    
xhat0 = sysAbs.states(:,j);
xhatsim = [xhat0];
uhat = pol(:,j);

disp([' Satisfaction probability at xhat0 = ', num2str(satProp(j))])

N = 40;

u = uhat+K*(x0-xhat0);
for i = 1:N
    w = sysLTI.mu + sqrt(sysLTI.sigma).*randn(1);
    % compute next state
    xnext = sysLTI.A*xsim(:,i)+sysLTI.B*u+sysLTI.Bw*w;
    xsim = [xsim, xnext];

    % stop if you have reached the parking region or violated the spec
    if xnext>min(P1.V) & xnext<max(P1.V)
        disp(['Reached parking area ','after ', num2str(i), ' time instances']);
        break;
    elseif xnext>min(P2.V) & xnext<max(P2.V)
        disp(['Violated specification']);
        break;
    end

    % find next abstract state, by looking at the maximizer in R wrt the value
    % function. 
  
    inR = rel.inR(xsim(:,end),sysAbs.states);
    indices_valid = indexing(inR);
    [value_max, index_aux] = max(satProp(inR)); 
    j = indices_valid(index_aux); % find maximizing index of abstract state
    xhatnext = sysAbs.states(:,j);

    % compute next input
    uhat = [pol(:,j)];
    u = uhat+K*(xnext-xhatnext);
    
    if i==N
        disp(['Did not reach parking after ', num2str(N), ' time instance.'])
    end
end


%% Show results
figure;
plot(sysAbs.states,satProp)
xlabel('x')
title('Satisfaction Probability')

ypoints = 1:length(xsim);
Plot_sysLTI(sysLTI)
hold on
plot(xsim,ypoints)
ylim([-1 length(xsim)])