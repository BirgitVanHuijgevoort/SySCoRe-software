clear all
close all
Install

%% Specify system parameters and regions
% LTI systems of the form
% x(t+1) = Ax(t) + Bu(t) + Bw w(t)
% y(t) = Cx(t) + Du(t)

% Define parameters
Ns = .1;
n = 1:3;
Ns_sequence = [1 Ns.^(n)./factorial(n)];
c_sequence = Ns_sequence.*0;
c_sequence(1)=1;

% Define dynamics
dim = length(Ns_sequence)-1;
A = toeplitz(c_sequence(1:end-1),Ns_sequence(1:end-1));
B = Ns_sequence(end-1:-1:1)';
C = eye(dim);
D = zeros(dim,1);
Bw = .1*eye(dim);  % 10x larger than ARCH

% Specify mean and variance of disturbance w(t)
%%% In ARCH report they use mu = [0;0] and sigma = eye(2), so not the
%%% size of dimension....
mu = zeros(dim,1);
sigma = eye(dim);

% Save all system parameters (incl. Bw) into a struct
sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);

%% Spaces and sets
% Take state space equal to safe set
xil = -12;   % Lowerbound x1
xiu = 12;   % Upperbound x1
bounds = repmat({[xil,xiu]}, n(end),1);
sysLTI.X = Polyhedron(combvec(bounds{:})');

% Single input 
ul = [-1];   % Lowerbound input u
uu = [1];     % Upperbound input u
sysLTI.U = Polyhedron([ul(1),uu(1)]');

% Specify regions for the specification
% (using same format as used in function inpolygon)
p1 = combvec(bounds{:});
P1 = Polyhedron(p1');

pil = -8;   % Lowerbound x1
piu = 8;   % Upperbound x1
boundsp = repmat({[pil,piu]}, n(end),1);
P2 = Polyhedron(combvec(boundsp{:})');

% Define the regions and the atomic propositions
sysLTI.regions = [P1, P2];
sysLTI.AP = {'p1', 'p2'};
%sysLTI.regions = [P2];
%sysLTI.AP = {'p2'};

%% Synthesize scLTL formula  
%%% use LTL2BA and check if determinstic and accepting state with loop with
% one (i.e., sink state).
% input: (sc)LTL formula and atomic propositions (see readme in folder
% LTL2BA)
% output: struct DFA containing (among others) the transitions

formula = '(p1 & X p1 & X X p1 & X X X p1 & X X X X p2) ';
%formula = '((p1 | p2) & X (p1 | p2) & X X (p1 | p2) & X X X (p1 | p2) & X X X X p2) ';
%formula = 'X X X X p2';

[DFA] = TranslateSpec(formula, sysLTI.AP);

%% Construct abstract model
disp('start gridding')
tic
ula = 1*ul;   % part of input for actuation (lowerbound)
uua = 1*uu;
ulf = ul-ula;   % part of input for feedback (lowerbound)
uuf = uu-uua;

lu = 5;
uhat = combvec(linspace(ula(1),uua(1),lu));

l = 100*ones(1,dim);  % number of grid cells in x1-, x2- and x3-direction
tol = 10^-4;

sysAbs = Gridding(sysLTI,uhat,l,tol,'TensorComputation',true, 'TensorToolbox', 'tensortoolbox');

% Label output space
label = zeros(1,prod(l));
[label(1:prod(l))] = deal(1);
inP1 = inpolygon(sysAbs.states(1,:),sysAbs.states(2,:),P1.V(:,1)',P1.V(:,2)');
inP2 = inpolygon(sysAbs.states(1,:),sysAbs.states(2,:),P2.V(:,1)',P2.V(:,2)');
[label(inP1)] = deal(3);
[label(inP2)] = deal(4);    % all states in p2 are also in p1, so we have DFA.act(4)
sysAbs.labels = label;

disp(['---> finished gridding in ', num2str(toc), ' seconds.'])
%% Compute delta based on epsilon
disp('start computing eps delta');tic;

epsilon = 2;    % should be larger than vector beta!

%[epsilonBounds] = ComputeEpsilonBounds(sysLTI,mu,sigma,sysAbs.beta)

[delta, D_m, K] = ComputeDelta(epsilon,sysLTI,sysLTI.mu,sysLTI.sigma,sysAbs.beta);

disp(['delta = ', num2str(delta), ', epsilon = ', num2str(epsilon) ])
rel = SimRel(epsilon,delta,D_m);

disp(['---> finished computing eps delta in ', num2str(toc), ' seconds'] )

%% Synthesize controller
rel.NonDetLabels  = NonDeterministicLabelling(sysAbs.states, sysLTI.regions, rel, 'Efficient', sysAbs);

tic
disp('Start computing Robust policy')
N=5; % Horizon of specification
[satProb,pol]  = SynthesizeRobustController(sysAbs,DFA, rel, N, true);
disp(['---> finished computing robust policy in ', num2str(toc), ' seconds.'] )

%% Plot value function
if dim == 2
    figure;
    [X1hat, X2hat] = ndgrid(sysAbs.hx{1},sysAbs.hx{2});
    surf(X1hat, X2hat,reshape(satProb,l(1),l(2)),'EdgeColor','interp')
    xlabel('x_1')
    ylabel('x_2')
    title('Value function')
end

%% Simulation
if dim == 3
    x0 = [-5;0;6];
    xsim = x0;
    indexing = 1:length(sysAbs.states);

    if P2.contains(x0)
        disp('Initial state is inside target area')
    elseif sysLTI.X.contains(x0)
        disp('Initial state is inside safe area')
    else
        disp('Initial state is not in safe area... please specify a different initial state')
    end

    % find initial abstract state
    diff = abs(x0.*ones(size(sysAbs.states))-sysAbs.states);
    inR = (([1 1 1]*((D_m^.5*diff).^2)).^.5)<=epsilon;
    indices_valid = indexing(inR);
    [value_max, index_aux] = max(satProb(inR)); 
    j = indices_valid(index_aux); % find maximizing index of abstract state    
    xhat0 = sysAbs.states(:,j);
    xhatsim = [xhat0];
    uhat = pol(:,j);

    disp([' Satisfaction probability at xhat0 = ', num2str(satProb(j))])

    N = 6;

    u = uhat+K*(x0-xhat0);
    for i = 1:N
        w1 = sysLTI.mu(1) + sqrt(sysLTI.sigma(1,1)).*randn(1);
        w2 = sysLTI.mu(2) + sqrt(sysLTI.sigma(2,2)).*randn(1);
        w3 = sysLTI.mu(3) + sqrt(sysLTI.sigma(3,3)).*randn(1);
        w = [w1;w2;w3];
        % compute next state
        xnext = sysLTI.A*xsim(:,i)+sysLTI.B*u+sysLTI.Bw*w;
        xsim = [xsim, xnext];

        if P2.contains(xnext)
            disp('State is inside target area')
        elseif sysLTI.X.contains(xnext)
            disp('State is inside safe area')
        else
            disp('Specification violated!')
        end

        % find next abstract state, by looking at the maximizer in R wrt the value
        % function. 
        inR = rel.inR(xsim(:,end),sysAbs.states);
        indices_valid = indexing(inR);
        [value_max, index_aux] = max(satProb(inR)); 
        j = indices_valid(index_aux); % find maximizing index of abstract state
        xhatnext = sysAbs.states(:,j);

        % compute next input
        uhat = [pol(:,j)];
        u = uhat+K*(xnext-xhatnext);
    end
end