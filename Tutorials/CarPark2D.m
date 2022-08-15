clc
clear
close all

tStart = tic;

%% Specify system parameters and regions
% LTI systems of the form
% x(t+1) = Ax(t) + Bu(t) + Bw w(t)
% y(t) = Cx(t) + Du(t)

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

sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);
 
% Bounds on state space 
x1l = -10;   % Lowerbound x1
x1u = 10;   % Upperbound x1
x2l = -10;   % Lowerbound x2
x2u = 10;   % Upperbound x2
sysLTI.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');

% Bounds on  and input space
ul = [-1;-1];   % Lowerbound input u
uu = [1;1];     % Upperbound input u
sysLTI.U = Polyhedron(combvec([ul(1),uu(1)],[ul(2),uu(2)])');

% Specify regions for the specification
p1x = [4 4 10 10 4];    % x1-coordinates
p1y = [0 -4 -4 0 0];    % x2-coordinates
p1 = [p1x; p1y];         % parking region
P1 = Polyhedron(p1');

p2x = [4 4 10 10 4];    % x1-coordinates
p2y = [0 4 4 0 0];      % x2-coordinates
p2 = [p2x; p2y];        % avoid region
P2 = Polyhedron(p2');

sysLTI.regions = [P1;P2]; % regions that get specific atomic propositions
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

lu = 3;
uhat = combvec(linspace(ula(1),uua(1),lu),linspace(ula(2),uua(2),lu));

l = [200, 200];  % number of grid cells in x1- and x2-direction
tol=10^-6;
sysAbs = Gridding(sysLTI,uhat,l,tol,'TensorComputation',true);

% Save some extra system parameters into struct
sysAbs.orig = sysLTI;

label = zeros(1,prod(l));
[label(1:prod(l))] = deal(1);
inP1 = inpolygon(sysAbs.states(1,:),sysAbs.states(2,:),p1x,p1y);
inP2 = inpolygon(sysAbs.states(1,:),sysAbs.states(2,:),p2x,p2y);
[label(inP1)] = deal(3);
[label(inP2)] = deal(2);
sysAbs.labels = label;

tGridEnd = toc(tGridStart);
disp(['---> finished gridding in ', num2str(tGridEnd), ' seconds.'])

%% Quantify similarity
disp('Quantify similarity');
tSimStart = tic;

epsilon = 1.005;

[rel, K] = QuantifySim(sysLTI, sysAbs, epsilon, sysLTI.mu,sysLTI.sigma,sysAbs.beta);

tSimEnd = toc(tSimStart);

disp(['delta = ', num2str(rel.delta), ', epsilon = ', num2str(rel.epsilon) ])
disp(['---> finished quantifying similarity in ', num2str(tSimEnd), ' seconds'] )

%% Synthesize controller
disp('start computing robust controller')
N = 150;     % time horizon

[satProp,pol] = SynthesizeRobustController(sysAbs, DFA, rel, N, true);

tEnd = toc(tStart);
disp(['Total runtime = ', mat2str(tEnd)])

%% Start simulation
x0 = [-4;-5];
xsim = x0;
indexing = 1:length(sysAbs.states);

% find initial abstract state
diff = abs(x0.*ones(size(sysAbs.states))-sysAbs.states);
inR = (([1 1]*((rel.R^.5*diff).^2)).^.5)<=epsilon;
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
    w1 = sysLTI.mu(1) + sqrt(sysLTI.sigma(1,1)).*randn(1);
    w2 = sysLTI.mu(2) + sqrt(sysLTI.sigma(2,2)).*randn(1);
    w = [w1;w2];
    % compute next state
    xnext = sysLTI.A*xsim(:,i)+sysLTI.B*u+sysLTI.Bw*w;
    xsim = [xsim, xnext];

    % stop if you have reached the parking region
    if inpolygon(xnext(1),xnext(2),p1x,p1y)
        disp(['Reached parking area ','after ', num2str(i), ' time instances']);
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
        disp(['Did not reach parking after ', N, ' time instance.'])
    end
end


%% Show results
figure;
[X1hat, X2hat] = ndgrid(sysAbs.hx{1},sysAbs.hx{2});

surf(X1hat, X2hat,reshape(satProp,l(1),l(2)),'EdgeColor','interp')
xlabel('x_1')
ylabel('x_2')
title('Value function')

figure(10);
plot(p1x,p1y,'g+-','LineWidth',2)
axis equal
hold on
plot(p2x,p2y,'r+-','LineWidth',2)
xlim([x1l x1u])
ylim([x2l x2u])
xlabel('x_1')
ylabel('x_2')
plot(xsim(1,:),xsim(2,:))
plot(xsim(1,:),xsim(2,:),'rx')
