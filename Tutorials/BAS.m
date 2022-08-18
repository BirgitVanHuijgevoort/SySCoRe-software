%%% Building_automation case study
% sysLTI -> original model
% sysLTIr -> reduced order model
%
%% Case study: Model with large number of continuous variables
% author: Birgit
% author of 7D model: Nathalie Cauchi
% -------------------------------------------------------
% x_c[k+1] = A_cx_c[k] + B_cu_c[k] +F_cd_c[k] + Q_c
% y_c[k]   = [1 0 0 0 0 0 0]
% x_c = [T_z1 T_z2 T_w5 T_w6 T_w2 T_w3 T_w7]^T
% u_c = T_sa
% d_c =[T_out T_hall CO2_1 CO2_2 T_rw,r1 T_rw,r2]^T
% -------------------------------------------------------
% -------------------------------------------------------

clear;close all;

run Install.m

tStart = tic;

% Load parameters needed to build model
BASmodel = load('BASmodel.mat');
Z1m = BASmodel.Z1m;
Qc = BASmodel.Qc;

%% Specify parameters of 7-dimensional model
A = Z1m.A;
B = Z1m.B;
C = Z1m.C;
D = zeros(1,1);
Bw = Z1m.F(:,1:6);
dim = Z1m.dim;

mu = [9;15;500;500;35;35]; % mean of disturbance
sigma = [1, zeros(1,5); 0, 1, zeros(1,4); 0, 0 , 100, zeros(1,3);
                0,0,0,100,zeros(1,2); zeros(1,4), 5, 0;
                zeros(1,5), 5];% variance of disturbance

% Transform model, such that w comes from Gaussian distribution with mean 0
% and variance eye().
Bw = Bw*sigma^(1/2);
Qc = Qc+Bw*mu;
mu = zeros(6,1);
sigma = eye(6);

% input bounds
ul = 15;
uu = 30;

% Work with system with state x-xss and output y-C*xss
sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);

% Remove Qc by looking at steady state solution
uss = (uu-ul)/2+ul;
uss = 15;
xss = -inv(A-eye(dim))*(B*uss+Qc);

% New bounds on input space (u-uss)
ul = ul-uss;     % Lowerbound input u
uu = uu-uss;     % Upperbound input u
sysLTI.U = Polyhedron(combvec([ul(1),uu(1)])');

% Specify division of input space for actuation and feedback
ula = -0.6*((uu-ul)/2)+((uu-ul)/2)+ul;   % part of input for actuation R*uhat, with R=1
uua = 0.6*((uu-ul)/2)+((uu-ul)/2)+ul;
ulf = -0.175*((uu-ul)/2);   % part of input for feedback K(x-P*xhat)
uuf = 0.175*((uu-ul)/2);

% New bounds on state space (x-xss)
x1l = 19.5-xss(1);   % Lowerbound x1
x1u = 20.5-xss(1);   % Upperbound x1
x2l = 19-xss(2);   % Lowerbound x2
x2u = 22-xss(2);   % Upperbound x2
x3l = 18-xss(3);   % Lowerbound x3
x3u = 22-xss(3);   % Upperbound x3
x4l = 18-xss(4);   % Lowerbound x4
x4u = 22-xss(4);   % Upperbound x4
x5l = 18-xss(5);   % Lowerbound x5
x5u = 22-xss(5);   % Upperbound x5
x6l = 18-xss(6);   % Lowerbound x6
x6u = 22-xss(6);   % Upperbound x6
x7l = 18-xss(7);   % Lowerbound x7
x7u = 22-xss(7);   % Upperbound x7
Xbound = [x1l, x1u; x2l, x2u; x3l, x3u; x4l, x4u; x5l, x5u; x6l, x6u; x7l x7u];
sysLTI.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u],[x3l,x3u],[x4l,x4u],[x5l,x5u],[x6l,x6u],[x7l,x7u])');

% Specify regions for the specification
P1 = Polyhedron(combvec([x1l,x1u],[-100,100],[-100,100],[-100,100],[-100,100],[-100,100],[-100,100])');

%% Reduced order model
f = 0.098;
dimr = 2;
[sysLTIr,F] = ModelReduction(sysLTI,dimr,f);

% Make sure that the reduced model has first state as output (C=[1 0])
if sysLTIr.C(1)~= 1
    % flip outputs
    T2  = [0 1; 1 0];
    
    % set to 1
    T=eye(length(sysLTIr.C));
    T(length(sysLTIr.C),length(sysLTIr.C))= sysLTIr.C(end);
    
    % state transformation
    sysLTIr.A = T2*T*sysLTIr.A*inv(T)*inv(T2);
    sysLTIr.C = sysLTIr.C*inv(T)*inv(T2);
    BBw= T2*T*[sysLTIr.B, sysLTIr.Bw];
    sysLTIr.B = BBw(:,1);
    sysLTIr.Bw = BBw(:,2:end);
end

% Compute P and Q
[P, Q, R] = ComputeProjection(sysLTI,sysLTIr);

% Bounds on state space
x1l = x1l;   % Lowerbound x1
x1u = x1u;   % Upperbound x1

% Compute state space bounds of reduced-order model
Ax = [];
bx = [];
for i = 1:2^dim
    xv = sysLTI.X.V(i,:)';

    for j = 1:dim
        if P(j,2)>=0
            Axi(2*(j-1)+1,:) = [0 -P(j,2)];
            Axi(2*j,:) = [P(j,2) 0];
            bxi(2*(j-1)+1,1) = [-xv(j)+P(j,1)*x1u];
            bxi(2*j,1) = [xv(j)-P(j,1)*x1l];
        else
            Axi(2*(j-1)+1,:) = [0 P(j,2)];
            Axi(2*j,:) = [-P(j,2) 0];
            bxi(2*(j-1)+1,1) = [xv(j)-P(j,1)*x1u];
            bxi(2*j,1) = [-xv(j)+P(j,1)*x1l];
        end
    end

    Ax = [Ax; Axi];
    bx = [bx; bxi];
end

xhatBound = linprog([-1;1], Ax,bx);
if ~isempty(xhatBound)
    x2l = xhatBound(1);   
    x2u = xhatBound(2);   
end

sysLTIr.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');

sysLTIr.U = sysLTI.U;

%% Specify output region
% Specify regions for the specification
TSP = 20-xss(1);               % temperature set point
p1x = [19.5-xss(1) 20.5-xss(1)  20.5-xss(1) 19.5-xss(1) 19.5-xss(1)];    % x1-coordinates
p1y = [-100 -100 100 100 -100];    % x2-coordinates
p1 = [p1x; p1y];        % goal region
P1r = Polyhedron(p1x(1:2)');

sysLTIr.regions = [P1r]; % regions that get specific atomic propositions
sysLTIr.AP = {'p1'}; % with the corresponding atomic propositions

%% Synthesize scLTL formula (or input DFA yourself)
N = 6;
formula = '(p1 & X p1 & X X p1 & X X X p1 & X X X p1 & X X X X p1 & X X X X X p1)';
[DFA] = TranslateSpec(formula, sysLTIr.AP);

%% Reduce State space
[~, output_set] = IncreaseDecreasePolytope(P1r, 0.1);
[sysLTIr,~] = ReduceX(sysLTIr, [ula(1),uua(1)], output_set, 'invariance', 5);

%% Construct abstract model by gridding it
disp('start gridding');tic

lu = 3;
uhat = combvec(linspace(ula(1),uua(1),lu));

l = [3000*3000];  % number of grid cells in x1- and x2-direction
tol=10^-6;
sysAbs = Gridding(sysLTIr,uhat,l,tol,'TensorComputation',true);
l = sysAbs.l;% load updated l
% Save some extra system parameters into struct
sysAbs.orig = sysLTIr;

label = zeros(1,prod(l));
[label(1:prod(l))] = deal(1);
inP1 = inpolygon(sysAbs.states(1,:),sysAbs.states(2,:),p1x,p1y);
[label(inP1)] = deal(2);
sysAbs.labels = label;

toc
disp(['---> finished gridding in ', num2str(toc), ' seconds.'])
%% Compute delta based on epsilon

% Compute polytope 
beta = sysAbs.beta;
Uhat = Polyhedron(sysAbs.inputs');

Wlb = sysLTIr.mu-3*sum(sysLTIr.sigma,2);
Wub = sysLTIr.mu+3*sum(sysLTIr.sigma,2);
Wset = Polyhedron('lb', Wlb, 'ub', Wub);

% Compute additional error, by truncating the disturbance
onemindel = mvncdf(Wlb,Wub,mu,sigma);
del_trunc = 1-onemindel;

Z = (sysLTI.B*R-P*sysLTIr.B)*Uhat+(sysLTI.Bw-P*sysLTIr.Bw)*Wset;
Zred = Z;
Zred = Z.minVRep();

% [OPTIONAL] This function computes the bounds for epsilon (takes a lot of time)
% epsilonBounds = [eps_max,eps_min]
% [epsilonBounds] = ComputeEpsilonBoundsMOR(sysLTI,sysLTIr,mu,sigma,uuf,Zred,P)
epsilonBounds = [0.2413, 0.0520]; 

epsilon_1 = epsilonBounds(1);
if epsilon_1>=epsilonBounds(2) && epsilon_1<=epsilonBounds(1)
    disp('Feasible epsilon chosen')
elseif epsilon_1>epsilonBounds(1)
    disp('Epsilon is larger then necessary')
elseif epsilon_1<=epsilonBounds(2)
    disp('Infeasible epsilon chosen')
end

disp('start computing simulation relation');tic

% Compute MOR simulation relation
[delta_1, D_1, K_1, F_1] = ComputeDelta_intPQRK(epsilon_1,sysLTI,sysLTIr,mu,sigma,uuf,Zred,P);
del = delta_1;
delta_1 = delta_1+del_trunc;

% Compute gridding simulation relation
epsilon_2 = max(0.35-epsilon_1,0);
[delta_2, D_2, K_2] = ComputeDelta(epsilon_2,sysLTIr,sysLTIr.mu,sysLTIr.sigma,beta);

delta = delta_1+delta_2;
epsilon = epsilon_1+epsilon_2;

rel_1 = SimRel(epsilon_1,delta_1,D_1);
rel_2 = SimRel(epsilon_2,delta_2,D_2);
rel = rel_1.Combine(rel_2,sysLTIr.X);
disp(['delta = ', num2str(rel.delta), ', epsilon = ', num2str(rel.epsilon) ])

toc
disp('---> finished computing simulation relation')

rel.NonDetLabels  = NonDeterministicLabelling(sysAbs.outputs, sysLTIr.regions, rel);
%rel.NonDetLabels  = NonDeterministicLabelling(sysAbs.outputs, sysLTIr.regions, rel, 'Efficient', sysAbs);


%% Synthesize controller
disp('---> Start Value iteration')
[V,pol]  = SynthesizeRobustController(sysAbs,DFA,rel,N,false);

q = DFA.S0;
% Fix Value at initial q0 based on labeling
for i = 1:size(sysAbs.labels,2)
    if sysAbs.labels(i) == 2
        V(q,i) = V(q,i);
    else
        V(q,i) = 0;
    end
end

% Determine computation time
tEnd = toc(tStart);

X1 = reshape(sysAbs.states(1,:),l)+xss(1);
X2 = reshape(sysAbs.states(2,:),l)+xss(2);

VC = reshape(V(q,:),l);

% Plot satisfaction probability
figure;
surf(X1(1:15:end,1:15:end),X2(1:15:end,1:15:end),VC(1:15:end,1:15:end),'EdgeColor','interp')
xlabel('x_{r1}', 'FontSize', 16)
ylabel('x_{r2}', 'FontSize', 16)
view(2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

%% Start simulation
YSIM = [];
USIM = [];
[i,j] = max(V(DFA.S0,:));
xr0 = sysAbs.states(:,j);
x0 = P*xr0;

% run the simulation 6 times
for l = 1:6 
    xsim = [x0];
    xrsim = [xr0];

    q = DFA.S0;
        if P1r.contains(sysLTIr.C*xr0)
            label = 2;
        else
            label = 1;
        end
    q_next = DFA.trans(q,label);

    % find initial abstract state, look at xr0
    diff = abs(xr0.*ones(1,size(sysAbs.states,2))-sysAbs.states);
    dis = sqrt(sum(diff.^2,1));
    [~,j] = min(dis);

    xhat0 = sysAbs.states(:,j);
    xhatsim = [xhat0];

    % satisfaction probability of this initial state
    SatProp = V(q,j);
    disp(['Value function at x0 equals ', mat2str(SatProp,4),])

    uhat = pol(:,j, q_next);
    indexing = 1:length(sysAbs.states);

    ursim = [];
    usim = [];
    q = q_next;
    for i = 1:N
        disp(['x = ', mat2str(xsim(:,i)+xss,4), ',    q^+ = ',num2str(q_next)])    % compute next continuous state

        if q ==DFA.F
            disp("satisfied specification")
            break;
        elseif q == DFA.sink
            disp("failed specification")
            break;
        end
        w = mvnrnd(sysLTIr.mu,sysLTIr.sigma)';

        u = R*uhat + Q*xrsim(:,i) + K_1*(xsim(:,i)-P*xrsim(:,i));
        ur = uhat;
        ursim = [ursim, ur];
        usim = [usim, u];

        xnext = sysLTI.A*xsim(:,i)+sysLTI.B*u+sysLTI.Bw*w;
        xsim = [xsim, xnext];

        wr = w+F_1*(xsim(:,i)-P*xrsim(:,i));
        xrnext = sysLTIr.A*xrsim(:,i)+sysLTIr.B*ur+sysLTIr.Bw*wr;
        xrsim = [xrsim, xrnext];

        if P1.contains(xsim(:,end)) 
            label = 2;
        elseif xsim(1,end)>=19.5-xss(1) && xsim(1,end)<=20.5-xss(1)
            label = 2;
        else
            label = 1;
        end
        q_next = DFA.trans(q,label);

        % find next abstract state, by looking at the maximizer in R wrt the value
        % function. 

        diff2 = abs(xrsim(:,end).*ones(1,size(sysAbs.states,2))-sysAbs.states);
        inR2 = (([1 1]*((D_2^.5*diff2).^2)).^.5)<=epsilon_2;
        inR = inR2;
        indices_valid = indexing(inR);
        V_q = V(q_next,:);
        [value_max, index_aux] = max(V_q(inR));
        j = indices_valid(index_aux); % find maximizing index of abstract state
        xhatnext = sysAbs.states(:,j);
        xhatsim = [xhatsim, xhatnext];

        uhat = [pol(:,j,q_next)];

        q = q_next;

    end
    
    % compute output
    y = C*xsim+C*xss;
    usim_ss = usim+uss;
    
    YSIM = [YSIM; y];
    USIM = [USIM; usim_ss];
    
end

%% plot results

figure;
subplot(2,1,1)
plot(y(1,:))
hold on
plot(y(1,:),'bo')
plot(1:N,19.5*ones(1,N),'r--')
plot(1:N,20.5*ones(1,N),'r--')
title('State evolution')

subplot(2,1,2)
plot(usim+uss);
hold on
plot(usim+uss,'bo')
title('Input')

k = 0:5;
figure;
plot(k,YSIM')
hold on
plot(0:N-1,19.5*ones(1,N),'r--')
plot(0:N-1,20.5*ones(1,N),'r--')
title('State evolution')

%% Display execution time 
% end time is before simulation part
disp(['Total runtime = ', mat2str(tEnd)])

%% Verify input bound
Qx = Q*sysLTIr.X;
Qx = Qx.minVRep();

u_lb = ula+ulf+min(Qx.V);
u_ub = uua+uuf+max(Qx.V);

% use only (abstract) states with V>0
[~,j] = find(V(DFA.S0,:)>0);
XX = sysAbs.states(:,j);

u_lb2 = ula+ulf+min(Q*XX);
u_ub2 = uua+uuf+max(Q*XX);

disp(['Input is between ', mat2str(u_lb2+uss,4), ' and ', mat2str(u_ub2+uss,4)])