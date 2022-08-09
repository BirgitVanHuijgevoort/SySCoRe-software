clear all
clc
close all
%% Van der Pol tutorial
run Install.m

tStart = tic;
disp(' Van der Pol Oscillator tutorial ')

% Set value for epsilon
epsilon = 0.1;

% Load model into sysNonLin
Vanderpol_dist

% Bounds on state space 
x1l = -4;   % Lowerbound x1
x1u = 4;   % Upperbound x1
x2l = -4;   % Lowerbound x2
x2u = 4;   % Upperbound x2
sysNonLin.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');

% Bounds on  and input space
ul = [-1];   % Lowerbound input u
uu = [1];     % Upperbound input u
sysNonLin.U = Polyhedron([ul(1),uu(1)]');

% Specify division of input space for actuation and feedback
% if PWAsystem is stable
ula = 0.6*ul;   % part of input for actuation (lowerbound)
uua = 0.6*uu;
ulf = ul-ula;   % part of input for feedback (lowerbound)
uuf = uu-uua;
% if PWAsystem is unstable
ula_2 = 0.6*ul;  
uua_2 = 0.6*uu;
ulf_2 = ul-ula_2;   
uuf_2 = uu-uua_2;

% Stay in bounded region
P1 = Polyhedron(combvec([-3,3],[-3,3])');
%P1 = Polyhedron(combvec([-5,5],[-5,5])');

% Reach in bounded region
%P2 = Polyhedron(combvec([2,3],[-1,1])');
P2 = Polyhedron(combvec([-1.2,-0.9],[-2.9,-2])');

sysNonLin.regions = [P1;P2]; % regions that get specific atomic propositions
sysNonLin.AP = {'p1', 'p2'}; % with the corresponding atomic propositions

Plot_sysLTI(sysNonLin)

% Compute symbolic expression for Taylor expansion and remainder 
[T,R] = ComputeTaylorAndRemainder(sysNonLin);

%% 1. Synthesize scLTL formula (or input DFA yourself)
%%% use LTL2BA and check if determinstic and accepting state with loop with 1.
% input: (sc)LTL formula and atomic propositions (see readme in folder
% LTL2BA)
% output: struct DFA containing (among other) the transitions

formula = '(p1 U p2)';  % p1 = safe region, p2 = target region
% formula should use atomic propositions in sysLTI.AP. 

% Make sure your current folder is the main SySCoRe folder
[DFA] = TranslateSpec(formula,sysNonLin.AP);

%% 2. Construct abstract model

%% Abstraction part 1 PWA approximation 

showDynamics = 0;   % Set to 1 if you want to see the dynamics of the VanderPol oscillator and its PWA approximation
showAnalysis = 1;

if showDynamics
    % plot dynamics of autonomous system
    Xt = zeros(2,100);
    x1 = -1;
    x2 = 1.5;
    Xt(:,1) = [x1;x2];
    for t = 1:500
        x1 = Xt(1,t);
        x2 = Xt(2,t);
        x1next = x1+x2*param.tau;
        x2next = x2+(-x1+(1-x1^2)*x2)*param.tau;

        Xt(:,t+1) = [x1next;x2next];
    end

    figure(10);
    plot(Xt(1,1),Xt(2,1),'r*')
    hold on
    plot(Xt(1,:),Xt(2,:))
    xlabel('x1')
    ylabel('x2')
    title('Dynamics of nominal autonomous system')

    % quiver plot of dynamics
    [X1,X2] = meshgrid(x1l:.25:x1u, x2l:.25:x2u);
    U = zeros(size(X1));
    V = zeros(size(X2));
    for i = 1:size(X1,1)
        for j = 1:size(X1,2)
            x1 = X1(i,j);
            x2 = X2(i,j);
            x1next = x1+x2*param.tau;
            x2next = x2+(-x1+(1-x1^2)*x2)*param.tau;

            u = x1next-x1;
            v = x2next-x2;

            U(i,j) = u/sqrt(u^2+v^2);
            V(i,j) = v/sqrt(u^2+v^2);
        end
    end       
    quiver(X1,X2,U,V,0.5, 'color',[0 0 1]);
end

%%% Define partitions and dynamics of PWA approximation
% N = [36,36]; % number of grid points in each direction @ todo put this back
N = [4,4]; % number of grid points in each direction

[sysPWA] = PWAapproximation(sysNonLin,T,N,showDynamics);
sysPWA.X = sysNonLin.X;
sysPWA.U = sysNonLin.U;
sysPWA.regions = sysNonLin.regions;
sysPWA.AP = sysNonLin.AP;

% Change division of input depending on stability of PWA system
for k = 1:sysPWA.Np
    if sum(abs(eig(sysPWA.Partition(k).Dynamics.A))>=1)>0   % system is unstable
        sysPWA.Partition(k).Uu = max(abs(ulf_2),abs(uuf_2));
    else
        sysPWA.Partition(k).Uu = max(abs(ulf),abs(uuf));
    end
end

%%% Quantify difference 

% Only works for f(x) being a polynomial function.
% p = polynomialDegree(sysNonLin.fsym);
for k = 1:sysPWA.Np
    sysPWA.Partition(k).K = ComputeTaylorAccuracy(sysPWA.Partition(k),R);
end

%%% Analysis of K per partition
if showAnalysis
    figure(10);
    for j = 1:sysPWA.Np
        xc = sysPWA.Partition(j).Polyhedron.chebyCenter;
        K = sysPWA.Partition(j).K;
        text(xc.x(1),xc.x(2), num2str(K.V(1,2),2))
    end
end

%% Abstraction part 2 gridding
disp('start gridding');tic
% input: 
% output: 
lu = 5;
%uhat = linspace(ul(1),uu(1),lu);
uhat = linspace(ula(1),uua(1),lu);
[~,idx]=sort(abs(uhat));
uhat=uhat(idx);
% l = [450, 450];  % number of grid cells in x1- and x2-direction
l = [100, 100];
%l = [20, 20];
tol = 10^-8;

% Grid system with efficient nonlinear gridding
sysAbs = GridSpace_nonlin_tensor_v2(sysPWA,uhat,l,tol, 'tensortoolbox');

% Save some extra system parameters into struct
sysAbs.orig = sysNonLin;

disp('start labelling');tic

label = zeros(1,prod(l));
[label(1:prod(l))] = deal(1);
index_labels = sysNonLin.regions.contains(sysAbs.states);
[label(index_labels(1,:))] = deal(3);
[label(index_labels(2,:))] = deal(4);
sysAbs.labels = label;

toc
disp('---> finished gridding')

disp('determine partition for each abstract state');tic
% Determine for each abstract state in which partition it lies
for k = 1:sysPWA.Np
     sysAbs.Partition(find(sysPWA.Partition(k).Polyhedron.contains(sysAbs.states))) = k;
end
toc

% Limit input uhat by changing P_det --> VERIFY IF CORRECT?!?
RemoveInput = sysAbs.inputs>uua_2 | sysAbs.inputs<ula_2;
[~,RemoveInput] = find(RemoveInput);
NrOfStates = size(sysAbs.states,2);
for i = 1:length(RemoveInput)
    Col = [(RemoveInput(i)-1)*NrOfStates+1:RemoveInput(i)*NrOfStates];
    sysAbs.P.P_det(:,Col) = sparse(1,1);
end


%% 3. Compute delta based on epsilon
disp('start computing eps delta');tic;

% rewrite to be able to use parfor instead of for
for i = 1:sysPWA.Np
    Kset(i) = sysPWA.Partition(i).K;
    Dynamics(i) = sysPWA.Partition(i).Dynamics;
    Uu(:,i) = sysPWA.Partition(i).Uu;
end

% Compute epsilon bounds for biggest Bset! (hence for partition with
% largest K)
%Kmax = 0;
%for i = 1:sysPWA.Np
%    K = max(max(Kset(i).V));
%    if K>Kmax
%        Kmax = K;
%        index = i;
%    end
%end

mu = sysPWA.mu;
sigma = sysPWA.sigma;
Beta = sysAbs.beta;
%Bset = plus(Kset(index),Beta);
%Bset = minVRep(Bset);
%[epsilonBounds] = ComputeEpsilonBounds(Dynamics(index),mu,sigma,Bset,2,Uu(:,index),'Uniform')

f_delta = ones(1,sysPWA.Np);
DD = [];
KfK = [];
parfor i = 1:sysPWA.Np
%for i = 1:sysPWA.Np
    Bset = plus(Kset(i),Beta);
    Bset = minVRep(Bset);
    
    %@Sofie
    [delta, D, Kf] = ComputeDelta(epsilon,Dynamics(i),mu,sigma,Bset,2,Uu(:,i),'Uniform');
    disp(i)
    disp(delta)
%     [delta, D, Kf] = ComputeDelta(epsilon,Dynamics(i),mu,sigma,Bset,2,Uu(:,i));
    f_delta(1,i) = delta;
    DD = [DD, D];
    KfK = [KfK; Kf];
end
for i = 1:sysPWA.Np
    sysPWA.Partition(i).delta = f_delta(1,i);
    sysPWA.Partition(i).D = DD(:,(i-1)*dim+1:i*dim);
    sysPWA.Partition(i).Kf = KfK(i,:);
    sysPWA.Partition(i).rel = SimRel(epsilon,sysPWA.Partition(i).delta,eye(2));
end

disp([', epsilon = ', num2str(epsilon), 'delta = ', num2str(f_delta)])

rel = SimRel(epsilon,f_delta,eye(2));

rel.NonDetLabels  = NonDeterministicLabelling(sysAbs.outputs, sysNonLin.regions, rel);
toc    
disp('---> finished computing eps delta')

%%% Analysis of delta per partition
if showDynamics
    figure(10);
    for j = 1:sysPWA.Np
        xc = sysPWA.Partition(j).Polyhedron.chebyCenter;
        text(xc.x(1),xc.x(2)-0.2, num2str(f_delta(j),2))
        %text(xc.x(1),xc.x(2)-0.25, num2str(Kset(j).V(1,2),2))
    end
end
if showAnalysis    
    X = [];
    Y = [];
    F_delta = [];
    KK = [];
    for i = 1:size(sysAbs.states,2)
        x = sysAbs.states(1,i);
        y = sysAbs.states(2,i);
        delta = f_delta(sysAbs.Partition(:,i));
        K = sysPWA.Partition(sysAbs.Partition(:,i)).K.V(1,2);

        X = [X x];
        Y = [Y y];
        F_delta = [F_delta, delta];
        KK = [KK, K];
    end

    X = reshape(X,sqrt(size(sysAbs.states,2)),sqrt(size(sysAbs.states,2)));
    Y = reshape(Y,sqrt(size(sysAbs.states,2)),sqrt(size(sysAbs.states,2)));
    F_delta = reshape(F_delta,sqrt(size(sysAbs.states,2)),sqrt(size(sysAbs.states,2)));
    KK = reshape(KK,sqrt(size(sysAbs.states,2)),sqrt(size(sysAbs.states,2)));

    figure;
    surf(X,Y,F_delta,'EdgeColor','interp');
    %zlim([0 1]);
    xlabel('\hat{x}_1', 'FontSize', 20)
    ylabel('\hat{x}_2', 'FontSize', 20)
    zlabel('f_{\delta}(\hat{x})','FontSize', 20)
    zticks([0.01 0.02])
    xlim([-4,4]);
    ylim([-4,4]);
    %ax = gca;
    %ax.FontSize = 16; 
    
    figure;
    surf(X,Y,KK,'EdgeColor','interp');
    xlabel('x1')
    ylabel('x2')
    zlabel('\kappa_{i,max}')
end


%% Synthesize controller
disp('start computing robust controller')

N = 40;     % time horizon

% Get delta in correct shape for SynthesizeRobustController
rel.delta = rel.delta(sysAbs.Partition)';

[satProp,pol] = SynthesizeRobustController(sysAbs, DFA, rel, N, true);

[X1hat, X2hat] = ndgrid(sysAbs.hx{1},sysAbs.hx{2});

figure;
surf(X1hat, X2hat,reshape(satProp,l(1),l(2)),'EdgeColor','interp')
xlabel('x_1', 'FontSize', 20)
ylabel('x_2', 'FontSize', 20)
xlim([-4,4]);
ylim([-4,4]);
colorbar

[satProp,pol] = SynthesizeRobustController(sysAbs,DFA, rel, N, false);

disp('---> finished computing robust controller')

tEnd = toc(tStart);

%% Display execution time 
% end time is before simulation part
disp(['Total runtime = ', mat2str(tEnd)])

%% Simulation 
x0 = [-1;1.5];

xsim = [x0];

% first step in DFA
q = DFA.S0;
if P2.contains(x0)
    label = 2;
elseif P1.contains(x0)
    label = 3;
else
    label = 1;
end
q_next = DFA.trans(q,label);

% initial abstract state
diff = x0-sysAbs.states;
dis = sqrt(sum(diff.^2,1));
[~,j] = min(dis);
xhat0 = sysAbs.states(:,j);
xhatsim = [xhat0];

% satisfaction probability of initial state
SatProp = satProp(q,j);
disp(['Value function at x0 equals ', mat2str(SatProp,4),])

% initial input
uhat = pol(:,j, q_next);

% run simulation 
usim = [];
uhatsim = [];
indexing = 1:length(sysAbs.states);
q=q_next;
for i = 1:N
    disp(['x = ', mat2str(xsim(:,i),4), ',    q^+ = ',num2str(q_next)])    % show next continuous state

    if q == DFA.F
        disp("satisfied specification")
        break;
    elseif q == DFA.sink
        disp("failed specification")
        break;
    end
    w = mvnrnd(sysNonLin.mu,sysNonLin.sigma)';

    Pr = sysAbs.Partition(j);     % Partition nr of abstract state
    
    Kf = sysPWA.Partition(Pr).Kf;
    u = uhat+Kf*(xsim(:,end)-xhatsim(:,end));
    usim = [usim, u];
    
    A = sysPWA.Partition(Pr).Dynamics.A;
    B = sysPWA.Partition(Pr).Dynamics.B;
    a = sysPWA.Partition(Pr).Dynamics.a;
    Bw = sysPWA.Partition(Pr).Dynamics.Bw;
    xnext = A*xsim(:,end)+B*u+Bw*w;

    xsim = [xsim,xnext];

    if P2.contains(xsim(:,end))
        label = 2;
    elseif P1.contains(xsim(:,end))
        label = 3;
    else
        label = 1;
    end
    q_next = DFA.trans(q,label);

    % find next abstract state, by looking at the maximizer in R wrt the value function. 
    D = sysPWA.Partition(Pr).D;
    diff2 = xsim(:,end)-sysAbs.states;
    inR = (([1 1]*((D^.5*diff2).^2)).^.5)<=epsilon;
    indices_valid = indexing(inR);
    satProp_q = satProp(q_next,:);
    [value_max, index_aux] = max(satProp_q(inR));
    j = indices_valid(index_aux); % find maximizing index of abstract state
    xhatnext = sysAbs.states(:,j);
    xhatsim = [xhatsim, xhatnext];

    uhat = [pol(:,j,q_next)];
    uhatsim = [uhatsim, uhat];
    
    q = q_next;
    
    if size(j,2)==0
        q = DFA.sink;
    end
    
end

%% Show results
figure;
plot_x = plot(sysNonLin.X);
set(plot_x, 'FaceColor','None');
plot_y = plot(sysNonLin.C*sysNonLin.X);
set(plot_y, 'FaceColor','None');
hold on 
for i =1:length(sysNonLin.regions)
    plot_x = plot(sysNonLin.regions(i));
    set(plot_x, 'FaceColor','None');
    hold on
    xc = sysNonLin.regions(i).chebyCenter;
    text(xc.x(1),xc.x(2), sysNonLin.AP{i})
end
xlim([x1l x1u])
ylim([x2l x2u])
xlabel('x_1')
ylabel('x_2')
plot(xsim(1,:),xsim(2,:))
plot(xsim(1,:),xsim(2,:),'rx')
