function [sysPWA] = PWAapprox(sysNonLin,T,N,NL)
% PWAapproximation computes a PWA approximation of the system.
% Here, the deterministic part of the system should be of the form x(t+1) =
% f(x,u)
%
% Outputs:
% sysPWA = piecewise-affine system that approximation the nonlinear system
% sysNonLin
%
% Inputs:
% sysNonLin = nonlinear system
% T = symbolic expression for the linear taylor approximation
% N = number of partitions in each direction [N1 N2 N3 ...]
% NL = vector with NL(i) = 1 if i-th function in sysNonLin.fsym(i) is
% nonlinear
%
% Copyright 2022 Birgit van Huijgevoort b.c.v.huijgevoort@tue.nl

%% Initialization

nx = sysNonLin.dim;
nu = sysNonLin.U.Dim;

% Check if partitions are square, if not, make square
if ~all(N == N(1))
    warning('Current implementation requires square partitions for the PWA approximation, however, you have not supplied this. I will change this for you.')
    Nmax = max(N);
    N = Nmax.*ones(1,nx+nu);
end

% Check if only nonlinear wrt x or only wrt u
LinearInput = 0;
LinearState = 0; 
if all(NL(:,2) == 0) % linear wrt input!
    LinearInput = 1;
    % remove input partitioning from N
    if size(N,2) > nx
        warning('The systems function is affine wrt the input, so no need to partition the input space. I changed this for you.')
        N = N(1,1:nx);
    end
end
if all(NL(:,1) == 0) % linear wrt states!
    LinearState = 1;
    % remove input partitioning from N
    if size(N,2) > nu
        warning('The systems function is affine wrt the states, so no need to partition the state space. I changed this for you.')
        N = N(1,nx+1:nx+nu);
    end
end

%% Create partitioning of the state space and input space

if LinearInput
    XUbounds = [min(sysNonLin.X.V)' max(sysNonLin.X.V)'];
elseif LinearState
    XUbounds = [min(sysNonLin.U.V)' max(sysNonLin.U.V)'];
else
    Xbounds = [min(sysNonLin.X.V)' max(sysNonLin.X.V)'];
    Ubounds = [min(sysNonLin.U.V)' max(sysNonLin.U.V)'];
    XUbounds = [Xbounds;Ubounds];
end

% % Create points of combined state and input space (if necessary)
for i = 1:size(N,2)
    points = linspace(XUbounds(i,1),XUbounds(i,2),N(i));
    Xsp(:,i) = points;
end

% Create square partitions of specified points
v = ones(1,size(XUbounds,1));
ready = false;
iC = 1;
while ~ready
    
    for i = 1:size(N,2)
        tmp{i} = Xsp(v(i):v(i)+1,i)';
    end
    
    Partition(iC).Polyhedron = Polyhedron([combvec(tmp{:})']);
    iC = iC+1;
    
    % Shift the index vector:
    ready = true;       % Claim temporarily that the WHILE loop is ready
    for k = 1:size(N,2)
        v(k) = v(k) + 1;
        if v(k) <= N(k)-1
            ready = false;
            break;         % v(k) increased successfully, leave "for k" loop
        end
        v(k) = 1;         % v(k) reached nx+nu: reset it, proceed with next v
    end
end

Np = size(Partition,2);   % number of partitions

%% Define dynamics (using taylor expansion)
x = sym('x', [nx 1]);
u = sym('u', [nu 1]);

if all(all(NL == 0))
    disp('This system is already linear!')
else
    [i, ~] = find(NL == 1);
    i = unique(i);
    fprintf('Performing piecewise-affine approximation for state update of state %d', i)
end

fsym = sysNonLin.fsym;

if all(all(NL == 0))   % all functions are linear, so no need to perform PWA approximation, just rewrite dynamics to matrix form.
    for k = 1:Np
        for i = 1:length(fsym)
            for j = 1:nx
                tmp = coeffs(fsym(i), x(j)); 
                if has(tmp(end),x) == 0
                    A(i,j) = tmp(end);
                end 
            end
            for j = 1:nu
                tmp = coeffs(fsym(i), u(j));
                if has(tmp(end),u) == 0
                    B(i,j) = tmp(end);
                end
            end
            tmp = coeffs(fsym(i),'All');
            a(i,1) = tmp(end);
        end

        Partition(k).Dynamics.A = A;
        Partition(k).Dynamics.a = a;
        Partition(k).Dynamics.B = B;

        Partition(k).Dynamics.C = sysNonLin.C;  % seems redundant, but necessary for ComputeDelta function
        Partition(k).Dynamics.D = 0;
        Partition(k).Dynamics.Bw = sysNonLin.Bw;  % seems redundant, but necessary for ComputeDelta function
        Partition(k).Dynamics.dim = sysNonLin.dim;
    end
else
    for k = 1:Np
        xc = Partition(k).Polyhedron.chebyCenter;
        Partition(k).Dynamics.A = T.A(xc.x');
        Partition(k).Dynamics.a = T.a(xc.x');   
        Partition(k).Dynamics.B = T.B(xc.x');

        Partition(k).Dynamics.C = sysNonLin.C;  % seems redundant, but necessary for ComputeDelta function
        Partition(k).Dynamics.D = 0;
        Partition(k).Dynamics.Bw = sysNonLin.Bw;  % seems redundant, but necessary for ComputeDelta function
        Partition(k).Dynamics.dim = sysNonLin.dim;
    end
end
C = sysNonLin.C;    % Assume that this part is linear!
Bw = sysNonLin.Bw;
mu = sysNonLin.mu;
sigma = sysNonLin.sigma;

sysPWA = PWAModel(Partition,C,Bw,mu,sigma);
end