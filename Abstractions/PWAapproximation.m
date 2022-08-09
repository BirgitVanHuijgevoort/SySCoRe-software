function [sysPWA] = PWAapproximation(sysNonLin,B,T,N,NL)
%Function that creates a PWA approximation of the system.
% Here, the deterministic part of the system should be of the form x(t+1) = f(x)x + Bu 

n = sysNonLin.dim;

%% Create partitioning of the state space
% Xbounds = zeros(n,2); % Xbounds contains lower and upper-bound of state space
% Xsp = cell(1,n);
% for i = 1:n
%     Xbounds(i,1) = min(sysNonLin.X.V(:,i));
%     Xbounds(i,2) = max(sysNonLin.X.V(:,i));
% 
%     points = linspace(Xbounds(i,1),Xbounds(i,2),N(i));
%     Xsp(1,i) = {points};
% end
% 
% Xpoints = combvec(Xsp{:})';
% 
% [X1,X2] = ndgrid(Xpoints(:,1)',Xpoints(:,2)');
% 
% % Write something similar as ndgrid yourself. Use repmat, reshape to get
% % grid points?!?!
% 
% % ------ Write to nD! --------
% % N has dimension n
% % Xpoints has n columns
% % for i = 1:N(1)-1
% %     for j = 1:N(2)-1
% %         Partition((i-1)*(N(1)-1)+j).Polyhedron = Polyhedron([combvec(Xpoints(j:j+1,1), Xpoints(i:i+1,2))]');
% %     end
% % end
% 
% for i = 1:prod(N)
%     Partition(i).Polyhedron = Polyhedron('lb', [X1(i,i); X1(i,i+1)], 'ub', [X1(i+1,i); X1(i+1,i+1)]);
% end
% 
% Np = size(Partition,2);   % number of partitions

Xbounds = [min(sysNonLin.X.V)' max(sysNonLin.X.V)'];
for i = 1:n
    points = linspace(Xbounds(i,1),Xbounds(i,2),N(i));
    Xsp(:,i) = points;
end

v = ones(1,n);
ready = false;
iC = 1;
while ~ready
    
    for i = 1:n
        tmp{i} = Xsp(v(i):v(i)+1,i)';
    end
    
    Partition(iC).Polyhedron = Polyhedron([combvec(tmp{:})']);
    iC = iC+1;
    
    % Shift the index vector:
    ready = true;       % Claim temporarily that the WHILE loop is ready
    for k = 1:n
        v(k) = v(k) + 1;
        if v(k) <= N(k)-1
            ready = false;
            break;         % v(k) increased successfully, leave "for k" loop
        end
        v(k) = 1;         % v(k) reached nx: reset it, proceed with next v
    end
end

Np = size(Partition,2);   % number of partitions

%% Define dynamics (using taylor expansion)
x = sym('x', [n 1]);

if all(NL == 0)
    disp('This system is already linear!')
else
    i = find(NL == 1);
    fprintf('Performing piecewise-affine approximation for state update of state %d', i)
end

fsym = sysNonLin.fsym;

if all(NL==0)   % both functions are linear, so no need to perform PWA approximation, just rewrite dynamics to matrix form.
    for k = 1:Np
        for i = 1:length(fsym)
            for j = 1:n
                tmp = coeffs(fsym(i), x(j)); 
                if has(tmp(end),x) == 0
                    A(i,j) = tmp(end);
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
        
        Partition(k).Dynamics.B = B;
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