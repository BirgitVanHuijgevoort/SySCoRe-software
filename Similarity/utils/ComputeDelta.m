function [delta, Dmin, Kmin] = ComputeDelta(epsilon,sysLTI,mu,sigma,Bset,varargin)
% Written by: Birgit van Huijgevoort
% This code is based on the method developed in the paper "Similarity
% quantification for linear stochastic systems as a set-theoretic control
% problem" by B.C. van Huijgevoort and S. Haesaert

% inputs:
% epsilon - row vector with values of epsilon
% sysLTI - LTI systems according to LinModel.m (or one partition of
% PWAModel.m) in Model classes
% mu - mean of Gaussian disturbance of original model
% sigma - variance of Gaussian disturbance of original model
% Bset - set of distubances in system x_Delta as a polyhedron

% outputs:
% delta: row vector with values of delta corresponding to epsilon vector
% Kmin: K in interface function u=uhat+K(x-xhat)
% Dmin: Weighting matrix in simulation relation R:= {||x-xhat||_D \leq \epsilon}

% Options:
% 'interface' - specify the interface function by following 'interface'
% with 0 or 1 and with an upperbound on the part of the input for K(x-xhat)
% 0. (default) u=uhat, 1. u=uhat+K(x-xhat)
%
% 'distr' - 'Gaussian' if a Gaussian distribution (default) is considered 
%               'Uniform' if a uniform distribution is considered.

%% 

%if sigma~=1
%    error('only sigma=1 allowed')
%end
%if mu~=0
%    error('only mu=0 allowed')
%end

% default values if unspecified
alpha = 1;
interfaceK = 0;
uuf = zeros(size(sysLTI.B,2),1);
distrUni = 0;

for i = 1:length(varargin)
    % try to find 'alpha'
    if strcmp(varargin{i},'alpha')
        alpha = varargin{i+1};
    end
end

for i = 1:length(varargin)
    % try to find 'interface'
    if strcmp(varargin{i},'interface')
        interfaceK = varargin{i+1};
        if interfaceK
            uuf = varargin{i+2};
        else
            uuf = zeros(size(sysLTI.B,2),1);
        end
    end
end

for i = 1:length(varargin)
    % try to find distr
    if strcmp(varargin{i},'distr')
        distribution = varargin{i+1};
        if distribution == 'Uniform'
            distrUni = 1;
        else
            distrUni = 0;
        end
    end
end

%%
% system parameters
A = sysLTI.A;
B = sysLTI.B;
C = sysLTI.C;
Bw = sysLTI.Bw;
dim = sysLTI.dim;

disturbSize = size(Bw, 2);

deleps = []; % Delta-epsilon pairs
Del = []; % Delta values
lambda_v = linspace(0.9999, 0, 100);
Lambda = []; % Corresponding lambda values
R = []; % Corresponding values of r

% Predefine matrices
L = sdpvar(disturbSize, dim);
Dinv = sdpvar(dim, dim, 'symmetric');
M1 = [  Dinv, Dinv * C'; % First LMI
        C * Dinv, eye(size(C, 1))];
r_sq = sdpvar(1); % r^2
if interfaceK
    Q = sdpvar(size(B,2),dim,'full');
else
    Q = zeros(size(B,2),dim);
end

% Solver settings
options = sdpsettings('verbose', 0, 'solver', 'mosek', 'cachesolvers', 1);
% Use either solver mosek or sedumi
% Allow the optimizer to cache the solver in order to reduce runtime
% The used solver will be stored in a persistent variable!

% Loop to compute delta for multiple epsilon values
for i = 1:length(epsilon)
    eps = epsilon(i);
    eps_sq = eps.^2;
    DelMin = 1; % Set default value for minimum delta

    % Line search
    for k = 1:length(lambda_v)
        lambda = lambda_v(k);
        
        M2 = [(1/eps_sq) * Dinv, L'; % Second LMI
            L, r_sq * eye(size(L, 1))];
        
        obj = r_sq; % Objective function
        
        if interfaceK
            M2b =   [   (1/eps^2)*Dinv, Q';
                        Q, uuf^2*eye(size(Q,1))];
            LMI = [(r_sq>=0); (Dinv>=0); (M1>=0); (M2>=0); (M2b>=0)]; % alt: r_sq>=1e-6
        else
            LMI = [(r_sq>=0); (Dinv>=0); (M1>=0); (M2>=0)]; % alt: r_sq>=1e-6
        end
        
        % Loop over vertices beta_l
        for l = 1:size(Bset.V,1)
                beta = Bset.V(l, :)';

                M3 =    [lambda*Dinv,zeros(dim,1), Dinv*A'+Q'*B'+L'*Bw';
                        zeros(1,dim), (1-lambda)*eps^2, beta';
                        A*Dinv+B*Q+Bw*L, beta, Dinv];
                    
                LMI = [LMI; (M3>=0)];
        end

        % Solve LMI
        sol = optimize(LMI, obj, options);

        if sol.problem==0 
            % Solution found
            r = sqrt(round(double(value(r_sq)), 6));
            if distrUni
                del = abs(1-(1-0.5*abs(r\sqrt(disturbSize)))^disturbSize);
            else
                del = abs(1 - 2 * normcdf(-r/2, 0, 1));
            end
            D = value(Dinv)^-1;
            K = value(Q)*D;
            if DelMin > del % update the value of D and del in case the
                % optimization improves upon the current best values
                DelMin = del;
                Dmin = value(Dinv)^-1;
                Kmin = K;
            end
        elseif sol.problem == -2 || sol.problem == -3 || sol.problem == -4 || sol.problem == -5 || sol.problem == -9 || sol.problem == -11
            % Problem
            disp('WARNING: there seems to be a problem with the solver ')
            display(sol.info)
            yalmiperror(sol.problem)
            del = 1;
            r = 1000;
            D = eye(dim);
            K = zeros(1,dim);
            break;
        else
            % Infeasible problem
            del = 1;
            r = 1000;
            D = eye(dim);
            K = zeros(1,dim);
        end
        Del = [Del, del];
        Lambda = [Lambda, lambda];
        R = [R r];
    end

    deleps = [deleps; DelMin, eps];

end

delta = deleps(:, 1)';

end
