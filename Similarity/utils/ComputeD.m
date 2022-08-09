function [D, K] = ComputeD(eps,sysPWA,mu,sigma,Bset_all,index,varargin)
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
% varargin{1} - number (1 or 2) that determines the interface function.
% varargin{2} - upperbound on the part of the input for K(x-xhat)
% varargin{3} - 'Gaussian' if a Gaussian distribution (default) is considered 
%               'Uniform' if a uniform distribution is considered.

% Interface function: 1. (default) u=uhat, 2. u=uhat+K(x-xhat)

% output:
% delta: row vector with values of delta corresponding to epsilon vector

%% 

%if sigma~=1
%    error('only sigma=1 allowed')
%end
%if mu~=0
%    error('only mu=0 allowed')
%end

interfaceK = 1;

if nargin>=7
   if varargin{1} == 2
        interfaceK = 1;
        uuf = varargin{2};
   elseif varargin{1} == 1
       interfaceK = 0;
       uuf = zeros(size(sysLTI.B,2),1);
   else 
       error(['6th input argument should be 1 or 2 depending ' ...
           'on the desired interface function'])
   end
else
  interfaceK = 0;
  uuf = zeros(size(sysLTI.B,2),1);
end
%%
dim = sysPWA.dim;
Dinv = sdpvar(dim, dim, 'symmetric');
LMI2 = [];

%%
for i = 1:length(index)
    PN = index(i);
    Bset = Bset_all(PN);
    
    % system parameters
    A = sysPWA.Partition(PN).Dynamics.A;
    B = sysPWA.Partition(PN).Dynamics.B;
    C = sysPWA.Partition(PN).Dynamics.C;
    Bw = sysPWA.Partition(PN).Dynamics.Bw;

    disturbSize = size(Bw, 2);

    lambda_v = linspace(0.9999, 0, 100);
    Lambda = []; % Corresponding lambda values
    R = []; % Corresponding values of r

    % Predefine matrices
    L = sdpvar(disturbSize, dim);
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
            
            
            del = abs(1 - 2 * normcdf(-r/2, 0, 1));
            
            if DelMin > del % update the value of D and del in case the
                % optimization improves upon the current best values
                DelMin = del;
                LambdaMin(i) = lambda;
                LMIMin = LMI;
            end
        elseif sol.problem == -2 || sol.problem == -3 || sol.problem == -4 || sol.problem == -5 || sol.problem == -9 || sol.problem == -11
            % Problem
            disp('WARNING: there seems to be a problem with the solver ')
            display(sol.info)
            yalmiperror(sol.problem)
            del = 1;
            r = 1000;
            break;
        else
            % Infeasible problem
            del = 1;
            r = 1000;
        end
    end

    assert(exist("LMIMin", "var"), "No feasible solution to LMIs found.")
    LMI2 = [LMI2; LMIMin];
  
end

sol = optimize(LMI2, obj, options);
D = value(Dinv)^-1;
K = value(Q)*D;
r = sqrt(round(double(value(r_sq)), 6));
delta = abs(1 - 2 * normcdf(-r/2, 0, 1));

end
