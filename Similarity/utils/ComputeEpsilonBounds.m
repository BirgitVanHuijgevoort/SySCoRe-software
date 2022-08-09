function [epsilonBounds] = ComputeEpsilonBounds(sysLTI,mu,sigma,Bset,varargin)
% Written by: Birgit van Huijgevoort

% inputs:
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

% Assumes interface function: u=uhat

%%

if nargin>=5
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

if nargin>=7
   if strcmp(varargin{3},'Uniform')
        distrUni = 1;
   else 
       error(['6th input argument should be Uniform or Gaussian depending ' ...
           'on the considered probability distribution of the disturbance'])
   end
else
  distrUni = 0;
end

%% 
% system parameters
A = sysLTI.A;
B = sysLTI.B;
C = sysLTI.C;
Bw = sysLTI.Bw;
dim = sysLTI.dim;

disturbSize = size(Bw, 2);

deleps = [];
Eps = [];
lambda_v = linspace(0.9999,0,100);

delta = [0, 0.999999];

options = sdpsettings('verbose',0,'solver','mosek', 'cachesolvers', 1); % mosek or sedumi

for i = 1:length(delta) 
    del = delta(i);
    Lambda = [];
    DD = [];
    % Compute r
    if distrUni
        r = 2*sqrt(disturbSize)*(1-nthroot((1-del),disturbSize));
    else
        gamma = norminv((1-del)/2,0,1);
        r = abs(2*gamma);
    end
    EpsMax = 1e-12;
    for k = 1:length(lambda_v)
        % Predefine matrices
        L = sdpvar(disturbSize,dim);
        Dinv = sdpvar(dim,dim,'symmetric');
        M1 =    [Dinv, Dinv*C';
                C*Dinv, eye(size(C,1))];
        if interfaceK
            Q = sdpvar(size(B,2),dim,'full');
        else
            Q = zeros(size(B,2),dim);
        end
        lambda = lambda_v(k);
        eps_sq = sdpvar(1);
        
        M2 =    [r^2*Dinv, L';
                L, eps_sq*eye(size(L,1))];
            
        obj = -eps_sq;
        
        if interfaceK
            M2b = [ uuf^2*Dinv, Q';
                    Q, eps_sq*eye(size(Q,1))];
            LMI = [(eps_sq>=0);(Dinv>=0);(M1>=0);(M2>=0);(M2b>=0)];
        else
            LMI = [(eps_sq>=0);(Dinv>=0);(M1>=0);(M2>=0)];
        end
            
        for l = 1:size(Bset.V,1)
            z = Bset.V(l,:)';
            
            M3 =   [lambda*Dinv,zeros(dim,1), Dinv*A'+Q'*B'+L'*Bw';
                    zeros(1,dim), (1-lambda)*eps_sq, -eps_sq*z';
                    A*Dinv+B*Q+Bw*L, -eps_sq*z, Dinv];
            LMI = [LMI; (M3>=0)];
        end
            
        sol = optimize(LMI,obj,options);
        if sol.problem==0
        %if  all([(round((1/(LengthBeta^2))-(value(eps_sq)),6)>=0);(round(eig(value(M1)),6)>=0);(round(eig(value(M2)),6)>=0);(round(eig(value(M31)),6)>=0);(round(eig(value(M32)),6)>=0);(round(eig(value(M33)),6)>=0);(round(eig(value(M34)),6)>=0)])
            eps_sq = value(eps_sq);
            eps_sq = double(eps_sq);
            
            Dinv = value(Dinv);
            D = inv(Dinv);
            K = value(Q)*D;

            if eps_sq>EpsMax  % maximize eps_sq
                EpsMax = eps_sq;
                D = value(Dinv)^-1;
                K = value(Q)*D;
            end
        else
            eps_sq=1000;
            D = eye(dim);
            K = zeros(1,dim);
        end
    end
    
    eps = 1/(sqrt(EpsMax));
    Eps = [Eps, eps];
    deleps = [deleps; del eps];
    
    if (i>2 && round(deleps(i,2)-deleps(i-1,2),6)==0)
        break;
    end
  
end

epsilonBounds = deleps(:,2)';

end
