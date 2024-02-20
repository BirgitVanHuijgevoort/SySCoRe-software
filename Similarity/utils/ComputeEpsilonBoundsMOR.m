function [epsilonBounds] = ComputeEpsilonBoundsMOR(sysLTI,sysLTIr,mu,sigma,uuf,Z)
% Written by: Birgit van Huijgevoort

% inputs:
% sys:
% ...

% output:
% epsilon: row vector with values of epsilon_max and epsilon_min

% Assumes interface function: u=R*uhat+Q*xhat+K*(x-P*xhat)

%% 
% system parameters
A = sysLTI.A;
B = sysLTI.B;
C = sysLTI.C;
Bw = sysLTI.Bw;
dim = sysLTI.dim;

Bwr = sysLTIr.Bw;

if isequal(mu,zeros(size(Bwr,2),1)) && isequal(sigma,eye(size(Bwr,2)))
else
    error('mu different from zero and sigma different from identity has not been implemented yet') 
end

% abstraction parameters
Zset = Z.V; 

deleps = [];
Eps = [];
lambda_v = linspace(0.9999,0,100);

delta = [0, 0.999999];

for i = 1:length(delta) 
    del = delta(i);
    Lambda = [];
    DD = [];
    % Compute r
    gamma = norminv((1-del)/2,0,1);
    r = abs(2*gamma);
    EpsMax = 1e-12;
    for k = 1:length(lambda_v)
        lambda = lambda_v(k);
        eps_sq = sdpvar(1);
        L = sdpvar(size(Bw,2),dim);
        E = sdpvar(size(B,2),dim,'full');
        Dinv = sdpvar(dim,dim,'symmetric');
        
        M1 =    [Dinv, Dinv*C';
                C*Dinv, eye(size(C,1))];
        M2 =    [r^2*Dinv, L';
                L, eps_sq*eye(size(L,1))];
        M2b =   [uuf^2*Dinv, E';
                E, eps_sq*eye(size(E,1))];
            
        obj = -eps_sq;
        LMI = [(eps_sq>=0);(Dinv>=0);(M1>=0);(M2>=0);(M2b>=0)];
        
        for j = 1:size(Zset,1)
            z = Zset(j,:)';
            
            M3 =   [lambda*Dinv,zeros(dim,1), Dinv*A'+E'*B'+L'*Bw';
                    zeros(1,dim), (1-lambda)*eps_sq, -eps_sq*z';
                    A*Dinv+B*E+Bw*L, -eps_sq*z, Dinv];
            LMI = [LMI; (M3>=0)];
        end
            
        options = sdpsettings('verbose',0,'solver','mosek'); % mosek or sedumi
        sol = optimize(LMI,obj,options);
        if sol.problem==0
        %if  all([(round((1/(LengthBeta^2))-(value(eps_sq)),6)>=0);(round(eig(value(M1)),6)>=0);(round(eig(value(M2)),6)>=0);(round(eig(value(M31)),6)>=0);(round(eig(value(M32)),6)>=0);(round(eig(value(M33)),6)>=0);(round(eig(value(M34)),6)>=0)])
            eps_sq = value(eps_sq);
            eps_sq = double(eps_sq);
            
            Dinv = value(Dinv);
            D = inv(Dinv);
            K = value(E)*D;

            if eps_sq>EpsMax  % maximize eps_sq
                EpsMax = eps_sq;
                %D = value(Dinv)^-1;
            end
        else
            eps_sq=1000;
            D = eye(dim);
            K = eye(dim);
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
