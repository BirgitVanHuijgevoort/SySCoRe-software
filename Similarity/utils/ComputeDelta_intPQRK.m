function [delta, Dmin, Kmin, Fmin] = ComputeDelta_intPQRK(epsilon,sysLTI,sysLTIr,mu,sigma,uuf,Z,P)
% Written by: Birgit van Huijgevoort
% Difference wrt ComputeDelta_intPQRK: +Bw*what-Bw*what instead of PBrw*...
% This code is based on the method developed in the paper "Similarity
% quantification for linear stochastic systems as a set-theoretic control
% problem" by B.C. van Huijgevoort and S. Haesaert

% inputs:
% epsilon: row vector with values of epsilon
% sys:
% ...

% Interface function: u=R*uhat+Q*xhat+K*(x-P*xhat)

% output:
% delta: row vector with values of delta corresponding to epsilon vector

% This code only works if the disturbance has a Gaussian distibrution with
% mean 0 and variance identity!

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
    
%if length(gridSize)==1
%    gridSize = gridSize*ones(size(dim,1));
%end

% abstraction parameters
Zset = Z.V; 

deleps = [];
Del = [];
lambda_v = linspace(0.9999,0,100);
Lambda = [];
R = [];
Dmin = eye(dim); % default values in case a solution does not exist
Kmin = ones(size(B,2),dim);

options = sdpsettings('verbose',0,'solver','mosek', 'cachesolvers', 1); % mosek or sedumi

for i = 1:length(epsilon)
    eps = epsilon(i);
    DelMin = 1;
    for k = 1:length(lambda_v)
        lambda = lambda_v(k);
        L = sdpvar(size(Bw,2),dim);
        E = sdpvar(size(B,2),dim,'full');
        Dinv = sdpvar(dim,dim,'symmetric');
        r_sq = sdpvar(1);
        M1 =    [Dinv, Dinv*C';
                C*Dinv, eye(size(C,1))];
        M2 =    [(1/eps^2)*Dinv, L';
                L, r_sq*eye(size(L,1))];
        M2b =   [(1/eps^2)*Dinv, E';
                E, uuf^2*eye(size(E,1))];
       
        obj = r_sq;
        LMI = [(r_sq>=0);(Dinv>=0);(M1>=0);(M2>=0);(M2b>=0)];   
        
        for j = 1:size(Zset,1)
                z = Zset(j,:)';
  
                M3 =    [lambda*Dinv,zeros(dim,1), Dinv*A'+E'*B'+L'*Bw';
                        zeros(1,dim), (1-lambda)*eps^2, z';
                        A*Dinv+B*E+Bw*L, z, Dinv];

                LMI = [LMI;(M3>=0)];
        end

        sol = optimize(LMI,obj,options);
        if sol.problem==0 
            r_sq = double(value(r_sq));
            r_sq = round(r_sq,6);
            r = sqrt(r_sq);
            
            del = abs(1-2*normcdf(-r/2,0,1));
            D = value(Dinv)^-1;
            K = value(E)*D;
            F = value(L)*D;
            if DelMin > del % update the value of D and del in case the
                % optimization improves upon the current best values
                DelMin = del;
                Dmin = value(Dinv)^-1;
                Kmin = value(E)*D;
                Fmin = value(L)*D;
            end
        elseif sol.problem == -2 || sol.problem == -3 || sol.problem == -4 || sol.problem == -5 || sol.problem == -9 || sol.problem == -11
            disp('WARNING: there seems to be a problem with the solver ')
            display(sol.info)
            yalmiperror(sol.problem)
            del = 1;
            r = 1000;
            D = eye(dim);
            K = eye(dim);
            F = eye(dim);
            break;
        else
            del = 1;
            r = 1000;
            D = eye(dim);
            K = eye(dim);
            F = eye(dim);
        end
        Del = [Del, del];
        Lambda = [Lambda, lambda];
        R = [R r];
    end

    deleps = [deleps; DelMin, eps];

end

delta = deleps(:,1)';


end
