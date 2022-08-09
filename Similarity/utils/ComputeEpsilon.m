function [epsilon] = ComputeEpsilon(delta,sys,Bw,mu,sigma,gridSize)
% Written by: Birgit van Huijgevoort

% inputs:
% delta: row vector with values of delta
% sys:
% ...

% output:
% epsilon: row vector with values of epsilon corresponding to delta vector

% TODO: 
% change gridSize from scalar to vector?
% choose different interface function --> adjust LMIs accordingly
% add warning messages (if delta is too big or too small?, if there is no solution?)

%% 
% system parameters
A = sys.A;
C = sys.C;

% abstraction parameters
% beta is beta_max
Bset = [gridSize/2, gridSize/2, -gridSize/2, -gridSize/2; -gridSize/2, gridSize/2, -gridSize/2, gridSize/2]; 
beta1 = Bset(:,1);
beta2 = Bset(:,2);
beta3 = Bset(:,3);
beta4 = Bset(:,4);

LengthBeta = sqrt(beta1(1)^2+beta1(2)^2);

epsGam = [];
L = [];
deleps = [];
Lameps = [];
Lambda = [];
delD = [];
lambda_v = linspace(0.9999,0,200);

for i = 1:length(delta) 
    del = delta(i);
    Eps = [];
    Lambda = [];
    DD = [];
    % Compute r
    gamma = norminv((1-del)/2,mu,sigma);
    r = abs(2*gamma);
 
    for k = 1:length(lambda_v)
        lambda = lambda_v(k);
        L = sdpvar(2,2);
        Dinv = sdpvar(2,2,'symmetric');
        eps_sq = sdpvar(1);
        M1 =    [Dinv, Dinv*C';
                C*Dinv, eye(2)];
        M2 =    [r^2*Dinv, L';
                L, eps_sq*eye(2)];
        M31 =   [lambda*Dinv,zeros(2,1), Dinv*A'+L'*Bw';
                zeros(1,2), (1-lambda)*eps_sq, -eps_sq*beta1';
                A*Dinv+Bw*L, -eps_sq*beta1, Dinv];
        M32 =   [lambda*Dinv,zeros(2,1), Dinv*A'+L'*Bw';
                zeros(1,2), (1-lambda)*eps_sq, -eps_sq*beta2';
                A*Dinv+Bw*L, -eps_sq*beta2, Dinv];
        M33 =   [lambda*Dinv,zeros(2,1), Dinv*A'+L'*Bw';
                zeros(1,2), (1-lambda)*eps_sq, -eps_sq*beta3';
                A*Dinv+Bw*L, -eps_sq*beta3, Dinv];
        M34 =   [lambda*Dinv,zeros(2,1), Dinv*A'+L'*Bw';
                zeros(1,2), (1-lambda)*eps_sq, -eps_sq*beta4';
                A*Dinv+Bw*L, -eps_sq*beta4, Dinv];
            
        obj = -eps_sq;
        LMI = [(eps_sq>=0);((1/(LengthBeta^2))-eps_sq>=0);(Dinv>=0);(M1>=0);(M2>=0);(M31>=0);(M32>=0);(M33>=0);(M34>=0)]; %eps>length(beta) 
        options = sdpsettings('verbose',0,'solver','mosek'); % mosek or sedumi
        sol = optimize(LMI,obj,options);
        
        %if sol.problem==0
        if  all([(round((1/(LengthBeta^2))-(value(eps_sq)),6)>=0);(round(eig(value(M1)),6)>=0);(round(eig(value(M2)),6)>=0);(round(eig(value(M31)),6)>=0);(round(eig(value(M32)),6)>=0);(round(eig(value(M33)),6)>=0);(round(eig(value(M34)),6)>=0)])
            eps_sq = value(eps_sq);
            eps_sq = double(eps_sq);
            
            Dinv = value(Dinv);
            D = inv(Dinv);
            
            Eps = [Eps, eps_sq];
            Lambda = [Lambda, lambda];
            DD = [DD, D];
            if length(Eps)>1
                if Eps(end)<Eps(end-1)  % maximize eps_sq
                    Eps(end) = [];
                    Lambda(end) = [];
                    DD(:,end-1:end);
                    break
                end
            end
        else
            eps_sq=[];
        end
    end
    
    eps = 1/(sqrt(Eps(end)));
    D = DD(:,end-1:end);
    deleps = [deleps; del eps];
    %Lameps = [Lameps; Lambda(end), eps];
    %delD = [delD; del*eye(2), D];
    

    if (i>2 && round(deleps(i,2)-deleps(i-1,2),6)==0)
        break;
    end

    
end

epsilon = deleps(:,2)';

end
