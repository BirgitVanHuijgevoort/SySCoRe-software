function [sysAbs] = GridSpace_nd(sys, uhat,l, tol)
% [P,X1hat,X2hat,XhatSpace,gridSize] = GridSpace(sys, uhat,l) 
%
% Inputs: 
% sys   = LTI systems with fields A and B 
% uhat  = Finite set of inputs,  example uhat = combvec(linspace(ul(1),uu(1),3),linspace(ul(2),uu(2),3));
% l     = Number of finite states in each dimension  [l1 l2 l3 ...]
%
% 
% Outputs:
% sysabs with
% P = matrix describing the transition probabilities.
% P = [Pu1, Pu2, ...] is a concatenation of P-matrices for each uhat
% P(i,j, 1:l) is the probability of going from state i to state j with
% input uhat(:,1)
%
%% 
try
    A = sys.A;
catch
    error('sys.A does not exist')
end
try
    B = sys.B;
catch 
    error('sys.B does not exist') 
end
try
    Bw = sys.Bw;
catch
    error('sys.Bw does not exist')
end
try
    sigma = sys.sigma;
catch
    error('sys.sigma does not exist')
end
try
    dim = sys.dim;
catch
    error('sys.dim does not exist')
end
try
    mu = sys.mu;
catch
    mu = zeros(dim,1);
end
try
    sys.X.computeVRep;
    Xl = min(sys.X.V);
    Xu = max(sys.X.V);
catch 
     error('sys.X does not exist or is not a Polyhedron')
end 
%if Bw ~= eye(dim) 
%    error('Bw different from identity has not been implemented yet') 
%end
if dim > 2
    error('Dimension different from 2 has not been implemented yet') 
end 
if length(l)==1
    l = l*ones(1,dim);
end
if length(sigma)==1
    sigma = sigma*ones(1,dim);
end
if length(mu)==1
    mu = mu*ones(1,dim);
end

%==================== Transform state space ========================
if isequal(Bw,eye(dim))
    Uz = eye(dim);
    sigma_z = diag(sigma);
    mu_z = eye(dim)*sys.mu;
else
    Sigma = sys.Bw*sigma*sys.Bw';
    [Uz,Sz,~] = svd(Sigma);
    sigma_z = diag(Sz);
    mu_z = Uz'* sys.Bw*sys.mu;
end

if sum(abs(mu_z))>0
warning('Implementation does not hold for mu not equal to zero')
end

% zstates = Uz'*xstates
% xstates = Uz*zstates
try
Z = Uz'*sys.X;
catch 
     error('sys.X does not exist or is not a Polyhedron')
end

Z.computeVRep;
Zl = min(Z.V);
Zu = max(Z.V);

%===================== Marginals for uniform grid==========================
% width of single slot
gridSize = (Xu-Xl)./l;   

% selection of representative points in each dimension 
% using a cell function
% hx = Xl+gridSize/2:gridSize:Xu;
hx  = arrayfun( @(xl,xu,gridsize)  xl+gridsize/2:gridsize:xu, Xl,Xu,gridSize,'UniformOutput',false);

% intervals 
% mx = Xl:gridSize:Xl+l*gridSize ;        
mx  = arrayfun( @(xl,xu,gridsize, l_i)  xl:gridsize:xl+l_i*gridsize, Xl,Xu,gridSize,l,'UniformOutput',false);

XhatSpace = combvec(hx{:});

P = zeros(prod(l),prod(l),length(uhat)); 
for k = 1:length(uhat)
    for index = 1:length(XhatSpace)
            xhat =XhatSpace(:,index);
            pij_ = 1;
            for d_index = 1:dim
               cp =  diff(normcdf(mx{d_index}, A(d_index,:)*xhat+B(d_index,:)*uhat(:,k)+mu_z(d_index), sigma_z(d_index)));
               cp(cp<tol) = 0; % equivalent to (t1(2:length(t1))-t1(1:length(t1)-1));
               pij_ = reshape(pij_'*cp,1,[]);
            end
            P(index,:,k) = pij_;
    end
end
P = TransitionProbability(P);
states = XhatSpace;
beta = Polyhedron((diag(gridSize)*(ff2n(dim)-0.5)')');
sysAbs = MDP_model(P,hx,states,beta,sys);
sysAbs.inputs= uhat;          




%==========================================================================