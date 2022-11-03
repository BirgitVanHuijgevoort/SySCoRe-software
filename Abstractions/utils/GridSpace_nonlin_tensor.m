function [sysAbs] = GridSpace_nonlin_tensor(sys, uhat,l, tol, tensortool)
% [sysAbs] = GridSpace(sys, uhat,l, tol, tensortool) 
%
% Inputs: 
% sys   = LTI systems with fields A and B 
% uhat  = Finite set of inputs,  example uhat = combvec(linspace(ul(1),uu(1),3),linspace(ul(2),uu(2),3));
% l     = Number of finite states in each dimension  [l1 l2 l3 ...]
% tol  =  tolerance for truncating to 0
% tensortool = "tensortoolbox" or "tensorlab" or "2d"
%  
% Outputs:
% P = "matrix" describing the transition probabilities. 
% P(i,j,k) is the probability of going from state i to state j with
% input uhat(:,k)
% output 2: abstract state space
%
%% 

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

if nargin < 5 
    tensortool = "tensortoolbox"; 
end

    
if length(l)==1
    uniform=0;
else
    uniform =1;
end
if length(sigma)==1
    sigma = diag(sigma*ones(1,dim));
elseif size(sigma,1) == 1
    sigma = diag(sigma);
elseif size(sigma,2) == 1
    sigma = diag(sigma);
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


if ~uniform
Zsize = (Zu-Zl);
lx = floor(sqrt(l*Zsize(1)/Zsize(2)));
ly = floor(l/lx);
l=[lx ly];
    
end
 %===================== Compute uniform grid==========================
% width of single slot
gridSize = (Zu-Zl)./l; 
% selection of representative points in each dimension 
% using a cell function
% hx = Xl+gridSize/2:gridSize:Xu;
hz = arrayfun(@(zl, zu, gridsize) zl+gridsize/2:gridsize:zu, ...
    Zl, Zu, gridSize, 'UniformOutput', false);
ZhatSpace = combvec(hz{:});

XhatSpace = Uz*ZhatSpace;
%===================== Compute deterministic probability matrix==========================
% P_det is going to be a sparse matrix. 
nXhatSpace = size(XhatSpace, 2);
nUhat = size(uhat, 2);
index_total = 1:length(ZhatSpace);
j_indices = [];
i_indices = [];
for k = 1:nUhat % for each control action
    z_n = Uz'*sys.f_det(XhatSpace, repmat(uhat(:, k), 1, nXhatSpace));
    z_n_ind = ones(dim, 1) + floor(diag(gridSize.^-1)*(z_n-Zl'));
    indices = min( ...
        [1 <= z_n_ind; ...
        l' >= z_n_ind], ...
        [], 1);
    zi_indices= arrayfun(@(i) z_n_ind(i, indices), ...
        1:dim, 'UniformOutput', false);
    i_indices = [i_indices, sub2ind(l, zi_indices{:})];
    j_indices = [j_indices, index_total(indices) + (k-1)*nXhatSpace];
end
% Probability of travelling to state i under state-input pair j for
% deterministic transition
P_det = sparse(i_indices, j_indices, ones(size(i_indices)), ...
    nXhatSpace, nXhatSpace*nUhat);

%===================== Compute stochastic probability matrix==========================
% P_stoch
% compute extended grid (2x the size) around the origin. 
mx2  = arrayfun(@(gridsize,l_i) -gridsize*0.5:gridsize:gridsize*(l_i-0.5), ...
    gridSize, l, 'UniformOutput', false);


P = cell(dim,1);
for k = 1:nUhat
    for d_index = 1:dim
        cp =  diff(normcdf(mx2{d_index},  mu_z(d_index), sigma_z(d_index)));
        cp(cp<tol) = 0; % equivalent to (t1(2:length(t1))-t1(1:length(t1)-1));
        P{d_index} = toeplitz(cp);
    end
end

if strcmp(tensortool,'2d')
   P = TensorTransitionProbability(l,P_det,P{:});
elseif strcmp(tensortool,'tensorlab' )
   P = TensorTransitionProbability_tensorlab(l,P_det,P{:});
elseif strcmp(tensortool,'tensortoolbox')
   P = TensorTransitionProbability_tensortoolbox(l,P_det,P{:});
end

beta = Uz*Polyhedron((diag(2*gridSize)*(ff2n(dim)-0.5)')'); 
% !!!! gridsize is twice the size that you would have with a direct gridding. 
states = XhatSpace;
sysAbs = MDP_model(P,hz,states,beta, sys);
sysAbs.zstates = ZhatSpace;
sysAbs.outputs = sys.C*XhatSpace;
sysAbs.inputs = uhat; 
sysAbs.l = l;
sysAbs.outputmap = sys.C*Uz; 



